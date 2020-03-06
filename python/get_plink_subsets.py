#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 09:17:42 2020

Get PLINK subsets for clumping

@author: nbaya
"""

import hail as hl
import argparse

hl.init(log='/tmp/hail.log')

bucket = 'gs://ukb-diverse-pops'
REFERENCE_GENOME = 'GRCh37'
CHROMOSOMES = list(map(str, range(1, 23))) + ['X', 'XY']
POPS = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']

MIN_CASES = 50
MIN_CASES_ALL = 100
MIN_CASES_EUR = 100

ldprune_wd = f'{bucket}/ld_prune' 

pop_dict = {'AFR': 6637, # dict with sample counts for each population 
            'AMR': 982,
            'CSA': 8876,
            'EAS': 2709,
            'EUR': 420542,
            'MID': 1599}


#rev_pop_dict = {}
#for k, v in sorted(list(pop_dict.items()), key=lambda x:x[0].lower(), reverse=True):
#    rev_pop_dict[k] = v

# from https://github.com/atgu/ukbb_pan_ancestry/blob/master/resources/generic.py
#   - get_hq_samples()
#   - get_pruned_tsv_path()
#   - get_age_sex_tsv_path()
#   - get_covariates_ht_path()
#   - get_covariates()

def get_hq_samples():
    ht = hl.import_table(f'{bucket}/misc/ukb31063_samples_qc_FULL.txt', no_header=True)
    drop_samples = hl.import_table(f'{bucket}/misc/ukb31063.withdrawn_samples_20190321.txt', no_header=True, key='f0')
    ht = ht.key_by(s=ht.f0).drop('f0')
    return ht.filter(hl.is_missing(drop_samples[ht.s]))


def get_pruned_tsv_path():
    return f'{bucket}/pca/ukb_diverse_pops_pruned.tsv.bgz'


def get_age_sex_tsv_path():
    return f'{bucket}/Phenotypes/uk_round2_allSamples_phenos_phesant.6148_5.tsv.gz'


def get_covariates_ht_path():
    return f'{bucket}/pca/all_pops_non_eur_pruned_within_pop_pc_covs.ht'


def get_covariates(key_type = hl.str):
    ht = hl.read_table(get_covariates_ht_path())
    return ht.key_by(s=key_type(ht.s))

# from https://github.com/atgu/ukbb_pan_ancestry/blob/master/resources/genotypes.py
#   - get_sample_file()
#   - get_ukb_imputed_data()
#   - get_filtered_mt()
#   - get_ukb_af_ht_path()

ukb_imputed_bgen_path = 'gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chr{}_v3.bgen'
ukb_imputed_info_path = 'gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_mfi_chr{}_v3.txt'
ukb_imputed_info_ht_path = f'{bucket}/imputed/ukb_mfi_v3.ht'

def get_sample_file(chromosome: str = '1'):
    if chromosome not in ('X', 'XY'):
        chromosome = 'autosomes'
    elif not chromosome.startswith('chr'):
        chromosome = f'chr{chromosome}'
    return f'gs://ukb31063/ukb31063.{chromosome}.sample'

def get_ukb_imputed_data(chromosome: str = '1', variant_list: hl.Table = None, entry_fields = ('GP', )):
    if chromosome == 'all':
        chromosome = '{' + ','.join(map(str, range(1, 23))) + '}'
    add_args = {}
    if variant_list is not None:
        add_args['variants'] = variant_list
    return hl.import_bgen(ukb_imputed_bgen_path.format(chromosome), entry_fields=entry_fields,
                          sample_file=get_sample_file(chromosome), **add_args)

def get_mt(chrom: str = 'all', pop: str = 'all', imputed: bool = True, min_mac: int = 20, 
                    entry_fields = ('GP', ), not_pop: bool = False):
    if imputed:
        ht = hl.read_table(get_ukb_af_ht_path())
        if pop == 'all':
            ht = ht.filter(hl.any(lambda x: ht.af[x] * ht.an[x] >= min_mac, hl.literal(POPS)))
        else:
            ht = ht.filter(ht.af[pop] * ht.an[pop] >= min_mac)
        mt = get_ukb_imputed_data(chrom, variant_list=ht, entry_fields=entry_fields)
    else:
        mt = hl.read_matrix_table('gs://ukb31063/ukb31063.genotype.mt')

    covariates_ht = get_covariates()
    hq_samples_ht = get_hq_samples()
    # TODO: confirm that this is correct set
    mt = mt.annotate_cols(**covariates_ht[mt.s])
    mt = mt.filter_cols(hl.is_defined(mt.pop) & hl.is_defined(hq_samples_ht[mt.s]))

    return mt

def get_ukb_af_ht_path(with_x = True):
    return f'{bucket}/imputed/ukb_frequencies{"_with_x" if with_x else ""}.ht'



def get_pop_prop_dict(pop_dict: dict, pop: str, not_pop: bool = False) -> (dict, int):
    r'''
    Get population proportions in `pop_dict` for a population `pop`.
    If `not_pop`=True, this gets population proportions for every population
    except population `pop`.
    '''
    if not_pop:
        not_pop_dict = {k:v for k,v in pop_dict.items() if k!=pop}
        n_total = sum(not_pop_dict.values())
        pop_prop_dict = {k: v/n_total for k,v in not_pop_dict.items()}
    else:
        pop_prop_dict = {pop: 1.0}
        n_total = pop_dict[pop]
    return pop_prop_dict, n_total

def get_subset(mt_pop, pop_dict: dict, pop: str, n_max: int, not_pop: bool = False):
    r'''
    Get Hail table sample of max size = `n_max` for population `pop`.
    If `not_pop`=True, this gets the sample subset with every population except
    population `pop`.
    '''
    pop_prop_dict, n_total = get_pop_prop_dict(pop_dict=pop_dict, 
                                               pop=pop, 
                                               not_pop=not_pop)

    limiting_pop = min(pop_prop_dict, key=pop_prop_dict.get)
    n_sample = int(min(pop_dict[limiting_pop]/pop_prop_dict[limiting_pop], n_max))
    if n_sample != n_max:
        print(f'Using sample size of {n_sample} instead of {n_max} due to limiting population size in {limiting_pop}')
    print({k:v*n_sample for k,v in pop_prop_dict.items()})
        
    cols = mt_pop.cols()
    if not not_pop and n_sample == pop_dict[pop]: # if sampling a single population `pop` and n_sample is the same as the population's size
        ht_sample = cols
    else:
        cols = cols.annotate(tmp_rand = hl.rand_norm())
        cols = cols.order_by('tmp_rand')
        cols = cols.add_index(name = 'rand_idx')
        ht_sample = cols.filter(cols.rand_idx<n_sample)
        ht_sample = ht_sample.drop('tmp_rand','rand_idx')
    ht_sample = ht_sample.key_by('s')
    ht_sample = ht_sample.select('pop') # keyed by 's', thus the two remaining fields are 'pop' and 's'
    
    return ht_sample
    
def to_plink(pop: str, mt, ht_sample, not_pop: bool = False, overwrite=False):
    r'''
    Exports matrix table to PLINK2 files
    '''
    assert 'GT' in mt.entry, "mt must have 'GT' as an entry field"
    assert mt.GT.dtype==hl.tcall, "entry field 'GT' must be of type `Call`"
    bfile_path = f'{ldprune_wd}/subsets/{"not_" if not_pop else ""}{pop}'
    
    if not overwrite and not all([hl.hadoop_exists(f'{bfile_path}.{suffix}') for suffix in ['bed','bim','fam']]):
        print(f'\nPLINK files already exist for {"not_" if not_pop else ""}{pop}')
    else:
        mt_sample = mt.filter_cols(hl.is_defined(ht_sample[mt.s]))
        mt_sample = mt_sample.annotate_rows(varid = mt_sample.locus+mt_sample.alleles[0]+mt_sample.alleles[1])
        hl.export_plink(dataset = mt_sample, 
                        output = bfile_path, 
                        ind_id = mt_sample.s,
                        varid = mt_sample.rsid) # varid used to be rsid
    

if __name__=='__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--pop', default=None, type=str, help='population to use')
    parser.add_argument('--overwrite_plink', default=False, action='store_true', help='whether to overwrite existing PLINK files')
    args = parser.parse_args()
    
    n_max = 5000 # maximum number of samples in subset (equal to final sample size if there are sufficient samples for each population)
    not_pop = True
    
    
    if args.pop is None:
        pops = [pop for pop, pop_ct in pop_dict.items()]
    else:
        pops = [args.pop]
    
    print(
    f'''
    pops: {pops}
    not_pop: {not_pop}
    overwrite_plink: {args.overwrite_plink}
    ''')

    mt0 = get_mt(chrom='all', entry_fields = ('GT',)) # default entry_fields will be 'GP', we need 'GT' for exporting to PLINK
 
    ## get mt of population (or every population but that population if `not_pop`=True)
    for pop in pops:
        if not_pop:
            mt_pop = mt0.filter_cols(mt0.pop != pop)
        else:    
            mt_pop = mt0.filter_cols(mt0.pop == pop)
        print(f'\n\nmt sample ct ({pop}, not_pop={not_pop}): {mt_pop.count_cols()}\n\n')
        

    for pop in pops:
        ht_sample_path = f'{ldprune_wd}/subsets/{"not_" if not_pop else ""}{pop}.ht'
        
        if hl.hadoop_exists(f'{ht_sample_path}/_SUCCESS'):
            ht_sample = hl.read_table(ht_sample_path)
            print(f'... Subset ht already exists for pop={pop}, not_pop={not_pop} ...')
            print(f'\n\nSubset ht sample ct: {ht_sample.count()}\n\n')
        else:
            print(f'... getting subset (pop={pop}, not_pop={not_pop}) ...')
            
            ht_sample = get_subset(mt_pop = mt_pop,
                                   pop_dict = pop_dict, 
                                   pop = pop, 
                                   n_max = n_max, 
                                   not_pop = not_pop)
            
            ht_sample_ct = ht_sample.count()
            print(f'\n\nht_sample_ct: {ht_sample_ct}\n\n')
            ht_sample = ht_sample.checkpoint(ht_sample_path)
        
        print(f'... exporting to plink (pop={pop}, not_pop={not_pop}) ...')
        to_plink(pop = pop,
                 mt = mt_pop,
                 ht_sample = ht_sample,
                 not_pop = not_pop,
                 overwrite=args.overwrite_plink)
