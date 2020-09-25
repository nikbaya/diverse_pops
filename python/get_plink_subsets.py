#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 09:17:42 2020

Get PLINK subsets for clumping

@author: nbaya
"""

import hail as hl
import argparse
from ukbb_pan_ancestry.resources.genotypes import get_filtered_mt

hl.init(log='/tmp/hail.log')

bucket = 'gs://ukb-diverse-pops'
REFERENCE_GENOME = 'GRCh37'
CHROMOSOMES = list(map(str, range(1, 23))) + ['X', 'XY']
POPS = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']

MIN_CASES = 50
MIN_CASES_ALL = 100
MIN_CASES_EUR = 100

pop_dict = {
    'AFR': 6637, # dict with sample counts for each population 
    'AMR': 982,
    'CSA': 8876,
    'EAS': 2709,
    'EUR': 420542,
    'MID': 1599
            }


def get_filtered_not_pop_mt(pops: list,
                            chrom: str = 'all',
                            imputed: bool = True,
                            min_mac: int = 20,
                            entry_fields=('GP',),
                            filter_mac_instead_of_ac: bool = False):
    r'''
    Wraps `get_filtered_mt()` from ukbb_pan_ancestry.resources.genotypes
    This filters to 
    '''
    assert len(pops)>0 and pops.issubset(POPS)
    
    kwargs = {
        'pop':('all' if len(pops)>1 else pops[0]),
        'imputed':imputed,
        'min_mac':min_mac,
        'entry_fields':entry_fields,
        'filter_mac_instead_of_ac':filter_mac_instead_of_ac
        }
    
    mt = get_filtered_mt(chrom=chrom, **kwargs)
    
    if len(pops)>1:
        mt = mt.filter_cols(hl.set(pops).contains(mt.pop))
        
    return mt


def get_pop_prop_dict(pop_dict: dict, pops: list) -> (dict, int):
    r'''
    Get population proportions in `pop_dict` for a list of populations `pops`
    '''
    tmp_pop_dict = {pop:n_pop for pop,n_pop in pop_dict.items() if pop in pops}
    n_total = sum(tmp_pop_dict.values())
    pop_prop_dict = {k: v/n_total for k,v in tmp_pop_dict.items()}
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
    if not not_pop and pop!='ALL_POPS' and n_sample == pop_dict[pop]: # if sampling a single population `pop` and n_sample is the same as the population's size. WARNING: Order of conditional statements matters.
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
    
def to_plink(pop: str, 
             subsets_dir, 
             mt, 
             ht_sample, 
             not_pop: bool = False, 
             chr_x_only: bool = False,
             export_varid: bool = True,
             overwrite=False):
    r'''
    Exports matrix table to PLINK2 files
    NOTE: These files will need to split up by chromosome before plink_clump.py
    can be run. 
    '''
    assert 'GT' in mt.entry, "mt must have 'GT' as an entry field"
    assert mt.GT.dtype==hl.tcall, "entry field 'GT' must be of type `Call`"
    print(chr_x_only)
    bfile_path = f'{subsets_dir}/{"not_" if not_pop else ""}{pop}/{"not_" if not_pop else ""}{pop}'
    if chr_x_only: bfile_path += '.chrX'
    
    if not overwrite and all([hl.hadoop_exists(f'{bfile_path}.{suffix}') for suffix in ['bed','bim','fam']]):
        print(f'\nPLINK files already exist for {"not_" if not_pop else ""}{pop}')
        print(bfile_path)
    else:
        print(f'Saving to bfile prefix {bfile_path}')
        mt_sample = mt.annotate_rows(varid = hl.str(mt.locus)+':'+mt.alleles[0]+':'+mt.alleles[1])
        mt_sample = mt_sample.filter_cols(hl.is_defined(ht_sample[mt_sample.s]))
        hl.export_plink(dataset = mt_sample, 
                        output = bfile_path, 
                        ind_id = mt_sample.s,
                        varid = mt_sample.varid) # varid used to be rsid
def export_varid(args):
    n_max = 5000
    subsets_dir = f'{bucket}/ld_prune/subsets_{round(n_max/1e3)}k' 
    
    mt = get_filtered_not_pop_mt(chrom='all' if not args.chr_x_only else 'X', 
                                 entry_fields=('GT',)) # default entry_fields will be 'GP', we need 'GT' for exporting to PLINK
    
    mt_sample = mt.annotate_rows(chrom = mt.locus.contig,
                                 pos = mt.locus.position,
                                 varid = hl.str(mt.locus)+':'+mt.alleles[0]+':'+mt.alleles[1])
    
    rows = mt_sample.rows()
    rows = rows.key_by()
    rows = rows.select('chrom','pos','varid')
    rows.export(f'{subsets_dir}/varid.txt',delimiter=' ')
        
def main(args):
    n_max = 5000 # maximum number of samples in subset (equal to final sample size if there are sufficient samples for each population)
    not_pop = args.not_pop
    
    subsets_dir = f'{bucket}/ld_prune/subsets_{round(n_max/1e3)}k' 
    
    pop = args.pop.upper()
    assert pop in POPS+['ALL_POPS'], f'Invalid population passed pop flag: {pop}'
    
    print(
    f'''
    pop: {pop}
    not_pop: {not_pop}
    chr_x_only: {args.chr_x_only}
    overwrite_plink: {args.overwrite_plink}
    ''')

    mt_pop = get_filtered_not_pop_mt(chrom='all' if not args.chr_x_only else 'X',  # chrom='all' includes autosomes and chrX
                                     pop='all' if args.pop=='ALL_POPS' else pop,
                                     entry_fields=('GT',),
                                     not_pop=not_pop) # default entry_fields will be 'GP', we need 'GT' for exporting to PLINK

    print(f'\n\nmt sample ct ({pop}, not_pop={not_pop}): {mt_pop.count_cols()}\n\n')
    print(f'\n\nmt variant ct ({pop}, not_pop={not_pop}): {mt_pop.count_rows()}\n\n')
    
    ht_sample_path = f'{subsets_dir}/{"not_" if not_pop else ""}{pop}.ht'
    print(ht_sample_path)
    if hl.hadoop_exists(f'{ht_sample_path}/_SUCCESS'):
        ht_sample = hl.read_table(ht_sample_path)
        print(f'... Subset ht already exists for pop={pop}, not_pop={not_pop} ...')
        print(f'\n\nSubset ht sample ct: {ht_sample.count()}\n\n')
    else:
        pass
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
             subsets_dir=subsets_dir,
             mt = mt_pop,
             ht_sample = ht_sample,
             not_pop = not_pop,
             chr_x_only=args.chr_x_only,
             overwrite=args.overwrite_plink)
    

if __name__=='__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--pop', default=None, type=str, help='population to use')
    parser.add_argument('--not_pop', default=False, action='store_true', help='if true, get leave-one-out pop set without `pop`')
    parser.add_argument('--chr_x_only', default=False, action='store_true', help='Only export chr X')
    parser.add_argument('--overwrite_plink', default=False, action='store_true', help='whether to overwrite existing PLINK files')
    parser.add_argument('--export_varid', default=False, action='store_true', help='export varids')
    args = parser.parse_args()
    
    if args.export_varid:
        export_varid(args=args)
    else:
        main(args=args)