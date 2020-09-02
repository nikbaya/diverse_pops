#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 17:52:37 2020

@author: nbaya
"""

from hail.linalg import BlockMatrix
import hail as hl
from ukbb_pan_ancestry import *
from collections import namedtuple

P_THRESHOLDS = {'s1': 5e-8, 's2': 1e-6, 's3': 1e-4, 's4': 1e-3, 's5': 1e-2, 's6': .05, 's7': .1, 's8': .2, 's9': .5, 's10': 1.}

scratch_dir = 'gs://ukbb-diverse-temp-30day/nb-scratch'

# NOTE: These methods are hardcoded because they had some issues being loaded from ukbb_pan_ancestry 

def get_clumping_results_path(pop: str = 'full', high_quality_only: bool = False, 
                              not_pop: bool = True):
    mt_name = f'{"not_" if not_pop else ""}{pop}.mt' if pop != 'full' else 'full_clump_results.mt'
    return f'{bucket}/ld_prune/clump_results{"_high_quality" if high_quality_only else ""}/{mt_name}'


def get_prs_mt_path(high_quality: bool = True):
    return f'{bucket}/prs/all_combos_prs{"" if high_quality else "_raw"}.mt'


def get_clump_sumstats_bm_path(high_quality: bool = True):
    return f'{temp_bucket}/prs/clumped_sumstats{"" if high_quality else "_raw"}.bm'


def get_clump_sumstats_col_ht_path(high_quality: bool = True):
    return f'{temp_bucket}/prs/clumped_sumstats{"" if high_quality else "_raw"}.cols.ht'

def get_prs_bm_path(high_quality: bool = True):
    return f'{temp_bucket}/prs/prs{"" if high_quality else "_raw"}.bm'


def get_prs_assess_ht_path(high_quality: bool = True):
    return f'{temp_bucket}/prs/assess_prs{"" if high_quality else "_raw"}.ht'

def get_filtered_mt_with_x(pop: str = 'all'):
    mt = get_filtered_mt(pop=pop, entry_fields=('dosage',))
    mt_x = get_filtered_mt('X', pop=pop, entry_fields=('dosage',))
    # assert mt.s.collect() == mt_x.s.collect()
    mt = mt.union_rows(mt_x)
    return mt

def separate_results_mt_by_pop(mt, col_field = 'pheno_data', entry_field = 'summary_stats', skip_drop: bool = False):
    mt = mt.annotate_cols(col_array=hl.zip_with_index(mt[col_field])).explode_cols('col_array')
    mt = mt.transmute_cols(pop_index=mt.col_array[0], **{col_field: mt.col_array[1]})
    mt = mt.annotate_entries(**{entry_field: mt[entry_field][mt.pop_index]})
    if not skip_drop:
        mt = mt.drop('pop_index')
    return mt

def explode_by_p_threshold(mt):
    mt = mt.annotate_cols(p_threshold=hl.literal(list(P_THRESHOLDS.items()))).explode_cols('p_threshold')
    mt = mt.transmute_cols(p_threshold_name=mt.p_threshold[0], p_threshold=mt.p_threshold[1])
    return mt

def get_test_genotypes_bm(chrom, genotype_bm_path):
    
    meta_mt = hl.read_matrix_table(get_meta_analysis_results_path())
#   

    if chrom=='all':
        mt = get_filtered_mt_with_x()
    else:
        mt = get_filtered_mt(chrom=chrom, entry_fields=('dosage',))
        
    mt = mt.filter_rows(hl.is_defined(meta_mt.rows()[mt.row_key]))
    
#    if not hl.hadoop_is_file(f'{genotype_samples_ht_path}/_SUCCESS'):
#        samples = mt.s.take(10)
#        mt = mt.filter_cols(hl.literal(samples).contains(mt.s))
#        mt = mt.key_cols_by(userId=hl.int32(mt.s))
#        mt.cols().add_index().write(genotype_samples_ht_path, overwrite=True)
#    else:
#        samples_ht = hl.read_table(genotype_samples_ht_path)
    controls = hl.read_table(f'{scratch_dir}/genotype_samples_n10.ht')
    cases = hl.read_table(f'{scratch_dir}/genotype_samples_n10_cases.ht')
    samples_ht = cases.union(controls)
    mt = mt.filter_cols(hl.is_defined(samples_ht[hl.int32(mt.s)]))

    mt = mt.key_cols_by(userId=hl.int32(mt.s))
    print(mt.count())
    
    mt = mt.select_cols().select_rows()
    mt = mt.repartition(1000)
    BlockMatrix.write_from_entry_expr(mt.dosage, genotype_bm_path, overwrite=True)
    
def make_sumstats_bm(sumstats_bm_path, high_quality):
    meta_mt = hl.read_matrix_table(get_meta_analysis_results_path())
    clump_mt = hl.read_matrix_table(get_clumping_results_path(high_quality_only=high_quality)).rename({'pop': 'clump_pops'})
    mt = all_axis_join(meta_mt, clump_mt)
    mt = separate_results_mt_by_pop(mt, 'clump_pops', 'plink_clump', skip_drop=True)
    mt = separate_results_mt_by_pop(mt, 'meta_analysis_data', 'meta_analysis', skip_drop=True)
    mt = mt.filter_cols(mt.meta_analysis_data.pop == mt.clump_pops)
    mt = explode_by_p_threshold(mt).unfilter_entries()
    
    mt = mt.filter_cols((mt.description=='Type 2 diabetes')&(mt.p_threshold==1))
    
    BlockMatrix.write_from_entry_expr(
                hl.or_else(mt.meta_analysis.BETA * hl.is_defined(mt.plink_clump.TOTAL) * hl.int(mt.meta_analysis.Pvalue < mt.p_threshold), 0.0),
                sumstats_bm_path, overwrite=True)

def get_test_genotypes_mt(chrom, genotype_samples_ht_path, genotype_mt_path, cases_only):
    meta_mt = hl.read_matrix_table(get_meta_analysis_results_path())
#   
    if chrom=='all':
        mt = get_filtered_mt_with_x()
    else:
        mt = get_filtered_mt(chrom=chrom, entry_fields=('dosage',))
        
    if status=='cases':
        t2d_ht = hl.read_table(f'gs://ukbb-diverse-temp-30day/nb-scratch/t2d.ht/')
        t2d_ht = t2d_ht.filter(t2d_ht.both_sexes==1)
        t2d_ht = t2d_ht.key_by('userId')
        mt = mt.filter_cols(hl.is_defined(t2d_ht[hl.int32(mt.s)]))
        
    mt = mt.filter_rows(hl.is_defined(meta_mt.rows()[mt.row_key]))
    
    if not hl.hadoop_is_file(f'{genotype_samples_ht_path}/_SUCCESS'):
        samples = mt.s.take(10)
        mt = mt.filter_cols(hl.literal(samples).contains(mt.s))
        mt = mt.key_cols_by(userId=hl.int32(mt.s))
        mt.cols().add_index().write(genotype_samples_ht_path, overwrite=True)
    else:
        samples_ht = hl.read_table(genotype_samples_ht_path)
        samples = samples_ht.s.collect()
        mt = mt.filter_cols(hl.literal(samples).contains(mt.s))
        mt = mt.key_cols_by(userId=hl.int32(mt.s))
    
    mt = mt.select_entries('dosage')
    mt = mt.select_rows()
    mt = mt.select_cols()
    mt = mt.repartition(10)
    mt.write(genotype_mt_path)
    
def compute_test_prs_bm(genotype_bm_path, prs_bm_path, args):
    sumstats_bm = BlockMatrix.read(get_clump_sumstats_bm_path(args.high_quality))
    genotype_bm = BlockMatrix.read(genotype_bm_path)
    prs_bm: BlockMatrix = genotype_bm.T @ sumstats_bm
    prs_bm.write(prs_bm_path, args.overwrite)
    
def repartition(genotype_mt_path1, genotype_mt_path2):
    mt = hl.read_matrix_table(genotype_mt_path1)
    mt = mt.repartition(10)
    mt = mt.write(genotype_mt_path2)
    


def compute_prs_mt(genotype_mt_path, prs_mt_path):
    scratch_dir = 'gs://ukbb-diverse-temp-30day/nb-scratch'

    clumped = hl.read_table('gs://ukb-diverse-pops/ld_prune/results_high_quality/not_AMR/phecode-250.2-both_sexes/clump_results.ht/')
    sumstats = hl.import_table('gs://ukb-diverse-pops/sumstats_flat_files/phecode-250.2-both_sexes.tsv.bgz',
                               impute=True)
    sumstats = sumstats.annotate(locus = hl.locus(sumstats.chr, sumstats.pos),
                                 alleles = hl.array([sumstats.ref, sumstats.alt]))
    sumstats = sumstats.key_by('locus','alleles')
    sumstats.describe()
#    mt = hl.read_matrix_table(genotype_mt_path) # read genotype mt subset
    
    # get full genotype mt
    meta_mt = hl.read_matrix_table(get_meta_analysis_results_path())
    mt = get_filtered_mt_with_x()        
    mt = mt.filter_rows(hl.is_defined(meta_mt.rows()[mt.row_key]))
    mt = mt.select_entries('dosage')
    mt = mt.select_rows()
    mt = mt.select_cols()
    
    mt = mt.annotate_rows(beta = hl.if_else(hl.is_defined(clumped[mt.row_key]), sumstats[mt.row_key].beta_meta, 0))
    mt = mt.annotate_cols(score = hl.agg.sum(mt.beta*mt.dosage))
    mt_cols = mt.cols()
    mt_cols = mt_cols.repartition(1000)
    mt_cols.write(f'{scratch_dir}/prs_all_samples.ht')
#    mt.write(prs_mt_path)

def main():
     # args
    high_quality=True
    overwrite=False
    struct = namedtuple("struct", "high_quality overwrite")
    args = struct(high_quality=high_quality, overwrite=overwrite)
    
    hl.init(default_reference='GRCh37', log='/prs.log',
            spark_conf={'spark.hadoop.fs.gs.requester.pays.mode': 'AUTO', 
                        'spark.hadoop.fs.gs.requester.pays.project.id': 'ukbb-diversepops-neale'})
    
#    cases_only = True
    
    chrom='all'
#    genotype_samples_ht_path = f'{scratch_dir}/genotype_samples_n10{"_cases" if cases_only else ""}.ht'
    genotype_bm_path = f'{scratch_dir}/genotype_bm_n10.bm'

    get_test_genotypes_bm(chrom=chrom,
                          genotype_bm_path=genotype_bm_path)
   
#    genotype_mt_path1 = f'{scratch_dir}/genotype_bm_n10{f"_chr{chrom}" if chrom!="all" else ""}{"_cases" if cases_only else ""}.mt'
#    genotype_mt_path2 = f'{scratch_dir}/genotype_mt_n10{f"_chr{chrom}" if chrom!="all" else ""}{"_cases" if cases_only else ""}.mt'

#    get_test_genotypes_mt(chrom=chrom,
#                          genotype_samples_ht_path=genotype_samples_ht_path, 
#                          genotype_mt_path=genotype_mt_path1,
#                          cases_only=cases_only)

#    repartition(genotype_mt_path1, genotype_mt_path2)
    
#    sumstats_bm_path = f'{scratch_dir}/sumstats_t2d.bm'
#    make_sumstats_bm(sumstats_bm_path, high_quality)
    
#    prs_bm_path = f'{scratch_dir}/prs_n10.bm'
#    prs_mt_path = f'{scratch_dir}/prs_n10{"_cases" if cases_only else ""}.mt'
#    
#    compute_test_prs_bm(genotype_bm_path, prs_mt_path, args)


##    prs_mt_path = f'{scratch_dir}/prs_n10{"_cases" if cases_only else ""}.mt'
#    genotype_mt_path = f'{scratch_dir}/genotype_all_samples.mt'
#    prs_mt_path = f'{scratch_dir}/prs_all_samples.mt'
#    compute_prs_mt(genotype_mt_path, prs_mt_path)
    
    

if __name__=='__main__':
    main()
   
    
    