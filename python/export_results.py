#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 21:32:56 2020

Export flat file summary statistics from matrix tables

@author: nbaya
"""

import argparse
import hail as hl
from itertools import combinations
from time import time

bucket = 'gs://ukb-diverse-pops'
ldprune_dir = f'{bucket}/ld_prune'

def get_variant_results_path(pop: str):
    return f'{bucket}/combined_results/results_{pop}.mt'


def export_results(num_pops, batch_size=256):
    r'''
    `num_pops`: exact number of populations for which phenotype is defined
    '''
    hl.init(default_reference='GRCh38', log='/tmp/export_entries_by_col.log')
    mt0 = hl.read_matrix_table(get_variant_results_path(pop='full'))
    meta_mt0 = hl.read_matrix_table('gs://ukb-diverse-pops/combined_results/meta_analysis.mt')
    
    mt0 = mt0.annotate_cols(pheno_id = (mt0.trait_type+'-'+
                            mt0.phenocode+'-'+
                            mt0.pheno_sex+
                            hl.if_else(hl.len(mt0.coding)>0, '-'+mt0.coding, '')+
                            hl.if_else(hl.len(mt0.modifier)>0, '-'+mt0.modifier, '')
                            ).replace(' ','_').replace('/','_'))
    mt0 = mt0.annotate_rows(chr = mt0.locus.contig,
                            pos = mt0.locus.position,
                            ref = mt0.alleles[0],
                            alt = mt0.alleles[1])
    
    all_pops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']
    
    trait_types_to_run = ['continuous','biomarkers','categorical','phecode', 'icd10', 'prescriptions'] # list of which trait_type to run
    pop_sets = [set(i) for i in list(combinations(all_pops, num_pops))] # list of exact set of pops for which phenotype is defined
    
    
    # fields specific to each category of trait
    quant_meta_fields = ['AF_Allele2']
    quant_fields = ['AF_Allele2']
    
    binary_meta_fields = ['AF_Cases','AF_Controls']
    binary_fields = ['AF.Cases','AF.Controls']
    
    # dictionaries for renaming fields
    quant_meta_field_rename_dict = {'AF_Allele2':'af_meta',
                              'BETA':'beta_meta',
                              'SE':'se_meta',
                              'Pvalue':'pval_meta',
                              'Pvalue_het':'pval_heterogeneity'}
    quant_field_rename_dict = {'AF_Allele2':'af',
                         'BETA':'beta',
                         'SE':'se',
                         'Pvalue':'pval',
                         'low_confidence':'low_confidence'} # decided on this implementation to make later code cleaner
    
    binary_meta_field_rename_dict = {'BETA':'beta_meta',
                                     'SE':'se_meta',
                                     'Pvalue':'pval_meta',
                                     'AF_Cases':'af_cases_meta',
                                     'AF_Controls':'af_controls_meta',
                                     'Pvalue_het':'pval_heterogeneity'}
    binary_field_rename_dict = {'AF.Cases':'af_cases',
                                'AF.Controls':'af_controls',
                                'BETA':'beta',
                                'SE':'se',
                                'Pvalue':'pval',
                                'low_confidence':'low_confidence'} # decided on this implementation to make later code cleaner
    
    all_quant_trait_types = {'continuous','biomarkers'}
    all_binary_trait_types = {'categorical','phecode', 'icd10', 'prescriptions'}
    
    quant_trait_types = all_quant_trait_types.intersection(trait_types_to_run) # get list of quant trait types to run
    binary_trait_types = all_binary_trait_types.intersection(trait_types_to_run) # get list of binary trait types to run
    error_trait_types = set(trait_types_to_run).difference(quant_trait_types.union(binary_trait_types))
    assert len(error_trait_types)==0, f'ERROR: The following trait_types are invalid: {error_trait_types}'
        
    for trait_category, trait_types in [('quant', quant_trait_types), ('binary', binary_trait_types)]:
        if len(trait_types)==0: #if no traits in trait_types list
            continue
        print(f'{trait_category} trait types to run: {trait_types}')
        
        if trait_category == 'quant':
            meta_fields = quant_meta_fields
            fields = quant_fields
            meta_field_rename_dict = quant_meta_field_rename_dict
            field_rename_dict = quant_field_rename_dict
        elif trait_category == 'binary':
            meta_fields = binary_meta_fields
            fields = binary_fields
            meta_field_rename_dict = binary_meta_field_rename_dict
            field_rename_dict = binary_field_rename_dict
    
        meta_fields += ['BETA','SE','Pvalue','Pvalue_het']
        fields += ['BETA','SE','Pvalue','low_confidence']
            
        for pop_set in pop_sets:    
            start = time()
            
            mt1 = mt0.filter_cols((hl.literal(trait_types).contains(mt0.trait_type))&
                                  (hl.set(mt0.pheno_data.pop)==hl.literal(pop_set)))
            
            col_ct = mt1.count_cols()
            if col_ct==0:
                print(f'Skipping {trait_types},{sorted(pop_set)}, no phenotypes found')
                continue
            
            pop_list = sorted(pop_set)
            
            pop_idx_list = [pop_list.index(pop) for pop in all_pops if pop in pop_list]
            
            annotate_dict = {}
            
            keyed_mt = meta_mt0[mt1.row_key,mt1.col_key]
            for field in meta_fields: # NOTE: Meta-analysis columns go before per-population columns
        #        annotate_dict.update({f'{meta_field_rename_dict[field]}': hl.float64(hl.format('%.3e', keyed_mt.meta_analysis[field][0]))})
                field_expr = keyed_mt.meta_analysis[field][0]
                annotate_dict.update({f'{meta_field_rename_dict[field]}': hl.if_else(hl.is_nan(field_expr),
                                                                          hl.str(field_expr),
                                                                          hl.format('%.3e', field_expr))})
        
            for field in fields:
                for pop_idx in pop_idx_list:
        #            annotate_dict.update({f'{field_rename_dict[field]}_{pops[pop_idx]}': hl.format('%.3e', mt1.summary_stats[field][pop_idx])})
                    field_expr = mt1.summary_stats[field][pop_idx]
                    annotate_dict.update({f'{field_rename_dict[field]}_{all_pops[pop_idx]}': hl.if_else(hl.is_nan(field_expr),
                                                                                             hl.str(field_expr),
                                                                                             hl.str(field_expr) if field=='low_confidence' else hl.format('%.3e', field_expr))})
            
            mt2 = mt1.annotate_entries(**annotate_dict)
            
            mt2 = mt2.filter_cols(mt2.coding != 'zekavat_20200409')
            mt2 = mt2.key_cols_by('pheno_id')
            mt2 = mt2.key_rows_by().drop('locus','alleles','gene','annotation','summary_stats')
            
            batch_idx = 1
            export_out = f'{ldprune_dir}/release/{trait_category}/{"-".join(pop_list)}_batch{batch_idx}'
            while hl.hadoop_is_dir(export_out):
                batch_idx += 1
                export_out = f'{ldprune_dir}/release/{trait_category}/{"-".join(pop_list)}_batch{batch_idx}'
            print(f'\nExporting {col_ct} phenos to: {export_out}\n')
            hl.experimental.export_entries_by_col(mt = mt2,
                                                  path = export_out,
                                                  bgzip = True,
                                                  batch_size = batch_size,
                                                  use_string_key_as_file_name = True,
                                                  header_json_in_file = False)
            end = time()
            print(f'\nExport complete for:\n{trait_types}\n{pop_list}\ntime: {round((end-start)/3600,2)} hrs')

def make_pheno_manifest():
    # mt0 = hl.read_matrix_table('gs://ukb-diverse-pops/combined_results/results_full.mt')
    # mt0_cols = mt0.cols()
    ht0 = hl.read_table('gs://ukb-diverse-pops-public/sumstats_qc_analysis/lambda/lambdas_full.ht')
    
    annotate_dict = {}
    for field in ['n_cases','n_controls','heritability','lambda_gc']:
        for pop in ['AFR','AMR','CSA','EAS','EUR','MID']:
            new_field = field if field!='heritability' else 'saige_heritability' # new field name (only applicable to saige heritability)
            idx = ht0.pheno_data.pop.index(pop)
            field_expr = ht0.pheno_data[field]
            annotate_dict.update({f'{new_field}_{pop}': hl.if_else(hl.is_nan(idx),
                                                               hl.null(field_expr[0].dtype),
                                                               field_expr[idx])})
    annotate_dict.update({'filename':(ht0.trait_type+'-'+
                                     ht0.phenocode+'-'+
                                     ht0.pheno_sex+
                                     hl.if_else(hl.len(ht0.coding)>0, '-'+ht0.coding, '')+
                                     hl.if_else(hl.len(ht0.modifier)>0, '-'+ht0.modifier, '')+
                                     '.tsv.bgz'
                                     ).replace(' ','_').replace('/','_')})
    ht1 = ht0.annotate(**annotate_dict)
    ht1 = ht1.drop('pheno_data')
    ht1.export(f'{ldprune_dir}/phenotype_manifest.tsv')
    
def make_variant_manifest():
    # mt0 = hl.read_matrix_table('gs://ukb-diverse-pops/combined_results/results_full.mt')
    # TODO: annotate on the data from get_gene_intervals_path()
    pass


if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--num_pops', type=int, help='number of defined pops (options: 6,5,...,2,1)')
    args = parser.parse_args()

    export_results(num_pops=args.num_pops, 
                   batch_size=256)