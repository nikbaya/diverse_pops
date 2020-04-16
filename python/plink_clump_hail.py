#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 15:09:04 2020

Hail script for clumping GWAS results with PLINK

@author: nbaya
"""

import argparse
import hail as hl
import sys

bucket = 'gs://ukb-diverse-pops'
ldprune_dir = f'{bucket}/ld_prune'


PILOT_PHENOTYPES = set(map(lambda x: (x, 'irnt', 'continuous'), {'50', '699', '23104'})).union(
    set(map(lambda x: (*x, 'categorical'), {('20004', '1095'), ('20004', '1479')}))).union(
    set(map(lambda x: (x, '', 'icd'), {'K519', 'K509', 'E109', 'E119', 'J459', 'I251',
                                       'K51', 'K50', 'E10', 'E11', 'J45', 'I25'}))).union(
    set(map(lambda x: (x, 'icd10', 'icd_all'), {'K519', 'K509', 'E109', 'E119', 'J459', 'I251',
                                       'K51', 'K50', 'E10', 'E11', 'J45', 'I25'}))).union(
    set(map(lambda x: (x, 'both_sexes', 'phecode'), {'401', '411'}))).union(
    {('1717', '1717', 'continuous'),
     ('random', 'random', 'continuous'),
     ('random', 'random_strat', 'continuous'),
     ('whr', 'whr', 'continuous'),
     ('1747', '4', 'categorical'),
     ('30040', 'irnt', 'continuous'),
     ('30890', '30890', 'biomarkers'),
     ('HMG CoA reductase inhibitor|statin', '', 'prescriptions')})

def ht_to_tsv(args):
    r'''
    Convert Hail table of variant results to a tsv 
    '''
    add_args = {}
    if args.n_threads is not None:
        add_args['master'] = f'local[{args.n_threads}]'
    hl.init(default_reference='GRCh38', log='/tmp/ht_to_tsv.log', **add_args)
    ht = hl.read_table(args.input_file)
    print(ht.describe())
    ht.export(args.output_file)

    
def tsv_to_ht(args):
    r'''
    Convert tsv to a Hail table
    '''
    print(sys.version)
    add_args = {}
    if args.n_threads is not None:
        add_args['master'] = f'local[{args.n_threads}]'
    hl.init(default_reference='GRCh38', log='/tmp/tsv_to_ht.log', **add_args)
    ht = hl.import_table(args.input_file, no_header=True)
    ht = ht.rename({'f0':'contig',
                    'f1':'F',
                    'f2':'varid',
                    'f3':'pos',
                    'f4':'P',
                    'f5':'TOTAL',
                    'f6':'NSIG',
                    'f7':'S05',
                    'f8':'S01',
                    'f9':'S001',
                    'f10':'S0001',
                    'f11':'SP2'})
    print(ht.describe())
    ht.write(args.output_file, overwrite=args.overwrite)

def write_meta_sumstats(args):
    r'''
    Extracting entries from full meta-analysis mt. Run separately, not in pipeline
    '''
    add_args = {}
    if args.n_threads is not None:
        add_args['master'] = f'local[{args.n_threads}]'
    hl.init(default_reference='GRCh38', log='/tmp/write_meta_sumstats.log', **add_args)
    mt0 = hl.read_matrix_table('gs://ukb-diverse-pops/combined_results/meta_analysis.mt')
    
    mt = mt0.filter_cols((hl.len(mt0.meta_analysis_data[0].pop)==6)|
            (hl.literal(PILOT_PHENOTYPES).contains((mt0.pheno, mt0.coding, mt0.trait_type))))
  # get phenos with results for all 6 populations
    mt2 = mt.select_cols().select_rows()
    mt2 = mt2.select_entries(Pvalue = mt2.meta_analysis.Pvalue)
#    hl.experimental.export_entries_by_col(mt=mt2, path=, batch_size: int = 256, bgzip: bool = True, header_json_in_file: bool = True)
    ht1 = mt2.entries()
    ht1 = ht1.key_by('pheno','trait_type','coding','locus','alleles')
    ht1.describe()
    ht1 = ht1.annotate(varid = hl.str(ht1.locus)+'_'+ht1.alleles[0]+'_'+ht1.alleles[1])
    out = f'{ldprune_dir}/meta_analysis.all_pops.ht'
    print(f'\n...Writing meta-analysis sumstats for all pops...\nFrom: {args.input_file}\nTo:{out}')
    ht1.write(out, overwrite=True)

def get_meta_sumstats(args):
    print(ldprune_dir)
    print(f'\n...Getting meta-analysis sumstats...\nUsing: {args.input_file}\n')
    pops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']
    add_args = {}
    if args.n_threads is not None:
        add_args['master'] = f'local[{args.n_threads}]'
    hl.init(default_reference='GRCh38', log='/tmp/mt_to_tsv.log', **add_args)
    ht = hl.read_table(args.input_file)
    ht_pheno = ht.filter((ht.pheno==args.pheno)&
                         (ht.coding==args.coding)&
                         (ht.trait_type==args.trait_type))
    loo_meta_idx = pops.index(args.pop)+1 # get index of leave-one-out meta-analysis result for population `pop`
    ht_pheno = ht_pheno.key_by('varid') # key by varid, then select P so that those are the only fields when exporting to tsv
    ht_pheno = ht_pheno.select(P=ht_pheno.Pvalue[loo_meta_idx])
    ht_pheno = ht_pheno.filter(hl.is_defined(ht_pheno.P))
    ht_pheno.export(args.output_file)
    
    
def test_get_meta_sumstats(args):
    print(ldprune_dir)
    print(f'\n...Getting meta-analysis sumstats...\nUsing: {args.input_file}\n')
#    pops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']
    add_args = {}
    if args.n_threads is not None:
        add_args['master'] = f'local[{args.n_threads}]'
    hl.init(default_reference='GRCh38', log='/tmp/mt_to_tsv.log', **add_args)
    ht = hl.read_table(args.input_file)
    ht_pheno = ht.filter((ht.pheno==args.pheno)&
                         (ht.coding==args.coding)&
                         (ht.trait_type==args.trait_type))
    ht_pheno = ht_pheno.key_by('varid') # key by varid, then select P so that those are the only fields when exporting to tsv
#    loo_meta_idx = pops.index(args.pop)+1 # get index of leave-one-out meta-analysis result for population `pop`
#    ht_pheno = ht_pheno.select(P=ht_pheno.Pvalue[loo_meta_idx])
#    ht_pheno = ht_pheno.filter(hl.is_defined(ht_pheno.P))
    ht_pheno = ht_pheno.annotate(b = hl.rand_norm(0,0.1), # SNP effect
                                 se = hl.rand_norm(0,0.1)**2, # standard error of SNP effect
                                 A1 = ht_pheno.varid.split('_')[-1], # effect alele (A2 in UKB)
                                 A2 = ht_pheno.varid.split('_')[-2], # non-effect alele (A1 in UKB)
                                 af = hl.rand_unif(0.05, 0.95),
                                 N = 360000) # frequency of A1 allele
    ht_pheno = ht_pheno.annotate(P = 2*hl.pnorm(-hl.abs(ht_pheno.b/ht_pheno.se)))
    ht_pheno.export(args.output_file)
    
    

def export_cont_results():
    hl.init(default_reference='GRCh38', log='/tmp/export_entries_by_col.log')
    mt0 = hl.read_matrix_table('gs://ukb-diverse-pops/combined_results/results_full.mt')
    meta_mt0 = hl.read_matrix_table('gs://ukb-diverse-pops/combined_results/meta_analysis.mt')
    
    mt0 = mt0.annotate_cols(pheno_id = mt0.trait_type+'-'+mt0.pheno+'-'+mt0.coding)
    mt0 = mt0.annotate_rows(chr = mt0.locus.contig,
                            pos = mt0.locus.position,
                            ref = mt0.alleles[0],
                            alt = mt0.alleles[1])
    
#    mt1 = mt0.filter_cols((hl.literal(['biomarkers']).contains(mt0.trait_type))&(hl.len(mt0.pheno_data.pop)==6))    
    mt1 = mt0.filter_cols((hl.literal(['continuous']).contains(mt0.trait_type))&(hl.len(mt0.pheno_data.pop)==6))    
##    mt1 = mt0.filter_cols(((hl.literal(['biomarkers']).contains(mt0.trait_type))&(hl.len(mt0.pheno_data.pop)==5))) 
#    mt1 = mt0.filter_cols((hl.literal(['biomarkers','continuous']).contains(mt0.trait_type))&(
#                            hl.len(mt0.pheno_data.pop)==6)&(
#                            hl.literal(PILOT_PHENOTYPES).contains((mt0.pheno, mt0.coding, mt0.trait_type))))  
    
#    mt1 = mt0.filter_cols(((hl.literal(['continuous']).contains(mt0.trait_type))&(hl.len(mt0.pheno_data.pop)==3)))
#    pheno_id_list = mt1.pheno_id.take(5)
#    mt1 = mt1.filter_cols((hl.literal(pheno_id_list).contains(mt1.pheno_id)))
    
    pops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']
    pop_list = mt1.pheno_data.pop.take(1)[0]
    
    pop_idx_list = [pop_list.index(pop) for pop in pops if pop in pop_list]
    print(f'pop_idx_list: {pop_idx_list}')
    
    
    field_rename_dict = {'AF_Allele2':'AF',
                         'BETA':'beta',
                         'SE':'se',
                         'Pvalue':'pval',
                        }
    meta_field_rename_dict = {'BETA':'beta_meta',
                              'SE':'se_meta',
                              'Pvalue':'pval_meta',
                              'Q':'Q',
                              'Pvalue_het':'pval_het',
                              'AF_Allele2':'AF_meta'}
    
    annotate_dict = {}
    
    for field in ['AF_Allele2','BETA','SE','Pvalue']:
        for pop_idx in pop_idx_list:
            annotate_dict.update({f'{field_rename_dict[field]}_{pops[pop_idx]}': hl.if_else(hl.is_nan(mt1.summary_stats[field][pop_idx]),
                                                                                        hl.str(mt1.summary_stats[field][pop_idx]),
                                                                                        hl.format('%.3e', mt1.summary_stats[field][pop_idx]))})
#        for pop_idx in range(len(pops)):
##                annotate_dict.update({f'{field_rename_dict[field]}_{pops[pop_idx]}': mt1.summary_stats[field][pop_idx] if pop_idx in pop_idx_list else hl.null(hl.tfloat64)})
##            annotate_dict.update({f'{field_rename_dict[field]}_{pops[pop_idx]}': hl.float64(hl.format('%.3e', mt1.summary_stats[field][pop_idx]) if pop_idx in pop_idx_list else hl.null(hl.tstr))})
#                annotate_dict.update({f'{field_rename_dict[field]}_{pops[pop_idx]}': hl.if_else(hl.is_nan(mt1.summary_stats[field][pop_idx]),
#                                                                                        hl.str(mt1.summary_stats[field][pop_idx]),
#                                                                                        hl.format('%.3e', mt1.summary_stats[field][pop_idx])) if pop_idx in pop_idx_list else hl.null(hl.tfloat64)})

    
    keyed_mt = meta_mt0[(mt1.locus,mt1.alleles),(mt1.pheno, mt1.coding, mt1.trait_type)]
    for field in ['AF_Allele2','BETA','SE','Pvalue','Pvalue_het']:
        # TODO: Make this generalizable to any combination of populations for meta-analysis
#            annotate_dict.update({f'{meta_field_rename_dict[field]}': keyed_mt.meta_analysis[field][0]})
#            annotate_dict.update({f'{meta_field_rename_dict[field]}': hl.float64(hl.format('%.3e', keyed_mt.meta_analysis[field][0]))})
            annotate_dict.update({f'{meta_field_rename_dict[field]}': hl.if_else(hl.is_nan(keyed_mt.meta_analysis[field][0]),
                                                                                hl.str(keyed_mt.meta_analysis[field][0]),
                                                                                hl.format('%.3e', keyed_mt.meta_analysis[field][0]))})
    
    mt2 = mt1.annotate_entries(**annotate_dict)
    
    mt2 = mt2.key_cols_by('pheno_id')
    mt2 = mt2.key_rows_by().drop('locus','alleles','gene','annotation','summary_stats')
    
    batch_idx = 1
    export_out = f'{ldprune_dir}/release/batch{batch_idx}'
    while hl.hadoop_is_dir(export_out):
        batch_idx += 1
        export_out = f'{ldprune_dir}/release/batch{batch_idx}'
    print(f'\nSaving to: {export_out}\n')
    hl.experimental.export_entries_by_col(mt = mt2,
                                          path = export_out,
                                          bgzip = True,
                                          batch_size = 256,
                                          use_string_key_as_file_name = True,
                                          header_json_in_file = False)

def export_binary_results():
    hl.init(default_reference='GRCh38', log='/tmp/export_entries_by_col.log')
    mt0 = hl.read_matrix_table('gs://ukb-diverse-pops/combined_results/results_full.mt')
    meta_mt0 = hl.read_matrix_table('gs://ukb-diverse-pops/combined_results/meta_analysis.mt')
    
    
    mt0 = mt0.annotate_cols(pheno_id = mt0.trait_type+'-'+mt0.pheno+'-'+mt0.coding)
    mt0 = mt0.annotate_rows(chr = mt0.locus.contig,
                            pos = mt0.locus.position,
                            ref = mt0.alleles[0],
                            alt = mt0.alleles[1])
    
    mt1 = mt0.filter_cols((hl.literal(['icd_all','prescriptions']).contains(mt0.trait_type))&(hl.len(mt0.pheno_data.pop)==6))
    
    pops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']
    pop_list = mt1.pheno_data.pop.take(1)[0]
    
    pop_idx_list = [pop_list.index(pop) for pop in pops if pop in pop_list]
    print(f'pop_idx_list: {pop_idx_list}')
    
    
    field_rename_dict = {'AF.Cases':'AF_cases',
                         'AF.Controls':'AF_controls',
                         'BETA':'beta',
                         'SE':'se',
                         'Pvalue':'pval',
                        }
    meta_field_rename_dict = {'BETA':'beta_meta',
                              'SE':'se_meta',
                              'Pvalue':'pval_meta',
                              'AF_Cases':'AF_cases_meta',
                              'AF_Controls':'AF_controls_meta',
                              'Q':'Q',
                              'Pvalue_het':'pval_het'}
    
    annotate_dict = {}
    
    for field in ['AF.Cases','AF.Controls','BETA','SE','Pvalue']:
        for pop_idx in pop_idx_list:
            annotate_dict.update({f'{field_rename_dict[field]}_{pops[pop_idx]}': hl.if_else(hl.is_nan(mt1.summary_stats[field][pop_idx]),
                                                                                        hl.str(mt1.summary_stats[field][pop_idx]),
                                                                                        hl.format('%.3e', mt1.summary_stats[field][pop_idx]))})
    
    keyed_mt = meta_mt0[(mt1.locus,mt1.alleles),(mt1.pheno, mt1.coding, mt1.trait_type)]
    for field in ['AF_Cases','AF_Controls','BETA','SE','Pvalue','Pvalue_het']:
            annotate_dict.update({f'{meta_field_rename_dict[field]}': hl.if_else(hl.is_nan(keyed_mt.meta_analysis[field][0]),
                                                                                hl.str(keyed_mt.meta_analysis[field][0]),
                                                                                hl.format('%.3e', keyed_mt.meta_analysis[field][0]))})
    
    mt2 = mt1.annotate_entries(**annotate_dict)
    
    mt2 = mt2.key_cols_by('pheno_id')
    mt2 = mt2.key_rows_by().drop('locus','alleles','gene','annotation','summary_stats')
    
    batch_idx = 1
    export_out = f'{ldprune_dir}/release/batch{batch_idx}'
    while hl.hadoop_is_dir(export_out):
        batch_idx += 1
        export_out = f'{ldprune_dir}/release/batch{batch_idx}'
    print(f'\nSaving to: {export_out}\n')
    hl.experimental.export_entries_by_col(mt = mt2,
                                          path = export_out,
                                          bgzip = True,
                                          batch_size = 256,
                                          use_string_key_as_file_name = True,
                                          header_json_in_file = False)

    
def export_pop_pheno_pairs():
    hl.init(default_reference='GRCh38', log='/tmp/export_entries_by_col.log')
    mt0 = hl.read_matrix_table('gs://ukb-diverse-pops/combined_results/results_full.mt')
#    meta_mt0 = hl.read_matrix_table('gs://ukb-diverse-pops/combined_results/meta_analysis.mt')
    mt0 = mt0.annotate_cols(pheno_id = mt0.trait_type+'-'+mt0.pheno+'-'+mt0.coding)
    mt0 = mt0.annotate_rows(variant = mt0.locus.contig+':'+hl.str(mt0.locus.position)+':'+mt0.alleles[0]+':'+mt0.alleles[1])

    mt = mt0.filter_cols((hl.literal(['biomarkers']).contains(mt0.trait_type))&(hl.len(mt0.pheno_data.pop)==6))    

    pop_list = mt.pheno_data.pop.take(1)[0]
    
    ldprune_dir = f'{bucket}/ld_prune'
    export_out = f'{ldprune_dir}/release-test'
    
    pheno_id_list = mt.pheno_id.collect()
    
    for pop in pop_list:
        ht_pop = hl.read_table(f'gs://ukb-diverse-pops/combined_results/results_{pop}.ht')
        rename_dict = {'AF_Allele2':'AF',
                       'BETA':'beta',
                       'SE':'se',
                       'Pvalue':'pval'}
        ht_pop = ht_pop.rename(rename_dict)
        ht_pop = ht_pop.annotate(pheno_id = ht_pop.trait_type+'-'+ht_pop.pheno+'-'+ht_pop.coding)
        ht_pop = ht_pop.annotate(variant = ht_pop.locus.contig+':'+hl.str(ht_pop.locus.position)+':'+ht_pop.alleles[0]+':'+ht_pop.alleles[1])
        for pheno_id in pheno_id_list:
            ht_pheno = ht_pop.filter(ht_pop.pheno_id==pheno_id)    
            ht_pheno = ht_pheno.key_by('variant').select('AF','beta','se','pval')
            if ht_pheno.count()>0:
                ht_pheno.export(f'{export_out}/{pop}.{pheno_id}.tsv.bgz')
        
def export_loo():
    meta_mt0 = hl.read_matrix_table('gs://ukb-diverse-pops/combined_results/meta_analysis.mt')
    
    meta_mt0 = meta_mt0.annotate_cols(pheno_id = meta_mt0.trait_type+'-'+meta_mt0.pheno+'-'+meta_mt0.coding)
    
    meta_mt0 = meta_mt0.annotate_rows(SNP = meta_mt0.locus.contig+':'+hl.str(meta_mt0.locus.position)+':'+meta_mt0.alleles[0]+':'+meta_mt0.alleles[1],
                                      A1 = meta_mt0.alleles[1],
                                      A2 = meta_mt0.alleles[0])
    
    meta_field_rename_dict = {'BETA':'b',
                              'SE':'se',
                              'Pvalue':'p',
                              'AF_Allele2':'freq'}
    
    loo = meta_mt0.filter_cols(hl.len(meta_mt0.pheno_data.pop)==5)
    for pop in ['AFR','EUR']: #['AFR','AMR','CSA','EAS','EUR','MID']:
        print(f'not_{pop}')
        loo_pop = loo.filter_cols(~loo.pheno_data.pop.contains(pop))
        annotate_dict = {}
        for field in ['AF_Allele2','BETA','SE','Pvalue']:
            annotate_dict.update({meta_field_rename_dict[field]: loo_pop.meta_analysis[field][0]}) # first entry in array always is the meta-analysis of 5 populations
        loo_pop = loo_pop.annotate_entries(**annotate_dict)
    
        loo_pop = loo_pop.key_cols_by('pheno_id')
        loo_pop = loo_pop.key_rows_by().drop('locus','alleles','gene','annotation','meta_analysis')
        
        batch_idx = 1
        export_out = f'{ldprune_dir}/release/batch{batch_idx}'
        while hl.hadoop_is_dir(export_out):
            batch_idx += 1
            export_out = f'{ldprune_dir}/release/batch{batch_idx}'
        print(f'\nSaving to: {export_out}\n')
        hl.experimental.export_entries_by_col(mt = loo_pop,
                                              path = export_out,
                                              bgzip = True,
                                              batch_size = 256,
                                              use_string_key_as_file_name = True,
                                              header_json_in_file = False)
        
#    all_pops = meta_mt0.filter_cols(hl.len(meta_mt0.pheno_data.pop)==6)
#    for pop_idx, pop in enumerate(['AFR','AMR','CSA','EAS','EUR','MID'],1):
#        if pop in ['AFR','EUR']:
#            annotate_dict = {}
#            for field in ['AF_Allele2','BETA','SE','Pvalue']:
#                annotate_dict.update({meta_field_rename_dict[field]: all_pops.meta_analysis[field][pop_idx]}) # 2nd entry corresponds to not_AFR, 3rd corresponds to not_AMR, etc.
#            loo_pop = all_pops.annotate_entries(**annotate_dict)
#            loo_pop = loo_pop.key_cols_by('pheno_id')
#            loo_pop = loo_pop.key_rows_by().drop('locus','alleles','gene','annotation','meta_analysis')
#            batch_idx = 1
#            export_out = f'{ldprune_dir}/release/batch{batch_idx}'
#            while hl.hadoop_is_dir(export_out):
#                batch_idx += 1
#                export_out = f'{ldprune_dir}/release/batch{batch_idx}'
#            print(f'\nSaving to: {export_out}\n')
#            hl.experimental.export_entries_by_col(mt = loo_pop,
#                                                  path = export_out,
#                                                  bgzip = True,
#                                                  batch_size = 256,
#                                                  use_string_key_as_file_name = True,
#                                                  header_json_in_file = False)

def make_pheno_master():
    mt0 = hl.read_matrix_table('gs://ukb-diverse-pops/combined_results/results_full.mt')
    mt0_cols = mt0.cols()
    annotate_dict = {}
    for field in ['n_cases','n_controls','saige_heritability']:
        for pop in ['AFR','AMR','CSA','EAS','EUR','MID']:
            annotate_dict.update({f'{field}_{pop}': hl.if_else(hl.is_nan(mt0_cols.pheno_data.pop.index(pop)),
                                                               hl.null(mt0_cols.pheno_data[field][0].dtype),
                                                               mt0_cols.pheno_data[field][mt0_cols.pheno_data.pop.index(pop)])})
    # print(annotate_dict)
    mt0_cols = mt0_cols.annotate(**annotate_dict)
    mt0_cols = mt0_cols.drop('pheno_data')
    mt0_cols.export(f'{ldprune_dir}/phenotype_master.tsv.bgz')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--n_threads', help='number of threads')
    parser.add_argument('--input_file', help='Input file of variant results')
    parser.add_argument('--output_file', help='Output file of variant results')
    parser.add_argument('--pop', type=str, help='Population to be left out')
    parser.add_argument('--pheno', type=str, help='Phenotype of meta-analyzed sumstats')
    parser.add_argument('--coding', type=str, help='coding of meta-analyzed sumstats')
    parser.add_argument('--trait_type', type=str, help='trait_type of  of meta-analyzed sumstats')
    parser.add_argument('--ht_to_tsv', action='store_true')
    parser.add_argument('--tsv_to_ht', action='store_true')
    parser.add_argument('--write_meta_sumstats', action='store_true')
    parser.add_argument('--get_meta_sumstats', action='store_true')
    parser.add_argument('--test_get_meta_sumstats', action='store_true')
    parser.add_argument('--export_cont_results', action='store_true')
    parser.add_argument('--export_binary_results', action='store_true')
    parser.add_argument('--export_pop_pheno_pairs', action='store_true')
    parser.add_argument('--export_loo', action='store_true')
    parser.add_argument('--overwrite', default=False, action='store_true', help='overwrite existing files')
    args = parser.parse_args()

    if args.ht_to_tsv:
        ht_to_tsv(args)
    elif args.tsv_to_ht:
        tsv_to_ht(args)
    elif args.get_meta_sumstats:
        get_meta_sumstats(args)
    elif args.test_get_meta_sumstats:
        test_get_meta_sumstats(args)
    elif args.write_meta_sumstats:
        write_meta_sumstats(args)
    elif args.export_cont_results:
        export_cont_results()
    elif args.export_binary_results:
        export_binary_results()
    elif args.export_pop_pheno_pairs:
        export_pop_pheno_pairs()
    elif args.export_loo:
        export_loo()

        

        