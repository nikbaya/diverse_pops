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
import ukb_common
#from ukbb_pan_ancestry import get_pheno_manifest_path
from ukb_common import mwzj_hts_by_tree

bucket = 'gs://ukb-diverse-pops'
ldprune_dir = f'{bucket}/ld_prune'


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
    ht = ht.filter(ht.contig!='') # added because there were some files with extra lines with empty strings
    ht = ht.key_by(locus = hl.locus(contig=ht.contig, 
                                    pos=hl.int(ht.pos),
                                    reference_genome='GRCh37'),
                   alleles = hl.array([ht.varid.split(':')[2],
                                       ht.varid.split(':')[3]])
    )
    ht = ht.annotate_globals(trait_type = args.trait_type,
                             phenocode = args.phenocode,
                             pheno_sex = args.pheno_sex,
                             coding = args.coding,
                             modifier = args.modifier)
    ht = ht.drop('contig','varid','pos')
#    print(ht.describe())
#    ht.select().show()
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
    
#    mt = mt0.filter_cols((hl.len(mt0.meta_analysis_data[0].pop)==6)|
#            (hl.literal(PILOT_PHENOTYPES).contains((mt0.pheno, mt0.coding, mt0.trait_type))))
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
    r'''
    '''
    print(ldprune_dir)
    print(f'\n...Getting meta-analysis sumstats...')
    all_pops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']
    pop = args.pop
    add_args = {}
    if args.n_threads is not None:
        add_args['master'] = f'local[{args.n_threads}]'
    hl.init(default_reference='GRCh38', log='/tmp/mt_to_tsv.log', **add_args)
    meta_mt0 = hl.read_matrix_table('gs://ukb-diverse-pops/combined_results/meta_analysis.mt')
    trait_type, phenocode, pheno_sex, coding, modifier = "biomarkers", "30600", "both_sexes", "30600", ""
    meta_mt1 = meta_mt0.filter_cols((meta_mt0.trait_type==trait_type)&
                                    (meta_mt0.phenocode==phenocode)&
                                    (meta_mt0.pheno_sex==pheno_sex)&
                                    (meta_mt0.coding==coding)&
                                    (meta_mt0.modifier==modifier))
    
    req_pop_list = [p for p in all_pops if p is not pop]
    meta_mt1 = meta_mt1.annotate_cols(idx = meta_mt1.meta_analysis_data.pop.index(hl.literal(req_pop_list))) # get index of which meta-analysis is the leave-on-out for current pop
    loo_mt = meta_mt1.filter_cols(hl.is_defined(meta_mt1.idx))
    loo_ht0 = loo_mt.entries()
    loo_ht0 = loo_ht0.key_by()
    loo_ht1 = loo_ht0.select(pval = loo_ht0.meta_analysis.Pvalue[loo_ht0.idx],
                             SNP = loo_ht0.locus.contig+':'+hl.str(loo_ht0.locus.position)+':'+loo_ht0.alleles[0]+':'+loo_ht0.alleles[1])
    loo_ht1 = loo_ht1.filter(hl.is_defined(loo_ht1.pval))
    loo_ht1.export(args.output_file)
   
        
def export_ma_format(batch_size=256):
    r'''
    Export columns for .ma format (A1, A2, freq, beta, se, N) for select phenotypes
    '''
    meta_mt0 = hl.read_matrix_table('gs://ukb-diverse-pops/combined_results/meta_analysis.mt')
    
    highprev = hl.import_table(f'{ldprune_dir}/joined_ukbb_lancet_age_high_prev.tsv', impute=True)
    highprev = highprev.annotate(pheno = highprev.code.replace('_irnt',''))
    pheno_list = highprev.pheno.collect()
    pheno_list = [p for p in pheno_list if p is not None]
    meta_mt0 = meta_mt0.filter_cols(hl.literal(pheno_list).contains(meta_mt0.pheno))

    meta_mt0 = meta_mt0.annotate_cols(pheno_id = (meta_mt0.trait_type+'-'+
                                      meta_mt0.phenocode+'-'+
                                      meta_mt0.pheno_sex+
                                      hl.if_else(hl.len(meta_mt0.coding)>0, '-'+meta_mt0.coding, '')+
                                      hl.if_else(hl.len(meta_mt0.modifier)>0, '-'+meta_mt0.modifier, '')
                                      ).replace(' ','_').replace('/','_'))
    
    meta_mt0 = meta_mt0.annotate_rows(SNP = meta_mt0.locus.contig+':'+hl.str(meta_mt0.locus.position)+':'+meta_mt0.alleles[0]+':'+meta_mt0.alleles[1],
                                      A1 = meta_mt0.alleles[1], # .ma format requires A1 = effect allele, which in this case is A2 for UKB GWAS
                                      A2 = meta_mt0.alleles[0])

    meta_field_rename_dict = {'BETA':'b',
                          'SE':'se',
                          'Pvalue':'p',
                          'AF_Allele2':'freq',
                          'N':'N'}
    
    all_pops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']

    for pop in ['AFR','EUR']: #['AFR','AMR','CSA','EAS','EUR','MID']:
        print(f'not_{pop}')

        req_pop_list = [p for p in all_pops if p is not pop]

        loo_pop = meta_mt0.annotate_cols(idx = meta_mt0.meta_analysis_data.pop.index(hl.literal(req_pop_list))) # get index of which meta-analysis is the leave-on-out for current pop
        loo_pop = loo_pop.filter_cols(hl.is_defined(loo_pop.idx))
        
        annotate_dict = {}
        for field in ['AF_Allele2','BETA','SE','Pvalue','N']:
            annotate_dict.update({meta_field_rename_dict[field]: loo_pop.meta_analysis[field][loo_pop.idx]}) 
        batch_idx = 1

    export_out = f'{ldprune_dir}/loo/not_{pop}/batch{batch_idx}'
    while hl.hadoop_is_dir(export_out):
        batch_idx += 1
        export_out = f'{ldprune_dir}/loo/not_{pop}/batch{batch_idx}'
    checkpoint_path = f'gs://ukbb-diverse-temp-30day/loo/not_{pop}/batch{batch_idx}.mt'
#        print(f'\nCheckpointing to: {checkpoint_path}\n')
    loo_pop = loo_pop.checkpoint(checkpoint_path,
                                 _read_if_exists=True,
                                 overwrite=True)
    loo_pop = loo_pop.filter_entries(hl.is_defined(loo_pop.b))
    print(f'\nExporting to: {export_out}\n')
    hl.experimental.export_entries_by_col(mt = loo_pop,
                                          path = export_out,
                                          bgzip = True,
                                          batch_size = batch_size,
                                          use_string_key_as_file_name = True,
                                          header_json_in_file = False)
    

#def mwzj_hts_by_tree(all_hts, temp_dir, globals_for_col_key, 
#                           debug=False, inner_mode = 'overwrite', repartition_final: int = None):
#    r'''
#    Adapted from ukb_common mwzj_hts_by_tree()
#    Uses read_clump_ht() instead of read_table()
#    '''
#    chunk_size = int(len(all_hts) ** 0.5) + 1
#    outer_hts = []
#    
#    checkpoint_kwargs = {inner_mode: True}
#    if repartition_final is not None:
#        intervals = ukb_common.get_n_even_intervals(repartition_final)
#        checkpoint_kwargs['_intervals'] = intervals
#    
#    if debug: print(f'Running chunk size {chunk_size}...')
#    for i in range(chunk_size):
#        if i * chunk_size >= len(all_hts): break
#        hts = all_hts[i * chunk_size:(i + 1) * chunk_size]
#        if debug: print(f'Going from {i * chunk_size} to {(i + 1) * chunk_size} ({len(hts)} HTs)...')
#        try:
#            if isinstance(hts[0], str):
#                hts = list(map(lambda x: hl.read_table(x), hts))
#            ht = hl.Table.multi_way_zip_join(hts, 'row_field_name', 'global_field_name')
#        except:
#            if debug:
#                print(f'problem in range {i * chunk_size}-{i * chunk_size + chunk_size}')
#                _ = [ht.describe() for ht in hts]
#            raise
#        outer_hts.append(ht.checkpoint(f'{temp_dir}/temp_output_{i}.ht', _read_if_exists=True, **checkpoint_kwargs))
#    ht = hl.Table.multi_way_zip_join(outer_hts, 'row_field_name_outer', 'global_field_name_outer')
#    ht = ht.transmute(inner_row=hl.flatmap(lambda i:
#                                           hl.cond(hl.is_missing(ht.row_field_name_outer[i].row_field_name),
#                                                   hl.range(0, hl.len(ht.global_field_name_outer[i].global_field_name))
#                                                   .map(lambda _: hl.null(ht.row_field_name_outer[i].row_field_name.dtype.element_type)),
#                                                   ht.row_field_name_outer[i].row_field_name),
#                                           hl.range(hl.len(ht.global_field_name_outer))))
#    ht = ht.transmute_globals(inner_global=hl.flatmap(lambda x: x.global_field_name, ht.global_field_name_outer))
#    mt = ht._unlocalize_entries('inner_row', 'inner_global', globals_for_col_key)
#    return mt

def resume_mwzj(temp_dir, globals_for_col_key):
    ls = hl.hadoop_ls(temp_dir)
    paths = [x['path'] for x in ls if 'temp_output' in x['path'] ]
    chunk_size = len(paths)
    outer_hts = []
    for i in range(chunk_size):
        outer_hts.append(hl.read_table(f'{temp_dir}/temp_output_{i}.ht'))
    ht = hl.Table.multi_way_zip_join(outer_hts, 'row_field_name_outer', 'global_field_name_outer')
    ht = ht.transmute(inner_row=hl.flatmap(lambda i:
                                           hl.cond(hl.is_missing(ht.row_field_name_outer[i].row_field_name),
                                                   hl.range(0, hl.len(ht.global_field_name_outer[i].global_field_name))
                                                   .map(lambda _: hl.null(ht.row_field_name_outer[i].row_field_name.dtype.element_type)),
                                                   ht.row_field_name_outer[i].row_field_name),
                                           hl.range(hl.len(ht.global_field_name_outer))))
    ht = ht.transmute_globals(inner_global=hl.flatmap(lambda x: x.global_field_name, ht.global_field_name_outer))
    mt = ht._unlocalize_entries('inner_row', 'inner_global', globals_for_col_key)
    return mt

def join_clump_results(pop):
    r'''
    Wrapper for mwzj
    '''
    pop = pop.upper()
    pheno_manifest = hl.import_table(
    #            get_pheno_manifest_path(), 
            'gs://ukb-diverse-pops/ld_prune/phenotype_manifest.tsv.bgz', # hardcoded path to avoid having to change user-pays
            impute=True, 
            key=ukb_common.PHENO_KEY_FIELDS)
    pheno_manifest = pheno_manifest.annotate(pheno_id = pheno_manifest.filename.replace('.tsv.bgz',''))
    
    clump_results_dir = f'{ldprune_dir}/results/not_{pop}'
    ls = hl.hadoop_ls(f'{clump_results_dir}/*')
    all_hts = [x['path'] for x in ls if 'clump_results.ht' in x['path']]
    
    temp_dir = f'gs://ukbb-diverse-temp-30day/nb-temp/{pop}'
    globals_for_col_key = ukb_common.PHENO_KEY_FIELDS
    mt = mwzj_hts_by_tree(all_hts=all_hts,
                         temp_dir=temp_dir,
                         globals_for_col_key=globals_for_col_key)
#    mt = resume_mwzj(temp_dir=temp_dir, # NOTE: only use if all the temp hts have been created
#                     globals_for_col_key=globals_for_col_key)
    mt.write(f'{ldprune_dir}/clump_results/not_{pop}.mt')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--n_threads', help='number of threads')
    parser.add_argument('--input_file', help='Input file of variant results')
    parser.add_argument('--output_file', help='Output file of variant results')
    parser.add_argument('--pop', type=str, help='Population to be left out')
    parser.add_argument('--trait_type', type=str, help='trait_type in meta-analyzed sumstats')
    parser.add_argument('--phenocode', type=str, help='phenocode in meta-analyzed sumstats')
    parser.add_argument('--pheno_sex', type=str, help='pheno_sex in meta-analyzed sumstats')
    parser.add_argument('--coding', type=str, default='', help='coding in meta-analyzed sumstats')
    parser.add_argument('--modifier', type=str, default='', help='modifier in meta-analyzed sumstats')
    parser.add_argument('--ht_to_tsv', action='store_true')
    parser.add_argument('--tsv_to_ht', action='store_true')
    parser.add_argument('--write_meta_sumstats', action='store_true')
    parser.add_argument('--get_meta_sumstats', action='store_true')
    parser.add_argument('--test_get_meta_sumstats', action='store_true')
    parser.add_argument('--export_pop_pheno_pairs', action='store_true')    
    parser.add_argument('--join_clump_results', action='store_true')    
    parser.add_argument('--batch_size', type=int, default=256, help='max number of phenotypes per batch for export_entries_by_col')
    parser.add_argument('--overwrite', default=False, action='store_true', help='overwrite existing files')
    args = parser.parse_args()
    
#    try:        
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
    elif args.join_clump_results:
        join_clump_results(pop=args.pop)

#    except:
#        hl.copy_log('gs://ukbb-diverse-temp-30day/nb_hail.log')
        

        