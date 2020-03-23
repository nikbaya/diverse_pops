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
ldprune_wd = f'{bucket}/ld_prune'

def ht_to_tsv(args):
    r'''
    Convert Hail table of variant results to a tsv 
    '''
    add_args = {}
    if args.n_threads is not None:
        add_args['master'] = f'local[{args.n_threads}]'
    hl.init(default_reference='GRCh38', log='/tmp/ht_to_tsv.log', **add_args)
    ht = hl.read_table(args.input_file)
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
    mt = hl.read_matrix_table('gs://ukb-diverse-pops/combined_results/meta_analysis.mt')
    mt = mt.filter_cols(hl.len(mt.meta_analysis_data[0].pop)==6)  # get phenos with results for all 6 populations
    mt2 = mt.select_cols().select_rows()
    mt2 = mt2.select_entries(Pvalue = mt2.meta_analysis.Pvalue)
    ht1 = mt2.entries()
    ht1 = ht1.key_by('pheno','trait_type','coding','locus','alleles')
    ht1.describe()
    ht1 = ht1.annotate(varid = hl.str(ht1.locus)+'_'+ht1.alleles[0]+'_'+ht1.alleles[1])
    out = f'{ldprune_wd}/meta_analysis.all_pops.ht'
    print(f'\n...Writing meta-analysis sumstats for all pops...\nFrom: {args.input_file}\nTo:{out}')
    ht1.write(out, overwrite=True)

def get_meta_sumstats(args):
    print(ldprune_wd)
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
    parser.add_argument('--overwrite', default=False, action='store_true', help='overwrite existing files')
    args = parser.parse_args()

    if args.ht_to_tsv:
        ht_to_tsv(args)
    elif args.tsv_to_ht:
        tsv_to_ht(args)
    elif args.get_meta_sumstats:
        get_meta_sumstats(args)
    elif args.write_meta_sumstats:
        write_meta_sumstats(args)

        