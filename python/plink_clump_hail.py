#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 15:09:04 2020

Hail script for clumping GWAS results with PLINK

@author: nbaya
"""

import argparse
import hail as hl


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
    Convert tsv of variant results to a Hail table
    '''
    add_args = {}
    if args.n_threads is not None:
        add_args['master'] = f'local[{args.n_threads}]'
    hl.init(default_reference='GRCh38', log='/tmp/tsv_to_ht.log', **add_args)
    ht = hl.import_table(args.input_file)
    ht.write(args.output_file)

def get_meta_sumstats(args):
    add_args = {}
    if args.n_threads is not None:
        add_args['master'] = f'local[{args.n_threads}]'
    hl.init(default_reference='GRCh38', log='/tmp/mt_to_tsv.log', **add_args)
    mt = hl.read_matrix_table(args.input_file)
    mt = mt.filter_cols(hl.len(mt.meta_analysis_data[0].pop)==6)  # get phenos with results for all 6 populations
    mt_pheno = mt.filter_cols((mt.pheno==args.pheno)&
                              (mt.coding==args.coding)&
                              (mt.trait_type==args.trait_type))
    ht_pheno = mt_pheno.key_cols_by('pheno').select_entries('meta_analysis').make_table()
    pops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']
    loo_meta_idx = pops.index(args.pop)+1 # get index of leave-one-out meta-analysis result for population `pop`
    ht_pheno = ht_pheno.select(Pvalue=ht_pheno[f'{args.pheno}.meta_analysis'].Pvalue[loo_meta_idx])
    ht_pheno = ht_pheno.filter(hl.is_defined(ht_pheno.Pvalue))
    ht_pheno = ht_pheno.annotate(varid = ht_pheno.locus+'_'+ht_pheno.alleles[0]+'_'+ht_pheno.alleles[1]) # format: "chr:pos_ref_alt"
    ht_pheno.key_by('varid')
    ht_pheno.select('Pvalue').write(args.output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--n_threads', help='Overwrite everything')
    parser.add_argument('--input_file', help='Input file of variant results')
    parser.add_argument('--output_file', help='Output file of variant results')
    parser.add_argument('--pop', type=str, help='Population to be left out')
    parser.add_argument('--pheno', type=str, help='Phenotype of meta-analyzed sumstats')
    parser.add_argument('--coding', type=str, help='coding of meta-analyzed sumstats')
    parser.add_argument('--trait_type', type=str, help='trait_type of  of meta-analyzed sumstats')
    parser.add_argument('--ht_to_tsv', action='store_true')
    parser.add_argument('--tsv_to_ht', action='store_true')
    parser.add_argument('--get_meta_sumstats', action='store_true')
    args = parser.parse_args()

    if args.ht_to_tsv:
        ht_to_tsv(args)
    elif args.tsv_to_ht:
        tsv_to_ht(args)
    elif args.get_meta_sumstats:
        get_meta_sumstats(args)

        