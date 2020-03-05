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
    Convert Hail table of variant results to a gzipped tsv 
    '''
    add_args = {}
    if args.n_threads is not None:
        add_args['master'] = f'local[{args.n_threads}]'
    hl.init(default_reference='GRCh38', log='/tmp/plink_clump_hail.log', **add_args)
    ht = hl.read_table(args.input_file)
    ht.export(args.output_file)
    
def tsv_to_ht(args):
    r'''
    Convert tsv of variant results to a Hail table
    '''
    add_args = {}
    if args.n_threads is not None:
        add_args['master'] = f'local[{args.n_threads}]'
    hl.init(default_reference='GRCh38', log='/tmp/plink_clump_hail.log', **add_args)
    ht = hl.import_table(args.input_file)
    ht.write(args.output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--n_threads', help='Overwrite everything')
    parser.add_argument('--input_file', help='Input file of variant results')
    parser.add_argument('--output_file', help='Output file of variant results')
    parser.add_argument('--from_ht', action='store_true')
    args = parser.parse_args()

    if args.from_ht:
        ht_to_tsv(args)
    else:
        tsv_to_ht(args)

        