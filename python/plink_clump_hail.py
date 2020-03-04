#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 15:09:04 2020

Hail script for Clumping GWAS results with PLINK

@author: nbaya
"""

import argparse
import hail as hl


def ht_to_tsv(args):
    r'''
    Convert Hail table of variant results to a gzipped tsv for use in the pipeline
    '''
    add_args = {}
    if args.n_threads is not None:
        add_args['master'] = f'local[{args.n_threads}]'
    hl.init(default_reference='GRCh38', log='/tmp/plink_clump_hail.log', **add_args)

    ht = hl.read_table(args.input_file)
    
    
    
    ht.export(args.output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--n_threads', help='Overwrite everything')
    parser.add_argument('--input_file', help='Input Hail table of variant results')
    parser.add_argument('--output_file', help='Output tsv (gzipped) of variant results')
    args = parser.parse_args()

    ht_to_tsv(args)