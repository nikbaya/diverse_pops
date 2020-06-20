#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 11:10:23 2020

Fixes issue of A1 in varid being minor allele instead of A1 in UKB.

@author: nbaya
"""

import argparse
import hail as hl

#wd = 'file:///home/nbaya'
#cloud_wd = f'gs://ukb-diverse-pops/ld_prune/subsets_5k'
ldprune_dir = 'gs://ukb-diverse-pops/ld_prune'

def main(args):
    variants = hl.import_table(f'{ldprune_dir}/variants.txt', delimiter=' ')
    variants = variants.annotate(alleles = hl.set([variants.ref, variants.alt]))
    variants = variants.key_by('chrom','pos','alleles')
    
    
    all_pops = [
            'AFR', 
            'AMR', 
            'CSA', 
            'EAS', 
            'EUR', 
            'MID'
            ]

    pops = [args.pop.upper()] if args.pop!=None else all_pops
    
    for pop in pops:
        bim0 = hl.import_table(f'{ldprune_dir}/subsets_5k/not_{pop}/not_{pop}.bim', 
                              delimiter=' ', 
                              no_header=True)
        bim0 = bim0.rename({'f0':'chrom',
                          'f3':'pos',
                          'f4':'A1_bim',
                          'f5':'A2_bim'})
        bim0 = bim0.add_index()

        bim = bim0.annotate(alleles = hl.set([bim0.A1_bim, bim0.A2_bim]))
        bim = bim.key_by('chrom','pos','alleles')
        
        bim = bim.annotate(ref = variants[bim.key].ref,
                           alt = variants[bim.key].alt)
        bim = bim.annotate(SNP = bim.chrom+':'+bim.pos+':'+bim.ref+':'+bim.alt)
        bim = bim.key_by()
        bim = bim.order_by(bim.idx)
        bim = bim.select('chrom','SNP','f2','pos','A1_bim','A2_bim')
        bim.export(f'{ldprune_dir}/subsets_5k/not_{pop}/not_{pop}.bim_new', 
                   header=False, 
                   delimiter=' ')

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pop', default=None, type=str, help='population to use')
    args = parser.parse_args()

    main(args=args)