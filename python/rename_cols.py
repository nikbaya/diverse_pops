#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 13:26:26 2020

Rename columns of exported flat files

@author: nbaya
"""

import subprocess
from itertools import combinations
import hail as hl
hl.init()

for num_pops in range(1,6)[::-1]:
    
    all_pops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']
    pop_sets = [set(i) for i in list(combinations(all_pops, num_pops))] # list of exact set of pops for which phenotype is defined

    incorrrect_pops = all_pops[:num_pops] # incorrect column suffixes

    for trait_category in ['quant','binary']:
        
        for pop_set in pop_sets:
            
            pop_list = sorted(pop_set) # correct column suffixes
            
            bucket = f'gs://ukb-diverse-pops/ld_prune/release/{trait_category}/{"-".join(pop_list)}_batch1'
            print(bucket)
            
            if hl.hadoop_is_dir(bucket):
                subprocess.call(['gsutil','-m','cp',f'{bucket}/*bgz', './'])
                assert False
#                break
