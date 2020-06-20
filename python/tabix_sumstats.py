#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 17:33:14 2020

Tabix sumstats tsvs using Hail Batch

@author: nbaya
"""

import hail as hl
import hailtop.batch as hb

bucket = 'gs://ukb-diverse-pops'
ldprune_dir = f'{bucket}/ld_prune'
scratch_dir = f'gs://ukbb-diverse-temp-30day/nb-tabix-scratch'

def get_ss_path_list(sumstats_dir):
    ss_list = hl.hadoop_ls(sumstats_dir)
    ss_path_list = [x['path'] for x in ss_list if 'bgz' in x['path']]
    print(f'\nNumber of sumstats files: {len(ss_path_list)}\n')
    return ss_path_list
    
def tabix(b, ss_path, out_dir):
    r'''
    tabix's a bgz file with gcloud path `path` using Batch `b`
    '''
    fname = ss_path.split('/')[-1]
    f = b.read_input(ss_path)
    j = b.new_job(name=fname.split('.')[0])
    j.command(f'tabix -s 1 -b 2 -e 2 -c chr {f}') # treat header (which begins with "chr") as a comment
    j.command(f'mv {f}.tbi {j.ofile}')
    b.write_output(j.ofile, f'{out_dir}/{fname}.tbi')

if __name__=="__main__":
    hl.init(log='/Users/nbaya/Downloads/tabix_sumstats.log')
    backend = hb.ServiceBackend(billing_project='ukb_diverse_pops',
                                bucket='ukbb-diverse-temp-30day/nb-batch-tmp')
    
    b = hb.Batch(name='tabix', backend=backend,
                 default_image='gcr.io/ukbb-diversepops-neale/nbaya_tabix:latest',
                 default_storage='100M', # works with 2G
                 default_cpu=1)
    
#    sumstats_dir = f'{bucket}/sumstats_flat_files'
#    sumstats_dir = f'{ldprune_dir}/export_results/update'
#    sumstats_dir = f'{ldprune_dir}/loo/sumstats/batch1'
    sumstats_dir = f'{ldprune_dir}/variant_qc'
    print(f'\nUsing sumstats from {sumstats_dir}')
    
    ss_path_list = get_ss_path_list(sumstats_dir=sumstats_dir)
    
    out_dir = f'{sumstats_dir}_tabix'
    print(f'\nSaving tabix files to {out_dir}\n')
    
    for ss_path in ss_path_list:
        tabix(b=b, 
              ss_path=ss_path, 
              out_dir=out_dir)
        
    b.run(open=True)
        