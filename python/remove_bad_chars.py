#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 14:36:36 2020

Remove spaces from filenames of existing files

@author: nbaya
"""

from hail.utils import hadoop_ls
import subprocess

release_dir = 'gs://ukb-diverse-pops/ld_prune/release'

# list of buckets containing files to be changed
buckets= ['quant_batch2',
          'binary_batch3']


for bucket in buckets:
    ls = hadoop_ls(f'{release_dir}/{bucket}')
    
    paths = [f['path'] for f in ls if (' ' in f['path'])] # paths of files with ' ' in path
    
    assert ~any([f['is_dir'] for f in ls]), 'WARNING: Directories exist within this bucket'
    
    for path in paths:
        cmd = ['gsutil', 'mv', path, path.replace(" ","_")]
#        print(cmd)
        subprocess.run(cmd)