#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 13:06:01 2020

Flatten Dropbox manifest

@author: nbaya
"""

import pandas as pd

wd = '/Users/nbaya/Documents/lab/diverse_pops'

df = pd.read_csv(f'{wd}/UKBB_Pan_Populations-Manifest-20200615-manifest_info.tsv',
                 sep='\t')

bgz = df[~df.File.str.contains('tbi')]
bgz = bgz.rename(columns={'File':'filename'})
tbi = df[df.File.str.contains('tbi')]
tbi['filename'] = tbi.File.str.replace('.tbi','')

merged = bgz.merge(tbi, on='filename')

