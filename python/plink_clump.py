#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 09:35:41 2020

Clumping GWAS results with PLINK

@author: nbaya
"""

import hailtop.pipeline as hp

ldprune_wd = 'gs://ukb-diverse-pops/ld_prune/'



def plink_clump(phen, pop, not_pop=False):
    

    backend = hp.BatchBackend(billing_project='ukb_diverse_pops')
    p = hp.Pipeline(name='plink_clump', backend=backend,
                          default_image='gcr.io/ukbb-diversepops-neale/nbaya_plink:latest',
                          default_storage='500Mi', default_cpu=8)
    


def local_plink_clump():
    backend = hp.LocalBackend(gsa_key_file='/Users/nbaya/.hail/ukb-diverse-pops.json')
    p = hp.Pipeline(name='plink_clump', backend=backend)

    
    create = p.new_task(name='create-dummy')
    create.declare_resource_group(bfile={'bed': '{root}.bed',
                                         'bim': '{root}.bim',
                                         'fam': '{root}.fam'})
    
#create.command(f'plink --dummy 10 100 --make-bed --out {create.bfile}')
create.command(f'/Users/nbaya/Documents/lab/smiles/smiles/bash/plink/plink --dummy 10 100 --make-bed --out {create.bfile}')

p.run(open=True) 




#bfile_root = f'{ldprune_wd}/subsets/{"not_" if not_pop else ""}{pop}'
bfile_root = '/Users/nbaya/Documents/lab/smiles/data/arabidopsis1'
bfile = p.read_input_group(bed=f'{bfile_root}.bed',
                           bim=f'{bfile_root}.bim',
                           fam=f'{bfile_root}.fam')

wc_bim = p.new_task(name='wc-bim')
wc_bim.command(f'wc -l {bfile.bim}')
p.run()

    clump = p.new_task(name='clump')
    
    clump.declare_resource_group(bfile={'bed': f'{root}.bed',
                                        'bim': f'{root}.bim',
                                        'fam': f'{root}.fam'})
    clump.command(f'plink --dummy 10 100 --make-bed --out {clump.bfile}')
    p.run() 
    
    pass