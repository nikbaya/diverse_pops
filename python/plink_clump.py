#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 09:35:41 2020

Clumping GWAS results with PLINK

@author: nbaya
"""

import hailtop.pipeline as hp

bucket = 'gs://ukb-diverse-pops'
ldprune_wd = f'{bucket}/ld_prune'

def get_variant_results_path(pop: str, extension: str = 'ht'):
    return f'{bucket}/combined_results/results_{pop}.{extension}'



def declare_rg(t, name):
    t.declare_resource_group(**{name: {'bed': '{root}.bed',
                                       'bim': '{root}.bim',
                                       'fam': '{root}.fam'}})

def read_plink_input_group(p, bfile_path):
    r'''Reads input group of PLINK files into pipeline'''
    p.read_input_group(**dict((suffix, f'{bfile_path}.{suffix}') for suffix in ['bed','bim','fam']))


def plink_clump():
backend = hp.BatchBackend(billing_project='ukb_diverse_pops')
p = hp.Pipeline(name='plink_clump', backend=backend,
                      default_image='gcr.io/ukbb-diversepops-neale/nbaya_plink:latest',
                      default_storage='500Mi', default_cpu=8)

pop = 'AFR'
    
## read ref LD PLINK files with samples from `pop` excluded
bfile_path = f'{ldprune_wd}/subsets/not_{pop}'
bfile = read_plink_input_group(p=p, bfile_path=bfile_path)

# TODO: Make loop over phenotypes

pheno='1717'
trait_type='continuous'
coding='1717'

## convert variant results from ht to tsv
t1 = p.new_task(name='to-tsv')
    
hail_script='/Users/nbaya/Documents/lab/diverse_pops/python/plink_clump_hail.py'
ss_ht_path = f'{bucket}/results/result/{pop}/{trait_type}-{pheno}-{coding}/variant_results.ht' # path to sumstats ht

t1.command(
f"""
PYTHONPATH=$PYTHONPATH:/ 
python3 {hail_script}
--input_file {ss_path}
--output_file {t1.ofile}
""".replace('\n', ' '))

## run plink clumping
t2 = p.new_task(name='clump')

t2.command(
f"""
plink \
--bfile {bfile} \
--clump {t2.ofile} \ 
--clump-field p.value.NA \
--clump-p1 1 \
--clump-p2 1 \
--clump-r2 0.1 \
--clump-kb 500 \
--out {t2.ofile}
""")

## gzip output
t3 = p.new_task(name='gzip')
t3.command(
f"""
column -t {t2.ofile} | gzip > {t3.ofile}
""")

pruned_ss_path = f'{ldprune_wd}/pruned/{pop}/refld_not{pop}/{trait_type}-{pheno}-{coding}.tsv.gz'
p.write_output(t3.ofile, pruned_ss_path)

p.run(open=True)

def local_plink_clump():
    pass
#    pop = 'AFR'
#    
#    backend = hp.LocalBackend(gsa_key_file='/Users/nbaya/.hail/ukb-diverse-pops.json')
#    p = hp.Pipeline(name='plink_clump', backend=backend)
#    
#    bfile_path = f'{ldprune_wd}/subsets/not_{pop}'
#    bfile = read_plink_input_group(p=p, bfile_path=bfile_path)
#    
#    
#    t1 = p.new_task(name='task1')
#    
#    trait_type='continuous'
#    pheno='1717'
#    coding='1717'
#    
#    ss_path = f'gs://ukb-diverse-pops/results/result/{pop}/{trait_type}-{pheno}-{coding}/variant_results.ht' # path to sumstats ht
#    hail_script='/Users/nbaya/Documents/lab/diverse_pops/python/plink_clump_hail.py'
#    
#    t1.command(
#    f"""
#    PYTHONPATH=$PYTHONPATH:/ PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=24g pyspark-shell"
#    python3 {hail_script}
#    --input_file {ss_path}
#    --output_file {t1.ofile}
#    """)
#    
#    t2 = p.new_task(name='task2')
#    
#    plink_local = '/Users/nbaya/Documents/lab/smiles/smiles/bash/plink/plink'
#    t2.command(f"""
#               {plink_local} \
#               --bfile {bfile} \
#               --clump --dummy 10 100 --make-bed --out {create.bfile1}
#               """)
    
#    p.write_output()



if __name__=="__main__":
    
    plink_clump
    