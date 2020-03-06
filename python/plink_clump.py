#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 09:35:41 2020

Clumping GWAS results with PLINK

@author: nbaya
"""

import hail as hl
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
    return p.read_input_group(**dict((suffix, f'{bfile_path}.{suffix}') for suffix in ['bed','bim','fam']))

def get_pheno_list(pop: str):
    r'''
    Returns list of phenotypes for population `pop`.
    '''
#    pheno_list_path = f'{ldprune_wd}/pheno_lists/pheno_list_{pop}.ht'
#    pheno_list_path = f'{ldprune_wd}/pheno_lists/pheno_list_intersection.ht' # intersection of 921 traits with results for all pops
#
#    if hl.hadoop_is_file(f'{pheno_list_path}/_SUCCESS'):
#        pheno_list_ht = hl.read_table(pheno_list_path) # phenotypes from population `pop`
#    else:
#        raise Exception(f'pheno list ht does not exist for pop {pop}')
#    
#    assert list(pheno_list_ht.key)==['pheno','coding','trait_type'], f'key ordering is incorrect for pop {pop}'
#    
#    pheno_list = list(zip(*[pheno_list_ht[f].collect() for f in pheno_list_ht.key])) # list of tuples
    
    pheno_list = [('100001', 'irnt', 'continuous'),
                  ('100002', 'irnt', 'continuous')] # dummy list for testing
    
    return pheno_list


def run_clumping(p, pop, bfile, pheno, coding, trait_type, hail_script):
    
    task_suffix = f'{pop}-{pheno}-{coding}-{trait_type}'
    
    t1 = p.new_task(name=f'get_meta_ss_{task_suffix}')
    t1 = t1.image('gcr.io/ukbb-diversepops-neale/hail_utils:3.2')
    t1.storage('4G')
    t1.memory('8G')
        
    meta_mt_path = f'{bucket}/combined_results/meta_analysis.mt' # path to sumstats ht
    
    t1.command(
    f"""
    PYTHONPATH=$PYTHONPATH:/ 
    python3 {hail_script}
    --input_file {meta_mt_path}
    --pop {pop}
    --pheno {pheno}
    --coding {coding}
    --trait_type {trait_type}
    --output_file {t1.ofile}
    --get_meta_sumstats
    """.replace('\n', ' '))
    
    
    ## run plink clumping
    t2 = p.new_task(name=f'plink_clump_{task_suffix}')
    t2.depends_on(t1)
    t2.storage('40G').memory('2G')
    
    # TODO: change default storage & memory?
    
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
    
#    ## gzip output
#    t3 = p.new_task(name=f'gzip_{task_suffix}')
#    t3.depends_on(t2)
#    
#    t3.command(
#    f"""
#    column -t {t2.ofile} | gzip > {t3.ofile}
#    """)
#    
#    pruned_ss_path = f'{ldprune_wd}/pruned/{pop}/refld_not{pop}/{trait_type}-{pheno}-{coding}.tsv.gz'
#    p.write_output(t3.ofile, pruned_ss_path)
    
    ## to ht
    

    

def main():
hl.init(default_reference='GRCh38')

pops = ['AFR']

backend = hp.BatchBackend(billing_project='ukb_diverse_pops')
p = hp.Pipeline(name='test-clump', backend=backend,
                      default_image='gcr.io/ukbb-diversepops-neale/nbaya_plink:latest',
                      default_storage='500Mi', default_cpu=8)

## download hail script to VM
hail_script = p.read_input(f'{ldprune_wd}/scripts/python/plink_clump_hail.py')

for pop in pops:   
    pheno_list = get_pheno_list(pop)
    
    ## read ref ld plink files 
#bfile_path = f'{ldprune_wd}/subsets/not_{pop}' # samples from `pop` are excluded from plink files
bfile_path = f'{ldprune_wd}/subsets/not_{pop}.chr22' # dummy bfile for testing
bfile = read_plink_input_group(p=p, bfile_path=bfile_path)

for pheno, coding, trait_type in pheno_list:
pheno, coding, trait_type = pheno_list[0]
run_clumping(p=p, 
             pop=pop, 
             bfile=bfile,
             pheno=pheno, 
             coding=coding, 
             trait_type=trait_type,
             hail_script=hail_script)

p.run(open=True)

    


# TODO: Make loop over phenotypes


## convert variant results from ht to tsv

p.run(open=True)

def local_main():

    pop = 'AFR'
    
    backend = hp.LocalBackend(gsa_key_file='/Users/nbaya/.hail/ukb-diverse-pops.json')
    p = hp.Pipeline(name='plink_clump', backend=backend)
    
    bfile_path = f'{ldprune_wd}/subsets/not_{pop}'
    bfile = read_plink_input_group(p=p, bfile_path=bfile_path)
    
    
    t1 = p.new_task(name='task1')
    
    trait_type='continuous'
    pheno='1717'
    coding='1717'
    
    ss_path = f'gs://ukb-diverse-pops/results/result/{pop}/{trait_type}-{pheno}-{coding}/variant_results.ht' # path to sumstats ht
    hail_script='/Users/nbaya/Documents/lab/diverse_pops/python/plink_clump_hail.py'
    
    t1.command(
    f"""
    PYTHONPATH=$PYTHONPATH:/ PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=24g pyspark-shell"
    python3 {hail_script}
    --input_file {ss_path}
    --output_file {t1.ofile}
    """)
    
    t2 = p.new_task(name='task2')
    
    plink_local = '/Users/nbaya/Documents/lab/smiles/smiles/bash/plink/plink'
    t2.command(f"""
               {plink_local} \
               --bfile {bfile} \
               --clump --dummy 10 100 --make-bed --out {create.bfile1}
               """)
    
    p.write_output()



if __name__=="__main__":
    
    main()
    