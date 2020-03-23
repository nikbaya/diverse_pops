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

def read_plink_input_group_chrom(p, subset, chrom):
    r'''
    Reads input group of PLINK files into pipeline.
    NOTE: This format allows for use of a single fam (not specific to a chromosome)
    '''
    prefix = f'{ldprune_wd}/subsets/{subset}/{subset}'
    return p.read_input_group(bed=f'{prefix}.chr{chrom}.bed',
                              bim=f'{prefix}.chr{chrom}.bim',
                              fam=f'{prefix}.fam')

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
    
    pheno_list = [('100001', 'irnt', 'continuous')]#,
#                  ('100002', 'irnt', 'continuous')] # dummy list for testing
    
    return pheno_list


def run_clumping(p, pop, pheno, coding, trait_type, hail_script):
    task_suffix = f'not_{pop}-{trait_type}-{pheno}-{coding}'
    attributes_dict = {'pop': pop,
                       'trait_type': trait_type,
                       'pheno': pheno,
                       'coding': coding}    
    n_threads = 8
    
    output_dir = f'{ldprune_wd}/results/not_{pop}/{trait_type}-{pheno}-{coding}'
    output_txt = f'{output_dir}/clumped_results-test.txt'

    if not hl.hadoop_is_file(output_txt):
        t1 = p.new_task(name=f'get_meta_ss_{task_suffix}')
        t1.attributes.update(attributes_dict)
        t1 = t1.image('gcr.io/ukbb-diversepops-neale/hail_utils:3.4')
        t1.storage('1G') # prev: 8
        t1.memory('100M') # prev: 4
        t1.cpu(n_threads)
            
    #    meta_ht_path = f'{ldprune_wd}/meta_analysis.all_pops_old.ht' # path to old meta-analyzed sumstats ht
        meta_ht_path = f'{ldprune_wd}/meta_analysis.all_pops.ht' # path to meta-analyzed sumstats ht
        
        t1.command(
        f"""
        PYTHONPATH=$PYTHONPATH:/ PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=4g --conf spark.executor.memory=24g pyspark-shell"
        python3 {hail_script}
        --input_file {meta_ht_path}
        --pop {pop}
        --pheno {pheno}
        --coding {coding}
        --trait_type {trait_type}
        --output_file {t1.ofile}
        --get_meta_sumstats
        --n_threads {n_threads}
        """.replace('\n', ' '))
        
        clumping_tasks = []
        
        ## run plink clumping
        for chrom in range(1,23):
            ## read ref ld plink files 
            bfile = read_plink_input_group_chrom(p=p, 
                                                 subset=f'not_{pop}',
                                                 chrom=chrom)
            
            t2 = p.new_task(name=f'plink_clump_{task_suffix}.chr{chrom}')
            t2.attributes.update(attributes_dict)
            t2.depends_on(t1)
            t2.storage('5G')
            t2.memory('1G')
            
            t2.command(' '.join(['set','-ex']))
            t2.command(' '.join(['plink',
                                 '--bfile', str(bfile),
                                 '--clump', str(t1.ofile),
                                 '--clump-field P',
                                 '--clump-snp-field varid',
                                 '--clump-p1 1',
                                 '--clump-p2 1',
                                 '--clump-r2 0.1',
                                 '--clump-kb 500',
                                 '--chr', str(chrom),
                                 '--out',f'{t2.ofile}_tmp']))
            t2.command(' '.join(['awk',"'{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12 }'","OFS='\t'",f'{t2.ofile}_tmp.clumped',
                                 '|','tail -n+2','>',str(t2.ofile)]))
            clumping_tasks.append(t2)
    
        t3 = p.new_task(name=f'plink_clump_sink_{task_suffix}')
        t3.depends_on(*clumping_tasks)
        t3.command(f'cat {" ".join([t.ofile for t in clumping_tasks])} > {t3.ofile}')
        p.write_output(t3.ofile, output_txt)
    
    ## import as hail table and save
    t4 = p.new_task(name=f'tsv_to_ht_{task_suffix}')
    t4.attributes.update(attributes_dict)
    t4 = t4.image('gcr.io/ukbb-diversepops-neale/hail_utils:3.4')
    t4.storage('1G')
    t4.memory('100M')
    t4.cpu(n_threads)
    if not hl.hadoop_is_file(output_txt):
        t4.depends_on(t3)
    
    output_ht = f'{output_dir}/clump_results.ht'
    
    t4.command(' '.join(['set','-ex']))
    t4.command(' '.join(['PYTHONPATH=$PYTHONPATH:/',
                        'PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=4g --conf spark.executor.memory=24g pyspark-shell"']))
    t4.command(' '.join(['python3', str(hail_script),
                         '--input_file', output_txt,
                         '--pop', pop,
                         '--pheno', pheno,
                         '--coding', coding,
                         '--trait_type', trait_type,
                         '--output_file', output_ht,
                         '--tsv_to_ht',
                         '--n_threads', str(n_threads),
                         '--overwrite']))    

def main():
    pops = ['EUR'] #'AFR']#, 
    
    backend = hp.BatchBackend(billing_project='ukb_diverse_pops')
    p = hp.Pipeline(name='test-clump', backend=backend,
                    default_image='gcr.io/ukbb-diversepops-neale/nbaya_plink:0.1',
                    default_storage='500Mi', default_cpu=8)
    
    ## download hail script to VM
    hail_script = p.read_input(f'{ldprune_wd}/scripts/python/plink_clump_hail.py')
    
    for pop in pops:   
        pheno_list = get_pheno_list(pop)
    
    
        for pheno, coding, trait_type in pheno_list:
            run_clumping(p=p, 
                         pop=pop, 
                         pheno=pheno, 
                         coding=coding, 
                         trait_type=trait_type,
                         hail_script=hail_script)
    
    p.run(open=True)


if __name__=="__main__":
    
    hl.init(default_reference='GRCh38')
    main()
    