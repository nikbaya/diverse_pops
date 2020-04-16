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
subsets_dir = f'{ldprune_wd}/subsets_50k'


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
    prefix = f'{subsets_dir}/{subset}/{subset}'
    return p.read_input_group(bed=f'{prefix}.hm3.chr{chrom}.bed',
                              bim=f'{prefix}.hm3.chr{chrom}.bim',
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
    
    pheno_list = [('100002', 'irnt', 'continuous')]#,
#                  ('100001', 'irnt', 'continuous')] # dummy list for testing
    
    return pheno_list
    

def run_method(p, pop, pheno, coding, trait_type, hail_script, output_txt,
                output_ht, get_meta_ss, method):
    r'''
    Runs either PLINK clump (method = 'clump') or SBayesR (method = 'sbayesr')
    '''
    task_suffix = f'not_{pop}-{trait_type}-{pheno}-{coding}'
    # TODO: if method = 'sbayesr' check if LD matrix has already been calculated
    
    n_threads = 8
    
    tasks = []
        
    ## run plink clumping
    for chrom in range(22,23):
        ## read ref ld plink files 
        bfile = read_plink_input_group_chrom(p=p, 
                                             subset=f'not_{pop}',
                                             chrom=chrom)
        
        get_betas = p.new_task(name=f'{method}_{task_suffix}.chr{chrom}')
        # TODO: change image to include GCTB if running SBayesR?
        get_betas.depends_on(get_meta_ss)
        get_betas.storage('5G')
        get_betas.cpu(1) # plink clump cannot multithread
        
        get_betas.command(' '.join(['set','-ex']))
        if method == 'clump':
            get_betas.memory('18G')
            get_betas.command(' '.join(['plink',
                                        '--bfile', str(bfile),
                                        '--clump', str(get_meta_ss.ofile),
                                        '--clump-field P',
                                        '--clump-snp-field varid',
                                        '--clump-p1 1',
                                        '--clump-p2 1',
                                        '--clump-r2 0.1',
                                        '--clump-kb 500',
                                        '--chr', str(chrom),
                                        '--out',f'{get_betas.ofile}_tmp']))
            get_betas.command(' '.join(['awk',"'{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12 }'","OFS='\t'",f'{get_betas.ofile}_tmp.clumped',
                                 '|','tail -n+2','>',str(get_betas.ofile)]))
        elif method == 'sbayesr':
            ldm_type = 'sparse' # options: full, sparse
            ldm_path = f'{subsets_dir}/not_{pop}/ldm/chr{chrom}.ldm.{ldm_type}'
            if hl.hadoop_is_file(ldm_path):
                ldm = p.read_input(ldm_path)
            else:
                make_ldm = p.new_task(name=f'make_{ldm_type}_ldm_{task_suffix}.chr{chrom}')
                make_ldm.memory('60G')
                make_ldm.command(' '.join(['wget',
                                        'https://cnsgenomics.com/software/gctb/download/gctb_2.0_Linux.zip',
                                        '-P', '~/']))
                make_ldm.command(' '.join(['unzip','~/gctb_2.0_Linux.zip','-d','~/']))
                make_ldm.command(' '.join(['ls','-ltrR','~/']))
                make_ldm.command(' '.join(['mv','~/gctb_2.0_Linux/gctb','/usr/local/bin/']))
                make_ldm.command(' '.join(['plink',
                                           '--bfile', str(bfile),
                                           '--maf 0.0000000001',
                                           '--make-bed',
                                           '--out', f'{make_ldm.ofile}_tmp1']))
                make_ldm.command(' '.join(['gctb',
                                            '--bfile', f'{make_ldm.ofile}_tmp1',
#                                            '--snp 1-1000',
                                            f'--make-{ldm_type}-ldm', 
                                            '--out',f'{make_ldm.ofile}_tmp2']))
                # TODO: use both .bin and .info files
                make_ldm.command(' '.join(['mv',f'{make_ldm.ofile}_tmp2.ldm.{ldm_type}', str(make_ldm.ofile)]))
                p.write_output(make_ldm.ofile, ldm_path)
                ldm = make_ldm.ofile
            get_betas.command(' '.join(['wget',
                                        'https://cnsgenomics.com/software/gctb/download/gctb_2.0_Linux.zip',
                                        '-P', '~/']))
            get_betas.memory('18G')
            get_betas.command(' '.join(['unzip','~/gctb_2.0_Linux.zip','-d','~/']))
            get_betas.command(' '.join(['ls','-ltrR','~/']))
            get_betas.command(' '.join(['mv','~/gctb_2.0_Linux/gctb','/usr/local/bin/']))
            get_betas.command(' '.join(['gctb',
                                        '--sbayes R', 
                                        '--ldm', f'{ldm}',
                                        '--pi 0.95,0.02,0.02,0.01',
                                        '--gamma 0.0,0.01,0.1,1',
                                        '--gwas-summary', str(get_meta_ss.ofile),
                                        '--chain-length 10000',
                                        '--burn-in 2000',
                                        '--out-freq 10',
                                        '--out',f'{get_betas.ofile}_tmp']))
            get_betas.command(' '.join(['head',f'{get_betas.ofile}_tmp.snpRes']))
            get_betas.command(' '.join(['mv',f'{get_betas.ofile}_tmp.snpRes', str(get_betas.ofile)]))
            
        tasks.append(get_betas)

    get_betas_sink = p.new_task(name=f'{method}_sink_{task_suffix}')
    get_betas_sink.command(f'cat {" ".join([t.ofile for t in tasks])} > {get_betas_sink.ofile}') # this task implicitly depends on the chromosome scatter tasks
    p.write_output(get_betas_sink.ofile, output_txt)
    
    ## import as hail table and save
    tsv_to_ht = p.new_task(name=f'{method}_to_ht_{task_suffix}')
    tsv_to_ht = tsv_to_ht.image('gcr.io/ukbb-diversepops-neale/hail_utils:3.4')
    tsv_to_ht.storage('1G')
    tsv_to_ht.memory('100M')
    tsv_to_ht.cpu(n_threads)
    tsv_to_ht.depends_on(get_betas_sink)
    
    tsv_to_ht.command(' '.join(['set','-ex']))
    tsv_to_ht.command(' '.join(['PYTHONPATH=$PYTHONPATH:/',
                        'PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=4g --conf spark.executor.memory=24g pyspark-shell"']))
    tsv_to_ht.command(' '.join(['python3', str(hail_script),
                         '--input_file', output_txt,
                         '--pop', pop,
                         '--pheno', pheno,
                         '--coding', coding,
                         '--trait_type', trait_type,
                         '--output_file', output_ht,
                         '--tsv_to_ht',
                         '--n_threads', str(n_threads),
                         '--overwrite']))    

def get_adj_betas(p, pop, pheno, coding, trait_type, hail_script):
    r'''
    wrapper method for both PLINK clumping and SBayesR
    '''
    
    task_suffix = f'not_{pop}-{trait_type}-{pheno}-{coding}'
    attributes_dict = {'pop': pop,
                       'trait_type': trait_type,
                       'pheno': pheno,
                       'coding': coding}    
    n_threads = 8
    
#    output_dir = f'{ldprune_wd}/results/not_{pop}/{trait_type}-{pheno}-{coding}'
    output_dir = f'{ldprune_wd}/test-results/not_{pop}/{trait_type}-{pheno}-{coding}' # for testing

#    clump_output_txt = f'{output_dir}/clumped_results-test.txt' # PLINK clump output txt file
    clump_output_ht = f'{output_dir}/clump_results-test.ht' # PLINK clump output hail table
    
    sbayesr_output_txt = f'{output_dir}/sbayesr_results-test.txt' # SBayesR output txt file
    sbayesr_output_ht = f'{output_dir}/sbayesr_results-test.ht' # SBayesR output hail table
    
    if not hl.hadoop_is_file(f'{clump_output_ht}/_SUCCESS') or not hl.hadoop_is_file(f'{sbayesr_output_ht}/_SUCCESS'):
                
        get_meta_ss = p.new_task(name=f'get_meta_ss_{task_suffix}')
        get_meta_ss.attributes.update(attributes_dict)
        get_meta_ss = get_meta_ss.image('gcr.io/ukbb-diversepops-neale/hail_utils:3.4')
        get_meta_ss.storage('5G') # prev: 8
        get_meta_ss.memory('100M') # prev: 4
        get_meta_ss.cpu(n_threads)
            
        meta_ht_path = f'{ldprune_wd}/meta_analysis.all_pops.ht' # path to meta-analyzed sumstats ht
                
        get_meta_ss.command(' '.join(['PYTHONPATH=$PYTHONPATH:/',
                        'PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=4g --conf spark.executor.memory=24g pyspark-shell"']))
        get_meta_ss.command(' '.join(['python3', str(hail_script),
                             '--input_file', meta_ht_path,
                             '--pop', pop,
                             '--pheno', pheno,
                             '--coding', coding,
                             '--trait_type', trait_type,
                             '--output_file', str(get_meta_ss.ofile),
#                             '--get_meta_sumstats',                      
                             '--test_get_meta_sumstats', # NOTE: Only for testing
                             '--n_threads', str(n_threads)]))
    
#    if not hl.hadoop_is_file(f'{clump_output_ht}/_SUCCESS'):
#        run_method(p=p, 
#                   pop=pop, 
#                   pheno=pheno, 
#                   coding=coding, 
#                   trait_type=trait_type, 
#                   hail_script=hail_script, 
#                   output_txt=clump_output_txt,
#                   output_ht=clump_output_ht,
#                   get_meta_ss=get_meta_ss,
#                   method='clump')
    if not hl.hadoop_is_file(f'{sbayesr_output_ht}/_SUCCESS'):
        run_method(p=p, 
                   pop=pop, 
                   pheno=pheno, 
                   coding=coding, 
                   trait_type=trait_type, 
                   hail_script=hail_script, 
                   output_txt=sbayesr_output_txt,
                   output_ht=sbayesr_output_ht,
                   get_meta_ss=get_meta_ss,
                   method='sbayesr')
        
        

def main():
    pops = ['EUR'] # 'AFR']#
    
    backend = hp.BatchBackend(billing_project='ukb_diverse_pops')
    p = hp.Pipeline(name='test-clump', backend=backend,
                    default_image='gcr.io/ukbb-diversepops-neale/nbaya_plink:0.1',
                    default_storage='500Mi', default_cpu=8)
    
    ## download hail script to VM
    hail_script = p.read_input(f'{ldprune_wd}/scripts/python/plink_clump_hail.py')
    
    for pop in pops:   
        pheno_list = get_pheno_list(pop)
    
    
        for pheno, coding, trait_type in pheno_list:
            get_adj_betas(p=p, 
                          pop=pop, 
                          pheno=pheno, 
                          coding=coding, 
                          trait_type=trait_type,
                          hail_script=hail_script)
    
    p.run(open=True)
    
    backend.close()


if __name__=="__main__":
    
    hl.init(default_reference='GRCh38')
    main()
    
    