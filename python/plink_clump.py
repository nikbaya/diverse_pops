#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 09:35:41 2020

Clumping GWAS results with PLINK

@author: nbaya
"""

import hail as hl
import hailtop.batch as hb
import random 

bucket = 'gs://ukb-diverse-pops'
ldprune_dir = f'{bucket}/ld_prune'

all_pops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']


def declare_rg(t, name):
    t.declare_resource_group(**{name: {'bed': '{root}.bed',
                                       'bim': '{root}.bim',
                                       'fam': '{root}.fam'}})

def read_plink_input_group_chrom(p, method, subset, chrom):
    r'''
    Reads input group of PLINK files into Batch.
    NOTE: This format allows for use of a single fam (not specific to a chromosome)
    '''
    assert method in {'clump','sbayesr'}
    
    if method == 'clump':
        subsets_dir = f'{ldprune_dir}/subsets_5k'
        prefix = f'{subsets_dir}/{subset}/{subset}'
        return p.read_input_group(bed=f'{prefix}.chr{chrom}.bed',
                                  bim=f'{prefix}.chr{chrom}.bim',
                                  fam=f'{prefix}.fam')
    elif method == 'sbayesr':
        subsets_dir = f'{ldprune_dir}/subsets_50k'
        prefix = f'{subsets_dir}/{subset}/{subset}'
        return p.read_input_group(bed=f'{prefix}.hm3.chr{chrom}.bed',
                                  bim=f'{prefix}.hm3.chr{chrom}.bim',
                                  fam=f'{prefix}.fam')

def get_pheno_list(pheno_manifest, pop: str):
    r'''
    Returns list of phenotypes for population `pop`.
    '''
    
    pheno_manifest = pheno_manifest.filter(~pheno_manifest.pops.contains(pop)) # get 5 pop traits
    
    pheno_list = list(zip(pheno_manifest.trait_type.collect(),
                          pheno_manifest.phenocode.collect(), 
                          pheno_manifest.pheno_sex.collect(),
                          pheno_manifest.coding.collect(),
                          pheno_manifest.modifier.collect()))

    seed = 1
    random.seed(a=seed)
    random.shuffle(pheno_list)
#    
    n_traits = 2
    print(f'\nWARNING: For testing purposes, only using first {n_traits} traits (random seed: {seed})\n')
    pheno_list = pheno_list[:n_traits]
    
#    pheno_list = [
##            ('biomarkers', '30600', 'both_sexes', '', 'irnt'), # 6 pop LOO quant
##            ('categorical', '100240', 'both_sexes', '100240', ''), # 6 pop LOO binary
#            ('continuous', '1408', 'both_sexes', '', ''),# 5 pop not_AFR quant
#            ('icd10','Z37','both_sexes','',''), # 5 pop not_EUR binary
##            ('categorical', '1150', 'both_sexes', '3', ''), # 5 pop not_AFR binary
#            ]
    print(pheno_list)
    print(f'\nNumber of phenotypes for not_{pop}: {len(pheno_list)}')
    
    return pheno_list

def make_pheno_id(trait_type, phenocode, pheno_sex, coding, modifier):
    pheno_id = (trait_type+'-'+
                phenocode+'-'+
                pheno_sex+
                ('-'+coding if len(coding)>0 else '')+
                ('-'+modifier if len(modifier)>0 else '')).replace(' ','_').replace('/','_')
    return pheno_id    

def get_meta_sumstats(p, pop, pheno_id, method):
    assert method in {'clump','sbayesr'}
    filename=f'{pheno_id}.tsv.bgz'
    
    loo_dir ='gs://ukb-diverse-pops/ld_prune/loo/sumstats/batch1' 
    if hl.hadoop_is_file(f'{loo_dir}/{filename}'): # 6 pop LOO
        get_meta_ss = p.new_job(name=f'get_meta_ss_{pop}-{pheno_id}')
        get_meta_ss.storage('2G')
        get_meta_ss.cpu(1)
        get_meta_ss.command(' '.join(['set','-ex']))
        meta_ss = p.read_input(f'{loo_dir}/{filename}')
        if method=='clump':
            get_meta_ss.command(' '.join(['gunzip','-c',str(meta_ss),'|',
                                          'cut',f'-f1,{2+all_pops.index(pop)}','|',
                                          'sed',f"'s/pval_not_{pop}/P/g'",'|',
                                              'awk',"""'$2!="NA" {print}'""",
                                          '>',get_meta_ss.ofile]))
        elif method=='sbayesr':
            pass
    else:    
        assert trait_type in ['continuous','biomarkers','categorical','phecode', 'icd10', 'prescriptions']
        trait_category = 'quant' if trait_type in ['continuous','biomarkers'] else 'binary'
        five_pop_dir = f'gs://ukb-diverse-pops/sumstats_flat_files'
        tabix_dir = f'gs://ukb-diverse-pops/sumstats_flat_files_tabix'
        if hl.hadoop_is_file(f'{five_pop_dir}/{filename}'):
            get_meta_ss = p.new_job(name=f'get_meta_ss_{pop}-{pheno_id}')
            get_meta_ss = get_meta_ss.image('gcr.io/ukbb-diversepops-neale/nbaya_plink:latest')
            get_meta_ss.storage('2G')
            get_meta_ss.cpu(1)
            get_meta_ss.command(' '.join(['set','-ex']))
            meta_ss = p.read_input(f'{five_pop_dir}/{filename}')
            if method=='clump':
                pval_col_idx = 8 if trait_category == 'quant' else 9 # due to additional AF columns in binary traits, pvalue column location may change
                get_meta_ss.command(' '.join(['gunzip','-c',str(meta_ss),'|',
                                              'awk',f"""'{{print $1=$1":"$2":"$3":"$4, $2=${pval_col_idx}}}'""",'|',
                                              'sed','-e',"""'s/pval_meta/P/g'""",
                                              '-e',"""'s/chr:pos:ref:alt/SNP/g'""",'|',
                                              'awk',"""'$2!="NA" {print}'""",
                                              '>',get_meta_ss.ofile]))
            elif method=='sbayesr':
                assert False, "sbayesr is not an option"
        else:
            print(f'ERROR: No sumstats for {pheno_id} for not_{pop}')
            return None            
    return get_meta_ss.ofile

def get_sumstats_v2(p, pop, pheno_id, method, start_chrom=1, end_chrom=22):
    assert method in {'clump','sbayesr'}
    filename=f'{pheno_id}.tsv.bgz'
    trait_type = pheno_id.split('-')[0]
    trait_category = 'quant' if trait_type in ['continuous','biomarkers'] else 'binary'

    ss_dir = f'gs://ukb-diverse-pops/sumstats_flat_files'
    tabix_dir = f'gs://ukb-diverse-pops/sumstats_flat_files_tabix'
    meta_ss = p.read_input(f'{ss_dir}/{filename}')
    tabix = p.read_input(f'{tabix_dir}/{filename}.tbi')

    get_ss = p.new_job(name=f'get_ss_{pheno_id}')
    get_ss = get_ss.image('gcr.io/ukbb-diversepops-neale/nbaya_tabix:latest')
    get_ss.storage('1G')
    get_ss.cpu(1)
    get_ss.command(' '.join(['set','-ex']))
#    bgz_fname = f'{get_ss.ofile}.bgz'
#    tbi_fname = f'{get_ss.ofile}.bgz.tbi'
    bgz_fname = f'{get_ss.ofile}.bgz'
    tbi_fname = f'{get_ss.ofile}.bgz.tbi'
    get_ss.command(' '.join(['mv',meta_ss,bgz_fname]))
    get_ss.command(' '.join(['mv',tabix,tbi_fname]))
    pval_col_idx = 8 if trait_category == 'quant' else 9 # due to additional AF columns in binary traits, pvalue column location may change
    get_ss.command('\n'.join(
            f'''
            tabix -h {bgz_fname} {chrom} | \\
            awk '{{print $1=$1":"$2":"$3":"$4, $2=${pval_col_idx}}}' | \\
            sed -e 's/pval_meta/P/g' \\
                -e 's/chr:pos:ref:alt/SNP/g' | \\
            awk '$2!="NA" {{print}}' > {get_ss[f"ofile_{chrom}"]}
            '''
            for chrom in range(start_chrom, end_chrom+1)
            )
    )
    ss_dict = {
            chrom:get_ss[f'ofile_{chrom}']
            for chrom in range(start_chrom,end_chrom+1)
    }
    return ss_dict
    
def get_ss_chrom(p, pheno_id, chrom):
    filename = f'{pheno_id}.tsv.bgz'
    trait_type = pheno_id.split('-')[0]
    trait_category = 'quant' if trait_type in ['continuous','biomarkers'] else 'binary'

    ss_dir = f'gs://ukb-diverse-pops/sumstats_flat_files'
    tabix_dir = f'gs://ukb-diverse-pops/sumstats_flat_files_tabix'
    meta_ss = p.read_input(f'{ss_dir}/{filename}')
    tabix = p.read_input(f'{tabix_dir}/{filename}.tbi')

    get_ss_chrom = p.new_job(name=f'get_ss_{pheno_id}_chr{chrom}')
    get_ss_chrom = get_ss_chrom.image('gcr.io/ukbb-diversepops-neale/nbaya_tabix:latest')
    get_ss_chrom.storage('2G')
    get_ss_chrom.cpu(1)
    get_ss_chrom.command(' '.join(['set','-ex']))
    pval_col_idx = 8 if trait_category == 'quant' else 9 # due to additional AF columns in binary traits, pvalue column location may change
    bgz_fname = f'{get_ss_chrom.ofile}.tsv.bgz'
    tbi_fname = f'{get_ss_chrom.ofile}.tsv.bgz.tbi'
    get_ss_chrom.command(' '.join(['mv',meta_ss,bgz_fname]))
    get_ss_chrom.command(' '.join(['mv',tabix,tbi_fname]))
    get_ss_chrom.command(' '.join(['tabix','-h',bgz_fname,str(chrom),'|',
                                   'awk',f"""'{{print $1=$1":"$2":"$3":"$4, $2=${pval_col_idx}}}'""",'|',
                                   'sed','-e',"""'s/pval_meta/P/g'""",
                                   '-e',"""'s/chr:pos:ref:alt/SNP/g'""",'|',
                                   'awk',"""'$2!="NA" {print}'""",
                                   '>',get_ss_chrom.ofile]))
    return get_ss_chrom.ofile
    
def get_adj_betas(p, pop, pheno_id, hail_script):
    r'''
    wrapper method for both PLINK clumping and SBayesR
    '''
                
#    task_suffix = f'not_{pop}-{trait_type}-{phenocode}-{coding}'
#    attributes_dict = {'pop': pop,
#                       'trait_type': trait_type,
#                       'pheno': phenocode,
#                       'coding': coding}    
#    n_threads = 8
    
#    output_dir = f'{ldprune_dir}/results/not_{pop}/{trait_type}-{pheno}-{coding}'
    temp_bucket = 'gs://ukbb-diverse-temp-30day'
    output_dir = f'{temp_bucket}/results-test/not_{pop}/{pheno_id}' # for testing

    clump_output_txt = f'{output_dir}/clumped_results-test.txt' # PLINK clump output txt file
    clump_output_ht = f'{output_dir}/clump_results-test.ht' # PLINK clump output hail table
    
    sbayesr_output_txt = f'{output_dir}/sbayesr_results-test.txt' # SBayesR output txt file
    sbayesr_output_ht = f'{output_dir}/sbayesr_results-test.ht' # SBayesR output hail table
    
    test_run = True
    
    if test_run:
        print('\n\nWARNING: Test run may overwrite existing results!\n')
    
    if not hl.hadoop_is_file(f'{clump_output_ht}/_SUCCESS') or test_run:
#        ss = get_meta_sumstats(p=p,
#                               hail_script=hail_script,
#                               pop=pop, 
#                               trait_type=trait_type, 
#                               phenocode=phenocode, 
#                               pheno_sex=pheno_sex, 
#                               coding=coding,
#                               modifier=modifier,
#                               method='clump')
        ss_dict = get_sumstats_v2(p=p,
                                  pop=pop, 
                                  pheno_id=pheno_id,
                                  method='clump')
        
        if ss_dict!=None:
            
            run_method(p=p, 
                       pop=pop, 
                       pheno_id=pheno_id,  
                       hail_script=hail_script, 
                       output_txt=clump_output_txt,
                       output_ht=clump_output_ht,
                       ss_dict=ss_dict,
                       method='clump')
        
#    if not hl.hadoop_is_file(f'{sbayesr_output_ht}/_SUCCESS'):
#            
#        meta_ss = p.read_input(f'{ldprune_dir}/loo/not_AFR/batch2/biomarkers-30740-30740.tsv.bgz')
#
#        run_method(p=p, 
#                   pop=pop, 
#                   pheno=pheno, 
#                   coding=coding, 
#                   trait_type=trait_type, 
#                   hail_script=hail_script, 
#                   output_txt=sbayesr_output_txt,
#                   output_ht=sbayesr_output_ht,
#                   ss=meta_ss,
#                   method='sbayesr')
        
        
def run_method(p, pop, pheno_id, hail_script, output_txt, output_ht, ss_dict, method):
    r'''
    Runs either PLINK clump (method = 'clump') or SBayesR (method = 'sbayesr')
    '''
    assert method in {'clump','sbayesr'}
    
    task_suffix = f'not_{pop}-{pheno_id}'
    # TODO: if method = 'sbayesr' check if LD matrix has already been calculated
    
    tasks = []
        
    ## run plink clumping
    for chrom, ss_chrom in ss_dict.items():
        ## read ref ld plink files 
        bfile = read_plink_input_group_chrom(p=p, 
                                             method=method,
                                             subset=f'not_{pop}',
                                             chrom=chrom)
        
#        ss_chrom = get_ss_chrom(p=p, 
#                          pheno_id=pheno_id, 
#                          chrom=chrom)
        
        get_betas = p.new_job(name=f'{method}_{task_suffix}_chr{chrom}')
        
        # TODO: change image to include GCTB if running SBayesR?
        get_betas.storage('5G')
        get_betas.cpu(1) # plink clump cannot multithread
        
        get_betas.command(' '.join(['set','-ex']))
        
        if method == 'clump':
#            clump_memory = -15*(chrom-1)+400 # Memory requested for PLINK clumping in MB. equation: -15*(chrom-1) + 500 is based on 400 MB for chr 1, 80 MB for chr 22
            clump_memory = 3.75 # in GB
            get_betas.memory(clump_memory) # default: 30G
#            get_betas.command(f' gunzip -c {ss} | grep "^{chrom}:\|SNP" |'+"awk '{ print $1,$7 }' | head")
#            get_betas.command(f'head {bfile.bim}')
#            get_betas.command(' '.join([f'gunzip -c {ss} | grep "^{chrom}:\|SNP" > {get_betas.ofile}_ss']))
            
            get_betas.command(' '.join(['plink',
                                        '--bfile', str(bfile),
                                        '--memory',str(clump_memory*1000), # memory in MB
                                        '--threads','1', # explicitly set threads to 1
                                        '--clump', ss_chrom,
#                                        '--clump', '<(','grep',f'"^{chrom}:\|SNP"', str(ss), ')',
                                        '--clump-field P',
                                        '--clump-snp-field SNP',
                                        '--clump-p1 1',
                                        '--clump-p2 1',
                                        '--clump-r2 0.1',
                                        '--clump-kb 500',
                                        '--chr', str(chrom),
                                        '--out',f'{get_betas.ofile}_tmp']))
            get_betas.command(' '.join(['awk',"'{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12 }'",
                                        "OFS='\t'",f'{get_betas.ofile}_tmp.clumped','2>','/dev/null',
                                 '|','tail -n+2','>',str(get_betas.ofile)]))
        elif method == 'sbayesr':
            ldm_type = 'full' # options: full, sparse
            ldm_path = f'{ldprune_dir}/subsets_50k/not_{pop}/ldm/not_{pop}.hm3.chr{chrom}.maf_gt_0.ldm.{ldm_type}'
#            ldm_path = f'{ldprune_dir}/subsets_50k/not_{pop}/ldm/not_{pop}.hm3.chr{chrom}.maf_gt_0.chisq_5.ldm.sparse'
            
            if hl.hadoop_is_file(f'{ldm_path}.info') and hl.hadoop_is_file(f'{ldm_path}.bin'):
                ldm = p.read_input_group(info=f'{ldm_path}.info',
                                         bin=f'{ldm_path}.bin')
            else:
                make_ldm = p.new_job(name=f'make_{ldm_type}_ldm_{task_suffix}.chr{chrom}')
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
                
            get_betas.declare_resource_group(out={'log':'{root}.log',
                                                  'snpRes':'{root}.snpRes',
                                                  'parRes':'{root}.parRes',
                                                  'mcmcsamples.SnpEffects':'{root}.mcmcsamples.SnpEffects',
                                                  'mcmcsamples.Par':'{root}.mcmcsamples.Par'})
            get_betas.command(' '.join(['wget',
                                        'https://cnsgenomics.com/software/gctb/download/gctb_2.0_Linux.zip',
                                        '-P', '~/']))
            get_betas.memory('18G')
            get_betas.command(' '.join(['unzip','~/gctb_2.0_Linux.zip','-d','~/']))
            get_betas.command(' '.join(['ls','-ltrR','~/']))
            get_betas.command(' '.join(['mv','~/gctb_2.0_Linux/gctb','/usr/local/bin/']))
            get_betas.command(' '.join(['gctb',
                                        '--sbayes R', 
                                        '--ldm', str(ldm),
                                        '--pi 0.95,0.02,0.02,0.01',
                                        '--gamma 0.0,0.01,0.1,1',
                                        '--gwas-summary', f' <( gunzip -c {ss} | grep -v "NA" )',
                                        '--chain-length 10000',
                                        '--burn-in 2000',
                                        '--out-freq 10',
                                        '--out',f'{get_betas.out}']))
            get_betas.command(' '.join(['head',f'{get_betas.out}.snpRes']))
            get_betas.command(' '.join(['mv',f'{get_betas.out}.snpRes', str(get_betas.ofile)]))
            
        tasks.append(get_betas)

    get_betas_sink = p.new_job(name=f'{method}_sink_{task_suffix}')
    get_betas_sink.command(f'cat {" ".join([t.ofile for t in tasks])} > {get_betas_sink.ofile}') # this task implicitly depends on the chromosome scatter tasks
    p.write_output(get_betas_sink.ofile, output_txt)
    
    ## import as hail table and save
    n_threads = 8
    tsv_to_ht = p.new_job(name=f'{method}_to_ht_{task_suffix}')
    tsv_to_ht = tsv_to_ht.image('gcr.io/ukbb-diversepops-neale/hail_utils:3.4')
    tsv_to_ht.storage('1G')
    tsv_to_ht.memory('100M')
    tsv_to_ht.cpu(n_threads)
    tsv_to_ht.depends_on(get_betas_sink)
    
    tsv_to_ht.command(' '.join(['set','-ex']))
    tsv_to_ht.command(' '.join(['PYTHONPATH=$PYTHONPATH:/',
                        'PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=4g --conf spark.executor.memory=24g pyspark-shell"']))
    tsv_to_ht.command(' '.join(['python3', str(hail_script),
                         '--input_file', f'"{output_txt}"', # output_txt must be doubly enclosed by quotes needed for files with "|" in their pheno_id
                         '--tsv_to_ht',
                         '--pop', pop,
                         '--output_file', f'"{output_ht}"', # output_ht must be doubly enclosed by quotes needed for files with "|" in their pheno_id
                         '--n_threads', str(n_threads),
                         '--overwrite']))    

def main():    
    backend = hb.ServiceBackend(billing_project='ukb_diverse_pops',
                                   bucket='ukbb-diverse-temp-30day/nb-batch-tmp')
#    backend = batch.LocalBackend(tmp_dir='/tmp/batch/')
    
    p = hb.batch.Batch(name='test-clump', backend=backend,
                          default_image='gcr.io/ukbb-diversepops-neale/nbaya_plink:0.1',
                          default_storage='500Mi', default_cpu=8)
    
    ## download hail script to VM
    hail_script = p.read_input(f'{ldprune_dir}/scripts/python/plink_clump_hail.py')
    
    ## get phenotype manifest
    pheno_manifest = hl.import_table(f'{bucket}/combined_results/phenotype_manifest.tsv.bgz',
                                     impute=True)
    pheno_manifest = pheno_manifest.filter(pheno_manifest.num_pops>=5) # necessary for LOO clumping
    
    for pop in all_pops:   
        pheno_list = get_pheno_list(pheno_manifest=pheno_manifest,
                                    pop=pop)
    
        for trait_type, phenocode, pheno_sex, coding, modifier in pheno_list:
            pheno_id = make_pheno_id(trait_type, phenocode, pheno_sex, coding, modifier)
            get_adj_betas(p=p, 
                          pop=pop, 
                          pheno_id=pheno_id, 
                          hail_script=hail_script)
    p.run(open=True)
#    if type(backend)==batch.ServiceBackend:
#        print('running')
#    else:
#        p.run(verbose=True,
#              delete_scratch_on_exit=True)
    
    backend.close()


if __name__=="__main__":
    
    hl.init(default_reference='GRCh38')
    main()
    
    