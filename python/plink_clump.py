#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 09:35:41 2020

Clumping GWAS results with PLINK

@author: nbaya
"""

import hail as hl
from hailtop import batch


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

def get_pheno_list(pop: str):
    r'''
    Returns list of phenotypes for population `pop`.
    '''
##    pheno_list_path = f'{ldprune_dir}/pheno_lists/pheno_list_{pop}.ht'
##    pheno_list_path = f'{ldprune_dir}/pheno_lists/pheno_list_intersection.ht' # intersection of 921 traits with results for all pops
##
##    if hl.hadoop_is_file(f'{pheno_list_path}/_SUCCESS'):
##        pheno_list_ht = hl.read_table(pheno_list_path) # phenotypes from population `pop`
##    else:
##        raise Exception(f'pheno list ht does not exist for pop {pop}')
##    
##    assert list(pheno_list_ht.key)==['pheno','coding','trait_type'], f'key ordering is incorrect for pop {pop}'
##    
#a
##    pheno_list = list(zip(*[pheno_list_ht[f].collect() for f in pheno_list_ht.key])) # list of tuples
##    
##    pheno_list = [('100420','100420','continuous')]
##    [('100002', 'irnt', 'continuous')]#,
##                  ('100001', 'irnt', 'continuous')] # dummy list for testing
##    pheno_list = [("biomarkers", "30740", "both_sexes", "30740", "")]
##    pheno_list = [("categorical", "1150", "both_sexes", "", "3")]
   
#    ht = hl.import_table('gs://ukb-diverse-pops/ld_prune/release/phenotype_manifest.tsv',
#                         impute=True)
#    
#    pheno_list = list(zip(ht.trait_type.collect(),
#                          ht.phenocode.collect(), 
#                          ht.pheno_sex.collect(),
#                          ht.coding.collect(),
#                          ht.modifier.collect()))
#    
#    print(f'\nWARNING: For testing purposes, only using first 10 traits\n')
#    pheno_list = pheno_list[:10]
#    
    pheno_list = [
#            ('biomarkers', '30600', 'both_sexes', '', 'irnt'), # 6 pop LOO quant
#            ('categorical', '100240', 'both_sexes', '100240', ''), # 6 pop LOO binary
            ('continuous', '1408', 'both_sexes', '', ''),# 5 pop not_AFR quant
            ('icd10','Z37','both_sexes','',''), # 5 pop not_EUR binary
#            ('categorical', '1150', 'both_sexes', '3', ''), # 5 pop not_AFR binary
            ]
    print(pheno_list)
    
    return pheno_list

def make_pheno_id(trait_type, phenocode, pheno_sex, coding, modifier):
    pheno_id = (trait_type+'-'+
                phenocode+'-'+
                pheno_sex+
                ('-'+coding if len(coding)>0 else '')+
                ('-'+modifier if len(modifier)>0 else ''))
    return pheno_id    

def get_meta_sumstats(p, hail_script, pop, trait_type, phenocode, pheno_sex, 
                      coding, modifier, method):
    assert method in {'clump','sbayesr'}
    pheno_id = make_pheno_id(trait_type, phenocode, pheno_sex, coding, modifier)
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
        pop_list = [p for p in all_pops if p!=pop]
        five_pop_dir = f'gs://ukb-diverse-pops/ld_prune/release/{trait_category}/{"-".join(pop_list)}_batch1'
        if hl.hadoop_is_file(f'{five_pop_dir}/{filename}'):
            get_meta_ss = p.new_job(name=f'get_meta_ss_{pop}-{pheno_id}')
            get_meta_ss.command(' '.join(['set','-ex']))
            get_meta_ss.storage('2G')
            get_meta_ss.cpu(1)
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
                pass
        else:
            print(f'ERROR: No sumstats for {pheno_id} for not_{pop}')
            return None            
    return get_meta_ss.ofile


def get_adj_betas(p, pop, trait_type, phenocode, pheno_sex, coding, modifier, hail_script):
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
    output_dir = f'{temp_bucket}/results-test/not_{pop}/{make_pheno_id(trait_type, phenocode, pheno_sex, coding, modifier)}' # for testing

    clump_output_txt = f'{output_dir}/clumped_results-test.txt' # PLINK clump output txt file
    clump_output_ht = f'{output_dir}/clump_results-test.ht' # PLINK clump output hail table
    
    sbayesr_output_txt = f'{output_dir}/sbayesr_results-test.txt' # SBayesR output txt file
    sbayesr_output_ht = f'{output_dir}/sbayesr_results-test.ht' # SBayesR output hail table
    
    test_run = True
    
    if test_run:
        print('\n\nWARNING: Test run may overwrite existing results!\n')
    
    if not hl.hadoop_is_file(f'{clump_output_ht}/_SUCCESS') or test_run:
        ss = get_meta_sumstats(p=p,
                               hail_script=hail_script,
                               pop=pop, 
                               trait_type=trait_type, 
                               phenocode=phenocode, 
                               pheno_sex=pheno_sex, 
                               coding=coding,
                               modifier=modifier,
                               method='clump')
        
        if ss!=None:
            
#            head_ss = p.new_job(name='head_ss')
#            head_ss.command(' '.join(['set','-ex']))
##            head_ss.command(' '.join(['gunzip -c', str(ss), '|', 'head']))
#            head_ss.command(' '.join(['head',str(ss)]))
            
            run_method(p=p, 
                       pop=pop, 
                       trait_type=trait_type,
                       phenocode=phenocode,
                       pheno_sex=pheno_sex,
                       coding=coding,
                       modifier=modifier,  
                       hail_script=hail_script, 
                       output_txt=clump_output_txt,
                       output_ht=clump_output_ht,
                       ss=ss,
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
        
        
def run_method(p, pop, trait_type, phenocode, pheno_sex, coding, modifier, 
               hail_script, output_txt, output_ht, ss, method):
    r'''
    Runs either PLINK clump (method = 'clump') or SBayesR (method = 'sbayesr')
    '''
    assert method in {'clump','sbayesr'}
    
    task_suffix = f'not_{pop}-{make_pheno_id(trait_type, phenocode, pheno_sex, coding, modifier)}'
    # TODO: if method = 'sbayesr' check if LD matrix has already been calculated
    
    n_threads = 8
    
    tasks = []
        
    ## run plink clumping
    for chrom in range(22,23):
        ## read ref ld plink files 
        bfile = read_plink_input_group_chrom(p=p, 
                                             method=method,
                                             subset=f'not_{pop}',
                                             chrom=chrom)
        
        get_betas = p.new_job(name=f'{method}_{task_suffix}.chr{chrom}')
        
        # TODO: change image to include GCTB if running SBayesR?
        get_betas.storage('5G')
        get_betas.cpu(1) # plink clump cannot multithread
        
        get_betas.command(' '.join(['set','-ex']))
        
        if method == 'clump':
            clump_memory = -15*(chrom-1)+400 # Memory requested for PLINK clumping in MB. equation: -15*(chrom-1) + 400 is based on 400 MB for chr 1, 80 MB for chr 22
            get_betas.memory(f'{clump_memory}M') # default: 30G
#            get_betas.command(f' gunzip -c {ss} | grep "^{chrom}:\|SNP" |'+"awk '{ print $1,$7 }' | head")
#            get_betas.command(f'head {bfile.bim}')
#            get_betas.command(' '.join([f'gunzip -c {ss} | grep "^{chrom}:\|SNP" > {get_betas.ofile}_ss']))
            
            get_betas.command(' '.join(['plink',
                                        '--bfile', str(bfile),
                                        '--memory',str(clump_memory),
#                                        '--clump', f'{get_betas.ofile}_ss',
                                        '--clump', '<(','grep',f'"^{chrom}:\|SNP"', str(ss), ')',
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
#    tsv_to_ht = p.new_job(name=f'{method}_to_ht_{task_suffix}')
#    tsv_to_ht = tsv_to_ht.image('gcr.io/ukbb-diversepops-neale/hail_utils:3.4')
#    tsv_to_ht.storage('1G')
#    tsv_to_ht.memory('100M')
#    tsv_to_ht.cpu(n_threads)
#    tsv_to_ht.depends_on(get_betas_sink)
#    
#    tsv_to_ht.command(' '.join(['set','-ex']))
#    tsv_to_ht.command(' '.join(['PYTHONPATH=$PYTHONPATH:/',
#                        'PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=4g --conf spark.executor.memory=24g pyspark-shell"']))
#    tsv_to_ht.command(' '.join(['python3', str(hail_script),
#                         '--input_file', output_txt,
#                         '--tsv_to_ht',
#                         '--pop', pop,
#                         '--trait_type', trait_type,
#                         '--phenocode', phenocode,
#                         '--pheno_sex', pheno_sex,
#                         '--output_file', output_ht,
#                         '--n_threads', str(n_threads),
#                         '--overwrite']+
#                         (['--coding',coding] if coding != '' else [])+
#                         (['--modifier',modifier] if modifier != '' else [])))    

def main():
    pops = ['AFR','EUR'] # 'AFR']#
    
    backend = batch.ServiceBackend(billing_project='ukb_diverse_pops')
#    backend = batch.LocalBackend(tmp_dir='/tmp/batch/')
    
    p = batch.batch.Batch(name='test-clump', backend=backend,
                          default_image='gcr.io/ukbb-diversepops-neale/nbaya_plink:0.1',
                          default_storage='500Mi', default_cpu=8)
    
    ## download hail script to VM
    hail_script = p.read_input(f'{ldprune_dir}/scripts/python/plink_clump_hail.py')
    
    ## get phenotype master table
#    pheno_master = p.read_input(f'{ldprune_dir}/release/phenotype_master.tsv')
    
    for pop in pops:   
        pheno_list = get_pheno_list(pop)
    
        for trait_type, phenocode, pheno_sex, coding, modifier in pheno_list:
            
            get_adj_betas(p=p, 
                          pop=pop, 
                          trait_type=trait_type,
                          phenocode=phenocode,
                          pheno_sex=pheno_sex,
                          coding=coding,
                          modifier=modifier, 
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
    
    