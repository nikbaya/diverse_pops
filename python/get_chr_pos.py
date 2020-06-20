#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 13:50:20 2020

Extract chr and pos from SNP ID "chr:pos:a1:a2" string in order to tabix

@author: nbaya
"""

import hail as hl
import hailtop.batch as hb

bucket = 'gs://ukb-diverse-pops'
ldprune_dir = f'{bucket}/ld_prune'
loo_sumstats_dir = f'{ldprune_dir}/loo/sumstats/batch1'    
#scratch_dir = f'gs://ukbb-diverse-temp-30day/nb-tabix-scratch'
out_dir = f'{ldprune_dir}/loo/sumstats/batch2'

def get_paths():
    pheno_manifest = hl.import_table(f'{ldprune_dir}/phenotype_manifest.tsv.bgz', impute=True)
    pheno_manifest = pheno_manifest.filter(pheno_manifest.num_pops==6)
    filenames = pheno_manifest.filename.collect()
    all_files = hl.hadoop_ls(loo_sumstats_dir)
    print(len(all_files))
    all_paths = list(map(lambda f: f['path'], all_files))
    print(len(all_paths))
    filename_path_dict = dict(zip([path.split('/')[-1] for path in all_paths],all_paths))
    paths = [filename_path_dict[f] for f in filenames if f in filename_path_dict and filename_path_dict[f] in all_paths]
    print(len(paths))
    return paths
    

def annotate_chr_pos(b, path, use_tabix=False):
    pheno_id = path.replace('.tsv.bgz','').split('/')[-1]
    
    ss = b.read_input(path=path)
    
    j = b.new_job(f'get_chr_pos_{pheno_id}')
    j.command('set -ex')
    j.command(' '
            f'''
            gunzip -c {ss} | \\
            sed -e 's/SNP/chr:pos:a1:a2/g' \\
                -e 's/:/\t/g'  | \\
            awk 'BEGIN {{ OFS = "\t" }} {{ print $1=$1,$2=$2,$3=$3,$4=$4,$5=$1":"$2":"$3":"$4,$5,$6,$7,$8,$9,$10 }}' | \\
            sed 's/chr:pos:a1:a2/SNP/g' | \\
            bgzip > {j["bgz"]}{".tsv.bgz" if use_tabix else ""}
            '''
            )
    
    if use_tabix:
        j.command(
                f'''
                tabix -s 1 -b 2 -e 2 -c chr {j["bgz"]}.tsv.bgz
                ''')
        j.command(f'mv {j["bgz"]}.tsv.bgz {j["bgz"]}')
        j.command(f'mv {j["bgz"]}.tsv.bgz.tbi {j["tbi"]}')
        
        b.write_output(j['tbi'], f'{out_dir}_tabix/{pheno_id}.tsv.bgz.tbi')
        
    b.write_output(j['bgz'], f'{out_dir}/{pheno_id}.tsv.bgz')
    
    
def main():
    
    use_tabix = True
    
    hl.init(log='/Users/nbaya/Downloads/get_chr_pos.log')
    backend = hb.ServiceBackend(billing_project='ukb_diverse_pops',
                                bucket='ukbb-diverse-temp-30day/nb-batch-tmp')
    
    b = hb.Batch(name='get_chr_pos', backend=backend,
                 default_image='gcr.io/ukbb-diversepops-neale/nbaya_tabix:latest',
                 default_storage='2G', default_cpu=1)

    
    paths = get_paths()
    
    for path in paths:
        print(path)
        annotate_chr_pos(b=b,
                         path=path,
                         use_tabix=use_tabix)
    
    b.run(open=True)
    
    backend.close()
    
if __name__=='__main__':
        
    main()
    