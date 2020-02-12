#!/usr/bin/env python3

__author__ = 'nbaya'

import time
#from gnomad_hail import *
#from shlex import quote as shq
#from ukb_common.utils.saige_pipeline import *
#from ukb_exomes import *
import hail as hl

bucket = 'gs://ukbb-pharma-exome-analysis'
temp_bucket = 'gs://ukbb-pharma-exome-analysis-temp'
results_dir = f'{bucket}/finngen/result'
final_results_ht = f'{bucket}/finngen/variant_results.ht'
vep_path = f'{bucket}/finngen/all_variants_vep.ht'


def main(args):
    hl.init(default_reference='GRCh38', log='/load_results.log')
    start_time = time.time()
    all_phenos_ht = hl.import_table('gs://finngen-public-data-r2/summary_stats/r2_manifest.tsv', impute=True)
    # all_phenos_ht = all_phenos_ht.annotate(code=all_phenos_ht.phenocode.split('_', 2)[0])
    all_phenos = all_phenos_ht.collect()

    backend = pipeline.BatchBackend(billing_project='ukb_pharma')
    # backend = pipeline.LocalBackend(gsa_key_file='/Users/konradk/.hail/ukb-diverse-pops.json')
    p = pipeline.Pipeline(name='finngen_load', backend=backend,
                          default_image='gcr.io/ukbb-exome-pharma/hail_utils:3.3',
                          default_storage='500Mi', default_cpu=8)

    tasks = []
    for i, pheno in enumerate(all_phenos):
        variant_results_ht_path = f'{results_dir}/ht/{pheno.phenocode}.ht'
        if not args.overwrite_results and hl.hadoop_exists(f'{variant_results_ht_path.replace(".ht", ".mt")}/_SUCCESS'):
            continue
        t: pipeline.pipeline.Task = p.new_task(name='load_pheno', attributes={'pheno': pheno.phenocode}).cpu(args.n_threads)
        t.command(f"""
        PYTHONPATH=$PYTHONPATH:/ PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=24g pyspark-shell"
        python3 /ukb_exomes/hail/load_finngen_results_hail.py
        --input_file {pheno.path_bucket} --n_threads {args.n_threads}
        --load_single --vep_path {vep_path}
        --additional_dict {shq(json.dumps(dict(pheno)))}
        --output_ht {variant_results_ht_path}
        --output_mt {variant_results_ht_path.replace('.ht', '.mt')}
        --overwrite
        """.replace('\n', ' '))
        tasks.append(t)
        if args.limit and i == args.limit:
            break

    t: pipeline.pipeline.Task = p.new_task(name='combine').cpu(args.n_threads)

    t.depends_on(*tasks)
    t.command(f"""
    PYTHONPATH=$PYTHONPATH:/ PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=4g --conf spark.executor.memory=24g pyspark-shell"
    python3 /ukb_exomes/hail/load_finngen_results_hail.py --combine_all
    --input_directory {results_dir}/ht
    --output_ht {final_results_ht}
    --output_mt {final_results_ht.replace('.ht', '.mt')}
    --overwrite --n_threads {args.n_threads}
    """.replace('\n', ' '))

    logger.info(f'Setup took: {time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))}')
    logger.info(f'Submitting: {get_tasks_from_pipeline(p)}')
    p.run(dry_run=args.dry_run, verbose=True, delete_scratch_on_exit=False)
    logger.info(f'Finished: {get_tasks_from_pipeline(p)}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--limit', help='Overwrite everything', type=int)
    parser.add_argument('--dry_run', help='Overwrite everything', action='store_true')
    parser.add_argument('--overwrite_results', help='Overwrite everything', action='store_true')
    parser.add_argument('--n_threads', help='Overwrite everything', default=8, type=int)
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@nbaya')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)