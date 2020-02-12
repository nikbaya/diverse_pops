from ukb_common import *
from ukb_exomes import *

temp_bucket = 'gs://ukbb-pharma-exome-analysis-temp'
temp_mt_path = f'{temp_bucket}/finngen.mt'


def main(args):
    add_args = {}
    if args.n_threads is not None:
        add_args['master'] = f'local[{args.n_threads}]'
    hl.init(default_reference='GRCh38', log='/load_finngen.log', **add_args)

    if args.load_single:
        ht = hl.import_table(args.input_file, impute=True, force_bgz=True, min_partitions=100).rename({'#chrom': 'chrom'})
        ht = ht.transmute(locus=hl.locus('chr' + ht.chrom, ht.pos), alleles=[ht.ref, ht.alt]).key_by('locus', 'alleles')
        ht = ht.transmute(Pvalue=ht.pval).annotate_globals(**json.loads(args.additional_dict))
        ht = ht.annotate(**get_vep_formatted_data(args.vep_path)[ht.key])
        ht = ht.checkpoint(args.output_ht, overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        ht = ht.select_globals().annotate(**json.loads(args.additional_dict))
        mt = ht.to_matrix_table(['locus', 'alleles'], ['phenocode'],
                                ['rsids', 'nearest_genes', 'gene', 'annotation'],
                                ['category', 'name', 'n_cases', 'n_controls'])
        mt.checkpoint(args.output_mt, overwrite=args.overwrite, _read_if_exists=not args.overwrite)

    if args.combine_all:
        # all_hts = list(filter(lambda y: y.endswith('.ht'), map(lambda x: x['path'], hl.hadoop_ls(args.input_directory))))
        # print(f'Got {len(all_hts)} HTs...')
        # mt = mwzj_hts_by_tree(all_hts, temp_bucket + '/finngen', ['phenocode'], debug=True)
        # mt.checkpoint(temp_mt_path, overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        mt = hl.read_matrix_table(temp_mt_path)
        mt.naive_coalesce(5000).write(args.output_mt, args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_single', help='Overwrite everything', action='store_true')
    parser.add_argument('--combine_all', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_ht', help='Overwrite everything', action='store_true')
    parser.add_argument('--n_threads', help='Overwrite everything')
    parser.add_argument('--input_file', help='Overwrite everything')
    parser.add_argument('--vep_path', help='Overwrite everything')
    parser.add_argument('--additional_dict', help='Overwrite everything')
    parser.add_argument('--output_ht', help='Overwrite everything')
    parser.add_argument('--output_mt', help='Overwrite everything')
    parser.add_argument('--input_directory', help='Overwrite everything')
    args = parser.parse_args()

    main(args)