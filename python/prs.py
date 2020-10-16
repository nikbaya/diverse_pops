#!/usr/bin/env python3

__author__ = 'konradk'

import argparse
from hail.linalg import BlockMatrix
# from gnomad.utils import slack
from ukbb_pan_ancestry import *


def main(args):
    hl.init(default_reference='GRCh37', log='/prs.log',
            spark_conf={'spark.hadoop.fs.gs.requester.pays.mode': 'AUTO', 'spark.hadoop.fs.gs.requester.pays.project.id': 'ukbb-diversepops-neale'})

    if args.prepare_sumstats_matrix:
        # get meta mt and separate by pop combo
        meta_mt = hl.read_matrix_table(get_meta_analysis_results_path())
        meta_mt = separate_results_mt_by_pop(meta_mt, 'meta_analysis_data', 'meta_analysis')
        meta_mt = meta_mt.annotate_cols(clump_pops=meta_mt.meta_analysis_data.pop)
        meta_mt = meta_mt.key_cols_by('clump_pops', *meta_mt.col_key)
        
        # get sumstats mt and separate by pop combo
        ss_mt = get_final_sumstats_mt_for_export()
        ss_mt = separate_results_mt_by_pop(ss_mt, 'pheno_data', 'summary_stats')
        ss_mt = ss_mt.annotate_cols(clump_pops=hl.array([ss_mt.pheno_data.pop]))
        ss_mt = ss_mt.key_cols_by(*meta_mt.col_key)
        
        # join meta results and sumstats mt
        # NOTE: union_cols() requires the same entry fields schema
        meta_mt = meta_mt.select_entries(BETA = meta_mt.meta_analysis.BETA).select_cols().select_rows()
        ss_mt = ss_mt.select_entries(BETA = ss_mt.summary_stats.BETA).select_cols().select_rows()
        mt = meta_mt.union_cols(ss_mt)
        
        # filter to distinct cols
        # NOTE: distinct_by_col() does not allow a col key of type `list`
        mt = mt.annotate_cols(clump_pops_str = hl.delimit(mt.clump_pops)).key_cols_by('clump_pops_str', *[k for k in mt.col_key if k!='clump_pops']).distinct_by_col()
        mt = mt.distinct_by_col()
        
        # ensure that betas are not missing
        ss_mt = ss_mt.annotate_cols(clump_pops_str = hl.delimit(ss_mt.clump_pops)).key_cols_by('clump_pops_str', *[k for k in ss_mt.col_key if k!='clump_pops'])
        mt = mt.annotate_entries(BETA = hl.or_else(mt.BETA, ss_mt[mt.row_key, mt.col_key].BETA))
        
        # read clump mt and separate by pop combo
        clump_mt = hl.read_matrix_table(get_clumping_results_path(high_quality=args.high_quality))
        clump_mt = separate_results_mt_by_pop(clump_mt, 'clump_pops', 'plink_clump', skip_drop=True)
        clump_mt = clump_mt.annotate_cols(clump_pops_str = hl.delimit(clump_mt.clump_pops))
        clump_mt = clump_mt.drop('clump_pops').key_cols_by(*mt.col_key)
        
        # join sumstats/meta-analysis with clump mt
        mt = all_axis_join(mt, clump_mt)
        
        # filter to cols that only exist in clump_mt (pop_index is a field from clump_mt)
        mt = mt.filter_cols(hl.is_defined(mt.pop_index))
        
        mt = explode_by_p_threshold(mt).unfilter_entries()
        # Write pheno data for later use
        mt.add_col_index('idx').key_cols_by('idx').cols().write(get_clump_sumstats_col_ht_path(args.high_quality), args.overwrite)
        BlockMatrix.write_from_entry_expr(
            hl.or_else(mt.meta_analysis.BETA * hl.is_defined(mt.plink_clump.TOTAL) * hl.int(mt.meta_analysis.Pvalue < mt.p_threshold), 0.0),
            get_clump_sumstats_bm_path(args.high_quality), args.overwrite)
        # 2020-06-25 01:49:32 Hail: INFO: Wrote all 7078 blocks of 28987534 x 3530 matrix with block size 4096.
        # If clump_mt is significantly smaller than meta_mt, consider putting that on the left of the join,
        # then filter the genotype matrix to only those SNPs (pilot would go from 28.9M -> 21.2M)

    if args.prepare_genotype_matrix:
        meta_mt = hl.read_matrix_table(get_meta_analysis_results_path())
        mt = get_filtered_mt_with_x()
        mt = mt.filter_rows(hl.is_defined(meta_mt.rows()[mt.row_key]))
        # Write sample data for later use
        mt = mt.key_cols_by(userId=hl.int32(mt.s))
        mt.cols().add_index().write(genotype_samples_ht_path, args.overwrite)
        BlockMatrix.write_from_entry_expr(mt.dosage, genotype_bm_path, args.overwrite)
        # 2020-06-25 19:18:14 Hail: INFO: Wrote all 764424 blocks of 28987534 x 441345 matrix with block size 4096.

    if args.compute_prs:
        sumstats_bm = BlockMatrix.read(get_clump_sumstats_bm_path(args.high_quality))
        genotype_bm = BlockMatrix.read(genotype_bm_path)
        prs_bm: BlockMatrix = genotype_bm.T @ sumstats_bm
        prs_bm.write(get_prs_bm_path(args.high_quality), args.overwrite)
        # Pilot of 353 phenos (3530 columns), ~58 hours
        # 2020-06-25 19:18:14 Hail: INFO: Wrote all 764424 blocks of 28987534 x 441345 matrix with block size 4096.
        # 2020-06-28 05:27:54 Hail: INFO: wrote matrix with 441345 rows and 3530 columns as 108 blocks of size 4096

    if args.create_prs_mt:
        prs_bm = BlockMatrix.read(get_prs_bm_path(args.high_quality))
        pheno_ht = hl.read_table(get_clump_sumstats_col_ht_path(args.high_quality)).key_by('idx')
        samples_ht = hl.read_table(genotype_samples_ht_path).key_by('idx')
        mt = BlockMatrix.to_matrix_table_row_major(prs_bm, n_partitions=200).rename({'element': 'score'})
        mt = mt.annotate_cols(**pheno_ht[mt.col_key]).key_cols_by(*PHENO_KEY_FIELDS)
        mt = mt.annotate_rows(**samples_ht[mt.row_key]).key_rows_by('userId')
        mt.write(get_prs_mt_path(args.high_quality), args.overwrite)

    if args.assess_prs:
        prs_mt = hl.read_matrix_table(get_prs_mt_path(args.high_quality))
        pheno_mt = get_ukb_pheno_mt()  # TODO: fix all phenos to new keying scheme
        pheno_mt = pheno_mt.key_cols_by(
            **pheno_mt.col_key.annotate(modifier=hl.if_else(pheno_mt.trait_type == "biomarkers", "irnt", pheno_mt.modifier)))
        mt = prs_mt.annotate_entries(**pheno_mt[prs_mt.row_key, prs_mt.col_key])
        mt = mt.filter_cols(mt.description == 'Type 2 diabetes')
        import gnomad.utils.file_utils
        gnomad.utils.file_utils.select_primitives_from_ht(mt.entries()).export(f'{temp_bucket}/prs/prs_test.txt.bgz')
        # mt = mt.annotate_cols(prs_corr=hl.agg.linreg(mt.both_sexes, [1.0, mt.score]))
        # ht = mt.cols().checkpoint(get_prs_assess_ht_path(args.high_quality), args.overwrite, _read_if_exists=not args.overwrite)
        # print(ht.aggregate(hl.agg.fraction(ht.prs_corr.p_value[1] < 0.05)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--prepare_sumstats_matrix', help='Prepare summary stats blockmatrix', action='store_true')
    parser.add_argument('--prepare_genotype_matrix', help='Prepare genotype blockmatrix', action='store_true')
    parser.add_argument('--compute_prs', help='Compute PRS', action='store_true')
    parser.add_argument('--high_quality', help='Overwrite everything', action='store_true')
    parser.add_argument('--create_prs_mt', help='Convert PRS blockmatrix to MT', action='store_true')
    parser.add_argument('--assess_prs', help='Assess PRS performance', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    # if args.slack_channel:
    #     from slack_token_pkg.slack_creds import slack_token
    #     with slack.slack_notifications(slack_token, args.slack_channel):
    main(args)
