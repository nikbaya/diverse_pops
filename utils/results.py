
import hail as hl

bucket = 'gs://ukb-diverse-pops'
public_bucket = 'gs://ukb-diverse-pops-public/'

def get_gene_intervals_path(reference: str = 'GRCh37'):
    return f'{public_bucket}/misc/gene_intervals_{reference}.ht'


def get_variant_results_path(pop: str, extension: str = 'mt'):
    return f'{bucket}/combined_results/results_{pop}.{extension}'


def get_variant_results_qc_path(extension: str = 'ht'):
    return f'{bucket}/combined_results/full_variant_qc_metrics.{extension}'


def get_meta_analysis_results_path(extension: str = 'mt'):
    return f'{bucket}/combined_results/meta_analysis.{extension}'


def get_phenotype_results_qc_path(extension: str = 'ht'):
    return f'{bucket}/combined_results/full_phenotype_qc_metrics.{extension}'


def get_analysis_data_path(subdir: str, dataset: str, pop: str, extension: str = 'txt.bgz'):
    return f'{public_bucket}/sumstats_qc_analysis/{subdir}/{dataset}_{pop}.{extension}'


def get_results_timing_tsv_path(timing_type: str, pop: str = ''):
    check_timing_type(timing_type)

    pop = f'_{pop}' if timing_type == 'saige' else ''

    return f'{bucket}/results/misc/timings_{timing_type}{pop}.txt'


def get_results_timing_ht_path(timing_type: str):
    check_timing_type(timing_type)
    return f'{bucket}/results/misc/timings_{timing_type}.ht'


def get_heritability_txt_path(from_date: str = None):
    return f'{bucket}/results/misc/all_heritabilities{"_" + from_date if from_date else ""}.txt'


def annotate_nearest_gene(t, add_contig: bool = False, add_only_gene_symbols_as_str: bool = False, loc: str = 'nearest_genes'):
    intervals_ht = hl.read_table(get_gene_intervals_path())
    if add_contig:
        intervals_ht = intervals_ht.annotate(contig=intervals_ht.interval.start.contig)
    annotation = intervals_ht.index(t.locus, all_matches=True)
    if add_only_gene_symbols_as_str:
        annotation = hl.delimit(annotation.gene_name)
    if loc: annotation = {loc: annotation}
    return t.annotate_rows(**annotation) if isinstance(t, hl.MatrixTable) else t.annotate(**annotation)


def filter_lambda_gc(lambda_gc):
    return (lambda_gc > 0.8) & (lambda_gc < 1.2)


def load_final_sumstats_mt(filter_phenos: bool = True, filter_variants: bool = True,
                           filter_sumstats: bool = True, separate_columns_by_pop: bool = True,
                           annotate_with_nearest_gene: bool = True):
    mt = hl.read_matrix_table(get_variant_results_path('full', 'mt'))
    variant_qual_ht = hl.read_table(get_variant_results_qc_path())
    mt = mt.annotate_rows(**variant_qual_ht[mt.row_key])
    pheno_qual_ht = hl.read_table(get_analysis_data_path('lambda', 'lambdas', 'full', 'ht'))
    mt = mt.annotate_cols(**pheno_qual_ht[mt.col_key])

    if filter_phenos:
        keep_phenos = hl.zip_with_index(mt.pheno_data).filter(
            lambda x: filter_lambda_gc(x[1].lambda_gc))

        mt = mt.annotate_cols(
            pheno_indices=keep_phenos.map(lambda x: x[0]),
            pheno_data=keep_phenos.map(lambda x: x[1]))
        mt = mt.annotate_entries(
            summary_stats=hl.zip_with_index(mt.summary_stats).filter(
                lambda x: mt.pheno_indices.contains(x[0])).map(lambda x: x[1])
        )
        mt = mt.filter_cols(hl.len(mt.pheno_data) > 0)

    if filter_sumstats:
        mt = mt.annotate_entries(summary_stats=mt.summary_stats.map(
            lambda x: hl.or_missing(~x.low_confidence, x)
        ))
        mt = mt.filter_entries(~mt.summary_stats.all(lambda x: hl.is_missing(x.Pvalue)))

    if filter_variants:
        mt = mt.filter_rows(mt.high_quality)

    if annotate_with_nearest_gene:
        mt = annotate_nearest_gene(mt)

    if separate_columns_by_pop:
        mt = separate_results_mt_by_pop(mt)

    return mt


def separate_results_mt_by_pop(mt):
    mt = mt.annotate_cols(pheno_data=hl.zip_with_index(mt.pheno_data)).explode_cols('pheno_data')
    mt = mt.annotate_cols(pop_index=mt.pheno_data[0], pheno_data=mt.pheno_data[1])
    mt = mt.annotate_entries(summary_stats=mt.summary_stats[mt.pop_index]).drop('pop_index')
    return mt
