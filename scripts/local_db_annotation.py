from typing import List
from pybedtools import BedTool

import pandas as pd
import argparse
import sys
import glob


def parse_args(args) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Annotate CNVs regarding their frequency on an in-house CNV db')
    parser.add_argument('-b', '--baseline-path', type=str, required=True, help='unix style pathname pattern')
    # parser.add_argument('-r', '--regex', type=str, required=True, help='regex for sample ids')
    parser.add_argument('-t', '--test-path', type=str, required=True, help='test sample filename')
    parser.add_argument('-c', '--cnv-targets', type=str, required=True, help='')
    parser.add_argument('-a', '--annotation-output', type=str, required=True, help='annotation output filename')
    parser.add_argument('-g', '--gene-output', type=str, required=True, help='gene output filename')
    parser.add_argument('-f', '--min-overlap-fraction', type=float, required=False, default=0.5,
                        help='minimum overlap required as a fraction of CNVs in test sample')

    return parser.parse_args(args)


def get_sortable_chrom(chrom: str) -> str:
    chrom = chrom.strip('chr')
    if chrom.isdigit():
        return chrom if int(chrom) > 9 else f'0{chrom}'
    else:
        return chrom


def create_bed(filenames: List) -> BedTool:
    c = ['chrom', 'chromStart', 'chromEnd', 'call', 'cnv id']
    df_list = [pd.read_table(filename, usecols=c) for filename in filenames]
    data = pd.concat(df_list)
    data['sortable chrom'] = data['chrom'].apply(lambda x: get_sortable_chrom(x))
    data = data.sort_values(by=['sortable chrom', 'chromStart'])
    return BedTool(data[c].to_string(header=False, index=False, index_names=False), from_string=True)


if __name__ == '__main__':
    parsed_args = parse_args(sys.argv[1:])

    min_overlap_fraction = parsed_args.min_overlap_fraction

    columns = ['chrom', 'start', 'end', 'call', 'cnv id']
    cov_columns = ['# overlaps', 'covered bases', 'length', 'fraction']

    targets = pd.read_table(parsed_args.cnv_targets, header=None,
                            usecols=[0, 1, 2, 3], names=['chrom', 'start', 'end', 'gene'])

    baseline_files = glob.glob(parsed_args.baseline_path)
    grouped_baseline_bed = create_bed(baseline_files)
    baseline_beds = [create_bed([filename]) for filename in baseline_files]
    print(f'Loaded {len(baseline_files)} baseline files...')

    sample_bed = create_bed([parsed_args.test_path])

    # compute coverage considering all cnvs in the baseline
    unified_coverage = sample_bed.coverage(grouped_baseline_bed,
                                           f=min_overlap_fraction).to_dataframe(disable_auto_names=True,
                                                                                names=columns + cov_columns)
    # compute coverage separately
    coverage_beds = [sample_bed.coverage(bed, f=min_overlap_fraction) for bed in baseline_beds]
    all_coverages = pd.concat([bed.to_dataframe(disable_auto_names=True,
                                                names=columns + cov_columns) for bed in coverage_beds if bed])
    # compute overlap frequency in pop
    unified_coverage['# samples'] = unified_coverage['cnv id'].apply(
        lambda cnv_id: all_coverages[(all_coverages['cnv id'] == cnv_id) &
                                     (all_coverages['fraction'] > 0)]['cnv id'].count())

    unified_coverage['frequency (local db)'] = unified_coverage['# samples'] / len(baseline_files) * 100
    unified_coverage['% overlapped bases (local db)'] = unified_coverage['fraction'] * 100

    # get unique regions on CNV
    unique = sample_bed.subtract(grouped_baseline_bed)

    unified_coverage['# non-overlapped bases (local db)'] = \
        unified_coverage['length'] - unified_coverage['covered bases']

    genes_in_unique = BedTool(targets.to_string(header=False, index=False, index_names=False),
                              from_string=True).intersect(unique, wb=True) \
        .to_dataframe(disable_auto_names=True,
                      names=['chrom_A', 'start_A', 'end_A', 'gene'] + columns)
    genes_in_unique['length'] = genes_in_unique['end_A'] - genes_in_unique['start_A']

    unified_coverage['genes in non-overlapped bases (local db)'] = unified_coverage['cnv id'].apply(
        lambda cnv_id: ','.join(genes_in_unique[genes_in_unique['cnv id'] == cnv_id]['gene'].unique()))
    unified_coverage['genes in non-overlapped bases (local db)'] = \
        unified_coverage['genes in non-overlapped bases (local db)'].apply(lambda genes: genes if genes else 'none')

    cov_data: pd.DataFrame = unified_coverage[['cnv id',
                                               'frequency (local db)',
                                               '% overlapped bases (local db)',
                                               '# non-overlapped bases (local db)',
                                               'genes in non-overlapped bases (local db)']]

    sample_df = pd.read_table(parsed_args.test_path)
    annotated = sample_df.merge(cov_data, on='cnv id')

    annotated.to_csv(parsed_args.annotation_output, sep='\t', index=False)

    gene_data: pd.DataFrame = genes_in_unique[['chrom_A', 'start_A', 'end_A', 'length', 'gene', 'call', 'cnv id']]
    gene_data.columns = ['chrom', 'start', 'end', 'length', 'gene', 'call', 'cnv id']
    gene_data.to_csv(parsed_args.gene_output, sep='\t', index=False)
