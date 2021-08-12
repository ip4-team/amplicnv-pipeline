from typing import List

import pandas as pd
import argparse
import sys


def parse_args(args) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Add CNV classification (ClassifyCNV) status')
    parser.add_argument('-s', '--cnv-stats-path', type=str, required=True, help='CNV stats filename')
    parser.add_argument('-c', '--cnv-classification', type=str, required=True, help='CNV classification filename')
    parser.add_argument('-g', '--gene-list', type=str, help='target genes filename')
    parser.add_argument('-o', '--output', type=str, required=True, help='output filename')

    return parser.parse_args(args)


def get_classification(row: pd.Series, data: pd.DataFrame) -> str:
    c: pd.Series = data[(data['Chromosome'] == row['chrom']) &
                        (data['Start'] == row['chromStart']) &
                        (data['End'] == row['chromEnd'])]['Classification']
    assert len(c) == 1, 'only one CNV is allowed at the same position!'

    return c.values[0]


def get_gene_info(cnv_genes: str, gene_data: List) -> str:
    cnv_genes = cnv_genes.split(',')
    target_genes_in_cnv = [gene for gene in gene_data if gene in cnv_genes]
    return ','.join(target_genes_in_cnv) if target_genes_in_cnv else 'none'


if __name__ == '__main__':
    parsed_args = parse_args(sys.argv[1:])

    cnv_stats = pd.read_table(parsed_args.cnv_stats_path)
    class_data = pd.read_table(parsed_args.cnv_classification)

    cnv_stats['classification'] = cnv_stats.apply(lambda row: get_classification(row, class_data), axis=1)

    if parsed_args.gene_list:
        genes = list(pd.read_table(parsed_args.gene_list, header=None).iloc[:, 0].unique())
        cnv_stats['target genes'] = cnv_stats['genes'].apply(lambda cnv_genes: get_gene_info(cnv_genes, genes))

    cnv_stats.to_csv(parsed_args.output, sep='\t', index=False)
