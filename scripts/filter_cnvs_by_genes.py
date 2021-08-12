import pandas as pd
import argparse
import sys


def parse_args(args) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Filter CNVs file by genes')
    parser.add_argument('-i', '--input', type=str, required=True, help='input file')
    parser.add_argument('-o', '--output', type=str, required=True, help='output file')
    parser.add_argument('-g', '--gene-list', type=str, required=True, help='gene list file')
    return parser.parse_args(args)


if __name__ == '__main__':
    parsed_args = parse_args(sys.argv[1:])
    gene_list = pd.read_csv(parsed_args.gene_list, header=None, names=['symbol'])
    cnv_data = pd.read_csv(parsed_args.input, sep='\t')
    filtered = cnv_data[
        cnv_data['genes'].apply(lambda x: sum([gene in gene_list['symbol'].unique() for gene in x.split(',')]) > 0)]
    filtered.to_csv(parsed_args.output, sep='\t', index=False)
