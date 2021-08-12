import pandas as pd
import argparse
import sys


def parse_args(args) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Filter targets file by genes')
    parser.add_argument('-i', '--input', type=str, required=True, help='input file')
    parser.add_argument('-o', '--output', type=str, required=True, help='output file')
    parser.add_argument('-g', '--gene-list', type=str, required=True, help='gene list file')
    return parser.parse_args(args)


if __name__ == '__main__':
    parsed_args = parse_args(sys.argv[1:])
    gene_list = pd.read_csv(parsed_args.gene_list, header=None, names=['symbol'])
    targets = pd.read_table(parsed_args.input, index_col=False, header=0)
    filtered = targets[targets['gene'].isin(gene_list['symbol'])]
    filtered.to_csv(parsed_args.output, sep='\t', index=False)


