import pandas as pd
import argparse
import sys


def parse_args(args) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Create a file to be used as input for PhenogramViz (Cytoscape)')
    parser.add_argument('-i', '--input', type=str, required=True, help='input file')
    parser.add_argument('-o', '--output', type=str, required=True, help='output file')
    parser.add_argument('-f', '--filtered-output', type=str, required=True, help='filtered output file')
    parser.add_argument('-s', '--min-score', type=float, default=0.76, help='min score value')
    parser.add_argument('-p', '--p-value', type=float, default=0.05, help='adjusted p-value cutoff')
    return parser.parse_args(args)


def get_classify_cnv_columns(data: pd.DataFrame) -> pd.DataFrame:
    blocks_vep = data.loc[:, ['chrom', 'chromStart', 'chromEnd']]
    blocks_vep['chrom'] = blocks_vep['chrom'].apply(lambda x: x.strip('chr'))
    blocks_vep['call'] = data['call'].apply(lambda x: 'DUP' if x == 'gain' else 'DEL')
    return blocks_vep


if __name__ == '__main__':
    parsed_args = parse_args(sys.argv[1:])
    df = pd.read_csv(parsed_args.input, sep='\t')
    get_classify_cnv_columns(df).to_csv(parsed_args.output, sep='\t', index=False, header=False)

    filtered = df[(df['score'] > parsed_args.min_score) &
                  (df['adjusted p-value (Bonferroni)'] < parsed_args.p_value)]
    get_classify_cnv_columns(filtered).to_csv(parsed_args.filtered_output, sep='\t', index=False, header=False)
