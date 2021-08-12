import pandas as pd
import argparse
import sys


def parse_args(args) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Create a file to be used as input for Ensembl Variant Effect '
                                                 'Predictor (VEP)')
    parser.add_argument('-i', '--input', type=str, required=True, help='input file')
    parser.add_argument('-o', '--output', type=str, required=True, help='output file')
    return parser.parse_args(args)


if __name__ == '__main__':
    parsed_args = parse_args(sys.argv[1:])

    df = pd.read_csv(parsed_args.input, sep='\t')

    blocks_vep = df.loc[:, ['chrom', 'chromStart', 'chromEnd']]
    blocks_vep['chrom'] = blocks_vep['chrom'].apply(lambda x: x.strip('chr'))
    blocks_vep['call'] = df['call'].apply(lambda x: 'DUP' if x == 'gain' else 'DEL')
    blocks_vep['strand'] = blocks_vep['chrom'].apply(lambda x: '.')
    blocks_vep['id'] = df['cnv id']

    blocks_vep.to_csv(parsed_args.output, sep='\t', index=False, header=False)
