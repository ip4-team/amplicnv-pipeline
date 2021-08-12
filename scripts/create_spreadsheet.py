import pandas as pd
import argparse
import sys


def parse_args(args) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Create a Excel spreadsheet from input files')
    parser.add_argument('-c', '--cnvs', type=str, required=True, help='CNVs filename')
    parser.add_argument('-t', '--targets', type=str, required=True, help='targets filename')
    parser.add_argument('-f', '--files', type=str, action='append',
                        required=False, default=[], help='other files')
    parser.add_argument('-o', '--output', type=str, required=True, help='output filename')

    return parser.parse_args(args)


if __name__ == '__main__':
    parsed_args = parse_args(sys.argv[1:])
    files = [parsed_args.cnvs, parsed_args.targets] + parsed_args.files
    tables = [pd.read_table(filename, index_col=False, header=0) for filename in files]

    with pd.ExcelWriter(parsed_args.output) as writer:
        for i, table in enumerate(tables):
            sheet_name = files[i].split('/')[-1].split('.')[0]
            table.to_excel(writer, sheet_name=sheet_name, index=False)
