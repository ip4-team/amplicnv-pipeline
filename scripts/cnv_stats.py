import pandas as pd
import argparse
import sys
from scipy import stats


def parse_args(args) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Annotate CNVs regarding their quality and reliability.')
    parser.add_argument('-i', '--input', type=str, required=True, help='input file')
    parser.add_argument('-o', '--output', type=str, required=True, help='output file')
    return parser.parse_args(args)


if __name__ == '__main__':
    parsed_args = parse_args(sys.argv[1:])

    df = pd.read_table(parsed_args.input, index_col=False, header=0)

    df['cnv id'] = df['cnv id'].apply(lambda x: x.split('; '))

    df = df.explode('cnv id')

    grouped_by_cnv = df.groupby('cnv id')
    blocks = grouped_by_cnv[['chrom', 'chromStart']].min()
    blocks['chromEnd'] = grouped_by_cnv[['chromEnd']].max()
    blocks['length'] = blocks['chromEnd'] - blocks['chromStart']
    blocks['genes'] = grouped_by_cnv[['gene']].agg(lambda x: ','.join(set(sorted(x))))
    blocks['# targets'] = grouped_by_cnv[['ratio']].count()

    blocks['% targets > 0.76 map'] = df[df['mappability'] > 0.76].groupby('cnv id')[[
        'ratio']].count() / grouped_by_cnv[['ratio']].count()
    blocks['% targets > 0.76 map'] = blocks['% targets > 0.76 map'].fillna(0).round(2)

    blocks['% median (baseline) > 25'] = df[df['# median (baseline)'] > 25].groupby('cnv id')[[
        'ratio']].count() / grouped_by_cnv[['ratio']].count()
    blocks['% median (baseline) > 25'] = blocks['% median (baseline) > 25'].fillna(0).round(2)

    na_dropped = df.dropna(subset=['% mean (std)'])
    blocks['% std/mean < 0.5'] = na_dropped[na_dropped['% mean (std)'] < 50].groupby('cnv id')[[
        'ratio']].count() / na_dropped.groupby('cnv id')[['ratio']].count()
    blocks['% std/mean < 0.5'] = blocks['% std/mean < 0.5'].fillna(0).round(2)

    na_dropped = df.dropna(subset=['% median (mad)'])
    blocks['% mad/median < 0.4'] = na_dropped[na_dropped['% median (mad)'] < 40].groupby('cnv id')[['ratio']].count() / \
                                   na_dropped.groupby('cnv id')[['ratio']].count()
    blocks['% mad/median < 0.4'] = blocks['% mad/median < 0.4'].fillna(0).round(2)

    blocks['score'] = ((blocks['% targets > 0.76 map'] +
                        blocks['% median (baseline) > 25'] +
                        blocks['% mad/median < 0.4']) / 3).round(2)

    blocks.loc[:, ['% targets > 0.76 map',
                   '% median (baseline) > 25',
                   '% std/mean < 0.5',
                   '% mad/median < 0.4']] = blocks[
                                                    ['% targets > 0.76 map',
                                                     '% median (baseline) > 25',
                                                     '% std/mean < 0.5',
                                                     '% mad/median < 0.4']
                                                ] * 100

    blocks['sum # reads norm test'] = grouped_by_cnv[['# reads norm (test)']].apply(sum)
    blocks['sum # median norm (baseline)'] = grouped_by_cnv[['# median norm (baseline)']].apply(sum)

    blocks['cnv ratio'] = blocks['sum # reads norm test'] / blocks['sum # median norm (baseline)']

    blocks['call'] = blocks['cnv ratio'].apply(lambda x: 'gain' if x > 1 else 'loss')

    cnv_ids = df['cnv id'].unique()
    for cnv_id in cnv_ids:
        cnvs = df[df['cnv id'] == cnv_id]
        blocks.loc[cnv_id, 'p-value'] = stats.wilcoxon(cnvs['# median norm (baseline)'],
                                                       cnvs['# reads norm (test)'])[1]

    blocks['adjusted p-value (Bonferroni)'] = blocks['p-value'] * len(cnv_ids)

    blocks = blocks.sort_values(by=['score'], ascending=False)

    blocks.to_csv(parsed_args.output, sep='\t')
