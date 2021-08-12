import argparse
import sys

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as mticker


def get_state(ratio, q1, q3, iqr) -> str:
    if ratio < q1 - 1.5 * iqr:
        return 'loss'
    if ratio > q3 + 1.5 * iqr:
        return 'gain'
    return 'neutral'


def get_mappability(map_data: pd.DataFrame, row: pd.Series) -> float:
    chrom = row['chrom']
    chrom_start = row['chromStart']
    chrom_end = row['chromEnd']
    t_map = map_data.query('(chrom == @chrom) & (chromStart <= @chrom_start) & (chromEnd >= @chrom_end)')['mappability']
    return float(t_map.round(2))


def get_target_quality(row: pd.Series) -> str:
    is_dispersion_high = row['mad/median'] >= 0.4
    is_mappability_bad = row['mappability'] < 0.76
    is_low_rd = row['median'] < 25
    if is_low_rd:
        return 'low baseline read-depth'
    if is_mappability_bad and is_mappability_bad:
        return 'poor map & high var'
    if is_mappability_bad:
        return 'poor mappability'
    if is_dispersion_high:
        return 'high variability'
    return 'stable'


def plot_gene(selected_gene: str, nrr_table: pd.DataFrame, map_table: pd.DataFrame,
              q1: float, q3: float, iqr: float, output_path: str) -> None:
    selected = nrr_table[nrr_table['gene'] == selected_gene].copy()
    selected['mappability'] = selected.apply(lambda row: get_mappability(map_table, row), axis=1)
    selected['copy state'] = selected['ratio'].apply(
        lambda r: get_state(r, q1, q3, iqr))
    selected['target quality'] = selected.apply(lambda row: get_target_quality(row), axis=1)
    selected['x'] = [f'T{i}' for i in range(1, len(selected) + 1)]
    selected['2x ratio'] = selected['ratio'].round(1) * 2

    with sns.axes_style('darkgrid'):
        if len(selected) > 18:
            plt.figure(figsize=(12, 4))

        palette = {'neutral': 'royalblue',
                   'loss': 'orangered',
                   'gain': 'forestgreen'}
        markers = {'stable': 'o',
                   'poor mappability': 'X',
                   'high variability': 'P',
                   'poor map & high var': 'd',
                   'low baseline read-depth': 'v'}

        ax = sns.scatterplot(y=selected['2x ratio'], x=selected['x'],
                             hue=selected['copy state'], palette=palette,
                             style=selected['target quality'], markers=markers)
        ax.set(ylim=(-0.1, 4))
        ax.axhline((q1 - 1.5 * iqr) * 2, ls='--', linewidth=1, alpha=.7)
        ax.axhline((q3 + 1.5 * iqr) * 2, ls='--', linewidth=1, alpha=.7)
        ax.set_ylabel('copy state', fontsize=13)
        ax.set_xlabel('targets', fontsize=13)
        ticks_loc = ax.get_yticks().tolist()
        ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
        ax.set_yticklabels(['{:.0f}'.format(i) if i in [0, 1.0, 2.0, 3.0, 4.0] else '' for i in ticks_loc])
        if len(selected) > 15:
            plt.xticks(rotation=45)

        copy_legend_elements = [
            Line2D([0], [0], label='neutral', marker='o', markersize=8, color='#eaeaf2',
                   markerfacecolor=palette['neutral']),
            Line2D([0], [0], label='loss', marker='o', markersize=8, color='#eaeaf2',
                   markerfacecolor=palette['loss']),
            Line2D([0], [0], label='gain', marker='o', markersize=8, color='#eaeaf2',
                   markerfacecolor=palette['gain']),
        ]
        ax.set_title(r'gene: $\it' + selected_gene + '$', fontsize=14)
        legend = plt.legend(title=r'$\bf copy \ state$', title_fontsize='large', handles=copy_legend_elements,
                            markerscale=1, fontsize=13, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        ax.add_artist(legend)

        target_legend_elements = [
            Line2D([0], [0], label='stable', marker=markers['stable'],
                   markersize=8, color='#eaeaf2', markerfacecolor='dimgray'),
            Line2D([0], [0], label='poor mappability', marker=markers['poor mappability'],
                   markersize=8, color='#eaeaf2', markerfacecolor='dimgray'),
            Line2D([0], [0], label='high variability', marker=markers['high variability'],
                   markersize=8, color='#eaeaf2', markerfacecolor='dimgray'),
            Line2D([0], [0], label='poor map & high var', marker=markers['poor map & high var'],
                   markersize=8, color='#eaeaf2', markerfacecolor='dimgray'),
            Line2D([0], [0], label='low baseline read-depth', marker=markers['low baseline read-depth'],
                   markersize=8, color='#eaeaf2', markerfacecolor='dimgray')
        ]
        markers_legend = plt.legend(title=r'$\bf target\ quality$', title_fontsize='large',
                                    handles=target_legend_elements,
                                    markerscale=1, fontsize=13, bbox_to_anchor=(1.05, 0.5), loc=2, borderaxespad=0.)
        ax.add_artist(markers_legend)

        plt.savefig(f'{output_path}/{selected_gene}.png', dpi=300, bbox_inches='tight')


def parse_args(args) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Filter targets file by genes')
    parser.add_argument('-i', '--targets', type=str, required=True, help='input targets file (nrrtest.tsv)')
    parser.add_argument('-m', '--targets-map', type=str, required=True, help='input targets mappability file (BED)')
    parser.add_argument('-o', '--output-dir', type=str, required=True, help='output dir')
    parser.add_argument('-g', '--gene-list', type=str, required=True, help='gene list file')
    return parser.parse_args(args)


if __name__ == '__main__':
    parsed_args = parse_args(sys.argv[1:])
    gene_list = pd.read_csv(parsed_args.gene_list, header=None, names=['symbol'])
    table_index = ['chrom', 'chromStart', 'chromEnd']
    mappability_table = pd.read_table(parsed_args.targets_map, names=table_index + ['mappability'])
    target_table = pd.read_table(parsed_args.targets).sort_values(by=table_index)
    target_table['mad/median'] = (target_table['mad'] / target_table['norm_median']).round(2)

    # compute Q1, Q3, IQR for ratios
    ratios = target_table[(target_table['median'] >= 25) &
                          (target_table['ratio'] > 0.8) &
                          (target_table['ratio'] < 1.2)]['ratio']
    Q1 = ratios.quantile(0.25)
    Q3 = ratios.quantile(0.75)
    IQR = Q3 - Q1

    for gene in gene_list['symbol']:
        plot_gene(gene, target_table, mappability_table, Q1, Q3, IQR, parsed_args.output_dir)
