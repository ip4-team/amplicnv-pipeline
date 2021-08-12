import pandas as pd


def create_cfg(bed_file: str, sample_df: pd.DataFrame, baseline_df: pd.DataFrame,
               test_sample_id: str, rd_size: int, rd_step: int, filter_opt: str,
               use_vars: bool, output_dir: str):

    output_file = f'{test_sample_id}{filter_opt}-{rd_size}-{rd_step}.cfg'
    with open(output_dir + output_file, 'w') as file:
        file.write('[bed]\n')
        file.write(f'\tbedfile = {bed_file}\n')

        file.write('[sample]\n')
        file.write(f'\tcovfile = {sample_df.loc[test_sample_id, "cov_file"]}\n')
        if use_vars:
            file.write(f'\tvcffile = {sample_df.loc[test_sample_id, "var_file"]}\n')

        file.write('[baseline]\n')
        baseline_cov = '\n\t\t'.join(baseline_df['cov_file'].tolist())
        file.write(f'\tcovfiles =\n\t\t{baseline_cov}')
        if use_vars:
            baseline_vcf = '\n\t\t'.join(baseline_df['var_file'].tolist())
            file.write(f'\n\tvcffiles =\n\t\t{baseline_vcf}')

        file.write('\n[targtest]\n')
        params = '\n\t'.join([
            f'size = {rd_size}',
            f'step = {rd_step}',
            f'metric = IQR',
            f'interval_range = 1.5'
        ])
        file.write(f'\t{params}')

        file.write('\n[output]\n')
        file.write(f'\tpath = {output_dir}/{output_file.strip(".cfg")}\n'.replace('//', '/'))
