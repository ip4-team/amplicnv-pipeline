import pandas as pd
import glob
import re


def get_directory(path: str):
    return f'{path}/'.replace('//', '/')


def load_sample_info(tsv_file):
    return pd.read_table(tsv_file, usecols=[0, 1], names=['id', 'sex'], header=None).set_index('id')


def get_test_sample_ids():
    return list(test_samples.index)


def load_samples_table(info, path):
    directory = get_directory(path)
    cov_files = glob.glob(f'{directory}*.cov') + glob.glob(f'{directory}*.xls')
    var_files = glob.glob(f'{directory}*.vcf.gz') if USE_VARIANTS else []

    info.loc[:, 'cov_file'] = info.apply(
        lambda row: list(filter(re.compile(f'{row.name}').findall, cov_files))[0].split('/')[-1],
        axis=1, result_type='expand')

    if USE_VARIANTS:
        info.loc[:, 'var_file'] = info.apply(
            lambda row: list(filter(re.compile(f'{row.name}').findall, var_files))[0].split('/')[-1],
            axis=1, result_type='expand')

    return info


def get_cov(wildcards):
    test_dir = get_directory(config['test']['dir'])
    baseline_dir = get_directory(config['baseline']['dir'])
    if wildcards.sample in test_samples.index:
        return test_dir + test_samples.loc[wildcards.sample, 'cov_file']
    else:
        return baseline_dir + baseline_samples.loc[wildcards.sample, 'cov_file']


def get_baseline(wildcards):
    baseline_dir = get_directory(config['baseline']['dir'])
    baseline_df = build_df(
        baseline_samples[baseline_samples['sex'] == test_samples.loc[wildcards.sample, 'sex']],
        baseline_dir)
    baseline_df = baseline_df[~baseline_df.index.isin([wildcards.sample])]
    if USE_VARIANTS:
        return baseline_df.loc[:, ['cov_file', 'var_file']].values.flatten().tolist()
    else:
        return baseline_df.loc[:, ['cov_file']].values.flatten().tolist()


def get_test(wildcards):
    test_dir = get_directory(config['test']['dir'])
    test_df = build_df(test_samples, test_dir)
    if USE_VARIANTS:
        return test_df.loc[wildcards.sample, ['cov_file', 'var_file']].tolist()
    else:
        return test_df.loc[wildcards.sample, ['cov_file']].tolist()


def build_df(data, data_dir, sex = None):
    directory = get_directory(data_dir)

    if sex:
        data = data[data['sex'] == sex]

    if USE_VARIANTS:
        if config['filter']['skip']:
            cov_full = data['cov_file'].apply(lambda f: directory + f)
            vcf_full = data['var_file'].apply(lambda f: directory + f)
            return pd.DataFrame({'cov_file': cov_full, 'var_file': vcf_full})
        else:
            cov_full = map(lambda index: f'{OUTPUT_DIR}preprocessed/{index}_filtered.tsv',
                list(data.index))
            vcf_full = data['var_file'].apply(lambda f: directory + f)
            return pd.DataFrame({'cov_file': cov_full, 'var_file': vcf_full})
    else:
        if config['filter']['skip']:
            cov_full = data['cov_file'].apply(lambda f: directory + f)
            return pd.DataFrame({'cov_file': cov_full})
        else:
            cov_full = map(lambda index: f'{OUTPUT_DIR}preprocessed/{index}_filtered.tsv',
                list(data.index))
            return pd.DataFrame({'cov_file': cov_full})




def get_sample_cov(wildcards):
    skip_map = config['filter']['skip']
    filename = test_samples.loc[wildcards.sample, 'cov_file'] if skip_map else f'{wildcards.sample}_filtered.tsv.bed'
    return OUTPUT_DIR + f'cnv-analysis/{wildcards.sample}{FILTER_OPT}-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/bed/{filename}'


