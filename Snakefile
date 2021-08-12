from snakemake.utils import min_version
from utils import cfg_creator


# set minimum snakemake version
min_version('5.3.0')


# load configuration file
configfile: 'config-example.yaml'

# setup
OUTPUT_DIR = f'{config["output"]["dir"]}/'.replace('//', '/')

wildcard_constraints:
    sample="EMPTY-MOCK-\d+"


FILTER_OPT = '' if config['filter']['skip'] else '_filtered'
RD_WINDOW_SIZE = config['detection']['read-depth']['window-size']
RD_WINDOW_STEP = config['detection']['read-depth']['window-step']
USE_VARIANTS = config['detection']['variants']['use']
LOCAL_ANNOTATION = not config['local-db-annotation']['skip']


include: 'rules/common.smk'

# load test samples
test_samples = load_samples_table(load_sample_info(config['test']['samples']), config['test']['dir'])
print(f'Test samples loaded: {len(test_samples)}')

# load baseline samples
baseline_samples = load_samples_table(load_sample_info(config['baseline']['samples']), config['baseline']['dir'])
print(f'Baseline samples loaded: {len(baseline_samples)}')

include: 'rules/preprocess.smk'
include: 'rules/cnvcall.smk'
include: 'rules/annotation.smk'
include: 'rules/create_inputs.smk'
include: 'rules/filter.smk'
include: 'rules/classify.smk'
include: 'rules/local_db_annotation.smk'
include: 'rules/create_spreadsheet.smk'
include: 'rules/create_plots.smk'

rule all:
    input:
        cnv_stats=expand(rules.annotate_stats.output, sample=get_test_sample_ids()),
        vep_input=expand(rules.vep_input.output, sample=get_test_sample_ids()),
        phenogramviz_input=expand(rules.phenogramviz_input.output, sample=get_test_sample_ids()),
        cnvs_filtered_by_gene=expand(rules.filter_cnvs.output, sample=get_test_sample_ids()),
        cnv_targets_filtered_by_gene=expand(rules.filter_cnv_targets.output, sample=get_test_sample_ids()),
        all_targets_filtered_by_gene=expand(rules.filter_all_targets.output, sample=get_test_sample_ids()),
        classifed_cnvs=expand(rules.class_annotation.output, sample=get_test_sample_ids()),
        annotated_cnvs=expand(rules.localdb_annotation.output.cnvs, sample=get_test_sample_ids()) if LOCAL_ANNOTATION else [],
        spreadsheet=expand(rules.spreadsheet.output, sample=get_test_sample_ids()),
        plots=expand(rules.plots_by_genes.output, sample=get_test_sample_ids())
