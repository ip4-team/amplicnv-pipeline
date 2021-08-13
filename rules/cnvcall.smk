rule create_cfg:
    input:
        bed_file=config['exome']['bed'],
        test_sample=get_test,
        baseline_samples=get_baseline
    output:
        OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}' + ".cfg"
    run:
        cfg_creator.create_cfg(
            input.bed_file,
            build_df(test_samples, config['test']['dir']),
            build_df(baseline_samples[~baseline_samples.index.isin([wildcards.sample])],
                config['baseline']['dir'], sex=test_samples.loc[wildcards.sample, 'sex']),
            wildcards.sample,
            RD_WINDOW_SIZE, RD_WINDOW_STEP,
            FILTER_OPT, USE_VARIANTS,
            OUTPUT_DIR + 'cnv-analysis/')


rule detect_cnv:
    input:
        rules.create_cfg.output
    output:
        OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/tables/bam/nrrtest.tsv'
    log:
        OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}' + ".log"
    shell:
        "amplicnv detectcfg --cfg-file {input} &> {log}"

