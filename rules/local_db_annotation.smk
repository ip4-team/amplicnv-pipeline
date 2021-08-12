rule localdb_annotation:
    input:
        test=rules.class_annotation.output,
        targets=rules.annotate_map.output
    output:
        cnvs=OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +
             "{sample}_cnvs-annotated.tsv",
        genes=OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +
              "{sample}_non-overlapping-genes.tsv"
    params:
        baseline_regex=config['local-db-annotation']['baseline-path-regex'],
        min_overlap_fraction=config['local-db-annotation']['min-overlap-fraction']
    shell:
        "python scripts/local_db_annotation.py -b '{params.baseline_regex}' -t {input.test} " +
        "-c {input.targets} -a {output.cnvs} -g {output.genes} -f {params.min_overlap_fraction}"