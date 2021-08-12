rule filter_cnv_targets:
    input:
        targets=rules.annotate_map.output,
        gene_list=config['gene-filter']['list']
    output:
        OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +
        "filtered/{sample}_cnv-targets-filtered-by-genes.tsv"
    shell:
        "python scripts/filter_targets_by_genes.py -i {input.targets} -g {input.gene_list} -o {output}"


rule filter_cnvs:
    input:
        cnvs=rules.annotate_stats.output,
        gene_list=config['gene-filter']['list']
    output:
        OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +
        "filtered/{sample}_cnvs-filtered-by-genes.tsv"
    shell:
        "python scripts/filter_cnvs_by_genes.py -i {input.cnvs} -g {input.gene_list} -o {output}"


rule filter_all_targets:
    input:
        targets=rules.detect_cnv.output,
        gene_list=config['gene-filter']['list']
    output:
        OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +
        "filtered/{sample}_all-targets-filtered-by-genes.tsv"
    shell:
        "python scripts/filter_targets_by_genes.py -i {input.targets} -g {input.gene_list} -o {output}"






