rule vep_input:
    input:
        rules.annotate_stats.output
    output:
        OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +
        "{sample}_vep-input.tsv"
    shell:
        "python scripts/create_vep_input.py -i {input} -o {output}"


rule phenogramviz_input:
    input:
        rules.annotate_stats.output
    output:
        all=OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +
            "{sample}_phenogramviz-input.txt",
        filtered=OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +
                 "{sample}_phenogramviz-input-filtered.txt"
    shell:
        "python scripts/create_phenogramviz_input.py -i {input} -o {output.all} -f {output.filtered}"


rule classify_cnv_input:
    input:
        rules.annotate_stats.output
    output:
        all=OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +
            "{sample}_classifyCNV-input.txt",
        filtered=OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +
                 "{sample}_classifyCNV-input-filtered.txt"
    shell:
        "python scripts/create_classifyCNV_input.py -i {input} -o {output.all} -f {output.filtered}"
