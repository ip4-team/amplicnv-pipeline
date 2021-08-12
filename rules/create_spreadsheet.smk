def get_cnv_file(wildcards):
    if LOCAL_ANNOTATION:
        return OUTPUT_DIR + f'cnv-analysis/{wildcards.sample}' +\
               f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +\
               f'{wildcards.sample}_cnvs-annotated.tsv'
    else:
        return OUTPUT_DIR + f'cnv-analysis/{wildcards.sample}' +\
               f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +\
               f'{wildcards.sample}_cnvs-classified.tsv'


rule spreadsheet:
    input:
        cnvs=get_cnv_file,
        targets=OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' +
                f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +
                "{sample}_cnv-targets.tsv",
        cnvs_filtered=OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' +
                      f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +
                      "filtered/{sample}_cnvs-filtered-by-genes.tsv",
        targets_filtered=OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' +
                         f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +
                         "filtered/{sample}_cnv-targets-filtered-by-genes.tsv"
    output:
        OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/' +
        "{sample}.xlsx"
    shell:
        "python scripts/create_spreadsheet.py -c {input.cnvs} -t {input.targets} " +
        "-f {input.cnvs_filtered} -f {input.targets_filtered} " +
        "-o {output}"
