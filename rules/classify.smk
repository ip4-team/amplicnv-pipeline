CLASSIFY_CNV_VERSION = "1.1.0"
GENOME_BUILD = config['genome']['build']

rule install_classify_cnv:
    output:
        directory("ClassifyCNV-" + CLASSIFY_CNV_VERSION)
    shell:
        "bash scripts/download_classifyCNV.sh " + CLASSIFY_CNV_VERSION


rule classify_cnv:
    input:
        tool_dir=rules.install_classify_cnv.output,
        cnvs=OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +
             "{sample}_classifyCNV-input.txt"
    output:
        directory(OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/classified/' +
        "{sample}_ClassifyCNV-output")
    shell:
        "python ClassifyCNV-{CLASSIFY_CNV_VERSION}/ClassifyCNV.py --infile {input.cnvs} --GenomeBuild {GENOME_BUILD} " +
        "--outdir {output}"


rule class_annotation:
    input:
        cnvs=rules.annotate_stats.output,
        class_dir=rules.classify_cnv.output,
        gene_list=config['gene-filter']['list']
    output:
        OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +
        "{sample}_cnvs-classified.tsv"
    shell:
        "python scripts/add_classification.py -s {input.cnvs} -c {input.class_dir}/Scoresheet.txt -g {input.gene_list} -o {output}"

