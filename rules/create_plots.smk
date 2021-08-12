rule plots_by_genes:
    input:
        targets=rules.detect_cnv.output,
        exome_mappability=config['annotation']['mappability'],
        gene_list=config['gene-filter']['list']
    output:
        directory(OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' +
                  f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/plots/genes')
    shell:
        "mkdir {output} && python scripts/create_plots_by_genes.py \
        --targets {input.targets} \
        --targets-map  {input.exome_mappability} \
        --gene-list {input.gene_list} \
        --output-dir {output}"
