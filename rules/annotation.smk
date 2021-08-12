INTERSECT_FRACTION = 0.99 if config['filter']['skip'] else 1
header = 'chrom\tchromStart\tchromEnd\tgene\t' \
         '# reads (test)\t# mean (baseline)\t# median (baseline)\t' \
         '# reads norm (test)\t# mean norm (baseline)\t' \
         'std\t% mean (std)\t' \
         '# median norm (baseline)\tiqr\tmad\t% median (mad)\t' \
         'cnv id\tratio\tmappability'

rule annotate_map:
    input:
        exome_mappability=config['annotation']['mappability'],
        cnv_bed=rules.detect_cnv.output
    output:
        OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +
        "{sample}_cnv-targets.tsv"
    params:
        bed_input=get_sample_cov,
        slop="-l 1 -r 0",
        intersect=f'-f {INTERSECT_FRACTION} -r -wa -wb',
        header=header
    shell:
        "bedtools intersect -a {params.bed_input}.bed -b {input.exome_mappability} {params.intersect} | " +
        "awk -v FS='\t' -v OFS='\t' '{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $21}}' > {output} && " +
        "sed -i '1i {params.header}' {output}"


rule annotate_stats:
    input:
        rules.annotate_map.output
    output:
        OUTPUT_DIR + "cnv-analysis/{sample}" + f'{FILTER_OPT}' + f'-{RD_WINDOW_SIZE}-{RD_WINDOW_STEP}/annotated/' +
        "{sample}_cnv-stats.tsv"
    shell:
        "python scripts/cnv_stats.py -i {input} -o {output}"
