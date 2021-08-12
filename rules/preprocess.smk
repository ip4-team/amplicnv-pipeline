rule strip_header:
    input:
        get_cov
    output:
        OUTPUT_DIR + "preprocessed/{sample}_headerless.tsv"
    shell:
        "awk 'NR > 1 {{print}}' < {input} > {output}"


rule sort:
    input:
        rules.strip_header.output
    output:
        OUTPUT_DIR + "preprocessed/{sample}_headerless_sorted.tsv"
    shell:
        "sort -k 1,1 -k2,2n {input} > {output}"


rule slop:
    input:
        sizes=config["genome"]["index"],
        tsv=rules.sort.output
    output:
        OUTPUT_DIR + "preprocessed/{sample}_hless_sorted_0based.tsv"
    params:
        "-l 1 -r 0"
    shell:
        "bedtools slop -i {input.tsv} -g {input.sizes} {params} > {output}"


rule subtract:
    input:
        bed=config['filter']['mappability'],
        tsv=rules.slop.output
    output:
        OUTPUT_DIR + "preprocessed/{sample}_filtered.tsv"
    params:
        "-f 1 -r"
    shell:
        "bedtools subtract -a {input.tsv} -b {input.bed} {params} > {output}"

