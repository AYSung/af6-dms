SAMPLES = [
    '1_ko1a_p0',
    '2_ko1a_p3',
    '3_ko1a_p6',
    '4_ko1b_p0',
    '5_ko1b_p3',
    '6_ko1b_p6',
    '7_ko2a_p0',
    '8_ko2a_p3',
    '9_ko2a_p6',
    '10_ko2b_p0',
    '11_ko2b_p3',
    '12_ko2b_p6',
    '13_ko2c_p0',
    '14_ko2c_p3',
    '15_ko2c_p6',
]

rule all:
    input:
        expand("variant-counts/{sample}.counts.txt", sample=SAMPLES)

rule filter_aligned:
    input:
        r1="raw-data/{sample_prefix}.R1.fastq.gz",
        r2="raw-data/{sample_prefix}.R2.fastq.gz"
    threads:
        workflow.cores
    params:
        "-x bowtie2-indices/NDUFAF6_endo/ndufaf6_endo_cdna_reference -t"
    output:
        "{output_dir}/{sample_prefix}.R1.fastq.gz",
        "{output_dir}/{sample_prefix}.R2.fastq.gz"
    log:
        "{output_dir}/{sample_prefix}.summary.txt"
    shell:
        "bowtie2 -p {threads} {params} -1 {input.r1} -2 {input.r2} --al-conc-gz {wildcards.output_dir}/{wildcards.sample_prefix}.R%.fastq.gz 1> /dev/null 2> {log}"

rule merge_pairs:
    input:
        "filtered-data/{sample_prefix}.R1.fastq.gz",
        "filtered-data/{sample_prefix}.R2.fastq.gz"
    params:
        "-M 150 -m 20 -z -o {sample_prefix} -d {output_dir}"
    threads:
        workflow.cores
    output:
        "{output_dir}/{sample_prefix}.extendedFrags.fastq.gz"
    log:
        "{output_dir}/{sample_prefix}.flash.txt"
    shell:
        "flash -t {threads} {input} {params} 2>&1 | tee {log}"

rule filter_wt:
    input:
        "merged-data/{sample_prefix}.extendedFrags.fastq.gz"
    output:
        "variant-data/{sample_prefix}.variants.fastq"
    shell:
        "python filter_variants.py {input}"


rule gzip_variants:
    input:
        "variant-data/{sample_prefix}.variants.fastq"
    output:
        "variant-data/{sample_prefix}.variants.fastq.gz"
    shell:
        "gzip {input}"


rule orient_reads:
    input:
        "variant-data/{sample_prefix}.variants.fastq.gz"
    output:
        "oriented-data/{sample_prefix}.oriented.fasta"
    shell:
        "python orient_variants.py {input}"

rule gzip_oriented_reads:
    input:
        "oriented-data/{sample_prefix}.oriented.fasta"
    output:
        "oriented-data/{sample_prefix}.oriented.fasta.gz"
    shell:
        "gzip {input}"

rule count_variants:
    input:
        "oriented-data/{sample_prefix}.oriented.fasta.gz"
    output:
        "variant-counts/{sample_prefix}.counts.txt"
    shell:
        "python count_variants.py {input}"