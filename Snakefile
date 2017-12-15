import os

workdir: '/mnt/large/treetesting'

LOG = 'vedro-pipeline.log'

VEDRO = '/common/contrib/bin/vedro-0.2.0'
VEDRO_TREE_BUILD = '/common/contrib/bin/vedro-tree-build-0.2.0'
VEDRO_BUILD = '/common/contrib/bin/vedro-build-0.2.0'
VEDRO_COLLAPSE = '/common/contrib/bin/vedro-collapse-0.2.0'

REFERENCE_FASTA = 'reference/all_seqs.fasta'
NCBI_TAXDUMP = 'reference/taxdump.tar.gz'
CHUNKS = list(range(0, 5))

INPUTS = [
'bact_SRR1270040.fastq.gz',
'virs_ERR862184.fastq.gz',
'bact_SRR1284500.fastq.gz',
'virs_ERR885561.fastq.gz']

SAMPLES = [os.path.join('samples', f) for f in INPUTS]

rule report:
    log: LOG
    input:
        "bowtie2/alignments_of_interest.sam",
        "informative_hits.txt"
    output:
        "report.pdf"
    shell:
        "touch {output[0]}"

rule align_informatives_against_genomes_of_interest:
    log: LOG
    input:
        "intermediate/informative_reads.fasta",
        "bowtie2/of_interest.bt2"
    output:
        "bowtie2/alignments_of_interest.sam"
    threads: 1
    shell:
        "bowtie2 -f --very-sensitive --end-to-end --all --no-unal --no-sq --no-head"
        "--omit-sec-seq --threads {threads} -x {input[1]} -U {input[0]} -S {output[0]}"

rule extract_informative_reads:
    log: LOG
    input:
        "informative_hits.txt",
        "intermediate/prepared_reads.fasta"
    output:
        "intermediate/informative_reads.fasta"
    shell:
        "touch {output}"

rule build_aligner_index_of_genomes_of_interest:
    log: LOG
    input:
        "bowtie2/genomes_of_interest.fasta"
    output:
        "bowtie2/of_interest.bt2"
    shell:
        "bowtie2-build -f {input} {output}"

rule extract_genomes:
    log: LOG
    input:
        "informative_hits.txt",
        expand("reference/{chunkid}.fasta", chunkid=CHUNKS)
    output: "bowtie2/genomes_of_interest.fasta"
    shell:
        "cat {input} > {output[0]}"

rule identify_informatives:
    log: LOG
    input:
        "intermediate/collapsed_vedro_results.txt",
        "indexes/index.tax"
    output:
        "informative_hits.txt"
    shell:
        "echo \"{input}\" > {output[0]}"

rule collapse_vedro:
    log: LOG
    input: expand("intermediate/{chunkid}vedresults.txt", chunkid=CHUNKS)
    output: "intermediate/collapsed_vedro_results.txt"
    shell: VEDRO_COLLAPSE + "-v -o {output} {input} 2> {log}"

rule bin_reads:
    log: LOG
    input: reads="intermediate/prepared_reads.fasta", index="indexes/chunk.{chunkid}.index"
    output: "intermediate/{chunkid}vedresults.txt"
    threads: 16
    shell:
        VEDRO + '-v -e 3 -f {input.reads} -i {input.index} -t {threads} -r {output} 2> {log}'

rule qc_reads:
    log: LOG
    input: expand("{sample}", sample=SAMPLES)
    output: "intermediate/prepared_reads.fasta"
    threads: 1
    shell:
        "gunzip -cd {input} | "
        "fastq_quality_filter -q 25 -p 80 | "
        "fastx_trimmer -f 5 -l 55 | "
        "fastq_to_fasta -n -r | "
        "fastx_collapser -o {output} 2> {log}"

rule build_index:
    log: LOG
    input: fasta="reference/{chunkid}.fasta"
    output: "indexes/chunk.{chunkid}.index"
    threads: 1
    shell: VEDRO_BUILD + ' -v -f {input.fasta} -i {output} 2> {log}'

rule build_taxonomy_tree:
    log: LOG
    output: "indexes/index.tax"
    input: ref=REFERENCE_FASTA, dump=NCBI_TAXDUMP
    shell: VEDRO_TREE_BUILD + "-v -f {input.ref} -i {output} -d {input.dump} 2> {log}"

rule reference_chunks:
    log: LOG
    output: expand("reference/{chunkid}.fasta", chunkid=CHUNKS)
    shell: "true"
