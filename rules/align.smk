#!/usr/bin/env python

configfile: "config.yaml"
FASTQDIR = config["resources"]["fastqdir"]
GENOMEDIR = config["resources"]["genomedir"]
## Parsing output directory
OUTDIR = config["resources"]["outdir"]

REF_GENOMES, = glob_wildcards(GENOMEDIR + "/{ref}_ref.fa")
ALIGNERS=["ngmlr", "minimap2"]
MOVIES, = glob_wildcards(FASTQDIR + "/{movie}.fastq.gz")
SVCALLER=["sniffles", "pbsv", "svim"]


# Parameter: sample_name
sample = "sv_sample01"
if "sample_name" in config:
    sample = config['sample_name']


checkpoint ngmlr_aln:
    input:
        fastq = FASTQDIR + "/{movie}.fastq.gz",
        ref = GENOMEDIR + "/{ref}_ref.fa"
    output:
        OUTDIR + "/{ref}/ngmlr/alignment/aln_{movie}.bam"
    params:
        ngmlr = config["tools"]["NGMLR"],
        sample = config['sample_name']
    threads:
        20
    envmodules:
        "gcc/4.8.2",
        "samtools/1.8"
    shell:
        """
        {params.ngmlr} -t {threads} -r {input.ref} -q {input.fastq} --rg-id {wildcards.movie} --rg-sm {params.sample} |
 \
                samtools sort -@ {threads} -T /cluster/work/pausch/temp_scratch/fang/ -o {output} -
        """

checkpoint pbmm2_aln:
    input:
        fastq = FASTQDIR + "/{movie}.fastq.gz",
        ref = GENOMEDIR + "/{ref}_ref.fa"
    output:
        temp(OUTDIR + "/{ref}/minimap2/alignment/aln_{movie}.bam")
    threads:
        10
    params:
        sample = config['sample_name']
    conda:
        "../envs/pacbio.yml"
    shell:
        """

        pbmm2 align {input.ref} {input.fastq} --rg '@RG\tID:{wildcards.movie}\tSM:{params.sample}' -j {threads} | \
                samtools calmd -bS - {input.ref} | samtools sort -@ {threads} -T $TMPDIR > {output}
        """

rule merge_bam:
    input:
        lambda wildcards: expand(OUTDIR + "/{{ref}}/{{aligner}}/alignment/aln_{movie}.bam",movie=MOVIES)
    output:
        BAM = protected(OUTDIR + "/{ref}/{aligner}/alignment/{sample}_merge.bam"),
        BAI = protected(OUTDIR + "/{ref}/{aligner}/alignment/{sample}_merge.bam.bai")
    output:
        BAM = protected(OUTDIR + "/{ref}/{aligner}/alignment/{sample}_merge.bam"),
        BAI = protected(OUTDIR + "/{ref}/{aligner}/alignment/{sample}_merge.bam.bai")
    envmodules:
        "gcc/4.8.2",
        "samtools/1.8"
    threads: 10
    shell:
        """
        samtools merge -@ {threads} - {input} | \
                samtools sort -@ {threads} -T $TMPDIR > {output.BAM} && samtools index -@ {threads} {output.BAM}

        """
