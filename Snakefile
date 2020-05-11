#!/usr/bin/env python
import glob

include: "rules/align.smk"
include: "rules/sv_callers.smk"
include: "rules/qc.smk"


# Deal with the lsf logfile
import os
if not os.path.exists("loglsf"):
    os.makedirs("loglsf")

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

#target rule
rule all:
    input:
        expand(OUTDIR + "/{ref}/{aligner}/sv_calls/svim_{sample}/variants.vcf",ref=REF_GENOMES,aligner=ALIGNERS,sample=sample),
        expand(OUTDIR+"/{ref}/{aligner}/sv_calls/{sample}_pbsv.vcf",ref=REF_GENOMES, aligner=ALIGNERS, sample=sample),
        expand(OUTDIR+"/{ref}/{aligner}/sv_calls/{sample}_sniffles.vcf",ref=REF_GENOMES, aligner=ALIGNERS, sample=sample),
        expand(OUTDIR+"/{ref}/{aligner}/sv_calls/{sample}_svim.vcf",ref=REF_GENOMES, aligner=ALIGNERS, sample=sample),
        expand(OUTDIR+ "/{ref}/{aligner}/sv_calls/{sample}_{svcaller}_filtered.vcf",ref=REF_GENOMES, aligner=ALIGNERS,sample=sample,svcaller=SVCALLER),
        expand(OUTDIR + "/{ref}/sv_combined/high_conf_{ref}.vcf",ref=REF_GENOMES),
        expand(OUTDIR+ "/{ref}/{aligner}/alignment/{sample}_bam.stats",ref=REF_GENOMES, aligner=ALIGNERS, sample=sample),
        expand(OUTDIR+ "/{ref}/{aligner}/qc_{sample}_bam",ref=REF_GENOMES, aligner=ALIGNERS,sample=sample),
        expand(OUTDIR+ "/{ref}/{aligner}/alignment/{sample}_coverage_bam.png",ref=REF_GENOMES, aligner=ALIGNERS,sample=sample)
