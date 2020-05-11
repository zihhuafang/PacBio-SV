rule aln_stats:
    input:
        OUTDIR + "/{ref}/{aligner}/alignment/{sample}_merge.bam"
    output:
        OUTDIR + "/{ref}/{aligner}/alignment/{sample}_bam.stats"
    envmodules:
        "gcc/4.8.5",
        "python_cpu/3.6.4"
    shell:
        """
        export PATH=$PATH:/cluster/work/pausch/group_bin/anaconda3/envs/surpyvor
        ../scripts/alignment_stats.py {input} --output {output}

        """

rule nanoplot_qc:
    input:
         OUTDIR + "/{ref}/{aligner}/alignment/{sample}_merge.bam"
    output:
         DIR = directory(OUTDIR + "/{ref}/{aligner}/qc_{sample}_bam")
    params:
         sample = config['sample_name']
    threads:
         10
    conda:
        "../envs/NanoPack.yml"   
    shell:
        "NanoPlot -t {threads} --bam {input} --raw -o {output.DIR} -p {params.sample}_{wildcards.aligner} --N50 --title {params.sample}_{wildcards.aligner} "


rule plot_coverage:
    input:
        OUTDIR + "/{ref}/{aligner}/alignment/{sample}_merge.bam"
    output:
        OUTDIR + "/{ref}/{aligner}/alignment/{sample}_coverage_bam.png"
    conda:
        "../envs/deeptools.yml"
    shell:
        """
        plotCoverage -b {input} -o {output} --minMappingQuality 20
        """
