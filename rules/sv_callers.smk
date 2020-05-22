checkpoint sniffles:
    input:
        OUTDIR + "/{ref}/{aligner}/alignment/{sample}_merge.bam"
    output:
        OUTDIR + "/{ref}/{aligner}/sv_calls/{sample}_sniffles.vcf"
    params:
        sniffles= config["tools"]["SNIFFLES"],
        min_mq = 20,
        read_support = 5,
        min_length = 50
    envmodules:
        "intel/18.0.1"
    shell:
        """
        {params.sniffles} -m {input} -v {output} -s {params.read_support} -q {params.min_mq} -l {params.min_length} --genotype --report_read_strands
        """

checkpoint pbsv_signature:
    input:
        OUTDIR + "/{ref}/{aligner}/alignment/{sample}_merge.bam"
    output:
        temp(OUTDIR + "/{ref}/{aligner}/sv_calls/{sample}_pbsv.svsig.gz")
    conda:
        "../envs/pacbio.yml"
    params:
        map_Q = 20
    shell:
        """
        pbsv discover -q {params.map_Q} {input} {output}
        
        """

rule pbsv_call_sv:
    input:
        svsig = OUTDIR + "/{ref}/{aligner}/sv_calls/{sample}_pbsv.svsig.gz",
        ref = GENOMEDIR + "/{ref}_ref.fa"
    output:
        OUTDIR + "/{ref}/{aligner}/sv_calls/{sample}_pbsv.vcf"
    conda:
        "../envs/pacbio.yml"
    threads:
        10
    shell:
        """
        pbsv call -j {threads} {input.ref} {input.svsig} {output}

        """

checkpoint svim_call:
    input:
        BAM = OUTDIR + "/{ref}/{aligner}/alignment/{sample}_merge.bam",
        ref = GENOMEDIR + "/{ref}_ref.fa"
    output:
        OUTDIR + "/{ref}/{aligner}/sv_calls/svim_{sample}/variants.vcf"
    params:
    conda:
        "../envs/svim.yml"
    shell:
        """
        svim alignment --sample {wildcards.sample} --min_mapq 20 --min_sv_size 50 \
         {wildcards.ref}/{wildcards.aligner}/sv_calls/svim_{wildcards.sample}/ {input.BAM} {input.ref}
        """

rule svim_filter:
    input:
        OUTDIR + "/{ref}/{aligner}/sv_calls/svim_{sample}/variants.vcf"
    output:
        vcf = OUTDIR + "/{ref}/{aligner}/sv_calls/{sample}_svim.vcf"
    params:
        filtered = OUTDIR + "/{ref}/{aligner}/sv_calls/svim_{sample}/variants_filtered.vcf.gz"
    shell:
        """
        cat {input}  \
        | sed 's/INS:NOVEL/INS/g' \
        | sed 's/DUP:INT/DUP/g' \
        | sed 's/DUP:TANDEM/DUP/g' \
        | awk '{{ if($1 ~ /^#/) {{ print $0 }} \
         else {{ if($6>5) {{ print $0 }} }} }}' > {output.vcf}

         #cat {input}  \
         #| awk '{{ if($1 ~ /^#/) {{ print $0 }} \
         #else {{ if($6>5) {{ print $0 }} }} }}' \
         #| bgzip > {params.filtered}

         #tabix {params.filtered}

         """

rule filter_vcf:
    input:
        VCF = OUTDIR + "/{ref}/{aligner}/sv_calls/{sample}_{svcaller}.vcf",
        BED = GENOMEDIR + "/{ref}_chr.bed"
    output:
        OUTDIR + "/{ref}/{aligner}/sv_calls/{sample}_{svcaller}_filtered.vcf"
    params:
        bcftools = config["tools"]["BCFTOOLS"]
    shell:
        """
        {params.bcftools} view -T {input.BED} {input.VCF} \
                | {params.bcftools} view -f PASS \
                | {params.bcftools} sort > {output}

        """

rule merge_sv_vcf:
    input:
        expand(OUTDIR + "/{{ref}}/{aligner}/sv_calls/{sample}_{svcaller}_filtered.vcf", aligner=ALIGNERS, sample=sample, svcaller=SVCALLER)
    output:
        vcf= temp(OUTDIR + "/{ref}/sv_combined/tmp.vcf"),
        fofn = OUTDIR + "/{ref}/sv_combined/samples.fofn",
        slist= OUTDIR + "/{ref}/sv_combined/samples_name.txt"
    params:
        survivor = config["tools"]["SURVIVOR"]
    shell:
        """
        ls {input} > {output.fofn} ;
        awk -F '/|_' '{{print $14"_"$19}}' {output.fofn} > {output.slist};
        {params.survivor} merge {output.fofn} 1000 1 1 1 0 50 {output.vcf}
        """

rule high_conf:
    input:
        vcf= OUTDIR + "/{ref}/sv_combined/tmp.vcf",
        slist= OUTDIR + "/{ref}/sv_combined/samples_name.txt"
    output:
        rehead_vcf= OUTDIR + "/{ref}/sv_combined/sv_merge.vcf",
        high_conf= OUTDIR + "/{ref}/sv_combined/high_conf_{ref}.vcf"
    params:
        bcftools = config["tools"]["BCFTOOLS"]
    shell:
        """
        {params.bcftools} reheader -s {input.slist} {input.vcf} -o {output.rehead_vcf} &&
        {params.bcftools} filter -i 'INFO/SUPP = "6"' -o {output.high_conf} {output.rehead_vcf}
        """
