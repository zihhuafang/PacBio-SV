#!/bin/env bash
cache_dir='/cluster/work/pausch/fang/vep_99'

#ref='UCD'
ref='Angus'

cd ${cache_dir}

#copy file to cache directory


#singularity exec -H ${cache_dir}:/cluster/home/fangzi \
#	-B /cluster/work/pausch/fang/vep_99:/opt/vep/.vep \
#	ensembl-vep.simg vep \
#	--species bos_taurus \
#        --cache --dir /opt/vep/.vep \
#	--format vcf \
#	--offline --vcf --force_overwrite \
#	--overlaps --sift b --max_sv_size 30343437 \
#	-fasta ./Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz \
#	--hgvs --symbol \
#	--input_file high_conf_${ref}.vcf \
 #       --output_file high_conf_${ref}_vep.vcf \
#	--stats_text high_conf_${ref}_vep.stats

singularity exec -H ${cache_dir}:/cluster/home/fangzi \
        -B /cluster/work/pausch/fang/vep_99:/opt/vep/.vep \
        ensembl-vep.simg vep \
        --species bos_taurus_hybrid \
        --cache --dir /opt/vep/.vep \
        --format vcf \
        --offline --vcf --force_overwrite \
        --overlaps --max_sv_size 30343437 \
        -fasta ./Bos_taurus_hybrid.UOA_Angus_1.dna.toplevel.fa.gz \
        --hgvs --symbol \
        --input_file high_conf_${ref}.vcf \
        --output_file high_conf_${ref}_vep.vcf \
        --stats_text high_conf_${ref}_vep.stats

