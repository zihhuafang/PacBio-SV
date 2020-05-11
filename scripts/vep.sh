#!/bin/env bash

module purge
module load gcc/4.8.2 perl/5.18.4 samtools/1.8
perl-init
perl -e 'use DBI'
PERL5LIB=${PERL5LIB}:/cluster/work/pausch/meenu/vcf_UCD_2019_05/ensembl-vep/biodbhts/lib

vep='/cluster/work/pausch/meenu/vcf_UCD_2019_05/ensembl-vep/vep'
fd='/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/fang/pacbio_snake/'
ref_fd='/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/fang/pacbio_snake/ref/'


for ref in UCD Angus
do
${vep} -i ${fd}${ref}/sv_combined/high_conf_${ref}.vcf \
--gtf ${ref_fd}${ref}_1.99.chr_sorted.gtf.gz \
--fasta ${ref_fd}${ref}_ref.fa \
--hgvs --symbol --overlaps --vcf --force \
-o ${fd}${ref}/sv_combined/high_conf_${ref}_vep.vcf \
--stats_text ${fd}${ref}/sv_combined/high_conf_${ref}_vep.stats
done
