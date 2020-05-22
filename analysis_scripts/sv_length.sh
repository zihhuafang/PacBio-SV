#!/bin/bash

Angus_vcf='/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/fang/pacbio_snake/Angus/sv_combined/sv_type/'
UCD_vcf='/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/fang/pacbio_snake/UCD/sv_combined/sv_type/'

for sv in DEL DUP INS TRA
do
cat ${Angus_vcf}high_conf_Angus_${sv}.vcf | grep -v '#' | awk -v sv="${sv}" '{split($8,a,";"); for(i in a) { split(a[i],b,"="); if (b[1] == "SVLEN") printf(b[2] "\t");  } print sv }' > ${Angus_vcf}${sv}_length.txt
done


for sv in DEL DUP INS INV TRA
do
cat ${UCD_vcf}high_conf_UCD_${sv}.vcf | grep -v '#' | awk -v sv="${sv}" '{split($8,a,";"); for(i in a) { split(a[i],b,"="); if (b[1] == "SVLEN") printf(b[2] "\t");  } print sv }' > ${UCD_vcf}${sv}_length.txt
done
