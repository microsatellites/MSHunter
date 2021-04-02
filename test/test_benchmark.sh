#!/bin/sh

output_dir="/mnt/project/ProjectMSI/MSCalling/note/py/test_mshunter/"
/home/pengjia/miniconda3/envs/default/bin/python /mnt/project/ProjectMSI/MSCalling/note/py/MSHunter/src/main.py benchmark \
  -i /media/pengjia/ZS6Y_Feng_1/NA12878/1000G/1000G_CCS_contig/HG00732-hap1.bam \
  -o ${output_dir}ilm_test_hap1 \
  -m /mnt/project/ProjectMSI/MSCalling/Data/ref/microsatellite/GRCh38_l5_m100_20200528.no.list.mshunter_chr1 \
  -r /mnt/project/REF/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  -t 3 -b 1000 -d True --tech contig
/home/pengjia/miniconda3/envs/default/bin/python /mnt/project/ProjectMSI/MSCalling/note/py/MSHunter/src/main.py benchmark \
  -i /media/pengjia/ZS6Y_Feng_1/NA12878/1000G/1000G_CCS_contig/HG00732-hap2.bam \
  -o ${output_dir}ilm_test_hap2 \
  -m /mnt/project/ProjectMSI/MSCalling/Data/ref/microsatellite/GRCh38_l5_m100_20200528.no.list.mshunter_chr1 \
  -r /mnt/project/REF/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  -t 3 -b 1000 -d True --tech contig

/home/pengjia/miniconda3/envs/default/bin/python /mnt/project/ProjectMSI/MSCalling/note/py/MSHunter/src/main.py benchmark_merge \
  -1 ${output_dir}ilm_test_hap1.vcf.gz -2 ${output_dir}ilm_test_hap1.vcf.gz \
  -o ilm_test.vcf.gz -d True
