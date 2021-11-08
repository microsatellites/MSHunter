#!/bin/bash

output_dir="/home/DATA/ProjectSTR/codes/test/ILM/"
path_microsatellite_regions="/home/DATA/ProjectSTR/data/microsatellite_regions/GRCh38_l5_m100_20200528.no.list_chr8"
path_ref="/home/pengjia/REF/GRCh38_full_analysis_set_plus_decoy_hla/genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
#/home/pengjia/miniconda3/envs/default/bin/python /home/DATA/ProjectSTR/codes/MSHunter/mshunter/main.py genotype \
#  -i /home/DATA/ProjectSTR/codes/test/CCS/HG00733.CCS.GRCh38.bam \
#  -o ${output_dir}ccs \
#  -m ${path_microsatellite_regions} \
#  -r ${path_ref} \
#  -t 1 -b 1000 -d True -tech ccs --only_simple  True
#


/home/pengjia/miniconda3/envs/default/bin/python /home/DATA/ProjectSTR/codes/MSHunter/src/main.py genotype \
  -i /home/DATA/ProjectSTR/codes/test/ILM/HG00733.ILM.GRCh38.bam \
  -o ${output_dir}ilm \
  -m ${path_microsatellite_regions} \
  -r ${path_ref} \
  -t 1 -b 1000 -tech ccs --only_simple  True


