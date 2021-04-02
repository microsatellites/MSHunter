!/bin/bash
output_dir="/mnt/project/ProjectMSI/MSCalling/note/py/test_mshunter/"
/home/pengjia/miniconda3/envs/default/bin/python /mnt/project/ProjectMSI/MSCalling/note/py/MSHunter/src/main.py genotype \
  -i /media/pengjia/ZS6Y_Feng_1/GIAB/NA12878/HG001_GRCh38.haplotag.RTG.trio.bam \
  -o ${output_dir}ccs \
  -m /mnt/project/ProjectMSI/MSCalling/Data/ref/microsatellite/GRCh38_l5_m100_20200528.no.list.mshunter_chr1 \
  -r /mnt/project/REF/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  -t 1 -b 1000 -d True -tech ccs --only_simple  False

