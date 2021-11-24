# MShunter

MSHunter is still in development. Please report bugs via [GitHub issues](https://github.com/microsatellites/MSHunter/issues/new).

If you have any question about this software, please contact with Peng Jia (pengjia@stu.xjtu.edu.cn)

# Installation

```shell
git clone https://github.com/microsatellites/MSHunter.git
cd MSHunter
mamba env create -f ./mshunter_env.yaml
conda activate mshunter_env
python setup.py
mshunter -h 
```

# Quik start

```shell 
 mshunter genotype  -m GRCh38_noSimpleRepeat_chr1.list  -i HG00731.ILM.GRCh38.bam -r GRCh38_full_analysis_set_plus_decoy_hla.fa -o output_test -tech ilm  -t 20 
```

parameters: 
* -m: microsatellite files, generated by scan module of msisensor-pro ([see example](https://raw.githubusercontent.com/microsatellites/MSHunter/master/test/data/GRCh38_noSimpleRepeat_chr1.list))
* -i aligned bam file or cram file
* -r reference file 
* -o output directory 
* -tech sequencing technology 
* -t threads used
