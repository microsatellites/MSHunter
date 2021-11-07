#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : genotype.py
# Author : Peng Jia
# Date   : 2020.07.13
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""
from src.ngs import *
from src.ccs import *


def genotype_init(args):
    """
    argument procress
    """
    paras = {}
    paras["input"] = args.input[0]
    paras["output"] = args.output[0]
    paras["microsatellite"] = args.microsatellite[0]
    paras["reference"] = args.reference[0]
    paras["microsatellite_region_format"] = args.microsatellite_region_format[0]
    paras["tech"] = args.technology[0]
    # paras["hap"] = args.haplotype_bam[0]
    paras["prefix_len"] = args.prefix_len[0]
    paras["suffix_len"] = args.suffix_len[0]
    paras["debug"] = True if args.debug[0].lower() == "true" else False
    paras["only_homopolymer"] = True if str(args.only_homopolymers[0]).lower() == "true" else False
    paras["only_simple"] = True if args.only_simple[0].lower() == "true" else False
    paras["using_phasing_info"] = args.using_phasing_info[0]
    paras["minimum_support_reads"] = args.minimum_support_reads[0]
    paras["minimum_mapping_quality"] = args.minimum_mapping_quality[0]
    # paras["allow_mismatch"] = args.allow_mismatch[0]
    paras["threads"] = args.threads[0]
    paras["batch"] = args.batch[0]
    paras["minimum_phasing_reads"] = args.minimum_phasing_reads[0]
    paras["min_allele_fraction"] = args.min_allele_fraction[0]
    paras["sequencing_error"] = args.sequencing_error[0]
    paras["maximum_distance_of_two_complex_events"] = args.maximum_distance_of_two_complex_events[0]
    paras["sample"] = args.sample[0]
    paras["ranges_of_repeat_times"] = {}

    for i in args.minimum_repeat_times[0].split(";"):
        unitRange, repeatRange = i.split(":")
        if "-" in unitRange:
            unitStart, unitEnd = tuple(map(int, unitRange.split("-")))
        else:
            unitStart = int(unitRange)
            unitEnd = unitStart
        repeatStart = int(repeatRange)
        # print(unitStart,unitEnd,"  ",repeatStart, repeatEnd)
        for ur in range(unitStart, unitEnd + 1):
            if ur not in paras["ranges_of_repeat_times"]:
                paras["ranges_of_repeat_times"][ur] = {}
            paras["ranges_of_repeat_times"][ur]["min"] = repeatStart
        for i in args.maximum_repeat_times[0].split(";"):
            # print(i)
            unitRange, repeatRange = i.split(":")
            if "-" in unitRange:
                unitStart, unitEnd = tuple(map(int, unitRange.split("-")))
            else:
                unitStart = int(unitRange)
                unitEnd = unitStart
            repeatStart = int(repeatRange)
            # print(unitStart,unitEnd,"  ",repeatStart, repeatEnd)
            for ur in range(unitStart, unitEnd + 1):
                if ur not in paras["ranges_of_repeat_times"]:
                    paras["ranges_of_repeat_times"][ur] = {}
                paras["ranges_of_repeat_times"][ur]["max"] = repeatStart
    error_stat = False
    if os.path.exists(paras["input"]):
        logger.info("The input file is : " + paras["input"] + ".")
    else:
        logger.error('The input file ' + paras["input"] + ' is not exist, please check again')
        error_stat = True
    if os.path.isfile(paras["microsatellite"]):
        logger.info("The microsatellites file  is : " + paras["microsatellite"])
    else:
        logger.error('The microsatellites file ' + paras["microsatellite"] + ' is not exist, please check again')
        error_stat = True
    if os.path.isfile(paras["reference"]):
        logger.info("The reference file is : " + paras["reference"] + ".")
    else:
        paras["reference"] = "" if paras["reference"] == "." else paras["reference"]
        logger.error('The reference file ' + paras["reference"] + ' is not exist, please check again')
        error_stat = True
    if paras["input"][-4:] == "cram":
        paras["input_format"] = "cram"
        cramfile = pysam.AlignmentFile(paras["input"], mode="rb", reference_filename=paras["reference"])
        if not cramfile.has_index():
            logger.info("Build index for the input cram ...")
            pysam.index(paras["input"])
        cramfile.close()
    else:
        paras["input_format"] = "bam"
        bamfile = pysam.AlignmentFile(paras["input"], mode="rb")
        if not bamfile.has_index():
            logger.info("Build index for the input bam ...")
            pysam.index(paras["input"])
        bamfile.close()

    if not os.path.exists(paras["output"]):
        pass
    else:
        if paras["debug"]:
            pass
        else:
            logger.error('The output ' + paras["output"] +
                         ' is still exist! in case of overwrite files in this workspace, '
                         'please check your script!')
            error_stat = True
    if error_stat:
        return False
    logger.info("The output is : " + paras["output"] + ".")
    output_path = paras["output"]
    output_path = output_path if output_path[-1] == "/" else output_path + "/"

    if not os.path.exists(output_path):
        os.makedirs(output_path)
    paras["output"] = output_path
    input_path = paras["input"]
    input_path = input_path[:-1] if input_path[-1] == "/" else input_path
    if paras["sample"] == "default":
        case = input_path.split("/")[-1].strip(".bam")
        case = case.strip(".cram")
    else:
        case = paras["sample"]
    paras["output_pre"] = paras["output"] + case + ".pre.vcf.gz"
    paras["output_tmp"] = paras["output"] + case + "_tmp"
    if not os.path.exists(paras["output_tmp"]):
        os.makedirs(paras["output_tmp"])
    paras["output_model"] = paras["output"] + case + ".model"
    paras["output_micro"] = paras["output"] + case + "_micro.vcf.gz"
    paras["output_indel"] = paras["output"] + case + "_indel.vcf.gz"
    paras["output_snv"] = paras["output"] + case + "_snv.vcf.gz"
    paras["output_complex"] = paras["output"] + case + "_complex.vcf.gz"
    set_value("case", case)
    set_value("paras", paras)
    # print(paras)
    return True


def genotype(parase):
    if not genotype_init(parase):
        logger.error("Genotype init ERROR!")
        return -1
    paras = get_value("paras")
    if get_value("paras")["tech"] == "ilm":
        genotype_ngs(paras)
    elif get_value("paras")["tech"] == "ccs":
        genotype_ccs(paras)
