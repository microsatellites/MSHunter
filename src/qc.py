#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : qc.py
# Author : Peng Jia
# Date   : 2020.08.28
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""
import os
import pysam
from src.units import *
from src.Microsatellite import Microsatellite
from multiprocess.pool import Pool


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
    paras["hap"] = args.haplotype_bam[0]
    paras["prefix_len"] = args.prefix_len[0]
    paras["suffix_len"] = args.suffix_len[0]
    paras["debug"] = args.debug[0]
    paras["only_homopolymer"] = args.only_homopolymers[0]
    paras["minimum_support_reads"] = args.minimum_support_reads[0]
    paras["minimum_mapping_quality"] = args.minimum_mapping_quality[0]
    # paras["allow_mismatch"] = args.allow_mismatch[0]
    paras["threads"] = args.threads[0]
    paras["batch"] = args.batch[0]
    paras["minimum_phasing_reads"] = args.minimum_phasing_reads[0]
    paras["min_allele_fraction"] = args.min_allele_fraction[0]
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
    case = input_path.split("/")[-1].strip(".bam")
    case = case.strip(".cram")
    paras["output_dis"] = paras["output"] + case + ".dis.vcf.gz"
    paras["output_tmp"] = paras["output"] + case + "_tmp"
    if not os.path.exists(paras["output_tmp"]):
        os.makedirs(paras["output_tmp"])
    paras["output_model"] = paras["output"] + case + ".model"
    paras["output_call"] = paras["output"] + case + ".vcf.gz"
    set_value("case", case)
    set_value("paras", paras)
    return True


def process_one_ms(info):
    ms = Microsatellite(info)
    # dis = ms.get_dis()
    dis, quals = ms.get_dis_qual()

    if len(dis) > 0:
        dis_str = "|".join([":".join([str(a), str(b)]) for a, b in dis.items()])
        dis_mean, dis_std = dis_stats(dis, ["mean", 'std'])
        # prefix_str = "".join(map(int2str, list(quals["prefix"])))
        # suffix_str = "".join(map(int2str, list(quals["suffix"])))
        # ms_str = "".join(map(int2str, quals["ms"]))
        return ms.ms_id, ms.repeat_unit, ms.repeat_times, dis_mean, dis_std, dis_str, dis, \
               quals
    else:
        return None, None, None, None


def pre_stat(paras, df_microsatellites):
    # reference=paras["reference"]
    path_pre_stat = paras["output"].rstrip("/") + "/" + get_value("case") + ".stat"
    path_pre_stat_tmp = paras["output_tmp"].rstrip("/") + "/" + get_value("case") + ".stat"
    file_all_stat = open(path_pre_stat, "w")
    file_all_stat.write("\t".join(["repeat_unit_length", "repeat_times", "num_forward", "num_reversed",
                                   "this_repeat_mean_mean", "this_repeat_mean_std",
                                   "this_repeat_std_mean", "this_repeat_std_std",
                                   "forward_prefix", "forward_ms", "forward_suffix",
                                   "reversed_prefix", "reversed_ms", "reversed_suffix"]) + "\n")

    df_microsatellites_download_sample = microsatellites_sampling(df_microsatellites, paras)

    for repeat_unit, info in df_microsatellites_download_sample.items():
        for repeat_times, ms_infos in info.items():
            logger.info("Processing   repeat unit: " + str(repeat_unit) + " repeat times: " + str(repeat_times))
            infos = []
            for id, info in ms_infos.iterrows():
                info["reference"] = paras["reference"]
                info["prefix_len"] = paras["prefix_len"]
                info["suffix_len"] = paras["suffix_len"]
                infos.append(info)
            pool = Pool(processes=paras["threads"])
            res_infos = pool.map(process_one_ms, infos)
            pool.close()
            pool.join()
            suffix_str = "." + str(repeat_unit) + "." + str(repeat_times)
            file = open(path_pre_stat_tmp + suffix_str + ".repeat", "w")
            this_repeat_means = []
            this_repeat_stds = []
            num_forward = 0
            num_reversed = 0
            prefix_forward = []
            suffix_forward = []
            ms_forward = []
            prefix_reversed = []
            suffix_reversed = []
            ms_reversed = []
            for res in res_infos:
                if None not in res:
                    file.write("\t".join(map(str, res[:-2]))+"\n")
                    this_repeat_means.append(res[3])
                    this_repeat_stds.append(res[4])
                    prefix_forward.extend(res[-1]["prefix_forward"])
                    suffix_forward.extend(res[-1]["suffix_forward"])
                    ms_forward.extend(res[-1]["ms_forward"])
                    prefix_reversed.extend(res[-1]["prefix_reversed"])
                    suffix_reversed.extend(res[-1]["suffix_reversed"])
                    ms_reversed.extend(res[-1]["ms_reversed"])
                    num_forward += res[-1]["num_forward"]
                    num_reversed += res[-1]["num_reversed"]

            file.close()
            if num_forward + num_reversed < 2: continue
            this_repeat_mean_mean = np.mean(this_repeat_means)
            this_repeat_mean_std = np.std(this_repeat_means)
            this_repeat_std_mean = np.mean(this_repeat_stds)
            this_repeat_std_std = np.std(this_repeat_stds)
            pd.concat([pd.DataFrame([np.nanmean(np.array(prefix_forward), axis=0)]),
                       pd.DataFrame([np.nanmean(np.array(ms_forward), axis=0)]),
                       pd.DataFrame([np.nanmean(np.array(suffix_forward), axis=0)])
                       ], axis=1, ).to_csv(path_pre_stat_tmp + suffix_str + ".forward.qual")
            pd.concat([pd.DataFrame([np.nanmean(np.array(prefix_reversed), axis=0)]),
                       pd.DataFrame([np.nanmean(np.array(ms_reversed), axis=0)]),
                       pd.DataFrame([np.nanmean(np.array(suffix_reversed), axis=0)])
                       ], axis=1, ).to_csv(path_pre_stat_tmp + suffix_str + ".reversed.qual")
            forward_prefix = np.nanmean(prefix_forward)
            forward_ms = np.nanmean(ms_forward)
            forward_suffix = np.nanmean(suffix_forward)

            reversed_prefix = np.nanmean(prefix_reversed)
            reversed_ms = np.nanmean(ms_reversed)
            reversed_suffix = np.nanmean(suffix_reversed)
            this_info_list = list(map(str, [repeat_unit, repeat_times, num_forward, num_reversed,
                                            this_repeat_mean_mean, this_repeat_mean_std,
                                            this_repeat_std_mean, this_repeat_std_std,
                                            forward_prefix, forward_ms, forward_suffix,
                                            reversed_prefix, reversed_ms, reversed_suffix
                                            ]))
            file_all_stat.write("\t".join(this_info_list) + "\n")
    file_all_stat.close()
    return


def microsatellites_sampling(df_microsatellite, paras, sample_num=1000):
    df_microsatellite_downsample = {}
    repeat_ranges = paras["ranges_of_repeat_times"]
    for repeat_unit_length in sorted(repeat_ranges.keys()):
        df_microsatellite_repeat_unit = df_microsatellite[df_microsatellite["motifLen"] == repeat_unit_length]
        if len(df_microsatellite_repeat_unit) > 0 and repeat_unit_length not in df_microsatellite_downsample:
            df_microsatellite_downsample[repeat_unit_length] = {}
        for repeat_times in range(repeat_ranges[repeat_unit_length]["min"], repeat_ranges[repeat_unit_length]["max"]):
            df_microsatellite_repeat_times = df_microsatellite_repeat_unit[
                df_microsatellite_repeat_unit["repeatTimes"] == repeat_times]
            if len(df_microsatellite_repeat_times) > sample_num:
                df_microsatellite_repeat_times = df_microsatellite_repeat_times.sample(sample_num)
            if len(df_microsatellite_repeat_times) > 0:
                df_microsatellite_downsample[repeat_unit_length][repeat_times] = df_microsatellite_repeat_times
    return df_microsatellite_downsample


def qc(paras):
    if not genotype_init(paras):
        logger.error("QC init ERROR!")
        return
    paras = get_value("paras")
    df_microsatellites = load_microsatellites(paras)
    pre_stat(paras, df_microsatellites)
