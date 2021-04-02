#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : estimate_error.py
# Author : Peng Jia
# Date   : 2020.09.03
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""
import pysam
import yaml
# from src.global_dict import *
from src.units import *


def class_dis_by_motif(paras):
    path_pre = paras["output_pre"]
    path_dis_parameter = paras["output_tmp"]
    min_support_reads = paras["minimum_support_reads"]
    logger.info("Scanning the distribution file of microsatellite!")
    vcf_file = pysam.VariantFile(path_pre)
    files_motif = {}
    record_num = 0
    for rec in vcf_file.fetch():
        record_num += 1
        record_info = rec.info
        motif = record_info["motif"]
        support_reads = record_info["depth"]
        if support_reads > min_support_reads:
            if motif not in files_motif:
                files_motif[motif] = pysam.VariantFile(path_dis_parameter + "/tmp_motif_" + motif + ".vcf.gz", 'w',
                                                       header=vcf_file.header)
            files_motif[motif].write(rec)

    motif_list = []
    for motif in files_motif:
        files_motif[motif].close()
        pysam.tabix_index(path_dis_parameter + "/tmp_motif_" + motif + ".vcf.gz", force=True, preset="vcf")

        motif_list.append(motif)
    set_value("motif_list", motif_list)


def get_homo_normal_dis(motif_dis_tmp, max_repeat, motif_len):
    homo_list = {}
    for homo in motif_dis_tmp:
        tmp_dis = {}
        for pot in range(0, max_repeat + 1):
            tmp_dis[pot] = 0
        read_num = len(motif_dis_tmp[homo])
        read_num=0
        for site in motif_dis_tmp[homo]:
            read_num+=sum(site.values())
        # print(motif_dis_tmp[homo])
        for oner in motif_dis_tmp[homo]:
            for pot in oner:
                tmp_dis[pot] += oner[pot]
        homo_list[homo] = {}
        for pot in tmp_dis:
            if tmp_dis[pot] > 0:
                homo_list[homo][pot] = round(tmp_dis[pot] / read_num, 6)

    for homo in range(0, max_repeat + 1):
        if homo not in homo_list:
            homo_list[homo] = {homo * motif_len: 1}
        x, y = [], []
        for i in homo_list[homo]:
            x.append(i)
            y.append(homo_list[homo][i])
        # print(sum(y))
    return homo_list


def get_one_motif_prosess(paras, motif):
    path_dis_parameter = paras["output_tmp"]
    motif_len = len(motif)
    motif_dis_tmp = {}
    vcffile = pysam.VariantFile(path_dis_parameter + "/tmp_motif_" + motif + ".vcf.gz", "r")
    # path_dis_parameter + "tmp_motif_" + motif + ".bcf
    max_repeat = 0
    repeat_list = []
    for rec in vcffile.fetch():
        recoed_info = rec.info
        # repeatTimes = int(rec.info["repeat_times"])
        # repeat_length = int(recoed_info["ref_repeat_length"])
        repeat_times = int(recoed_info["repeat_times"])
        dis_list = [list(map(int, i.split(":"))) for i in recoed_info["dis"].split("|")]
        this_max_repeat = max([i[0] for i in dis_list])
        max_repeat = max_repeat if max_repeat >= this_max_repeat else this_max_repeat
        support_reads = int(rec.info["depth"])
        dis_list_normal = {}
        for i in dis_list:
            dis_list_normal[i[0]] = i[1] / support_reads
        repeat_list.append(repeat_times)

        if repeat_times not in motif_dis_tmp:
            motif_dis_tmp[repeat_times] = []
        thistmp = motif_dis_tmp[repeat_times]
        thistmp.append(dis_list_normal)
        motif_dis_tmp[repeat_times] = thistmp
    max_repeat = max_repeat // motif_len + 1
    homo_list = get_homo_normal_dis(motif_dis_tmp, max_repeat, motif_len=motif_len)
    mixture = {}
    for first in range(0, max_repeat + 1):
        for second in range(0, max_repeat + 1):
            if first <= second:
                first_dis = homo_list[first]
                second_dis = homo_list[second]
                thismaxture = {}
                for rp in set([i for i in first_dis] + [j for j in second_dis]):
                    if rp not in first_dis:
                        first_dis[rp] = 0
                    if rp not in second_dis:
                        second_dis[rp] = 0
                for rp in first_dis:
                    thismaxture[rp] = round((first_dis[rp] + second_dis[rp]) / 2, 6)
                mixture[first * 1000 + second] = remove_zero_dict(thismaxture)
    with open(path_dis_parameter + "/tmp_motif_" + motif + ".model", "w") as f:
        yaml.dump({"mixture": mixture}, f)
    return {"mixture": mixture, "max_repeat": max_repeat}


def estimate_error():
    args = get_value("paras")
    class_dis_by_motif(args)
    motif_list = get_value("motif_list")
    model = {}
    for motif in motif_list:
        logger.info("Build error model for motif: " + motif)
        model[motif] = get_one_motif_prosess(args, motif)
    with open(args["output_model"], "w") as f:
        yaml.dump(model, f)
    set_value("model", model)
    return model


if __name__ == "__main__":
    pass
