#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : tmp.py
# Author : Peng Jia
# Date   : 2020.10.26
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""
model_all = get_value("model")
if self.repeat_unit not in model_all:
    self.check = False
    self.check_status.append("No_error_model")
    self.model_stat = False
    self.ref_str_ms = pysam.FastaFile(self.reference).fetch(self.chrom, self.start, self.end + 1)
    self.alt = "."
    return False
hap1_repeat_length = 0
hap2_repeat_length = 0
model = model_all[self.repeat_unit]["mixture"]
max_repeat = model_all[self.repeat_unit]["max_repeat"]
del model_all
self.model_stat = True
# ms_dis = {}
# ms_dis_strand = {True: {}, False: {}}
# ms_dis_hap = {0: {}, 1: {}, 2: {}}

# for read_id, read_info in self.muts.items():
#     strand = read_info.strand
#     hap = read_info.hap
#     # mut_dict_by_pos = {}
#     # for mut in read_info.mismatches:
#     #     if mut[0] not in mut_dict_by_pos:
#     #         mut_dict_by_pos[mut[0]] = []
#     #     mut_dict_by_pos[mut[0]].append([mut[0], "SNV", read_id, hap, strand, mut])
#     # for mut in read_info.insertions:
#     #     if mut[0] not in mut_dict_by_pos:
#     #         mut_dict_by_pos[mut[0]] = []
#     #     mut_dict_by_pos[mut[0]].append([mut[0], "INS", read_id, hap, strand, mut])
#     # for mut in read_info.deletions:
#     #     if mut[0] not in mut_dict_by_pos:
#     #         mut_dict_by_pos[mut[0]] = []
#     #     mut_dict_by_pos[mut[0]].append([mut[0], "DEL", read_id, hap, strand, mut])
#     if read_info.repeat_length not in ms_dis:
#         ms_dis[read_info.repeat_length] = 1
#     else:
#         ms_dis[read_info.repeat_length] += 1
#
#     if read_info.repeat_length not in ms_dis_hap[hap]:
#         ms_dis_hap[hap][read_info.repeat_length] = 0
#     ms_dis_hap[hap][read_info.repeat_length] += 1
#     if read_info.repeat_length not in ms_dis_strand[strand]:
#         ms_dis_strand[strand][read_info.repeat_length] = 0
#     ms_dis_strand[strand][read_info.repeat_length] += 1
#
# self.ms_dis = ms_dis
# self.ms_dis_hap0 = ms_dis_hap[0]
# self.ms_dis_hap1 = ms_dis_hap[1]
# self.ms_dis_hap2 = ms_dis_hap[2]
# self.support_hap0 = len(ms_dis_hap[0])
# self.support_hap1 = len(ms_dis_hap[1])
# self.support_hap2 = len(ms_dis_hap[2])
# self.ms_dis_forward = ms_dis_strand[True]
# self.ms_dis_reversed = ms_dis_strand[False]
# if self.depth - self.support_hap0 < self.depth * 0.5:  # TODO add in command parameters
#     self.reads_phased = False
# if abs(self.support_hap1 - self.support_hap2) > self.depth * 0.4:  # TODO add in input arguments
#     self.reads_phased = False
# elif self.support_hap0 > self.depth * 0.5:  # TODO add in input arguments
#     self.reads_phased = False
# else:
#     self.reads_phased = True
min_repeat = max([min(self.ms_dis.keys()) - 1, 1])
max_repeat = min([max(self.ms_dis.keys()) + 1, max_repeat])

def get_by_model():
    if self.reads_phased:
        model_id_list = [first * 1000 + first for first in range(min_repeat, max_repeat + 1)]
        ms_dis_hap1_normal = {k: v / self.support_hap1 for k, v in self.ms_dis_hap1.items()}
        ms_dis_hap2_normal = {k: v / self.support_hap2 for k, v in self.ms_dis_hap2.items()}
        distance_dict_hap1 = {}
        distance_dict_hap2 = {}
        for model_id in model_id_list:
            distance_dict_hap1[model_id] = get_disdistance(ms_dis_hap1_normal, model[model_id])
            distance_dict_hap2[model_id] = get_disdistance(ms_dis_hap2_normal, model[model_id])
        distance_tuple = sorted(distance_dict_hap1.items(), key=lambda kv: (kv[1], kv[0]))
        hap1_repeat_length = distance_tuple[0][0] // 1000
        if len(distance_dict_hap1) > 1:
            first_two_distance_ratio_hap1 = (distance_tuple[1][1] - distance_tuple[0][1]) / (
                    distance_tuple[0][1] + 0.000001)
            # first_two_distance_hap1 = distance_tuple[1][1] - distance_tuple[0][1]
            # distance_hap1 = distance_tuple[0][1]
            qual_hap1 = first_two_distance_ratio_hap1 * self.support_hap1
        # else:
        #     qual_hap1 = 0

        distance_tuple = sorted(distance_dict_hap2.items(), key=lambda kv: (kv[1], kv[0]))
        hap2_repeat_length = distance_tuple[0][0] // 1000
        if len(distance_dict_hap2) > 1:
            first_two_distance_ratio_hap2 = (distance_tuple[1][1] - distance_tuple[0][1]) / (
                    distance_tuple[0][1] + 0.000001)
            # first_two_distance_hap2 = distance_tuple[1][1] - distance_tuple[0][1]
            # distance_hap2 = distance_tuple[0][1]
            qual_hap2 = first_two_distance_ratio_hap2 * self.support_hap2
        # else:
        #     qual_hap2 = 0
        self.qual_ms = round(0.5 * (qual_hap1 + qual_hap2), 2)
        self.qual_ms_hap1 = round(qual_hap1, 2)
        self.qual_ms_hap2 = round(qual_hap2, 2)
    else:
        model_id_list = []
        for first in range(min_repeat, max_repeat + 1):
            for second in range(min_repeat, max_repeat + 1):
                if first <= second:
                    model_id_list.append(first * 1000 + second)
        ms_dis_normal = {k: v / self.depth for k, v in self.ms_dis.items()}
        distance_dict = {}
        for model_id in model_id_list:
            distance_dict[model_id] = get_disdistance(ms_dis_normal, model[model_id])
        distance_tuple = sorted(distance_dict.items(), key=lambda kv: (kv[1], kv[0]))
        hap1_repeat_length = distance_tuple[0][0] // 1000
        hap2_repeat_length = distance_tuple[0][0] % 1000
        if len(distance_dict) > 1:
            firsttwoDistanceRatio = (distance_tuple[1][1] - distance_tuple[0][1]) / (
                    distance_tuple[0][1] + 0.000001)
            # first_two_distance = distance_tuple[1][1] - distance_tuple[0][1]
            # distance = distance_tuple[0][1]
            self.qual_ms = round(firsttwoDistanceRatio * self.support_reads, 2)