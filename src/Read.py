#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : Read.py
# Author : Peng Jia
# Date   : 2020.08.05
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""
from src.units import *
import pysam
import numpy as np


class Read:
    """
    Description: Read
    """

    def __init__(self, read_id, alignment, reference, chrom, tech=""):
        self.chrom = chrom
        self.read_name = alignment.query_name
        self.align_start = alignment.reference_start
        self.align_end = alignment.reference_end
        self.this_read_str = alignment.query_sequence.upper()
        self.tech = tech
        # print(alignment)
        if tech == "contig":
            self.this_read_quals = []
        else:
            self.this_read_quals = "".join([chr(i + 33) for i in alignment.query_qualities])
        self.strand = False if alignment.is_reverse else True  # True: forward False: reverse
        self.this_read_list = []
        self.this_quals_list = []
        self.this_ref_str = ""
        self.this_ref_list = []
        self.read_id = read_id
        self.reference = reference
        self.support_microsatellites = []
        # self.alignment = alignment
        self.microsatellites = {}
        self.direction = False if alignment.is_reverse else True
        self.hap = int(alignment.get_tag("HP")) if alignment.has_tag("HP") else 0
        # print(self.hap)
        self.cigartuples = alignment.cigartuples
        self.mut_info = {}
        self.repeat_lengths = {}

    # def get_microsatellite_detail(self, ms_info):
    #     self. = ms_info

    def get_read_str(self):
        pass
        # TODO:  To be optimized

        self.this_ref_str = pysam.FastaFile(self.reference).fetch(self.chrom, start=self.align_start,
                                                                  end=self.align_end).upper()
        self.this_ref_list = list(self.this_ref_str)
        sub_read_str = []
        sub_read_quals = []

        # read_mut_info = ReadInfo()
        # read_mut_info.direction = self.direction
        # read_mut_info.hap = self.hap
        read_pos = 0
        for cigartuple in self.cigartuples:
            if cigartuple[0] in [0, 7, 8]:  # 0 : M : match or mishmatch ; 7: :=:match; 8:X:mismatch
                match_read = list(self.this_read_str[read_pos:read_pos + cigartuple[1]])
                match_quals = list(self.this_read_quals[read_pos:read_pos + cigartuple[1]])
                sub_read_str.extend(match_read)
                sub_read_quals.extend(match_quals)
                read_pos += cigartuple[1]
            elif cigartuple[0] in [1, 4, 5]:  # 1:I:inserion ;4:S:soft clip 5:H:hardclip
                if cigartuple[0] == 1:
                    if len(sub_read_str) == 0:
                        continue
                    else:
                        sub_read_str[-1] += self.this_read_str[read_pos:read_pos + cigartuple[1]]
                        if self.tech == "contig": continue
                        sub_read_quals[-1] += self.this_read_quals[read_pos:read_pos + cigartuple[1]]
                elif cigartuple[0] == 5:
                    continue
                read_pos += cigartuple[1]
            elif cigartuple[0] in [2, 3]:  # 2:D; 3:N: skip region of reference
                sub_read_str.extend([""] * cigartuple[1])
                if self.tech == "contig": continue
                sub_read_quals.extend([""] * cigartuple[1])
            else:
                return -1
        self.this_read_list = sub_read_str
        # self.this_read_str = ""
        self.this_quals_list = sub_read_quals
        self.this_read_quals = ""

    def get_repeat_length(self, ms_start, ms_end):  # give a start and end
        query_repeat_length = len(
            "".join(self.this_read_list[ms_start - 1 - self.align_start:ms_end - self.align_start - 1]))
        return query_repeat_length

    def get_repeat_length_all_ms(self):  # return all microsatellite covered
        self.microsatellites_num = len(self.microsatellites)
        repeat_lengths = {}
        for ms_id, ms in self.microsatellites.items():
            ms_start = ms.start
            ms_end = ms.end
            query_repeat_length = len(
                "".join(self.this_read_list[ms_start - 1 - self.align_start:ms_end - self.align_start - 1]))
            repeat_lengths[ms_id] = query_repeat_length
        self.repeat_lengths = repeat_lengths

    def get_quals(self, q_start, q_end):
        quals_list = self.this_quals_list[q_start - self.align_start - 1:q_end - self.align_start - 1]
        return quals_list, self.strand

    def get_ms_info_one_read(self):
        self.microsatellites_num = len(self.microsatellites)
        # print(self.read_id, self.microsatellites_num,len(self.support_microsatellites))
        read_muts = {}
        repeat_lengths = {}
        for ms_id, ms in self.microsatellites.items():
            ms_start = ms.start
            ms_end = ms.end
            ms_start_pre = ms.start_pre
            ms_end_suf = ms.end_suf
            query_repeat_length = len(
                "".join(self.this_read_list[ms_start - 1 - self.align_start:ms_end - self.align_start - 1]))
            prefix = self.this_read_list[ms_start_pre - self.align_start:ms_start - self.align_start]
            suffix = self.this_read_list[ms_end - self.align_start:ms_end_suf - self.align_start]
            ms_content = self.this_read_list[ms_start - self.align_start:ms_end - self.align_start]
            mismatches = []
            deletions = []
            insertions = []
            pos_based_info = {}
            pos_deletion = []
            ref_pos = ms_start_pre - 2
            for pot in range(ms_start_pre - self.align_start - 1, ms_end_suf - self.align_start):
                ref_pos += 1
                this_read_base = self.this_read_list[pot]
                this_ref_base = self.this_ref_list[pot]
                if this_read_base == this_ref_base:
                    continue
                else:
                    if ref_pos < ms_start_pre:
                        band = [1]
                    elif ref_pos < ms_start:
                        band = [2]
                    elif ref_pos < ms_end:
                        band = [3]
                    elif ref_pos < ms_end_suf:
                        band = [4]
                    else:
                        band = [5]
                    this_read_base_len = len(this_read_base)
                    this_ref_base_len = len(this_ref_base)
                    if this_read_base_len == this_ref_base_len:
                        mismatches.append([ref_pos, this_read_base, band])
                        pos_based_info[ref_pos] = ["SNV", this_ref_base, this_read_base, band]
                    else:
                        if this_read_base_len < this_ref_base_len:
                            deletions.append([ref_pos, this_read_base, band])
                            # pos_based_info[ref_pos] = ["DEL", this_ref_base, ""]
                            pos_deletion.append([ref_pos, this_ref_base, band[0]])
                        else:
                            if ref_pos == ms_start - 1 and this_read_base_len > ms.repeat_unit_len:
                                # TODO
                                band = [3]
                            insertions.append([ref_pos, this_read_base, band])
                            if band[0] != 3:
                                pos_based_info[ref_pos] = ["INS", "", this_read_base[1:], band]

            if len(pos_deletion) > 0:
                deletion_start = pos_deletion[0][0]
                deletion_len = 1
                del_str = pos_deletion[0][1]
                band = [pos_deletion[0][2]]
                for i in range(1, len(pos_deletion)):
                    if pos_deletion[i - 1][0] + 1 == pos_deletion[i][0]:
                        deletion_len += 1
                        del_str += pos_deletion[i][1]
                        band.append(pos_deletion[i][2])
                    else:
                        pos_based_info[deletion_start] = ["DEL", del_str, "", set(band)]

                        deletion_len = 1
                        deletion_start = pos_deletion[i][0]
                        del_str = pos_deletion[i][1]
                        band = [pos_deletion[i][2]]
                # print(del_str,deletion_start,band)
                if len(set(band)) > 1 or list(set(band))[0] != 3:
                    pos_based_info[deletion_start] = ["DEL", del_str, "", set(band)]
            read_str = self.this_read_list[ms_start_pre - 1 - self.align_start:ms_end_suf + 2 - self.align_start]
            read_muts[ms_id] = Read_Mutation(self.read_id, repeat_length=query_repeat_length, strand=self.strand,
                                             hap=self.hap,
                                             mismatches=mismatches, deletions=deletions, insertions=insertions,
                                             prefix=prefix, suffix=suffix, ms_content=ms_content,
                                             pos_based_info=pos_based_info,
                                             read_str=read_str
                                             )
            repeat_lengths[ms_id] = query_repeat_length
        self.repeat_lengths = repeat_lengths
        self.mut_info = read_muts
