#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : benchmark.py
# Author : Peng Jia
# Date   : 2020.07.13
# Email  : pengjia@stu.xjtu.edu.cn
# Description: gnotyping microsatellite by Next generation data (for Illumina and
               Pacbio CCS data.)
=============================================================================="""
import os
import re
import collections
import pysam
import multiprocessing
# from mshunter.global_dict import *
from mshunter.units import *


class ReadInfo:
    """
    Description: Mutation
    pos_del_prefix: deletion position in prefix
    pos_del_ms:  deletion position in ms region
    pos_del_suffix: deletion position in suffix
    pos_ins_prefix: insertion position in prefix
    pos_ins_suffix: insertion position in suffix
    pos_ins_ms: insertion position in ms region
    pos_snp_prefix: mismatch position in prefix
    pos_snp_ms: mismatch position in ms regions
    pos_snp_suffix: mismatch position in suffix
    var_type_prefix: variant type in prefix
    var_type_suffix: variant type in suffix
    var_type_ms: variant type in ms region
    var_type: variant type in all region
    var_list: variant type list
    """

    def __init__(self):
        self.var_prefix = []
        self.var_suffix = []
        self.var_ms = []
        self.var_type_prefix = []
        self.var_type_suffix = []
        self.var_type_ms = []
        self.var_type = ""
        self.var_type_list = []
        self.deletion_ms = []
        self.insertion_ms = []
        self.mismatch_ms = []
        self.deletion_prefix = []
        self.insertion_prefix = []
        self.mismatch_prefix = []
        self.deletion_suffix = []
        self.insertion_suffix = []
        self.mismatch_suffix = []
        self.direction = True  # True read is forward, False: read is reversed
        self.rpl = 0  # repeat length
        self.read_name = ""
        self.read_str = ""
        self.read_list = []
        self.del_span = "None"  # all ,left,right,
        self.hap = 0  # 0: unknow , 1, 2

    def prefix_var_type(self):
        for var in self.var_prefix:
            if len(var[1]) == len(var[2]):
                self.var_type_prefix.append("SNP")
            elif len(var[1]) > len(var[2]):
                self.var_type_prefix.append("INS")
            elif len(var[1]) < len(var[2]):
                self.var_type_prefix.append("DEL")

    def suffix_var_type(self):
        for var in self.var_suffix:
            if len(var[1]) == len(var[2]):
                self.var_type_suffix.append("SNP")
            elif len(var[1]) > len(var[2]):
                self.var_type_suffix.append("INS")
            elif len(var[1]) < len(var[2]):
                self.var_type_suffix.append("DEL")

    def ms_var_type(self, offset):

        if offset < 0:
            self.var_type_ms.append("DEL")
        if offset > 0:
            self.var_type_ms.append("INS")

        for var in self.var_ms:
            if len(var[1]) == len(var[2]):
                self.var_type_ms.append("SNP")

    def comput(self, offset):
        self.ms_var_type(offset)
        self.prefix_var_type()
        self.suffix_var_type()
        self.var_type_list = self.var_type_prefix + self.var_type_ms + self.var_type_suffix
        var_num = len(self.var_type_list)
        if var_num > 1:
            self.var_type = "Complex"
        elif var_num == 1:
            self.var_type = self.var_type_list[0]
        else:
            self.var_type = "None"

    def show_info(self):
        # if len()
        print("read_name", self.read_name)
        print("read_str", self.read_str)
        print("deletion", self.deletion)
        print("insertion", self.insertion)
        print("mismatch", self.mismatch)
        print("direction", self.direction)
        print("rpl", self.rpl)
        print()


class Microsatellite:
    """
    @chrom : chromsome
    @pos_start: start position of microsatellite
    @motif: repeat unit of microsatellite
    @motif_len = repeat unit length of microsatellite
    @repeat_times =  the repeat times of microsatellite
    @pos_end :end position of microsatellite
    @bam_path: bam file path
    @reference_path: reference path
    @start_pre: start pos of this call, default 10bp upstream the microsatellite
    @end_suf:  start pos of this call, default 10bp downstream the microsatellite
    @ref_repeat_length: the reference repeat length
    @prefix_len: how many bps to detect upstream
    @suffix_len = how many bps to detect downstream
    @repeat_length_dis: repeat length distribution
    @query_repeat_length： the repeat length supported by the most reads
    dis_stat: True, at least one read covered this microsatellite
    more_than_one_alleles: more than one reads covered this microsatellite, and have more alleles
    more_than_one_alleles_ms: more than one reads covered this microsatellite, and have more alleles in microsatellite region
    check: True, could be used as benchmark
    check_stats: Why do not be applied as benchmark
    mismatch: position of mishmatch
    ms_var_type: variant type in ms region (DEL,INS,SNP)
    up_var_type: variant type in upstream of ms region (DEL,INS,SNP)
    down_var_type: variant type in downstream of ms region (DEL,INS,SNP)
    microsatellite_id: microsatellite_id (chrom_pos)
    ref_str: reference string
    alt_str: alternative string
    allele: allele in this haplotype
    mut_start: mutation start
    mut_end: mutation end
    mut_type: mutation type and details
    """

    def __init__(self,
                 chrom,
                 pos_start,
                 pos_end,
                 motif,
                 motif_len,
                 repeat_times,
                 bam_path,
                 reference_path,
                 prefix_len=0,
                 suffix_len=0,
                 min_mapping_qual=1,
                 phasing=False,
                 min_allele_fraction=0.2,
                 minimum_support_reads=2
                 ):
        """
        @param chrom: chromsome
        @param pos_start: start position
        @param pos_end: end position
        @param motif: repeat unit
        @param motif_len: repeat unit length
        @param repeat_times: repeat times
        @param bam_path: bam path
        @param reference_path: reference path
        @prefix_len: {prefix_len} bps upstream of microsatellite to analysis
        @suffix_len: {suffix_len} bps downstream of microsatellite to analysis
        """
        self.chrom = str(chrom)
        self.pos_start = pos_start
        self.motif = motif
        self.motif_len = motif_len
        self.repeat_times = repeat_times
        self.pos_end = pos_end
        self.bam_path = bam_path
        self.reference_path = reference_path
        self.ref_repeat_length = repeat_times * motif_len
        self.mirosatellite_id = chrom + "_" + str(pos_start)
        self.prefix_len = prefix_len
        self.suffix_len = suffix_len
        self.start_pre = pos_start - self.prefix_len
        self.end_suf = pos_end + self.suffix_len
        self.mut_start = pos_start - 1
        self.mut_end = pos_end
        self.repeat_length_dis = {}
        self.query_repeat_length = self.ref_repeat_length
        self.dis_stat = False
        self.more_than_one_alleles = False
        self.more_than_one_alleles_ms = False
        self.check = True
        self.check_stats = []
        self.ref_str = "."
        self.alt_str = "."
        self.allele = 0
        self.mut_type = None
        self.min_mapping_qual = min_mapping_qual
        self.reads_info = {}
        self.phasing = phasing
        self.support_reads = 0
        self.dis = {}
        self.dis_hap = {}
        self.hap_info_mismatch = {}
        self.direction_info_mismatch = {}
        self.min_allele_fraction = min_allele_fraction
        self.minimum_support_reads = minimum_support_reads

    def get_reads_alignment(self, bam_file):
        alignment_list = [alignment for alignment in bam_file.fetch(self.chrom, self.pos_start - 1, self.pos_end + 1)]
        reads_com = []
        for alignment in alignment_list:
            if alignment.is_unmapped or alignment.is_duplicate or alignment.is_secondary:
                continue
            # if alignment.reference_start > self.pos_start or alignment.reference_end < self.pos_end:
            if alignment.reference_start > self.start_pre or alignment.reference_end < self.end_suf:
                continue
            if alignment.mapping_quality < self.min_mapping_qual: continue
            reads_com.append(alignment)
        return reads_com

    def get_dis(self):
        """
        Description: get the distribution of the microsateliite repeat length
        """
        bam_file = pysam.AlignmentFile(self.bam_path, mode="rb", reference_filename=self.reference_path)
        repeat_length_dict = {}
        reads = self.get_reads_alignment(bam_file)
        bam_file.close()
        for alignment in reads:
            repeat_length = self.get_repeat_length(alignment)
            if repeat_length < 0: continue
            if repeat_length not in repeat_length_dict:
                repeat_length_dict[repeat_length] = 0
            repeat_length_dict[repeat_length] += 1
        self.repeat_length_dis = repeat_length_dict
        self.allele = len(repeat_length_dict)
        if self.allele < 1:

            self.check = False
            self.check_stats.append("No_read_covered")
        else:
            self.dis_stat = True
            # print(repeat_length_dict)
            self.query_repeat_length = get_max_support_index(repeat_length_dict)
            if self.allele > 1:
                self.more_than_one_alleles = True
                self.more_than_one_alleles_ms = True
                self.check = False
                self.check_stats.append("More_alleles_in_MS")
        # print(self.query_repeat_length, self.repeat_length_dis, self.ref_repeat_length, "000")

        if not self.check:
            return -1
        else:
            return 1

    def get_reads_info(self):
        bam_file = pysam.AlignmentFile(self.bam_path, mode="rb", reference_filename=self.reference_path)
        repeat_length_dict = {}
        reads = self.get_reads_alignment(bam_file)
        bam_file.close()
        self.ref_str_ms = pysam.FastaFile(self.reference_path).fetch(self.chrom, start=self.start_pre,
                                                                     end=self.end_suf)
        reads_info = {}
        for alignment in reads:
            one_read_info = self.get_repeat_info_from_one_read(alignment)
            reads_info[one_read_info.read_name] = one_read_info
        self.reads_info = reads_info

    def remove_low_fraction(self):
        # for
        candidate_snp = {}

        for pos, pos_info in self.direction_info_mismatch.items():
            read_num = 0
            all_check = True
            reads_num_debug = []
            for dir in ["F", "R"]:
                this_dir_num = len(pos_info[dir]) if dir in pos_info else 0
                read_num += this_dir_num
                reads_num_debug.append([dir, this_dir_num])
                if all_check and this_dir_num < self.min_expect_support_half:
                    all_check = False
            if all_check and read_num >= self.min_expect_support:
                candidate_snp[pos] = read_num
                reads_num_debug.append(["all", read_num])
                print(reads_num_debug)
                print(candidate_snp, self.min_expect_support, self.min_expect_support_half)
                print()

        pass

    def genotype_sits(self):
        pass

    def phasing(self):
        pass

    def genotype_microsatellite(self):
        self.support_reads = len(self.reads_info)
        if self.support_reads < get_value("paras")["minimum_support_reads"]:
            return False
        dis_hap = {0: {}, 1: {}, 2: {}}
        # hap_info_mismatch = {0: {}, 1: {}, 2: {}}

        hap_info_mismatch = {}
        # direction_info_mismatch = {"F": {}, "R": {}}
        direction_info_mismatch = {}
        for read_name, read in self.reads_info.items():
            # print(read.read_list)
            repeat_length = read.rpl
            if repeat_length < 0:
                self.support_reads -= 1
                continue
            hap = read.hap
            direction = "F" if read.direction else "R"
            for mis_info in read.mismatch_prefix + read.mismatch_ms + read.mismatch_suffix:
                if mis_info[0] not in hap_info_mismatch:
                    hap_info_mismatch[mis_info[0]] = {}
                if hap not in hap_info_mismatch[mis_info[0]]:
                    hap_info_mismatch[mis_info[0]][hap] = []

                hap_info_mismatch[mis_info[0]][hap] += [read_name]

                if mis_info[0] not in direction_info_mismatch:
                    direction_info_mismatch[mis_info[0]] = {}
                if direction not in direction_info_mismatch[mis_info[0]]:
                    direction_info_mismatch[mis_info[0]][direction] = []
                direction_info_mismatch[mis_info[0]][direction] += [read_name]
            if repeat_length not in dis_hap[hap]:
                dis_hap[hap][repeat_length] = 0
            dis_hap[hap][repeat_length] += 1
        dis_all = dis_sum([dis_hap[0], dis_hap[1], dis_hap[2]])
        self.dis = dis_all
        self.dis_hap = dis_hap
        self.hap_info_mismatch = hap_info_mismatch
        self.direction_info_mismatch = direction_info_mismatch
        self.min_expect_support = max(2, self.minimum_support_reads, int(self.support_reads * self.min_allele_fraction))
        self.min_expect_support_half = max(1, self.minimum_support_reads * 7 // 20, self.min_expect_support * 7 // 20)
        self.remove_low_fraction()
        # self.remove_strand_bias()
        # snps_pos

        # print(self.dis)
        # TODO

    # print()

    def get_pileup_info(self):
        """
        Description: get the detail mutational information of upstream and downstream
        """
        mut = ReadInfo()
        left_pos = self.end_suf
        right_pos = self.start_pre
        alt_str = []
        fa_file = pysam.FastaFile(self.reference_path)

        self.ref_str = fa_file.fetch(self.chrom, self.start_pre, self.end_suf)

        if self.check:
            bam_file = pysam.AlignmentFile(self.bam_path, mode="rb", reference_filename=self.reference_path)
            self.ref_str = fa_file.fetch(self.chrom, self.start_pre, self.end_suf)
            pos = self.start_pre
            segment_pos = 0

            for pot in bam_file.pileup(self.chrom,
                                       self.start_pre,
                                       self.end_suf,
                                       truncate=True,
                                       fastafile=fa_file):
                pot_alleles = list(set(map(lambda x: x.upper(),
                                           pot.get_query_sequences(mark_matches=True,
                                                                   mark_ends=False,
                                                                   add_indels=True))))
                pot_alleles = collections.Counter(pot_alleles).most_common()
                if pot_alleles[0][0][0] in [",", "."]:
                    alt_str.append(self.ref_str[segment_pos])
                elif pot_alleles[0][0][0] in ["*"]:
                    alt_str.append("")
                elif pot_alleles[0][0][0] in ["A", "G", "C", "T"]:
                    left_pos = min(left_pos, pos)
                    right_pos = max(right_pos, pos)
                    alt_str.append(pot_alleles[0][0][0])
                    if pos < self.pos_start:
                        mut.var_prefix.append([pos, pot_alleles[0][0][0], self.ref_str[segment_pos]])
                    elif pos >= self.pos_end:
                        mut.var_suffix.append([pos, pot_alleles[0][0][0], self.ref_str[segment_pos]])
                    else:
                        mut.var_ms.append([pos, pot_alleles[0][0][0], self.ref_str[segment_pos]])
                if len(pot_alleles[0][0]) > 1:
                    if pot_alleles[0][0][1] == "+":
                        left_pos = min(left_pos, pos)
                        right_pos = max(right_pos, pos)
                        p = re.compile("[\\+\\-][0-9]+")
                        indel_f = p.findall(pot_alleles[0][0])[0]
                        indel_len = int(indel_f[1:])
                        indel_str = pot_alleles[0][0][-indel_len:]
                        alt_str[-1] = alt_str[-1] + indel_str
                        if pos < self.pos_start - 1:
                            mut.var_prefix.append([pos, indel_str, ""])
                        elif pos >= self.pos_end:
                            mut.var_suffix.append([pos, indel_str, ""])
                        else:
                            mut.var_ms.append([pos, indel_str, ""])
                    else:
                        p = re.compile("[\\+\\-][0-9]+")
                        indel_f = p.findall(pot_alleles[0][0])[0]
                        indel_len = int(indel_f[1:])
                        indel_str = pot_alleles[0][0][-indel_len:]
                        del_end = pos + 1 + indel_len
                        left_pos = min(left_pos, pos + 1)
                        right_pos = max(right_pos, del_end)
                        if del_end < self.pos_start:
                            mut.var_prefix.append([pos + 1, '', indel_str])
                        elif pos + 1 >= self.pos_end:
                            mut.var_suffix.append([pos + 1, '', indel_str])
                        else:

                            if del_end < self.pos_end or pos + 1 >= self.pos_start:
                                mut.var_ms.append([pos + 1, '', indel_str])
                            else:
                                self.check = False
                                self.check_stats.append("DEL_Span")
                            # else:
                            #     mut.var_ms.append([pos + 1, '', indel_str[0:self.pos_end - (pos + 1)]])  # TODO
                            #     mut.var_suffix.append([self.pos_end, '', indel_str[self.pos_end - (pos + 1):]])  # TODO

                pos += 1
                segment_pos += 1
            bam_file.close()

            self.mut_start = min(left_pos, self.pos_start - 1)
            self.mut_end = max(right_pos, self.pos_end)
            self.alt_str = "".join(alt_str[self.mut_start - self.start_pre:self.mut_end - self.start_pre])
        else:
            pass
        mut.comput(self.query_repeat_length - self.ref_repeat_length)
        self.ref_str = self.ref_str[self.mut_start - self.start_pre:self.mut_end - self.start_pre]
        self.mut_type = mut
        # if len(self.alt_str) < 1:
        #     print(self.ref_str)
        #     print(self.alt_str)
        #     print(self.pos_start)
        #     print(self.mut_start - self.start_pre, self.mut_end - self.start_pre)
        #     print(self.check_stats)
        #     print()
        fa_file.close()

    def pos_convert_ref2read(self, ref_block: list, read_block: list, pos: int, direction="start") -> tuple:
        """
        Description: get the read position according the reference position and cigar staring
        @param direction:  start of end
        @param ref_block:
        @param read_block:
        @param start:
        @param end:
        @return:
        """
        if direction == "start":
            block_index = 0
            for sub_ref_block in ref_block:
                if sub_ref_block[0] <= pos <= sub_ref_block[1]:
                    # print(sub_ref_block,pos)
                    break
                block_index += 1

            if pos == ref_block[block_index][1]:  # M|D M|I D|M
                read_pos = read_block[block_index][1]
                if ref_block[block_index][2] == 2:
                    self.start_pre = min([self.start_pre, ref_block[block_index][0] - 1])
                return pos, read_pos

            if ref_block[block_index][2] == 0:  # match
                read_pos = read_block[block_index][0] + (pos - ref_block[block_index][0])
                return pos, read_pos
            elif ref_block[block_index][2] == 2:  # deletion
                # pos = pos
                read_pos = read_block[block_index][1]
                self.start_pre = min([ref_block[block_index][0] - 1, self.start_pre])

                # pos = ref_block[block_index - 1][1] - 1
                # read_pos = read_block[block_index - 1][1] - 1
                return pos, read_pos
            else:
                pos = ref_block[block_index - 1][1] - 1
                read_pos = read_block[block_index - 1][1] - 1
                return pos, read_pos

        if direction == "end":
            block_index = len(ref_block) - 1
            for sub_ref_block in ref_block[::-1]:
                if sub_ref_block[0] <= pos <= sub_ref_block[1]:
                    # print(sub_ref_block,pos)
                    break
                block_index = block_index - 1
            if pos == ref_block[block_index][0]:  # D|M I|M M|D
                read_pos = read_block[block_index][0]
                if ref_block[block_index][2] == 2:
                    self.end_suf = max([self.end_suf, ref_block[block_index][1] + 1])
                return pos, read_pos

            if ref_block[block_index][2] == 0:
                read_pos = read_block[block_index][0] + (pos - ref_block[block_index][0])
                return pos, read_pos
            elif ref_block[block_index][2] == 2:
                # pos = pos
                read_pos = read_block[block_index][0]
                self.end_suf = max([self.end_suf, ref_block[block_index][1] + 1])
                # pos = ref_block[block_index + 1][0] + 1
                # read_pos = read_block[block_index + 1][0] + 1
                return pos, read_pos
            else:
                pos = pos + 1
                read_pos = read_block[block_index - 1][0] - 1
                return pos, read_pos

    def get_repeat_info_from_one_read(self, alignment):
        # TODO:  To be optimized
        align_start = alignment.reference_start
        align_end = alignment.reference_end
        this_ref_str = pysam.FastaFile(self.reference_path).fetch(self.chrom, start=align_start, end=align_end)
        this_read_str = alignment.query_sequence
        sub_read_str = []
        sub_ref_str = []
        ref_block = []
        read_block = []
        read_mut_info = ReadInfo()
        read_mut_info.direction = False if alignment.is_reverse else True
        if alignment.has_tag("HP"):
            read_mut_info.hap = int(alignment.get_tag("HP"))
        read_pos = 0
        ref_pos = align_start
        ref_pos_str = 0
        for cigartuple in alignment.cigartuples:
            if cigartuple[0] in [0, 7, 8]:  # 0 : M : match or mishmatch ; 7: :=:match; 8:X:mismatch
                ref_block.append((ref_pos, ref_pos + cigartuple[1], 0))
                read_block.append((read_pos, read_pos + cigartuple[1], 0))
                match_ref = list(this_ref_str[ref_pos_str:ref_pos_str + cigartuple[1]])
                match_read = list(this_read_str[read_pos:read_pos + cigartuple[1]])
                sub_ref_str.extend(match_ref)
                sub_read_str.extend(match_read)
                if ref_pos + cigartuple[1] < self.start_pre or ref_pos > self.end_suf:
                    pass
                else:
                    pos = ref_pos - 1
                    for ref_pot, read_pot in zip(match_ref, match_read):
                        pos += 1
                        if pos < self.start_pre or pos > self.end_suf: continue
                        if read_pot == "N" or ref_pot == "N": continue
                        if read_pot == ref_pot: continue
                        if pos < self.pos_start:
                            read_mut_info.mismatch_prefix.append([pos, ref_pot, read_pot])
                        elif pos >= self.pos_end:
                            read_mut_info.mismatch_suffix.append([pos, ref_pot, read_pot])
                        else:
                            read_mut_info.mismatch_ms.append([pos, ref_pot, read_pot])
                read_pos += cigartuple[1]
                ref_pos += cigartuple[1]
                ref_pos_str += cigartuple[1]
            elif cigartuple[0] in [1, 4, 5]:  # 1:I:inserion ;4:S:soft clip 5:H:hardclip
                ref_block.append((ref_pos, ref_pos + 0, 1))
                read_block.append((read_pos, read_pos + cigartuple[1], 1))
                if cigartuple[0] == 1:
                    # print(read_pos)
                    # print(alignment.cigartuples)
                    sub_read_str[-1] += this_read_str[read_pos:read_pos + cigartuple[1]]
                    pos = ref_pos - 1
                    if ref_pos + cigartuple[1] < self.start_pre or ref_pos > self.end_suf:
                        pass
                    else:
                        if pos < self.start_pre or pos > self.end_suf: continue
                        if pos < self.pos_start - 1:
                            read_mut_info.insertion_prefix.append([pos, sub_ref_str[-1], sub_read_str[-1]])
                        elif pos >= self.pos_end:
                            read_mut_info.insertion_suffix.append([pos, sub_ref_str[-1], sub_read_str[-1]])
                        else:
                            read_mut_info.insertion_ms.append([pos, sub_ref_str[-1], sub_read_str[-1]])
                read_pos += cigartuple[1]
            elif cigartuple[0] in [2, ]:  # 2:D; 3:N: skip region of reference
                ref_block.append((ref_pos, ref_pos + cigartuple[1], 2))
                read_block.append((read_pos, read_pos, 2))
                delete_str = this_ref_str[ref_pos_str:ref_pos_str + cigartuple[1]]
                sub_ref_str.extend(list(delete_str))
                sub_read_str.extend([""] * cigartuple[1])
                pos = ref_pos
                end_pos = ref_pos + cigartuple[1]
                if ref_pos < self.start_pre:  # ref_pos + cigartuple[1]
                    if end_pos < self.start_pre:
                        pass
                    elif end_pos < self.end_suf:
                        read_mut_info.del_span = "Left"
                    else:
                        read_mut_info.del_span = "All"
                elif ref_pos > self.end_suf:
                    pass
                else:
                    if end_pos >= self.end_suf:
                        read_mut_info.del_span = "Right"
                    else:
                        # if pos < self.start_pre or pos > self.end_suf: continue
                        if pos < self.pos_start:
                            read_mut_info.deletion_prefix.append([pos, delete_str, ""])
                        elif pos >= self.pos_end:
                            read_mut_info.deletion_suffix.append([pos, delete_str, ""])
                        else:
                            read_mut_info.deletion_ms.append([pos, delete_str, ""])
                    # if
                ref_pos += cigartuple[1]
                ref_pos_str += cigartuple[1]
            else:
                return -1
        this_repeat_length = len("".join(sub_read_str[self.pos_start - 1 - align_start:self.pos_end - align_start - 1]))
        read_mut_info.rpl = this_repeat_length
        read_mut_info.read_name = alignment.query_name + "_" + ("1" if alignment.is_read1 else "2")
        read_mut_info.read_str = "".join(sub_read_str[self.start_pre - align_start:self.end_suf - align_start])
        read_mut_info.read_list = sub_read_str[self.start_pre - align_start:self.end_suf - align_start]
        return read_mut_info

    def get_repeat_length(self, alignment):
        """
        Description: get the repeat length according to the aligned read.
        """
        align_start = alignment.reference_start
        ref_block = []
        read_block = []
        read_pos = 0
        ref_pos = align_start
        for cigartuple in alignment.cigartuples:
            if cigartuple[0] in [0, 7, 8]:  # 0 : M : match or mishmatch ; 7: :=:match; 8:X:mismatch
                ref_block.append((ref_pos, ref_pos + cigartuple[1], 0))
                read_block.append((read_pos, read_pos + cigartuple[1], 0))
                read_pos += cigartuple[1]
                ref_pos += cigartuple[1]
            elif cigartuple[0] in [1, 4, 5]:  # 1:I:inserion ;4:S:soft clip 5:H:hardclip
                ref_block.append((ref_pos, ref_pos + 0, 1))
                read_block.append((read_pos, read_pos + cigartuple[1], 1))
                read_pos += cigartuple[1]
            elif cigartuple[0] in [2, ]:  # 2:D; 3:N: skip region of reference
                ref_block.append((ref_pos, ref_pos + cigartuple[1], 2))
                read_block.append((read_pos, read_pos, 2))
                ref_pos += cigartuple[1]
            else:
                return -1

        if self.pos_start >= ref_block[-1][1] or self.pos_start <= ref_block[0][0]:
            return -1
        if self.pos_end >= ref_block[-1][1] or self.pos_end <= ref_block[0][0]:
            return -1

        ref_start, read_start = self.pos_convert_ref2read(ref_block, read_block, self.pos_start, direction="start")
        ref_end, read_end = self.pos_convert_ref2read(ref_block, read_block, self.pos_end, direction="end")
        rpt = self.ref_repeat_length + ((read_end - read_start) - (ref_end - ref_start))
        return rpt


def ngs_process_one_ms_site(microsatellite):
    # msDetail.get_dis()
    # if msDetail.check:
    # msDetail.get_pileup_info()
    microsatellite.get_reads_info()
    microsatellite.genotype_microsatellite()
    # print()
    return microsatellite


def ngs_write_vcf_init(outputpath):
    outputfile = pysam.VariantFile(outputpath, "w")
    bam_file = pysam.AlignmentFile(get_value("paras")["input"], "rb")
    contigs = bam_file.references
    contigsLen = bam_file.lengths
    chromList = get_value("chrom_list")
    contigs_len_dict = {}
    sortedContig = []
    for contig, length in zip(contigs, contigsLen):
        if contig in chromList:
            sortedContig.append(contig)
            contigs_len_dict[contig] = length
    for contig in sortedContig:
        outputfile.header.add_line(
            "##contig=<ID={chrom},length={length}>".format(chrom=contig, length=contigs_len_dict[contig]))
    set_value("contigs_info", contigs_len_dict)
    outputfile.header.add_line('##INFO=<ID=chrom,Number=1,Type=String,Description="Chromosome">')
    outputfile.header.add_line('##INFO=<ID=pos,Number=1,Type=Integer,Description="Position">')
    outputfile.header.add_line('##INFO=<ID=ms_start,Number=1,Type=Integer,Description='
                               '"Position start of microsatellite">')
    outputfile.header.add_line('##INFO=<ID=ms_end,Number=1,Type=Integer,Description='
                               '"Position end of microsatellite">')
    outputfile.header.add_line('##INFO=<ID=motif,Number=1,Type=String,Description="Repeat unit">')
    outputfile.header.add_line('##INFO=<ID=repeat_times,Number=1,Type=Integer,Description='
                               '"Repeat times of motif in reference">')
    outputfile.header.add_line('##INFO=<ID=motif_len,Number=1,Type=Integer,Description="Repeat unit length">')
    outputfile.header.add_line('##INFO=<ID=ref_repeat_length,Number=1,Type=Integer,Description='
                               '"length of microsatellite">')
    outputfile.header.add_line('##INFO=<ID=start_pre,Number=1,Type=Integer,Description="Start position for analysis">')
    outputfile.header.add_line('##INFO=<ID=end_suf,Number=1,Type=Integer,Description="End position for analysis">')
    outputfile.header.add_line('##INFO=<ID=mut_start,Number=1,Type=Integer,Description="Start position of mutaiton">')
    outputfile.header.add_line('##INFO=<ID=mut_end,Number=1,Type=Integer,Description="End position of mution">')
    outputfile.header.add_line('##INFO=<ID=query_repeat_length,Number=1,Type=Integer,Description='
                               '"Evaluation repeat length of microsatellite">')
    outputfile.header.add_line('##INFO=<ID=dis_stat,Number=1,Type=String,Description='
                               '"True,the distribution is available">')
    outputfile.header.add_line('##INFO=<ID=check,Number=1,Type=String,Description='
                               '"True,the site is available for benchmark">')
    outputfile.header.add_line('##INFO=<ID=check_stats,Number=1,Type=String,Description='
                               '"Why this site is not available for benchmark">')
    outputfile.header.add_line('##INFO=<ID=allele,Number=1,Type=Integer,Description="Allele number in this site">')
    outputfile.header.add_line('##INFO=<ID=dis,Number=1,Type=String,Description='
                               'Distribution of repeat length>')
    outputfile.header.add_line('##INFO=<ID=var_type,Number=1,Type=String,Description='
                               '"Variant typeComplex, SNP, DEL, INS">')
    outputfile.header.add_line('##INFO=<ID=var_type_list,Number=1,Type=String,Description='
                               '"Variant type, Complex, SNP, DEL, INS">')
    outputfile.header.add_line('##INFO=<ID=var_detail,Number=1,Type=String,Description='
                               '"Variant Details, prefix|ms|suffix, prefix: record1:record2, record1: pos-alt.ref">')
    return outputfile


def ngs_write_vcf_close(outputfile):
    outputfile.close()
    pysam.tabix_index(get_value("paras")["output_dis"], force=True, preset="vcf")


def ngs_write_vcf(outputfile, dataList):
    for msDetail in dataList:
        vcfrec = outputfile.new_record()
        # print("infoKey",vcfrec.info.keys())
        vcfrec.contig = msDetail.chrom
        # vcfrec.stop = pos + msDetail.repeat_times * len(msDetail.motif)
        vcfrec.pos = msDetail.mut_start
        vcfrec.ref = msDetail.ref_str
        vcfrec.alts = (msDetail.alt_str,) if msDetail.alt_str != "" else ("N",)
        vcfrec.id = msDetail.mirosatellite_id
        vcfrec.stop = msDetail.pos_end
        vcfrec.info["ms_start"] = msDetail.pos_start
        vcfrec.info["ms_end"] = msDetail.pos_end
        vcfrec.info["motif"] = msDetail.motif
        vcfrec.info["repeat_times"] = msDetail.repeat_times
        vcfrec.info["motif_len"] = msDetail.motif_len
        vcfrec.info["ref_repeat_length"] = msDetail.ref_repeat_length
        vcfrec.info["start_pre"] = msDetail.start_pre
        vcfrec.info["end_suf"] = msDetail.end_suf
        vcfrec.info["mut_start"] = msDetail.mut_start
        vcfrec.info["mut_end"] = msDetail.mut_end
        vcfrec.info["query_repeat_length"] = msDetail.query_repeat_length
        vcfrec.info["dis_stat"] = str(msDetail.dis_stat)
        vcfrec.info["check"] = str(msDetail.check)
        vcfrec.info["check_stats"] = "|".join(msDetail.check_stats)
        vcfrec.info["dis"] = "|".join(
            [str(key) + ":" + str(value) for key, value in msDetail.repeat_length_dis.items()])
        vcfrec.info["allele"] = msDetail.allele
        # if msDetail.check:
        vcfrec.info["var_type"] = msDetail.mut_type.var_type
        vcfrec.info["var_type_list"] = '|'.join([":".join(msDetail.mut_type.var_type_prefix),
                                                 ":".join(msDetail.mut_type.var_type_ms),
                                                 ":".join(msDetail.mut_type.var_type_suffix)])
        # print(msDetail.mut_type.var_prefix,msDetail.mut_type.var_ms,msDetail.mut_type.var_suffix)
        vcfrec.info["var_detail"] = \
            "!".join([
                "|".join([":".join([str(one[0]), one[1], one[2]]) for one in msDetail.mut_type.var_prefix]),
                "|".join([":".join([str(one[0]), one[1], one[2]]) for one in msDetail.mut_type.var_ms]),
                "|".join([":".join([str(one[0]), one[1], one[2]]) for one in msDetail.mut_type.var_suffix]),
            ])
        outputfile.write(vcfrec)


def ngs_multi_run(thread, datalist):
    # pool = multiprocessing.Pool(processes=thread)
    # result_list = pool.map(ngs_process_one_ms_site, datalist)
    # pool.close()
    # pool.join()
    result_list = []
    for ms in datalist:
        result_list.append(ngs_process_one_ms_site(ms))

    # print("input",len(datalist))
    # print("output",len(result_list))

    return result_list


def genotype_ngs(paras):
    args = get_value("paras")
    dis = args["output_dis"]
    thread = args["threads"]
    batch = args["batch"]
    min_mapping_quailty = args["minimum_mapping_quality"]
    min_allele_fraction = args["min_allele_fraction"]
    minimum_support_reads = args["minimum_support_reads"]

    outputfile = ngs_write_vcf_init(dis)
    contigs_info = get_value("contigs_info")
    df_microsatellites = load_microsatellites(args)
    curentMSNum = 0
    tmp_window = []
    ms_number = get_value("ms_number")
    prefix_len = args["prefix_len"]
    suffix_len = args["suffix_len"]
    for ms_id, info in df_microsatellites.iterrows():
        curentMSNum += 1
        chrom = info["chr"]
        if chrom not in contigs_info:
            continue
        pos_start = int(info["pos"])
        repeat_times = int(info["repeatTimes"])
        motif = info["motif"]
        phasing = args["hap"]
        motif_len = len(motif)
        pos_end = pos_start + motif_len * repeat_times
        this_ms_bm = Microsatellite(chrom=chrom,
                                    pos_start=pos_start,
                                    pos_end=pos_end,
                                    motif=info["motif"],
                                    motif_len=motif_len,
                                    repeat_times=repeat_times,
                                    bam_path=args["input"],
                                    reference_path=args["reference"],
                                    prefix_len=prefix_len,
                                    suffix_len=suffix_len,
                                    min_mapping_qual=min_mapping_quailty,
                                    phasing=phasing,
                                    min_allele_fraction=min_allele_fraction,
                                    minimum_support_reads=minimum_support_reads,
                                    )
        tmp_window.append(this_ms_bm)

        if curentMSNum % (batch * thread) == 0:
            print("[INFO] Build Benchmark: Total", ms_number, "microsatelite, processing:",
                  curentMSNum - batch * thread + 1,
                  "-", curentMSNum, "(" + str(round(curentMSNum / ms_number * 100, 2)) + "%)")
            result_list = ngs_multi_run(thread=thread, datalist=tmp_window)
            tmp_window = []
            # ngs_write_vcf(outputfile, result_list)
            # ngs_write_vcf(outputfile, result_list)
            # ngs_write_vcf_close(outputfile)

    print("[INFO] Build Benchmark: Total", ms_number, "microsatelite, finish all!")
    result_list = ngs_multi_run(thread=thread, datalist=tmp_window)
    ngs_write_vcf(outputfile, result_list)
    ngs_write_vcf_close(outputfile)
