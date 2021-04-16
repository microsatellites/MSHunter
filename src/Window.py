#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : Window.py
# Author : Peng Jia
# Date   : 2020.08.04
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""
import multiprocessing
from src.global_dict import *
from src.Microsatellite import Microsatellite
from src.Read import Read
import pysam


class Window:
    """
    Description: class Window
    Member variables:
        contig: the contig of this window
        paras: parameters of this program, containing command input and some default values
        ms_list: information of microsatellite in this window
        bam_path: bam file path
        win_start: start position of this window
        win_end: end position of this window
        reads: Read Object of this window; dict: read_id:Read
        reads_num: reads numbers
        microsatellites: Microsatellite Object of this window; dict: ms_id:Microsatellite
        microsatellites_id: list of ms_id
        vcf_recs: list of vcf records infomation in this window
    Member methods:

    """

    def __init__(self, ms_info_list, tech=""):
        self.tech = tech
        self.contig = ms_info_list[0]["chr"]
        self.paras = get_value("paras")
        self.ms_list = ms_info_list
        self.bam_path = self.paras["input"]
        # self.threads = self.paras["threads"]
        self.win_start = ms_info_list[0]["pos"] - self.paras["prefix_len"]
        self.win_end = ms_info_list[-1]["pos"] + ms_info_list[-1]["repeatTimes"] * ms_info_list[-1]["motifLen"] + \
                       self.paras["suffix_len"]
        self.reads = {}
        self.reads_num = 0
        self.microsatellites = {}
        self.microsatellites_id = [it["chr"] + "_" + str(it["pos"]) for it in ms_info_list]
        self.vcf_recs = []

    def init_microsatellites(self, only_simple=True):
        """
        Description: Microsatellite class init for this window
        Returns:
        """
        microsatellites = []
        for ms in self.ms_list:
            microsatellites.append(Microsatellite(ms, only_simple=only_simple))
        self.microsatellites = {ms_info.ms_id: ms_info for ms_info in microsatellites}

    def init_reads(self):
        """
        Description: 1. Extract all read covered the microsatellite regions in this window.
                     2. Read Init
        Returns:
        """
        reads = {}
        sam_file = pysam.AlignmentFile(self.paras["input"], mode="rb", reference_filename=self.paras["reference"])
        # pysam.AlignmentFile(self.bam_path, mode="rb", reference_filename=self.reference_path)
        # TODO optimize
        for ms_id, ms_info in self.microsatellites.items():
            for alignment in sam_file.fetch(ms_info.chrom, ms_info.start_pre - 1, ms_info.end_suf + 1):
                if alignment.is_unmapped or alignment.is_duplicate or alignment.is_secondary or alignment.is_supplementary:
                    continue
                if alignment.reference_start > ms_info.start_pre - 1 or alignment.reference_end < ms_info.end_suf + 1:
                    continue
                if len(alignment.query_sequence) < 2:
                    continue
                read_id = alignment.query_name + "_" + str(alignment.reference_start)
                if read_id not in reads:
                    reads[read_id] = Read(read_id=read_id,
                                          chrom=self.contig,
                                          alignment=alignment,
                                          reference=self.paras["reference"],
                                          tech=self.tech)
                if ms_info.ms_id not in reads[read_id].support_microsatellites:
                    reads[read_id].support_microsatellites.append(ms_info.ms_id)
        self.reads = reads
        self.reads_num = len(self.reads)

    def get_reads_dis(self):
        """
        Description: get repeat length distribution
        Returns:
        """
        result_list = []
        for read in self.reads.values():
            read.microsatellites = {ms_id: self.microsatellites[ms_id] for ms_id in read.support_microsatellites}
            read.get_read_str()
            read.get_repeat_length_all_ms()
            result_list.append(read)
        self.reads = {read.read_id: read for read in result_list}

    def get_reads_info(self, only_simple=True):
        """
        Description: get mutation of each read and repeat length destribution
        Returns:
        """
        result_list = []
        for read in self.reads.values():
            read.microsatellites = {ms_id: self.microsatellites[ms_id] for ms_id in read.support_microsatellites}
            read.get_read_str()
            if only_simple:
                read.get_repeat_length_all_ms()
            else:
                read.get_ms_info_one_read()
            result_list.append(read)
        self.reads = {read.read_id: read for read in result_list}

    def merge_reads_repeat_length_distribution(self):
        microsatellites_dict = {ms_id: {} for ms_id in self.microsatellites}
        for read_id, read in self.reads.items():
            stand = read.strand
            hap = read.hap
            for ms_id, ms_read_repeat_length in read.repeat_lengths.items():
                microsatellites_dict[ms_id][read_id] = [ms_read_repeat_length, stand, hap]
        self.reads = {}
        for ms_id, ms_read_repeat_length_info in microsatellites_dict.items():
            self.microsatellites[ms_id].set_read_dis_info(ms_read_repeat_length_info)

    def merge_reads_info(self):
        microsatellites_dict = {ms_id: {} for ms_id in self.microsatellites}
        for read_id, read in self.reads.items():
            # print("read_id",read_id)
            for ms_id in read.microsatellites:
                # print("ms_id", ms_id)
                microsatellites_dict[ms_id][read_id] = read
        self.reads = {}
        for ms_id, reads_info in microsatellites_dict.items():
            self.microsatellites[ms_id].set_reads_info(reads_info)

    def merge_muts_info(self, only_simple=True):
        # logger.info("\tMerge microsatellites infomation from different reads... ")
        microsatellites_dict_dis = {ms_id: {} for ms_id in self.microsatellites}
        microsatellites_dict_mut = {ms_id: {} for ms_id in self.microsatellites}
        for read_id, read in self.reads.items():
            stand = read.strand
            hap = read.hap
            for ms_id, ms_read_repeat_length in read.repeat_lengths.items():
                microsatellites_dict_dis[ms_id][read_id] = [ms_read_repeat_length, stand, hap]
            if not only_simple:
                for ms_id, ms_read_mut in read.mut_info.items():
                    microsatellites_dict_mut[ms_id][read_id] = ms_read_mut
        self.reads = {}
        for ms_id, ms_read_repeat_length_info in microsatellites_dict_dis.items():
            self.microsatellites[ms_id].set_read_dis_info(ms_read_repeat_length_info)
            if not only_simple:
                self.microsatellites[ms_id].set_muts_info(microsatellites_dict_mut[ms_id])
        # microsatellites_dict = {ms_id: {} for ms_id in self.microsatellites}
        # for read_id, read in self.reads.items():
        #     for ms_id, ms_read_mut in read.mut_info.items():
        #         microsatellites_dict[ms_id][read_id] = ms_read_mut
        # self.reads = {}  # release memory
        # for ms_id, reads_info in microsatellites_dict.items():
        #     self.microsatellites[ms_id].set_muts_info(reads_info)

    # def genotype_one_microsatellite_ccs_contig(self, microsatellite):
    #     # microsatellite.get_dis()
    #     microsatellite.one_hap_genotype()
    #     return microsatellite

    def genotype_one_microsatellite_ccs(self, microsatellite):
        # microsatellite.get_dis()

        return microsatellite

    # def genotype_microsatellite_ccs_contig(self):
    #     microsatellites = []
    #     for microsatellite in self.microsatellites.values():
    #         microsatellites.append(self.genotype_one_microsatellite_ccs_contig(microsatellite))
    #     self.microsatellites = {ms.ms_id: ms for ms in microsatellites}

    def call_variants(self):
        microsatellites = []
        for microsatellite in self.microsatellites.values():
            if self.paras["only_simple"]:
                microsatellite.call_micro()
            else:
                microsatellite.call_micro_and_other()
            # microsatellite.remove_noise()
            # microsatellite.ccs_genotype()
            microsatellites.append(microsatellite)
        self.microsatellites = {ms.ms_id: ms for ms in microsatellites}

    # def write_to_vcf_ccs_contig(self, file_output):
    #     # logger.info("\tWrite to vcf ... ")
    #     recs = []
    #     for ms_id in self.microsatellites_id:
    #         # print(ms_id)
    #         # print(self.microsatellites)
    #         ms = self.microsatellites[ms_id]
    #         vcfrec = file_output.new_record()
    #         # print("infoKey",vcfrec.info.keys())
    #         vcfrec.contig = ms.chrom
    #         # vcfrec.stop = pos + ms.repeat_times * len(ms.motif)
    #         vcfrec.pos = ms.mut_start
    #         vcfrec.ref = ms.ref_str
    #         vcfrec.alts = (ms.alt_str,) if ms.alt_str != "" else (".",)
    #         vcfrec.id = ms.ms_id
    #         vcfrec.stop = ms.mut_end
    #         vcfrec.info["ms_start"] = ms.start
    #         vcfrec.info["ms_end"] = ms.end
    #         vcfrec.info["motif"] = ms.repeat_unit
    #         vcfrec.info["repeat_times"] = ms.repeat_times
    #         vcfrec.info["motif_len"] = ms.repeat_unit_len
    #         vcfrec.info["ref_repeat_length"] = ms.repeat_len
    #         vcfrec.info["start_pre"] = ms.start_pre
    #         vcfrec.info["end_suf"] = ms.end_suf
    #         vcfrec.info["mut_start"] = ms.mut_start
    #         vcfrec.info["mut_end"] = ms.mut_end
    #
    #         vcfrec.info["query_repeat_length"] = ms.query_repeat_length
    #         vcfrec.info["dis_stat"] = str(ms.dis_stat)
    #         vcfrec.info["check"] = str(ms.check)
    #         vcfrec.info["check_stats"] = "|".join(ms.check_status)
    #         vcfrec.info["dis"] = "|".join(
    #             [str(key) + ":" + str(value) for key, value in ms.ms_dis.items()])
    #         # vcfrec.info["allele"] = ms.allele
    #         # if ms.check:
    #         vcfrec.info["var_type"] = ms.mut.var_type
    #         vcfrec.info["var_type_list"] = ms.mut.var_type_detail
    #         # print(ms.mut_type.var_prefix,ms.mut_type.var_ms,ms.mut_type.var_suffix)
    #         vcfrec.info["var_detail"] = ms.mut.var_detail_str
    #         # file_output.write(vcfrec)
    #         recs.append(vcfrec)
    #     return recs

    def write_to_vcf_call_variants_complex(self, file_output):
        recs = []
        for ms_id in self.microsatellites_id:
            # print(ms_id)
            # print(self.microsatellites)
            ms = self.microsatellites[ms_id]
            print(get_value("case"))
            print(ms.format_GT)
            print("ALT", ms.alt_str, ms.alt)
            print(ms.ref_str)

            vcfrec = file_output.new_record()
            vcfrec.contig = ms.chrom
            vcfrec.stop = ms.start + ms.repeat_times * ms.repeat_unit_len
            vcfrec.pos = ms.start
            vcfrec.ref = ms.ref_str
            # vcfrec.alts = (ms.alt,) if ms.alt != "" else ("N",)
            vcfrec.alts = ms.alt
            vcfrec.id = ms.ms_id
            vcfrec.stop = ms.end
            vcfrec.info["ms_start"] = ms.start
            vcfrec.info["ms_end"] = ms.end
            vcfrec.info["motif"] = ms.repeat_unit
            vcfrec.info["repeat_times"] = ms.repeat_times
            vcfrec.info["motif_len"] = ms.repeat_unit_len
            vcfrec.info["ref_repeat_length"] = ms.repeat_len
            vcfrec.info["query_repeat_length"] = ms.query_repeat_length
            vcfrec.info["depth"] = ms.depth
            vcfrec.info["dis_stat"] = str(ms.dis_stat)
            vcfrec.info["dis"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis.items()]
            )
            vcfrec.info["dis_hap0"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap0.items()]
            )
            vcfrec.info["dis_hap1"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap1.items()]
            )
            vcfrec.info["dis_hap2"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap2.items()]
            )
            vcfrec.info["dis_forward"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_forward.items()]
            )
            vcfrec.info["dis_reversed"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_reversed.items()]
            )
            vcfrec.info["Quality"] = "|".join(map(str, [ms.qual_ms, ms.qual_ms_hap1, ms.qual_ms_hap2]))
            # vcfrec.info["Distance"] = mscall.distance
            # vcfrec.qual = round(mscall.qual, 6)
            # vcfrec.info["Alleles"] = ms.alleles
            # print(mscall.format_DP)

            vcfrec.samples[get_value("case")]["GT"] = ms.format_GT
            vcfrec.samples[get_value("case")]["DP"] = ms.format_DP
            vcfrec.samples[get_value("case")]["QL"] = ms.format_QL
            vcfrec.samples[get_value("case")]["AL"] = ms.format_AL
            vcfrec.samples[get_value("case")].phased = ms.reads_phased

            # recs.append(vcfrec)
            recs.append(vcfrec)
        return recs

    def write_to_vcf_call_variants_micro(self, file_output):
        recs = []
        for ms_id in self.microsatellites_id:
            # print(ms_id)
            # print(self.microsatellites)
            ms = self.microsatellites[ms_id]
            # print(get_value("case"))
            # print(ms.format_GT)
            # print("ALT", ms.alt_str, ms.alt)
            # print(ms.ref_str)

            vcfrec = file_output.new_record()
            vcfrec.contig = ms.chrom
            vcfrec.stop = ms.start + ms.repeat_times * ms.repeat_unit_len
            vcfrec.pos = ms.start
            # vcfrec.ref = ms.ref_str_ms
            vcfrec.ref=str(ms.repeat_times)+"["+ms.repeat_unit+"]"
            # vcfrec.alts = (ms.alt,) if ms.alt != "" else ("N",)
            if ms.report_micro:
                vcfrec.alts = ms.alt_ms
            vcfrec.id = ms.ms_id
            vcfrec.stop = ms.end
            vcfrec.info["ms_start"] = ms.start
            vcfrec.info["ms_end"] = ms.end
            vcfrec.info["motif"] = ms.repeat_unit
            vcfrec.info["repeat_times"] = ms.repeat_times
            vcfrec.info["motif_len"] = ms.repeat_unit_len
            vcfrec.info["ref_repeat_length"] = ms.repeat_len
            vcfrec.info["query_repeat_length"] = ms.query_repeat_length
            vcfrec.info["depth"] = ms.depth
            vcfrec.info["dis_stat"] = str(ms.dis_stat)
            vcfrec.info["dis"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis.items()]
            )
            vcfrec.info["dis_hap0"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap0.items()]
            )
            vcfrec.info["dis_hap1"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap1.items()]
            )
            vcfrec.info["dis_hap2"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap2.items()]
            )
            vcfrec.info["dis_forward"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_forward.items()]
            )
            vcfrec.info["dis_reversed"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_reversed.items()]
            )
            vcfrec.info["Quality"] = "|".join(map(str, [ms.qual_ms, ms.qual_ms_hap1, ms.qual_ms_hap2]))
            # vcfrec.info["Distance"] = mscall.distance
            # vcfrec.qual = round(mscall.qual, 6)
            # vcfrec.info["Alleles"] = ms.alleles
            # print(mscall.format_DP)
            if ms.report_micro:
                vcfrec.samples[get_value("case")]["GT"] = ms.format_GT_ms
                vcfrec.samples[get_value("case")]["DP"] = ms.format_DP_ms
                vcfrec.samples[get_value("case")]["QL"] = ms.format_QL_ms
                vcfrec.samples[get_value("case")]["AL"] = ms.format_AL_ms
                vcfrec.samples[get_value("case")].phased = ms.reads_phased
            recs.append(vcfrec)
        return recs

    def write_to_vcf_call_variants_indel(self, file_output):
        recs = []
        for ms_id in self.microsatellites_id:
            # print(ms_id)
            # print(self.microsatellites)
            ms = self.microsatellites[ms_id]
            print(get_value("case"))
            print(ms.format_GT)
            print("ALT", ms.alt_str, ms.alt)
            print(ms.ref_str)
            vcfrec = file_output.new_record()
            vcfrec.contig = ms.chrom
            vcfrec.stop = ms.start + ms.repeat_times * ms.repeat_unit_len
            vcfrec.pos = ms.start
            vcfrec.ref = ms.ref_str
            # vcfrec.alts = (ms.alt,) if ms.alt != "" else ("N",)
            vcfrec.alts = ms.alt
            vcfrec.id = ms.ms_id
            vcfrec.stop = ms.end
            vcfrec.info["ms_start"] = ms.start
            vcfrec.info["ms_end"] = ms.end
            vcfrec.info["motif"] = ms.repeat_unit
            vcfrec.info["repeat_times"] = ms.repeat_times
            vcfrec.info["motif_len"] = ms.repeat_unit_len
            vcfrec.info["ref_repeat_length"] = ms.repeat_len
            vcfrec.info["query_repeat_length"] = ms.query_repeat_length
            vcfrec.info["depth"] = ms.depth
            vcfrec.info["dis_stat"] = str(ms.dis_stat)
            vcfrec.info["dis"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis.items()]
            )
            vcfrec.info["dis_hap0"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap0.items()]
            )
            vcfrec.info["dis_hap1"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap1.items()]
            )
            vcfrec.info["dis_hap2"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap2.items()]
            )
            vcfrec.info["dis_forward"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_forward.items()]
            )
            vcfrec.info["dis_reversed"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_reversed.items()]
            )
            vcfrec.info["Quality"] = "|".join(map(str, [ms.qual_ms, ms.qual_ms_hap1, ms.qual_ms_hap2]))
            # vcfrec.info["Distance"] = mscall.distance
            # vcfrec.qual = round(mscall.qual, 6)
            # vcfrec.info["Alleles"] = ms.alleles
            # print(mscall.format_DP)

            vcfrec.samples[get_value("case")]["GT"] = ms.format_GT
            vcfrec.samples[get_value("case")]["DP"] = ms.format_DP
            vcfrec.samples[get_value("case")]["QL"] = ms.format_QL
            vcfrec.samples[get_value("case")]["AL"] = ms.format_AL
            vcfrec.samples[get_value("case")].phased = ms.reads_phased

            # recs.append(vcfrec)
            recs.append(vcfrec)
        return recs

    def write_to_vcf_call_variants_snv(self, file_output):
        recs = []
        for ms_id in self.microsatellites_id:
            # print(ms_id)
            # print(self.microsatellites)
            ms = self.microsatellites[ms_id]
            print(get_value("case"))
            print(ms.format_GT)
            print("ALT", ms.alt_str, ms.alt)
            print(ms.ref_str)

            vcfrec = file_output.new_record()
            vcfrec.contig = ms.chrom
            vcfrec.stop = ms.start + ms.repeat_times * ms.repeat_unit_len
            vcfrec.pos = ms.start
            vcfrec.ref = ms.ref_str
            # vcfrec.alts = (ms.alt,) if ms.alt != "" else ("N",)
            vcfrec.alts = ms.alt
            vcfrec.id = ms.ms_id
            vcfrec.stop = ms.end
            vcfrec.info["ms_start"] = ms.start
            vcfrec.info["ms_end"] = ms.end
            vcfrec.info["motif"] = ms.repeat_unit
            vcfrec.info["repeat_times"] = ms.repeat_times
            vcfrec.info["motif_len"] = ms.repeat_unit_len
            vcfrec.info["ref_repeat_length"] = ms.repeat_len
            vcfrec.info["query_repeat_length"] = ms.query_repeat_length
            vcfrec.info["depth"] = ms.depth
            vcfrec.info["dis_stat"] = str(ms.dis_stat)
            vcfrec.info["dis"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis.items()]
            )
            vcfrec.info["dis_hap0"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap0.items()]
            )
            vcfrec.info["dis_hap1"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap1.items()]
            )
            vcfrec.info["dis_hap2"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap2.items()]
            )
            vcfrec.info["dis_forward"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_forward.items()]
            )
            vcfrec.info["dis_reversed"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_reversed.items()]
            )
            vcfrec.info["Quality"] = "|".join(map(str, [ms.qual_ms, ms.qual_ms_hap1, ms.qual_ms_hap2]))
            # vcfrec.info["Distance"] = mscall.distance
            # vcfrec.qual = round(mscall.qual, 6)
            # vcfrec.info["Alleles"] = ms.alleles
            # print(mscall.format_DP)

            vcfrec.samples[get_value("case")]["GT"] = ms.format_GT
            vcfrec.samples[get_value("case")]["DP"] = ms.format_DP
            vcfrec.samples[get_value("case")]["QL"] = ms.format_QL
            vcfrec.samples[get_value("case")]["AL"] = ms.format_AL
            vcfrec.samples[get_value("case")].phased = ms.reads_phased

            # recs.append(vcfrec)
            recs.append(vcfrec)
        return recs

    def write_to_vcf_call_variants(self, file_output):
        recs = []
        for ms_id in self.microsatellites_id:
            # print(ms_id)
            # print(self.microsatellites)
            ms = self.microsatellites[ms_id]
            print(get_value("case"))
            print(ms.format_GT)
            print("ALT", ms.alt_str, ms.alt)
            print(ms.ref_str)

            vcfrec = file_output.new_record()
            vcfrec.contig = ms.chrom
            vcfrec.stop = ms.start + ms.repeat_times * ms.repeat_unit_len
            vcfrec.pos = ms.start
            vcfrec.ref = ms.ref_str
            # vcfrec.alts = (ms.alt,) if ms.alt != "" else ("N",)
            vcfrec.alts = ms.alt
            vcfrec.id = ms.ms_id
            vcfrec.stop = ms.end
            vcfrec.info["ms_start"] = ms.start
            vcfrec.info["ms_end"] = ms.end
            vcfrec.info["motif"] = ms.repeat_unit
            vcfrec.info["repeat_times"] = ms.repeat_times
            vcfrec.info["motif_len"] = ms.repeat_unit_len
            vcfrec.info["ref_repeat_length"] = ms.repeat_len
            vcfrec.info["query_repeat_length"] = ms.query_repeat_length
            vcfrec.info["depth"] = ms.depth
            vcfrec.info["dis_stat"] = str(ms.dis_stat)
            vcfrec.info["dis"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis.items()]
            )
            vcfrec.info["dis_hap0"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap0.items()]
            )
            vcfrec.info["dis_hap1"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap1.items()]
            )
            vcfrec.info["dis_hap2"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap2.items()]
            )
            vcfrec.info["dis_forward"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_forward.items()]
            )
            vcfrec.info["dis_reversed"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_reversed.items()]
            )
            vcfrec.info["Quality"] = "|".join(map(str, [ms.qual_ms, ms.qual_ms_hap1, ms.qual_ms_hap2]))
            # vcfrec.info["Distance"] = mscall.distance
            # vcfrec.qual = round(mscall.qual, 6)
            # vcfrec.info["Alleles"] = ms.alleles
            # print(mscall.format_DP)

            vcfrec.samples[get_value("case")]["GT"] = ms.format_GT
            vcfrec.samples[get_value("case")]["DP"] = ms.format_DP
            vcfrec.samples[get_value("case")]["QL"] = ms.format_QL
            vcfrec.samples[get_value("case")]["AL"] = ms.format_AL
            vcfrec.samples[get_value("case")].phased = ms.reads_phased
            # recs.append(vcfrec)
            recs.append(vcfrec)
        return recs

    def write_to_vcf_pre_stat(self, file_output):
        # logger.info("\tWrite to vcf ... ")
        recs = []
        for ms_id in self.microsatellites_id:
            # print(ms_id)
            # print(self.microsatellites)
            ms = self.microsatellites[ms_id]
            vcfrec = file_output.new_record()
            # print("infoKey",vcfrec.info.keys())
            vcfrec.contig = ms.chrom
            # vcfrec.stop = pos + ms.repeat_times * len(ms.motif)
            vcfrec.pos = ms.start
            # vcfrec.ref = ms.ref_str
            vcfrec.ref = "."
            vcfrec.alts = (ms.alt_str,) if ms.alt_str != "" else (".",)
            vcfrec.id = ms.ms_id
            vcfrec.stop = ms.end
            vcfrec.info["ms_start"] = ms.start
            vcfrec.info["ms_end"] = ms.end
            vcfrec.info["motif"] = ms.repeat_unit
            vcfrec.info["repeat_times"] = ms.repeat_times
            vcfrec.info["motif_len"] = ms.repeat_unit_len
            vcfrec.info["ref_repeat_length"] = ms.repeat_len
            vcfrec.info["query_repeat_length"] = ms.query_repeat_length
            vcfrec.info["depth"] = ms.depth
            vcfrec.info["dis_stat"] = str(ms.dis_stat)
            vcfrec.info["dis"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis.items()]
            )
            vcfrec.info["dis_hap0"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap0.items()]
            )
            vcfrec.info["dis_hap1"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap1.items()]
            )
            vcfrec.info["dis_hap2"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap2.items()]
            )
            vcfrec.info["dis_forward"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_forward.items()]
            )
            vcfrec.info["dis_reversed"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_reversed.items()]
            )
            recs.append(vcfrec)
        return recs

    # def run_window_benchmark(self):
    #     """
    #     For CCS contig benchmark
    #     Returns:
    #     """
    #     self.init_microsatellites()  # 并行
    #     self.init_reads()  # 扫描read 确实其对应的 MS
    #     self.get_reads_info()  # 处理read 并行
    #     self.merge_muts_info()  # 合并read信息为MS信息
    #     self.genotype_microsatellite_ccs_contig()  # 变异检测 并行
    #     # self.write_to_vcf_ccs_contig(file_output)  # 一条一条写入

    def run_window_pre_stat(self):
        """
        For ngs/ccs pre_stat
        Returns:
        """
        self.init_microsatellites()  # 并行
        self.init_reads()  # 扫描read 确实其对应的 MS
        self.get_reads_dis()  # 处理read 并行
        self.merge_reads_repeat_length_distribution()  # 合并read信息为MS信息

    def run_window_call_variant(self):
        """
        For ngs/ccs variant calling
        Returns:

        """
        self.init_microsatellites(only_simple=self.paras["only_simple"])  # 并行
        self.init_reads()  # 扫描read 确实其对应的 MS
        self.get_reads_info(only_simple=self.paras["only_simple"])
        self.merge_muts_info(only_simple=self.paras["only_simple"])  # 合并read信息为MS信息
        self.call_variants()
