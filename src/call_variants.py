#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : call_variants.py
# Author : Peng Jia
# Date   : 2020.09.03
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""

import os
import re
import collections
import pysam
import multiprocessing
from src.global_dict import *
from src.units import *
from src.Window import Window


def call_variant_write_vcf_init_snv(outputpath):
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
    outputfile.header.add_sample(get_value("case"))
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
    outputfile.header.add_line('##INFO=<ID=allele,Number=1,Type=Integer,Description="Allele number in this site">')
    outputfile.header.add_line('##INFO=<ID=Quality,Number=1,Type=String,Description="Variant Quality">')
    outputfile.header.add_line('##INFO=<ID=dis,Number=1,Type=String,Description='
                               'Distribution of repeat length>')
    outputfile.header.add_line('##INFO=<ID=depth,Number=1,Type=Integer,Description='
                               'Support Reads>')
    outputfile.header.add_line('##INFO=<ID=dis_hap0,Number=1,Type=String,Description='
                               'Distribution of repeat length from unphased reads>')
    outputfile.header.add_line('##INFO=<ID=dis_hap1,Number=1,Type=String,Description='
                               'Distribution of repeat length from hap 1>')
    outputfile.header.add_line('##INFO=<ID=dis_hap2,Number=1,Type=String,Description='
                               'Distribution of repeat length from hap 2>')
    outputfile.header.add_line('##INFO=<ID=dis_forward,Number=1,Type=String,Description='
                               'Distribution of repeat length forward reads>')
    outputfile.header.add_line('##INFO=<ID=dis_reversed,Number=1,Type=String,Description='
                               'Distribution of repeat length  reversed read>')
    outputfile.header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    outputfile.header.add_line('##FORMAT=<ID=DP,Number=1,Type=String,Description="Allele Depth">')
    outputfile.header.add_line('##FORMAT=<ID=AL,Number=1,Type=String,Description="Allele">')
    outputfile.header.add_line('##FORMAT=<ID=QL,Number=1,Type=String,Description="Allele Quality">')

    return outputfile


def call_variant_write_vcf_init_indel(outputpath):
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
    outputfile.header.add_sample(get_value("case"))
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
    outputfile.header.add_line('##INFO=<ID=allele,Number=1,Type=Integer,Description="Allele number in this site">')
    outputfile.header.add_line('##INFO=<ID=Quality,Number=1,Type=String,Description="Variant Quality">')
    outputfile.header.add_line('##INFO=<ID=dis,Number=1,Type=String,Description='
                               'Distribution of repeat length>')
    outputfile.header.add_line('##INFO=<ID=depth,Number=1,Type=Integer,Description='
                               'Support Reads>')
    outputfile.header.add_line('##INFO=<ID=dis_hap0,Number=1,Type=String,Description='
                               'Distribution of repeat length from unphased reads>')
    outputfile.header.add_line('##INFO=<ID=dis_hap1,Number=1,Type=String,Description='
                               'Distribution of repeat length from hap 1>')
    outputfile.header.add_line('##INFO=<ID=dis_hap2,Number=1,Type=String,Description='
                               'Distribution of repeat length from hap 2>')
    outputfile.header.add_line('##INFO=<ID=dis_forward,Number=1,Type=String,Description='
                               'Distribution of repeat length forward reads>')
    outputfile.header.add_line('##INFO=<ID=dis_reversed,Number=1,Type=String,Description='
                               'Distribution of repeat length  reversed read>')
    outputfile.header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    outputfile.header.add_line('##FORMAT=<ID=DP,Number=1,Type=String,Description="Allele Depth">')
    outputfile.header.add_line('##FORMAT=<ID=AL,Number=1,Type=String,Description="Allele">')
    outputfile.header.add_line('##FORMAT=<ID=QL,Number=1,Type=String,Description="Allele Quality">')
    return outputfile


def call_variant_write_vcf_init_complex(outputpath):
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
    outputfile.header.add_sample(get_value("case"))
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
    outputfile.header.add_line('##INFO=<ID=allele,Number=1,Type=Integer,Description="Allele number in this site">')
    outputfile.header.add_line('##INFO=<ID=Quality,Number=1,Type=String,Description="Variant Quality">')
    outputfile.header.add_line('##INFO=<ID=dis,Number=1,Type=String,Description='
                               'Distribution of repeat length>')
    outputfile.header.add_line('##INFO=<ID=depth,Number=1,Type=Integer,Description='
                               'Support Reads>')
    outputfile.header.add_line('##INFO=<ID=dis_hap0,Number=1,Type=String,Description='
                               'Distribution of repeat length from unphased reads>')
    outputfile.header.add_line('##INFO=<ID=dis_hap1,Number=1,Type=String,Description='
                               'Distribution of repeat length from hap 1>')
    outputfile.header.add_line('##INFO=<ID=dis_hap2,Number=1,Type=String,Description='
                               'Distribution of repeat length from hap 2>')
    outputfile.header.add_line('##INFO=<ID=dis_forward,Number=1,Type=String,Description='
                               'Distribution of repeat length forward reads>')
    outputfile.header.add_line('##INFO=<ID=dis_reversed,Number=1,Type=String,Description='
                               'Distribution of repeat length  reversed read>')
    outputfile.header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    outputfile.header.add_line('##FORMAT=<ID=DP,Number=1,Type=String,Description="Allele Depth">')
    outputfile.header.add_line('##FORMAT=<ID=AL,Number=1,Type=String,Description="Allele">')
    outputfile.header.add_line('##FORMAT=<ID=QL,Number=1,Type=String,Description="Allele Quality">')
    return outputfile


def call_variant_write_vcf_init(outputpath):
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
    outputfile.header.add_sample(get_value("case"))
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
    outputfile.header.add_line('##INFO=<ID=allele,Number=1,Type=Integer,Description="Allele number in this site">')
    outputfile.header.add_line('##INFO=<ID=Quality,Number=1,Type=String,Description="Variant Quality">')
    outputfile.header.add_line('##INFO=<ID=dis,Number=1,Type=String,Description='
                               'Distribution of repeat length>')
    outputfile.header.add_line('##INFO=<ID=depth,Number=1,Type=Integer,Description='
                               'Support Reads>')
    outputfile.header.add_line('##INFO=<ID=dis_hap0,Number=1,Type=String,Description='
                               'Distribution of repeat length from unphased reads>')
    outputfile.header.add_line('##INFO=<ID=dis_hap1,Number=1,Type=String,Description='
                               'Distribution of repeat length from hap 1>')
    outputfile.header.add_line('##INFO=<ID=dis_hap2,Number=1,Type=String,Description='
                               'Distribution of repeat length from hap 2>')
    outputfile.header.add_line('##INFO=<ID=dis_forward,Number=1,Type=String,Description='
                               'Distribution of repeat length forward reads>')
    outputfile.header.add_line('##INFO=<ID=dis_reversed,Number=1,Type=String,Description='
                               'Distribution of repeat length  reversed read>')
    outputfile.header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    outputfile.header.add_line('##FORMAT=<ID=DP,Number=1,Type=String,Description="Allele Depth">')
    outputfile.header.add_line('##FORMAT=<ID=AL,Number=1,Type=String,Description="Allele">')
    outputfile.header.add_line('##FORMAT=<ID=QL,Number=1,Type=String,Description="Allele Quality">')

    return outputfile


def call_variant_write_vcf_close(outputfile):
    args = get_value("paras")
    if args["only_simple"]:
        outputfile["micro"].close()
        pysam.tabix_index(args["output_micro"], force=True, preset="vcf")

    else:
        outputfile["indel"].close()
        pysam.tabix_index(args["output_indel"], force=True, preset="vcf")
        outputfile["snv"].close()
        pysam.tabix_index(args["output_snv"], force=True, preset="vcf")
        outputfile["complex"].close()
        pysam.tabix_index(args["output_complex"], force=True, preset="vcf")
        outputfile["micro"].close()
        pysam.tabix_index(args["output_micro"], force=True, preset="vcf")

    # outputfile.close()
    # pysam.tabix_index(get_value("paras")["output_pre"], force=True, preset="vcf")


def run_one_window(win_info):
    window = Window(win_info)
    window.run_window_call_variant()
    return window


def run_window_mul(windows, args, file_output):
    contig = windows[0][0]["chr"]
    start = windows[0][0]["pos"]
    end = windows[-1][-1]["pos"]
    num = 0
    for win in windows:
        num += len(win)
    logger.info("--------------------------------------------------------------------------------")
    logger.info("Processing " + contig + ":" + str(start) + "-" + str(end))
    logger.info("No. of Microsatellites in window: " + str(num))

    pool = multiprocessing.Pool(processes=args["threads"])
    windows = pool.map(run_one_window, windows)
    pool.close()
    pool.join()
    # windows_res=[]
    # for win in windows:
    #     windows_res.append(run_one_window(win))
    # windows=windows_res
    # win_recs=[]
    if args["only_simple"]:
        for win in windows:
            for rec in win.write_to_vcf_call_variants_micro(file_output["micro"]):
                file_output["micro"].write(rec)
            # for rec in win.write_to_vcf_call_variants_indel(file_output["indel"]):
            #     file_output["indel"].write(rec)
    else:
        for win in windows:
            for rec in win.write_to_vcf_call_variants_micro(file_output["micro"]):
                file_output["micro"].write(rec)
            # for rec in win.write_to_vcf_call_variants_indel(file_output["indel"]):
            #     file_output["indel"].write(rec)
            # for rec in win.write_to_vcf_call_variants_snv(file_output["snv"]):
            #     file_output["snv"].write(rec)
            # for rec in win.write_to_vcf_call_variants_indel(file_output["complex"]):
            #     file_output["complex"].write(rec)

    logger.info("Total Microsatellites: " + str(args["ms_num"]))
    logger.info("Finished Microsatellites: " + str(args["current_num"]) +
                " (" + str(round(args["current_num"] / args["ms_num"] * 100, 2)) + "%)")


def call_variants(df_microsatellites):
    args = get_value("paras")
    output_file_dict = {}
    if args["only_simple"]:
        output_micro = call_variant_write_vcf_init_indel(args["output_micro"])
        output_file_dict["micro"] = output_micro
    else:
        output_micro = call_variant_write_vcf_init_indel(args["output_micro"])
        output_indel = call_variant_write_vcf_init_indel(args["output_indel"])
        output_snv = call_variant_write_vcf_init_snv(args["output_snv"])
        output_complex = call_variant_write_vcf_init_complex(args["output_complex"])
        output_file_dict["indel"] = output_indel
        output_file_dict["snv"] = output_snv
        output_file_dict["complex"] = output_complex
        output_file_dict["micro"] = output_micro
    contigs_info = get_value("contigs_info")
    if args["debug"]:
        locis_num = 10000
        df_microsatellites = df_microsatellites.iloc[106000:locis_num + 106000, :]
    args["ms_num"] = len(df_microsatellites)
    total_current_num = 0
    for contig, contig_len in contigs_info.items():
        logger.info("--------------------------------------------------------------------------------")
        logger.info("Call variant: Processing " + contig + "...")
        this_contig_microsatellite = df_microsatellites[df_microsatellites["chr"] == contig].sort_values("pos")
        window_ms = []
        ms_num = 0
        win_num = 0
        window_sub = []
        for ms_id, info in this_contig_microsatellite.iterrows():
            ms_num += 1
            total_current_num += 1
            info["prefix_len"] = args["prefix_len"]
            info["suffix_len"] = args["suffix_len"]
            info["reference"] = args["reference"]
            window_sub.append(info)
            if ms_num % (args["batch"]) == 0:
                window_ms.append(window_sub)
                window_sub = []
                win_num += 1
                if win_num % args["threads"] == 0:
                    args["current_num"] = total_current_num
                    run_window_mul(window_ms, args, file_output=output_file_dict)
                    window_ms = []
        if len(window_ms) > 0:
            ms_num = 0
            for win in window_ms:
                ms_num += len(win)
            item_num = ms_num // args["threads"] + 1
            window_ms_tmp = []
            window_sub = []
            num = 0
            for win in window_ms:
                for ms in win:
                    num += 1
                    window_sub.append(ms)
                    if num % item_num == 0:
                        window_ms_tmp.append(window_sub)
                        window_sub = []
            if len(window_sub) > 0:
                window_ms_tmp.append(window_sub)
            del window_sub, window_ms
            run_window_mul(window_ms_tmp, args, file_output=output_file_dict)
            del window_ms_tmp
    call_variant_write_vcf_close(output_file_dict)
