#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : benchmark_merge.py
# Author : Peng Jia
# Date   : 2020.07.16
# Email  : pengjia@stu.xjtu.edu.cn
# Description: Merge two haplotype calling result
=============================================================================="""
from src.global_dict import *
import os
import pysam


def benchmark_merge_init(args):
    """
        argument procress
        """
    paras = {}
    paras["hap1"] = args.hap1[0]
    paras["hap2"] = args.hap2[0]
    paras["output"] = args.output[0]
    paras["sample"] = args.sample[0]
    paras["debug"] = args.debug[0]
    error_stat = False
    if os.path.exists(paras["hap1"]):
        print("[INFO] The haplotype 1 microsatellite calling result is : '" + paras["hap1"] + "'.")
    else:
        print('[ERROR] The aplotype 1 microsatellite calling result '
              + paras["input"] + ' is not exist, please check again')
        error_stat = True
    if os.path.exists(paras["hap2"]):
        print("[INFO] The haplotype 2 microsatellite calling result is : '" + paras["hap2"] + "'.")
    else:
        print('[ERROR] The aplotype 2 microsatellite calling result '
              + paras["input"] + ' is not exist, please check again')
        error_stat = True
    if not os.path.exists(paras["output"]):
        print("[INFO] The output is : " + paras["output"] + ".")
    else:
        print(
            '[ERROR] The output ' + paras["output"] +
            ' is still exist! in case of overwrite files in this workspace, '
            'please check your script!')
        if not paras["debug"]:
            error_stat = True
    set_value("paras", paras)
    if error_stat: return False
    return True


def merge():
    paras = get_value("paras")
    hap1 = pysam.VariantFile(paras["hap1"])
    hap2 = pysam.VariantFile(paras["hap2"])
    sample = paras["sample"]
    mergerd_file = pysam.VariantFile(paras["output"], "w", header=hap1.header)
    # outputfile.header.add_line('##INFO=<ID=MutationType,Number=1,Type=String,Description="Mutation type">')
    mergerd_file.header.add_line('##INFO=<ID=dis1,Number=1,Type=String,Description='
                                 'Repeat length distribution of hap1>')
    mergerd_file.header.add_line('##INFO=<ID=dis2,Number=1,Type=String,Description='
                                 'Repeat length distribution of hap2>')
    mergerd_file.header.add_line('##INFO=<ID=var_type_list_hap1,Number=1,Type=String,Description='
                                 '"Variant type, Complex, SNP, DEL, INS">')
    mergerd_file.header.add_line('##INFO=<ID=var_detail_hap1,Number=1,Type=String,Description='
                                 '"Variant Details, rec1|rec2 rec1: pos:type:content">')
    mergerd_file.header.add_line('##INFO=<ID=var_type_list_hap2,Number=1,Type=String,Description='
                                 '"Variant type, Complex, SNP, DEL, INS">')
    mergerd_file.header.add_line('##INFO=<ID=var_detail_hap2,Number=1,Type=String,Description='
                                 '"Variant Details, rec1|rec2 rec1: pos:type:content">')
    mergerd_file.header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    mergerd_file.header.add_line('##FORMAT=<ID=AL,Number=1,Type=String,Description="Observed repeat length">')
    mergerd_file.header.add_line('##FORMAT=<ID=start_pre,Number=1,Type=String,'
                                 'Description="Start position of this haplotype analysis!">')
    mergerd_file.header.add_line('##FORMAT=<ID=end_suf,Number=1,Type=String,'
                                 'Description="End position of this haplotype analysis!">')
    mergerd_file.header.add_line('##FORMAT=<ID=mut_start,Number=1,Type=String,Description="Mutation start position">')
    mergerd_file.header.add_line('##FORMAT=<ID=mut_end,Number=1,Type=String,Description="Mutation end position">')
    mergerd_file.header.add_line('##FORMAT=<ID=mut_type,Number=1,Type=String,Description="Mutation type">')
    mergerd_file.header.add_sample(sample)
    for hap1_rec, hap2_rec in zip(hap1.fetch(), hap2.fetch()):
        if hap1_rec.info["check"] == "False" or hap2_rec.info["check"] == "False":
            continue
        merged_rec = mergerd_file.new_record()
        merged_rec.contig = hap1_rec.contig
        merged_rec.id = hap1_rec.id
        mut_start = min(hap1_rec.pos, hap2_rec.pos)
        mut_end = max(hap1_rec.info["mut_end"], hap2_rec.info["mut_end"])
        merged_rec.pos = mut_start
        merged_rec.stop = mut_end
        if "None" == hap1_rec.info["var_type"]:
            hap1_rec.alts = [hap1_rec.ref]
        if "None" == hap2_rec.info["var_type"]:
            hap2_rec.alts = [hap2_rec.ref]
        ref_prefix = ""
        hap1_prefix = ""
        hap2_prefix = ""
        ref_suffix = ""
        hap1_suffix = ""
        hap2_suffix = ""
        if hap1_rec.info["mut_start"] > mut_start:
            ref_prefix = hap2_rec.ref[0:hap1_rec.info["mut_start"] - mut_start]
            hap1_prefix = hap2_rec.ref[0:hap1_rec.info["mut_start"] - mut_start]
        if hap2_rec.info["mut_start"] > mut_start:
            # ref_prefix = hap1_rec.ref[0:hap2_rec.info["mut_start"] - mut_start]
            hap2_prefix = hap1_rec.ref[0:hap2_rec.info["mut_start"] - mut_start]
        if hap1_rec.info["mut_end"] < mut_end:
            ref_suffix = hap2_rec.ref[(hap1_rec.info["mut_end"] - mut_end):]
            hap1_suffix = hap2_rec.ref[(hap1_rec.info["mut_end"] - mut_end):]
        if hap2_rec.info["mut_end"] < mut_end:
            # ref_suffix=hap1_rec.ref[-(hap2_rec.info["mut_end"]-mut_end):]
            hap2_suffix = hap1_rec.ref[(hap2_rec.info["mut_end"] - mut_end):]
        merged_rec.ref = ref_prefix + hap1_rec.ref + ref_suffix
        hap1_alt = hap1_prefix + hap1_rec.alts[0] + hap1_suffix
        hap2_alt = hap2_prefix + hap2_rec.alts[0] + hap2_suffix
        alts = [merged_rec.ref]
        if hap1_alt not in alts:
            alts.append(hap1_alt)
        if hap2_alt not in alts:
            alts.append(hap2_alt)
        gt = [alts.index(hap1_alt), alts.index(hap2_alt)]
        if len(alts) == 1:
            alts.append(".")
        merged_rec.alts = alts[1:]
        merged_rec.info["ms_start"] = hap1_rec.info["ms_start"]
        merged_rec.info["ms_end"] = hap1_rec.info["ms_end"]
        merged_rec.info["motif"] = hap1_rec.info["motif"]
        merged_rec.info["repeat_times"] = hap1_rec.info["repeat_times"]
        merged_rec.info["motif_len"] = hap1_rec.info["motif_len"]
        merged_rec.info["ref_repeat_length"] = hap1_rec.info["ref_repeat_length"]
        merged_rec.info["start_pre"] = min([hap1_rec.info["start_pre"], hap2_rec.info["start_pre"]])
        merged_rec.info["end_suf"] = max([hap1_rec.info["end_suf"], hap2_rec.info["end_suf"]])
        merged_rec.info["mut_start"] = mut_start
        merged_rec.info["mut_end"] = mut_end
        merged_rec.info["dis_stat"] = "|".join([hap1_rec.info["dis_stat"], hap2_rec.info["dis_stat"]])
        check = True if hap1_rec.info["check"] == "True" and hap2_rec.info["check"] == "True" else False
        merged_rec.info["check"] = "|".join([str(check), hap1_rec.info["check"], hap2_rec.info["check"]])
        merged_rec.info["check_stats"] = "|".join([str(check), hap1_rec.info["check"], hap2_rec.info["check"]])
        merged_rec.info["dis1"] = hap1_rec.info["dis"]
        merged_rec.info["dis2"] = hap2_rec.info["dis"]
        merged_rec.info["var_type"] = "|".join([hap1_rec.info["var_type"], hap2_rec.info["var_type"]])

        merged_rec.info["var_type_list_hap1"] = hap1_rec.info["var_type_list"]
        merged_rec.info["var_type_list_hap2"] = hap2_rec.info["var_type_list"]
        merged_rec.info["var_detail"] = hap1_rec.info["var_detail"]
        merged_rec.info["var_detail"] = hap2_rec.info["var_detail"]
        al = list(map(str, [hap1_rec.info["query_repeat_length"], hap2_rec.info["query_repeat_length"]]))
        start_pre = list(map(str, [hap1_rec.info["start_pre"], hap2_rec.info["start_pre"]]))
        end_suf = list(map(str, [hap1_rec.info["end_suf"], hap2_rec.info["end_suf"]]))
        mut_start = list(map(str, [hap1_rec.info["mut_start"], hap2_rec.info["mut_start"]]))
        mut_end = list(map(str, [hap1_rec.info["mut_end"], hap2_rec.info["mut_end"]]))
        mut_type = list(map(str, [hap1_rec.info["var_type"], hap2_rec.info["var_type"]]))
        merged_rec.samples[sample]["GT"] = tuple(gt)
        merged_rec.samples[sample]["AL"] = "|".join(al)
        merged_rec.samples[sample]["start_pre"] = "|".join(start_pre)
        merged_rec.samples[sample]["end_suf"] = "|".join(end_suf)
        merged_rec.samples[sample]["mut_start"] = "|".join(mut_start)
        merged_rec.samples[sample]["mut_end"] = "|".join(mut_end)
        merged_rec.samples[sample]["mut_type"] = "|".join(mut_type)
        merged_rec.samples[sample].phased = True
        if merged_rec.info["check"]:
            mergerd_file.write(merged_rec)

    mergerd_file.close()
    pysam.tabix_index(paras["output"], force=True, preset="vcf")

    print("[INFO] Finished")


def benchmark_merge(parase):
    if not benchmark_merge_init(parase):
        return -1
    merge()
