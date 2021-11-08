#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : call.py
# Author : Peng Jia
# Date   : 2020.07.13
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""

from mshunter.units import *
from mshunter.global_dict import *
import pysam
import multiprocessing


class MSCall:
    chr_id = ""
    pos = ""
    ref = ""
    info = {}
    dis = {}
    dis_norm = {}
    minAllele = 1
    maxAllele = 100
    modelStat = False
    model = {}
    distance = 0
    distance_dict = {}
    firstAllels = ""
    secondAllels = ""
    qual = 0
    filter = ""
    precision = "NoReadSpan"  # LowCov, Fuzzy, High
    format_GT = (0, 0)
    format_AL = "/".join(["0", "0"])
    format_DP = "/".join(["0", "0", "0"])
    format_QL = "/".join(["0", "0", "0"])
    alleles = "0:0"
    alt = (".",)
    first_two_distance = 0
    first_two_distance_hap1 = 0
    first_two_distance_hap2 = 0
    qual_hap1 = 0
    qual_hap2 = 0
    dis_hap1_normal = {}
    dis_hap2_normal = {}
    distance_hap1 = 0
    distance_hap2 = 0

    def __init__(self, chr_id, pos, info, sample):
        self.support_reads_hap1 = 0
        self.support_reads_hap2 = 0
        self.support_reads = 0
        self.info = info
        self.chr_id = chr_id
        self.pos = pos
        self.sample = sample
        self.disStat = False if info["disStat"] == "False" else True
        self.ref = str(info["repeatTimes"]) + "[" + info["motif"] + "]"
        self.phased = True if self.info["Phased"] == "True" else False

    def getDisdistance2(self, dict1, dict2):
        dictkey = list(dict1.keys()) + list(dict2.keys())
        for key in dictkey:
            if key not in dict1:
                dict1[key] = 0
            if key not in dict2:
                dict2[key] = 0
        sum = 0
        for key in dictkey:
            err = dict1[key] - dict2[key]
            sum += (err * err)
        return round(np.sqrt(sum), 6)

    def getDisdistance(self, dict1, dict2):
        dictkey = list(dict1.keys()) + list(dict2.keys())
        list1 = []
        list2 = []
        for key in dictkey:
            if key not in dict1:
                list1.append(0)
            else:
                list1.append(dict1[key])
            if key not in dict2:
                list2.append(0)
            else:
                list2.append(dict2[key])
        return np.linalg.norm(np.array(list1) - np.array(list2))

    def modelpre(self):
        if not self.disStat:
            return
        else:
            model = get_value("model")
            motif = self.info["motif"]
            if motif not in model:
                # self.modelStat = False
                return
            self.modelStat = True
            (self.support_reads, self.support_reads_hap1, self.support_reads_hap2) = list(
                map(int, self.info["support_reads"].split("|")))
            maxRepeat = model[motif]["maxRepeat"]

            if self.phased:
                dis_normal = {}
                this_dis = [list(map(int, i.split("-"))) for i in self.info["dis"].split("|")[1].split(":")]
                for i in this_dis:
                    dis_normal[i[0]] = i[1] / self.support_reads_hap1
                self.dis_hap1_normal = dis_normal
                # print(self.info)
                dis_normal = {}
                this_dis = [list(map(int, i.split("-"))) for i in self.info["dis"].split("|")[2].split(":")]
                for i in this_dis:
                    dis_normal[i[0]] = i[1] / self.support_reads_hap2
                self.dis_hap2_normal = dis_normal

                self.minAllele = max([min(self.dis_hap1_normal.keys()) - 1, min(self.dis_hap2_normal.keys()) - 1, 1])
                self.maxAllele = min(
                    [max(self.dis_hap2_normal.keys()) + 1, max(self.dis_hap2_normal.keys()) + 1, maxRepeat])
                model_id_list = []
                for first in range(self.minAllele, self.maxAllele + 1):
                    model_id_list.append(first * 1000 + first)
                thismodel = {}
                for model_id in model_id_list:
                    thismodel[model_id] = model[motif]["maxture"][model_id]
                self.model = thismodel
            else:
                dis_normal = {}
                this_dis = [list(map(int, i.split("-"))) for i in self.info["dis"].split("|")[0].split(":")]
                for i in this_dis:
                    dis_normal[i[0]] = i[1] / self.support_reads
                self.dis_norm = dis_normal
                self.minAllele = max([min(dis_normal.keys()) - 1, 1])
                self.maxAllele = min([max(dis_normal.keys()) + 1, maxRepeat])
                model_id_list = []
                for first in range(self.minAllele, self.maxAllele + 1):
                    for second in range(self.minAllele, self.maxAllele + 1):
                        if first <= second:
                            model_id_list.append(first * 1000 + second)
                thismodel = {}
                for model_id in model_id_list:
                    thismodel[model_id] = model[motif]["maxture"][model_id]
                self.model = thismodel

    def get_allele_from_hap(self):

        distance_dict = {}
        for model_id in self.model:
            distance_dict[model_id] = get_disdistance(self.dis_hap1_normal, self.model[model_id])
        distance_tuple = sorted(distance_dict.items(), key=lambda kv: (kv[1], kv[0]))
        first_1 = distance_tuple[0][0] // 1000
        first_two_distance_ratio_hap1 = (distance_tuple[1][1] - distance_tuple[0][1]) / (
                distance_tuple[0][1] + 0.000001)
        self.first_two_distance_hap1 = distance_tuple[1][1] - distance_tuple[0][1]
        self.distance_hap1 = distance_tuple[0][1]
        self.qual_hap1 = first_two_distance_ratio_hap1 * self.support_reads_hap1

        distance_dict = {}
        for model_id in self.model:
            distance_dict[model_id] = getDisdistance(self.dis_hap2_normal, self.model[model_id])
        self.model = {}
        distance_tuple = sorted(distance_dict.items(), key=lambda kv: (kv[1], kv[0]))
        first_2 = distance_tuple[0][0] // 1000
        first_two_distance_ratio_hap2 = (distance_tuple[1][1] - distance_tuple[0][1]) / (
                distance_tuple[0][1] + 0.000001)
        self.first_two_distance_hap2 = distance_tuple[1][1] - distance_tuple[0][1]
        self.distance_hap2 = distance_tuple[0][1]
        self.qual_hap2 = first_two_distance_ratio_hap2 * self.support_reads_hap1

        if self.info["lowSupport"] == "True":
            self.precision = "LowCov"
            self.filter = "LowCov"
        elif first_two_distance_ratio_hap1 < 0.3:
            self.precision = "Fuzzy"
            self.filter = "Fuzzy"
        else:
            self.precision = "High"
            self.filter = "PASS"
        self.firstAllels = first_1
        self.secondAllels = first_2
        self.qual = (self.qual_hap1 + self.qual_hap2) * 0.5

        if first_1 == first_2:
            if first_1 == self.info["repeatTimes"]:
                self.format_GT = (0, 0)
                # self.alt = (".")
            else:
                self.format_GT = (1, 1)
                self.alt = (str(first_1) + "[" + self.info["motif"] + "]",)
        else:
            if self.info["repeatTimes"] in [first_1, first_2]:
                self.format_GT = (0, 1)
                if first_1 == self.info["repeatTimes"]:
                    self.alt = (str(first_2) + "[" + self.info["motif"] + "]",)
                else:
                    self.alt = (str(first_1) + "[" + self.info["motif"] + "]",)
            else:
                self.format_GT = (1, 2)
                if first_1 < first_2:
                    self.alt = (str(first_1) + "[" + self.info["motif"] + "]",
                                str(first_2) + "[" + self.info["motif"] + "]")
                else:
                    self.alt = (str(first_2) + "[" + self.info["motif"] + "]",
                                str(first_1) + "[" + self.info["motif"] + "]")
        self.format_AL = "/".join(list(map(str, [first_1, first_2])))
        self.format_DP = "/".join(
            list(map(str, [self.support_reads, self.support_reads_hap1, self.support_reads_hap2])))
        self.format_QL = "/".join(list(map(str, [self.qual, self.qual_hap1, self.qual_hap2])))

    def get_alleles_from_diploid(self):
        distance_dict = {}
        for model_id in self.model:
            distance_dict[model_id] = get_disdistance(self.dis_norm, self.model[model_id])
        self.model = {}
        distance_tuple = sorted(distance_dict.items(), key=lambda kv: (kv[1], kv[0]))
        first_1 = distance_tuple[0][0] // 1000
        first_2 = distance_tuple[0][0] % 1000
        firsttwoDistanceRatio = (distance_tuple[1][1] - distance_tuple[0][1]) / (distance_tuple[0][1] + 0.000001)
        self.first_two_distance = distance_tuple[1][1] - distance_tuple[0][1]
        self.distance = distance_tuple[0][1]
        qual = firsttwoDistanceRatio * self.support_reads
        if self.info["lowSupport"] == "True":
            self.precision = "LowCov"
            self.filter = "LowCov"
        elif firsttwoDistanceRatio < 0.3:
            self.precision = "Fuzzy"
            self.filter = "Fuzzy"
        else:
            self.precision = "High"
            self.filter = "PASS"
        self.firstAllels = first_1
        self.secondAllels = first_2
        self.qual = qual
        if first_1 == first_2:
            if first_1 == self.info["repeatTimes"]:
                self.format_GT = (0, 0)
                # self.alt = (".")
            else:
                self.format_GT = (1, 1)
                self.alt = (str(first_1) + "[" + self.info["motif"] + "]",)
        else:
            if self.info["repeatTimes"] in [first_1, first_2]:
                self.format_GT = (0, 1)
                if first_1 == self.info["repeatTimes"]:
                    self.alt = (str(first_2) + "[" + self.info["motif"] + "]",)
                else:
                    self.alt = (str(first_1) + "[" + self.info["motif"] + "]",)
            else:
                self.format_GT = (1, 2)
                if first_1 < first_2:
                    self.alt = (str(first_1) + "[" + self.info["motif"] + "]",
                                str(first_2) + "[" + self.info["motif"] + "]")
                else:
                    self.alt = (str(first_2) + "[" + self.info["motif"] + "]",
                                str(first_1) + "[" + self.info["motif"] + "]")
        if first_1 <= first_2:
            self.format_AL = "/".join(list(map(str, [first_1, first_2])))
        else:
            self.format_AL = "/".join(list(map(str, [first_2, first_1])))
        # self.format_AL = ":".join(list(map(str, [first_1, first_2])))
        self.format_DP = "/".join(
            list(map(str, [self.support_reads, 0, 0])))
        self.format_QL = "/".join(list(map(str, [self.qual, 0, 0])))

    def mscall(self):
        qual = 0
        if not self.disStat:
            self.qual = qual
            self.precision = "NoReadSpan"
            self.filter = "NoReadSpan"
            return
        if (len(self.model) < 2) or (not self.modelStat):
            self.qual = qual
            self.precision = "NoAvailableModel"
            self.filter = "NoAvailableModel"
            return

        if self.phased:
            self.get_allele_from_hap()
        else:
            self.get_alleles_from_diploid()


def call_one_ms(msCall):
    msCall.modelpre()
    msCall.mscall()
    return msCall


def multicallMS(mscall_list, outputfile, thread=4):
    # pool = multiprocessing.Pool(thread)
    # result_list = pool.map(call_one_ms, mscall_list)
    # pool.close()
    # pool.join()
    result_list = []
    for i in mscall_list:
        result_list.append(call_one_ms(i))

    write_call2vcf(result_list, outputfile)
    return result_list


def write_call2vcf_init():
    paras = get_value("paras")
    inputpath = paras["output_dis"]
    outputpath = paras["output_call"]
    inputfile = pysam.VariantFile(inputpath, )
    outputfile = pysam.VariantFile(outputpath, "w", header=inputfile.header)
    outputfile.header.add_line(
        '##INFO=<ID=FirstAlleles,Number=1,Type=String,Description="The first allele type of this point">')
    outputfile.header.add_line(
        '##INFO=<ID=SecondAlleles,Number=1,Type=String,Description="The second allele type of this point">')
    outputfile.header.add_line('##INFO=<ID=Quality,Number=1,Type=Float,Description="Genotype quality">')
    outputfile.header.add_line(
        '##INFO=<ID=Distance,Number=1,Type=Float,Description="Distance between two distribution.">')
    outputfile.header.add_line(
        '##INFO=<ID=FirstTwoDistance,Number=1,Type=String,Description="Distance between two distribution.">')
    outputfile.header.add_line('##INFO=<ID=Precision,Number=1,Type=String,Description="Genotype quality level">')
    outputfile.header.add_line('##INFO=<ID=Alleles,Number=1,Type=String,Description="Alleles">')
    outputfile.header.add_line('##INFO=<ID=Mutation,Number=1,Type=String,Description="Mutation">')
    outputfile.header.add_line('##INFO=<ID=MutationType,Number=1,Type=String,Description="Mutation type">')
    # outputfile.header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    outputfile.header.add_line('##FORMAT=<ID=DP,Number=1,Type=String,Description="Allele Depth">')
    outputfile.header.add_line('##FORMAT=<ID=AL,Number=1,Type=String,Description="Allele">')
    outputfile.header.add_line('##FORMAT=<ID=QL,Number=1,Type=String,Description="Allele Quality">')
    outputfile.header.add_line('##INFO=<ID=Type,Number=1,Type=String,Description="Mutation type">')
    return outputfile


def write_call2vcf(mscall_list, outputfile):
    for mscall in mscall_list:
        vcfrec = outputfile.new_record()
        vcfrec.contig = mscall.chr_id
        vcfrec.pos = mscall.pos
        vcfrec.ref = mscall.ref
        vcfrec.alts = mscall.alt
        vcfrec.info["chrom"] = mscall.chr_id
        vcfrec.info["pos"] = mscall.info["pos"]
        vcfrec.info["Start"] = mscall.info["Start"]
        vcfrec.info["End"] = mscall.info["End"]
        vcfrec.info["motif"] = mscall.info["motif"]
        vcfrec.info["repeatTimes"] = mscall.info["repeatTimes"]
        vcfrec.info["prefix"] = mscall.info["prefix"]
        vcfrec.info["suffix"] = mscall.info["suffix"]
        vcfrec.info["depth"] = mscall.info["depth"]
        vcfrec.info["support_reads"] = mscall.info["support_reads"]
        vcfrec.info["dis"] = mscall.info["dis"]
        vcfrec.info["proD"] = mscall.info["proD"]
        vcfrec.info["proI"] = mscall.info["proI"]
        vcfrec.info["lowSupport"] = str(mscall.info["lowSupport"])
        vcfrec.info["disStat"] = str(mscall.info["disStat"])
        vcfrec.info["Precision"] = mscall.precision
        vcfrec.info["Quality"] = round(mscall.qual, 6)
        if not mscall.phased:
            vcfrec.info["FirstTwoDistance"] = str(mscall.first_two_distance)
        else:
            vcfrec.info["FirstTwoDistance"] = "|".join(
                [str(mscall.first_two_distance_hap1), str(mscall.first_two_distance_hap2)])
        vcfrec.info["Distance"] = mscall.distance
        vcfrec.qual = round(mscall.qual, 6)
        vcfrec.info["Alleles"] = mscall.alleles
        # print(mscall.format_DP)
        vcfrec.samples[get_value("case")]["GT"] = mscall.format_GT
        vcfrec.samples[get_value("case")]["DP"] = mscall.format_DP
        vcfrec.samples[get_value("case")]["QL"] = mscall.format_QL
        vcfrec.samples[get_value("case")]["AL"] = mscall.format_AL
        vcfrec.samples[get_value("case")].phased = mscall.phased
        outputfile.write(vcfrec)
    return outputfile


def write_call2vcf_close(outputfile):
    outputfile.close()
    pysam.tabix_index(get_value("paras")["output_call"], force=True, preset="vcf")


def call():
    paras = get_value("paras")
    thread = paras["threads"]
    batch = paras["batch"]
    disvcfpath = paras["output_dis"]
    vcffile = pysam.VariantFile(disvcfpath, "rb")
    curentMSNum = 0
    tmpWindow = []
    ms_number = get_value("ms_number")
    outputfile = write_call2vcf_init()
    sampleID = get_value("case")
    for rec in vcffile.fetch():
        msCall = MSCall(chr_id=rec.chrom, pos=rec.pos, info=dict(rec.info), sample=sampleID)
        tmpWindow.append(msCall)
        curentMSNum += 1
        # if curentMSNum > 10000 and paras["debug"]:
        #     break
        if curentMSNum % (batch * thread) == 0:
            print("[INFO] MS calling: Total", ms_number, "microsatelite, processing:", curentMSNum - batch * thread + 1,
                  "-", curentMSNum, "(" + str(round(curentMSNum / (0.1 + ms_number) * 100, 2)) + "%)")
            multicallMS(tmpWindow, outputfile=outputfile, thread=thread)
            tmpWindow = []
    multicallMS(tmpWindow, outputfile=outputfile, thread=thread)
    print("[INFO] MS calling: Total", ms_number, "microsatelite, finish all!")
