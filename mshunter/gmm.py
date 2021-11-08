#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : gmm.py
# Author : Peng Jia
# Date   : 2020.10.26
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""

import numpy as np
from sklearn import mixture
from collections import Counter


# from scipy import stats
# dis = {8: 13}


def get_more_times(repeat_list):
    """
    Args:
        repeat_list (list): repeat length of microsatellite site, such as [10,10,9,9,9,9,8,5,5,5,5,5]
    Returns: central value of the input list

    """
    counts = Counter(repeat_list)
    f = sorted(zip(counts.values(), counts.keys()))
    # print(f)
    if len(f) == 1:
        return f[-1][-1]
    else:
        # print(abs(f[-1][0] - f[-2][0]) / len(repeat_list))
        if abs(f[-1][0] - f[-2][0]) / len(repeat_list) < 0.3:  # TODO add on command input
            # print("ldldldl")
            return int(round(np.mean(repeat_list)))
        else:
            return f[-1][-1]


def get_repeat_gmm(dis, target=1):
    """
    Args:
        dis (dict):distribution of microsatellite repeat length
        target (int): the haplotype num 1/2
    Returns:
        the eliminated repeat length of microsatellite

    """
    support = sum(dis.values())
    cluster = len(dis)
    if cluster == 1:
        genotype = list(dis.keys()) if target == 1 else list(dis.keys()) * 2
        qual = support
        return {"genotype": genotype, "qual": qual}
    repeat_list = []
    for k, v in dis.items():
        repeat_list.extend([k] * v)
    repeat_list = np.array(repeat_list).reshape(-1, 1)
    dpgmm = mixture.BayesianGaussianMixture(n_components=cluster,
                                            covariance_type='full',
                                            tol=0.0001,
                                            max_iter=400, ).fit(repeat_list)
    pre_dis = {}
    pre_num = {}
    # print(dpgmm.means_)
    for k, v in dis.items():
        k_pre = dpgmm.predict(np.array([[k]]))
        # print(k, k_pre)
        if k_pre[0] not in pre_dis:
            pre_dis[k_pre[0]] = []
            pre_num[k_pre[0]] = 0
        pre_dis[k_pre[0]].extend([k] * v)
        pre_num[k_pre[0]] += v
    m = sorted(pre_num.keys(), key=(lambda x: pre_num[x]))
    if target == 1:
        # genotype = [int(round(np.mean(pre_dis[m[-1]])))]
        # genotype = [get_more_times(pre_dis[m[-1]])]
        if len(m) > 1:
            genotype = [get_more_times(pre_dis[m[-1]] + pre_dis[m[-2]])]
            qual = (pre_num[m[-1]] + pre_num[m[-2]]) / (1 + np.std(pre_dis[m[-1]]))
        else:
            genotype = [get_more_times(pre_dis[m[-1]])]
            qual = (pre_num[m[-1]]) / (1 + np.std(pre_dis[m[-1]]))
    elif target == 2:
        if pre_num[m[-1]] > support * 0.7:
            # genotype = [int(round(np.mean(pre_dis[m[-1]])))] * 2
            genotype = [get_more_times(pre_dis[m[-1]])] * 2
            qual = support / (1 + np.std(pre_dis[m[-1]]))
        else:
            # hap1 = int(round(np.mean(pre_dis[m[-1]])))
            hap1 = get_more_times(pre_dis[m[-1]])
            qual1 = pre_num[m[-1]] / (1 + np.std(pre_dis[m[-1]]))
            # hap2 = int(round(np.mean(pre_dis[m[-2]])))
            hap2 = get_more_times(pre_dis[m[-2]])
            qual2 = pre_num[m[-2]] / (1 + np.std(pre_dis[m[-2]]))
            genotype = [hap1, hap2]
            qual = (qual1 + qual2) / 2
    return {"genotype": genotype, "qual": qual}


if __name__ == "__main__":
    dis = {8: 13, 9: 2, 7: 1, 10: 18}
    dis = {23: 1, 21: 3, 22: 3}
    dis = {27: 7, 26: 2, 30: 1, 28: 6, 29: 4}
    dis = {19: 8, 21: 5, 20: 5}
    print(get_repeat_gmm(dis, target=1))
