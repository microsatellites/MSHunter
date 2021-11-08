#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : ccs.py
# Author : Peng Jia
# Date   : 2020.08.14
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""
from mshunter.units import *
from mshunter.pre_stat import pre_stat
from mshunter.estimate_error import estimate_error
from mshunter.call_variants import call_variants

import yaml


def load_model(paras):
    model_path = paras["output_model"]
    stream = open(model_path, "r")
    model = yaml.load(stream.read(), Loader=yaml.FullLoader)
    set_value("model", model)


def genotype_ccs(paras):
    df_microsatellites = load_microsatellites(paras)
    # pre_stat(df_microsatellites) # TODO run this step in release version
    # estimate_error() # TODO run this step in release version
    # load_model(paras)  # TODO ignore run pre_stat and estimate process when debug
    call_variants(df_microsatellites)
