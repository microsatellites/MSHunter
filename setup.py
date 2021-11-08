#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: msisensor_rna-rna
# Script : setup.py.py
# Author : Peng Jia
# Date   : 2021.03.29
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""
from setuptools import setup
import os
import sys
curpath = os.path.abspath(os.path.dirname(sys.argv[0]))
sys.path.append(os.path.dirname(curpath))
from mshunter.global_dict import *

global_init()
setup(
    name=get_value("tools_name"),
    version=get_value("tools_version"),
    description=get_value("description"),
    url="https://github.com/microsatellites/MSHunter",
    author=get_value("author"),
    author_email=get_value("email"),
    license='Custom License',
    keywords='Microsatellite, Multiplatform, NGS, CCS/HiFi ',
    packages=['mshunter'],
    install_requires=[
        "pysam>=0.17.0", "matplotlib>=3.4.3", "numpy>=1.21.4", "scikit-learn>=1.0.1 ", 'pandas>=1.0',
        "pyyaml", "multiprocess"
    ],
    python_requires='>=3.6',
    entry_points={'console_scripts': ['mshunter = mshunter.main:main']},
    # scripts={"mshunter/main.py"},
)
