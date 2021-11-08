#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : main.py
# Author : Peng Jia
# Date   : 2020.07.13
# Email  : pengjia@stu.xjtu.edu.cn
# Description: main function and arguments processing
=============================================================================="""
import argparse
import os
import sys

curpath = os.path.abspath(os.path.dirname(sys.argv[0]))
sys.path.append(os.path.dirname(curpath))

from mshunter.genotype import *
from mshunter.qc import *
from mshunter.paras import args_process

logger.info(" ".join(sys.argv))


def main():
    """
    Main function.
    :return:
    """
    global_init()
    arg = args_process()

    if arg:
        parase = arg.parse_args()
        if parase.command == "genotype":
            genotype(parase)
            # genotype_ngs(parase)
        if parase.command == "qc":
            qc(parase)

        # if parase.command == "ngs":
        #     # genotype(parase)
        #     genotype_ngs(parase)


if __name__ == "__main__":
    main()
