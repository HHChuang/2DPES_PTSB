#!/usr/bin/env python
#################################################################
# Program:                                                      #
# 2018/08/22, Grace, Searching the path from near the shoulder  #
# (i.e. shift away from TS1) of a bifurcation reaction to TS2   #
# code: genModIRCstruc.f90                                      #
# 2019/04/16, Grace, Rewrite source code from Fortran to Python,#
# because of I/O issue cause the low efficiency.                #
#                                                               #
# Reference:                                                    #
#   JCP 1977, 66, 2153. 'The intrinsic reaction coordinate. An  #
#   ab initio calculation for HNC->HCN and H-+CH4-> CH4+H_.     #
#                                                               #
# Input:                                                        #
#       $1 = * .out (gaussian output file w/ force and          #
#       freq=hpmodes)                                           #
#       $2 = *.com (gaussian input file)                        #
# Output:                                                       #
#       *.inp (user defined name = $2)                          #
#################################################################

import numpy as np
import re  # regular
import sys


def main():



    # Main program
main()
