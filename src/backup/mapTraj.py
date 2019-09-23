#!/usr/bin/env python3
#########################################################
#                                                       #
# Input :                                               #
#       $1 = numerical PES                              #
#           format:                                     #
#                      # of atom                        #
#                       coord_x coord_y E               #
#                       #atom   x   y   z               #
#                                                       #
#       $2 = coord. of important points                 #
#           format:                                     #
#           $(name of point) $(x) $(y) $(z)             #
#                                                       #
# History:                                              #
# 2019/03/18, Grace, rewrite this code from Fortran into#
#   python code, because of the efficiency problem.     #
#########################################################

import sys
import os
import numpy as np

PATH_iniRP1 = '/Users/Grace/Google_Drive/Code/GitHub/ModifyIRC/run/H3CO/QCISD/Traj/RP1.reorder.rot./'
PATH_iniRP2 = '/Users/Grace/Google_Drive/Code/GitHub/ModifyIRC/run/H3CO/QCISD/Traj/RP2.reorder.rot./'
PATH_finRP1 = '/Users/Grace/Google_Drive/Code/GitHub/ModifyIRC/run/H3CO/QCISD/Traj/RP1.reorder.rot.coord/'
PATH_finRP2 = '/Users/Grace/Google_Drive/Code/GitHub/ModifyIRC/run/H3CO/QCISD/Traj/RP2.reorder.rot.coord/'


def totline(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    f.close()
    return i + 1


# def getInput():


# def main():
#     # 1. I/O


#     # Main program
# main()
