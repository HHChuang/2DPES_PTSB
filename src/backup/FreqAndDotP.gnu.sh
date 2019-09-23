#!/bin/bash
#########################################################
# Program                                               #
#   Plot the imaginary frequency (from negative         #
#   eigenvalue) and the corresponding dot product with  #
#   the selective mode (TS2's imaginary frequency mode  #
#   or one of the negative eigenvalue mode from its own #
#   Hessian).                                           #
#                                                       #
# I/O part                                              #
#   input:                                              #
#       $1 = $PATH/E.dat
#       $1 = $PATH/dot.dat                              #
#   output:                                             #
#       FreqAndDotP.eps                                 #
#                                                       #
# History                                               #
# 2018/07/11, Grace                                     #
#########################################################


gnuplot << EOF
set terminal postscript eps enhanced color size 10cm,15cm
set output 'EandFreqAndDotP.eps'

set lmargin 15
set rmargin 3
set multiplot layout 3,1 
set xrange [0:10]

# Energy ################################################
set bmargin 0 
unset xtics
unset key
set origin 0,0.73
set size 1,0.28
set ylabel 'Energy (Hartree)' font ',15'
p "$1" u 1:2 w line lw 5 lc 'black'  

# Frequency #############################################
set tmargin 0
set bmargin 0
unset xtics
unset key
set yrange [-300:100]
set origin 0,0.45
set size 1,0.28
set ylabel 'Frequency (cm^{-1})' font ',15' offset -2,0
p "$2" u 1:2 w linespoints title '1st. mode' lc 'black' pt 6#, \
    # "$2" u 1:4 w linespoints title '2nd. mode' lc 'black' pt 13, \
    # "$2" u 1:6 w linespoints title '3rd. mode' lc 'black' pt 1

# Dot product ###########################################
set tmargin 0
set bmargin 5
set xtics
set key
set yrange [-0.006:0.001]
set origin 0,0.1
set size 1,0.35
set ylabel 'Dot Product' font ',15' offset 1,0
set xlabel 'Intrinsic reaction coordinate' font ',15'
p "$2" u 1:3 w linespoints title '1st. mode' lc 'black' pt 6#, \
    # "$2" u 1:5 w linespoints title '2nd. mode' lc 'black' pt 13, \
    # "$2" u 1:7 w linespoints title '3rd. mode' lc 'black' pt 1
EOF