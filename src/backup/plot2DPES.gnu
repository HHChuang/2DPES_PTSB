#!/bin/bash

# contour
unset clabel
set contour base
set dgrid3d 20,20
set cntrparam levels incr -114.42,0.005,-114.32

unset key
unset xtics
unset ytics
set ztics font 'Helvetica,15' -114.44,0.04,-114.26

splot '39_E.dat' u 1:2:3 w lines