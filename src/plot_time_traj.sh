#!/bin/bash
# $1 = indix of NCHn system, e.g. n=1 -> NCH1

dir1='RP1.coord'
dir2='RP2.coord'
home=`pwd`

function main {
    cd $home
    rm -rf mix
    mkdir mix

    makelist $dir1 mix
    makelist $dir2 mix 
    
    cd $home/mix 
    RP1list=( $(cat RP1.coord.list) )
    RP2list=( $(cat RP2.coord.list) )
    xlim=(3.1 6.5 2.8 3.2 6.1)
    index=$(( $1 - 1 ))
    plot_time_traj ${xlim[$index]}

    # 
    mv RP1.eps NCH$1.RP1.eps 
    mv RP2.eps NCH$1.RP2.eps 
    mv RP1RP2.eps NCH$1.RP1RP2.eps 
    mv $home/mix/*.eps $home 
}

function makelist {
    # $1 = name of dir.
    # $2 = final dir.
    # output: $1.list 

    cd $home/$1
    ls | grep coord > list 
    for name in `cat list`
    do 
        nl $name | awk '{print $2,$1}' > new.$name
    done 
    awk '{print "new."$1}' list > $1.list 
    rm -f list

    mv new.* $home/$2 
    mv $1.list $home/$2 
}



function plot_time_traj {
    # $1 = max value of x
gnuplot << EOF
set term postscript eps enhanced color font 'Times-Roman,20' size 4,4
set output "RP1.eps"
unset key

RP1list="${RP1list[*]}"
RP2list="${RP2list[*]}"

set size square
set xrange [0:$1]
set yrange [20:200]
set ylabel 'fs'
set xlabel 'VRI region ({/Symbol @{\140\140\140}\326}amu Bohr)'
plot for [file in RP1list] file u 1:2 w line lw 3 lc rgb 'red'

set output "RP2.eps"
plot for [file in RP2list] file u 1:2 w line lw 3 lc rgb 'blue'

set output "RP1RP2.eps"
plot for [file in RP1list] file u 1:2 w line lw 3 lc rgb 'red', \
     for [file in RP2list] file u 1:2 w line lw 3 lc rgb 'blue'
EOF
}

main $1