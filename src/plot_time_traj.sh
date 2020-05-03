#!/bin/bash
# $1 = indix of NCHn system, e.g. n=1 -> NCH1

dir1='RP1.coord'
dir2='RP2.coord'
home=`pwd`

# maximum limit of x-range
xmaxList=(3.1 6.5 2.8 3.2 6.1) 
# green area
objiniList=(0.5 2.0 0.5 0.5 0.5)
objfinList=(1.0 5.5 1.0 1.5 4.0)

function main {
    cd $home
    rm -rf mix
    mkdir mix

    echo $dir1
    makelist $dir1 mix $1
    echo $dir2
    makelist $dir2 mix $1 
    
     
    cd $home/mix 
    RP1list=( $(cat RP1.coord.list) )
    RP2list=( $(cat RP2.coord.list) )
    plot_time_traj $1 

    # 
    mv RP1.eps NCH$1.RP1.eps 
    mv RP2.eps NCH$1.RP2.eps 
    mv RP1RP2.eps NCH$1.RP1RP2.eps 
    mv $home/mix/*.eps $home 
}

function makelist {
    # $1 = name of dir.
    # $2 = final dir.
    # $3 = indes of selected system
    # output: $1.list 

    cd $home/$1
    ls | grep coord > list 
    nTraj=$(wc -l list | awk '{print $1}')
    sumTraj=0
    for name in `cat list`
    do  
        TSS1toTSS2 $name # output: TSS1toTSS2$name
        nl TSS1toTSS2$name | awk '{print $2,$1}' > new.$name
        #
        countGreenT $3 new.$name # output: $singleTraj
        # echo 'new'.$name $singleTraj
        sumTraj=$(( $sumTraj + $singleTraj))
    done 
    # TODO:
    awk '{print "new."$1}' list > $1.list 
    avgTraj=$(echo "$sumTraj/$nTraj" | bc -l )
    echo $nTraj $avgTraj 

    rm -f list TSS1toTSS2*

    mv new.* $home/$2 
    mv $1.list $home/$2 
}

function countGreenT {
    # $1 = index of system 
    # $2 = name of traj
    # output: 
    # modify the value of $singleTraj

    index=$(( $1 - 1 ))
    objini=${objiniList[$index]}
    objfin=${objfinList[$index]}

    awk '{print $1}' $2 > testTraj
    i=0
    singleTraj=0
    while read LINE 
    do 
        i=$(( $i + 1 ))
        testXcoord=$(sed -n "$i,$i p" testTraj)
        if (( $(echo "$testXcoord > $objini" | bc -l) \
            && $(echo "$testXcoord < $objfin" | bc -l) ))
        then
            singleTraj=$(( $singleTraj + 1 )) 
        fi 
    done < testTraj
    rm -f testTraj
}

function TSS1toTSS2 {
    # $1 = name of the file

    # 1. Define the position of TSS1
      negSign=$(awk '{print $1}' $1 | grep -c '-')
      numL=$(($negSign + 1))

    # nl $1 | awk '{print $1,$2}' | grep '0.00000' \
    #     | grep -v '-' | awk '{print $1}' > x.dat 
    # nl $1 | awk '{print $1,$3}' | grep '0.00000' \
    #     | grep -v '-' | awk '{print $1}' > y.dat 
    # while read -r line
    # do 

    # done

    totL=$(wc -l $1 | awk '{print $1}')
    sed -n "$numL,$totL p" $1 > TSS1toTSS2$1 
}

function plot_time_traj {
    # $1 = index depends on selected system

    index=$(( $1 - 1 ))

    xmax=${xmaxList[$index]}
    objini=${objiniList[$index]}
    objfin=${objfinList[$index]}

gnuplot << EOF
# postscript does't support transparent
set term postscript eps enhanced color font 'Times-Roman,20' size 4,4
# set term png tranparent enhanced color font 'Times-Roman,20' size 4,4
set output "RP1.eps"
unset key

RP1list="${RP1list[*]}"
RP2list="${RP2list[*]}"

set size square
set xrange [0:$xmax]
# set yrange [20:200]
set yrange [0:100]
set ylabel 'fs'
set xlabel 'TSS1 {/Symbol \256} TSS2 ({/Symbol @{\140\140\140}\326}amu Bohr)'
set style rect fillcolor rgb "seagreen" fill transparent solid 0.5 noborder  
set obj rect from $objini, graph 0 to $objfin, graph 1 #front

plot for [file in RP1list] file u 1:2 w line lw 3 lc rgb 'red'

set output "RP2.eps"
plot for [file in RP2list] file u 1:2 w line lw 3 lc rgb 'blue'

set output "RP1RP2.eps"
plot for [file in RP1list] file u 1:2 w line lw 3 lc rgb 'red', \
     for [file in RP2list] file u 1:2 w line lw 3 lc rgb 'blue'
EOF
}

main $1