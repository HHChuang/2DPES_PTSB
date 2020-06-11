#!/bin/bash
# $1 = indix of NCHn system, e.g. n=1 -> NCH1

dir1='RP1.coord'
dir2='RP2.coord'
home=`pwd`
sysIndex=$1

# maximum limit of x-range
xmaxList=(3.1 6.5 2.8 3.2 6.1) 
# green area
objiniList=(0.5 2.0 0.5 0.5 0.5)
objfinList=(1.0 5.5 1.0 1.5 4.0)

function main {
    cd $home
    rm -rf mix
    mkdir mix

    makelist $dir1 mix # output: $dir1.list
    makelist $dir2 mix 
        
    # plot traj.
    cd $home/mix 
    RP1list=( $(cat RP1.coord.list) )
    RP2list=( $(cat RP2.coord.list) )
    plot_time_traj $sysIndex

    # calculate avgTime and SD
    # echo $dir1
    # calcStatistic $dir1.list 
    # echo $dir2
    # calcStatistic $dir2.list

    # collect figures
    mv RP1.eps NCH$sysIndex.RP1.eps 
    mv RP2.eps NCH$sysIndex.RP2.eps 
    mv RP1RP2.eps NCH$sysIndex.RP1RP2.eps 
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
        TSS1toTSS2 $name # output: TSS1toTSS2$name
        nl TSS1toTSS2$name | awk '{print $2,$1}' > new.$name
    done 
    awk '{print "new."$1}' list > $1.list 

    rm -f list TSS1toTSS2*
    mv new.* $home/$2 
    mv $1.list $home/$2 
}

function countGreenT {
    # $1 = index of system 
    # $2 = name of traj
    # output: 
    # modify the value of $singleTraj

    index=$(( $sysIndex - 1 ))
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
    objmid=$(echo "($objini + $objfin)/2.0" | bc -l )

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
set ylabel 'Time (fs)'
set xlabel 'TSS1 {/Symbol \256} TSS2 ({/Symbol @{\140\140\140}\326}amu Bohr)'
set style rect fillcolor rgb "seagreen" fill transparent solid 0.5 noborder  
set obj rect from $objini, graph 0 to $objfin, graph 1 #front
set x2tics ("TSS1" 0, "VRI" $objmid, "TSS2" $xmax)

plot for [file in RP1list] file u 1:2 w line lw 3 lc rgb 'red'

set output "RP2.eps"
plot for [file in RP2list] file u 1:2 w line lw 3 lc rgb 'blue'

set output "RP1RP2.eps"
plot for [file in RP1list] file u 1:2 w line lw 3 lc rgb 'red', \
     for [file in RP2list] file u 1:2 w line lw 3 lc rgb 'blue'
EOF
}

function calcStatistic {
    # $1 = name of list 
    # output (std-out): 
    # $avgTime; average time 
    # $SD; standard deviation

    nTraj=$(wc -l $1 | awk '{print $1}')
    sumTraj=0
    timeTraj=()
    timeind=0
    for name in `cat $1`
    do 
        countGreenT $sysIndex $name # output: $singleTraj
        timeTraj[$timeind]=$singleTraj
        timeind=$(( $timeind + 1 ))
        # echo $name $singleTraj
        sumTraj=$(( $sumTraj + $singleTraj))
    done 
    avgTraj=$(echo "$sumTraj/$nTraj" | bc -l )
    # echo ${timeTraj[@]}

    SD=0
    for ((i=0;i<$nTraj;i++))
    do 
        SD=$(echo "( ${timeTraj[$i]}*1.0 - $avgTraj )^2 + $SD" | bc -l)
    done 
    SD=$( echo "sqrt( $SD / ($nTraj-1) )"|bc -l )

    echo '$nTraj $avgTime $SD'
    printf "%d\t %.2f\t %.2f \n" $nTraj $avgTraj $SD
}
main 