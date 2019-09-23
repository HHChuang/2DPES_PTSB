# !/bin/bash
function collectRMSD(){
    Prog='/Users/Grace/Google_Drive/Code/GitHub/Rotation/src/kabsch.py'
    ref='/Users/Grace/Google_Drive/Code/GitHub/2DPES_PTSB/run/H3CO/tsPES.xyz'
    for namelist in RP1list RP2list
    do
        if [ $namelist == RP1list ] 
        then
            dir='RP1.ts'
            rawdata='RP1RMSD'
        else
            dir='RP2.ts'
            rawdata='RP2RMSD'
        fi

        # collect RMSD for all trajectory
        for name in `cat $namelist`  
        do 
            $Prog $dir/ts.$name $ref >> `echo $rawdata`
        done

        # count number of trajectory

        # before rotation and translation
        filename=$(echo $namelist | sed 's/list//g')
        rm -f $filename.old $filename.new
        for num in `cat $rawdata| awk '{print $1}' | sort -n | uniq `
        do
            traj=$(grep -c $num `echo $rawdata` )
            echo $num $traj >> $filename.old
        done 

        # after rotation and translation
        filename=$(echo $namelist | sed 's/list//g')
        for num in `cat $rawdata| awk '{print $2}' | sort -n | uniq`
        do
            traj=$(grep -c $num `echo $rawdata` )
            echo $num $traj >> $filename.new
        done
    done

    echo '#RP1old Traj RP2old Traj RP1new Traj RP2new Traj' > RMSD.dat
    paste RP1.old RP2.old RP1.new RP2.new >> RMSD.dat
    rm -f RP1.old RP2.old RP1.new RP2.new RP1RMSD RP2RMSD
}

function plotHis(){
    # $1 = RMSD.dat
    gnuplot << EOF
set terminal postscript eps color enhanced
set style fill solid

set xlabel 'RMSD'
set ylabel 'Amount of trajectory'
set size ratio 1

set output 'oldRMSD.eps'
set title 'Statistic of RMSD before rotation'
p "$1" u 1:2 w boxes lc rgb 'black' title 'RP1', \
    "$1" u 3:4 w boxes lc rgb 'red' title 'RP2'

set output 'newRMSD.eps'
set title 'Statistic of RMSD after rotation'
p "$1" u 5:6 w boxes lc rgb 'black' title 'RP1', \
    "$1" u 7:8 w boxes lc rgb 'red' title 'RP2'
EOF
}

# 1. collect RMSD datat from trajectory
collectRMSD #output: RMSD.dat and rotM.Traj*
# 2. draw the histogram of RMSD vs. # Traj
plotHis RMSD.dat #output: oldRMSD.eps and newRMSD.eps