#!/bin/bash
#########################################################################
# Program:                                                              #
#   Rotate and shift the structures of trajectories in order to be the  #
#   same orientation and same original of center of this 2D-PES.        #
#                                                                       #
# input:                                                                #
#   1. $(Grid)_Struc.xyz; structures of this surface                    #
#   2. trajectories in /RP1.reorder and /RP2.reorder                    #
#                                                                       #
# output:                                                               #
#   /RP1.reorder.rot and /RP2.reorder.rot                               #
#                                                                       #
# History:                                                              #
# 2018/12/27, Grace, PATH problem                                       #
# 2019/07/23, Grace, make IRC trajectories                              #
#########################################################################

# main program
function main(){
    # Set up PATHs
    rotProg='/Users/Grace/Google_Drive/Code/GitHub/2DPES_PTSB/src/rot.py'
    # first input, $1 #######################################################
    tssPESname='tss.59_PES.xyz' 
    #########################################################################
    tssPESpath='/Users/Grace/Google_Drive/Code/GitHub/2DPES_PTSB/run/NCH1'
    # second input, $2
    iniRP1=$tssPESpath/Traj/RP1.reorder/
    iniRP2=$tssPESpath/Traj/RP2.reorder/
    iniRP3=$tssPESpath/Traj/IRC.reorder/
    #########################################################################
    finRP1=$tssPESpath/Traj/RP1.reorder.rot/
    finRP2=$tssPESpath/Traj/RP2.reorder.rot/
    finRP3=$tssPESpath/Traj/IRC.reorder.rot/

    # echo 'RP1'
    # rot $iniRP1 $finRP1
    # echo 'RP2'
    # rot $iniRP2 $finRP2
    echo 'IRC'
    ircTraj $iniRP3
    rot $iniRP3 $finRP3
}

# extract TS structure from one trajectory
function getTSStraj(){
    # $1 = name of trajectory
    # $2 = name of ts = ts.$1
    NAtoms=$(head -n 1 $1 )
    grep -B 1 -A $NAtoms 'runpoint 1' $1 > $2
}

function rot(){
    # $1 = initial dir. #absoulate path 
    # $2 = final dir.
    rm -rf $2
    mkdir $2 
    cp -r $1 $2
    cp $tssPESpath/$tssPESname $2
    cd $2
    ls | grep reorder.Traj > list
    
    for name in `cat list`
    do 
        getTSStraj $name tss.$name # output: ts.$name
        $rotProg tss.$name $tssPESname $name # output: rot.$name 
    done
    rm -f reorder.Traj* list tss*
}

function ircTraj(){
    cd $1
    pwd
    ls | grep reorder > list
    for name in `cat list`
    do 
        line=$(grep -n ^0.0000 $name | cut -d ':' -f 1)
        sed -e "$line s/0.0000/0.0000 runpoint 1/" $name > tmp 
        mv tmp $name
    done
    rm -f list 
}

main