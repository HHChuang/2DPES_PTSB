#!/bin/bash

upfilehome=/Users/Grace/Google_Drive/Code/GitHub/ModifyIRC/run/Nitrene
# for upfile in NCH2 NCH3 NCH4 NCH5 NCH6 NCH7
for upfile in NCH4
do
    cd $upfilehome/$upfile
    home=$(pwd)
    # cp -r /Users/Grace/Google_Drive/Code/GitHub/ModifyIRC/aux/asymmetricPES/$upfile/traj/n*_modifiedAtomName .
    filelist=$(ls | grep modifiedAtomName)
    for file in `echo $filelist`
    do cd $home/$file
        for traj in `cat list`
        do
            /Users/Grace/Google_Drive/Code/GitHub/ModifyIRC/src/mapTrajectory $traj $home/35_E_Struc.dat
            # /Users/Grace/Google_Drive/Code/GitHub/ModifyIRC/src/mapTrajectory $traj $home/15_E_Struc.dat
        done
    done
done
