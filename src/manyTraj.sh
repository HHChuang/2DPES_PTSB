#!/bin/bash 
#########################################################################
# Program:																#
#	Classify trajectories into two groups form dir. /totTraj to two		#
#	sub-directory /RP1.reorder and /RP2.reorder.						#
#																		#
# Input; command-line arguments:										#
#	1. defR.dat															#
#	2. defProd1.dat														#
#	3. defProd2.dat														#
#																		#
# Output:																#
#	1. /RP1.reorder														#
#	2. /RP2.reorder														#
#																		#
# History:																#
# 2018/12/14, Grace, the previous one just dispear! fuck				#
# 2019/07/02, Grace, clear the set-up of Path variabls					#
#########################################################################

# Set-up path
	# $1 = defR.dat
	# $2 = defP1.dat 
	# $3 = defP2.dat

	countProg='/work/Grace/ProgDyn/Program/CountProd' # path of program
	inputPath=$(pwd)

	defR_path=$inputPath/$1
	defP1_path=$inputPath/$2
	defP2_path=$inputPath/$3
	colTraj=$inputPath/totTraj 
	RP1dir=$inputPath/RP1.reorder
	RP2dir=$inputPath/RP2.reorder

# change directory
	cd $colTraj
	ls | grep ^Traj'[0-9]' > trajlist.dat
	num=$(wc -l trajlist.dat | awk '{print $1}')
	cp $defR_path $defP1_path $defP2_path .

	rm -f record.dat reorder.Traj*
	for name in `cat trajlist.dat ` 
	do
		record=$($countProg $name $1 $2 $3)
		echo $name $record >> record.dat
		echo $name $record
	done
	rm -f trajlist.dat $1 $2 $3 

# count and group trajectories
	NN=$(grep -c 'none none' record.dat)
	RR=$(grep -c 'R R' record.dat)
	P1P1=$(grep -c 'P1 P1' record.dat)
	P2P2=$(grep -c 'P2 P2' record.dat)
	NR=$(grep -c 'none R' record.dat)
	NP1=$(grep -c 'none P1' record.dat)
	NP2=$(grep -c 'none P2' record.dat)
	P1P2=$(grep -c 'P1 P2' record.dat)
	RP1=$(grep -c 'R P1' record.dat)
	RP2=$(grep -c 'R P2' record.dat)
	tot=$(($NN+$RR+$P1P1+$P2P2+$NR+$NP1+$NP2+$P1P2+$RP1+$RP2))

# make RP1 RP2 list and then copy selected trajectories
	rm -rf $RP1dir $RP2dir
	mkdir $RP1dir $RP2dir
	grep 'R P1' record.dat | awk '{print $1}'  > RP1list.dat
	grep 'R P2' record.dat | awk '{print $1}' > RP2list.dat
	rm -f record.dat

	for name in `cat RP1list.dat`
	do 
		mv reorder.$name $RP1dir
	done
	for name in `cat RP2list.dat`
	do 
		mv reorder.$name $RP2dir
	done
	rm -f reorder.Traj*
	mv RP1list.dat RP2list.dat $inputPath

cat << EOF | tee $inputPath/analyseTraj.dat
----------------------------------------------------------------

    Analyse the initial/final point of all the trajectories 
    
output:
    files
        1. $home/RP1list.dat
        2. $home/RP2list.dat
    directories
        1. $RP1dir
        2. $RP2dir

std-out: 
total traj.: $tot

OtherOther: $NN
RR: $RR
P1P1: $P1P1
P2P2: $P2P2
OtherR: $NR
OtherP1: $NP1
OtherP2: $NP2
P1P2: $P1P2

RP1: $RP1
RP2: $RP2
----------------------------------------------------------------
EOF

 
