#!/bin/bash
#################################################################
# Program:														#
#	collect trajectories from dir. /DoneTraj to /totTraj		#
#																#
# History:														#
# 	2018/11/29, Grace											#
#	2019/07/02, Grace, set up more clear path variables			#
#################################################################

# Set-up path
	home=$(pwd)
	inidir='DoneTraj'
	findir='totTraj'
	iniPath=$home/$inidir
	finPath=$home/$findir

# Print out purpose
cat << EOF
Current directory is $home
Initial directory is $iniPath
Final directory is $finPath

Collect trajectory and then,
1. rename them
2. change atomic number from character to integer

EOF

# Initialize directory
	rm -rf $finPath
	mkdir $finPath

# Main program
	cd $iniPath
	ls | grep ^n > fileList.tmp
	num=0
	for filename in `cat fileList.tmp`
	do
		cd $iniPath/$filename
		ls | grep traj| grep '[1-9]' > trajList.tmp
		for trajname in `cat trajList.tmp`
		do 
			num=$(($num+1))
			cp $trajname Traj$num
			# change atomic number 
			sed -i 's/C/6/g' Traj$num
			sed -i 's/H/1/g' Traj$num
			sed -i 's/N/7/g' Traj$num
			sed -i 's/O/8/g' Traj$num	
			sed -i 's/F/9/g' Traj$num
			mv Traj$num $finPath
		done
		rm -f trajList.tmp
	done
	echo "Total amount of trajectories: $num"

	cd $iniPath
	rm -f fileList.tmp



