#!/bin/bash
#############################################################
# 2018/12/14, Grace											#
# 2019/03/13, Grace, change the structure by using functions#
# prequiste script											#
# 	1. getJob in /bin 										#
#############################################################

home='/work/Grace/Tantillo_Dynamic/H3CO/QCISD/ProgRaw'
doneTraj='/work/Grace/Tantillo_Dynamic/H3CO/QCISD/DoneTraj'
tempDir='/work/Grace/Tantillo_Dynamic/Program/templet' #do not change this line
inifile=1 # first number of file -> n$inifile
finalfile=40 # final number of file -> n$finalfile

function purpose(){
	echo '-------------------------------------------------------------'
	echo "Calculate quasi-classical trjaectories by Signlton's program."
	echo '-------------------------------------------------------------'
	echo ''
	echo '1. First time to execute this script: '
	echo '	check the path is correct or not!!'
	echo ''
	echo '2. Following check the status of trajectories: '
	echo "Move all sucessful directories to $doneTraj "
	echo "and resubmit fail dir."
	echo ''
}

function submit1st(){
	# $1 = initial file number, e.g. 1 -> n1
	# $2 = final file numner, eg. 10 -> n10
	cd $home
	for ((i=$1;i<=$2;i++))
	do
		cd $home
		cp -r $tmpDir $home/n$i
		cd $home/n$i
		qsub job
	done
}

function countTraj(){
	# output: $home/dirtraj.dat
	cd $home 
	ls | grep ^n > dirlist.dat
	tottraj=0
	rm -f dirtraj.dat
	for dir in `cat dirlist.dat`
	do
		cd $home/$dir 
		ntraj=$(ls | grep -c traj'[0-9]')
		tottraj=$(( $ntraj + $tottraj ))
		echo $dir $ntraj >> $home/dirtraj.dat
	done
	echo "Current total trajectories in this dir.: $tottraj"
	rm -f dirlist.dat 
}

function move2Done(){
	# output: 
	#	1. faillist.dat
	#	2. worklist.dat
	cd $home
	getJob > runlist.dat
	# grep fail jobs
	rm -f faillist.dat worklist.dat
	awk '{print $2}' dirtraj.dat | grep -v -n [1-9] \
		| cut -d ':' -f 1 > tmp
	for line in `cat tmp`
	do 
		sed -n "$line,$line p" dirtraj.dat | awk '{print $1}' >> faillist.dat
	done
	# grep sucessful jobs
	awk '{print $2}' dirtraj.dat | grep -n [1-9] \
		| cut -d ':' -f 1 > tmp
	for line in `cat tmp`
	do 
		sed -n "$line,$line p" dirtraj.dat | awk '{print $1}' >> worklist.dat
	done
	rm -f tmp

	for dir in `cat worklist.dat`
	do
		run=$(grep -c "$home/$dir" runlist.dat)
		[ $run == 0 ] && mv $home/$dir $doneTraj
	done
}

function resubmitFail(){
	for dir in `cat faillist.dat`
	do
		run=$(grep -c "$home/$dir" runlist.dat)
		if [ $run == 0 ] 
		then
			rm -rf $home/$dir
			cp -r $tempDir $home/$dir
			cd $home/$dir
			qsub job
			cd $home
		fi
	done
}

# main program 
function main(){
	# Step 1. stdout the purpose 
	purpose
	# Step 2. Initialize and for the first time to submit all jobs
	cd $home 
	rm -f runlist.dat faillist.dat worklist.dat dirlist.dat 
	submit1st $inifile $finalfile
	# Step 3. Count the total amount of trajectories
	# countTraj # output: $home/dirtraj.dat
	# Step 4. Move all sucessful file to DoneTraj
	# move2Done # output: $home/faillist.dat and $home/worklist.dat
	# Step 5. resubmit fail lsit
	# resubmitFail
	# rm -f runlist.dat faillist.dat worklist.dat dirtraj.dat dirlist.dat 
}

main