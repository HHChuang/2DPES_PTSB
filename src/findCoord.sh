#!/bin/bash
#########################################################################
# Program:                                                      `       #
#   Searching the corresponding coordinates of trajectories.            #
#   This script integrates the following programs:                      #
#   1. CountProd.f90                                                    #
#   2. MapTraj.f90                                                      #
#   3. rot.py                                                           #
#                                                                       #
# Input command-line arguments:                                         #
#   $1 = defR.dat                                                       #
#   $2 = defP1.dat                                                      #
#   $3 = defP2.dat                                                      #
#   $4 = *_E.dat; coordinate and energy of numerical PES                #
#   $5 = *_Struc.xyz; all the structures of numberical PES              #
#                                                                       #
# Output files:                                                         #
#   1. AnalysedTraj.dat                                                 #
#   2. SummaryTraj.dat                                                  #
#                                                                       #
# Input directories:                                                    #
#   1. /rawTraj; rawadata from ProgDyn program                          #
#                                                                       #
# Output directories:                                                   #
#   1. /totTraj; collecting trajectories into one directory             #
#   2. /RP1.coord: coordinate of trajectories                           #
#   3. /RP2.coord: coordinate of trajectories                           #
#                                                                       #
# History:                                                              #
# 2019/06/04, Grace, revised                                            #
# 2019/07/23, Grace, change coordinate from integer to floating points  #
# 2019/07/25, Grace, integrate lots functions and assign this script as #
#   the working horse.                                                  #
#                                                                       #
#########################################################################

function main(){
    # Input arguments:
    #   $1 = defR.dat; input argument for manyTraj()
    #   $2 = defP1.dat; input argument for manyTraj()
    #   $3 = defP2.dat; input argument for manyTraj()
    #   $4 = *_E.dat; input argument for adjustTraj()
    #   $5 = *_Struc.xyz; input argument for adjustTraj()
    # Output:
    #   1. $home/AnalysedTraj.dat 
    #   2. $home/SummaryTraj.dat 
    #   3. /RP1.coord 
    #   4. /RP2.coord

    # Step 1. Set up PATHs and important variables
        home=$(pwd)
        # PATHs of programs
            # 1. countProg: classify trajectories by extract the last structure,
            #   and then assign them by user-defined criteria. e.g. bond length 
            #   or bond angle of reactant or products. 
                countProg='/work/Grace/ProgDyn/Program/CountProd' 
            # 2. rotProg:
                rotProg='/work/Grace/ProgDyn/Program/rot.py' 
                # rotProg='/Users/Grace/Google_Drive/Code/GitHub/2DPES_PTSB/src'
            # 3. mapProg: 
                mapProg='/work/Grace/ProgDyn/Program/MapTraj'
                # mapProg='/Users/Grace/Google_Drive/Code/GitHub/2DPES_PTSB/src/MapTraj'

            rawTraj="$home/DoneTraj"
            totTraj="$home/totTraj"
            RP1dir="$home/RP1.coord"    
            RP2dir="$home/RP2.coord"
   
    # Step 2. Collect the rawadata trajectories (i.e. calculate from ProgDyn)
    #   to one directories and also rename them. Replace the atomic symbol to 
    #   integer and then standard output the total amount of trajectories. 
        #FIXME:
        ##############################
        # collectTraj $rawTraj $totTraj
        ############################## 
  
    # Step 3. Classify trajectories by calling program CountProd.f90, assign 
    #   them into several groups based on user-defined criteria. 
    #   The main two groups are RP1 and RP2.
        #FIXME:
        ##############################
        #  manyTraj $1 $2 $3 $totTraj
        ##############################
        # Output:
        #   1. AnalysedTraj.dat 
        #   2. SummaryTraj.dat 

    # Step 4. Calling rot.py to make rotational matrix based on the overlap of the
    #   first point of trajectory and the TSS1 of potential. Apply this matrix on 
    #   this trajectory. 
        #FIXME: 
        ############################################################
        echo ''
        echo 'Rotate trajectories, print RMSD befor/after.' 
        echo ''
        adjustTraj $totTraj AnalysedTraj.dat $4 $5 
        ############################################################
        # Output: rot.*
 
    # Step 5. Map trajectory; search their corresponding coordinate on this 
    #   numerical potential hypersurface.  Calling MapTraj.f90 to search coordinate.
        ############################################################
        echo ''
        echo 'Search the corresponding coordinates.'
        echo ''
        time mapTraj $totTraj AnalysedTraj.dat $4 $5 | tee stdMapTraj.dat 
        rm -rf $RP1dir $RP2dir
        mkdir $RP1dir 
        mkdir $RP2dir
        echo 'Output directories: RP1.coord and RP2.coord'
        moveTraj AnalysedTraj.dat $totTraj $RP1dir $RP2dir 
        ############################################################
        # Output: 
        #   1. stdMapTraj.dat 
        #   2. coord.* in /$RP1dir and /$RP2dir
}

function collectTraj(){
    # $1 = $rawTraj ; PATH of the initial directory
    # $2 = $totTraj ; PATH of the final directory 
    iniPath=$1
    finPath=$2 

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
}

function extractFinalPts(){
    # Input: 
    #   $1 = trajectory 
    # Output: 
    #   1. Pts1.xyz
    #   2. Pts2.xyz 

    NAtoms=$(head -n 1 $1)
    tail -n $(( $NAtoms + 2 )) $1 > Pts1.xyz

    tmpL=$(grep -n 'runpoint 1' $1 | tail -n 1 | cut -d ':' -f 1)
    finL=$(( $tmpL - 2 ))
    iniL=$(( $finL - $NAtoms -1  ))
    sed -n "$iniL,$finL p" $1 > Pts2.xyz
}

function manyTraj(){
    # Input arguments:
    #   $1 = defR.dat
    #   $2 = defP1.dat
    #   $3 = defP2.dat
    #   $4 = PATH of collecting trajectories 
    # Output:
    #   1. AnalysedTraj.dat 
    #   2. SummaryTraj.dat

    defR=$1
    defP1=$2
    defP2=$3
    colTraj=$4 

    # 1. Copy user-defined criteria to the directory which has 
    # all the collecting trajectories
        cd $home 
        cp $defR $defP1 $defP2 $colTraj
        cd $colTraj
        # ls | grep ^Traj'[0-9]' > trajlist.dat
        ls | grep ^Traj > trajlist.dat
        num=$(wc -l trajlist.dat | awk '{print $1}')

        rm -f AnalysedTraj.dat 
        for name in `cat trajlist.dat ` 
        do
            extractFinalPts $name # output: Pts1.xyz Pts2.xyz
            Pts1=$($countProg Pts1.xyz $defR $defP1 $defP2)
            Pts2=$($countProg Pts2.xyz $defR $defP1 $defP2)
            echo $name $Pts1 $Pts2 >> AnalysedTraj.dat 
            # echo $name $Pts1 $Pts2 # FIXME: uncommand for debug 
        done
        rm -f trajlist.dat $defR $defP1 $defP2 Pts1.xyz Pts2.xyz 

    # 2. Count and group trajectories
        NN=$(grep -c 'none none' AnalysedTraj.dat )
        RR=$(grep -c 'R R' AnalysedTraj.dat )
        P1P1=$(grep -c 'P1 P1' AnalysedTraj.dat )
        P2P2=$(grep -c 'P2 P2' AnalysedTraj.dat )
        NR=$(grep -c 'none R' AnalysedTraj.dat )
        NP1=$(grep -c 'none P1' AnalysedTraj.dat )
        NP2=$(grep -c 'none P2' AnalysedTraj.dat )
        P1P2=$(grep -c 'P1 P2' AnalysedTraj.dat )
        RP1=$(grep -c 'R P1' AnalysedTraj.dat )
        RP2=$(grep -c 'R P2' AnalysedTraj.dat )
        tot=$(( $NN + $RR + $P1P1 + $P2P2 + $NR + $NP1 + $NP2 + $P1P2 + $RP1 + $RP2 ))

cat << EOF | tee $home/SummaryTraj.dat 
----------------------------------------------------------------

    Analyse the initial/final point of all the trajectories 

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
    cp AnalysedTraj.dat $home 
}

function adjustTraj(){
    # Use rot.py to rotate molecules 
    # Input: 
    #   $1 = $totTraj 
    #   $2 = AnalysedTraj.dat
    #   $3 = *_E.dat 
    #   $4 = *_Struc.xyz 
    # output: 
    #   rot.* 
    
    cd $1 
    cp $home/$3 .
    cp $home/$4 .
    awk '{print $1}' $2 > list.dat 
    testfile=$(head -n 1 list.dat)
    NAtoms=$(head -n 1 $testfile)
    jobL=$(( $NAtoms + 2 ))

    for name in `cat list.dat`
    do 
        # 1. extract TSS from one trajectory
        grep -B 1 -A $NAtoms 'runpoint 1' $name \
            | tail -n $(($NAtoms +2 )) > tss.$name
        # 2. extract TSS from the numerical trajectory 
            order=$(grep -n '0.0000 0.0000' $3 | cut -d ':' -f 1)
             # FIXME: only for test system : H2O
            # order=$(grep -n '0.97587863 108.25046046' $3 | cut -d ':' -f 1)
            iniL=$(( 1 + $jobL * ($order -1) ))
            finL=$(( $jobL * $order ))
            sed -n "$iniL,$finL p" $4 > tss.PES.xyz
            $rotProg tss.$name tss.PES.xyz $name # output: rot.$name 
        done 
    rm -f tss.* list.dat tss.PES.xyz
}

function mapTraj(){
    # Use MapTraj.f90 to search coordinate
    # Input: 
    #   $1 = $totTraj: PATH of input trajectory; rot.*
    #   $2 = AnalysedTraj.dat 
    #   $3 = $*_E.dat: 
    #   $4 = $*_Struc.xyz: 
    # Output: 
    #   

    EngList=$3
    StrucList=$4

    cd $home 
    cp $2 $3 $4 $1
    cd $1 

    # 1. Extract correct name list for trajectories 
    awk '{print "rot."$1, $2, $3 }' $2 > List.dat 

    # 2. Search the corresponding coordinate 
    njobs=$(wc -l List.dat | awk '{print $1}')
    for (( i=1; i<=$njobs; i++))
    do 
        name=$(sed -n "$i,$i p" List.dat | awk '{print $1}')
        ReverseTraj=$(sed -n "$i,$i p" List.dat | awk '{print $2}')
        ForwardTraj=$(sed -n "$i,$i p" List.dat | awk '{print $3}')
        $mapProg $name $EngList $StrucList $ReverseTraj $ForwardTraj # output: coord.$name
    done
    # rm -f List.dat 
}

function moveTraj(){
    # Input: 
    #   $1 = AnalysedTraj.dat 
    #   $2 = $totTraj 
    #   $3 = $RP1dir 
    #   $4 = $RP2dir 

    cd $2 

    # Generate the filename list of RP1 and RP2 trajectories
    grep R $1 | grep P1 | awk '{print "coord.rot."$1 }' > RP1list.dat 
    grep R $1 | grep P2 | awk '{print "coord.rot."$1 }' > RP2list.dat 

    for name in `cat RP1list.dat`
    do 
        mv $name $3 
    done

    for name in `cat RP2list.dat`
    do 
        mv $name $4 
    done

    rm -f RP1list.dat RP2list.dat 
}

main $1 $2 $3 $4 $5