#!/bin/bash
#####################################################################
# 2019/03/15, Grace, pick the first and the last structures         #
# 2019/03/26, Grace, add the second arguments in order to generate  #
#   the smooth energy curve.                                        #
# 2019/03/27, Grace, make sure the direction starts from TSS to     #
#   minimum for energy profile as well as structures.               #
# 2019/03/28, Grace, the orientation of TSS is different from other #
#   IRC structures, use a python code to fix this problem.          #
#                                                                   #
# pre-requisited program                                            #
#   1. rot.py                                                       #
#                                                                   #
# $1 = potential energy curve                                       #
#       format:                                                     #
#               coord. energy                                       #
# $2 = a serious of structures                                      #
#       format:                                                     #
#               NAtoms                                              #
#               comment line                                        #
#               aromtN atomCoord(x,y,z)                             #
# $3 = amount of selected structures (integer)                      #
#####################################################################

NAtoms=$(head -n 1 $2)
totStruc=$(wc -l $1 | awk '{print $1}' )

# the order of $1 should start from transition-state structure
function checkOrder(){
    # $1 = potential energy curve                  
    if [ `awk '{print $1}' $1 | grep -c \-` -gt 1 ]
    then 
        sort -nr -k 1 $1 > tmp  
        mv -f tmp $1 
    fi
}

# 2019/03/28, fix the TSS orientation problem of the $2 
function orienTSS(){
    # $1 = a serious of structures    
    fileLine=$(( $NAtoms+2 ))
    sed -n "1,$fileLine p" $1 > IRC_TS.xyz
    sed -n "$(( $fileLine+1 )), $(( $fileLine*2 )) p" $1 > IRC_2.xyz
    # rotate and shift TSS; output: rot.IRC_TS.xyz
    python3 /Users/Grace/Google_Drive/Code/GitHub/2DPES_PTSB/src/rot.py IRC_TS.xyz IRC_2.xyz IRC_TS.xyz
    # replace the first struc. in IRC structure file
    totLine=$( wc -l $1 | awk '{print $1}')
    cat rot.IRC_TS.xyz > tmp
    sed -n "$(( $fileLine +1 )), $totLine p " $1 >> tmp 
    mv -f tmp $1
    rm -f IRC_*.xyz rot.IRC_TS.xyz
}

function purpose(){
    # $1 = amount of selected structures
    echo '------------------------------------------'
    echo "  From $totStruc structures, select $1.   "
    echo '  Output: '
    echo "      1. $3_$1"
    echo "      2. $3_$2"
    echo '------------------------------------------'
}

function CalcIntE(){
    # input:
    #   $1 = PEC
    #   $2 = $totStruc
    #   $3 = amount of selected structures (integer)     
    # output:
    #   $intE
    maxE=$(sort -n -k 2 $1 | tail -n 1 | awk '{print $2}')
    minE=$(sort -n -k 2 $1 | head -n 1 | awk '{print $2}')
    intE=$(echo " scale=10; ($maxE - $minE)/($3*1.0) " | bc )
}

function pick(){
    # input:
    #   $1 = PEC
    #   $2 = a serious of structures   
    #   $3 = amount of selected structures (integer)   
    # output:
    #   $3_$1
    #   $3_$2
    # setup for reference file
    # extract the first and the last point in the original file
    E0=$(head -n 1 $1 | awk '{print $2}')
    echo $E0 > refE.dat 
    for ((i=1;i<=$(($3-2));i++))
    do  
        echo " scale=10; $E0 - $i*$intE " | bc >> refE.dat
    done
    tail -n 1 $1 | awk '{print $2}' >> refE.dat
    # setup for output files
    output1=$(echo $3_$1)
    output2=$(echo $3_$2)
    # extract the first point
    head -n 1 $1 > $output1
    head -n  $(($NAtoms+2))  $2 > $output2
    # extract the points in between
    iniIndex=1    
    for ((i=2;i<=$(($3-1));i++))
    do
        refE=$(sed -n "$i,$i p" refE.dat)
        for ((j=$iniIndex;j<=$totStruc;j++))
        do
            tmpE=$(sed -n "$j,$j p" $1 | awk '{print $2}')
            delE=$(echo "scale=10; $tmpE - $refE " | bc)
            count=$(echo $delE | grep -c \-)
            if [ "$count" = "1" ]
            then
                coord=$(sed -n "$j,$j p" $1 | awk '{print $1}')
                # print out the selected point 
                echo $coord $tmpE >> $output1
                # extract the corresponding structure
                iniL=$(( 1+($j-1)*($NAtoms+2) ))
                finL=$(( $j*($NAtoms+2) ))
                sed -n "$iniL,$finL p" $2 >> $output2
                iniIndex=j+1
                break
            fi
        done
    done
    # extract the last point
    tail -n 1 $1  >> $output1
    tail -n $(($NAtoms+2)) $2 >> $output2
    rm -f refE.dat
}

function main(){
    # purpose $1 $2 $3 
    # checkOrder $1 # output: $1 
    orienTSS $2 # output: $2
    # CalcIntE $1 $totStruc $3 # output: $intE
    # pick $1 $2 $3 $intE # also use $intE, output: $3_$1 and $3_$2 
}

main $1 $2 $3 