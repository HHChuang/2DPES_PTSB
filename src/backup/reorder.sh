#!/bin/bash
#############################################################
# input:                                                    #
#   $1 = filename                                           #
#   filename = *.xyz                                        #
#       format:                                             #
#               NAtoms                                      #
#               comment                                     #
#               atomname x y z                              #
#                                                           #
# output:                                                   #
#   reorder.$(comment) (one file or many files)             #
#                                                           #
# History:                                                  #
# 2018/10/04, Grace                                         #
#############################################################

function get_index(){
    # $1 = NAtoms
    # output: list.dat
    rm -f list.dat
    for ((i=1;i<=$1;i++))
    do
        read -p "atom number $i is originally :" index
        num=$(($index+2))
        echo $num >> list.dat
    done
}

function do_reorder(){
    # $1 = NAtoms
    # $2 = list.dat
    # $3 = filename
    # output: reorder.$3
    head -n 2 $3 > reorder.$3
    for ((i=i;i<=$1;i++))
    do  
        num=$(sed -n "$i,$i p" $2)
        sed -n "$num,$num p" $3 >> reorder.$3
    done
}

# main program
NAtoms=$(head -n 1 $1)
totalline=$(wc -l $1 | awk '{print $1}')
round=$(( $totalline/($NAtoms+2) ))
get_index $NAtoms #output: list.dat
for ((i=1;i<=$round;i++))
do
    head -n $(( ($NAtoms+2)*$i )) $1 | tail -n $(($NAtoms+2)) > tmp
    filename=$(sed -n '2,2 p' tmp)
    mv tmp $filename.xyz
    do_reorder $NAtoms list.dat $filename.xyz
    rm -f $filename.xyz
done
rm -f list.dat