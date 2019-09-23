#!/bin/bash
#################################################################
# Program                                                       #
#   Analysis the reaction path from TS1 to Min1 (one of the     #
#   product in bifurcating reaction), and then using these      #
#   information to construct the reaction path from             #
#   vally-ridge-inflection (VRI) point to TS2.                  #
#                                                               #
# I/O part                                                      #
#     input:                                                    #
#       $1 = TS1_IRC.log                                        #
#     output:                                                   #
#       n$coord.com; a lot G09 input files                      #
#                                                               #
# Flowchart                                                     #
#         1. Define the range of x-axis                         #
#         2. Extract a serious of structures with corresponding #
#            coordinate as its file name.                       #
#         3. Calculate gradient and Hessian                     #
#         4. Generate Gaussian09 input files                    #
#                                                               #
# History                                                       #
#   2018/07/02, Grace                                           #
#################################################################

# Functions
# 1. Define the range of x-axis  
function defineXaxis(){
    echo '------------------------------------------'
    echo 'Step 1. Define x-axis'
    echo 'The x-axis is the forward IRC path of TS1'
    echo '------------------------------------------'
}
# 2. Extract a serious of structures with corresponding 
    #    coordinate as its file name.          
    #   input : 
    #           $1 = TS1_IRC.log
    #   output : 
    #           PES.dat
function getPES(){
    echo '------------------------------------------'
    echo 'Step 2. Extract structures along IRC'
    echo 'output : PES.dat'
    echo '------------------------------------------'

    nirc=$(grep -c 'NET REACTION COORDINATE UP TO THIS POINT' $1) 
    natoms=$(grep 'NAtoms' $1 | head -n 1 | awk '{print $2}')
    
    # paste TS info
    echo $natoms > PES.dat
    echo 0.00000 >> PES.dat 
    grep -A $((4+$natoms)) 'Input orientation' $1 \
    | head -n $((5+$natoms)) | tail -n $natoms \
    | awk '{print $2,$4,$5,$6}' >> PES.dat

    grep -n 'Point Number:' $1 | awk '{print $1}' \
        | sed 's/://g' > tmpStartLine.dat
    grep -n 'NET REACTION COORDINATE UP TO THIS POINT' $1 \
        | awk '{print $1}' | sed 's/://g' > tmpEndLine.dat
    grep 'NET REACTION COORDINATE UP TO THIS POINT' $1 \
        | awk '{print $9}' > tmpCoord.dat
    
    for ((a = 1; a <= $nirc ; a++))
    do 
        startLine=$(sed -n "$a,$a p" tmpStartLine.dat)
        endLine=$(sed -n "$a,$a p" tmpEndLine.dat)
        coord=$(sed -n "$a,$a p" tmpCoord.dat)
        # file structue of PES.dat 
            #   natoms    
            #   coord name
            #   coordinate ($atom x y z)
        echo $natoms >> PES.dat
        echo $coord >> PES.dat
        sed -n "$startLine,$endLine p" $1 \
            | grep -A $(($natoms+2)) 'Coordinates (Angstroms)'  \
            | tail -n $natoms | awk '{print $2,$4,$5,$6}' >> PES.dat
    done
    rm -f tmp*.dat 
}             

# 3. Calculate gradient and Hessian                             
# 4. Generate Gaussian09 input files                    
    # input :
    #       $1 = coord.dat
    #       $2 = struc.dat
    # output : 
    #       n$coord.com
function genInput(){
    coord=$(cat $1)
cat << EOF > n$coord.com
# $dft/$basis force freq=hpmodes

ModifyIRC

$charge $multi
EOF
cat $2 >> n$coord.com
echo '' >> n$coord.com
}


# Main program
    defineXaxis
    getPES $1 # output: PES.dat
    # assign variables
        natoms=$(head -n 1 PES.dat)
        nline=$(wc -l PES.dat | awk '{print $1}')
        nfile=$(( $nline/($natoms+2) ))
        basis=$(grep 'Standard basis' $1 | head -n 1 | awk '{print $3}')
        dft=$(grep 'SCF Done' $1 | head -n 1| awk '{print $3}' \
            | sed 's/E(//g' | sed 's/)//g')
    charge=$(grep 'Charge' $1 | head -n 1 | awk '{print $3}')
    multi=$(grep 'Multiplicity' $1 | head -n 1 | awk '{print $6}')

    echo '------------------------------------------'
    echo 'Main program: generate G09 input files'
    echo '------------------------------------------'
    for ((i=1;i<=$nfile;i++))
        do
        head -n $(($i*($natoms+2))) PES.dat | tail -n $(($natoms+2)) \
        | head -n 2 | tail -n 1 > coord.dat
        head -n $(( $i*($natoms+2) )) PES.dat | tail -n $natoms \
            > struc.dat
        genInput coord.dat struc.dat 
    done

    rm -f coord.dat struc.dat 