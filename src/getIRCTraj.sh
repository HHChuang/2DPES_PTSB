# !/bin/bash
#################################################################
# Program:                                                      #
#                                                               #
# pre-requist script:                                           #
#   1. getIRCstruc                                              #
#   2. rev1Dstruc                                               #
#   3. rot.py                                                   #
#                                                               #
# Input:                                                        #
#   $1 = *_F.log; IRC output file                               #
#   $2 = *_R.log; IRC output file                               #
#                                                               #
# History:                                                      #
#   2019/07/05, Grace                                           #
#################################################################

function main(){
    # $1 = TS*_F.log 
    # $2 = TS*_R.log 

    name_F=$(echo $1 | sed 's/.log//g')
    name_R=$(echo $2 | sed 's/.log//g')
    name=$(echo $name_F | sed 's/_F//g')

    # 1. Extract IRC structure for forward and reverse direction
    changeName $name_F # output: $name.xyz 
    changeName $name_R
    rev1Dstruc $name_R.xyz rev.$name_R.xyz 

    # 2. Extract TSS
    NAtoms=$(grep NAtoms $name_F.log | tail -n 1 \
        | awk '{print $2}')
    grep -A $(($NAtoms + 4 )) 'Input orientation:' $name_F.log \
        | tail -n $NAtoms | awk '{print $2,$4,$5,$6}' > TSS.xyz 
    
    # 3. make sure the TSS and the rest structure have the same 
    #   orientation. 
   
}

function changeName(){
    # $1 = TS*.log 
    name=$1
    nPt=$(grep 'Total number of gradient calculations:' $name.log \
        | awk '{print $6}' )
    getIRCstruc $name.log $name.xyz >/dev/null 2>&1
    grep -A $(( $nPt + 2 )) 'Summary of reaction path following' $name.log \
        | tail -n $nPt | awk '{print $2}' > relE.dat 
    E_TS=$(grep 'Energies reported relative to the TS energy of' $name.log \
        | awk '{print $9}' )

    NAtoms=$(head -n 1 $name.xyz)
    jobLine=$(( $NAtoms + 2 ))
    rm -f new.xyz 

    for ((i=1;i<=$nPt;i++))
    do
        relE=$(sed -n "$i,$i p" relE.dat)
        E=$(echo $E_TS + $relE | bc)

        echo $NAtoms >> new.xyz 
        echo $E >> new.xyz 
        head -n $(( $i * $jobLine )) $name.xyz | tail -n $NAtoms >> new.xyz 
    done

    mv new.xyz $name.xyz 
    rm -f relE.dat
}

main $1 $2 