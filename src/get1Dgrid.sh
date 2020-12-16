#!/bin/bash
#####################################################################
# Program:                                                          #
#                                                                   #
#   pre-requisited scripts:                                         #
#       1. checkGau                                                 #                           #                                                                   #
#       2. rev1Dstruc                                               #
#       3. rot.py                                                   #
#                                                                   #
# History:                                                          #
# 2019/03/27, Grace                                                 #
# 2019/05/31, Grace, add asymmetric case                            #
# 2019/06/10, Grace, rotate the TSS2 in x.xyz in order to align all #
#   the structures. Also, use this rotated TSS2 to rotate the IRC of#
#   TS2.                                                            #
#####################################################################

home=`pwd`

function main(){
    purpose

    # Step 1. check if it is a symmetric/asymmetric case,
    # and get the amount of grid.
        # if Artic1D.xyz exist, it is a asymmetric case
        # variables: $TS1_F_SAS, $Ngrid
        getTS1_F_SAS_Ngrid

    # Step 2. Extract energy profile for original cases.
        # list=( $TS1_F_SAS TS1_R TS2_F TS2_R)
        list=( $TS1_F_SAS )
        # echo ${list[@]}
        getPEC $list # output: E_scan.*.dat 
        Ngrid=20

    # Step 3. Extract energy profile for the cases
    # with selected grid points
        list=( $Ngrid\_$TS1_F_SAS \
                $Ngrid\_TS1_R \
                $Ngrid\_TS2_F \
                $Ngrid\_TS2_R
                )
        list=( $Ngrid\_$TS1_F_SAS)
        getPEC $list # output: E_scan.*.dat 

        # x.xyz: Combine TS1 forward direction and reverse direction 
        getXandY $TS1_F_SAS $Ngrid TS2_F TS2_R # output: x.xyz 
}

function purpose(){
    cat << EOF
---------------------------------------------------------------------------------
    Negative sign is replaced by 'n', since the file name cannot 
    start by '-'. 

    output: 
        1. x.xyz 
        2. E_scan.*.dat
---------------------------------------------------------------------------------
EOF
}

function getTS1_F_SAS_Ngrid(){
    [ -f Artic1D.xyz ]
    case $? in 
        1) # symmetric cases
            TS1_F_SAS=TS1_F
        ;;
        0) # asymmetric cases
            TS1_F_SAS=TS1_F_Artic1D
        ;;
    esac
    Ngrid=$(ls -d */| grep TS1_F| head -n 1 | cut -d '_' -f 1)
}

function getPEC(){
    # $1 = $list
    list=$1
    for file in `echo ${list[@]}`
    do 
        cd $home/$file
        totjob=$(ls | grep -c .com )
        checkGau all | tail -n $totjob | sort -n -k 1 > E_scan.$file.dat 
        mv E_scan.$file.dat $home
    done
}

function getXandY(){
    # input:
    #   $1 = $TS1_F_SAS 
    #   $2 = $Ngrid
    #   $3 = TS2_R
    #   $4 = TS2_F
    # output:
    #   1. x.xyz
    #   2. $NGrid\_TS2_R.xyz (modified)
    #   3. $NGrid\_TS2_F.xyz (modified)

    cd $home
    rev1Dstruc $2\_TS1_R.xyz $2\_rev_TS1_R.xyz 
    NAtoms=$(head -n 1 $2\_rev_TS1_R.xyz)
    nlines=$(wc -l $2\_rev_TS1_R.xyz| awk '{print $1}')
    removeTS=$(( $nlines-($NAtoms+2) ))
    head -n $removeTS $2\_rev_TS1_R.xyz  > x.xyz

    # Step 1. Rotate the last point of $Ngrid_$TS1_F_SAS.xyz
    tail -n $(($NAtoms + 2)) $2\_$1.xyz \
        | awk '{print $1,$2,$3,$4}' > tss2.x.xyz
    tail -n $((2 * ($NAtoms+2) )) $2\_$1.xyz| \
        head -n $(($NAtoms + 2 )) | awk '{print $1,$2,$3,$4}' \
        > tss2_plusOne.x.xyz

    ~/bin/rot.py tss2.x.xyz tss2_plusOne.x.xyz tss2.x.xyz 
    # output: rot.tss2.x.xyz
    
    mv rot.tss2.x.xyz tss2.x.xyz
    nlines=$(wc -l $2\_$1.xyz| awk '{print $1}')

    # reserve the name of coordinate of TSS2 such 
    # that only remove the $NAtoms line
    head -n $(($nlines - $NAtoms)) $2\_$1.xyz > tmp 
    
    tail -n $NAtoms tss2.x.xyz >> tmp 
    mv tmp $2\_$1.xyz 
    rm -f tss2_plusOne.x.xyz

    # Step 2. Combine the modified forward direction of IRC from TS1 to x.xyz
    cat $2\_$1.xyz >> x.xyz

    # Step 3. Align y from x 
    head -n $(($NAtoms + 2)) $2\_$3.xyz \
        | awk '{print $1,$2,$3,$4}'> tss2.$2\_$3.xyz
    ~/bin/rot.py tss2.$2\_$3.xyz tss2.x.xyz $2\_$3.xyz
    ~/bin/rot.py tss2.$2\_$3.xyz tss2.x.xyz $2\_$4.xyz
    mv rot.$2\_$3.xyz $2\_$3.xyz
    mv rot.$2\_$4.xyz $2\_$4.xyz
    rm -f tss2.x.xyz tss2.$2\_$3.xyz
}

main