#!/bin/bash
#########################################################
# History:                                              #
# 2019/03/29, Grace                                     #
# 2019/05/27, Grace, add asymmetric cases               #
#                                                       #
#########################################################

home=`pwd` #/home=***/1D

function main(){
    # userIO Ngrid
    # Ngrid*2-1 = amount of grids for 1D 
    Ngrid=30 # Default grid points: Ngrid=30
    purpose $Ngrid 

    cd $home # PATH=/1D

    # 1. Generate original x- and y-axes; the step 3 in README.md
        genXandY

    # 2. Select grid points
        genSelectXandY
}

function purpose(){
    # $1 = $ngrid
cat << EOF
---------------------------------------------------------------------------------
    Current dir.: $home
    
    For symmetric cases, the input files are: 
        1. header.dat
        2. TS1_F/R.log 
        3. TS2_F/R.log. 
    
    For asymmetric cases, the input files are:
        1. header.dat
        2. TS1_F/R.log
        3. TS2_F/R.log 
        4. \$(coordinate of select point).log 
        5. Artic1D.xyz 

    The amount of selected grids: $1
---------------------------------------------------------------------------------
EOF
}

function userIO(){
    echo ''
    read -p 'Selected amount of grids (integer)  : ' $1 
    echo ''
}

function runG16(){
    # $1 = list.dat
    for name in `cat $1`
    do  
        g16sub $name.com $name.log 
    done
}

function originalPEC(){
    list=$1
    for name in `echo ${list[@]}`
    do
        cd $home
        rm -rf $name
        mkdir $name 

        getIRCcurve $name.log E_$name.dat
        getIRCstruc $name.log $name.xyz
        cp $name.xyz $name
        cp header.dat $name
        cd $home/$name
        writeGauInpV $name.xyz header.dat
        ls | grep .com | sed 's/.com//g' > list.dat 
        runG16 list.dat 
        rm -f list.dat 
    done
}

function asymPEC(){
    # input: 
    #   $1 = TS1_F
    #   $2 = coordinate of select point
    #   $3 = Artic1D.xyz 
    #   $4 = TS1_F_Artic1D
    # output:
    #   TS1_F_Artic.xyz 

    cd $home 
    mkdir $4

    getIRCstruc $1.log $1.xyz 
    lineP=$(grep -n $2 $1.xyz | cut -d ':' -f 1)
    nline=$(($lineP - 2 ))
    head -n $nline $1.xyz > $4.xyz
    rename $3 $2 
    cat $3 >> $4.xyz
    cp $4.xyz $4 
    cp header.dat $4 

    cd $home/$4 #PATH=/1D/TS1_F_Artic1D
    writeGauInpV $4.xyz header.dat
    ls | grep .com | sed 's/.com//g' > list.dat 
    runG16 list.dat 
    rm -f list.dat
}

function collectAsymPEC(){
    # $1 = TS1_F_Artic1D
    
    [ -d $1 ] || echo "No directory $1, check function collectAsymPEC(). Stop program!"
    [ -d $1 ] || exit

    cd $home/$1 
    njobs=$(ls | grep -c log )
    checkGau all | tail -n $njobs | sort -n -k 1 > E_$1.dat
    cp E_$1.dat $home
}

function rename(){
    # $1 = *.xyz 
    # $2 = coordinate of selected point

    # From genModIRCstruc.f90, one step = 0.005 MW
    NAtoms=$(head -n 1 $1)
    nfile=$(( $NAtoms + 2 ))
    nline=$(wc -l $1 | awk '{print $1}')
    njobs=$(( $nline/$nfile ))
    rm -f rename.xyz rename.tmp 
    for ((i=1;i<=$njobs;i++))
    do 
        head -n $(($nfile*$i)) $1 | tail -n $NAtoms > rename.tmp 
        j=$(($i-1))
        coord=$( echo "$2 + 0.005*$j" | bc )
        echo $NAtoms >> rename.xyz 
        echo $coord >> rename.xyz 
        cat rename.tmp >> rename.xyz 
    done
    mv rename.xyz $1
    rm -f rename.tmp 
}

function selectGrid(){
    # $1 = ngrid
    # $2 = list 
    list=$2
    cd $home
    for name in `echo ${list[@]}`
    do
        selectIRCstruc.sh E_$name.dat $name.xyz $1
    done
}

function runSelectPEC(){
    list=$1
    for name in `echo ${list[@]}`
    do
        cd $home
        mkdir $name
        cd $home/$name
        cp ../$name.xyz .
        cp ../header.dat .
        writeGauInpV $name.xyz header.dat 
        ls | grep .com | sed 's/.com//g' > list.dat 
        runG16 list.dat 
        rm -f list.dat 
    done
}

function genXandY(){
    [ -f Artic1D.xyz ] # Artic1D.xyz is exist or not
    case $? in 
        1) # symmetric cases
            list=(TS1_F TS1_R TS2_F TS2_R)
            # list=(TS1_F)
            originalPEC $list
        ;;
        0) # asymmetric cases; Artic1D.xyz is exist 
            [ -f TS1_F_Artic1D.xyz ] 
            case $? in 
                1) # TS1_F_Artic1D.xyz is not exist
                    selectP=$(ls | grep log| grep -v F \
                        | grep -v R| sed 's/.log//g')
                    asymPEC TS1_F $selectP Artic1D.xyz TS1_F_Artic1D #output: TS1_F_Artic1D.xyz
                ;;
                0) # TS1_F_Artic1D.xyz exists
                    collectAsymPEC TS1_F_Artic1D #output: E_TS1_F_Artic1D.dat
                ;;
            esac 
            list=(TS1_R TS2_F TS2_R)
            for name in `echo ${list[@]}`
            do
                cd $home
                [ -d $name ] || originalPEC $name
            done
        ;;
    esac
}

function genSelectXandY(){
    [ -f Artic1D.xyz ] # Artic1D.xyz is exist or not
    case $? in 
        1) # symmetric cases
            FirstFile='TS1_F'
        ;;
        0) # asymmetric cases
            FirstFile='TS1_F_Artic1D'
        ;;
    esac
    list=($FirstFile TS1_R TS2_F TS2_R)
    # list=($FirstFile)
    selectGrid $Ngrid $list 
    list=($Ngrid\_$FirstFile \
          $Ngrid\_TS1_R
          $Ngrid\_TS2_F
          $Ngrid\_TS2_R
        )
    # list=($Ngrid\_$FirstFile)
    runSelectPEC $list
}

main 