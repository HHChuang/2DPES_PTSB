#!/bin/bash
#########################################################################
# Program:                                                      `       #
#   Searching the corresponding coordinates of trajectories.            #
#                                                                       #
# Pre-requist program:                                                  #
#   MapTraj.f90                                                         #
#                                                                       #
# History:                                                              #
# 2019/06/04, Grace, revised                                            #
# 2019/07/23, Grace, change coordinate from integer to floating points  #
#########################################################################

function main(){
    progPath='/Users/Grace/Google_Drive/Code/GitHub/2DPES_PTSB/src/MapTraj'
    sysPath='/Users/Grace/Google_Drive/Code/GitHub/2DPES_PTSB/run/H2O'
    PESdirPath="$sysPath/Traj"
    PES_dir='PESxyz' 
    tsPESPath=$sysPath
    # PES_Struc='59_Struc.G_-0.6.xyz'
    # PES_Energy='59_E.G_-0.6.dat'
    PES_Struc='3_Struc.xyz'
    PES_Energy='3_E.dat'

    iniDir1='RP1.reorder.rot'
    iniDir2='RP2.reorder.rot'
    iniDir3='IRC.reorder.rot'
    iniDir4='test.rot'
    finDir1='RP1.reorder.rot.coord'
    finDir2='RP2.reorder.rot.coord'
    finDir3='IRC.reorder.rot.coord'
    finDir4='test.rot.coord'

    # 2DPES: split xyz files 2D-PES; time consuming step
    # echo 'split xyz files 2D-PES'
    # time splitPESfiles $PESdirPath $PES_dir $tsPESPath $PES_Struc # output: 1. $name.xyz 2. orderList.dat
    time runMapping $PESdirPath $iniDir4 $finDir4 $PES_dir $sysPath $PES_Energy

    # RP1
    # echo 'RP1'
    # time runMapping $PESdirPath $iniDir1 $finDir1 $PES_dir $sysPath $PES_Energy

    # RP2
    # echo 'RP2'
    # time runMapping $PESdirPath $iniDir2 $finDir2 $PES_dir $sysPath $PES_Energy

    # IRC
    # echo 'IRC'
    # time runMapping $PESdirPath $iniDir3 $finDir3 $PES_dir $sysPath $PES_Energy
}

function runMapping(){
    # $1 = $PESdirPath
    # $2 = $iniDir
    # $3 = $finDir
    # $4 = $PES_dir
    # $5 = $sysPath
    # $6 = $PES_Energy #coord: floating point; for int2real()
    
    cd $1/$2
    cp * $1/$3
    cd $1/$4
    cp * $1/$3
    cp $5/$6 $1/$3
    

    cd $1/$3
    ls | grep rot | grep -v coord > nameList

    for name in `cat nameList`
    do
        splitTraj $name
        $progPath $name orderList.dat $6 #output: coord.$name 
        # int2real coord.$name orderList.dat $6 
        # rm -f $name
    done
    # rm -f nameList tss.* forward.* reverse.*
    # rm -f orderList.dat $6 *.xyz 
}

function splitPESfiles(){
    # $1 = $PESPath
    # $2 = $PES_dir
    # $3 = $tsPESPath
    # $4 = $PES_Struc
    # output: 
    #   1. $name.xyz 
    #   2. orderList.dat

    cd $1
    rm -rf $2
    mkdir $2
    cd $2
    cp $3/$4 .

    NAtoms=$(head -n 1 $4)
    oneStruc=$(($NAtoms+2))
    num1D=$(echo $4 | cut -d '_' -f 1)
    num2D=$(($num1D*$num1D))
    rm -f orderList.dat
    for((i=1;i<=$num2D;i++))
    do
        head -n $(($i*$oneStruc)) $4 | tail -n $oneStruc > tmp 
        x=$(sed -n '2,2 p' tmp | awk '{print $1}')
        y=$(sed -n '2,2 p' tmp | awk '{print $2}')
        E=$(sed -n '2,2 p' tmp | awk '{print $3}')
        name=$x\_$y
        cat tmp > $name.xyz 
        echo $x $y $E >> orderList.dat
    done
    rm -f tmp
}

function splitTraj(){
    # input:
    #   $1 = name of trajectory
    # output:
    #   1. tss.$1.xyz
    #   2. forward.$1.xyz
    #   3. reverse.$1.xyz

    NAtoms=$(head -n 1 $1)
    oneStruc=$(($NAtoms+2))
    grep -B 1 -A $NAtoms 'runpoint 1' $1 > tss.$1.xyz
    Ntss=$(grep -n 'runpoint 1' $1 | cut -d ':' -f 1)
    # Reverse 
    head -n $(($Ntss - 2 )) $1 > reverse.xyz
    Nline=$(wc -l reverse.xyz| awk '{print $1}')
    Njobs=$(( $Nline/$oneStruc ))
    for ((i=1;i<=$Njobs;i++))
    do 
        tail -n $(( $i*$oneStruc )) reverse.xyz \
            | head -n $oneStruc > reverse.$i.xyz 
    done

    # Forward
    nFor1=$(($Ntss + $NAtoms + 1))
    nFor2=$(wc -l $1 | awk '{print $1}')
    sed -n "$nFor1,$nFor2 p " $1 > forward.xyz
    Nline=$(wc -l forward.xyz| awk '{print $1}')
    Njobs=$(( $Nline/$oneStruc ))
    for ((i=1;i<=$Njobs;i++))
    do 
        head -n $(( $i*$oneStruc )) forward.xyz \
            | tail -n $oneStruc > forward.$i.xyz 
    done
    rm -f reverse.xyz forward.xyz
}

function int2real(){
    # $1 = name of trajectory 
    # $2 = orderList; integer
    # $3 = $PES_Energy; real 

    # search specific column but return full row 
    # https://stackoverflow.com/questions/27094678/linux-bash-script-how-to-search-on-a-column-but-return-full-row

    totline=$(wc -l $1 | awk '{print $1}')
    rm -f tmp 
    for ((i=1;i<=$totline;i++))
    do 
        buffer=$(sed -n "$i,$i p" $1)
        x=$(echo $buffer | awk '{print $1}')
        y=$(echo $buffer | awk '{print $2}')
        line=$(nl $2 | awk -v col=2 -v term="$x" 'toupper($col)==toupper(term)' \
            | awk -v col=3 -v term="$y" 'toupper($col)==toupper(term)' \
            | awk '{print $1}')
        sed -n "$line,$line p" $3 >> tmp 
    done
    mv tmp $1 
}
main