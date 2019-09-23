#!/bin/bash
#########################################################################
# Program:                                                              #
#                                                                       #
# Pre-requisted script:                                                 #
#           1. getCoord                                                 #
#                                                                       #
# Input:                                                                #
#           1. header.dat                                               #
#           2. x.xyz                                                    #
#           3. y_F.xyz                                                  #
#           4. y_R.xyz                                                  #
#                                                                       #
# Output :                                                              #
#           1. $Pts_E.dat                                               #
#           2. $Pts_Struc.xyz                                           #
#           format:                                                     #
#                       # of atom                                       #
#                       coord_x coord_y E                               #
#                       #atom   x   y   z                               #
#                                                                       #
#   History:                                                            #
# 2018/10/13, Grace                                                     #
# 2019/05/31, Grace, modify its structure                               #
# 2019/07/12, Grace, add two more input files in the command-line in    #
#   order to change the unit into mass-weighted coordinate.             #
# 2019/07/24, Grace, fix the bug in the reverse part of PES             #
# 2019/08/05, Grace, modify function getE() and Etmp() to extract       #
#   post-HF energy. Like MP2, CISD and CCSD.                            #
# 2019/09/19, Grace, fix the bu at Etmp_R.dat in getE() function        #
#########################################################################

home=`pwd`
file=(Forward Reverse)

function main (){
    # Input:  
    #   1. header.dat                                                                
    #   2. x.xyz                                                    
    #   3. y_F.xyz                                                  
    #   4. y_R.xyz                                                  

    cd $home 
    getPts $2 # output variable: $Pts, $NAtoms
    purpose $Pts

    # Extract energy profile 
    getE $1 $Pts $2 $3 $4 # output: Etmp_F.dat, Etmp_R.dat, $Pts_E.dat

    # Extract structure 
    getStruc $Pts # output: Stmp_F.xyz, Stmp_R.xyz, $Pts_Struc.xyz
}


function purpose(){
    # input: 
    #   $1 = $Pts
cat << EOF
--------------------------------------------------------------------------
    Current folder: $home
    It should have two sub-folders, "Forward" and "Reverse"

    Output: 
    1. $1_E.dat
    2. $1_Struc.xyz
--------------------------------------------------------------------------
EOF
}

function getPts(){
    # input: 
    #   $1 = x.xyz 
    # output:
    #  $Pts
    #  $NAtoms

    cd $home
    NAtoms=$(head -n 1 $1 )
    nlines=$(wc -l $1 | awk '{print $1}')
    Pts=$(( $nlines/($NAtoms+2) ))
}

function Etmp(){
    # input: 
    #   $1 = $Method
    #   $2 = name of directory
    #   $3 = $Pts 
    #   $4 = x.xyz 
    #   $5 = y_F.xyz or y_R.dat 
    #   $6 = output name; Etmp.dat 

    cp $4 $5 $2
    cd $home/$2

    Method=$1
    case $Method in 
        'SCF Done') # HF and DFT, but other methods may be included 
            grep 'SCF Done' *.log \
                | awk '{print $1,$6}' \
                | sed 's/.log://g' | sed 's/_/ /g' | sort -n -k 1 -k 2 > $6
        ;;
        'EUMP2') # post-HF method: MP2
            grep 'EUMP2' *.log \
                | awk '{print $1,$7}' \
                | sed 's/.log://g' | sed 's/_/ /g' | sort -n -k 1 -k 2 > $6
        ;;
        'E(Corr)')
            grep 'Wavefunction amplitudes converged. E(Corr)' *.log \
                | awk '{print $1,$6}' \
                | sed 's/.log://g' | sed 's/_/ /g' | sort -n -k 1 -k 2 > $6
        ;;
    esac
    # 2019/07/12, Grace, change unit from # of grid to mass-weighted coordinate
    NAtoms=$(head -n 1 $4)
    rm -f coord_x.dat 
    for ((i=1;i<=$3;i++))
    do
        lines=$(( 2 + ($NAtoms +2) * ($i -1) ))
        sed -n "$lines,$lines p" $4 >> coord_x.dat
    done
    y_half=$(( ($3 + 1)/2 ))
    rm -f coord_y.dat 
    for ((i=1;$i<=$y_half;i++))
    do
        lines=$(( 2 + ($NAtoms +2) * ($i -1) ))
        sed -n "$lines,$lines p" $5 >> coord_y.dat
    done  
    sed -i 's/n/-/g' coord_x.dat 
    sed -i 's/n//g' coord_y.dat
    rm -f coord.dat 
    for ((i=1;i<=$3;i++))
    do
        x=$(sed -n "$i,$i p" coord_x.dat)
        for ((j=1;j<=$y_half;j++))
        do  
            y=$(sed -n "$j,$j p" coord_y.dat)
            echo $x $y >> coord.dat 
        done
    done 
    awk '{print $3}' $6 > E_tmp
    paste coord.dat E_tmp | awk '{print $1,$2,$3}' > $6
    rm -f $4 $5 coord.dat coord_x.dat coord_y.dat E_tmp
    cp $6 $home
    cd  $home
}

function getE(){
    # input:
    #   $1 = header.dat  
    #   $2 = $Pts
    #   $3 = x.xyz 
    #   $4 = y_F.xyz 
    #   $5 = y_R.xyz 
    # output: 
    #   1. Etmp_F.dat 
    #   2. Etmp_R.dat 
    #   3. $Pts_E.dat 

    cd $home
    # Check the level of theory 
    nHF=$(grep -c HF $1 )
    nDFT=$(grep -v HF $1 | grep -v MP2 | \
        grep -v CISD | grep -v CCSD | wc -l )
    if [ "$nHF" == 1 -o "$nDFT" == 1 ]; then 
        Method='SCF Done'=1
    else
        # post-HF method; MP2, CISD
        nMP2=$(grep -v MP2 $1)
        if [ "$nMP2" == 1 ]; then 
            Method='EUMP2'
        else 
            Method='E(Corr)'
        fi
    fi
    
    # Extract energy profile
    Etmp $Method ${file[0]} $2 $3 $4 Etmp_F.dat 
    Etmp $Method ${file[1]} $2 $3 $5 Etmp_R.dat 

    # Remove redundant part; y=0, for Etemp_R.dat,
    # and add negative sign
    cd $home
    nlines=$(wc -l Etmp_R.dat | awk '{print $1}')
    sort -rn -k 2 Etmp_R.dat| head -n $(( $nlines - $2 )) \
        | awk '{print $1,-$2,$3}' > tmp
    mv tmp Etmp_R.dat 
    # 2019/07/24, Grace, since the coordinate is not the number 
    #   of grid point anymore, the sorting method msut be fixed.
    # oneDpts=$(awk '{print $2}' Etmp_R.dat | sort -n | uniq | wc -l )
    # nGroup=$(($nlines/$oneDpts))
    # restPart=$(($oneDpts-1)) #remove the central part; x-axis
    # rm -f tmp
    # for ((i=1;i<=$nGroup;i++))
    # do 
    #     head -n $(($i*$oneDpts)) Etmp_R.dat | tail -n $restPart >> tmp
    # done
    # awk '{print $1,-$2,$3}' tmp > Etmp_R.dat 
    # rm -f tmp
    
    # Combine forward part and reverse part
    cat Etmp_F.dat > $2\_E.dat
    cat Etmp_R.dat >> $2\_E.dat 
    sort -n -k 1 -k 2 $2\_E.dat > tmp
    mv tmp $2\_E.dat 
}

function Stmp(){
    # input: 
    #   $1 = $Pts
    #   $2 = name of directory
    #   $3 = output name; Stmp.xyz
    
    cd $home/$2
    rm -f $3
    ymax=$(( ($1 -1)/2 ))
    for ((x=1;x<=$1;x++))
    do
        for ((y=0;y<=$ymax;y++))
        do
            name=$(echo $x\_$y)
            ~/bin/getCoord $name.log $name.xyz 1> /dev/null
            cat $name.xyz >> $3
        done
    done
    cp $3 $home
}

function getStruc(){
    # input:
    #   $1 = $Pts
    # output:
    #   1. Stmp_F.xyz
    #   2. Stmp_R.xyz
    #   3. $Pts_Struc.xyz    

    # Extract structures
    Stmp $Pts ${file[0]} Stmp_F.xyz
    Stmp $Pts ${file[1]} Stmp_R.xyz 

    # Reorder structures 
    cd $home
    halfGroup=$(( ($Pts+1)/2 ))
    LineOneM=$(($NAtoms+2)) # number of line for one molecule
    LineHalfG=$(( $LineOneM*$halfGroup ))
    rm -f $Pts\_Struc.xyz 
    for ((i=1;i<=$Pts;i++))
    do
        # Reverse part 
        head -n $(( $i*$LineHalfG )) Stmp_R.xyz | tail -n $LineHalfG > HG.tmp
        for ((j=1;j<=$(($halfGroup-1));j++))
        do
            tail -n $(($j*$LineOneM)) HG.tmp| head -n $LineOneM > OM.tmp
            # add negative sign in the name of coordinate 
            X_Y=$(grep log OM.tmp | awk '{print $2}' \
                | sed 's/.log//g' | sed 's/_/ /g' | awk '{print $1,-$2}')
            E=$(grep log OM.tmp | awk '{print $4}')
            echo $NAtoms >> $Pts\_Struc.xyz 
            echo $X_Y $E >> $Pts\_Struc.xyz 
            tail -n $NAtoms OM.tmp >> $Pts\_Struc.xyz 
        done
        # Forward part
        head -n $(( $i*$LineHalfG )) Stmp_F.xyz | tail -n $LineHalfG > HG.tmp
        for ((j=1;j<=$halfGroup;j++))
        do
            head -n $(($j*$LineOneM)) HG.tmp | tail -n $LineOneM > OM.tmp
            X_Y=$(grep log OM.tmp | awk '{print $2}' \
                | sed 's/.log//g' | sed 's/_/ /g' | awk '{print $1,$2}')
            E=$(grep log OM.tmp | awk '{print $4}')
            echo $NAtoms >> $Pts\_Struc.xyz 
            echo $X_Y $E >> $Pts\_Struc.xyz 
            tail -n $NAtoms OM.tmp >> $Pts\_Struc.xyz 
        done
    done
    rm -f HG.tmp OM.tmp
}

main $1 $2 $3 $4
