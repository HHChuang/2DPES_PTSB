#!/bin/bash
#################################################################
# Program                                                       #
#   Extract one or all eigenvalue(s) and the corresponding      #
#   eigenvector(s).                                             #
#                                                               #
# I/O part                                                      #
#     input:                                                    #
#       $1 = 'filename' or 'all'                                #
#            filename = one specific file name                  #
#            all = all files in the current directory           #
#       $2 = '1' or 'all'                                       #
#           1 = the first mode                                  #
#           all = all negative eigenvalue mode                  #
#                                                               #
#     output:                                                   #
#       1. fileList.dat                                         #
#       2. *.Eig.dat                                            #
#           (format)                                            #
#           # of atoms ; 1st. row                               #
#           'eigenvalue' ; 2nd. row                             #
#           'eigenvector' ; the following rows                  #
#       3. *.G.dat                                              #
#           (format)                                            #
#           # of atoms ; 1st. row                               #
#           'gradient' ; the following rows                     #
#                                                               #
# History                                                       #
#   2018/07/06, Grace                                           #
#   2018/07/10, Grace, function extractAll                      #
#   2018/07/17, Grace, function extractG                        #
#   2018/09/28, Grace, debug of having positive eigenvalue case #
#################################################################

# Functions
function purpose(){
    echo '-------------------------------------------'
    echo 'The input files should have "force" and '
    echo '"freq=hpmodes"'
    echo ''
    echo 'Unit for gradient is (Hartree/Bohr) and the'
    echo 'eigenvalue of Hessian is normalized'
    echo ''
    echo 'output: 1. *.Eig.dat, 2. *.G.dat'
    echo '-------------------------------------------'
}
function userInput(){
    # input:
    #   $1 = file name or all
    # output:
    #   fileList.dat
    case "$1" in
        'all')
            echo 'Analyse all the files in the same directory' 
            ls | grep .log | sed 's/.log//g' | sed 's/n//g' \
                | sort -n -k 1 | sed 's/^/n/' > fileList.dat
        ;;
        *)
            echo "Analyse file:  $1"
            echo $1 | sed 's/.log//g' > fileList.dat
        ;;
    esac
}
function extractOne(){
    # input :
    #   $1 = filename
    # output :
    #   $name.Eig.dat
    NAtoms=$(grep 'NAtoms' $1.log | head -n 1 | awk '{print $2}')
    echo $NAtoms > $1.Eig.dat
    grep 'Frequencies ---' $1.log |head -n 1 \
        |awk '{print $3}' >> $1.Eig.dat
    grep -A $((3*$NAtoms + 4)) 'Frequencies ---' $1.log \
        | head -n $((3*$NAtoms + 5)) | tail -n $((3*$NAtoms)) \
        | awk '{print $4}' >> $1.Eig.dat
}
function extractAll(){
    # input :
    #   $1 = filename
    # output :
    #   $name.Eig.dat
    NAtoms=$(grep 'NAtoms' $1.log | head -n 1 | awk '{print $2}')
    nNegMode=$(grep 'Frequencies ---' $1.log \
        | awk '{print $3,$4,$5,$6,$7}' | grep -o '-' | wc -l) 
    nLine=$(($nNegMode/5))
    nrest=$(($nNegMode%5))
    
    [ $nNegMode == 0 ] && nLine=1 # all positive frequencies

    echo '--' > tmp
    #FIXME: may different for different cases
    grep -A $((3*$NAtoms + 4)) 'Frequencies ---' $1.log >> tmp
    # grep -A $((3*$NAtoms + 7)) 'Frequencies ---' $1.log >> tmp
    head -n $(( ($nLine + 1)*(3*$NAtoms + 6) )) tmp > tmp.dat
    # head -n $(( ($nLine + 1)*(3*$NAtoms + 9) )) tmp > tmp.dat
    rm -f $1.Eig.dat

    for ((a=1;a<=$nLine;a++))
    do
        head -n $(( $a*(3*$NAtoms +6) )) tmp.dat \
            | tail -n $((3*$Natoms +6)) > tmp
        for ((b=1;b<=5;b++))
        do
            echo $NAtoms >> $1.Eig.dat
            grep 'Frequencies ---' tmp \
                | awk "{print $`echo $(($a+2))`}" >> $1.Eig.dat
            tail -n $((3*$NAtoms)) tmp \
                | awk "{print $`echo $(($a+3))`}" >> $1.Eig.dat
        done
    done
    
    for ((a=1;a<=$nrest;a++))
    do
        echo $NAtoms >> $1.Eig.dat
        grep 'Frequencies ---' tmp.dat \
            | awk "{print $`echo $(($a+2))`}" >> $1.Eig.dat
        tail -n $((3*$NAtoms)) tmp.dat \
            | awk "{print $`echo $(($a+3))`}" >> $1.Eig.dat
    done
    rm -f tmp tmp.dat
}

function extractG(){
    # input : 
    #   $1 = filename
    # output: 
    #   $name.G.dat
    NAtoms=$(grep 'NAtoms' $1.log | head -n 1 | awk '{print $2}')
    grep -A $(($NAtoms+2)) 'Forces (Hartrees/Bohr)' $1.log \
        | tail -n $NAtoms | awk '{print $3,$4,$5}' > G.tmp
    echo $NAtoms > $1.G.dat
    for ((a=1;a<=$NAtoms;a++))
    do
        sed -n "$a,$a p" G.tmp | awk '{print $1}' >> $1.G.dat
        sed -n "$a,$a p" G.tmp | awk '{print $2}' >> $1.G.dat
        sed -n "$a,$a p" G.tmp | awk '{print $3}' >> $1.G.dat
    done
    rm -f G.tmp
}

# Main program
    purpose
    userInput $1 # output: fileList.dat
    for name in `cat fileList.dat`
    do
        case "$2" in    
            'all')
                extractG $name # output: $name.G.dat
                extractAll $name # output: $name.Eig.dat
            ;;
            *)
                extractG $name # output: $1.G.dat
                extractOne $name # output: $1.Eig.dat
            ;;
        esac
    done