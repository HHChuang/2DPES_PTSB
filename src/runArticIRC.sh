#!/bin/bash
#########################################################################
# Program :                                                             #
#   Interface to integrate Gaussian output, Gaussian 09 and             #
#   GenArticStruc.f90                                                   #
#                                                                       #
# Pre-requisted script/program :                                        #
#   1. getCoord                                                         # 
#   2. GenArticStruc.f90                                                #
#   3. checkGau                                                         #
#                                                                       #
# Input :                                                               #
#   1. *.log                                                            #
#       selected gaussian output file, it's a shoulder in the IRC       #
#       energy profile which is the initial point of the artifical      #
#       reaction coordinate.                                            #
#   2. TS2.log                                                          #
#       extract its eigenvector with one negative eigenvalue            #
#       need freq=hpmodes in the route section                          #
#                                                                       #
# Output :                                                              #
#   1. Artic1D_PEC.dat                                                  #
#   2. Artic1D.xyz                                                      #
#                                                                       #
# Output :                                                              #
#   new structures and their gaussian input/output files                #
#                                                                       #
# How to use it :                                                       #
# 1. copy the source code `genIRCstruc.f90` to the cluster              #
# 2. compile if to binary code                                          #
# 3. using this script to generate structure one-by-one                 #
#                                                                       #
# History :                                                             #
# 2019/05/17, Grace, add the function to automatically stop the program #
# by comparing the energy difference between the last artifical point   #
# and TSS2.                                                             #
# 2019/06/09, Grace, modify function getEandXYZ() since the orientation #
# of TSS2 is not the same to other structures.                          #
#########################################################################

# Main program 
function main(){
    # Record that the input argument is a gaussian output file
    # $1 = *.log; selected point near the VRI regioin

    # 0. stdout purpose and then generate necessary variables
        purpose $1 $2 
        echo 'First point' > Eff.dat 
        name=$(echo $1 | sed 's/.log//g')
        extractCoord $name # output: $1.xyz 
        NAtoms=$(head -n 1 $name.xyz)

    # 1. Re-calculate the cutted point in order to extract gradient 
        (time firstPts $name) 2>>Eff.dat # output: $name_F.log

    # 2. Iteration steps 
    #   Build the artifical reaction coordinate; genModIRCstruc.f90

        getRoute $name
        charge=$(cat charge.dat)
        multi=$(cat multi.dat)

        echo 'Iteration' >> Eff.dat 
        getEigV $NAtoms $2 # output: TS2EigV.dat 
        cp $name\_F.log 0_F.log 

        E_TS2=$(checkGau $2 | tail -n 1 | awk '{print $2}')
        E_error=0.0001 # 0.0001 hartree ~ less than 1 kcal/mol 

        # for i in {1..2} #test
        i=1 ###############
        while true 
        do
            echo Point $i >> Eff.dat
            former=$(($i-1))
            (time newPts $former\_F TS2EigV.dat) 2>>Eff.dat # output: newCoord_F.com 
            mv newCoord_F.com $i\_F.com
            genSGEfile $i\_F $i\_F 
            (time qsub -sync y job$i\_F.sge ) 2>>Eff.dat
            # check the energy to exit the loop
            E_pts=$(checkGau $i\_F.log | tail -n 1 | awk '{print $2}')
            dE=$(echo "scale=6; $E_pts - $E_TS2 " | bc)
            result=$(echo $dE'<'$E_error | bc -l) # 1=true, 0=false
            [ $result -eq 1 ] && break
            i=$(($i+1)) ###
        done 
        echo '' >> Eff.dat 
        echo "The energy difference between the last point and TS2 is $dE hartree." >> Eff.dat 

    # 3. Collect all the structures and the corresponding energy profile
    getEandXYZ $i $2 Artic1D_PEC.dat Artic1D.xyz # output: Artic1D_PEC.dat and Artic1D.xyz
}

function purpose(){
    # $1 = cutted point; *.log 
    # $2 = TS2.log 
    echo ''
    echo 'Generate a serious of structure along artifical '
    echo 'reaction coordinate.'
    echo ''
    echo "Intput: 1. $1, 2. $2"
    echo 'output: 1. Artic1D_PEC.dat, 2. Artic1D.xyz'
    echo ''
}

function extractCoord(){
    # input: $name
    # output: $name.xyz
    getCoord $1.log $1.xyz > /dev/null 2>&1 
}

function getRoute(){
    # input: 
    #   $1 = $name
    # output:
    #   1. route.dat 
    #   2. charge.dat 
    #   3. multi.dat 
    grep \# $1.log | head -n 1 > route.dat
    grep Charge $1.log | head -n 1 | awk '{print $3}' > charge.dat 
    grep Multi $1.log | head -n 1 | awk '{print $6}' > multi.dat 
}

function genPts(){
    # input: 
    # $1 = $name
    # $2 = route.dat 
    # $3 = $charge
    # $4 = $multi
    # output: 
    #   1. $name_F.com 
    charge=$3
    multi=$4

cat << EOF > $1_F.com 
%nprocshared=16
%mem=32GB
`cat $2`
force

calculate force

$charge $multi
`tail -n $NAtoms $1.xyz`

C H O Br 0
6-31g(d)
****
Rh 0
lanl2dz
****

Rh 0
lanl2dz

EOF
}

function getF(){
    # input:
    #   1. NAtoms
    #   2. $name
    # output: F.dat 
    jobLine=$(( $1 + 2 ))
    grep -A $jobLine 'Forces (Hartrees/Bohr)' $2.log \
        | tail -n $1 | awk '{print $3,$4,$5}' > F.dat
}

function getEigV(){
    # input:
    #   1. NAtoms
    #   2. TS2.log 
    # output: TS2EigV.dat 
    line1=$(( 4 + 3*$1 ))
    line2=$(( 3*$1 ))
    grep -A $line1 'Frequencies ---' $2 \
        | tail -n $line2 | awk '{print $4}' > TS2EigV.dat 
    # arrange 1D array to 2D matrix
    for ((i=1;i<=$1;i++))
    do
        head -n $(($i*3)) TS2EigV.dat | tail -n 3 > tmp.dat
        atom_x=$(sed -n '1,1 p' tmp.dat)
        atom_y=$(sed -n '2,2 p' tmp.dat)
        atom_z=$(sed -n '3,3 p' tmp.dat)
        echo $atom_x $atom_y $atom_z >> 2DTS2EigV.dat 
    done
    mv 2DTS2EigV.dat TS2EigV.dat 
    rm -f tmp.dat 
}

function genSGEfile(){
# $1 = file index
# $2 = filename
cat << EOF > job$1.sge
#!/bin/bash
export PATH=\$PATH:\$HOME/bin:/opt/util

### Default SGE Control
#$ -S /bin/sh -w w -j y -cwd 	### Run job through bash shell
#$ -j y							### Join stdout and stderr
#$ -l hostname='q50|q51'			### Resource control
#$ -cwd							### Use current working directory

### SGE Environment
echo '== SGE Environment =='
echo "Working directory is \$SGE_O_WORKDIR"
cd \$SGE_O_WORKDIR

echo 'Job starts'
echo "    Host: \$HOSTNAME"
echo '    Date:' `date`
echo 'Directory:' `pwd` 

### G16 Setup 
export g16root=/opt
export GAUSS_EXEDIR=\$g16root/g16
export GAUSS_SCRDIR=/scratch/$USER/g03
export LD_LIBRARY_PATH=/opt/gcc-8.2.0/lib64:$LD_LIBRARY_PATH:\$g16root/g16:\$g16root/g16/bsd
export PATH=$PATH:\$GAUSS_EXEDIR:\$g16root/g16/bsd

echo g16root is \$g16root
echo g16 scratch is \$GAUSS_SCRDIR

### Job Script
echo "Your job:"
time g16 <  $2.com  >  $2.log
EOF
}

function firstPts(){
    # $1 = $name
    name=$1
    getRoute $name # output: 1. route.dat, 2. charge.dat 3. multi.dat 
    charge=$(cat charge.dat)
    multi=$(cat multi.dat)
    genPts $name route.dat $charge $multi # output: $name_F.com 
        # output variable: $NAtoms
    genSGEfile 0_F $name\_F
    # genSGEfile 0_EigV $name\_EigV
    qsub -sync y job0_F.sge # output: $name_F.log
    # qsub -sync y job0_EigV.sge # output: $name_EigV.sge 
}

function newPts(){
    # $1 = $name
    # $2 = TS2EigV.dat 
    name=$1
    extractCoord $name # output: $name.xyz 
    getF $NAtoms $name # output: F.dat
    GenArticStruc $name.xyz F.dat $2 # output: newCoord.xyz 
    genPts newCoord route.dat $charge $multi # output: newCoord_F.com 
}

function getEandXYZ(){
    # $1 = integer; index of the last point 
    # $2 = TS2.log 
    # $3 = Artic1D_PEC.dat 
    # $4 = Artic1D.xyz

    extractCoord $1\_F # last structures of the artifical rxn. coord.
    TS2name=$(echo $2 | sed 's/.log//g' ) 
    extractCoord $TS2name 
    # Orientation of TSS2 is different, 2019/06/09, Grace
    # in order to execute rot.py, the format of TSS2 need to be modified
    awk '{print $1,$2,$3,$4}' $TS2name.xyz > test
    mv test $TS2name.xyz  

    rm -f $3 $4
    for ((a=0;a<=$1;a++))
    do 
        # get energy profile; output: Artic1D_PEC.dat 
        E=$(checkGau $a\_F.log | grep $a\_F | awk '{print $2}')
        echo $a $E >> $3 
        # get all artifical structures; output: Artic1D.xyz
        cat $a\_F.xyz >> $4
    done
    # rotate TSS2
    rot.py $TS2name.xyz $1\_F.xyz #output: rot.$TS2name.xyz
    mv rot.$TS2name.xyz $TS2name.xyz 

    E=$(checkGau $TS2name.log | grep $TS2name | awk '{print $2}')
    echo $TS2name $E >> $3  
    cat $TS2name.xyz >> $4
}

main $1 $2