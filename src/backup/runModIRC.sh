#!/bin/bash
#################################################################
# 1. copy the source code `genIRCstruc.f90` to the cluster      #
# 2. compile if to binary code                                  #
# 3. using this script to generate structure one-by-one         #
#################################################################

function genSGEfile(){
# $1 = filename
cat << EOF > job$1.sge
#!/bin/bash
export PATH=\$PATH:\$HOME/bin:/opt/util

### Default Control
#$ -S /bin/sh 					### Run job through bash shell
#$ -j y							### Join stdout and stderr
#$ -l hostname='q29|q30|q31|q32|q33|q34|q35|q36|q37|q38'  ### Resource control
#$ -cwd							### Use current working directory

### SGE Environment
echo == SGE Environment ==
echo Working directory is \$SGE_O_WORKDIR
cd \$SGE_O_WORKDIR

echo Job starts
echo "     Host: "\$HOSTNAME
echo "     Date: "\`date\`
echo "Directory: "\`pwd\`
### G09 Setup
#
export g09root=/ssd/\$USER/g09/a02
export GAUSS_EXEDIR=\$g09root/g09
export GAUSS_SCRDIR=/scratch/\$USER/g03
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:\$g09root/g09:\$g09root/g09/bsd
export PATH=\$PATH:\$GAUSS_EXEDIR:\$g09root/g09/bsd

echo g09root is \$g09root
echo g09 scratch is \$GAUSS_SCRDIR
echo
### Job Script
#
echo "Your job:"
time g09 <  $1.com  >  $1.log
EOF
}
./genModIRCstruc $1 1.com
genSGEfile 1
qsub -sync y job1.sge
for i in {2..500}
do 
    ./genModIRCstruc $(($i-1)).log $i.com
    genSGEfile $i
    qsub -sync y job$i.sge
done