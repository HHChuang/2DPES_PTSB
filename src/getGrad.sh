# !/bin/bash
#################################################################
# Program:                                                      #
#   Extract atomic number and forces from gaussian output files #
#                                                               #
# Output:                                                       #
#   1. a serious files name as $name.grad                       #
#   2. AN.dat                                                   #
#                                                               #
# History:                                                      #
#   2019/06/23, Grace                                           #
#################################################################

ls | grep .log | sed 's/.log//g' > list
name=$(head -n 1 list)
NAtoms=$(head -n 1 $name.xyz)
line=$(($NAtoms + 2))

tail -n $NAtoms $name.xyz | awk '{print $1}' > AN.dat

# gradient is equal to negative force
rm -f *.grad
for name in `cat list`
do 
    grep -A $line 'Forces (Hartrees/Bohr)' $name.log \
        | tail -n $NAtoms | awk '{print -$3, -$4, -$5}' > $name.grad
done
rm -f list