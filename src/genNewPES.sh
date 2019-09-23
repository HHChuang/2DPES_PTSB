# !/bin/bash
#################################################################
# Program:                                                      #
#   Generate one or more than one new PESs in G_*.xyz files     #
#                                                               #
# Pre-requiste script:                                          #
#   1. GenGradStruc.f90                                         #
#                                                               #
# Output:                                                       #
#   1. G_*.xyz                                                  #
#                                                               #
# History:                                                      #
#   2019/06/27, Grace                                           #
#################################################################


num1D=$(ls | grep log | cut -d '_' -f 1 \
    | sort -n | uniq | tail -n 1)
range_ds=$(ls| grep xyz | cut -d '_' -f 3 \
    | sed 's/.xyz//g' | sed '/^$/d' | sort -n | uniq )

modX_i=1
modX_f=$(($num1D - 1))
modY_i=1
modY_f=$(( ($num1D - 1) / 2 ))
X_i=$modX_i
X_f=$num1D
Y_i=$(($modY_i - 1))
Y_f=$modY_f
name=$(echo $X_i\_$Y_f.xyz)
NAtoms=$(head -n 1 $name)

for ds in `echo $range_ds`
do 
    rm -f G_$ds.xyz 
    for ((x=$X_i;x<=$X_f;x++))
    do
        for ((y=$Y_f;y>=$Y_i;y=y-1))
        do 
            if [ "$x" == "$X_f" -o "$y" == "$Y_i" ] #-o "$y" == "$Y_f" ] 
                then 
                    name=$x\_$y.xyz
                else 
                    name=$x\_$y\_$ds.xyz 
            fi
            echo $NAtoms >> G_$ds.xyz 
            echo $x\_$y >> G_$ds.xyz 
            cat $name | tail -n $NAtoms | awk '{print $1,$2,$3,$4}' >> G_$ds.xyz 
        done
    done
done

