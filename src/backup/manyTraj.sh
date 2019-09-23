# 2018/12/14, Grace, the previous one just dispear! fuck
countProg='/work/Grace/Tantillo_Dynamic/H3CO/CountProd'
defR='/work/Grace/Tantillo_Dynamic/H3CO/defR.dat'
defP1='/work/Grace/Tantillo_Dynamic/H3CO/defProd1.dat'
defP2='/work/Grace/Tantillo_Dynamic/H3CO/defProd2.dat'
home=$(pwd)
rawdata="$home/rawdata"

rm -rf RP1 RP2
mkdir RP1
mkdir RP2
RP1dir="$home/RP1"
RP2dir="$home/RP2"
echo $RP1dir

# change directory
cd $rawdata
ls | grep Traj'[0-9]' > trajlist.dat
rm -f record.dat
# for name in `cat trajlist.dat`
for name in Traj1
do
    record=$($countProg $$defR $defP1 $defP2)
    echo $name $record >> record.dat
done

# count and group trajectories
NN=$(grep -c 'none none' record.dat)
RR=$(grep -c 'R R' record.dat)
P1P1=$(grep -c 'P1 P1' record.dat)
P2P2=$(grep -c 'P2 P2' record.dat)
NR=$(grep -c 'none R' record.dat)
NP1=$(grep -c 'none P1' record.dat)
NP2=$(grep -c 'none P2' record.dat)
P1P2=$(grep -c 'P1 P2' record.dat)
RP1=$(grep -c 'R P1' record.dat)
RP2=$(grep -c 'R P2' record.dat)
tot=$(($NN+$RR+$P1P1+$P2P2+$NR+$NP1+$NP2+$P1P2+$RP1+$RP2))

# make RP1 RP2 list and then copy selected trajectories
grep 'R P1' record.dat | awk '{print $1}'  > RP1list.dat
grep 'R P2' record.dat | awk '{print $1}' > RP2list.dat
# rm -f record.dat
for name in `cat RP1`
do 
    echo $name $RP1dir
    cp $name $RP1dir
done
for name in `cat RP2`
do 
    cp $name $RP2dir
done
mv RP1list.dat RP2list.dat $home

cat << EOF
----------------------------------------------------------------

    Analyse the initial/final point of all the trajectories 
    
output:
    files
        1. $home/RP1list.dat
        2. $home/RP2list.dat
    directories
        1. $RP1dir
        2. $RP2dir

std-out: 
total traj.: $tot

OtherOther: $NN
RR: $RR
P1P1: $P1P1
P2P2: $P2P2
OtherR: $NR
OtherP1: $NP1
OtherP2: $NP2
P1P2: $P1P2

RP1: $RP1
RP2: $RP2
----------------------------------------------------------------
EOF

 