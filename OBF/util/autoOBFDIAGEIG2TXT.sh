#!/bin/bash 


if [ $# -le 3 ];then
	echo "Usage QEgrid_withGamma.sh nprint eig_dir output_file_postfix num_eig_file (do spin polarized if argument exists)"
	echo "Example $0 25 Out_unk/system.save SBAND_0.01 1000 0"
        exit
fi

#nv=`grep 'num_val' $1 | sed -e "s/=/ /g" | awk '{print $NF}'`
#nc=`grep 'num_cond' $1 | sed -e "s/=/ /g" | awk '{print $NF}'`
#ntot_e=`echo "$nv $nc" | awk '{printf "%d",$1+$2}'`
echo "Working on tag ${3}"
#echo $nv $nc $ntot_e
for i in $(seq 1 1 $4)   
do
	echo $i
awk -v ne=$1 '{for (i=1;i<=ne;i++) print $i}' $2/eig${i}.txt > EIG${i}_$3.txt
done

if ! [ -z $5 ]
then
	echo "Also doing dw spin"
	for i in $(seq 1 1 $4)
	do
		echo $i
		awk -v ne=$1 '{for (i=1;i<=ne;i++) print $i}' $2/eig${i}dw.txt > EIGdw${i}_$3.txt
	done
fi
