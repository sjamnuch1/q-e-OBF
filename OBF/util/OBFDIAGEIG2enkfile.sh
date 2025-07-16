#!/bin/bash

if [ $# -le 5 ];then
	echo "Usage $0 eig_directory number_of_kpoint ival fval icond fcond (spin)"
	echo "ival fval icond fcond are index of the eigenvalues"
	echo "They should be the same as on in brange.ipt"
	echo "Ex. metallic system  1 10 10 35"
	echo "Ex. insulator system 1 10 11 35"
	echo "If spin is given then will try to do dw spin too"
        exit
fi

nkpt=`echo $2 | awk '{printf "%d",$1}'`
ival=`echo $3 | awk '{printf "%d",$1}'`
fval=`echo $4 | awk '{printf "%d",$1}'`
icond=`echo $5 | awk '{printf "%d",$1}'`
fcond=`echo $6 | awk '{printf "%d",$1}'`

#echo $nkpt $nval $ncond $ntot

for i in $(seq 1 1 $nkpt )
do
	awk -v ival=$ival -v fval=$fval '{for (i=ival;i<=fval;i++) print $i}' $1/eig${i}.txt | awk '{a[FNR]=$1} END {for (i=1;i<=NR;i++) if ( i%3 == 0 || i == NR) printf "%17.14f\n",a[i] ; else printf "%17.14f ",a[i]}'
	awk -v icond=$icond -v fcond=$fcond '{for (i=icond;i<=fcond;i++) print $i}' $1/eig${i}.txt | awk '{a[FNR]=$1} END {for (i=1;i<=NR;i++) if ( i%3 ==0 || i == NR) printf "%17.14f\n",a[i]; else printf "%17.14f ",a[i]}'
done > enkfile

if ! [ -z $7 ]
then
	for i in $(seq 1 1 $nkpt )
	do
		awk -v ival=$ival -v fval=$fval '{for (i=ival;i<=fval;i++) print $i}' $1/eig${i}dw.txt | awk '{a[FNR]=$1} END {for (i=1;i<=NR;i++) if ( i%3 == 0 || i == NR) printf "%17.14f\n",a[i] ; else printf "%17.14f ",a[i]}'
		awk -v icond=$icond -v fcond=$fcond '{for (i=icond;i<=fcond;i++) print $i}' $1/eig${i}dw.txt | awk '{a[FNR]=$1} END {for (i=1;i<=NR;i++) if ( i%3 ==0 || i == NR) printf "%17.14f\n",a[i]; else printf "%17.14f ",a[i]}'
	done >> enkfile
fi


#for i in {1..1000..1}; do cat ../../../../8Gamma/Out_unk/system.save/eig${i}.txt | awk '{for (i=1;i<=25;i++) print $i}' | awk '{a[FNR]=$1} END {for (i=1;i<=NR;i++) if ( i%3 == 0 || i ==NR) printf "%17.14f\n",a[i] ; else printf "%17.14f ",a[i]}'; done > eigobf.txt
