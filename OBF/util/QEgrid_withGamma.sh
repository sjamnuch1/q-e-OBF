#!/bin/bash
#Write out kpoint list same way as QE but will include addition 7 gamma points to the list

if [ $# -ne 6 ];then
	echo "Usage QEgrid_withGamma.sh NK1 NK2 NK3 k1 k2 k3"
	exit
fi

NK1=$1
NK2=$2
NK3=$3
k1=$4
k2=$5
k3=$6

echo "K_POINTS crystal"
echo "$NK1 $NK2 $NK3" | awk '{printf "%d\n",$1*$2*$3+7}'

for i in $(seq 1 1 $NK1)
do
	for j in $(seq 1 1 $NK2)
	do
		for k in  $(seq 1 1 $NK3)
		do
			echo "$i $j $k" | awk -v NK1=$NK1 -v NK2=$NK2 -v NK3=$NK3 -v k1=$k1 -v k2=$k2 -v k3=$k3 '{printf "%8.5f %8.5f %8.5f %4.1f\n",($1-1)/NK1+k1/NK1,($2-1)/NK2+k2/NK2 ,($3-1)/NK3+k3/NK3 ,1}'
		done
	done
done
echo " 0.00000  1.00000  1.00000  1.0"
echo " 1.00000  0.00000  1.00000  1.0"
echo " 1.00000  1.00000  0.00000  1.0"
echo " 1.00000  0.00000  0.00000  1.0"
echo " 0.00000  1.00000  0.00000  1.0"
echo " 0.00000  0.00000  1.00000  1.0"
echo " 1.00000  1.00000  1.00000  1.0"
