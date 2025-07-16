#!/bin/bash 


if [ $# -ne 3 ];then
        echo "Usage QEgrid_withGamma.sh obf_basis.in eig.txt output_file_name"
        exit
fi


nv=`grep 'val' $1 | sed -e "s/=/ /g" | awk '{print $NF}'`
nc=`grep 'con' $1 | sed -e "s/=/ /g" | awk '{print $NF}'`
ntot_e=`echo "$nv $nc" | awk '{printf "%d",$1+$2}'`

#echo $nv $nc $ntot_e

awk -v ne=$ntot_e '{for (i=1;i<=ne;i++) print $i*27.2114}' $2 > $3
