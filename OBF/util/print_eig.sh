#!/bin/bash 


if [ $# -ne 2 ];then
        echo "Usage QEgrid_withGamma.sh eig.txt num_print"
        exit
fi


awk -v ne=$2 '{for (i=1;i<=ne;i++) print $i}' $1
