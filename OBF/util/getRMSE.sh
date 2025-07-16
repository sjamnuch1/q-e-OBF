#/bin/bash


if [ $# -le 5 ];then
	echo "Usage $0 nkpt ival fval icon fcon eig_suffix (do spin)"
	echo "sh    $0 1000 1 20 21 100 NBASIS100 0" 
        exit
fi
ncond=`echo " $5 $3" | awk '{printf "%d",$1-$2}'`
cdir=$(pwd)

if ! [ -z $7 ]
then	
        dft_name=DFTup
else
	dft_name=DFT
fi

echo "Working on $6 in $cdir"

echo "Tag RMSE(Ryd) RMSE(eV)" >> HEADER

echo "computing RMSE for all bands at each kpt"
while read i 
do
	paste ${dft_name}${i}.txt EIG${i}_${6}.txt | awk '{a+=($1-$2)**2} END{print sqrt(a/FNR)}'
done < List > RMSE_${6}.txt

echo "computing RMSE for occupied bands at each kpt"
while read i
do
        paste ${dft_name}${i}.txt EIG${i}_${6}.txt | head -n $3 | awk '{a+=($1-$2)**2} END{print sqrt(a/FNR)}'
done < List > RMSE_${6}_occ.txt

echo "computing RMSE for unoccupied bands at each kpt"
while read i
do
        paste ${dft_name}${i}.txt EIG${i}_${6}.txt | tail -n ${ncond} | awk  '{a+=($1-$2)**2} END{print sqrt(a/FNR)}'
done < List > RMSE_${6}_unocc.txt

echo "computing RMSE for unoccupied bands at all kpts"
while read i
do
        paste ${dft_name}${i}.txt EIG${i}_${6}.txt | tail -n ${ncond} | awk '{a+=($1-$2)**2} END{print sqrt(a/FNR)}'
done < List | awk -v j=${6} '{a+=($1-$2)**2} END{print j,sqrt(a/FNR),sqrt(a/FNR)*27.2114/2}' >> tmp_unocc

echo "computing RMSE for occupied bands at all kpts"
while read i
do
        paste ${dft_name}${i}.txt EIG${i}_${6}.txt | head -n ${3} | awk '{a+=($1-$2)**2} END{print sqrt(a/FNR)}'
done < List | awk -v j=${6} '{a+=($1-$2)**2} END{print j,sqrt(a/FNR),sqrt(a/FNR)*27.2114/2}' >> tmp_occ

echo "computing RMSE for all bands at all kpts"
while read i
do
        paste ${dft_name}${i}.txt EIG${i}_${6}.txt | awk '{a+=($1-$2)**2} END{print sqrt(a/FNR)}'
done < List | awk -v j=${6} '{a+=($1-$2)**2} END{print j,sqrt(a/FNR),sqrt(a/FNR)*27.2114/2}' >> tmp_all

if ! [ -z $7 ]
then
	dft_name=DFTdw
	echo "Working on down spin for $6 in $cdir"
        echo "Tag RMSE(Ryd) RMSE(eV)" >> HEADER

        echo "computing RMSE for all down spin bands at each kpt"
        while read i
          do
          paste ${dft_name}${i}.txt EIGdw${i}_${6}.txt | awk '{a+=($1-$2)**2} END{print sqrt(a/FNR)}'
done < List > dwRMSE_${6}.txt

echo "computing RMSE for occupied bands at each kpt"
while read i
do
        paste ${dft_name}${i}.txt EIGdw${i}_${6}.txt | head -n $3 | awk '{a+=($1-$2)**2} END{print sqrt(a/FNR)}'
done < List > dwRMSE_${6}_occ.txt

echo "computing RMSE for unoccupied bands at each kpt"
while read i
do
        paste ${dft_name}${i}.txt EIGdw${i}_${6}.txt | tail -n ${ncond} | awk  '{a+=($1-$2)**2} END{print sqrt(a/FNR)}'
done < List > dwRMSE_${6}_unocc.txt

echo "computing RMSE for unoccupied bands at all kpts"
while read i
do
        paste ${dft_name}${i}.txt EIGdw${i}_${6}.txt | tail -n ${ncond} | awk '{a+=($1-$2)**2} END{print sqrt(a/FNR)}'
done < List | awk -v j=${6} '{a+=($1-$2)**2} END{print j,sqrt(a/FNR),sqrt(a/FNR)*27.2114/2}' >> tmp_unocc_dw

echo "computing RMSE for occupied bands at all kpts"
while read i
do
        paste ${dft_name}${i}.txt EIGdw${i}_${6}.txt | head -n ${3} | awk '{a+=($1-$2)**2} END{print sqrt(a/FNR)}'
done < List | awk -v j=${6} '{a+=($1-$2)**2} END{print j,sqrt(a/FNR),sqrt(a/FNR)*27.2114/2}' >> tmp_occ_dw

echo "computing RMSE for all bands at all kpts"
while read i
do
        paste ${dft_name}${i}.txt EIGdw${i}_${6}.txt | awk '{a+=($1-$2)**2} END{print sqrt(a/FNR)}'
done < List | awk -v j=${6} '{a+=($1-$2)**2} END{print j,sqrt(a/FNR),sqrt(a/FNR)*27.2114/2}' >> tmp_all_dw

fi
