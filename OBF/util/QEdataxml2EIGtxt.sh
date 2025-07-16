#/bin/bash

  [ $# -ne 1 ] && echo usage: $0 data.xml
  [ $# -ne 1 ] && exit

lspin=`grep 'lsda' $1 | head -n 1 | grep 'true' -c`

grep 'eigenvalues' $1  -n | awk '{a[FNR]=$1} END{for(i=1;i<=FNR/2;i++) printf "%d %d\n", a[2*i-1]+1,a[2*i]-1}' > Line
numeig=`wc -l Line | awk '{print $1}'`

echo "Found $numeig kpts"

for i in $(seq 1 1 ${numeig})
do
	echo $i
j=`sed -n " $i p" Line | awk '{print $1}'`
k=`sed -n " $i p" Line | awk '{print $2}'`
sed -n " $j , $k p" $1 | awk '{for(i=1;i<=NF;i++) print $i}' | awk '{printf "%12.8f\n",$1*2}' > DFT${i}.txt
done
rm Line

#If we do spin polarized calculation, eigenvalues block will contain up->dw eigenvalues so we have
#nbnd*2 elements in the DFT.txt file.
#Here we divide them into up/dw file
if [ $lspin == 1 ]
then
nbnd=`wc -l DFT1.txt | awk '{print $1/2}'`
for i in $(seq 1 1 ${numeig})
do
	head -n ${nbnd} DFT${i}.txt > DFTup${i}.txt
	tail -n ${nbnd} DFT${i}.txt > DFTdw${i}.txt
	rm DFT${i}.txt
done
fi
