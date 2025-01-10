#!/bin/bash

dataFolder="data"

startSize=20
endSize=20

apdx="test"

rm gurobiTest.x

g++ -g -Wl,-rpath=/usr/local/gurobi/linux64/lib -I/usr/local/gurobi/linux64/include -std=c++20 gurobiBugTest.cpp -L /usr/local/gurobi/linux64/lib -lgurobi_g++8.5 -lgurobi110 -o gurobiTest.x

for((i=$startSize;i<=$endSize;i+=10))
do
for((j=6;j<=6;j+=1))
do
if [ "$j" -le 9 ]; then
		datafile="t"$i"_0"$j"_data.txt"
	else
		datafile="t"$i"_10_data.txt"
fi
outputFile=dw-cg_"$i"_"$j"_"$dataFolder"_"$apdx".out 

./gurobiTest.x /users4/deptadm/emis/ydong/"$dataFolder"/$datafile > $outputFile 2>&1 

wait
done
done





