#!/bin/bash

#dataFolder="newData_nodes_5-50"
dataFolder="data"
#dataFolder="sparse_data"
#dataFolder="sparse_data_20pct_of_original"

startSize=10
endSize=20

#****** it means use genetic algorithm (GA) calling Gurobi (GRB) to calculate fitness, and use GA as initial solution (init-Sol) to Danzig-Wolfe Column Generation algorithm
#apdx="GA-GRB-init-Sol"
#apdx="binary-y_last-lambda"
apdx="dominance"
#apdx="gurobi"
#apdx="test"

statfile="stat_"$startSize"_"$endSize"_"$dataFolder"_"$apdx".txt"
echo "" >  $statfile

#g++ -Wl,-rpath=/usr/local/gurobi/linux64/lib -I/usr/local/gurobi/linux64/include -std=c++20 MVBPMP_LR.cpp -L /usr/local/gurobi/linux64/lib -lgurobi_g++5.2 -lgurobi100 -o mvbpmp_lr.x

rm bpmp_dw-cg.x

g++ -g -Wl,-rpath=/usr/local/gurobi/linux64/lib -I/usr/local/gurobi/linux64/include -std=c++20 BPMP_DW-CG.cpp -L /usr/local/gurobi/linux64/lib -lgurobi_g++8.5 -lgurobi110 -o bpmp_dw-cg.x
#g++ -Wl,-rpath=/usr/local/gurobi/linux64/lib -I/usr/local/gurobi/linux64/include -std=c++20 BPMP_DW-CG_old.cpp -L /usr/local/gurobi/linux64/lib -lgurobi_g++8.5 -lgurobi110 -o bpmp_dw-cg.x
#g++ -Wl,-rpath=/usr/local/gurobi/linux64/lib -I/usr/local/gurobi/linux64/include -std=c++20 BPMP_DW-CG_oldestVersion.cpp -L /usr/local/gurobi/linux64/lib -lgurobi_g++8.5 -lgurobi110 -o bpmp_dw-cg.x


for((i=$startSize;i<=$endSize;i+=10))
do
for((j=1;j<=10;j+=1))
do
if [ "$j" -le 9 ]; then
		datafile="t"$i"_0"$j"_data.txt"
	else
		datafile="t"$i"_10_data.txt"
fi
outputFile=dw-cg_"$i"_"$j"_"$dataFolder"_"$apdx".out 
./bpmp_dw-cg.x /users4/deptadm/emis/ydong/"$dataFolder"/$datafile 3 /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_"$i"_3_1.txt > $outputFile 2>&1 
#echo "./bpmp_dw-cg.x /users4/deptadm/emis/ydong/"$dataFolder"/$datafile 3 /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_"$i"_3_1.txt" 
wait
echo "====== problem size $i  $j ======" >> $statfile
#tail -4 dw-cg_"$i"_"$j$apdx".out >> $statfile 
tail -6 $outputFile >> $statfile 
echo "" >> $statfile
done
done

#echo " " | mail -s "finished $i" -a ga_"$i"_output.txt -a pro_"$i".txt ydong




