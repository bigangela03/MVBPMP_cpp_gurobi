#g++ -Wl,-rpath=/usr/local/gurobi/linux64/lib -I/usr/local/gurobi/linux64/include -std=c++20 MVBPMP_LR.cpp -L /usr/local/gurobi/linux64/lib -lgurobi_g++5.2 -lgurobi100 -o mvbpmp_lr.x

#g++ -Wl,-rpath=/usr/local/gurobi/linux64/lib -I/usr/local/gurobi/linux64/include -std=c++20 BPMP_LR.cpp -L /usr/local/gurobi/linux64/lib -o bpmp_lr.x

g++ -Wl,-rpath=/usr/local/gurobi/linux64/lib -I/usr/local/gurobi/linux64/include -std=c++20 BPMP_LR.cpp -L /usr/local/gurobi/linux64/lib -lgurobi_g++8.5 -lgurobi110 -o bpmp_lr.x

#g++ -std=c++20 BPMP_LR.cpp -o bpmp_lr.x

#=========== select an instance to run ==========
./bpmp_lr.x /users4/deptadm/emis/ydong/data/t10_01_data.txt 3 /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_10_3_1.txt

#./bpmp_lr.x /users4/deptadm/emis/ydong/data/t10_02_data.txt 3 /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_10_3_1.txt

#./bpmp_lr.x /users4/deptadm/emis/ydong/data/t10_03_data.txt 3 /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_10_3_1.txt

#./bpmp_lr.x /users4/deptadm/emis/ydong/data/t10_04_data.txt 3 /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_10_3_1.txt

#./bpmp_lr.x /users4/deptadm/emis/ydong/data/t10_05_data.txt 3 /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_10_3_1.txt

#./bpmp_lr.x /users4/deptadm/emis/ydong/data/t10_07_data.txt 3 /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_10_3_1.txt

#./bpmp_lr.x /users4/deptadm/emis/ydong/data/t20_01_data.txt 3 /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_20_3_1.txt

#./bpmp_lr.x /users4/deptadm/emis/ydong/data/t20_06_data.txt 3 /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_20_3_1.txt

