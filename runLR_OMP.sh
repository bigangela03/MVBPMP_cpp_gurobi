#========== compile the source code ==========
#g++ -Wl,-rpath=/usr/local/gurobi/linux64/lib -I/usr/local/gurobi/linux64/include -std=c++20 -fopenmp MVBPMP_LR_openmp.cpp -L /usr/local/gurobi/linux64/lib -lgurobi_g++5.2 -lgurobi100 -o mvbpmp_lr_omp.x

g++ -Wl,-rpath=/usr/local/gurobi/linux64/lib -I/usr/local/gurobi/linux64/include -std=c++20 -fopenmp MVBPMP_LR_openmp.cpp -L /usr/local/gurobi/linux64/lib -lgurobi_g++8.5 -lgurobi110 -o mvbpmp_lr_omp.x

#=========== select an instance to run ==========
#./mvbpmp_lr_omp.x /users4/deptadm/emis/ydong/data/t10_01_data.txt 3 /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_10_3_1.txt

#./mvbpmp_lr_omp.x /users4/deptadm/emis/ydong/data/t10_02_data.txt 3 /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_10_3_1.txt

#./mvbpmp_lr_omp.x /users4/deptadm/emis/ydong/data/t10_03_data.txt 3 /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_10_3_1.txt

#./mvbpmp_lr_omp.x /users4/deptadm/emis/ydong/data/t10_04_data.txt 3 /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_10_3_1.txt

#./mvbpmp_lr_omp.x /users4/deptadm/emis/ydong/data/t10_05_data.txt 3 /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_10_3_1.txt

./mvbpmp_lr_omp.x /users4/deptadm/emis/ydong/data/t10_07_data.txt 3 /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_10_3_1.txt

#./mvbpmp_lr_omp.x /users4/deptadm/emis/ydong/data/t20_01_data.txt 3 /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_20_3_1.txt

#./mvbpmp_lr_omp.x /users4/deptadm/emis/ydong/data/t20_06_data.txt 3 /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_20_3_1.txt


