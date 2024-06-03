g++ -Wl,-rpath=/usr/local/gurobi/linux64/lib -I/usr/local/gurobi/linux64/include -std=c++20 MVBPMP_LR.cpp -L /usr/local/gurobi/linux64/lib -lgurobi_g++5.2 -lgurobi100 -o mvbpmp_lr.x
./mvbpmp_lr.x /users4/deptadm/emis/ydong/data/t20_06_data.txt /users4/deptadm/emis/ydong/mvbpmp-cpp-gurobi/vehicle_location_data/vehicle_data_20_3_1.txt


