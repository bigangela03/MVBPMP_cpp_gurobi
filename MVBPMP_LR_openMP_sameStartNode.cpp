#include "gurobi_c++.h"
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <iostream> //input and output
#include <fstream>	//read file
#include <sstream>
#include <algorithm>
#include <vector>
#include <climits>
#include <set>
#include <string>
#include <ctime>	//to measure CPU time
#include <chrono> //to measure run time
#include <time.h>

#include <iomanip>

#include "readData.h"

#include <omp.h>

#include "MVBPMP_common_functions_LR_OMP.h"

using namespace std;
using namespace std::chrono;
using namespace std::string_literals;
// to generate vehicle data file name

int numV = 3; // number of vehicles; will be overwritten when reading vehicle file
// number of nodes in instance
// in main(), we also use n to present number of nodes
// since it's easy to read variables with short subscript, like x[n][n][numV]
int NUM_NODES;

// when runing on HPC, please change the required cores in sbatch file accordingly
// for example, if NUM_THREADS_VEH = 9, and 3 vehicles in the instance, then the cores
// should be set up at least 9*3=27
int NUM_THREADS_VEH = 28; // number of threads for solving each vehicle's BPMP model

double INITIAL_U_COEFFICIENT;

double LR_gap_tolerance = 0.05;
double LR_complementarity_tolerance = 0.0001;
double LR_min_lamda = 0.0001;
double LR_miu_tolerance = 0.00000001;

// feed the LR solution to Gurobi and find the optimal solution
bool DO_STAGE_TWO = false;

bool PRINT_x_WHEN_INTEGER_SOL = false;
bool PRINT_y_WHEN_INTEGER_SOL = false;

// according to p174 on "Integer Programming" by Laurence A. Wolsey 1st edition
bool USE_LR_MULTIPLIER_TYPE_B = false;
bool USE_LR_MULTIPLIER_TYPE_C = true;

bool PRINT_VAR_VALUR = false;
bool PRINT_CONFLICT_PICKUP = false;
bool PRINT_LRmultiplierTypeB_update_process = false;

bool ADD_PREPROCESS = true;

bool ADD_POTENTIAL_ARCS = true;
double percentageOfPotArcs = 0;
bool RUN_IN_PARALLEL_OMP = true; // it will be overwritten by passed arguments
int TIME_LIMIT = 3600;					 // it will be overwritten by passed arguments

double bigM = 10000000;

// solution from solving MVBPMP in Gurobi
double ***solx_GRB = NULL;
double ***soly_GRB = NULL;
double ****solu_GRB = NULL;
double **sols_GRB = NULL;
double LBinGRB = -bigM;
double UBinGRB = bigM;
double LBzero = 0.0000000001;
double LBinLR = -bigM;
double UBinLR = bigM;

void updateLRmultiplierTypeC(double ***, double, double *, double **, double, double);
void printVar(double ***, double ***, double ****, int *);
void updateLRmultiplierTypeB(double ***, double **, double, int, int);
void reportTime(clock_t, auto);
// void storeBestLB(double ***, double ***, double ****, double **, double ***, double ***,
// 								 double ****, double **,int);
void storeBestLB(double ***, double ***, double ***, double ***, int);
// double updateUB(double);
// void printBestLB(double);
// void printBestUB(double);

// string
// itos(int i)
// {
// 	stringstream s;
// 	s << i;
// 	return s.str();
// }

struct Cargo
{
	int origin;
	int end;
	double potentialProfit;
};

// Custom comparison function to sort by profit in descending order
bool compareCargoByProfit(const Cargo &a, const Cargo &b)
{
	return a.potentialProfit > b.potentialProfit; // Sorts in descending order of age
}

class printIntSol : public GRBCallback
{
public:
	int n;
	int numV;
	GRBVar ***xv;
	GRBVar ***yv;
	GRBVar ****uv;
	GRBVar **sv;

	printIntSol(GRBVar ***xvars, GRBVar ***yvars, GRBVar ****uvars, GRBVar **svars, int nvar, int numVvar)
	{
		xv = xvars;
		yv = yvars;
		uv = uvars;
		sv = svars;
		n = nvar;
		numV = numVvar;
	}

protected:
	void
	callback()
	{
		try
		{
			if (where == GRB_CB_MIPSOL)
			{
				// Found an integer feasible solution
				double objTemp = getDoubleInfo(GRB_CB_MIPSOL_OBJ);

				if (objTemp > LBinGRB)
					LBinGRB = objTemp;

				printf("--------> GRB: OBJ = %lf.\n", objTemp);

				UBinGRB = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);
				printf("--------> GRB:  UB = %lf.\n", UBinGRB);

				int i, j, q;
				// double ***x = NULL;
				// double ***y = NULL;

				// x = new double **[n];
				// y = new double **[n];

				// for (i = 0; i < n; i++)
				// {
				// 	x[i] = new double *[n];
				// 	y[i] = new double *[n];
				// 	for (j = 0; j < n; j++)
				// 	{
				// 		x[i][j] = new double[numV];
				// 		y[i][j] = new double[numV];
				// 	}
				// }

				for (i = 0; i < n; i++)
				{
					sols_GRB[i] = getSolution(sv[i], numV);
					for (j = 0; j < n; j++)
					{
						solx_GRB[i][j] = getSolution(xv[i][j], numV);
						soly_GRB[i][j] = getSolution(yv[i][j], numV);
						for (int k = 0; k < n; k++)
							solu_GRB[i][j][k] = getSolution(uv[i][j][k], numV);

						for (q = 0; q < numV; q++)
						{
							if (PRINT_x_WHEN_INTEGER_SOL)
								if (solx_GRB[i][j][q] > 0.5)
									printf("x: %3d ->%3d (v%d)\n", i + 1, j + 1,
												 q + 1);
							if (PRINT_y_WHEN_INTEGER_SOL)
								if (soly_GRB[i][j][q] > 0.5)
									printf("y: %3d ->%3d (v%d)\n", i + 1, j + 1,
												 q + 1);
						}
					}
				}

				// for (i = 0; i < n; i++)
				// {
				// 	for (j = 0; j < n; j++)
				// 	{
				// 		delete[] x[i][j];
				// 		delete[] y[i][j];
				// 	}
				// 	delete[] x[i];
				// 	delete[] y[i];
				// }
				// delete[] x;
				// delete[] y;
			}
		}
		catch (GRBException e)
		{
			cout << "Error number: " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...)
		{
			cout << "Error during callback" << endl;
		}
	}
};

//====================================================================================
// void storeBestLB(double ***solx_d, double ***soly_d, double ****solu_d,
// 								 double **sols_d, double ***solx_best, double ***soly_best,
// 								 double ****solu_best, double **sols_best, int n)
// {
// 	int i, j, k, q;
// 	for (i = 0; i < n; i++)
// 		for (q = 0; q < numV; q++)
// 		{
// 			sols_best[i][q] = sols_d[i][q];
// 			for (j = 0; j < n; j++)
// 			{
// 				solx_best[i][j][q] = solx_d[i][j][q];
// 				soly_best[i][j][q] = soly_d[i][j][q];

// 				for (k = 0; k < n; k++)
// 					solu_best[i][j][k][q] = solu_d[i][j][k][q];
// 			}
// 		}
// }

void storeBestLB(double ***solx_d, double ***soly_d, double ***solx_best, double ***soly_best, int n)
{
	int i, j, k, q;
	for (i = 0; i < n; i++)
		for (q = 0; q < numV; q++)
		{
			// sols_best[i][q] = sols_d[i][q];
			for (j = 0; j < n; j++)
			{
				solx_best[i][j][q] = solx_d[i][j][q];
				soly_best[i][j][q] = soly_d[i][j][q];

				// for (k = 0; k < n; k++)
				// 	solu_best[i][j][k][q] = solu_d[i][j][k][q];
			}
		}
}

void reportTime(clock_t begin, auto beginWallClock)
{
	// clock() gives cpu time on Linux, and wall time on Windows.
	// link: https://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows
	clock_t end = clock();
	double computerTime = (double)(end - begin);
	double second = computerTime / CLOCKS_PER_SEC;
	printf("CPU time (on Linux): %lf computer time,  %lf seconds\n",
				 computerTime, second);

	// wallclock time
	auto endWallClock = high_resolution_clock::now();
	auto elapsedWallClock = duration_cast<std::chrono::nanoseconds>(endWallClock - beginWallClock);
	printf("Wall clock time: %.3f seconds.\n", elapsedWallClock.count() * 1e-9);
}

void updateLRmultiplierTypeB(double ***soly, double **LR_u, double LR_rou,
														 int countItr, int numArcs)
{
	int i, j, q;
	int n = NUM_NODES;
	double slack[n][n];
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		{
			double sum = 0;
			for (q = 0; q < numV; q++)
				sum += soly[i][j][q];
			slack[i][j] = 1 - sum;
		}
	double LR_miu = 1.0 / ((double)10 * numArcs * (countItr + 1));

	// for (i = 0; i < countItr; i++)
	// LR_miu *= LR_rou;

	cout << "updated LR_miu = " << LR_miu << endl;

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (PRINT_LRmultiplierTypeB_update_process)
				printf("old LR_u_%d_%d = %lf  ", i, j, LR_u[i][j]);

			LR_u[i][j] = LR_u[i][j] - LR_miu * slack[i][j];

			if (PRINT_LRmultiplierTypeB_update_process)
			{
				printf("slack_%d_%d = %lf  ", i, j, slack[i][j]);
				printf("after updating: LR_u_%d_%d = %lf  ", i, j, LR_u[i][j]);
			}

			if (LR_u[i][j] < 0)
				LR_u[i][j] = 0;
		}
		if (PRINT_LRmultiplierTypeB_update_process)
		{
			cout << endl;
		}
	}
}

void printVar(double ***solx, double ***soly, double ****solu, int *origin)
{
	int i, j, k, q;
	int n = NUM_NODES;
	// cout << "SOLUTION:" << endl;

	for (q = 0; q < numV; q++)
	{
		cout
				<< "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
				<< endl;
		printf("Vehicle %d starting at node %d:\n", q + 1, origin[q] + 1);

		printf("x=1:\n");
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				if (solx[i][j][q] > 0.99 && solx[i][j][q] < 1.01)
					printf("(%d,%d) ", i + 1, j + 1);
		printf("\n");

		printf("y=1:\n");
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				if (soly[i][j][q] > 0.99 && soly[i][j][q] < 1.01)
					printf("(%d,%d) ", i + 1, j + 1);

		printf("\n");
		cout << "u>0.000001:" << endl;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				for (k = 0; k < n; k++)
					if (solu[i][j][k][q] > 0.000001)
						printf("(%d,%d,%d) ", i + 1, j + 1, k + 1);
		printf("\n");
	}
	cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
			 << endl;
}

void updateLRmultiplierTypeC(double ***soly, double LR_lamda, double *LR_miu,
														 double **LR_u, double LB, double totalProfit)
{
	int i, j, q;
	// update LR multiplier LR_u
	double norm = 0;
	int n = NUM_NODES;
	double slack[n][n];
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		{
			double sum = 0;
			for (q = 0; q < numV; q++)
				sum += soly[i][j][q];
			slack[i][j] = 1 - sum;
			// printf ("slack=%lf %d %d\n", slack[i][j], i, j);
			norm += slack[i][j] * slack[i][j];
		}
	cout << "norm = " << norm << endl;

	*LR_miu = LR_lamda * (totalProfit - LB) / norm;
	cout << "updated LR_miu = " << *LR_miu << endl;

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		{
			LR_u[i][j] = LR_u[i][j] - (*LR_miu) * slack[i][j];
			if (LR_u[i][j] < 0)
				LR_u[i][j] = 0;
		}
}

int main(int argc, char *argv[])
{
	//********* only read graph info and vehicle info by arguments
	//********* doesn't go through the graph info in the data folder
	if (argc != 6)
	{
		cout
				<< "Usage: ./mvbpmp_openmp.x nodesDataNameAndPath numberOfVehicles(int) vehicleDataNameAndPath TIME_LIMIT(int) useOMP(or noOMP)"
				<< endl;
		cout << "useOMP means MVBPMP will be solved by Gurobi in Thread 1 in paralle with Lagrangian relaxation in Thread 0" << endl;
		cout << "noOMP means program will only run Lagrangian relaxation in Thread 0." << endl;
		return 1;
	}

	if (argc == 6)
	{
		cout << "Nodes Data: " << argv[1] << endl;
		cout << "Number of Vehicles: " << argv[2] << endl;
		cout << "Vehicles Data: " << argv[3] << endl;
		numV = stoi(argv[2]);
		TIME_LIMIT = stoi(argv[4]);

		string argument5 = argv[5];
		if (argument5 == "useOMP")
			RUN_IN_PARALLEL_OMP = true;
		else if (argument5 == "noOMP")
			RUN_IN_PARALLEL_OMP = false;
		else
		{
			cout << "ERROR: the 5th argument should be string useOMP or noOMP" << endl;
			exit(1);
		}
	}

	clock_t beginTime, endTimeOfLastIteration;

	auto beginWallClock = high_resolution_clock::now();
	auto endTimeOfLastIterationWallClock = high_resolution_clock::now();

	readData route;

	string filename = argv[1];
	printf("reading file %s ...\n", filename.c_str());
	route.readSingleFile(filename);

	// route.printStats();
	// route.printData();

	int n = route.numOfNode;
	NUM_NODES = route.numOfNode;

	double wt[n][n];
	double dis[n][n];
	for (int i = 1; i <= n; i++)
		for (int j = 1; j <= n; j++)
		{
			wt[i - 1][j - 1] = route.w[i][j];
			dis[i - 1][j - 1] = route.d[i][j];
		}
	double price = route.priceCharged;
	double cost = route.travelCost;
	double vw = route.vehicleWeight;
	double Q = route.totalCapacity;
	double disLimit = (double)route.DIS;

	//==============reading vehicle info==============
	string vehicleFileName = argv[3];

	printf("reading vehicle file: %s ...\n", vehicleFileName.c_str());
	// exit (1);
	route.readSingleVehicleFile(vehicleFileName, numV);

	// vehicle and origin always starts from 0
	int vehicle[numV];
	int origin[numV];

	printf("vehicles (start from 0):\n");
	for (int i = 0; i < numV; i++)
	{
		vehicle[i] = route.vehicle[i];
		printf("%d ", vehicle[i]);
	}
	printf("\norigins (start from 0):\n");
	for (int i = 0; i < numV; i++)
	{
		origin[i] = route.origin[i];
		printf("%d ", origin[i]);
	}
	printf("\n");

	if (numV < 2)
	{
		cout << "The input only has 1 vehicle." << endl;
		cout << "This is desgined for multiple vehicle problem. Quit running." << endl;
		exit(1);
	}

	for (int i = 1; i < numV; i++)
	{
		if (origin[i] != origin[i - 1])
		{
			cout << "Vehicles are supposed to have the same start nodes (node 1 (0 in code)). Quit running." << endl;
			exit(1);
		}
	}

	// int allVehiclesStartNode = -1;
	// if (startNodesSame)
	// 	allVehiclesStartNode = 0;

	//====== vehicles' depots are the same ======
	// vehicles' startNode might be different, so the startNode will
	// be defined for each vehicle in LR
	int endNode = n - 1;

	//============== pre-select some cargoes which can potentially make more profit ============

	// // by initialize int selectedPositiveProfitCargos[n][n] = {0}; we get some value in [0][x]=a strange number
	// // so we still initialize by each element;
	// // int selectedPositiveProfitCargos[n][n] = {0};
	// // int selectedPositiveProfitCargos[n][n];
	// int selectedPositiveProfitCargos[n][n];
	// for (int i = 0; i < n; i++)
	// 	for (int j = 0; j < n; j++)
	// 		selectedPositiveProfitCargos[i][j] = 0;

	// if (ADD_POTENTIAL_ARCS)
	// {
	// 	vector<Cargo> sortedCargos;

	// 	double potentialProfit[n][n];
	// 	int numPositivePotentialProfit = 0;

	// 	for (int i = 0; i < n; i++)
	// 		for (int j = 0; j < n; j++)
	// 		{
	// 			if ((i == j) || (i == (n - 1)) || (j == 0))
	// 				potentialProfit[i][j] = -1000000;
	// 			else
	// 				potentialProfit[i][j] = price * dis[i][j] * wt[i][j] - cost * dis[i][j] * wt[i][j] - cost * vw * dis[i][j];
	// 			// printf("profit(%d, %d) = %lf  dis:%lf  wt: %lf\n", i, j, potentialProfit[i][j], dis[i][j], wt[i][j]);

	// 			if (potentialProfit[i][j] > 0)
	// 				numPositivePotentialProfit++;

	// 			// struct Cargo newCargo = Cargo(i, j, potentialProfit[i][j]);
	// 			struct Cargo newCargo = {i, j, potentialProfit[i][j]};

	// 			sortedCargos.push_back(newCargo);
	// 		}

	// 	// Sort the vector using the custom comparison function
	// 	sort(sortedCargos.begin(), sortedCargos.end(), compareCargoByProfit);

	// 	// Print sorted order
	// 	// cout << "\nSorted by profit (descending):" << endl;
	// 	// for (const Cargo &c : sortedCargos)
	// 	// {
	// 	// 	// cout << c.origin << ", " << c.end << " " << c.potentialProfit << endl;
	// 	// 	int i = c.origin;
	// 	// 	int j = c.end;
	// 	// 	printf("profit(%d, %d) = %lf  dis:%lf  wt: %lf\n",
	// 	// 				 i, j, potentialProfit[i][j], dis[i][j], wt[i][j]);
	// 	// }

	// 	cout << "percentageOfPotArcs" << endl;
	// 	cout << "===> number of cargos selected from best potential profit cargos: " << ceil(numPositivePotentialProfit * 0.1) << endl;
	// 	for (int i = 0; i < ceil(numPositivePotentialProfit * percentageOfPotArcs); i++)
	// 	{
	// 		selectedPositiveProfitCargos[sortedCargos[i].origin][sortedCargos[i].end] = 1;
	// 	}

	// 	// Print sorted order
	// 	// cout << "\nselected cargos:" << endl;
	// 	// for (int i = 0; i < n; i++)
	// 	// 	for (int j = 0; j < n; j++)
	// 	// 		if (selectedPositiveProfitCargos[i][j] > 0)
	// 	// 		{
	// 	// 			printf("profit(%d, %d) = %lf  dis:%lf  wt: %lf\n", i, j, potentialProfit[i][j], dis[i][j], wt[i][j]);
	// 	// 		}
	// }

	//============== add preprocess ==============
	//===> check if visiting some arcs over distance limit

	vector<vector<vector<int>>>
			allVehiclesInaccNeighbors;
	vector<vector<vector<int>>> allVehiclesNodeNeighbors;

	if (ADD_PREPROCESS)
	{
		for (int k = 0; k < numV; k++)
		{
			vector<vector<int>> nodeNeighbors;
			vector<vector<int>> nodeInaccNeighbors;
			int startNode = origin[k];

			// initialize nodeNeighbors and nodeInaccNeighbors
			for (int i = 0; i < n; i++)
			{
				vector<int> neighbors;
				if (i != endNode)
					for (int j = 0; j < n; j++)
					{
						if (j != startNode)
							neighbors.push_back(j);
					}
				nodeNeighbors.push_back(neighbors);

				vector<int> neighbors2;
				nodeInaccNeighbors.push_back(neighbors2);
			}

			// find the arcs violate distance limit
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
				{
					if (i != startNode && i != endNode && j != endNode && j != startNode && i != j)
					{
						double distTemp = dis[startNode][i] + dis[i][j] + dis[j][endNode];
						if (distTemp > disLimit)
						{
							// update nodeInaccNeighbors
							nodeInaccNeighbors[i].push_back(j);

							// update nodeNeighbors, remove j from neighbors
							auto it = find(nodeNeighbors[i].begin(), nodeNeighbors[i].end(), j);
							if (it != nodeNeighbors[i].end())
								nodeNeighbors[i].erase(it);
						}
					}
				}
			allVehiclesInaccNeighbors.push_back(nodeInaccNeighbors);
			allVehiclesNodeNeighbors.push_back(nodeNeighbors);
		}

		// cout << "allVehiclesInaccNeighbors" << endl;
		// for (auto &oneVeh : allVehiclesInaccNeighbors)
		// {
		// 	cout << "vehicle " << endl;
		// 	for (auto &oneNodeNb : oneVeh)
		// 	{
		// 		for (auto &e : oneNodeNb)
		// 			cout << e << " ";
		// 		cout << endl;
		// 	}
		// }

		// cout << "allVehiclesNodeNeighbors" << endl;
		// for (auto &oneVeh : allVehiclesNodeNeighbors)
		// {
		// 	cout << "vehicle " << endl;
		// 	for (auto &oneNodeNb : oneVeh)
		// 	{
		// 		for (auto &e : oneNodeNb)
		// 			cout << e << " ";
		// 		cout << endl;
		// 	}
		// }
	}

	double UB = bigM;
	double LB = -bigM;

	// double ***solx_GRB = new double **[n]; // solution x from solving MVBPMP in Gurobi
	// double ***soly_GRB = new double **[n];
	// double ****solu_GRB = new double ***[n];
	// double **sols_GRB = new double *[n];

	solx_GRB = new double **[n]; // solution x from solving MVBPMP in Gurobi
	soly_GRB = new double **[n];
	solu_GRB = new double ***[n];
	sols_GRB = new double *[n];

	for (int i = 0; i < n; i++)
	{
		solx_GRB[i] = new double *[n];
		soly_GRB[i] = new double *[n];
		solu_GRB[i] = new double **[n];
		sols_GRB[i] = new double[numV];
		for (int j = 0; j < n; j++)
		{
			solx_GRB[i][j] = new double[numV];
			soly_GRB[i][j] = new double[numV];
			solu_GRB[i][j] = new double *[n];

			for (int k = 0; k < n; k++)
				solu_GRB[i][j][k] = new double[numV];
		}
	}

	double ***solx_best = new double **[n];
	double ***soly_best = new double **[n];
	double ****solu_best = new double ***[n];
	double **sols_best = new double *[n];

	for (int i = 0; i < n; i++)
	{
		solx_best[i] = new double *[n];
		soly_best[i] = new double *[n];
		solu_best[i] = new double **[n];
		sols_best[i] = new double[numV];
		for (int j = 0; j < n; j++)
		{
			solx_best[i][j] = new double[numV];
			soly_best[i][j] = new double[numV];
			solu_best[i][j] = new double *[n];

			for (int k = 0; k < n; k++)
				solu_best[i][j][k] = new double[numV];
		}
	}

	//============== start Lagrangian Relaxation and Gurobi in parallel ==============

	// figure out how many can run at once based on available cores
	const int max_threads = omp_get_max_threads();
	cout << "max_threads: " << max_threads << endl;
	int requiredNumCores = NUM_THREADS_VEH * 2 + 4;
	cout << "required number of cores (requiredNumCores): " << requiredNumCores << endl;

	if (max_threads < requiredNumCores)
	{
		cout << "Quit running: there are not enough threads available."
				 << endl;
		exit(1);
	}

	omp_set_nested(1); // Enable nested parallelism

	int numThreadsNeeded = 2;
	if (!RUN_IN_PARALLEL_OMP)
		numThreadsNeeded = 1;

// #pragma omp parallel num_threads(2) default(shared) // Outer parallel region with 2 threads
#pragma omp parallel num_threads(numThreadsNeeded) // Outer parallel region with 2 threads
	{
		int outer_thread_id = omp_get_thread_num();
		// cout << "outer_thread_id = " << outer_thread_id << endl;

		//============== start Lagrangian Relaxation ==============
		if (outer_thread_id == 0)
		{

			// double UB = bigM;
			// double LB = -bigM;

			double LR_rou = 0.4; // for LR multiplier type b

			double LR_lamda = 2;
			int LR_maxNumNoImprovement = 3;
			int LR_numNoImprovement = 0;

			bool findOptimalSolution = false;

			// double LR_u[n][n];	      //the lagrangian multiplier
			double **LR_u = new double *[n];
			for (int i = 0; i < n; i++)
				LR_u[i] = new double[n];

			// initialize LR_u
			double numArcs = n * (n - 1);
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					LR_u[i][j] = 1 / numArcs;

			double LR_lowestUpperBound = bigM;

			// store the solution of MVBPMP_LR in each iteration
			// double solx[n][n][numV];
			// double soly[n][n][numV];
			// double solu[n][n][n][numV];

			// solution x for LR dual with numV best solution from Gurobi
			double ***solx_numV_best_d = new double **[n];
			// EXAMPLE:
			//  solution1 s1, profit(s1)=1.783456, profit(s2)=1.783457, profit(s3)=1.782222
			//  according to "same profit" definition later (diff<=0.000001 (1E-6))
			//  s1 and s2 has same profit, so they are two different solutions
			//(routes or selected cargoes are different) with the same profit.
			//  so solx_opt_d stores the solution of s1, s2, and another s1,
			// since solx_opt_d only stores the solution of the best profit, so s3 is replacew ith s1
			double ***solx_opt_d = new double **[n];

			// when vehicles start from the same node, soly_numV_best_d stores the numV best solution for
			// BPMP from Gurobi
			double ***soly_numV_best_d = new double **[n];
			// then we select the solution with optimal profit and store in soly_opt_d
			double ***soly_opt_d = new double **[n];

			//!!!!!! _numV_best_ means it is numV best solution from Gurobi
			// since we only need it when calculation LB, and only x,y and theta are needed in LB part
			//  so for other vars, like u,s, solu_d and sols_d are actually solutions of numV best solutions

			double ****solu_d = new double ***[n];

			double soltheta_d[n][n][numV];

			// double sols_d[n][numV];
			double **sols_d = new double *[n];
			double profit_d[numV];

			for (int i = 0; i < n; i++)
			{
				solx_numV_best_d[i] = new double *[n];
				solx_opt_d[i] = new double *[n];
				soly_numV_best_d[i] = new double *[n];
				soly_opt_d[i] = new double *[n];
				solu_d[i] = new double **[n];
				sols_d[i] = new double[numV];

				// solx_best[i] = new double *[n];
				// soly_best[i] = new double *[n];
				// solu_best[i] = new double **[n];
				// sols_best[i] = new double[numV];
				for (int j = 0; j < n; j++)
				{
					solx_numV_best_d[i][j] = new double[numV];
					solx_opt_d[i][j] = new double[numV];
					soly_numV_best_d[i][j] = new double[numV];
					soly_opt_d[i][j] = new double[numV];
					solu_d[i][j] = new double *[n];

					// solx_best[i][j] = new double[numV];
					// soly_best[i][j] = new double[numV];
					// solu_best[i][j] = new double *[n];

					for (int k = 0; k < n; k++)
					{
						solu_d[i][j][k] = new double[numV];
						// solu_best[i][j][k] = new double[numV];
					}
				}
			}

			int maxItr = 1;
			int countItr = 0;

			cout << "===> LR_gap_tolerance = " << LR_gap_tolerance << endl;

			beginTime = clock();

			auto timelimitnano = std::chrono::nanoseconds(TIME_LIMIT);

			endTimeOfLastIterationWallClock = high_resolution_clock::now();
			while (((UB - LB) / LB > LR_gap_tolerance) || ((UB - LB) / LB < -LR_gap_tolerance))
			{
				auto endWallClock = high_resolution_clock::now();
				auto elapsedWallClock = duration_cast<std::chrono::nanoseconds>(endWallClock - beginWallClock);

				cout << "elapsedWallClock:" << elapsedWallClock.count() * 1e-9 << endl;
				// cout << "timelimitnano: " << timelimitnano.count() << endl;
				if (elapsedWallClock.count() * 1e-9 > timelimitnano.count())
					break;

				countItr++;
				// stop after a certain iterations
				// if (countItr > maxItr)
				// 	break;

				cout << "" << endl;
				cout << "============= Itr " << countItr << " ==============" << endl;
				cout << "LB = " << LB << endl;
				cout << "UB = " << UB << endl;

				//=======================> SOLVING BPMP <======================
				{

					// int nthreads = omp_get_num_threads();
					// cout << "NUM_THREADS=" << nthreads << endl;

					// start calling Gurobi
					int i, j, k;
					// int startNode = origin[vehicleIndex];
					int startNode = 0;

					GRBEnv *env = NULL;
					GRBVar **x = NULL;
					GRBVar **y = NULL;
					GRBVar s[n];
					GRBVar u[n][n][n];
					GRBVar theta[n][n];

					x = new GRBVar *[n];
					y = new GRBVar *[n];
					for (i = 0; i < n; i++)
					{
						x[i] = new GRBVar[n];
						y[i] = new GRBVar[n];
					}

					try
					{
						env = new GRBEnv();
						GRBModel model = GRBModel(*env);

						// Create binary decision variables
						for (i = 0; i < n; i++)
						{
							s[i] = model.addVar(0.0, n, 0.0, GRB_CONTINUOUS, "s_" + itos(i));
							// s[i] = model.addVar(0.0, n, 0.0, GRB_INTEGER, "s_" + itos(i));

							for (j = 0; j < n; j++)
							{
								x[i][j] = model.addVar(0.0, 1.0, 0, GRB_BINARY,
																			 "x_" + itos(i) + "_" + itos(j));
								y[i][j] = model.addVar(0.0, 1.0, 0, GRB_BINARY,
																			 "y_" + itos(i) + "_" + itos(j));
								theta[i][j] = model.addVar(
										0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
										"theta_" + itos(i) + "_" + itos(j));

								for (k = 0; k < n; k++)
								{
									string s = "u_" + itos(i) + "_" + itos(j) + "_" + itos(k);
									u[i][j][k] = model.addVar(0.0, GRB_INFINITY, 0.0,
																						GRB_CONTINUOUS, s);
								}
							}
						}

						// force some x, y variables to be zeros
						for (i = 0; i < n; i++)
						{
							x[i][i].set(GRB_DoubleAttr_UB, 0);
							y[i][i].set(GRB_DoubleAttr_UB, 0);
							x[i][startNode].set(GRB_DoubleAttr_UB, 0);
							y[i][startNode].set(GRB_DoubleAttr_UB, 0);
							x[n - 1][i].set(GRB_DoubleAttr_UB, 0);
							y[n - 1][i].set(GRB_DoubleAttr_UB, 0);
						}

						//====================AVOID DUPLICATE SOLUTIONS=========================
						{
							// when w(i,j)=0, force y(i,j)=0 to avoid y(i,j)=1, which is another feasible solution, but actually the same route and cargo selection
							for (i = 0; i < n; i++)
								for (j = 0; j < n; j++)
								{
									if (wt[i][j] < 0.000001)
									{
										// cout << "wt[" << i + 1 << "][" << j + 1 << "]<0.000001" << endl;
										// cout << "LR_u[" << i + 1 << "][" << j + 1 << "]=" << LR_u[i][j] << endl;
										y[i][j].set(GRB_DoubleAttr_UB, 0);
									}
								}

							// add constraint to force s[i] be a specific value, not vary and leads to many solutions
							s[0].set(GRB_DoubleAttr_LB, 1);
							s[0].set(GRB_DoubleAttr_UB, 1);
							GRBLinExpr exprAvoidDuplicateSolutions = 0.0;
							for (i = 0; i < n; i++)
								for (j = 0; j < n; j++)
									exprAvoidDuplicateSolutions += x[i][j];

							exprAvoidDuplicateSolutions += 1;
							exprAvoidDuplicateSolutions -= s[n - 1];
							model.addConstr(exprAvoidDuplicateSolutions >= 0, "avoidDuplicate_" + itos(i));

							for (i = 1; i < n - 1; i++) // s[n-1] <= sum x(i,j) as defined above
							{
								GRBLinExpr exprSzeroWhenNotVisited = 0.0;
								for (j = 0; j < n; j++)
									exprSzeroWhenNotVisited += x[i][j];

								exprSzeroWhenNotVisited *= n;
								// exprSzeroWhenNotVisited += 0.000001;
								exprSzeroWhenNotVisited -= s[i];
								model.addConstr(exprSzeroWhenNotVisited >= 0, "avoidSzeroWhenNotVisited_" + itos(i));
							}

						} //===END OF AVOIDING DUPLICATE SOLUTIONS===

						// force some triples variables to be zeros
						for (i = 0; i < n; i++)
							for (j = 0; j < n; j++)
							{
								u[n - 1][i][j].set(GRB_DoubleAttr_UB, 0);
								u[i][startNode][j].set(GRB_DoubleAttr_UB, 0);
								u[i][j][startNode].set(GRB_DoubleAttr_UB, 0);
								u[i][j][n - 1].set(GRB_DoubleAttr_UB, 0);
							}

						if (ADD_PREPROCESS)
						{
							for (int i = 0; i < n; i++)
							{
								vector<int> nbsTemp = allVehiclesInaccNeighbors[startNode][i];
								for (j = 0; j < nbsTemp.size(); j++)
								{
									x[i][nbsTemp[j]].set(GRB_DoubleAttr_UB, 0);
									y[i][nbsTemp[j]].set(GRB_DoubleAttr_UB, 0);
									for (k = 0; k < n; k++)
									{
										u[i][nbsTemp[j]][k].set(GRB_DoubleAttr_UB, 0);
										u[i][k][nbsTemp[j]].set(GRB_DoubleAttr_UB, 0);
										u[k][nbsTemp[j]][i].set(GRB_DoubleAttr_UB, 0);
									}
								}
							}

							// int numViolatedTriples = 0;
							for (int i = 0; i < n; i++)
								if (i != startNode && i != endNode)
									for (auto &k : allVehiclesNodeNeighbors[startNode][i])
										for (auto &j : allVehiclesNodeNeighbors[startNode][k])
											if (i != j)
												if (dis[startNode][i] + dis[i][k] + dis[k][j] + dis[j][endNode] > disLimit)
												{
													u[i][j][k].set(GRB_DoubleAttr_UB, 0);
													// numViolatedTriples++;
												}
							// cout << "numViolatedTriples=" << numViolatedTriples << endl;
						}

						{
							double yUpperBoundCopy[n][n];
							for (i = 0; i < n; i++)
								for (j = 0; j < n; j++)
									yUpperBoundCopy[i][j] = 1;

							int countTemp = 0;
							for (i = 0; i < n; i++)
								for (j = 0; j < n; j++)
									if (price * dis[i][j] * wt[i][j] - LR_u[i][j] <= 0)
									{
										y[i][j].set(GRB_DoubleAttr_UB, 0);
										yUpperBoundCopy[i][j] = 0;
										countTemp++;
									}
							cout << "===> number of y that is set to be zero because Lagrangian multiplier is non-positive: " << countTemp << endl;

							// if a node has no cargoes in or out, then corresponding x var is zero
							cout << "===> nodes  that won't be visited:" << endl;
							countTemp = 0;
							for (i = 0; i < n; i++)
							{
								double sumYub = 0; // ub means upper bound
								for (j = 0; j < n; j++)
								{
									double ubTemp1 = yUpperBoundCopy[i][j];
									double ubTemp2 = yUpperBoundCopy[j][i];
									sumYub += ubTemp1 + ubTemp2;

									if (ubTemp1 > 0.99 || ubTemp2 > 0.99)
										break;
								}
								if (sumYub < 0.01)
								{
									cout << i + 1 << " ";
									for (j = 0; j < n; j++)
									{
										x[i][j].set(GRB_DoubleAttr_UB, 0);
										x[j][i].set(GRB_DoubleAttr_UB, 0);
										countTemp += 2;
									}
								}
							}
							cout << endl;
							cout << "===>number of x set to zero: " << countTemp << endl;
						}

						//==============generate constraints in Gurobi================
						// vehicle goes out of vehicle's origin
						GRBLinExpr expr1 = 0.0;
						for (i = 0; i < n; i++)
							expr1 += x[startNode][i];
						model.addConstr(expr1 == 1, "origin");

						// vehicle goes back to node n
						GRBLinExpr expr2 = 0.0;
						for (i = 0; i < n; i++)
							expr2 += x[i][n - 1];
						model.addConstr(expr2 == 1, "destination");

						// flow conservation (the last node is the destination)
						for (k = 0; k < n - 1; k++)
						{
							if (k != startNode)
							{
								GRBLinExpr expr = 0;
								for (i = 0; i < n; i++)
									expr += x[i][k];
								for (j = 0; j < n; j++)
									expr -= x[k][j];
								model.addConstr(expr == 0,
																"flow_conservation_" + itos(k));
							}
						}

						// distance
						GRBLinExpr expr3 = 0.0;
						for (i = 0; i < n; i++)
							for (j = 0; j < n; j++)
								expr3 += dis[i][j] * x[i][j];
						model.addConstr(expr3 <= route.DIS, "distance");

						// node degree less than 1
						for (k = 0; k < n; k++)
						{
							GRBLinExpr expr = 0.0;
							for (i = 0; i < n - 1; i++)
								expr += x[i][k];
							model.addConstr(expr <= 1, "indegree_" + itos(k));
						}

						// subtour elimination
						for (i = 0; i < n; i++)
							for (j = 0; j < n; j++)
							{
								GRBLinExpr expr = 0.0;
								expr += s[i] - s[j] + (n - 1) * x[i][j] + (n - 3) * x[j][i];
								model.addConstr(expr <= n - 2,
																"s_" + itos(i) + "_" + itos(j));
							}

						// arc flow
						for (i = 0; i < n; i++)
							for (j = 0; j < n; j++)
							{
								GRBLinExpr expr = 0.0;
								expr += wt[i][j] * y[i][j] - theta[i][j];
								for (k = 0; k < n; k++)
									expr += u[i][k][j] + u[k][j][i] - u[i][j][k];
								model.addConstr(expr == 0,
																"flow_" + itos(i) + "_" + itos(j));
							}

						// arc flow upperbound
						for (i = 0; i < n; i++)
							for (j = 0; j < n; j++)
							{
								GRBLinExpr expr = 0.0;
								expr += theta[i][j] - Q * x[i][j];
								model.addConstr(expr <= 0,
																"flowBound_" + itos(i) + "_" + itos(j));
							}

						//==============set objective function in Gurobi================

						GRBLinExpr obj = 0.0;
						for (i = 0; i < n; i++)
							for (j = 0; j < n; j++)
							{
								obj += (price * dis[i][j] * wt[i][j] - LR_u[i][j]) * y[i][j];
								obj -= cost * dis[i][j] * theta[i][j];
								obj -= cost * vw * dis[i][j] * x[i][j];
							}
						model.setObjective(obj, GRB_MAXIMIZE);

						//==============set up time limit in Gurobi================
						auto endWallClock = high_resolution_clock::now();
						auto elapsedWallClock = duration_cast<std::chrono::nanoseconds>(endWallClock - beginWallClock);

						cout << "elapsedWallClock before starting Gurobi for BPMP:" << elapsedWallClock.count() * 1e-9 << endl;
						// cout << "timelimitnano: " << timelimitnano.count() << endl;
						if (elapsedWallClock.count() * 1e-9 > timelimitnano.count())
							break;

						double availableTime = timelimitnano.count() - elapsedWallClock.count() * 1e-9;
						cout << "availableTime = " << availableTime << endl;

						model.set(GRB_DoubleParam_TimeLimit, availableTime);

						// set up logfile
						// model.set("LogFile", "mvbpmp_" + itos(n) + "_" + itos(vehicleIndex) + ".log");
						// send the log to a file only
						// model.set(GRB_IntParam_LogToConsole, 0);
						// keep outputflag 1 to print out log
						// model.set(GRB_IntParam_OutputFlag, 1);
						model.set(GRB_IntParam_OutputFlag, 0);

						model.set(GRB_IntParam_PoolSearchMode, 2);

						model.set(GRB_IntParam_PoolSolutions, numV);

						// set the number of threads to required number of Threads
						model.set(GRB_IntParam_Threads, NUM_THREADS_VEH);

						//==============solve the model in Gurobi================
						model.optimize();

						// write model to file
						// model.write("BPMP_" + itos(n) + "_" + itos(vehicleIndex) + ".lp");
						// model.write("BPMP_" + itos(n) + ".lp");

						//==============Extract solution from Gurobi================
						// since we set up time limit for Gurobi, it might stop and return solution before getting optimal
						// but we only need the optimal solution to get UB
						int status = model.get(GRB_IntAttr_Status);
						if (status == GRB_OPTIMAL)
						{
							// Print number of solutions stored
							int nSolutions = model.get(GRB_IntAttr_SolCount);
							cout << "Number of solutions found: " << nSolutions << endl;

							if (nSolutions < numV)
							{
								cout << "nSolutions=" << nSolutions << endl;
								cout << "Do not get enough different solutions. Quit Running." << endl;
								exit(1);
							}

							double runtime = model.get(GRB_DoubleAttr_Runtime);
							cout << "LR-BPMPRuntime: " << runtime << " seconds" << endl;

							// Print objective values of solutions
							for (int e = 0; e < nSolutions; e++)
							{
								model.set(GRB_IntParam_SolutionNumber, e);
								cout << model.get(GRB_DoubleAttr_PoolObjVal) << " ";
								if (e % 10 == 9)
									cout << endl;
							}
							cout << endl;

							// print
							// double oldobj = -1;
							// double newobj = 0;
							for (int e = 0; e < nSolutions; e++)
							// for (int e = 0; e < 2; e++)
							{
								model.set(GRB_IntParam_SolutionNumber, e);
								profit_d[e] = model.get(GRB_DoubleAttr_PoolObjVal);

								// printf("=== %d (obj=%lf) ===\n", e + 1, profit_d[e]);
								// printf("------ print x=1 ------\n");
								// for (q = 0; q < numV; q++)
								for (i = 0; i < n - 1; i++)
									for (j = 0; j < n; j++)
									{
										double xval = x[i][j].get(GRB_DoubleAttr_Xn);
										if (xval > 0.99 && xval < 1.01)
											printf("%3d %3d  \n", i + 1, j + 1);
										solx_numV_best_d[i][j][e] = xval;
									}
								// printf("------ print y=1 ------\n");
								// for (q = 0; q < numV; q++)
								for (i = 0; i < n - 1; i++)
									for (j = 0; j < n; j++)
									{
										double yval = y[i][j].get(GRB_DoubleAttr_Xn);
										if (yval > 0.99 && yval < 1.01)
											printf("%3d %3d  (w=%.2lf)\n", i + 1, j + 1, wt[i][j]);
										soly_numV_best_d[i][j][e] = yval;
									}
								// printf("------ print theta>0.000001 ------\n");
								// for (q = 0; q < numV; q++)
								for (i = 0; i < n - 1; i++)
									for (j = 0; j < n; j++)
									{
										double thetaval = theta[i][j].get(
												GRB_DoubleAttr_Xn);
										// if (thetaval > 0.000001)
										// 	printf("%3d %3d  flow=%lf\n", i + 1, j + 1, thetaval);
										soltheta_d[i][j][e] = thetaval;
									}
								// printf("------ print u>0.000001 ------\n");
								// for (q = 0; q < numV; q++)
								{
									// int ogn = origin[q];
									for (i = 0; i < n - 1; i++)
										for (j = 0; j < n; j++)
											for (k = 0; k < n - 1; k++)
											// 	if (k != ogn)
											{
												double uval = u[i][j][k].get(GRB_DoubleAttr_Xn);
												// if (uval > 0.000001)
												// 	printf("%3d %3d %3d  u=%lf\n", i + 1, j + 1, k + 1, uval);
												solu_d[i][j][k][e] = uval;
											}
								}
								// printf("------ print s>0.000001 ------\n");
								for (i = 0; i < n; i++)
								{
									double sval = s[i].get(GRB_DoubleAttr_Xn);
									// if (sval > 0.000001)
									// 	printf("s%d=%.2lf  \n", i + 1, sval);
									sols_d[i][e] = sval;
								}
							}
						}
						else
						{
							cout << "WARNING: Gurobi stop before getting optimum because of time limit." << endl;
							cout << "Will not read solution." << endl;
							continue;
						}
					}
					catch (GRBException e)
					{
						cout << "Error number: " << e.getErrorCode() << endl;
						cout << e.getMessage() << endl;
					}
					catch (...)
					{
						cout << "Error during optimization" << endl;
					}

					for (i = 0; i < n; i++)
					{
						delete[] x[i];
						delete[] y[i];
					}
					delete env;
				} // end of openmp parallel running using Gurobi

				// collect profit from all scenarios
				// each scenario represent one vehicle's schedule
				double totalProfit = 0;
				// for (int i = 0; i < numV; i++)
				// {
				// 	totalProfit += profit_d[i];
				// 	cout << "profit for vehicle " << i << ": " << profit_d[i] << endl;
				// }

				// When the starting nodes are the same, profit_d[1] is the profit of the 2nd best, etc.
				// Actually, all vehicles have the same profit since the start nodes are the same.
				// but there might be multiple best solutions, so we need to check if the 2nd and 3rd best
				// are the same with the 1st best profit.

				int numOptimalProfitSol = 1;

				for (int i = 1; i < numV; i++)
				{ // we assume than if the next best profit is within range below, it is the same as the best one
					// for example, if p1=1.999999, p2=1.999998, then by scaling up by a factor of 2500
					// p1=1.999999*2500=$4999.9975, p2=999998*2500=$4999.995, the diff is less than 1 cent
					// so we assume these two profits are the same
					double pftTolerance = 0.000001;
					if (profit_d[i] >= profit_d[0] - pftTolerance && profit_d[i] <= profit_d[0] + pftTolerance)
						numOptimalProfitSol++;
				}

				cout << "===> numOptimalProfitSol=" << numOptimalProfitSol << endl;

				//====== get optimal solution y from numV best solutions
				// if only numOptimalProfitSol solutions are optimal, then for the rest, use the best solution
				// for example, if 3 vehicles, 2 different solutions have optimal proit, the 3rd one is a little bit worse
				// then the 3rd one use the 1st vehicle's solution as it's solution
				for (int q = 0; q < numV; q++)
				{
					for (int i = 0; i < n; i++)
						for (int j = 0; j < n; j++)
						{
							if (q < numOptimalProfitSol)
							{
								solx_opt_d[i][j][q] = solx_numV_best_d[i][j][q];
								soly_opt_d[i][j][q] = soly_numV_best_d[i][j][q];
							}
							else
							{
								solx_opt_d[i][j][q] = solx_numV_best_d[i][j][0];
								soly_opt_d[i][j][q] = soly_numV_best_d[i][j][0];
							}
						}
				}
				totalProfit = numV * profit_d[0];

				cout << "===> sum of OBJ VALUE: " << totalProfit << endl;

				// add the lagrangian multiplier part sum_LR_u back to obj
				// since this is constant and not considered in scenarios
				// sum_LR_u_times_y is for printing here and future use when solution
				// is feasible for the original problem
				double sum_LR_u = 0;
				double sum_LR_u_times_y = 0;
				bool findConflictPickupInOptSol = false;
				bool findConflictPickupInNumVehBestSol = false;
				for (int i = 0; i < n; i++)
					for (int j = 0; j < n; j++)
					{
						sum_LR_u += LR_u[i][j];

						double sum = 0;

						// here sum is for calculation of UBinLR, so only the best profit solution is used
						sum = 0;
						for (int q = 0; q < numV; q++)
							sum += soly_opt_d[i][j][q];

						if (sum > 1.9)
						{
							findConflictPickupInOptSol = true;
							if (PRINT_CONFLICT_PICKUP)
								printf("===> conflict %d -> %d\n", i + 1, j + 1);
						}
						sum_LR_u_times_y += LR_u[i][j] * sum;
					}

				cout << "sum_LR_u = " << sum_LR_u << endl;
				cout << "sum_LR_u_times_y = " << sum_LR_u_times_y << endl;

				totalProfit += sum_LR_u;
				// cout << "totalProfit in obj including LR = " << totalProfit << endl;
				printf("totalProfit in obj including LR =%lf\n", totalProfit);

				double profitInLRdual = totalProfit;

				if (profitInLRdual < UBinLR)
					UBinLR = profitInLRdual;

				cout << "===> UBinLR = " << UBinLR << endl;

				if (UBinLR >= LR_lowestUpperBound)
					LR_numNoImprovement++;
				else
				{
					LR_numNoImprovement = 0;
					LR_lowestUpperBound = UBinLR;
					cout << "find lower UBinLR" << endl;
				}

				printf("===> LR_numNoImprovement = %d\n", LR_numNoImprovement);

				//===> update LR parameters for the next iteration
				if (LR_numNoImprovement >= LR_maxNumNoImprovement)
				{
					LR_lamda = 0.5 * LR_lamda;
					LR_numNoImprovement = 0;
				}

				cout << "===> LR_lamda = " << LR_lamda << endl;

				if (LR_lamda <= LR_min_lamda)
				{
					cout << "===> LR_lamda is less than LR_min_damda " << LR_min_lamda
							 << endl;
					cout << "end loop";
					break;
				}

				for (int i = 0; i < n; i++)
				{
					for (int j = 0; j < n; j++)
					{
						// check if there are pickup conflicts in numV best solutions
						double sum = 0;
						for (int q = 0; q < numV; q++)
							sum += soly_numV_best_d[i][j][q];

						if (sum > 1.9)
						{
							findConflictPickupInNumVehBestSol = true;
							break;
						}
					}
					if (findConflictPickupInNumVehBestSol)
						break;
				}

				//===> if there is no request conflict
				//===> use the current solution without LR item as LB
				// if (!findConflictPickupInOptSol)
				if (!findConflictPickupInNumVehBestSol)
				{
					cout << "===> There is no request conflict in numV best solutions :)" << endl;

					double profitForLB;

					if (numOptimalProfitSol == numV)
					{
						// check if the LR_item is close to zero
						// if yes, then this is the optimal solution
						double LR_item = sum_LR_u - sum_LR_u_times_y;
						cout << "===> the LR items added in obj = " << LR_item << endl;

						if ((LR_item >= 0 && LR_item <= LR_complementarity_tolerance) || (LR_item < 0 && LR_item >= -LR_complementarity_tolerance))
						{
							LB = profitInLRdual - LR_item;

							LBinLR = LB;

							findOptimalSolution = true;

							cout << "===> the abs(LR item) added in obj is less than "
									 << LR_complementarity_tolerance << endl;
							cout << "===> the optimal solution is found :)" << endl;
							cout << "===> end loop." << endl;

							reportTime(endTimeOfLastIteration,
												 endTimeOfLastIterationWallClock);

							cout << "**************** THE OPTIMAL SOLUTION ****************"
									 << endl;
							printVar(solx_opt_d, soly_opt_d, solu_d, origin);

							// no need to store the solution since the optimal solution is found
							// and there is no need to feed solution to MVBPMP model to find optimal solution
							break;
						}
						else
						{
							// if LR_item is still big, then we need to keep iteration
							cout << "===> Use this solution to calculate LB." << endl;

							// since no request conflict
							// so the solution for x, y, u are feasible for the original model
							// but - u_(i,j)*( 1-sum_q y(i,j,q) ) was added to obj
							// so we need to recalculte total profit
							profitForLB = profitInLRdual;
							profitForLB += sum_LR_u_times_y;
							profitForLB -= sum_LR_u;
						}
					}
					else
					{
						if (numOptimalProfitSol > numV)
						{
							cout << "ERROR: numOptimalProfitSol should be alwasy <= numV! Exit." << endl;
							exit(1);
						}

						profitForLB = 0;
						for (int q = 0; q < numV; q++)
							profitForLB += profit_d[q];

						for (int i = 0; i < n; i++)
							for (int j = 0; j < n; j++)
							{
								double sum = 0;

								// here sum is for calculation of UBinLR, so only the best profit solution is used
								sum = 0;
								for (int q = 0; q < numV; q++)
									sum += soly_numV_best_d[i][j][q];

								profitForLB += LR_u[i][j] * sum;
							}
					}

					if (profitForLB > LBinLR)
						LBinLR = profitForLB;

					if (profitForLB > LBinGRB && profitForLB > LB)
					{
						LB = profitForLB;
						cout << "===> During LR loop, find a better LB from LR dual = " << LB << endl;
						LBinLR = profitForLB;
						// store the solution as the best LB
						// storeBestLB(solx_numV_best_d, soly_numV_best_d, solu_d, sols_d, solx_best, soly_best, solu_best, sols_best, n);
						storeBestLB(solx_numV_best_d, soly_numV_best_d, solx_best, soly_best, n);
					}
					else if (LBinGRB > profitForLB && LBinGRB > LB)
					{
						LB = LBinGRB;
						cout << "===> During LR loop, find a better LB from LB_GRB =" << LBinGRB << endl;

						// storeBestLB(solx_GRB, soly_GRB, solu_GRB, sols_GRB, solx_best, soly_best,solu_best, sols_best, n);
						storeBestLB(solx_GRB, soly_GRB, solx_best, soly_best, n);
					}

					double LR_miu;

					if (USE_LR_MULTIPLIER_TYPE_C)
						updateLRmultiplierTypeC(soly_opt_d, LR_lamda, &LR_miu, LR_u, LB,
																		profitInLRdual);

					// if (USE_LR_MULTIPLIER_TYPE_B)
					// 	updateLRmultiplierTypeB(soly_opt_d, LR_u, LR_rou, countItr, numArcs);

					// double LR_miu;
					// if (USE_LR_MULTIPLIER_TYPE_C)
					// 	updateLRmultiplierTypeC(soly_opt_d, LR_lamda, &LR_miu, LR_u, LB,
					// 													profitInLRdual);

					reportTime(endTimeOfLastIteration, endTimeOfLastIterationWallClock);
					endTimeOfLastIteration = clock();
					endTimeOfLastIterationWallClock = high_resolution_clock::now();

					// if (UB - LB < 0)
					// {
					// 	cout << "UB = " << UB << endl;
					// 	cout << "LB = " << LB << endl;
					// 	cout
					// 			<< "UB is less than LB, which is an error. Stop running and check code."
					// 			<< endl;
					// 	exit(1);
					// }

					if (LB == 0)
						LB = LBzero;

					// // if ((UB - LB) / LB <= LR_gap_tolerance)
					// if (((UB - LB) / LB <= LR_gap_tolerance) && ((UB - LB) / LB >= -LR_gap_tolerance))
					// {
					// 	cout << "=== will stop loop because gap " << (UB - LB) / LB
					// 			 << " <= LR_gap_tolerance  and >= -LR_gap_tolerance " << LR_gap_tolerance << endl;

					// 	// here if(profitForLB > LB), it is already checked before
					// 	// and the new LB is stored
					// 	// if(profitForLB <= LB), then no need to store the sol?_d solution as best LB
					// 	// so to conclude, no need to storeBestLB here
					// 	// storeBestLB (solx_d, soly_d, solu_d, sols_d, solx_best, soly_best, solu_best, sols_best);
					// }
					// // this condition is added after adding OpenMP feature.
					// // if no added, the code will calculate MVBPMP for new LB.
					// // better LB might be found, but not necessary since gap is already within gap tolerance
					// else
					// 	continue;

					if (LR_miu <= LR_miu_tolerance)
					{
						cout << "=== will stop loop because LR_miu < LR_miu_tolerance "
								 << LR_miu_tolerance << endl;

						// here if(profitForLB > LB), it is already checked before
						// and the new LB is stored
						// if(profitForLB <= LB), then no need to store the sol?_d solution as best LB
						// so to conclude, no need to storeBestLB here
						// storeBestLB (solx_d, soly_d, solu_d, sols_d, solx_best, soly_best, solu_best, sols_best);

						break;
					}
					// else
					// 	continue; // jump to the next iteration of while loop

					continue;
				}
				else
				{
					// if there is request conflict, then call gurobi to find LB

					//========================== calculate LB ========================//

					// 1. If we found requests picked up by multiple vehicles
					// for example, r_(i,j) picked up by multiple vehicles
					// then we let y_(i,j) undecided
					// 2. If r_(i,j) is picked up by only one vehicle, let y_(i,j)=1
					// 3. force other y_(i,j)=0
					// in this way, we try to reallocate the request

					//==============start calling Gurobi===============

					printf("===> start LB calculation <===\n");

					int i, j, k, q;
					int status, nSolutions;

					GRBEnv *env = NULL;
					// GRBVar x[n][n][numV];
					// GRBVar y[n][n][numV];
					GRBVar s[n][numV];
					GRBVar u[n][n][n][numV];
					GRBVar theta[n][n][numV];

					GRBVar ***x = NULL;
					GRBVar ***y = NULL;
					x = new GRBVar **[n];
					y = new GRBVar **[n];
					for (i = 0; i < n; i++)
					{
						x[i] = new GRBVar *[n];
						y[i] = new GRBVar *[n];
						for (j = 0; j < n; j++)
						{
							x[i][j] = new GRBVar[numV];
							y[i][j] = new GRBVar[numV];
						}
					}

					try
					{
						env = new GRBEnv();
						GRBModel model = GRBModel(*env);

						// Create binary decision variables
						for (q = 0; q < numV; q++)
						{
							for (i = 0; i < n; i++)
							{
								s[i][q] = model.addVar(0.0, n, 0.0, GRB_CONTINUOUS,
																			 "s_" + itos(i) + "_" + itos(q));
								for (j = 0; j < n; j++)
								{
									x[i][j][q] = model.addVar(
											0.0, 1.0, 0, GRB_BINARY,
											"x_" + itos(i) + "_" + itos(j) + "_" + itos(q));
									y[i][j][q] = model.addVar(
											0.0, 1.0, 0, GRB_BINARY,
											"y_" + itos(i) + "_" + itos(j) + "_" + itos(q));
									theta[i][j][q] = model.addVar(
											0.0,
											GRB_INFINITY,
											0.0,
											GRB_CONTINUOUS,
											"theta_" + itos(i) + "_" + itos(j) + "_" + itos(q));

									for (k = 0; k < n; k++)
									{
										string s = "u_" + itos(i) + "_" + itos(j) + "_" + itos(k) + "_" + itos(q);
										u[i][j][k][q] = model.addVar(0.0, GRB_INFINITY, 0.0,
																								 GRB_CONTINUOUS, s);
									}
								}
							}

							int ogn = origin[q];
							for (i = 0; i < n; i++)
							{
								x[i][i][q].set(GRB_DoubleAttr_UB, 0);
								y[i][i][q].set(GRB_DoubleAttr_UB, 0);
								x[i][ogn][q].set(GRB_DoubleAttr_UB, 0);
								y[i][ogn][q].set(GRB_DoubleAttr_UB, 0);
								x[n - 1][i][q].set(GRB_DoubleAttr_UB, 0);
								y[n - 1][i][q].set(GRB_DoubleAttr_UB, 0);
								for (j = 0; j < n; j++)
									if (wt[i][j] > -0.00000001 && wt[i][j] < 0.00000001)
										y[i][j][q].set(GRB_DoubleAttr_UB, 0);
							}

							if (ADD_PREPROCESS)
							{
								for (int i = 0; i < n; i++)
								{
									vector<int> nbsTemp = allVehiclesInaccNeighbors[q][i];
									for (j = 0; j < nbsTemp.size(); j++)
									{
										x[i][nbsTemp[j]][q].set(GRB_DoubleAttr_UB, 0);
										// y[i][nbsTemp[j]].set(GRB_DoubleAttr_UB, 0);
										for (k = 0; k < n; k++)
										{
											u[i][nbsTemp[j]][k][q].set(GRB_DoubleAttr_UB, 0);
											u[i][k][nbsTemp[j]][q].set(GRB_DoubleAttr_UB, 0);
											u[k][nbsTemp[j]][i][q].set(GRB_DoubleAttr_UB, 0);
										}
									}
								}

								// int numViolatedTriples = 0;
								for (int i = 0; i < n; i++)
									if (i != ogn && i != endNode)
										for (auto &k : allVehiclesNodeNeighbors[q][i])
											for (auto &j : allVehiclesNodeNeighbors[q][k])
												if (i != j)
													if (dis[ogn][i] + dis[i][k] + dis[k][j] + dis[j][endNode] > disLimit)
													{
														u[i][j][k][q].set(GRB_DoubleAttr_UB, 0);
														// numViolatedTriples++;
													}
								// cout << "numViolatedTriples=" << numViolatedTriples << endl;
							}
						}

						int arcCandidates[n][n][numV];

						double sols_potentialArcs[n][numV];
						for (q = 0; q < numV; q++)
							for (i = 0; i < n; i++)
								sols_potentialArcs[i][q] = sols_d[i][q];

						if (ADD_POTENTIAL_ARCS)
						{
							//===> if cargo (i,j) is picked up multiple times,
							// select node x so that (i,x) and (x,j) is within distance limit and weight limit
							// if cargo (h,i) is not accepeted, (j,k) accepted, then also consider (h,x),(x,j)
							// if cargo (j,k) is not accepeted, (h,i) accepted, then also consier (i,x),(x,k)
							// if (h,i) (j,k) are not accepted, also consier (h,x),(x,k)

							cout << "====== check potential arcs ======" << endl;

							vector<double> distSumTemp; // temp: temparoray use
							int addedArcsCount = 0;

							for (q = 0; q < numV; q++)
								distSumTemp.push_back(0);

							for (q = 0; q < numV; q++)
								for (i = 0; i < n; i++)
									for (j = 0; j < n; j++)
									{
										arcCandidates[i][j][q] = 0;
										distSumTemp[q] += dis[i][j] * solx_numV_best_d[i][j][q];
									}

							for (int q1 = 0; q1 < numV - 1; q1++)
								for (int q2 = q1 + 1; q2 < numV; q2++)
									for (i = 0; i < n; i++)
										for (j = 0; j < n; j++)
										{

											vector<int> startNodesSet;
											vector<int> endNodesSet;
											if (soly_numV_best_d[i][j][q1] + soly_numV_best_d[i][j][q2] > 1.9) // if there is conflict
											{
												// printf("check arc %d to %d by veh%d and veh%d\n", i + 1, j + 1, q1 + 1, q2 + 1);

												startNodesSet.push_back(i);
												endNodesSet.push_back(j);

												vector<double> availDist = {disLimit - distSumTemp[q1] + dis[i][j], disLimit - distSumTemp[q2] + dis[i][j]};

												vector<int> vehiclesTemp = {q1, q2};

												// for (auto &vehTemp : vehiclesTemp)
												for (int vehIndex = 0; vehIndex < vehiclesTemp.size(); vehIndex++)
												{

													int vehTemp = vehiclesTemp[vehIndex];

													// cout << "===> vehTemp " << vehTemp + 1 << endl;

													for (int k = 0; k < n; k++)
														if (solx_numV_best_d[k][i][vehTemp] > 0.9)
															if (soly_numV_best_d[k][i][vehTemp] < 0.1) // if the previous visited arc is not an accepted cargo
															{
																startNodesSet.push_back(k);
																availDist[vehIndex] = availDist[vehIndex] + dis[k][i];
															}

													for (int k = 0; k < n; k++)
														if (solx_numV_best_d[j][k][vehTemp] > 0.9)
															if (soly_numV_best_d[j][k][vehTemp] < 0.1)
															{
																endNodesSet.push_back(k);
																availDist[vehIndex] = availDist[vehIndex] + dis[j][k];
															}

													for (auto &k1 : startNodesSet)
														for (auto &k2 : endNodesSet)
															for (int k3 = 0; k3 < n - 1; k3++)
																if ((sols_d[k3][vehTemp] < 0.001) && (dis[k1][k3] + dis[k3][k2] <= availDist[vehIndex]) && (k1 != k3) && (k2 != k3))
																{
																	if (k1 != endNode && k2 != origin[vehTemp] && k3 != origin[vehTemp])
																	{
																		arcCandidates[k1][k3][vehTemp] = 1;
																		addedArcsCount++;
																		// printf("startnode k1=%d, endnode k2=%d, k3=%d, veh=%d, (k1,k3) added\n", k1 + 1, k2 + 1, k3 + 1, vehTemp + 1);
																		arcCandidates[k3][k2][vehTemp] = 1;
																		addedArcsCount++;
																		// printf("startnode k1=%d, endnode k2=%d, k3=%d, veh=%d, (k3,k2) added\n", k1 + 1, k2 + 1, k3 + 1, vehTemp + 1);

																		sols_potentialArcs[k3][vehTemp] = -1;
																	}
																}
												}
											}
										}

							cout << "===> number of added potential arcs addedArcsCount = " << addedArcsCount << endl;

							for (q = 0; q < numV; q++)
							{
								cout << "vehicle " << q + 1 << endl;
								for (i = 0; i < n; i++)
									for (j = 0; j < n; j++)
										if (arcCandidates[i][j][q] > 0.9)
											printf("%d-%d, ", i + 1, j + 1);
								cout << endl;
							}

							// ======>  preset x variabless <======
							// it cause infeasibility in some instances
							{
								// print visited nodes
								for (q = 0; q < numV; q++)
								{
									cout << "===> vehilce " << q + 1 << " visited nodes and potentical arc nodes:" << endl;
									for (i = 0; i < n; i++)
										if (sols_potentialArcs[i][q] > 0.9 || (sols_potentialArcs[i][q] < -0.9 && sols_potentialArcs[i][q] > -1.1))
											cout << i + 1 << " ";
									cout << endl;
								}

								// cout << "check nodes visit status" << endl;

								for (q = 0; q < numV; q++)
								{
									// cout << "vehicle " << q + 1 << endl;
									for (i = 0; i < n; i++)
									{

										// ANECDOTE and WARNING:
										// use double sumY, do NOT use int sumY
										// when use "int sumY", somehow, a value of soly_numV_best_d[i][j][q]=-0.00000, y is double
										// When I added up all elements for sumY as below, although one value of soly_numV_best_d=1
										// but by adding -0.00000, the sumY<1, and c++ just round it to 0 by change double to int

										// double sumY = 0;
										// for (j = 0; j < n; j++)
										// {
										// 	sumY += soly_numV_best_d[i][j][q] + soly_numV_best_d[j][i][q];
										// 	// sumY += selectedPositiveProfitCargos[i][j] + selectedPositiveProfitCargos[j][i];
										// 	sumY += arcCandidates[i][j][q] + arcCandidates[j][i][q];
										// }

										// if a node has no cargos in or out, and not on route then it won't be visited
										// if (sumY < 0.01 && nodesVisitStatus[i][q] < 0.1)
										if (sols_potentialArcs[i][q] > -0.01 && sols_potentialArcs[i][q] < 0.01)
											for (j = 0; j < n; j++)
												x[j][i][q].set(GRB_DoubleAttr_UB, 0.0);
										else if (sols_potentialArcs[i][q] > 0.9) // node i is selected in original solution
										{
											for (j = 0; j < n; j++)
												if (sols_potentialArcs[j][q] > -0.01 && sols_potentialArcs[j][q] < 0.1) // node j is not in candidate visit node list
													x[i][j][q].set(GRB_DoubleAttr_UB, 0.0);
												else if (sols_potentialArcs[j][q] > 0.9)
												{ // node j is in candidate visite node list
													// both i and j are selected in original solution,but i is visited after j
													if (sols_potentialArcs[i][q] > sols_potentialArcs[j][q])
														x[i][j][q].set(GRB_DoubleAttr_UB, 0.0);
												}
												// In some instances, we might have a route like 1-5-7-10, request(1,5) and (5,7) are selected
												// and they are picked up by multiple vehicles
												// 1-4,4-5, 1-6,6-5 is ok to insert between 1 and 5
												// 5-4,4-10 is ok to be inserted between 5,10 (since 7-10 has no cargo, we check both 5-7 and 5-10)
												// so node 4 can be inserted between multiple nodes pair, so we need arcCandidates[i][j][q]
												// to know which arcs are ok to be inserted
												// if we only use sols_potentialArcs[i][q], in such case, we can't record more than one insert point
												else if (sols_potentialArcs[j][q] < -0.9 && sols_potentialArcs[j][q] > -1.1)
												{ // sols_potentialArcs[j][q]=-1 means node j is a candidate node for vehicle q's route
													if (arcCandidates[i][j][q] == 0)
														x[i][j][q].set(GRB_DoubleAttr_UB, 0.0);
												}
										}
										else if (sols_potentialArcs[i][q] < -0.9 && sols_potentialArcs[i][q] > -1.1)
											for (j = 0; j < n; j++)
												if (arcCandidates[i][j][q] == 0)
													x[i][j][q].set(GRB_DoubleAttr_UB, 0.0);
									}
								}
							}
						}

						//===> preset all y vars LB and UB
						for (i = 0; i < n; i++)
							for (j = 0; j < n; j++)
							{
								double sumY = 0;
								for (q = 0; q < numV; q++)
									sumY += soly_numV_best_d[i][j][q];

								if (sumY > 1.9)
									for (q = 0; q < numV; q++)
									{
										//  if (soly_numV_best_d[i][j][q] > 0.9)
										//  {
										//  y[i][j][q].set (GRB_DoubleAttr_ScenNLB, 1.0);
										//  y[i][j][q].set (GRB_DoubleAttr_ScenNUB, 1.0);
										//  }
										//  else
										// use numV best solution from solving single vehicle problem above
										// as the guide for calculating LB
										// if ADD_POTENTIAL_ARCS is false, selectedPositiveProfitCargos[i][j] is all zeros

										// if the cargo is not selected in LR dual, then do not consider them in MVBPMP
										// defaul arcCandidates are all zeros
										// if arcCandidates[i][j]==1, it means y[i][j]can be considered in model, so no need to set it to be zero
										if (soly_numV_best_d[i][j][q] < 0.1 && arcCandidates[i][j][q] < 0.1)
										// if (electedPositiveProfitCargos[i][j] < 0.1)
										{
											y[i][j][q].set(GRB_DoubleAttr_UB, 0.0);
										}
									}
								else if (sumY > 0.9 && sumY < 1.1)
								{
									for (q = 0; q < numV; q++)
										if (soly_numV_best_d[i][j][q] > 0.9 && soly_numV_best_d[i][j][q] < 1.1)
											y[i][j][q].set(GRB_DoubleAttr_LB, 1);
										else
											y[i][j][q].set(GRB_DoubleAttr_UB, 0);
								}
								else if (sumY < 0.1)
									for (q = 0; q < numV; q++)
										y[i][j][q].set(GRB_DoubleAttr_UB, 0);
							}

						// set up constraints
						GRBConstr *vehOriginConstr = 0;
						GRBConstr *vehDestConstr = 0;
						vehOriginConstr = new GRBConstr[numV];
						vehDestConstr = new GRBConstr[numV];

						for (q = 0; q < numV; q++)
						{
							int ogn = origin[q];

							// vehicle goes out of origins
							GRBLinExpr expr1 = 0.0;
							for (i = 0; i < n; i++)
								expr1 += x[ogn][i][q];

							// model.addConstr (expr1 == 1, "origin_" + itos (q));
							vehOriginConstr[q] = model.addConstr(expr1 == 1,
																									 "origin_" + itos(q));

							// vehicle goes back to node n
							GRBLinExpr expr2 = 0.0;
							for (i = 0; i < n - 1; i++)
								expr2 += x[i][n - 1][q];
							// model.addConstr (expr2 == 1, "destination_" + itos (q));
							vehDestConstr[q] = model.addConstr(expr2 == 1,
																								 "destination_" + itos(q));

							// flow conservation
							for (int k = 0; k < n - 1; k++)
							{
								if (k != ogn)
								{
									GRBLinExpr expr = 0;
									for (i = 0; i < n - 1; i++)
										expr += x[i][k][q];

									// BE CAREFUL!
									// I used j=1 to start which exclues node 0
									// which cause that the optimal profit is lower!
									for (j = 0; j < n; j++)
										expr -= x[k][j][q];
									model.addConstr(
											expr == 0,
											"flow_conservation_" + itos(k) + "_" + itos(q));
								}
							}

							// distance
							GRBLinExpr expr3 = 0.0;
							for (i = 0; i < n - 1; i++)
								for (j = 0; j < n; j++)
									expr3 += dis[i][j] * x[i][j][q];
							model.addConstr(expr3 <= route.DIS, "distance_" + itos(q));

							// node degree less than 1
							for (int j = 0; j < n - 1; j++)
							{
								GRBLinExpr expr = 0.0;
								for (i = 0; i < n - 1; i++)
									expr += x[i][j][q];
								model.addConstr(expr <= 1,
																"indegree_" + itos(j) + "_" + itos(q));
							}

							// subtour elimination
							for (i = 0; i < n - 1; i++)
								for (j = 0; j < n; j++)
								{
									GRBLinExpr expr = 0.0;
									expr += s[i][q] - s[j][q] + (n - 1) * x[i][j][q] + (n - 3) * x[j][i][q];
									model.addConstr(
											expr <= n - 2,
											"subtour_" + itos(i) + "_" + itos(j) + "_" + itos(q));
								}

							// arc flow
							for (i = 0; i < n - 1; i++)
								for (j = 0; j < n; j++)
								{
									GRBLinExpr expr = 0.0;

									expr += wt[i][j] * y[i][j][q] - theta[i][j][q];
									for (k = 0; k < n - 1; k++)
									{
										if (k != ogn)
											expr += u[i][k][j][q] + u[k][j][i][q] - u[i][j][k][q];
										else
											expr += u[k][j][i][q];
									}
									// when k==n-1
									if (j != n - 1)
										expr += u[i][n - 1][j][q];
									model.addConstr(
											expr == 0,
											"flow_" + itos(i) + "_" + itos(j) + "_" + itos(q));

									/*
									 expr += wt[i][j] * y[i][j][q] - theta[i][j][q];
									 for (k = 0; k < n; k++)
									 expr += u[i][k][j][q] + u[k][j][i][q] - u[i][j][k][q];
									 */
								}

							// arc flow upperbound
							for (i = 0; i < n - 1; i++)
								for (j = 0; j < n; j++)
								{
									GRBLinExpr expr = 0.0;
									expr += theta[i][j][q] - Q * x[i][j][q];
									model.addConstr(
											expr <= 0,
											"flowBound_" + itos(i) + "_" + itos(j) + "_" + itos(q));
								}
						}

						// one vehicle for one cargo
						// if (!ADD_CUTS){}
						for (i = 0; i < n - 1; i++)
							for (j = 0; j < n; j++)
							{
								GRBLinExpr expr = 0.0;
								for (q = 0; q < numV; q++)
									expr += y[i][j][q];
								model.addConstr(expr <= 1,
																"one-one" + itos(i) + "_" + itos(j));
							}

						// set up objective
						GRBLinExpr obj = 0.0;
						for (q = 0; q < numV; q++)
						{
							for (i = 0; i < n - 1; i++)
								for (j = 0; j < n; j++)
								{
									obj += price * dis[i][j] * wt[i][j] * y[i][j][q];
									obj -= cost * dis[i][j] * theta[i][j][q];
									obj -= cost * vw * dis[i][j] * x[i][j][q];
								}
						}

						model.setObjective(obj, GRB_MAXIMIZE);

						model.set(GRB_IntParam_OutputFlag, 0);

						model.set(GRB_IntParam_Threads, NUM_THREADS_VEH);

						// Optimize model
						model.optimize();

						// write model to file
						// model.write ("MVBPMP.lp");

						// Status checking
						status = model.get(GRB_IntAttr_Status);
						if (status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE || status == GRB_UNBOUNDED)
						{
							cout << "The model cannot be solved "
									 << "because it is infeasible or unbounded" << endl;
							// return 1;
							// exit(1);
						}
						if (status != GRB_OPTIMAL)
						{
							cout << "Optimization was stopped with status " << status << endl;
							exit(1);
							// return 1;
						}

						//=====================> READ SOLUTION FROM GUROBI <============================
						// Extract solution
						if (model.get(GRB_IntAttr_SolCount) > 0)
						{

							double runtime = model.get(GRB_DoubleAttr_Runtime);
							cout << "LBcalculationInLRLoopRuntime: " << runtime << " seconds" << endl;

							double objtemp = model.get(GRB_DoubleAttr_ObjVal);
							cout << "OBJ: " << objtemp << endl;
							if (objtemp > LBinLR)
							{
								LBinLR = objtemp;
								cout << "===> LBinLR is updated by LR LB calculation." << endl;
							}
							// update LB
							if (objtemp > LB)
							{
								LB = objtemp;
								cout << "===> LB is updated by LR LB calculation." << endl;
								LBinLR = objtemp;

								// reading and printing solutions
								{

									//====== start defining solutions for x, y, u, and s ======
									double ***solx = new double **[n]; // solution x for LR dual
									double ***soly = new double **[n];
									double ****solu = new double ***[n];
									double **sols = new double *[n];

									for (int i = 0; i < n; i++)
									{
										solx[i] = new double *[n];
										soly[i] = new double *[n];
										solu[i] = new double **[n];
										sols[i] = new double[numV];

										for (int j = 0; j < n; j++)
										{
											solx[i][j] = new double[numV];
											soly[i][j] = new double[numV];
											solu[i][j] = new double *[n];

											for (int k = 0; k < n; k++)
												solu[i][j][k] = new double[numV];
										}
									}

									//====== start reading solutions for x, y, u, and s ======
									// solx[i][j] = model.get (GRB_DoubleAttr_X, x[i][j],numV);
									// soly[i][j] = model.get (GRB_DoubleAttr_X, y[i][j],numV);
									for (q = 0; q < numV; q++)
										for (i = 0; i < n; i++)
										{
											sols[i][q] = s[i][q].get(GRB_DoubleAttr_X);
											for (j = 0; j < n; j++)
											{
												solx[i][j][q] = x[i][j][q].get(GRB_DoubleAttr_X);
												soly[i][j][q] = y[i][j][q].get(GRB_DoubleAttr_X);
												for (k = 0; k < n; k++)
													solu[i][j][k][q] = u[i][j][k][q].get(
															GRB_DoubleAttr_X);
											}
										}

									// store the solution as the best LB
									// storeBestLB(solx, soly, solu, sols, solx_best, soly_best,solu_best, sols_best, n);
									storeBestLB(solx, soly, solx_best, soly_best, n);

									printf("Selected arcs: \n");
									for (q = 0; q < numV; q++)
										for (i = 0; i < n; i++)
											for (j = 0; j < n; j++)
											{
												if (solx[i][j][q] > 0.9)
													printf("%d -- %d (%d)\n", i + 1, j + 1, q + 1);
											}

									printf("Selected requests: \n");
									for (q = 0; q < numV; q++)
										for (i = 0; i < n; i++)
											for (j = 0; j < n; j++)
											{
												if (soly[i][j][q] > 0.9)
													printf("%d -- %d (%d)\n", i + 1, j + 1, q + 1);
											}

									for (i = 0; i < n; i++)
									{
										for (j = 0; j < n; j++)
										{
											for (k = 0; k < n; k++)
												delete[] solu[i][j][k];

											delete[] solx[i][j];
											delete[] soly[i][j];
											delete[] solu[i][j];
										}
										delete[] solx[i];
										delete[] soly[i];
										delete[] solu[i];
										delete[] sols[i];
									}
									delete[] solx;
									delete[] soly;
									delete[] solu;
									delete[] sols;
								}
							}
							else
								cout << "===> LR LB calculation does not find better LB." << endl;
						}
					}
					catch (GRBException e)
					{
						cout << "Error number: " << e.getErrorCode() << endl;
						cout << e.getMessage() << endl;
					}
					catch (...)
					{
						cout << "Error during optimization" << endl;
					}
					for (i = 0; i < n; i++)
					{
						for (j = 0; j < n; j++)
						{
							delete[] x[i][j];
							delete[] y[i][j];
						}
						delete[] x[i];
						delete[] y[i];
					}
					delete[] x;
					delete[] y;
					delete env;

					// if (USE_LR_MULTIPLIER_TYPE_B)
					// 	updateLRmultiplierTypeB(soly_opt_d, LR_u, LR_rou, countItr, numArcs);

					double LR_miu;

					// if (USE_LR_MULTIPLIER_TYPE_C)
					// 	updateLRmultiplierTypeC(soly_opt_d, LR_lamda, &LR_miu, LR_u, LB,
					// 													profitInLRdual);

					// add comparison to LB_GRB
					if (LBinGRB > LB)
					{
						LB = LBinGRB;
						// cout print 4.834166 to 4.83417. So I use printf to print more digits
						//  cout << "===> When calculating LB in LR loop, find a better LB from LB_GRB = " << LBinGRB << endl;
						printf("===> When calculating LB in LR loop, find a better LB from LB_GRB = %lf\n", LBinGRB);
						// storeBestLB(solx_GRB, soly_GRB, solu_GRB, sols_GRB, solx_best, soly_best, solu_best, sols_best, n);
						storeBestLB(solx_GRB, soly_GRB, solx_best, soly_best, n);
					}

					if (USE_LR_MULTIPLIER_TYPE_C)
						updateLRmultiplierTypeC(soly_opt_d, LR_lamda, &LR_miu, LR_u, LB,
																		profitInLRdual);

					if (LR_miu <= LR_miu_tolerance)
					{
						cout << "=== will stop loop because LR_miu < LR_miu_tolerance "
								 << LR_miu_tolerance << endl;
						break;
					}

					//======> update UB <======

					UB = updateUB(UB, UBinLR, UBinGRB);

					if (LB == 0)
						LB = LBzero;

					printf("===> LB UB gap = %lf \n", (UB - LB) / LB);
					// if ((UB - LB) / LB <= LR_gap_tolerance)
					if (((UB - LB) / LB <= LR_gap_tolerance) && ((UB - LB) / LB >= -LR_gap_tolerance))
						cout << "=== will stop loop because gap <= LR_gap_tolerance and >= -LR_gap_tolerance "
								 << LR_gap_tolerance << endl;

					reportTime(endTimeOfLastIteration, endTimeOfLastIterationWallClock);
					endTimeOfLastIteration = clock();
					endTimeOfLastIterationWallClock = high_resolution_clock::now();
				}
			} // end of while loop

			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					for (int k = 0; k < n; k++)
						delete[] solu_d[i][j][k];

					delete[] solx_numV_best_d[i][j];
					delete[] solx_opt_d[i][j];
					delete[] soly_numV_best_d[i][j];
					delete[] soly_opt_d[i][j];
					delete[] solu_d[i][j];
				}
				delete[] solx_numV_best_d[i];
				delete[] solx_opt_d[i];
				delete[] soly_numV_best_d[i];
				delete[] soly_opt_d[i];
				delete[] solu_d[i];
				delete[] sols_d[i];
			}
			delete[] solx_numV_best_d;
			delete[] solx_opt_d;
			delete[] soly_numV_best_d;
			delete[] soly_opt_d;
			delete[] solu_d;
			delete[] sols_d;

			if (findOptimalSolution)
				cout << "**************** THE OPTIMAL SOLUTION ****************" << endl;
			else
				cout << "**************** THE BEST LB SOLUTION ****************" << endl;

			printVar(solx_best, soly_best, solu_best, origin);

			cout << endl
					 << "===> The total time: " << endl;
			reportTime(beginTime, beginWallClock);

			clock_t end = clock();
			double second = (double)(end - beginTime) / CLOCKS_PER_SEC;

			auto endWallClock = high_resolution_clock::now();
			auto elapsedWallClock = duration_cast<std::chrono::nanoseconds>(endWallClock - beginWallClock);

			cout << "===> SUMMARY:" << endl;
			printf("SolutionTimeCPU = %lf\n", second);
			printf("SolutionTimeWallClock = %.3f seconds\n",
						 elapsedWallClock.count() * 1e-9);
			printf("SolutionLB = %lf\n", LB);
			printf("SolutionUB = %lf\n", UB);
			printf("SolutionGap = %lf\n", (UB - LB) / LB);

			printBestLB(LB, LBinLR, LBinGRB);
			printBestUB(UB, UBinLR, UBinGRB);

			exit(1);
		}
		else
		{
			// solveMVBPMPinParallel(UB, UBinLR, UBinGRB, LB, LBinLR, LBinGRB, solx_GRB, soly_GRB,
			// 											allVehiclesInaccNeighbors, allVehiclesNodeNeighbors);
			printf("===> solving MVBPMP with GRB <===\n");

			int i, j, k, q;
			int status, nSolutions;

			GRBEnv *env = NULL;
			// GRBVar x[n][n][numV];
			// GRBVar y[n][n][numV];
			// GRBVar s[n][numV];
			// GRBVar u[n][n][n][numV];
			GRBVar theta[n][n][numV];

			GRBVar ***x = NULL;
			GRBVar ***y = NULL;
			x = new GRBVar **[n];
			y = new GRBVar **[n];
			for (i = 0; i < n; i++)
			{
				x[i] = new GRBVar *[n];
				y[i] = new GRBVar *[n];
				for (j = 0; j < n; j++)
				{
					x[i][j] = new GRBVar[numV];
					y[i][j] = new GRBVar[numV];
				}
			}

			GRBVar ****u = new GRBVar ***[n];
			GRBVar **s = new GRBVar *[n];

			for (int i = 0; i < n; i++)
			{
				u[i] = new GRBVar **[n];
				s[i] = new GRBVar[numV];
				for (int j = 0; j < n; j++)
				{
					u[i][j] = new GRBVar *[n];
					for (int k = 0; k < n; k++)
						u[i][j][k] = new GRBVar[numV];
				}
			}

			try
			{

				env = new GRBEnv();
				GRBModel model = GRBModel(*env);

				// Create binary decision variables
				for (q = 0; q < numV; q++)
				{
					for (i = 0; i < n; i++)
					{
						s[i][q] = model.addVar(0.0, n, 0.0, GRB_CONTINUOUS,
																	 "s_" + itos(i) + "_" + itos(q));
						for (j = 0; j < n; j++)
						{
							x[i][j][q] = model.addVar(
									0.0, 1.0, 0, GRB_BINARY,
									"x_" + itos(i) + "_" + itos(j) + "_" + itos(q));
							y[i][j][q] = model.addVar(
									0.0, 1.0, 0, GRB_BINARY,
									"y_" + itos(i) + "_" + itos(j) + "_" + itos(q));
							theta[i][j][q] = model.addVar(
									0.0,
									GRB_INFINITY,
									0.0,
									GRB_CONTINUOUS,
									"theta_" + itos(i) + "_" + itos(j) + "_" + itos(q));

							for (k = 0; k < n; k++)
							{
								string s = "u_" + itos(i) + "_" + itos(j) + "_" + itos(k) + "_" + itos(q);
								u[i][j][k][q] = model.addVar(0.0, GRB_INFINITY, 0.0,
																						 GRB_CONTINUOUS, s);
							}
						}
					}

					int ogn = origin[q];
					for (i = 0; i < n; i++)
					{
						x[i][i][q].set(GRB_DoubleAttr_UB, 0);
						y[i][i][q].set(GRB_DoubleAttr_UB, 0);
						x[i][ogn][q].set(GRB_DoubleAttr_UB, 0);
						y[i][ogn][q].set(GRB_DoubleAttr_UB, 0);
						x[n - 1][i][q].set(GRB_DoubleAttr_UB, 0);
						y[n - 1][i][q].set(GRB_DoubleAttr_UB, 0);
						for (j = 0; j < n; j++)
							if (wt[i][j] > -0.00000001 && wt[i][j] < 0.00000001)
								y[i][j][q].set(GRB_DoubleAttr_UB, 0);
					}

					if (ADD_PREPROCESS)
					{
						for (int i = 0; i < n; i++)
						{
							vector<int> nbsTemp = allVehiclesInaccNeighbors[q][i];
							for (j = 0; j < nbsTemp.size(); j++)
							{
								x[i][nbsTemp[j]][q].set(GRB_DoubleAttr_UB, 0);
								// y[i][nbsTemp[j]].set(GRB_DoubleAttr_UB, 0);
								for (k = 0; k < n; k++)
								{
									u[i][nbsTemp[j]][k][q].set(GRB_DoubleAttr_UB, 0);
									u[i][k][nbsTemp[j]][q].set(GRB_DoubleAttr_UB, 0);
									u[k][nbsTemp[j]][i][q].set(GRB_DoubleAttr_UB, 0);
								}
							}
						}

						// int numViolatedTriples = 0;
						for (int i = 0; i < n; i++)
							if (i != ogn && i != endNode)
								for (auto &k : allVehiclesNodeNeighbors[q][i])
									for (auto &j : allVehiclesNodeNeighbors[q][k])
										if (i != j)
											if (dis[ogn][i] + dis[i][k] + dis[k][j] + dis[j][endNode] > disLimit)
											{
												u[i][j][k][q].set(GRB_DoubleAttr_UB, 0);
												// numViolatedTriples++;
											}
						// cout << "numViolatedTriples=" << numViolatedTriples << endl;
					}
				}

				// set up constraints
				GRBConstr *vehOriginConstr = 0;
				GRBConstr *vehDestConstr = 0;
				vehOriginConstr = new GRBConstr[numV];
				vehDestConstr = new GRBConstr[numV];

				for (q = 0; q < numV; q++)
				{
					int ogn = origin[q];

					// vehicle goes out of origins
					GRBLinExpr expr1 = 0.0;
					for (i = 0; i < n; i++)
						expr1 += x[ogn][i][q];

					// model.addConstr (expr1 == 1, "origin_" + itos (q));
					vehOriginConstr[q] = model.addConstr(expr1 == 1,
																							 "origin_" + itos(q));

					// vehicle goes back to node n
					GRBLinExpr expr2 = 0.0;
					for (i = 0; i < n - 1; i++)
						expr2 += x[i][n - 1][q];
					// model.addConstr (expr2 == 1, "destination_" + itos (q));
					vehDestConstr[q] = model.addConstr(expr2 == 1,
																						 "destination_" + itos(q));

					// flow conservation
					for (int k = 0; k < n - 1; k++)
					{
						if (k != ogn)
						{
							GRBLinExpr expr = 0;
							for (i = 0; i < n - 1; i++)
								expr += x[i][k][q];

							// BE CAREFUL!
							// I used j=1 to start which exclues node 0
							// which cause that the optimal profit is lower!
							for (j = 0; j < n; j++)
								expr -= x[k][j][q];
							model.addConstr(
									expr == 0,
									"flow_conservation_" + itos(k) + "_" + itos(q));
						}
					}

					// distance
					GRBLinExpr expr3 = 0.0;
					for (i = 0; i < n - 1; i++)
						for (j = 0; j < n; j++)
							expr3 += dis[i][j] * x[i][j][q];
					model.addConstr(expr3 <= route.DIS, "distance_" + itos(q));

					// node degree less than 1
					for (int j = 0; j < n - 1; j++)
					{
						GRBLinExpr expr = 0.0;
						for (i = 0; i < n - 1; i++)
							expr += x[i][j][q];
						model.addConstr(expr <= 1,
														"indegree_" + itos(j) + "_" + itos(q));
					}

					// subtour elimination
					for (i = 0; i < n - 1; i++)
						for (j = 0; j < n; j++)
						{
							GRBLinExpr expr = 0.0;
							expr += s[i][q] - s[j][q] + (n - 1) * x[i][j][q] + (n - 3) * x[j][i][q];
							model.addConstr(
									expr <= n - 2,
									"subtour_" + itos(i) + "_" + itos(j) + "_" + itos(q));
						}

					// arc flow
					for (i = 0; i < n - 1; i++)
						for (j = 0; j < n; j++)
						{
							GRBLinExpr expr = 0.0;

							expr += wt[i][j] * y[i][j][q] - theta[i][j][q];
							for (k = 0; k < n - 1; k++)
							{
								if (k != ogn)
									expr += u[i][k][j][q] + u[k][j][i][q] - u[i][j][k][q];
								else
									expr += u[k][j][i][q];
							}
							// when k==n-1
							if (j != n - 1)
								expr += u[i][n - 1][j][q];
							model.addConstr(
									expr == 0,
									"flow_" + itos(i) + "_" + itos(j) + "_" + itos(q));

							/*
							 expr += wt[i][j] * y[i][j][q] - theta[i][j][q];
							 for (k = 0; k < n; k++)
							 expr += u[i][k][j][q] + u[k][j][i][q] - u[i][j][k][q];
							 */
						}

					// arc flow upperbound
					for (i = 0; i < n - 1; i++)
						for (j = 0; j < n; j++)
						{
							GRBLinExpr expr = 0.0;
							expr += theta[i][j][q] - Q * x[i][j][q];
							model.addConstr(
									expr <= 0,
									"flowBound_" + itos(i) + "_" + itos(j) + "_" + itos(q));
						}

					// add redundant constraint to test if same multi obj value
					// caused by binary variable
					/*
					 for (i = 0; i < n - 1; i++)
					 for (j = 0; j < n; j++)
					 for (k = 0; k < n - 1; k++)
					 {
					 GRBLinExpr expr = 0.0;
					 if (k != ogn)
					 expr += u[i][j][k][q] - x[i][k][q];
					 model.addConstr (
					 expr <= 0,
					 "unique_triples_" + itos (i) + "_" + itos (j)
					 + "_" + itos (q));
					 }
					 */
				}

				// one vehicle for one cargo
				// if (!ADD_CUTS){}
				for (i = 0; i < n - 1; i++)
					for (j = 0; j < n; j++)
					{
						GRBLinExpr expr = 0.0;
						for (q = 0; q < numV; q++)
							expr += y[i][j][q];
						model.addConstr(expr <= 1,
														"one-one" + itos(i) + "_" + itos(j));
					}

				// set up objective
				GRBLinExpr obj = 0.0;
				for (q = 0; q < numV; q++)
				{
					for (i = 0; i < n - 1; i++)
						for (j = 0; j < n; j++)
						{
							obj += price * dis[i][j] * wt[i][j] * y[i][j][q];
							obj -= cost * dis[i][j] * theta[i][j][q];
							obj -= cost * vw * dis[i][j] * x[i][j][q];
						}
				}

				model.setObjective(obj, GRB_MAXIMIZE);

				// Set callback function
				printIntSol cb = printIntSol(x, y, u, s, n, numV);
				model.setCallback(&cb);

				model.set(GRB_IntParam_OutputFlag, 0);

				model.set(GRB_IntParam_Threads, NUM_THREADS_VEH);

				// Optimize model
				model.optimize();

				// write model to file
				// model.write ("MVBPMP.lp");

				// Status checking
				status = model.get(GRB_IntAttr_Status);
				if (status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE || status == GRB_UNBOUNDED)
				{
					cout << "The model cannot be solved "
							 << "because it is infeasible or unbounded" << endl;
					exit(1);
				}
				if (status != GRB_OPTIMAL)
				{
					cout << "Optimization was stopped with status " << status << endl;
					exit(1);
				}
				if (status == GRB_OPTIMAL)
				{
					// Gurobi found optimal solution before LR converges
					// print out optimal solution and stop
					cout << "Optimization was completed by Gurobi. Optimal solution found. " << endl;

					if (model.get(GRB_IntAttr_SolCount) > 0)
					{

						double runtime = model.get(GRB_DoubleAttr_Runtime);
						cout << "Gurobi in parallel Runtime: " << runtime << " seconds" << endl;

						double objtemp = model.get(GRB_DoubleAttr_ObjVal);
						UB = objtemp;
						LB = objtemp;

						//====== start defining solutions for x, y, u, and s ======
						double ***solx = new double **[n]; // solution x for LR dual
						double ***soly = new double **[n];
						double ****solu = new double ***[n];
						double **sols = new double *[n];

						for (int i = 0; i < n; i++)
						{
							solx[i] = new double *[n];
							soly[i] = new double *[n];
							solu[i] = new double **[n];
							sols[i] = new double[numV];

							for (int j = 0; j < n; j++)
							{
								solx[i][j] = new double[numV];
								soly[i][j] = new double[numV];
								solu[i][j] = new double *[n];

								for (int k = 0; k < n; k++)
									solu[i][j][k] = new double[numV];
							}
						}

						//====== start reading solutions for x, y, u, and s ======
						// solx[i][j] = model.get (GRB_DoubleAttr_X, x[i][j],numV);
						// soly[i][j] = model.get (GRB_DoubleAttr_X, y[i][j],numV);
						for (q = 0; q < numV; q++)
							for (i = 0; i < n; i++)
							{
								sols[i][q] = s[i][q].get(GRB_DoubleAttr_X);
								for (j = 0; j < n; j++)
								{
									solx[i][j][q] = x[i][j][q].get(GRB_DoubleAttr_X);
									soly[i][j][q] = y[i][j][q].get(GRB_DoubleAttr_X);
									for (k = 0; k < n; k++)
										solu[i][j][k][q] = u[i][j][k][q].get(
												GRB_DoubleAttr_X);
								}
							}

						// store the solution as the best LB
						// storeBestLB(solx, soly, solu, sols, solx_best, soly_best,solu_best, sols_best, n);
						storeBestLB(solx, soly, solx_best, soly_best, n);

						for (i = 0; i < n; i++)
						{
							for (j = 0; j < n; j++)
							{
								for (k = 0; k < n; k++)
									delete[] solu[i][j][k];

								delete[] solx[i][j];
								delete[] soly[i][j];
								delete[] solu[i][j];
							}
							delete[] solx[i];
							delete[] soly[i];
							delete[] solu[i];
							delete[] sols[i];
						}
						delete[] solx;
						delete[] soly;
						delete[] solu;
						delete[] sols;
					}
				}
			}
			catch (GRBException e)
			{
				cout << "Error number: " << e.getErrorCode() << endl;
				cout << e.getMessage() << endl;
			}
			catch (...)
			{
				cout << "Error during optimization" << endl;
			}

			delete env;

			cout << "**************** THE OPTIMAL SOLUTION ****************" << endl;
			printVar(solx_best, soly_best, solu_best, origin);
			cout << endl
					 << "===> The total time: " << endl;
			reportTime(beginTime, beginWallClock);

			clock_t end = clock();
			double second = (double)(end - beginTime) / CLOCKS_PER_SEC;

			auto endWallClock = high_resolution_clock::now();
			auto elapsedWallClock = duration_cast<std::chrono::nanoseconds>(endWallClock - beginWallClock);

			cout << "===> SUMMARY:" << endl;
			printf("SolutionTimeCPU = %lf\n", second);
			printf("SolutionTimeWallClock = %.3f seconds\n",
						 elapsedWallClock.count() * 1e-9);
			printf("SolutionLB = %lf\n", LB);
			printf("SolutionUB = %lf\n", UB);
			printf("SolutionGap = %lf\n", (UB - LB) / LB);
			// if program finish in this thread (thread 1), it means Gurobi finishes and find optimal solution before LR finishes
			// cout << "Best LB is found by GRB in parallel (find optimal solution)" << endl;
			cout << "BestLowerBoundSolver = GRB" << endl;
			cout << "BestUpperBoundSolver = GRB" << endl;
			exit(1);
			// return 0;
		}
	}

	// if (findOptimalSolution)
	// 	cout << "**************** THE OPTIMAL SOLUTION ****************" << endl;
	// else
	// 	cout << "**************** THE BEST LB SOLUTION ****************" << endl;

	// printVar(solx_best, soly_best, solu_best, origin);

	// cout << endl
	// 		 << "===> The total time: " << endl;
	// reportTime(beginTime, beginWallClock);

	// clock_t end = clock();
	// double second = (double)(end - beginTime) / CLOCKS_PER_SEC;

	// auto endWallClock = high_resolution_clock::now();
	// auto elapsedWallClock = duration_cast<std::chrono::nanoseconds>(endWallClock - beginWallClock);

	// cout << "===> SUMMARY:" << endl;
	// printf("SolutionTimeCPU = %lf\n", second);
	// printf("SolutionTimeWallClock = %.3f seconds\n",
	// 			 elapsedWallClock.count() * 1e-9);
	// printf("SolutionLB = %lf\n", LB);
	// printf("SolutionUB = %lf\n", UB);
	// printf("SolutionGap = %lf\n", (UB - LB) / LB);

	//*********************** START STAGE TWO ****************************
	//********************************************************************
	// the model is supposed to be the same with the MVBPMP in LR process
	// but 1. LR model has pre-setup of variables upperbound
	// this model in stage two doesn't, but only the initial solution from LR
	// 2. this model set up the timeLimit: LR + stage2 = 24 hours

	// if (DO_STAGE_TWO && !findOptimalSolution)
	// {
	// 	// this part is skipped
	// }

	return 0;
}

// double updateUB(double UB)
// {
// 	double bestValue;
// 	if (UB > UBinGRB && UBinLR > UBinGRB)
// 	{
// 		bestValue = UBinGRB;
// 		printf("===> UBinGRB = %lf is a better upper bound\n", UBinGRB);
// 	}
// 	else if (UB > UBinLR && UBinGRB > UBinLR)
// 	{
// 		bestValue = UBinLR;
// 		printf("===> UBinLR = %lf is a better upper bound\n", UBinLR);
// 	}
// 	else if (UB > UBinLR && UBinGRB > UBinLR - 0.00000001 && UBinGRB < UBinLR + 0.00000001)
// 	{
// 		bestValue = UBinGRB;
// 		printf("===> both UBinLR = %lf and UBinGRB = %lf are better upper bounds\n", UBinLR, UBinGRB);
// 	}
// 	else
// 	{
// 		cout << "Did not update UB." << endl;
// 		printf("UB=%lf, UBinLR=%lf, UBinGRB=%lf\n", UB, UBinLR, UBinGRB);
// 		bestValue = UB;
// 	}

// 	return bestValue;
// }

// void printBestLB(double LB)
// {
// 	if (LBinLR == LB && LBinGRB < LB)
// 		// cout << "Best LB is found by LR" << endl;
// 		cout << "BestLowerBoundSolver = LR" << endl;
// 	else if (LBinLR < LB && LBinGRB == LB)
// 		// cout << "Best LB is found by GRB in parallel" << endl;
// 		cout << "BestLowerBoundSolver = GRB" << endl;
// 	else if (LBinLR == LB && LBinGRB == LB)
// 		// cout << "Best LB is found by both LR and GRB in parallel" << endl;
// 		cout << "BestLowerBoundSolver = LR_GRB" << endl;
// 	else if (LBinLR > LB - 0.00000001 && LBinLR < LB + 0.00000001 && LBinGRB > LB - 0.00000001 && LBinGRB < LB + 0.00000001)
// 		cout << "BestLowerBoundSolver = LR_GRB" << endl;
// 	else
// 		printf("ERROR: LB=%lf, LBinLR=%lf, LBinGRB=%lf\n", LB, LBinLR, LBinGRB);
// }

// void printBestUB(double UB)
// {
// 	if (UBinLR == UB && UBinGRB > UB)
// 		// cout << "Best UB is found by LR" << endl;
// 		cout << "BestUpperBoundSolver = LR" << endl;
// 	else if (UBinLR > UB && UBinGRB == UB)
// 		// cout << "Best UB is found by GRB in parallel" << endl;
// 		cout << "BestUpperBoundSolver = GRB" << endl;
// 	else if (UBinLR == UB && UBinGRB == UB)
// 		// cout << "Best UB is found by both LR and GRB in parallel" << endl;
// 		cout << "BestUpperBoundSolver = LR_GRB" << endl;
// 	else if (UBinLR > UB - 0.00000001 && UBinLR < UB + 0.00000001 && UBinGRB > UB - 0.00000001 && UBinGRB < UB + 0.00000001)
// 		cout << "BestUpperBoundSolver = LR_GRB" << endl;
// 	else
// 		printf("ERROR: UB=%lf, UBinLR=%lf, UBinGRB=%lf\n", UB, UBinLR, UBinGRB);
// }