// This is an example from https://groups.google.com/g/gurobi/c/pkBNfu-iX0k

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// although we read vehicle data, we only consider one vehicle scenario whose route is from 0 to n-1 right now

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
#include <ctime>  //to measure CPU time
#include <chrono> //to measure run time
#include <time.h>

#include <iomanip>

#include "readData.h"
#include "ga.h"

using namespace std;
using namespace std::chrono;
using namespace std::string_literals;

bool PRINT_FOR_DEBUG = false;

bool FORCE_TRIANGLE_INEQUALITY = true;

bool toRemove_addInitialSolutionInPricingProblem = true;

// double bigM = 10000000;
double epsilon = numeric_limits<double>::epsilon();

// string
// itos(int i)
// {
// 	stringstream s;
// 	s << i;
// 	return s.str();
// }

void reportTime(clock_t, auto);
set<vector<int>> findTriangelInequalityViolation(double **, int);

int main(int argc, char *argv[])
{

	//********* only read graph info and vehicle info by arguments
	//********* doesn't go through the graph info in the data folder
	if (argc != 2)
	{
		cout
			// << "Usage: ./bpmp_dw-cg.x nodesDataNameAndPath numberOfVehicles vehicleDataNameAndPath"<< endl;
			<< "Usage: ./bpmp_dw-cg.x nodesDataNameAndPath"
			<< endl;
		return 1;
	}

	if (argc == 2)
	{
		cout << "Nodes Data: " << argv[1] << endl;
	}

	clock_t beginTime = clock();

	auto beginWallClock = high_resolution_clock::now();

	readData route;

	string filename = argv[1];
	printf("reading file %s ...\n", filename.c_str());
	route.readSingleFile(filename);

	// route.printStats();
	// route.printData ();

	int n = route.numOfNode;

	// double wt[n][n];
	// double dis_v[n][n];

	double **wt = new double *[n];
	double **dis_v = new double *[n];

	for (int i = 0; i < n; i++)
	{
		wt[i] = new double[n];
		dis_v[i] = new double[n];
	}

	for (int i = 1; i <= n; i++)
		for (int j = 1; j <= n; j++)
		{
			wt[i - 1][j - 1] = route.w[i][j];
			dis_v[i - 1][j - 1] = route.d[i][j];
		}
	double price = route.priceCharged;
	double cost = route.travelCost;
	double vw = route.vehicleWeight;
	double capacity = route.totalCapacity;

	//================> end of reading data <================

	if (FORCE_TRIANGLE_INEQUALITY)
	{
		cout << "==> FORCE_TRIANGLE_INEQUALITY: " << FORCE_TRIANGLE_INEQUALITY << endl;
		set<vector<int>> TIV = findTriangelInequalityViolation(dis_v, n);
		if (TIV.size() > 0)
		{
			for (auto &temp : TIV)
			{
				int tempSize = temp.size();
				if (tempSize != 3)
				{
					printf("ERROR: TIV element always have vectors of size 3.");
					exit(1);
				}
				// TIV:{i,j,k} d(i,j)+d(j,k)-d(i,k)<0
				double tempDiff = dis_v[temp[0]][temp[1]] + dis_v[temp[1]][temp[2]] - dis_v[temp[0]][temp[2]];
				dis_v[temp[0]][temp[2]] = dis_v[temp[0]][temp[2]] + tempDiff;
			}
		}
		TIV = findTriangelInequalityViolation(dis_v, n);
		printf("===> after modifying distance,TIV.size=%d\n", TIV.size());
	}

	double xCoeff[n][n] = {{1000000.000000, -0.677700, -0.760000, -0.115100, -1.294900, -0.814800, -0.190600, -0.723800, -0.387900, -0.829800, -0.146000, -0.873300, 0.420900, -0.980100, -1.429800, -0.642300, -0.690900, -0.150100, -0.627200, -1.000000},
						   {1000000.000000, 0.000000, -0.397600, -0.729100, -0.984700, -1.025900, -0.513000, -1.643900, -0.885500, -1.193600, -0.544800, -1.383000, -1.098400, -0.828200, -1.202400, -0.843800, -0.242900, -0.742100, -0.712900, -0.856300},
						   {1000000.000000, -0.397600, 0.000000, -1.267300, -0.597400, -0.739600, -1.661600, -0.999000, -1.790300, -0.945000, -0.616600, 0.741600, -1.528300, -0.435700, -0.788300, -0.546900, -0.637500, -1.122000, -0.432100, -0.461900},
						   {1000000.000000, -0.637900, -0.667300, 0.000000, -1.140200, -0.710500, -0.224400, -0.657300, -0.315800, -0.744100, -0.100700, -0.811800, -0.484100, -0.851800, -1.309200, -0.530000, -0.691600, -0.107300, -0.512100, -0.886000},
						   {1000000.000000, -1.108300, -0.597400, -1.140200, 0.000000, -1.099500, -2.622200, -1.713200, -1.098100, -1.027900, -1.986400, -1.077000, -1.769900, -0.307400, -0.229400, -0.737500, -1.316300, -1.297400, -0.690300, -0.278600},
						   {1000000.000000, -1.358500, -0.722700, -1.805200, -0.740700, 0.000000, -0.903200, -0.938400, -0.484000, -0.261900, -0.961900, -0.344400, -0.956400, -0.463200, -0.721600, -0.208200, -1.207700, -1.179500, -0.309700, -0.479100},
						   {1000000.000000, -0.513000, -0.710600, -0.224400, -1.224200, -0.920000, 0.000000, -0.881700, -0.540200, -0.963700, -0.141900, -1.057100, -1.878600, -0.962000, -1.364300, -0.708600, -0.501500, -0.315900, -0.661700, -0.996500},
						   {1000000.000000, -1.189100, -0.999000, -0.657300, -1.145500, -0.409200, -0.881700, 0.000000, -0.341500, -0.191400, -0.749000, -0.177500, -1.028600, -0.690600, -1.156100, -0.466100, -1.675400, -0.628900, -0.577100, -0.875300},
						   {1000000.000000, -0.885500, -0.777500, -0.315800, -1.095700, -0.497500, -0.540200, -0.341500, 0.000000, -0.444700, -0.452800, -0.485400, -0.480400, -0.788500, -1.147800, -0.378700, -1.250500, -0.237900, -0.435200, -0.818600},
						   {1000000.000000, -1.915000, -0.947400, -1.251700, -0.999000, -0.261900, -0.984100, -0.191400, -1.503100, 0.000000, -0.823500, -0.084900, -1.409200, -0.725100, -1.254100, -0.394800, -1.906300, -0.680600, -0.513100, -0.742100},
						   {1000000.000000, -0.544800, -0.616600, -0.100700, -1.127200, -0.762700, -0.141900, -0.749000, -0.408700, -0.823500, 0.000000, -0.878900, -0.565200, -0.873200, -1.284400, -0.571100, -0.640400, -0.206300, -0.533100, -0.885800},
						   {1000000.000000, -2.085200, -1.029600, -0.922000, -1.159800, -0.344400, -1.017100, 0.177500, -0.533600, -0.084900, -0.878900, 0.000000, -0.851100, -0.807500, -1.488400, -0.478000, -1.417300, -0.771500, -0.597500, -0.822500},
						   {1000000.000000, -1.181200, -1.145600, -0.414500, -1.567100, -0.956400, -1.873900, -0.673600, -0.480400, -0.847800, -0.560500, -0.851100, 0.000000, -1.254200, -1.656800, -0.856200, -1.857600, -0.401600, -0.895500, -1.286000},
						   {1000000.000000, -0.889400, -0.435700, -1.079800, -0.307400, -0.463200, -1.817600, -1.226500, -0.788500, -0.725100, -1.332500, -0.808700, -1.430600, 0.000000, -0.402400, -0.418000, -1.199700, -0.882200, -0.358700, -0.034000},
						   {1000000.000000, -1.202400, -0.788300, -1.537200, -0.229400, -2.449600, -2.221100, -2.430600, -1.153800, -2.982100, -1.743700, -1.765500, -1.625600, -0.402400, 0.000000, -0.769600, -1.428000, -1.428700, -0.964600, -0.367800},
						   {1000000.000000, -1.292600, -0.550500, -1.269200, -0.724300, -0.208200, -0.708600, -0.481800, -0.378700, -0.394800, -0.571100, -0.478000, -0.856200, -0.418000, -0.798400, 0.000000, -1.000000, -1.646200, -0.118400, -0.445800},
						   {1000000.000000, -0.242900, -0.673500, -1.063600, -1.227500, -1.285700, -0.501500, -1.997900, -1.514500, -2.179900, -0.962900, -2.419300, -1.091700, -1.061700, -1.464000, -1.274800, 0.000000, -1.051700, -0.904500, -1.091200},
						   {1000000.000000, -0.742100, -0.744300, -0.105800, -1.179800, -0.493100, -0.315900, -0.574900, -0.237900, -0.704600, -0.205800, -0.723300, -0.403000, -0.882200, -1.271500, -0.523000, -0.796100, 0.000000, -0.528900, -0.915700},
						   {1000000.000000, -1.045500, -0.432100, -0.512100, -0.663900, -0.309700, -0.816500, -0.577100, -0.370400, -0.513100, -0.532700, -0.596300, -1.071900, -0.358700, -0.742600, -0.118400, -1.852700, -0.528900, 0.000000, -0.390700},
						   {1000000.000000, 1000000.000000, 1000000.000000, 1000000.000000, 1000000.000000, 1000000.000000, 1000000.000000, 1000000.000000, 1000000.000000, 1000000.000000, 1000000.000000, 1000000.000000, 1000000.000000, 1000000.000000, 1000000.000000, 1000000.000000, 1000000.000000, 1000000.000000, 1000000.000000, 1000000.000000}};

	GRBEnv *env = 0;

	try
	{
		int i, j, k;

		env = new GRBEnv();

		{ //==> start

			if (PRINT_FOR_DEBUG)
				cout << "Start Pricing Problem \n";

			GRBModel PPmodel = GRBModel(*env);
			// PP.set (GRB_IntAttr_ModelSense, -1);

			GRBVar **x = nullptr;
			GRBVar s[n];

			x = new GRBVar *[n];
			for (i = 0; i < n; i++)
				x[i] = new GRBVar[n];

			for (i = 0; i < n; i++)
			{
				s[i] = PPmodel.addVar(0.0, n, 0.0, GRB_CONTINUOUS,
									  "s_" + itos(i));
				for (j = 0; j < n; j++)
					x[i][j] = PPmodel.addVar(0.0, 1.0, 0, GRB_BINARY,
											 "x_" + itos(i) + "_" + itos(j));
			}

			for (i = 0; i < n; i++)
			{
				x[i][i].set(GRB_DoubleAttr_UB, 0);
				x[i][0].set(GRB_DoubleAttr_UB, 0);
				x[n - 1][i].set(GRB_DoubleAttr_UB, 0);
			}

			if (PRINT_FOR_DEBUG)
				cout << "PP variables added! \n";

			// vehicle goes out of node 1
			GRBLinExpr expr1 = 0.0;
			for (i = 1; i < n; i++)
				expr1 += x[0][i];
			PPmodel.addConstr(expr1 == 1, "origin");

			// vehicle goes back to node n
			GRBLinExpr expr2 = 0.0;
			for (i = 0; i < n - 1; i++)
				expr2 += x[i][n - 1];
			PPmodel.addConstr(expr2 == 1, "destination");

			// flow conservation
			for (int k = 1; k < n - 1; k++)
			{
				GRBLinExpr expr = 0;
				for (i = 0; i < n - 1; i++)
					expr += x[i][k];
				for (j = 1; j < n; j++)
					expr -= x[k][j];
				PPmodel.addConstr(expr == 0, "flow_conservation_" + itos(k));
			}

			// distance
			GRBLinExpr expr3 = 0.0;
			for (i = 0; i < n - 1; i++)
				for (j = 1; j < n; j++)
					expr3 += dis_v[i][j] * x[i][j];
			PPmodel.addConstr(expr3 <= route.DIS, "distance");

			// node degree less than 1
			for (int k = 1; k < n - 1; k++)
			{
				GRBLinExpr expr = 0.0;
				for (i = 0; i < n - 1; i++)
					expr += x[i][k];
				PPmodel.addConstr(expr <= 1, "indegree_" + itos(k));
			}

			// subtour elimination
			for (i = 0; i < n - 1; i++)
				for (j = 1; j < n; j++)
				{
					GRBLinExpr expr = 0.0;
					expr += s[i] - s[j] + (n - 1) * x[i][j] + (n - 3) * x[j][i];
					PPmodel.addConstr(expr <= n - 2,
									  "s_" + itos(i) + "_" + itos(j));
				}

			if (PRINT_FOR_DEBUG)
				cout << "PP constraints added! \n";

			// set objective
			// GRBLinExpr obj = 0.0;
			// for (i = 0; i < n - 1; i++)
			// 	for (j = 1; j < n; j++)
			// 	{
			// 		obj -= cost * vw * dis_v[i][j] * x[i][j];
			// 		obj += capacity * pi[i][j] * x[i][j];
			// 	}
			// obj -= pi_lambda_sum;
			// PPmodel.setObjective(obj, GRB_MAXIMIZE);

			GRBLinExpr obj = 0.0;
			for (i = 0; i < n - 1; i++)
				for (j = 1; j < n; j++)
				{
					// obj += cost * vw * dis_v[i][j] * x[i][j];
					// obj -= capacity * pi[i][j] * x[i][j];
					obj += xCoeff[i][j] * x[i][j];
				}
			// obj += pi_lambda_sum;
			obj += 1.9994;
			PPmodel.setObjective(obj, GRB_MINIMIZE);

			if (toRemove_addInitialSolutionInPricingProblem)
			{
				for (i = 0; i < n - 1; i++)
					for (j = 1; j < n; j++)
						x[i][j].set(GRB_DoubleAttr_Start, 0.0);

				x[0][17].set(GRB_DoubleAttr_Start, 1);
				x[17][3].set(GRB_DoubleAttr_Start, 1);
				x[3][8].set(GRB_DoubleAttr_Start, 1);
				x[8][10].set(GRB_DoubleAttr_Start, 1);
				x[10][15].set(GRB_DoubleAttr_Start, 1);
				x[15][19].set(GRB_DoubleAttr_Start, 1);
			}
			if (PRINT_FOR_DEBUG)
				cout << "PP objective function is set up! \n";

			PPmodel.getEnv().set(GRB_IntParam_OutputFlag, 0); // Silent Mode
			PPmodel.set(GRB_DoubleParam_OptimalityTol, 1e-9); // Optimality tolerance

			PPmodel.write("DW_pricingModel.lp");
			PPmodel.optimize();

			double objValue_PP = PPmodel.get(GRB_DoubleAttr_ObjVal);
			objValue_PP = -objValue_PP;
			cout << "objective of PP problem (gurobi): " << objValue_PP << "\n";

			for (i = 0; i < n; i++)
				delete[] x[i];
			delete[] x;
			x = nullptr;

		} //==> end
	}
	catch (GRBException e)
	{
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...)
	{
		cout << "Exception during optimization" << endl;
	}

	delete env;
	env = nullptr;

	reportTime(beginTime, beginWallClock);

	return 0;
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

set<vector<int>> findTriangelInequalityViolation(double **dis_v, int n)
{
	set<vector<int>> TIV; // triangel inequality violation
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
			{
				double temp = dis_v[i][j] + dis_v[j][k] - dis_v[i][k];
				if (temp < -0.0000000001)
				{
					vector<int> nodesTemp;
					nodesTemp.push_back(i);
					nodesTemp.push_back(j);
					nodesTemp.push_back(k);
					TIV.insert(nodesTemp);
					printf("d(%d,%d)+d(%d,%d)-d(%d,%d)=%lf\n", i, j, j, k, i, k, temp);
				}
			}
	return TIV;
}