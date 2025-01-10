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
// #include "ga_objfunc_greedy.h"
#include "ga.h"
#include "dominance.h"
// #include "dominance_inaccesiblenode.h"

// #include <omp.h>

using namespace std;
using namespace std::chrono;
using namespace std::string_literals;
// to generate vehicle data file name

int numV = 3; // number of vehicles; will be overwritten when reading vehicle file
// number of nodes in instance
// in main(), we also use n to present number of nodes
// since it's easy to read variables with short subscript, like x[n][n][numV]
int NUM_NODES;

bool PRINT_FOR_DEBUG = false;

bool USE_GA_FOR_INIT_SOL = false; // the experiment shows that a good initial solution doesn't reduce the number of iterations of sub-problem
bool FORCE_TRIANGLE_INEQUALITY = true;
bool USE_DOMINACE_METHOD_OR_SUB_PROBLEM = true; // if it is false, Gurobi will be used for pricing problem

bool testOneRouteObjAndDist = false;
bool toRemove_addInitialSolutionInPricingProblem = false;

// double bigM = 10000000;
double epsilon = numeric_limits<double>::epsilon();

// itos is defined in ga.h as well
//  string
//  itos(int i)
//  {
//  	stringstream s;
//  	s << i;
//  	return s.str();
//  }

void reportTime(clock_t, auto);
set<vector<int>> findTriangelInequalityViolation(double **, int);

int main(int argc, char *argv[])
{

	//********* only read graph info and vehicle info by arguments
	//********* doesn't go through the graph info in the data folder
	if (argc != 4)
	{
		cout
			<< "Usage: ./bpmp_dw-cg.x nodesDataNameAndPath numberOfVehicles vehicleDataNameAndPath"
			<< endl;
		return 1;
	}

	if (argc == 4)
	{
		cout << "Nodes Data: " << argv[1] << endl;
		cout << "Number of Vehicles: " << argv[2] << endl;
		cout << "Vehicles Data: " << argv[3] << endl;
		numV = stoi(argv[2]);
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
	NUM_NODES = route.numOfNode;

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

	GRBEnv *env = 0;

	try
	{
		int i, j, k;
		// Patterns
		// double pat[][2];

		//****** newCol[n-1][j]=0 and newCol[i][0]=0
		double **newCol = new double *[n];
		for (i = 0; i < n; i++)
			newCol[i] = new double[n];

		for (i = 0; i < n; i++)
		{
			newCol[n - 1][i] = 0;
			newCol[i][0] = 0;

			for (j = 0; j < n; j++)
				newCol[i][j] = 0;
		}

		//======> initial route solution <======
		cout << "before using ga" << endl;
		if (USE_GA_FOR_INIT_SOL)
		{
			double disLimit = (double)route.DIS;

			// right now it only returns the selected route, no customers and profit info
			vector<int> finalRoute = runga(n, price, cost, capacity, disLimit, vw, wt, dis_v);

			cout << "========> selected route using genetic algorithm:" << endl;
			for (int i : finalRoute)
				cout << i << ",";
			cout << endl;

			for (i = 0; i < finalRoute.size() - 1; i++)
				newCol[finalRoute[i]][finalRoute[i + 1]] = 1;
		}
		else
		{
			newCol[0][n - 1] = 1;
		}
		//======> end of initial route selection <======

		if (PRINT_FOR_DEBUG)
		{
			cout << "=> initial patterns: \n";
			for (i = 0; i < n; i++)
			{
				for (j = 0; j < n; j++)
					cout << newCol[i][j] << "  ";
				cout << "" << endl;
			}
		}

		env = new GRBEnv();
		GRBModel Master = GRBModel(*env);

		//****** Variables of master problem
		//*** set up var lambda
		vector<GRBVar> lambda;
		lambda.push_back(
			Master.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
						  "lambda_" + itos(0)));
		if (PRINT_FOR_DEBUG)
			cout << "variable lambda is added \n";

		//*** set up var y and triples variable z
		GRBVar **y = nullptr;
		GRBVar ***z = nullptr;
		y = new GRBVar *[n];
		z = new GRBVar **[n];
		for (i = 0; i < n; i++)
		{
			y[i] = new GRBVar[n];
			z[i] = new GRBVar *[n];
			for (int j = 0; j < n; j++)
				z[i][j] = new GRBVar[n];
		}

		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
			{
				// y[i][j] = Master.addVar (0.0, 1.0, 0, GRB_BINARY, "y_" + itos (i) + "_" + itos (j));
				y[i][j] = Master.addVar(0.0, 1.0, 0, GRB_CONTINUOUS,
										"y_" + itos(i) + "_" + itos(j));
				for (k = 0; k < n; k++)
				{
					string s = "z_" + itos(i) + "_" + itos(j) + "_" + itos(k);
					//*** in DW-CG, we force u<=capacity (capacity is the vehicle's capacity)
					z[i][j][k] = Master.addVar(0.0, capacity, 0.0, GRB_CONTINUOUS, s);
				}
			}
		if (PRINT_FOR_DEBUG)
			cout << "variable y and z are added \n";

		for (i = 0; i < n; i++)
		{
			y[i][i].set(GRB_DoubleAttr_UB, 0);
			y[i][0].set(GRB_DoubleAttr_UB, 0);
			y[n - 1][i].set(GRB_DoubleAttr_UB, 0);

			for (j = 0; j < n; j++)
			{
				z[i][j][0].set(GRB_DoubleAttr_UB, 0);
				z[i][j][n - 1].set(GRB_DoubleAttr_UB, 0);
				z[i][0][j].set(GRB_DoubleAttr_UB, 0);
				z[n - 1][i][j].set(GRB_DoubleAttr_UB, 0);
			}
		}

		// Constraints of master problem!

		GRBConstr **arcFlowCapacity = new GRBConstr *[n];
		for (i = 0; i < n; i++)
			arcFlowCapacity[i] = new GRBConstr[n];

		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
			{
				GRBLinExpr expr = 0.0;

				expr += wt[i][j] * y[i][j];

				for (k = 0; k < n; k++)
					expr += z[i][k][j] + z[k][j][i] - z[i][j][k];

				expr -= capacity * newCol[i][j] * lambda[0];

				arcFlowCapacity[i][j] = Master.addConstr(
					expr <= 0, "arcFlowCap_" + itos(i) + "_" + itos(j));
			}

		GRBConstr *lambdaSum = new GRBConstr[1];
		GRBLinExpr expr = 0.0;
		expr += lambda[0];
		lambdaSum[0] = Master.addConstr(expr == 1, "lambdaSum");

		if (PRINT_FOR_DEBUG)
			cout << "constraints added! \n";

		// set objective
		GRBLinExpr obj = 0.0;
		for (i = 0; i < n; i++)
			for (j = 1; j < n; j++)
			{
				obj += (price - cost) * dis_v[i][j] * wt[i][j] * y[i][j];
				obj -= cost * vw * dis_v[i][j] * newCol[i][j] * lambda[0];
			}

		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				for (k = 0; k < n; k++)
					obj -= (dis_v[i][k] + dis_v[k][j] - dis_v[i][j]) * z[i][j][k];

		Master.setObjective(obj, GRB_MAXIMIZE);

		int itrNum = 1;

		double convergePoint = -100;
		double PPobj = -100;

		bool findSolutionInThisIteration;

		double **sol = new double *[n];
		for (i = 0; i < n; i++)
			sol[i] = new double[n];

		while (true)
		{
			printf("======== DW-CG Iteration %d ========\n", itrNum);
			//*** the model is empty after I use model.relax()
			// GRBModel relax = Master.relax ();
			// relax.write ("dw-cg-relax.lp");
			// relax.optimize ();

			findSolutionInThisIteration = false;

			Master.write("dw-cg-master.lp");
			Master.getEnv().set(GRB_IntParam_OutputFlag, 0); // Silent Mode
			Master.optimize();

			if (PRINT_FOR_DEBUG)
				cout << "Master solved! \n";

			double objValue_master = Master.get(GRB_DoubleAttr_ObjVal);
			cout << "objective of Master problem: " << objValue_master << "\n";

			double pi[n][n];

			for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
					pi[i][j] = arcFlowCapacity[i][j].get(GRB_DoubleAttr_Pi);

			if (PRINT_FOR_DEBUG)
			{
				cout << "==> dual values:" << endl;
				for (i = 0; i < n; i++)
				{
					for (j = 0; j < n; j++)
						cout << pi[i][j] << ", ";
					cout << endl;
				}
			}

			double pi_lambda_sum = lambdaSum[0].get(GRB_DoubleAttr_Pi);

			// ===> set up the coefficient of x variables in objective function
			vector<vector<double>> xCoeff(n, std::vector<double>(n));
			for (i = 0; i < n - 1; i++)
				for (j = 1; j < n; j++)
				{
					xCoeff[i][j] = 0; //===> for min obj problem
					xCoeff[i][j] += cost * vw * dis_v[i][j];
					xCoeff[i][j] -= capacity * pi[i][j];
				}
			// assigen big value to unaccesiable arc in case it is used in runDominance()
			for (j = 0; j < n; j++)
				xCoeff[n - 1][j] = 1000000;
			for (i = 0; i < n; i++)
				xCoeff[i][0] = 1000000;
			// ===> end of setup

			if (PRINT_FOR_DEBUG)
				cout << "pi_lambda dual " << pi_lambda_sum << endl;

			if (testOneRouteObjAndDist)
			{
				if (itrNum > 13)
					exit(1);
				if (itrNum == 13)
				{
					// printf("pi(0,3)=%lf\n", pi[0][3]);
					// printf("pi(3,17)=%lf\n", pi[3][17]);
					// printf("pi(17,7)=%lf\n", pi[17][4]);
					// printf("pi(7,8)=%lf\n", pi[4][14]);
					// printf("pi(13,19)=%lf\n", pi[14][19]);
					// printf("d(0,3)=%lf\n", dis_v[0][3]);
					// printf("d(3,17)=%lf\n", dis_v[3][17]);
					// printf("d(17,7)=%lf\n", dis_v[17][7]);
					// printf("d(8,13)=%lf\n", dis_v[8][13]);
					// printf("d(13,19)=%lf\n", dis_v[13][19]);

					printf("xCoeff(0,1)=%lf\n", xCoeff[0][1]);
					printf("xCoeff(1,4)=%lf\n", xCoeff[1][4]);
					printf("xCoeffv(4,9)=%lf\n", xCoeff[4][9]);
					printf("xCoeff(0,7)=%lf\n", xCoeff[0][7]);
					printf("xCoeff(7,1)=%lf\n", xCoeff[7][1]);
					printf("d(0,1)=%lf\n", dis_v[0][1]);
					printf("d(1,4)=%lf\n", dis_v[1][4]);
					printf("d(4,9)=%lf\n", dis_v[4][9]);
					printf("d(0,7)=%lf\n", dis_v[0][7]);
					printf("d(7,1)=%lf\n", dis_v[7][1]);

					// cout << "pi_lambda dual " << pi_lambda_sum << endl;
					//===> print coefficients of x variables
					printf("{{%lf,", xCoeff[0][0]);
					for (int j = 1; j < n - 1; j++)
						printf("%lf,", xCoeff[0][j]);
					printf("%lf},\n", xCoeff[0][n - 1]);
					for (int i = 1; i < n - 1; i++)
					{
						printf("{%lf,", xCoeff[i][0]);
						for (int j = 1; j < n - 1; j++)
							printf("%lf,", xCoeff[i][j]);
						printf("%lf},\n", xCoeff[i][n - 1]);
					}
					printf("{%lf,", xCoeff[n - 1][0]);
					for (int j = 1; j < n - 1; j++)
						printf("%lf,", xCoeff[n - 1][j]);
					printf("%lf}};\n", xCoeff[n - 1][n - 1]);
					//===> end of printing coeff

					// double disGurobi = dis_v[0][17] + dis_v[17][15] + dis_v[15][11] + dis_v[11][9] + dis_v[9][5] + dis_v[5][13] + dis_v[13][19];
					// double objGurobi = -xCoeff[0][17] - xCoeff[17][15] - xCoeff[15][11] - xCoeff[11][9] - xCoeff[9][5] - xCoeff[5][13] - xCoeff[13][19] - pi_lambda_sum;

					// double disDominance = dis_v[0][1] + dis_v[1][4] + dis_v[4][9];
					// double objDominance = -xCoeff[0][1] - xCoeff[1][4] - xCoeff[4][9] - pi_lambda_sum;

					// cout << "objGurobi=" << objGurobi << endl;
					// cout << "disGurobi=" << disGurobi << endl;
					// cout << "obj=" << objDominance << endl;
					// cout << "dis=" << disDominance << endl;
				}
			}

			//================> Pricing Problem <================
			// double **sol = new double *[n];
			double objValue_PP;

			if (USE_DOMINACE_METHOD_OR_SUB_PROBLEM)
			{

				//===> the objective is to minimize value, so we get negative of the function used for Gurobi
				// // double **xCoeff = new double *[n];

				// vector<vector<double>> xCoeff(n, std::vector<double>(n));

				// for (i = 0; i < n - 1; i++)
				// 	for (j = 1; j < n; j++)
				// 	{
				// 		xCoeff[i][j] = 0;
				// 		xCoeff[i][j] += cost * vw * dis_v[i][j];
				// 		xCoeff[i][j] -= capacity * pi[i][j];
				// 	}
				// // assigen big value to unaccesiable arc in case it is used in runDominance()
				// for (j = 0; j < n; j++)
				// 	xCoeff[n - 1][j] = 1000000;
				// for (i = 0; i < n; i++)
				// 	xCoeff[i][0] = 1000000;

				//===> !!!!!! selectedRoute index starts from 0 !!!!!!
				//===> for example, selectedRoute[0]=0, selectedRoute[1]=6, selectedRoute[2]=9
				//===> it means the route is 1->7->10
				// return the cost of the selected Route !!!
				double disLimit = (double)route.DIS;
				// double objValue_PP;
				// vector<int> selectedRoute = runDominance(n, dis_v, xCoeff, disLimit, &objValue_PP);
				vector<int> selectedRoute;
				runDominance(n, dis_v, xCoeff, disLimit, &objValue_PP, selectedRoute);

				objValue_PP += pi_lambda_sum;
				objValue_PP = -objValue_PP;
				printf("objective of Pricing problem (dominace): %lf\n", objValue_PP);

				// if (objValue_PP < 0.000001 && objValue_PP > -0.000001)
				// {
				// 	cout << "Optimum Found!" << endl;

				// 	convergePoint = objValue_master;
				// 	PPobj = objValue_PP;

				// 	printf("selected reoute: \n");
				// 	for (i = 0; i < selectedRoute.size() - 1; i++)
				// 		printf("%d->%d ", selectedRoute[i], selectedRoute[i + 1]);
				// 	printf("\n");

				// 	break;
				// }
				if (selectedRoute.size() > 1)
				{
					findSolutionInThisIteration = true;

					//===> pre-assign sol array value
					for (i = 0; i < n; i++)
						for (j = 0; j < n; j++)
							sol[i][j] = 0;

					//===> assign value from returned route to sol array
					for (i = 0; i < selectedRoute.size() - 1; i++)
						sol[selectedRoute[i]][selectedRoute[i + 1]] = 1;

					if (PRINT_FOR_DEBUG)
					{
						cout << "selected reoute" << endl;
						for (i = 0; i < selectedRoute.size() - 1; i++)
							cout << selectedRoute[i] << "->" << selectedRoute[i + 1] << ", ";
						cout << endl;
					}
				}
				else
				{
					printf("ERROR: after dominance method, there are at least two nodes in selected route!\n");
					exit(1);
				}
			}
			else
			{

				//===> when not using dominace method, use Gurobi to solve the MIP

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
				obj += pi_lambda_sum;
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

				// double objValue_PP = PPmodel.get(GRB_DoubleAttr_ObjVal);
				objValue_PP = PPmodel.get(GRB_DoubleAttr_ObjVal);

				objValue_PP = -objValue_PP;
				cout << "objective of PP problem (gurobi): " << objValue_PP << "\n";

				// if (objValue_PP < 0.000001 && objValue_PP > -0.000001)
				// {
				// 	cout << "Optimum Found!" << endl;
				// 	convergePoint = objValue_master;
				// 	PPobj = objValue_PP;
				// 	break;
				// }

				if (PPmodel.get(GRB_IntAttr_SolCount) > 0)
				{

					findSolutionInThisIteration = true;

					// Extract solution and create new column

					for (i = 0; i < n; i++)
						sol[i] = PPmodel.get(GRB_DoubleAttr_X, x[i], n);
				}

				for (i = 0; i < n; i++)
					delete[] x[i];
				delete[] x;
				x = nullptr;
			}

			//======> AFTER SOLVING PRICING PROBLEM AND READ ROUTE INFO, ADD NEW VAR TO MASTER PROBLEM <======
			if (findSolutionInThisIteration)
			{

				GRBColumn col;

				// make a new column in arc flow capacity constraint
				for (i = 0; i < n; i++)
					for (j = 0; j < n; j++)
						if (sol[i][j] > 0.999 && sol[i][j] < 1.001)
							col.addTerm(-capacity, arcFlowCapacity[i][j]);

				// if (PRINT_FOR_DEBUG)
				{
					printf("selected reoute: \n");
					for (i = 0; i < n; i++)
						for (j = 0; j < n; j++)
							if (sol[i][j] > 0.999 && sol[i][j] < 1.001)
								printf("%d->%d ", i, j);
					printf("\n");
				}

				// make a new column in lamdaSum constraint
				col.addTerm(1, lambdaSum[0]);

				// add the new column into model
				double newLambdaCoeff = 0;

				for (i = 0; i < n; i++)
					for (j = 1; j < n; j++)
						newLambdaCoeff -= cost * vw * dis_v[i][j] * sol[i][j];

				cout << "newLambdaCoeff=" << newLambdaCoeff << endl;

				lambda.push_back(
					Master.addVar(0.0, GRB_INFINITY, newLambdaCoeff,
								  GRB_CONTINUOUS, col,
								  "lambda_" + itos(itrNum)));
			}
			else
			{
				printf("ERROR: no route solution in pricing problem is found!\n");
				exit(1);
			}

			if (objValue_PP < 0.000001 && objValue_PP > -0.000001)
			{
				cout << "Optimum Found!" << endl;

				convergePoint = objValue_master;
				PPobj = objValue_PP;

				break;
			}

			// Master.write("AddingVar.lp");

			// if (itrNum >= 3)
			// exit (1);

			itrNum++;
		}

		for (i = 0; i < n; i++)
			delete[] sol[i];
		delete[] sol;
		sol = nullptr;

		if (convergePoint < -99)
			cout << "No Result" << endl;
		else
		{
			cout << "========================" << endl;
			cout << "convergePoint=" << convergePoint << endl;
			cout << "PPobj=" << PPobj << endl;
			cout << "DW-CG itrNum=" << itrNum << endl;
		}

		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				y[i][j].set(GRB_CharAttr_VType, GRB_BINARY);

		// Fix variable other lambda var to be 0 except the last added lambda
		// lambda[0] represents route from 0 to (n-1)
		for (i = 0; i < itrNum - 1; i++)
			lambda[i].set(GRB_DoubleAttr_UB, 0);

		Master.optimize();

		double objValue = Master.get(GRB_DoubleAttr_ObjVal);
		cout << "objective with binary y and the last added route: " << objValue
			 << "\n";

		/*
		 if (Master.get (GRB_IntAttr_SolCount) > 0)
		 {

		 for (i = 0; i < itrNum ; i++)
		 {
		 double sol_val = lambda[i].get (GRB_DoubleAttr_X);
		 //if (sol[i] > 0.999 && sol[i] < 1.001)
		 cout << "lambda[" << i << "]=" << sol_val << endl;
		 }
		 }
		 */
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