// This is an example from https://groups.google.com/g/gurobi/c/pkBNfu-iX0k

// In this version, I added BPMP LP before starting column generation.
// Then use the routes from LP solution which satisfies distance limit and has no cycles
// as the initial solution for master problem.
// NOTE: RIGHT NOW IT ONLY WORKS FOR INSTANCES WHOSE LP OPTIMUM SOLUTION DOESN'T INCLUDE CYCLES
// when there is cycle, the code will probably falls into endless cycles
// THIS IS NOT COMPLETED WORK. SO ONLY TEST ON FOR EXAMPLE t10_05, t30_01 etc.
// for t30_01, the running time doesn't descrease :(
// so I decided to drop this idea
// BTW, for t30_01, I returned 2 best routes from dominance method, but it seems in master problem only the best route has positive corresponding value in master problem.
// so maybe we only need to return the best route from dominace.

// added preprocess and force some y, z variables zero in master problem, and some x variables zero in pricing problem
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
// #include "dominance.h"
// #include "dominance_inab.h"
// #include "dominance_inab_faster.h"
#include "dominance_inab_faster2.h"

// #include "dominance_inab_shortlabel.h"

// #include "dominance_inab_v0.h"
// #include "dominance_inab_v1.h"

// #include "dominance_doubleB2D.h"
// #include "dominance_doubleB2D_vector.h"

// #include "dominance_inaccessiblenode.h"
// #include "dominance_bit_hashmap.h"
// #include "dominance_bit_struct.h"

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
// bool USE_DOMINACE_METHOD_FOR_SUB_PROBLEM = false; // if it is false, Gurobi will be used for pricing problem

bool toRemove_addInitialSolutionInPricingProblem = false;

double bigM = 10000000;

// epsilon() indicates the precision with which numbers close to 1 can be represented.
// min() gives the lower bound of representable positive numbers
// epsilon() is typically much larger than min();
//  double epsilon = numeric_limits<double>::epsilon();
double min_double = numeric_limits<double>::min();

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
void findOneRouteFromCurrentNode(int &, int, vector<int> &, double &, double **, int, int, int);

int main(int argc, char *argv[])
{

	string pricingMethod;
	bool ADD_PREPROCESS = true;
	bool USE_LP_SOL_AS_INIT_SOL_IN_MASTER = true; // solve BPMP LP, and then pick up the elementary routes in solution as the initial routes in master problem

	//===> for dominance method solving pricing problem
	// if it is true, finish dominance process when the first maxNumNegativeRoutes are found
	// if it is false, then return at most maxNumNegativeRoutes good routes after dominance method is completed
	bool returnFirstFoundGoodRoutes = false;
	int maxNumNegativeRoutes = 2;

	//********* only read graph info and vehicle info by arguments
	//********* doesn't go through the graph info in the data folder
	if (argc != 5)
	{
		cout
			<< "Usage: ./bpmp_dw-cg.x nodesDataNameAndPath numberOfVehicles vehicleDataNameAndPath pricingMethod(dominance or gurobi)"
			<< endl;
		return 1;
	}

	if (argc == 5)
	{
		cout << "Nodes Data: " << argv[1] << endl;
		cout << "Number of Vehicles: " << argv[2] << endl;
		cout << "Vehicles Data: " << argv[3] << endl;
		numV = stoi(argv[2]);
		pricingMethod = argv[4];
		if (pricingMethod != "dominance" && pricingMethod != "gurobi")
		{
			cout << "ERROR: pricingMethod should be one of them: dominance, gurobi" << endl;
			exit(1);
		}
	}

	// clock_t beginTime = clock();

	// auto beginWallClock = high_resolution_clock::now();

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
	double disLimit = (double)route.DIS;

	// //==============reading vehicle info==============
	// string vehicleFileName = argv[3];

	// printf("reading vehicle file: %s ...\n", vehicleFileName.c_str());
	// // exit (1);
	// route.readSingleVehicleFile(vehicleFileName, numV);

	// // vehicle and origin always starts from 0
	// int vehicle[numV];
	// int origin[numV];

	// printf("vehicles (start from 0):\n");
	// for (int i = 0; i < numV; i++)
	// {
	// 	vehicle[i] = route.vehicle[i];
	// 	printf("%d ", vehicle[i]);
	// }
	// printf("\norigins (start from 0):\n");
	// for (int i = 0; i < numV; i++)
	// {
	// 	origin[i] = route.origin[i];
	// 	printf("%d ", origin[i]);
	// }
	// printf("\n");

	//================> end of reading data <================

	int startingNode = 0;
	int endingNode = n - 1;

	// check if there is demand associated with each node
	bool findNonUsedNode = false;
	for (int i = 0; i < n; i++)
	{
		double nodeDemand = 0;
		if (i != startingNode && i != endingNode)
		{
			for (int j = 0; j < n; j++)
			{
				if (j != startingNode)
					nodeDemand += wt[i][j];
				if (j != endingNode)
					nodeDemand += wt[j][i];
			}

			if (nodeDemand < 10 * min_double)
			{
				cout << "ERROR: node " << i << " has no delivery request." << endl;
				findNonUsedNode = true;
			}
		}
	}
	// if (findNonUsedNode)
	// 	exit(1);

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
		if (TIV.size() > 0)
		{
			cout << "ERROR: there are some triangle inequality violations in distance data!" << endl;
			exit(1);
		}
	}

	// // for a node i, when sum_(j in N) x(i,j)=0, then y(i,j)=0 for any j, and z(i,j,k), z(j,i,k), and z(j,k,i)=0 in optimal solution
	// // in column generation, we keep adding the newly found route to the RHS of the master problem,
	// // so after each iteration of CG, remove the newly added visited nodes from unvisitedNodes set;
	// this idea is wrong for finding upper bound. It forces y z be zero according to the result from all previous iteration results, but it might be no zero in next iterations.
	// the idea of setting up var=0 is to make sure that it is definitely zero in the final solution, so in brand and bound tree, it is always zero.
	// so that it is equilalents to adding some branching at root node of branch and bound tree using column generation
	// unordered_set<int> unvisitedNodes;

	//===> check if visiting some arcs over distance limit
	vector<vector<double>> inaccessibleArcs;
	vector<vector<int>> nodeNeighbors; // it is not used right now since it makes dominance run slower
	vector<vector<int>> nodeInaccNeighbors;
	if (ADD_PREPROCESS)
	{
		// initialize nodeNeighbors and nodeInaccNeighbors
		for (int i = 0; i < n; i++)
		{
			vector<int> neighbors;
			if (i != endingNode)
				for (int j = 0; j < n; j++)
				{
					if (j != startingNode)
						neighbors.push_back(j);
				}
			nodeNeighbors.push_back(neighbors);

			vector<int> neighbors2;
			nodeInaccNeighbors.push_back(neighbors2);

			// if (i != 0 && i != (n - 1))
			// 	unvisitedNodes.insert(i);
		}

		// find the arcs violate distance limit
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				if (i != startingNode && i != endingNode && j != endingNode && j != startingNode && i != j)
				{
					double distTemp = dis_v[startingNode][i] + dis_v[i][j] + dis_v[j][endingNode];
					if (distTemp > disLimit)
					{
						vector<double> arcTemp = {(double)i, (double)j, distTemp};
						inaccessibleArcs.push_back(arcTemp);

						// update nodeInaccNeighbors
						nodeInaccNeighbors[i].push_back(j);

						// update nodeNeighbors
						auto it = find(nodeNeighbors[i].begin(), nodeNeighbors[i].end(), j);

						if (it != nodeNeighbors[i].end())
						{
							nodeNeighbors[i].erase(it);
						}
					}
				}
			}
		cout << "# of inaccessibleArcs = " << inaccessibleArcs.size() << endl;

		// for (auto nodePairs : inaccessibleArcs)
		// {
		// 	for (auto node : nodePairs)
		// 		cout << node << " ";
		// 	cout << endl;
		// }

		// cout << "===> nodeInaccNeighbors:" << endl;
		// for (int i = 0; i < n; i++)
		// {
		// 	cout << "node " << i << endl;
		// 	for (auto &neighbor : nodeInaccNeighbors[i])
		// 		cout << neighbor << " ";
		// 	cout << endl;
		// }

		// cout << "===> nodeNeighbors:" << endl;
		// for (int i = 0; i < n; i++)
		// {
		// 	cout << "node " << i << endl;
		// 	for (auto &neighbor : nodeNeighbors[i])
		// 		cout << neighbor << " ";
		// 	cout << endl;
		// }
	}

	clock_t beginTime = clock();

	auto beginWallClock = high_resolution_clock::now();

	GRBEnv *env = 0;

	try
	{
		int i, j, k;

		// //****** newCol[n-1][j]=0 and newCol[i][0]=0
		// double **newCol = new double *[n];
		// for (i = 0; i < n; i++)
		// 	newCol[i] = new double[n];

		// for (i = 0; i < n; i++)
		// 	for (j = 0; j < n; j++)
		// 		newCol[i][j] = 0;

		// set up a new GRB environment
		env = new GRBEnv();

		//======> initial route solution <======
		vector<vector<int>> selectedRoutesVec;
		if (USE_GA_FOR_INIT_SOL)
		{
			// double disLimit = (double)route.DIS;

			// right now it only returns the selected route, no customers and profit info
			// returns only one routes as initial solution
			vector<int> finalRoute = runga(n, price, cost, capacity, disLimit, vw, wt, dis_v);

			cout << "========> selected route using genetic algorithm:" << endl;
			for (int i : finalRoute)
				cout << i << ",";
			cout << endl;

			selectedRoutesVec.push_back(finalRoute);

			// for (i = 0; i < finalRoute.size() - 1; i++)
			// {
			// 	newCol[finalRoute[i]][finalRoute[i + 1]] = 1;
			// 	// unvisitedNodes.erase(finalRoute[i]);
			// }
		}
		else if (USE_LP_SOL_AS_INIT_SOL_IN_MASTER)
		{
			cout << "start LP" << endl;
			//===> returns multiple routes as initial solution
			vector<vector<int>> finalRoutes;

			GRBModel BPMPlp = GRBModel(*env);

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
					y[i][j] = BPMPlp.addVar(0.0, 1.0, 0, GRB_CONTINUOUS,
											"y_" + itos(i) + "_" + itos(j));
					for (k = 0; k < n; k++)
					{
						string s = "z_" + itos(i) + "_" + itos(j) + "_" + itos(k);
						//*** in DW-CG, we force u<=capacity (capacity is the vehicle's capacity)
						z[i][j][k] = BPMPlp.addVar(0.0, capacity, 0.0, GRB_CONTINUOUS, s);
					}
				}

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
			if (ADD_PREPROCESS)
			{
				for (int i = 0; i < n; i++)
					for (j = 0; j < nodeInaccNeighbors[i].size(); j++)
					{
						vector<int> nbsTemp = nodeInaccNeighbors[i];
						// in BPMP model, we also set x[i][nbsTemp[j]]=0, but in CG's master problem, there is no x variables.
						// x is only in pricing problem.
						y[i][nbsTemp[j]].set(GRB_DoubleAttr_UB, 0);
						for (k = 0; k < n; k++)
						{
							z[i][nbsTemp[j]][k].set(GRB_DoubleAttr_UB, 0);
							z[i][k][nbsTemp[j]].set(GRB_DoubleAttr_UB, 0);
							// z[k][i][nbsTemp[j]].set(GRB_DoubleAttr_UB, 0);
							z[k][nbsTemp[j]][i].set(GRB_DoubleAttr_UB, 0);
						}
					}

				for (int i = 0; i < n; i++)
					if (i != startingNode && i != endingNode)
						for (auto &k : nodeNeighbors[i])
							for (auto &j : nodeNeighbors[k])
								if (i != j)
									if (dis_v[startingNode][i] + dis_v[i][k] + dis_v[k][j] + dis_v[j][endingNode] > disLimit)
										z[i][j][k].set(GRB_DoubleAttr_UB, 0);

				// cout << "numViolatedTriples=" << numViolatedTriples << endl;
			}

			GRBVar **x = nullptr;
			GRBVar s[n];
			GRBVar theta[n][n];

			x = new GRBVar *[n];
			for (i = 0; i < n; i++)
				x[i] = new GRBVar[n];

			for (i = 0; i < n; i++)
			{
				s[i] = BPMPlp.addVar(0.0, n, 0.0, GRB_CONTINUOUS,
									 "s_" + itos(i));
				for (j = 0; j < n; j++)
				{
					x[i][j] = BPMPlp.addVar(0.0, 1.0, 0, GRB_CONTINUOUS,
											"x_" + itos(i) + "_" + itos(j));
					theta[i][j] = BPMPlp.addVar(
						0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
						"theta_" + itos(i) + "_" + itos(j));
				}
			}

			for (i = 0; i < n; i++)
			{
				x[i][i].set(GRB_DoubleAttr_UB, 0);
				x[i][0].set(GRB_DoubleAttr_UB, 0);
				x[n - 1][i].set(GRB_DoubleAttr_UB, 0);
				theta[i][i].set(GRB_DoubleAttr_UB, 0);
				theta[i][0].set(GRB_DoubleAttr_UB, 0);
				theta[n - 1][i].set(GRB_DoubleAttr_UB, 0);
			}
			if (ADD_PREPROCESS)
			{
				for (i = 0; i < n; i++)
				{
					vector<int> nbsTemp = nodeInaccNeighbors[i];
					for (j = 0; j < nbsTemp.size(); j++)
					{
						x[i][nbsTemp[j]].set(GRB_DoubleAttr_UB, 0);
					}
				}
			}

			// vehicle goes out of node 1
			GRBLinExpr expr1 = 0.0;
			for (i = 1; i < n; i++)
				expr1 += x[0][i];
			BPMPlp.addConstr(expr1 == 1, "origin");

			// vehicle goes back to node n
			GRBLinExpr expr2 = 0.0;
			for (i = 0; i < n - 1; i++)
				expr2 += x[i][n - 1];
			BPMPlp.addConstr(expr2 == 1, "destination");

			// flow conservation
			for (int k = 1; k < n - 1; k++)
			{
				GRBLinExpr expr = 0;
				for (i = 0; i < n - 1; i++)
					expr += x[i][k];
				for (j = 1; j < n; j++)
					expr -= x[k][j];
				BPMPlp.addConstr(expr == 0, "flow_conservation_" + itos(k));
			}

			// distance
			GRBLinExpr expr3 = 0.0;
			for (i = 0; i < n - 1; i++)
				for (j = 1; j < n; j++)
					expr3 += dis_v[i][j] * x[i][j];
			BPMPlp.addConstr(expr3 <= disLimit, "distance");

			// node degree less than 1
			for (int k = 1; k < n - 1; k++)
			{
				GRBLinExpr expr = 0.0;
				for (i = 0; i < n - 1; i++)
					expr += x[i][k];
				BPMPlp.addConstr(expr <= 1, "indegree_" + itos(k));
			}

			// subtour elimination
			for (i = 0; i < n - 1; i++)
				for (j = 1; j < n; j++)
				{
					GRBLinExpr expr = 0.0;
					expr += s[i] - s[j] + (n - 1) * x[i][j] + (n - 3) * x[j][i];
					BPMPlp.addConstr(expr <= n - 2,
									 "s_" + itos(i) + "_" + itos(j));
				}

			// arc flow
			for (i = 0; i < n - 1; i++)
				for (j = 1; j < n; j++)
				{
					GRBLinExpr expr = 0.0;
					expr += wt[i][j] * y[i][j] - theta[i][j];
					for (k = 0; k < n; k++)
						expr += z[i][k][j] + z[k][j][i] - z[i][j][k];

					BPMPlp.addConstr(expr == 0,
									 "flow_" + itos(i) + "_" + itos(j));
				}

			// arc flow upperbound
			for (i = 0; i < n - 1; i++)
				for (j = 1; j < n; j++)
				{
					GRBLinExpr expr = 0.0;
					expr += theta[i][j] - capacity * x[i][j];
					BPMPlp.addConstr(expr <= 0,
									 "flowBound_" + itos(i) + "_" + itos(j));
				}

			// set objective
			GRBLinExpr obj = 0.0;
			for (i = 0; i < n - 1; i++)
				for (j = 1; j < n; j++)
				{
					obj += price * dis_v[i][j] * wt[i][j] * y[i][j];
					obj -= cost * dis_v[i][j] * theta[i][j];
					obj -= cost * vw * dis_v[i][j] * x[i][j];
				}
			BPMPlp.setObjective(obj, GRB_MAXIMIZE);

			// Silent Mode
			// BPMPlp.getEnv().set(GRB_IntParam_OutputFlag, 0);
			BPMPlp.optimize();

			// Extract solution
			if (BPMPlp.get(GRB_IntAttr_SolCount) > 0)
			{
				double **sol = new double *[n];

				cout << "Selected requests: " << endl;
				for (i = 0; i < n; i++)
				{
					sol[i] = BPMPlp.get(GRB_DoubleAttr_X, y[i], n);
					for (j = 0; j < n; j++)
						if (sol[i][j] > min_double)
							printf("%d -- %d\n", i + 1, j + 1);
				}

				cout << "Selected arcs: " << endl;
				for (i = 0; i < n; i++)
				{
					sol[i] = BPMPlp.get(GRB_DoubleAttr_X, x[i], n);
					for (int j = 0; j < n; j++)
						if (sol[i][j] > min_double)
							printf("%d -- %d %lf\n", i + 1, j + 1, sol[i][j]);
				}

				int nextNode = -1;
				for (j = 0; j < n; j++)
				{
					if (j != startingNode && sol[startingNode][j] > min_double)
					{
						nextNode = j;
						break;
					}
				}
				cout << "nextNode = " << nextNode << endl;
				while (nextNode > 0)
				{
					cout << "------start a new route------" << endl;
					double routeFlow = bigM;
					double distance = 0;
					int count = 1;
					vector<int> oneRouteNodes;
					oneRouteNodes.push_back(startingNode);
					oneRouteNodes.push_back(nextNode);
					if (sol[startingNode][nextNode] < routeFlow)
						routeFlow = sol[startingNode][nextNode];

					cout << "routeFlow=" << routeFlow << endl;
					if (nextNode != endingNode)
						findOneRouteFromCurrentNode(count, nextNode, oneRouteNodes, routeFlow, sol, n, startingNode, endingNode);
					else
						sol[startingNode][endingNode] -= routeFlow;

					// check if there is cycle in this route first.
					// if no cycle, add to selectedRoutesVec
					//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					// work here !!!!!!!!!!!!!!!!!!!!!!!!;

					for (i = 0; i < oneRouteNodes.size() - 1; i++)
						distance += dis_v[oneRouteNodes[i]][oneRouteNodes[i + 1]];

					cout << "===> one found route:" << endl;
					for (i = 0; i < oneRouteNodes.size(); i++)
						cout << oneRouteNodes[i] << " ";
					cout << endl;
					cout << "route flow: " << routeFlow << endl;
					cout << "distance: " << distance << endl;

					//===> BE CAREFUL: ONLY ADD THE ROUTE WITHIN DISTANCE LIMIT
					// SOME ROUTES MIGHT HAVE OVER LIMIT DISTANCE
					// IF WE ADD IT TO MASTER PROBLEM, IT MIGHT CAUSE COLUMN GENERATION COVERGE POINT > LP OPTIMAL SOLUTION
					if (distance <= disLimit)
						selectedRoutesVec.push_back(oneRouteNodes);

					// remove the route flow from sol[][]
					for (i = 0; i < oneRouteNodes.size() - 1; i++)
						sol[oneRouteNodes[i]][oneRouteNodes[i + 1]] -= routeFlow;

					nextNode = -1;
					for (j = 0; j < n; j++)
					{
						if (j != startingNode && sol[startingNode][j] > min_double)
						{
							nextNode = j;
							break;
						}
					}
				}
				reportTime(beginTime, beginWallClock);

				for (i = 0; i < n; i++)
					delete[] sol[i];
				delete[] sol;

				double objValue = BPMPlp.get(GRB_DoubleAttr_ObjVal);
				cout << "obj value: " << objValue << endl;

				reportTime(beginTime, beginWallClock);
			}

			for (i = 0; i < n; i++)
			{
				delete[] x[i];
				delete[] y[i];
				for (int j = 0; j < n; j++)
					delete[] z[i][j];
				delete[] z[i];
			}
			delete[] x;
			delete[] y;
			delete[] z;

			BPMPlp.reset();
		}
		else
		{
			// newCol[0][n - 1] = 1;
			vector<int> oneRouteNodes;
			oneRouteNodes.push_back(startingNode);
			oneRouteNodes.push_back(endingNode);
			selectedRoutesVec.push_back(oneRouteNodes);
		}
		//======> end of initial route selection <======

		GRBModel Master = GRBModel(*env);

		// Turn off presolve
		// Master.set(GRB_IntParam_Presolve, 0);
		// Silent Mode
		Master.getEnv().set(GRB_IntParam_OutputFlag, 0);

		//****** Variables of master problem
		//*** set up var lambda
		vector<GRBVar> lambda;
		// lambda.push_back(
		// 	Master.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
		// 				  "lambda_" + itos(0)));
		// if (PRINT_FOR_DEBUG)
		// 	cout << "variable lambda is added \n";

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

		if (ADD_PREPROCESS)
		{
			for (int i = 0; i < n; i++)
				for (j = 0; j < nodeInaccNeighbors[i].size(); j++)
				{
					vector<int> nbsTemp = nodeInaccNeighbors[i];
					// in BPMP model, we also set x[i][nbsTemp[j]]=0, but in CG's master problem, there is no x variables.
					// x is only in pricing problem.
					y[i][nbsTemp[j]].set(GRB_DoubleAttr_UB, 0);
					for (k = 0; k < n; k++)
					{
						z[i][nbsTemp[j]][k].set(GRB_DoubleAttr_UB, 0);
						z[i][k][nbsTemp[j]].set(GRB_DoubleAttr_UB, 0);
						// z[k][i][nbsTemp[j]].set(GRB_DoubleAttr_UB, 0);
						z[k][nbsTemp[j]][i].set(GRB_DoubleAttr_UB, 0);
					}
				}

			int numViolatedTriples = 0;
			for (int i = 0; i < n; i++)
				if (i != startingNode && i != endingNode)
					for (auto &k : nodeNeighbors[i])
						for (auto &j : nodeNeighbors[k])
							if (i != j)
								if (dis_v[startingNode][i] + dis_v[i][k] + dis_v[k][j] + dis_v[j][endingNode] > disLimit)
								{
									z[i][j][k].set(GRB_DoubleAttr_UB, 0);
									numViolatedTriples++;
								}
			cout << "numViolatedTriples=" << numViolatedTriples << endl;
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

				// expr -= capacity * newCol[i][j] * lambda[0];

				arcFlowCapacity[i][j] = Master.addConstr(
					expr <= 0, "arcFlowCap_" + itos(i) + "_" + itos(j));
			}

		GRBConstr *lambdaSum = new GRBConstr[1];
		GRBLinExpr expr = 0.0;
		// expr += lambda[0];
		lambdaSum[0] = Master.addConstr(expr == 1, "lambdaSum");

		// set objective
		GRBLinExpr obj = 0.0;
		for (i = 0; i < n; i++)
			for (j = 1; j < n; j++)
			{
				obj += (price - cost) * dis_v[i][j] * wt[i][j] * y[i][j];
				// obj -= cost * vw * dis_v[i][j] * newCol[i][j] * lambda[0];
			}

		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				for (k = 0; k < n; k++)
					obj -= (dis_v[i][k] + dis_v[k][j] - dis_v[i][j]) * z[i][j][k];

		Master.setObjective(obj, GRB_MAXIMIZE);

		//===> add pre-selected routes from 1 GA, 2 BPMP LP, 3 startingNode->endingNode to the master problem's constraints and obj
		{
			int count = 1;
			for (auto &selectedRoute : selectedRoutesVec)
			{
				if (selectedRoute.size() > 1)
				{
					// if (PRINT_FOR_DEBUG)
					{
						cout << "======> selected reoute" << endl;
						for (i = 0; i < selectedRoute.size() - 1; i++)
							cout << selectedRoute[i] << "->" << selectedRoute[i + 1] << ", ";
						cout << endl;
					}

					// cout << "======> add initial routes from BPMP LP to master model" << endl;

					GRBColumn col;

					for (i = 0; i < selectedRoute.size() - 1; i++)
						col.addTerm(-capacity, arcFlowCapacity[selectedRoute[i]][selectedRoute[i + 1]]);

					// make a new column in lamdaSum constraint
					col.addTerm(1, lambdaSum[0]);

					// add the new column into model
					double newLambdaCoeff = 0;

					for (i = 0; i < selectedRoute.size() - 1; i++)
						newLambdaCoeff -= cost * vw * dis_v[selectedRoute[i]][selectedRoute[i + 1]];

					cout << "newLambdaCoeff=" << newLambdaCoeff << endl;

					lambda.push_back(
						Master.addVar(0.0, GRB_INFINITY, newLambdaCoeff, GRB_CONTINUOUS, col, "lambda_itr" + itos(0) + "_" + itos(count)));
					count++;
				}
				else
				{
					printf("ERROR: after dominance method, there are at least two nodes in selected route!\n");
					exit(1);
				}
			}
		}
		if (PRINT_FOR_DEBUG)
			cout << "constraints added! \n";

		int itrNum = 1;

		double convergePoint = -100;
		double PPobj = -100;

		bool findSolutionInThisIteration;

		vector<vector<vector<double>>> solNew;

		while (true)
		{
			printf("======== DW-CG Iteration %d ========\n", itrNum);
			//*** the model is empty after I use model.relax()
			// GRBModel relax = Master.relax ();
			// relax.write ("dw-cg-relax.lp");
			// relax.optimize ();

			findSolutionInThisIteration = false;

			Master.write("dw-cg-master.lp");

			Master.optimize();

			if (PRINT_FOR_DEBUG)
				cout << "Master solved! \n";

			double objValue_master = Master.get(GRB_DoubleAttr_ObjVal);
			cout << "objective of Master problem: " << objValue_master << "\n";

			// Extract solution and create new column
			double **solTemp = new double *[n];
			for (i = 0; i < n; i++)
				solTemp[i] = new double[n];
			for (i = 0; i < n; i++)
			{
				solTemp[i] = Master.get(GRB_DoubleAttr_X, y[i], n);
				for (j = 0; j < n; j++)
				{
					if (solTemp[i][j] > min_double)
					{
						cout << "y[" << i << "," << j << "]=" << solTemp[i][j] << ", ";
						cout << "w = " << wt[i][j] << ", ";
						cout << "wy=" << wt[i][j] * solTemp[i][j] << endl;
					}

					for (int k = 0; k < n; k++)
					{
						if (z[i][j][k].get(GRB_DoubleAttr_X) > min_double)
							cout << "z_" << i << " " << j << " " << k << " = " << z[i][j][k].get(GRB_DoubleAttr_X);
					}
				}
			}

			for (i = 0; i < lambda.size(); i++)
			{
				double temp = lambda[i].get(GRB_DoubleAttr_X);
				cout << "lambda" << i + 1 << "=" << temp << endl;
			}

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

			//================> Pricing Problem <================
			vector<double> objValue_PP_vec; // stores the objValue of each route, the cost is asending, such as, objValue_PP[0]=-1.9, objValue_PP[1] =-1.6, objValue_PP[2]=-1.2...

			if (pricingMethod == "dominance")
			{
				//===> !!!!!! selectedRoute index starts from 0 !!!!!!
				//===> for example, selectedRoute[0]=0, selectedRoute[1]=6, selectedRoute[2]=9
				//===> it means the route is 1->7->10
				// return the cost of the selected Route !!!
				vector<vector<int>> selectedRoutesVec;

				// runDominance(n, dis_v, xCoeff, disLimit, &objValue_PP, selectedRoute, startingNode, endingNode, nodeNeighbors);
				// runDominance(n, dis_v, xCoeff, disLimit, &objValue_PP, selectedRoute, startingNode, endingNode);
				// right now these arguments is only for dominance_inab_faster2.h. for other dominance variants, use the above function
				runDominance(n, dis_v, xCoeff, disLimit, objValue_PP_vec, selectedRoutesVec, returnFirstFoundGoodRoutes, maxNumNegativeRoutes, startingNode, endingNode);

				if (objValue_PP_vec.size() == 0)
				{
					cout << "ERROR: objValue_PP_vec should have at least one element." << endl;
					exit(1);
				}
				else
				{
					//===> pricing model of dominance doesn't include the constant. Add it now
					for (i = 0; i < objValue_PP_vec.size(); i++)
					{
						double objValue_PP = objValue_PP_vec[i];
						objValue_PP += pi_lambda_sum;
						objValue_PP = -objValue_PP;
						objValue_PP_vec[i] = objValue_PP;
					}
				}

				printf("best found objective of Pricing problem (dominace): %lf\n", objValue_PP_vec[0]);

				solNew.clear();
				for (auto &selectedRoute : selectedRoutesVec)
				{
					if (selectedRoute.size() > 1)
					{
						findSolutionInThisIteration = true;

						vector<vector<double>> oneRouteSol(n, vector<double>(n, 0)); // initialize nxn matrix with 0
						//===> assign value from returned route to sol array
						for (i = 0; i < selectedRoute.size() - 1; i++)
							oneRouteSol[selectedRoute[i]][selectedRoute[i + 1]] = 1;

						solNew.push_back(oneRouteSol);

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
			}
			// else
			// right now Gurobi returns only 1 route
			else if (pricingMethod == "gurobi")
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
				if (ADD_PREPROCESS)
				{
					for (i = 0; i < n; i++)
					{
						vector<int> nbsTemp = nodeInaccNeighbors[i];
						for (j = 0; j < nbsTemp.size(); j++)
						{
							x[i][nbsTemp[j]].set(GRB_DoubleAttr_UB, 0);
						}
					}
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
				PPmodel.set(GRB_DoubleParam_OptimalityTol, 1e-9); // Optimality tolerance (constraints)
				PPmodel.set(GRB_DoubleParam_MIPGap, 0);			  // set up mip gap for obj
				PPmodel.set(GRB_IntParam_Presolve, 0);			  // Turn off presolve

				PPmodel.write("DW_pricingModel.lp");
				PPmodel.optimize();

				double objValue_PP = PPmodel.get(GRB_DoubleAttr_ObjVal);

				objValue_PP = -objValue_PP;
				objValue_PP_vec.push_back(objValue_PP);
				cout << "objective of PP problem (gurobi): " << objValue_PP << "\n";

				if (PPmodel.get(GRB_IntAttr_SolCount) > 0)
				{

					findSolutionInThisIteration = true;

					vector<vector<double>> oneRouteSol(n, vector<double>(n, 0)); // initialize nxn matrix with 0

					// Extract solution and create new column
					double **sol = new double *[n];
					for (i = 0; i < n; i++)
						sol[i] = new double[n];
					for (i = 0; i < n; i++)
					{
						sol[i] = PPmodel.get(GRB_DoubleAttr_X, x[i], n);
						for (j = 0; j < n; j++)
							oneRouteSol[i][j] = sol[i][j];
					}

					solNew.push_back(oneRouteSol);
				}

				for (i = 0; i < n; i++)
					delete[] x[i];
				delete[] x;
				x = nullptr;
			}

			//======> AFTER SOLVING PRICING PROBLEM AND READ ROUTE INFO, ADD NEW VAR TO MASTER PROBLEM <======
			double objValue_PP = objValue_PP_vec[0]; // the best obj found
			if (objValue_PP < 0.000001 && objValue_PP > -0.000001)
			{
				cout << "Optimum Found!" << endl;

				convergePoint = objValue_master;
				PPobj = objValue_PP;

				break;
			}

			if (findSolutionInThisIteration)
			{
				cout << "======> assigning routes to master model" << endl;
				int count = 1;
				for (auto sol : solNew)
				{
					GRBColumn col;

					// make a new column in arc flow capacity constraint
					for (i = 0; i < n; i++)
						for (j = 0; j < n; j++)
							if (sol[i][j] > 0.999 && sol[i][j] < 1.001)
							{
								col.addTerm(-capacity, arcFlowCapacity[i][j]);
								// unvisitedNodes.erase(i);
							}
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
						Master.addVar(0.0, GRB_INFINITY, newLambdaCoeff, GRB_CONTINUOUS, col, "lambda_itr" + itos(itrNum) + "_" + itos(count)));
					count++;
				}

				// // --------if a node is not visited, then corresponding y and z are zeros----------
				// for (auto &i : unvisitedNodes)
				// {
				// 	for (j = 0; j < n; j++)
				// 	{
				// 		y[i][j].set(GRB_DoubleAttr_UB, 0);
				// 		y[j][i].set(GRB_DoubleAttr_UB, 0);
				// 		for (k = 0; k < n; k++)
				// 		{
				// 			z[i][j][k].set(GRB_DoubleAttr_UB, 0);
				// 			z[j][k][i].set(GRB_DoubleAttr_UB, 0);
				// 			z[j][i][k].set(GRB_DoubleAttr_UB, 0);
				// 		}
				// 	}
				// }
			}
			else
			{
				printf("ERROR: no route solution in pricing problem is found!\n");
				exit(1);
			}

			// Master.write("AddingVar.lp");

			// if (itrNum >= 3)
			// exit (1);

			itrNum++;
		}

		// for (i = 0; i < n; i++)
		// 	delete[] sol[i];
		// delete[] sol;
		// sol = nullptr;

		if (convergePoint < -99)
			cout << "No Result" << endl;
		else
		{
			cout << "========================" << endl;
			cout << "convergePoint=" << convergePoint << endl;
			cout << "PPobj=" << PPobj << endl;
			cout << "DW-CG itrNum=" << itrNum << endl;
		}

		reportTime(beginTime, beginWallClock);

		// use the last added route to solve BPMP to get lower bound

		if (false)
		{
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

	// reportTime(beginTime, beginWallClock);

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

void findOneRouteFromCurrentNode(int &count, int currentNode, vector<int> &oneRouteNodes, double &routeFlow, double **sol, int n, int startingNode, int endingNode)
{
	cout << "--- count " << count << "---" << endl;
	cout << "currentNode=" << currentNode << endl;
	int nextNode;
	bool findNextNode = false;
	for (int j = 0; j < n; j++)
	{
		if (j != startingNode && sol[currentNode][j] > min_double)
		{
			nextNode = j;
			findNextNode = true;
			count++;
			break;
		}
	}

	if (findNextNode == false)
	{
		cout << "ERROR: can't find next node in BPMP LP solution route construction." << endl;
		exit(1);
	}
	oneRouteNodes.push_back(nextNode);
	if (sol[currentNode][nextNode] < routeFlow)
		routeFlow = sol[currentNode][nextNode];

	cout << "nextNode=" << nextNode << endl;
	cout << "routeFlow=" << routeFlow << endl;

	if (count > 100)
	{
		cout << "ERROR: repeat too many times in findOneRouteFromCurrentNode() function" << endl;
		exit(1);
	}

	if (nextNode == endingNode)
		return;
	else
		findOneRouteFromCurrentNode(count, nextNode, oneRouteNodes, routeFlow, sol, n, startingNode, endingNode);
}
