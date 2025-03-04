// added preprocess and force some y, z variables zero in master problem, and some x variables zero in pricing problem

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
#include "readData.h"

#include <bits/stdc++.h> //to use unordered_set, which is hash set

using namespace std;
using namespace std::chrono;

// if need to turn off callback function as well, please go to main and remove (or comment) the following commands:
// (1) model.set(GRB_IntParam_LazyConstraints, 1);
// (2) printIntSol cb = printIntSol(x, y, u, ratioSet, n);
// (3) model.setCallback(&cb);
// on the contrary, if need to use ratio cut, then the above 3 commands need to be active
// it has to be done manually since if we use (2) and (3) in if(USE_RATIO_CUT), when it is true, the program has problems as I observed.

bool USE_RATIO_CUT = false;

bool USE_PROFITABLE_2_NODES_CYCLES = false;

bool PRINT_VAR_VALUE = true;

bool ADD_PREPROCESS = true;

#define EPSILON 0.00001
double capacity;

string
itos(int i)
{
	stringstream s;
	s << i;
	return s.str();
}

class printIntSol : public GRBCallback
{
public:
	GRBVar **xv;
	GRBVar **yv;
	GRBVar ***uv;
	unordered_set<string> ratioSet;
	int n;
	printIntSol(GRBVar **xvars, GRBVar **yvars, GRBVar ***uvars,
				unordered_set<string> passedSet, int nvar)
	{
		xv = xvars;
		yv = yvars;
		uv = uvars;
		ratioSet = passedSet;
		n = nvar;
	}

protected:
	void
	callback()
	{
		try
		{
			if (where == GRB_CB_MIPSOL)
			{

				int i, j, k;
				double **x = new double *[n];
				double **y = new double *[n];
				double ***u = new double **[n];

				for (i = 0; i < n; i++)
				{
					u[i] = new double *[n];

					for (j = 0; j < n; j++)
						u[i][j] = new double[n];
				}

				// printf("x vars: ");
				for (i = 0; i < n; i++)
					x[i] = getSolution(xv[i], n);

				for (i = 0; i < n; i++)
					y[i] = getSolution(yv[i], n);

				for (i = 0; i < n; i++)
					for (j = 0; j < n; j++)
						u[i][j] = getSolution(uv[i][j], n);

				if (PRINT_VAR_VALUE)
				{
					// Found an integer feasible solution
					printf("\n======> an integer solution found.\n");

					printf("x vars: ");
					for (i = 0; i < n; i++)
					{
						for (j = 0; j < n; j++)
						{
							if (x[i][j] > 0.5)
								printf("%d->%d ", i, j);
						}
					}
					cout << endl;

					printf("y vars: ");
					for (i = 0; i < n; i++)
					{
						for (j = 0; j < n; j++)
						{
							if (y[i][j] > 0.5)
								printf("%d->%d ", i, j);
						}
					}
					cout << endl;
					// printf("u vars: ");
					// for (i = 0; i < n; i++)
					// 	for (j = 0; j < n; j++)
					// 		for (k = 0; k < n; k++)
					// 			if (u[i][j][k] > EPSILON)
					// 				printf("u[%d,%d,%d] = %lf ", i, j, k, u[i][j][k]);
					// printf("\n");
				}

				// if (USE_RATIO_CUT)
				// {
				// 	for (i = 0; i < n; i++)
				// 		for (j = 0; j < n; j++)
				// 			for (k = 0; k < n; k++)
				// 			{
				// 				string key = itos(i) + "," + itos(j) + "," + itos(k);

				// 				if (u[i][j][k] > EPSILON && y[i][j] > 0.5 && (ratioSet.find(key) != ratioSet.end()))
				// 				{
				// 					// cout << key << " found" << endl;
				// 					// exit (1);

				// 					// u and y can't be positive at the same time
				// 					// u(ijk)<=(1-y(ij))*Q

				// 					GRBLinExpr expr = 0.0;
				// 					expr += uv[i][j][k] - (1 - yv[i][j]) * capacity;

				// 					addLazy(expr <= 0);
				// 					printf(
				// 						"*** price cost ratio cut is added for (%d %d %d)***\n",
				// 						i, j, k);
				// 				}
				// 			}
				// }

				for (i = 0; i < n; i++)
				{
					delete[] x[i];
					delete[] y[i];
					for (j = 0; j < n; j++)
						delete[] u[i][j];
					delete[] u[i];
				}
				delete[] x;
				delete[] y;
				delete[] u;
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

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		cout << "Usage: ./BPMP.x dataFolder" << endl;
		return 1;
	}

	cout << argc << endl;
	for (int i = 1; i < argc; ++i)
		cout << argv[i] << endl;

	readData route;
	auto it = filesystem::directory_iterator(argv[1]);
	for (const auto &entry : it)
	{
		filesystem::path path = entry.path();
		if (entry.is_regular_file())
		{
			string filename = path.string();
			route.readSingleFile(filename);

			// route.printStats();
			// route.printData();
			cout << "reading file: " << filename << endl;
			;

			int n = route.numOfNode;
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
			capacity = Q;

			// start calling Gurobi
			int i;

			GRBEnv *env = NULL;
			GRBVar **x = NULL;
			GRBVar **y = NULL;

			GRBVar s[n];
			// GRBVar u[n][n][n];
			GRBVar ***u = NULL;
			GRBVar theta[n][n];

			x = new GRBVar *[n];
			y = new GRBVar *[n];
			u = new GRBVar **[n];
			for (i = 0; i < n; i++)
			{
				x[i] = new GRBVar[n];
				y[i] = new GRBVar[n];
				u[i] = new GRBVar *[n];
				for (int j = 0; j < n; j++)
					u[i][j] = new GRBVar[n];
			}

			//===> add preprocessing
			//===> check if visiting some arcs over distance limit
			vector<vector<double>> inaccessibleArcs;
			vector<vector<int>> nodeNeighbors; // it is not used right now since it makes dominance run slower
			vector<vector<int>> nodeInaccNeighbors;
			if (ADD_PREPROCESS)
			{
				int startingNode = 0;
				int endingNode = n - 1;

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
				}

				// find the arcs violate distance limit
				for (int i = 0; i < n; i++)
					for (int j = 0; j < n; j++)
					{
						if (i != startingNode && i != endingNode && j != endingNode && j != startingNode && i != j)
						{
							double distTemp = dis[startingNode][i] + dis[i][j] + dis[j][endingNode];
							if (distTemp > (double)route.DIS)
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
			auto start = std::chrono::high_resolution_clock::now();

			try
			{
				int j, k;

				env = new GRBEnv();
				GRBModel model = GRBModel(*env);

				// Create binary decision variables
				for (i = 0; i < n; i++)
				{
					s[i] = model.addVar(0.0, n, 0.0, GRB_CONTINUOUS,
										"s_" + itos(i));
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

				for (i = 0; i < n; i++)
				{
					x[i][i].set(GRB_DoubleAttr_UB, 0);
					y[i][i].set(GRB_DoubleAttr_UB, 0);
					x[i][0].set(GRB_DoubleAttr_UB, 0);
					y[i][0].set(GRB_DoubleAttr_UB, 0);
					x[n - 1][i].set(GRB_DoubleAttr_UB, 0);
					y[n - 1][i].set(GRB_DoubleAttr_UB, 0);
				}

				if (ADD_PREPROCESS)
				{
					//===> preprocess 1
					for (i = 0; i < n; i++)
					{
						vector<int> nbsTemp = nodeInaccNeighbors[i];
						for (j = 0; j < nodeInaccNeighbors[i].size(); j++)
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

					//===> preprocess 2
					int numViolatedTriples = 0;
					int startingNode = 0;
					int endingNode = n - 1;
					for (int i = 0; i < n; i++)
						// if (i != startingNode && i != endingNode)
						if (i != 0 && i != n - 1)
							for (auto &k : nodeNeighbors[i])
								for (auto &j : nodeNeighbors[k])
									if (i != j)
										if (dis[startingNode][i] + dis[i][k] + dis[k][j] + dis[j][endingNode] > (double)route.DIS)
										{
											u[i][j][k].set(GRB_DoubleAttr_UB, 0);
											numViolatedTriples++;
										}
					cout << "numViolatedTriples=" << numViolatedTriples << endl;
				}

				// vehicle goes out of node 1
				GRBLinExpr expr1 = 0.0;
				for (i = 1; i < n; i++)
					expr1 += x[0][i];
				model.addConstr(expr1 == 1, "origin");

				// vehicle goes back to node n
				GRBLinExpr expr2 = 0.0;
				for (i = 0; i < n - 1; i++)
					expr2 += x[i][n - 1];
				model.addConstr(expr2 == 1, "destination");

				// flow conservation
				for (int k = 1; k < n - 1; k++)
				{
					GRBLinExpr expr = 0;
					for (i = 0; i < n - 1; i++)
						expr += x[i][k];
					for (j = 1; j < n; j++)
						expr -= x[k][j];
					model.addConstr(expr == 0, "flow_conservation_" + itos(k));
				}

				// distance
				GRBLinExpr expr3 = 0.0;
				for (i = 0; i < n - 1; i++)
					for (j = 1; j < n; j++)
						expr3 += dis[i][j] * x[i][j];
				model.addConstr(expr3 <= route.DIS, "distance");

				// node degree less than 1
				for (int k = 1; k < n - 1; k++)
				{
					GRBLinExpr expr = 0.0;
					for (i = 0; i < n - 1; i++)
						expr += x[i][k];
					model.addConstr(expr <= 1, "indegree_" + itos(k));
				}

				// subtour elimination
				for (i = 0; i < n - 1; i++)
					for (j = 1; j < n; j++)
					{
						GRBLinExpr expr = 0.0;
						expr += s[i] - s[j] + (n - 1) * x[i][j] + (n - 3) * x[j][i];
						model.addConstr(expr <= n - 2,
										"s_" + itos(i) + "_" + itos(j));
					}

				// arc flow
				for (i = 0; i < n - 1; i++)
					for (j = 1; j < n; j++)
					{
						GRBLinExpr expr = 0.0;
						expr += wt[i][j] * y[i][j] - theta[i][j];
						for (k = 1; k < n - 1; k++)
							expr += u[i][k][j] + u[k][j][i] - u[i][j][k];
						// when k==0
						if (i != 0)
							expr += u[0][j][i];
						// when k==n-1
						if (j != n - 1)
							expr += u[i][n - 1][j];
						model.addConstr(expr == 0,
										"flow_" + itos(i) + "_" + itos(j));
					}

				// arc flow upperbound
				for (i = 0; i < n - 1; i++)
					for (j = 1; j < n; j++)
					{
						GRBLinExpr expr = 0.0;
						expr += theta[i][j] - Q * x[i][j];
						model.addConstr(expr <= 0,
										"flowBound_" + itos(i) + "_" + itos(j));
					}

				// set objective
				GRBLinExpr obj = 0.0;
				for (i = 0; i < n - 1; i++)
					for (j = 1; j < n; j++)
					{
						obj += price * dis[i][j] * wt[i][j] * y[i][j];
						obj -= cost * dis[i][j] * theta[i][j];
						obj -= cost * vw * dis[i][j] * x[i][j];
					}
				model.setObjective(obj, GRB_MAXIMIZE);
				// model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
				printf("price = %lf\n", price);
				printf("cost = %lf\n", cost);
				printf("vw = %lf\n", vw);

				unordered_set<string> ratioSet;
				if (USE_RATIO_CUT)
				{
					// useRatioCuts
					// find triples that (d(ik) +d(kj))/d(ij) > price/cost

					double ratio;
					ratio = price / cost;
					// ratio = 2.8;
					cout << "ratio = " << ratio << endl;

					for (int i = 0; i < n; i++)
						for (int j = 0; j < n; j++)
							for (int k = 0; k < n; k++)
							{
								if (i != j && i != k && k != j)
								{
									if ((dis[i][k] + dis[k][j]) / dis[i][j] > ratio)
									{
										// u and y can't be positive at the same time
										// u(ijk)<=(1-y(ij))*Q

										GRBLinExpr expr = 0.0;
										expr += u[i][j][k] - (1 - y[i][j]) * Q;
										model.addConstr(
											expr <= 0,
											"ratio-cut_" + itos(i) + "_" + itos(j) + "_" + itos(k));

										// put the found triples indices into Hash set
										// cout << (dis[i][k] + dis[k][j]) / dis[i][j] << endl;
										string tripleIndex = itos(i) + "," + itos(j) + "," + itos(k);

										ratioSet.insert(tripleIndex);
									}
								}
							}
					cout << "ratioSet size = " << ratioSet.size() << endl;
					// exit (1);

					//  for (const auto &elem : ratioSet)
					//  {
					//  cout << elem << endl;
					//  }
				}

				if (USE_PROFITABLE_2_NODES_CYCLES)
				{
					for (int i = 0; i < n; i++)
						for (int j = i + 1; j < n; j++)
						{
							double subtourProfit = (price - cost) * (dis[i][j] * wt[i][j] + dis[j][i] * wt[j][i]) - cost * vw * (dis[i][j] + dis[j][i]);
							double subtourDist = dis[i][j] + dis[j][i];

							if (subtourProfit > 0 && subtourDist < route.DIS)
							{
								GRBLinExpr expr = 0.0;
								expr += x[i][j] + x[j][i];
								model.addConstr(
									expr <= 1,
									"2nodes_detour" + itos(i) + "_" + itos(j));

								printf("%d and %d subtour is added\n", i, j);
							}
						}
				}

				//****** Must set LazyConstraints parameter when using lazy constraints
				// 1 means use lazy constraint; 0 means not using lazy constraints
				// model.set(GRB_IntParam_LazyConstraints, 0);
				// model.set(GRB_IntParam_LazyConstraints, 1);

				//****** Set callback function
				// when USE_RATIO_CUT=false, ratioSet is a pointer to NULL
				// printIntSol cb = printIntSol(x, y, u, ratioSet, n);
				// model.setCallback(&cb);

				model.getEnv().set(GRB_IntParam_OutputFlag, 0); // set quiet mode

				// Optimize model
				model.optimize();

				// write model to file
				model.write("BPMP.lp");

				// Extract solution
				if (model.get(GRB_IntAttr_SolCount) > 0)
				{
					double **sol = new double *[n];
					printf("Selected arcs: \n");
					for (i = 0; i < n; i++)
					{
						sol[i] = model.get(GRB_DoubleAttr_X, x[i], n);
						for (int j = 0; j < n; j++)
							if (sol[i][j] > 0.9)
								printf("%d -- %d\n", i + 1, j + 1);
					}
					printf("Selected requests: \n");
					for (i = 0; i < n; i++)
					{
						sol[i] = model.get(GRB_DoubleAttr_X, y[i], n);
						for (int j = 0; j < n; j++)
							if (sol[i][j] > 0.9)
								printf("%d -- %d\n", i + 1, j + 1);
					}

					for (i = 0; i < n; i++)
						delete[] sol[i];
					delete[] sol;

					double objValue = model.get(GRB_DoubleAttr_ObjVal);
					cout << "obj value: " << objValue << endl;

					auto end = std::chrono::high_resolution_clock::now();
					chrono::duration<double> elapsed = end - start;
					cout << "Running time: " << elapsed.count() << " seconds" << endl;
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
				for (int j = 0; j < n; j++)
					delete[] u[i][j];
				delete[] u[i];
			}
			delete[] x;
			delete[] y;
			delete[] u;
			delete env;
		}
	}

	return 0;
}
