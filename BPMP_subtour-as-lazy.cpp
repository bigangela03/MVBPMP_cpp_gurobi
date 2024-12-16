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

bool USE_RATIO_CUT = false;
bool PRINT_VAR_VALUE = false;
bool PRINT_CYCLE_INFO = false;

bool ADD_TWO_NODES_DETOUR_ELIMINATION_IN_MASTER = true;

bool USE_PROFITABLE_2_NODES_CYCLES_IN_MASTER = true;

#define EPSILON 0.00001
double capacity;

// this version of BPMP reads all files in the target data folder
// so these global variables need to be initialized for each data file
int numCycles = 0;
vector<int> numNodesInCycles;

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
					printf("u vars: ");
					for (i = 0; i < n; i++)
						for (j = 0; j < n; j++)
							for (k = 0; k < n; k++)
								if (u[i][j][k] > EPSILON)
									printf("u[%d,%d,%d] = %lf", i, j, k, u[i][j][k]);
					printf("\n");
				}

				if (USE_RATIO_CUT)
				{
					for (i = 0; i < n; i++)
						for (j = 0; j < n; j++)
							for (k = 0; k < n; k++)
							{
								string key = itos(i) + "," + itos(j) + "," + itos(k);

								if (u[i][j][k] > EPSILON && y[i][j] > 0.5 && (ratioSet.find(key) != ratioSet.end()))
								{
									// cout << key << " found" << endl;
									// exit (1);

									// u and y can't be positive at the same time
									// u(ijk)<=(1-y(ij))*Q

									GRBLinExpr expr = 0.0;
									expr += uv[i][j][k] - (1 - yv[i][j]) * capacity;

									addLazy(expr <= 0);
									printf(
										"*** price cost ratio cut is added for (%d %d %d)***\n",
										i, j, k);
								}
							}
				}

				unordered_set<int> visitedNodesSet;
				bool findCycle = false;

				//========> Floyd-Warshall algorithm <========
				// Floyd-Warshall algorithm finds all pairs shortest path and can find negative cycle

				for (i = 0; i < n; i++)
					for (j = 0; j < n; j++)
					{
						if (x[i][j] > 0.5)
						{
							visitedNodesSet.insert(i);
							visitedNodesSet.insert(j);
						}
					}

				int dimVN = visitedNodesSet.size();

				int **VA = new int *[dimVN]; // VA: visited arcs
				int **minDist = new int *[dimVN];
				for (i = 0; i < dimVN; i++)
				{
					VA[i] = new int[dimVN];
					minDist[i] = new int[dimVN];
				}

				// the index in VA to the real node
				// for example, we found positive arcs 0 -> 9, 2 -> 5 -> 7 ->2
				// then visitedNodesSet = {0, 9, 2, 5, 7}
				// the corresponding index is [0, 1, 2, 3 ,4] in 2-dim array VA
				// so that VA[2][3] means arcs 2 -> 5 (arc 3 -> 6 in real world)
				int *I2N = new int[dimVN];

				i = 0;
				for (const auto &element : visitedNodesSet)
					I2N[i++] = element;

				for (i = 0; i < dimVN; i++)
					for (j = 0; j < dimVN; j++)
					{
						if (x[I2N[i]][I2N[j]] > 0.5)
						{
							VA[i][j] = -1;
							minDist[i][j] = -1;
						}
						else // let the arcs not selected with a big number
						// for small size under 100 nodes in graph, 10000 should be fine
						{
							VA[i][j] = 10000;
							minDist[i][j] = 10000;
						}
					}

				if (PRINT_CYCLE_INFO)
				{
					printf("============= find the detour =============\n");
					cout << "dimVN (dimention of visited nodes set)\n"
						 << dimVN << endl;

					cout << "===> visited notes:" << endl;
					for (const auto &element : visitedNodesSet)
						cout << element << " ";
					cout << endl;

					cout << "===> I2N set" << endl;
					for (i = 0; i < dimVN; i++)
						cout << I2N[i] << " ";
					cout << endl;

					cout << "===> minDist matrix:" << endl;
					for (i = 0; i < dimVN; i++)
					{
						for (j = 0; j < dimVN; j++)
						{
							cout << minDist[i][j] << " ";
						}
						cout << endl;
					}
				}

				// the main part of Floyd-Warshall algorithm
				int cycleIndex;
				for (k = 0; k < dimVN; k++)
				{
					for (i = 0; i < dimVN; i++)
						for (j = 0; j < dimVN; j++)
						{
							int min = minDist[i][k] + minDist[k][j];
							if (min > minDist[i][j])
								min = minDist[i][j];
							minDist[i][j] = min;
						}

					for (i = 0; i < dimVN; i++)
					{
						if (minDist[i][i] < 0)
						{
							// cout << "k=" << k << endl;
							// cout << "minDist = " << minDist[i][i] << endl;

							printf("===> a cycle is found from %d to %d!\n", I2N[i], I2N[i]);
							findCycle = true;
							cycleIndex = i;
						}
					}

					if (findCycle)
						break;
				}

				if (findCycle)
				{

					//========> find the arcs on the cycle <========
					vector<int> cyclePath;
					cyclePath.push_back(cycleIndex);
					int startingIndex = cycleIndex;
					int endingIndex = -1;
					int count = 0;
					for (j = 0; j < dimVN; j++)
					{
						if (VA[startingIndex][j] == -1)
						{
							cyclePath.push_back(j);
							count++;
							if (j == cycleIndex)
								break;
							//===> a cycle starts and ends with cycleIndex
							//===> so if the end node is not cycleIndex, start searching again
							startingIndex = j;
							j = -1; // so that after j++, j=0, search from 0 again
							if (count > dimVN)
							{
								cout << "ERROR: the program to find arcs in cycle doesn't work!" << endl;
								break;
							}
						}
					}

					//===> update global variables <===
					numCycles++;
					numNodesInCycles.push_back(cyclePath.size() - 1); // the cycle path starts and ends with the same index, i.e. cycleIndex
					//===> end of update <===

					if (PRINT_CYCLE_INFO)
					{
						cout << "cycle starting index is " << I2N[cycleIndex] << endl;

						cout << "the cycle path is " << endl;
						for (i = 0; i < cyclePath.size(); i++)
							cout << I2N[cyclePath[i]] << " ";
						cout << endl;
					}

					if ((ADD_TWO_NODES_DETOUR_ELIMINATION_IN_MASTER && (cyclePath.size() - 1 > 2)) || !ADD_TWO_NODES_DETOUR_ELIMINATION_IN_MASTER)
					{
						GRBLinExpr expr = 0.0;
						for (i = 0; i < cyclePath.size() - 1; i++)
						{
							expr += xv[I2N[cyclePath[i]]][I2N[cyclePath[i + 1]]];
							if (PRINT_CYCLE_INFO)
								printf("%d -> %d is added in subtour elimination constraint\n", I2N[cyclePath[i]], I2N[cyclePath[i + 1]]);
						}

						addLazy(expr <= cyclePath.size() - 2);
					}
				}

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

			//===> initialize global variables for cycle information
			numCycles = 0;
			numNodesInCycles.clear();
			//===> end of initialization

			string filename = path.string();
			route.readSingleFile(filename);

			// route.printStats();
			// route.printData();

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
				// for (i = 0; i < n - 1; i++)
				// 	for (j = 1; j < n; j++)
				// 	{
				// 		GRBLinExpr expr = 0.0;
				// 		expr += s[i] - s[j] + (n - 1) * x[i][j] + (n - 3) * x[j][i];
				// 		model.addConstr(expr <= n - 2,
				// 						"s_" + itos(i) + "_" + itos(j));
				// 	}

				if (ADD_TWO_NODES_DETOUR_ELIMINATION_IN_MASTER)
				{
					for (i = 0; i < n - 1; i++)
						for (j = i + 1; j < n; j++)
						{
							GRBLinExpr expr = 0.0;
							expr += x[i][j] + x[j][i];
							model.addConstr(expr <= 1,
											"detour-elimination_" + itos(i) + "_" + itos(j));
						}
				}

				if (USE_PROFITABLE_2_NODES_CYCLES_IN_MASTER)
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
					// ratio = price / cost;
					ratio = 2.8;
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
					/*
					 for (const auto &elem : ratioSet)
					 {
					 cout << elem << endl;
					 }
					 */
				}

				//****** Must set LazyConstraints parameter when using lazy constraints
				// 1 means use lazy constraint; 0 means not using lazy constraints
				// model.set(GRB_IntParam_LazyConstraints, 0);
				model.set(GRB_IntParam_LazyConstraints, 1);

				//****** Set callback function
				// when USE_RATIO_CUT=false, ratioSet is a pointer to NULL
				printIntSol cb = printIntSol(x, y, u, ratioSet, n);
				model.setCallback(&cb);

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
				}
				cout << "====================" << endl;
				cout << "number of cycles\n"
					 << numCycles << endl;

				cout << "number of nodes in each cycle:" << endl;

				for (i = 0; i < numNodesInCycles.size(); i++)
					cout << numNodesInCycles[i] << " ";
				cout << endl;
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
