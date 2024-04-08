#include <gurobi_c++.h>
#include <omp.h>
#include <sstream>
#include <iostream>

#include <cassert>
#include <cstdlib>
//this is to test on TSP with the same data
#include <cmath>
#include <sstream>

#include <iostream> //input and output
#include <fstream>  //read file
#include <sstream>
#include <algorithm>
#include <vector>
#include <climits>
#include <set>
#include <string>
#include <ctime>  //to measure CPU time
#include <chrono> //to measure run time
#include "readData.h"

using namespace std;
using namespace std::chrono;

string
itos (int i)
{
  stringstream s;
  s << i;
  return s.str ();
}

int
main (int argc, char *argv[])
{
  if (argc < 2)
    {
      cout << "Usage: ./test.x dataFolder" << endl;
      return 1;
    }

  cout << argc << endl;
  for (int i = 1; i < argc; ++i)
    cout << argv[i] << endl;

  readData route;
  auto it = filesystem::directory_iterator (argv[1]);
  for (const auto &entry : it)
    {
      filesystem::path path = entry.path ();
      if (entry.is_regular_file ())
	{
	  string filename = path.string ();
	  route.readSingleFile (filename);

	  //route.printStats();
	  //route.printData();

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

	  //start calling Gurobi

	  //this is to test on TSP with the same data

	  int NUM_THREADS = omp_get_num_threads ();
	  cout << "NUM_THREADS=" << NUM_THREADS << endl;

	  int numV = 3;
	  int requiredThreads = 8;

	 // omp_set_nested (1);

#pragma omp parallel num_threads(numV)
	    {
	      for (int ip = 0; ip < numV; ip++)
		{
		  if (omp_get_thread_num () == ip)
		    {
#pragma omp parallel num_threads(requiredThreads)
			{
			  int nthreads = omp_get_num_threads ();
			  cout << "NUM_THREADS=" << nthreads << endl;

			  int i, j, k;

			  GRBEnv *env = NULL;
			  GRBVar **x = NULL;
			  GRBVar s[n];

			  x = new GRBVar*[n];
			  for (i = 0; i < n; i++)
			    x[i] = new GRBVar[n];

			  try
			    {
			      env = new GRBEnv ();
			      GRBModel model = GRBModel (*env);

			      model.set (GRB_IntParam_OutputFlag, 0);

			      // Must set LazyConstraints parameter when using lazy constraints
			      //model.set (GRB_IntParam_LazyConstraints, 1);

			      // Create binary decision variables
			      for (i = 0; i < n; i++)
				{
				  s[i] = model.addVar (0.0, n, 0.0,
						       GRB_CONTINUOUS,
						       "s_" + itos (i));
				  for (j = 0; j < n; j++)
				    {
				      x[i][j] = model.addVar (
					  0.0, 1.0, 0, GRB_BINARY,
					  "x_" + itos (i) + "_" + itos (j));
				    }
				}

			      for (i = 0; i < n; i++)
				{
				  x[i][i].set (GRB_DoubleAttr_UB, 0);
				  x[i][0].set (GRB_DoubleAttr_UB, 0);
				  x[n - 1][i].set (GRB_DoubleAttr_UB, 0);

				}

			      //vehicle goes out of node 1
			      GRBLinExpr expr1 = 0.0;
			      for (i = 1; i < n; i++)
				expr1 += x[0][i];
			      model.addConstr (expr1 == 1, "origin");

			      //vehicle goes back to node n
			      GRBLinExpr expr2 = 0.0;
			      for (i = 0; i < n - 1; i++)
				expr2 += x[i][n - 1];
			      model.addConstr (expr2 == 1, "destination");

			      //flow conservation
			      for (int k = 1; k < n - 1; k++)
				{
				  GRBLinExpr expr = 0;
				  for (i = 0; i < n - 1; i++)
				    expr += x[i][k];
				  for (j = 1; j < n; j++)
				    expr -= x[k][j];
				  model.addConstr (
				      expr == 0,
				      "flow_conservation_" + itos (k));
				}

			      //node degree
			      for (int k = 1; k < n - 1; k++)
				{
				  GRBLinExpr expr = 0.0;
				  for (i = 0; i < n - 1; i++)
				    expr += x[i][k];
				  model.addConstr (expr == 1,
						   "indegree_" + itos (k));
				}

			      //subtour elimination
			      for (i = 0; i < n - 1; i++)
				for (j = 1; j < n; j++)
				  {
				    GRBLinExpr expr = 0.0;
				    expr += s[i] - s[j] + (n - 1) * x[i][j]
					+ (n - 3) * x[j][i];
				    model.addConstr (
					expr <= n - 2,
					"s_" + itos (i) + "_" + itos (j));
				  }

			      //set objective
			      GRBLinExpr obj = 0.0;
			      for (i = 0; i < n - 1; i++)
				for (j = 1; j < n; j++)
				  {
				    obj += dis[i][j] * x[i][j];
				  }
			      model.setObjective (obj, GRB_MINIMIZE);
			      //model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

			      // Optimize model
			      model.optimize ();

			      //write model to file
			      model.write ("TSP" + itos (ip) + ".lp");

			      // Extract solution
			      if (model.get (GRB_IntAttr_SolCount) > 0)
				{
				  double **sol = new double*[n];
				  printf ("Selected arcs: \n");
				  for (i = 0; i < n; i++)
				    {
				      sol[i] = model.get (GRB_DoubleAttr_X,
							  x[i], n);
				      for (int j = 0; j < n; j++)
					if (sol[i][j] > 0.9)
					  printf ("%d -- %d\n", i + 1, j + 1);
				    }

				  for (i = 0; i < n; i++)
				    delete[] sol[i];
				}
			    }
			  catch (GRBException e)
			    {
			      cout << "Error number: " << e.getErrorCode ()
				  << endl;
			      cout << e.getMessage () << endl;
			    }
			  catch (...)
			    {
			      cout << "Error during optimization" << endl;
			    }

			  for (i = 0; i < n; i++)
			    {
			      delete[] x[i];
			    }
			  delete env;

			}
		    }

		}
	    }
	}
    }

  return 0;
}
