#include <gurobi_c++.h>
#include <omp.h>
#include <sstream>
#include <iostream>

#include <cassert>
#include <cstdlib>
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

int numVehicles = 3; //number of vehicles; will be overwritten when reading vehicle file
int numNodes; //number of nodes in instance

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
  if (argc != 3)
    {
      cout << "Usage: ./openmpMVBPMP.x dataFolder numberOfVehicles" << endl;
      return 1;
    }

  cout << argc << endl;
  for (int i = 1; i < argc; ++i)
    cout << argv[i] << endl;

  numVehicles = stoi (argv[2]);

  readData route;
  auto it = filesystem::directory_iterator (argv[1]);
  for (const auto &entry : it)
    {
      filesystem::path path = entry.path ();
      if (entry.is_regular_file ())
	{

	  //==============reading nodes, price, cost and other info==============

	  string filename = path.string ();
	  printf ("reading file %s ...\n", filename.c_str ());
	  route.readSingleFile (filename);

	  //route.printStats();
	  //route.printData();

	  int n = route.numOfNode;

	  numNodes = route.numOfNode;
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

	  //==============reading vehicles' info (2 or 3 vehicles)==============

	  string vehicleFileName = "vehicle_location_data/vehicle_data_"
	      + to_string (n) + "_" + to_string (numVehicles) + "_1.txt";
	  printf ("reading vehicle file: %s ...\n", vehicleFileName.c_str ());
	  //exit (1);
	  route.readSingleVehicleFile (vehicleFileName, numVehicles);

	  //vehicle and origin array always starts from 0
	  int vehicle[numVehicles];
	  int origin[numVehicles];

	  printf ("vehicles (start from 0):\n");
	  for (int i = 0; i < numVehicles; i++)
	    {
	      vehicle[i] = route.vehicle[i];
	      printf ("%d ", vehicle[i]);
	    }
	  printf ("\norigins (start from 0):\n");
	  for (int i = 0; i < numVehicles; i++)
	    {
	      origin[i] = route.origin[i];
	      printf ("%d ", origin[i]);
	    }
	  printf ("\n");

	  //==============start calling Gurobi using openMP===============
	  omp_set_nested (1);
	  int requiredThreads = 8;

	  //require a master thread for each model
#pragma omp parallel num_threads(numVehicles)
	    {
	      for (int vehicleIndex = 0; vehicleIndex < numVehicles;
		  vehicleIndex++)
		{
		  if (omp_get_thread_num () == vehicleIndex)
		    {
		      //require a given number of threads for each model
#pragma omp parallel num_threads(requiredThreads)
			{
			  int nthreads = omp_get_num_threads ();
			  cout << "NUM_THREADS=" << nthreads << endl;

			  //start calling Gurobi
			  int i, j, k;
			  int vehicleOrigin = origin[vehicleIndex];

			  GRBEnv *env = NULL;
			  GRBVar **x = NULL;
			  GRBVar **y = NULL;
			  GRBVar s[n];
			  GRBVar u[n][n][n];
			  GRBVar theta[n][n];

			  x = new GRBVar*[n];
			  y = new GRBVar*[n];
			  for (i = 0; i < n; i++)
			    {
			      x[i] = new GRBVar[n];
			      y[i] = new GRBVar[n];
			    }

			  try
			    {
			      env = new GRBEnv ();
			      GRBModel model = GRBModel (*env);

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
				      y[i][j] = model.addVar (
					  0.0, 1.0, 0, GRB_BINARY,
					  "y_" + itos (i) + "_" + itos (j));
				      theta[i][j] = model.addVar (
					  0.0, GRB_INFINITY, 0.0,
					  GRB_CONTINUOUS,
					  "theta_" + itos (i) + "_" + itos (j));

				      for (k = 0; k < n; k++)
					{
					  string s = "u_" + itos (i) + "_"
					      + itos (j) + "_" + itos (k);
					  u[i][j][k] = model.addVar (
					      0.0, GRB_INFINITY, 0.0,
					      GRB_CONTINUOUS, s);
					}
				    }
				}

			      //force some x, y variables to be zeros
			      for (i = 0; i < n; i++)
				{
				  x[i][i].set (GRB_DoubleAttr_UB, 0);
				  y[i][i].set (GRB_DoubleAttr_UB, 0);
				  x[i][vehicleOrigin].set (GRB_DoubleAttr_UB,
							   0);
				  y[i][vehicleOrigin].set (GRB_DoubleAttr_UB,
							   0);
				  x[n - 1][i].set (GRB_DoubleAttr_UB, 0);
				  y[n - 1][i].set (GRB_DoubleAttr_UB, 0);
				}
			      //force some triples variables to be zeros
			      for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
				  {
				    u[n - 1][i][j].set (GRB_DoubleAttr_UB, 0);
				    u[i][vehicleOrigin][j].set (
					GRB_DoubleAttr_UB, 0);
				    u[i][j][vehicleOrigin].set (
					GRB_DoubleAttr_UB, 0);
				    u[i][j][n - 1].set (GRB_DoubleAttr_UB, 0);
				  }

			      //==============generate constraints in Gurobi================
			      //vehicle goes out of vehicle's origin
			      GRBLinExpr expr1 = 0.0;
			      for (i = 0; i < n; i++)
				expr1 += x[vehicleOrigin][i];
			      model.addConstr (expr1 == 1, "origin");

			      //vehicle goes back to node n
			      GRBLinExpr expr2 = 0.0;
			      for (i = 0; i < n; i++)
				expr2 += x[i][n - 1];
			      model.addConstr (expr2 == 1, "destination");

			      //flow conservation
			      for (k = 0; k < n; k++)
				{
				  if (k != vehicleOrigin)
				    {
				      GRBLinExpr expr = 0;
				      for (i = 0; i < n; i++)
					expr += x[i][k];
				      for (j = 0; j < n; j++)
					expr -= x[k][j];
				      model.addConstr (
					  expr == 0,
					  "flow_conservation_" + itos (k));
				    }
				}

			      //distance
			      GRBLinExpr expr3 = 0.0;
			      for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
				  expr3 += dis[i][j] * x[i][j];
			      model.addConstr (expr3 <= route.DIS, "distance");

			      //node degree less than 1
			      for (k = 0; k < n; k++)
				{
				  GRBLinExpr expr = 0.0;
				  for (i = 0; i < n - 1; i++)
				    expr += x[i][k];
				  model.addConstr (expr <= 1,
						   "indegree_" + itos (k));
				}

			      //subtour elimination
			      for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
				  {
				    GRBLinExpr expr = 0.0;
				    expr += s[i] - s[j] + (n - 1) * x[i][j]
					+ (n - 3) * x[j][i];
				    model.addConstr (
					expr <= n - 2,
					"s_" + itos (i) + "_" + itos (j));
				  }

			      //arc flow
			      for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
				  {
				    GRBLinExpr expr = 0.0;
				    expr += wt[i][j] * y[i][j] - theta[i][j];
				    for (k = 0; k < n; k++)
				      expr += u[i][k][j] + u[k][j][i]
					  - u[i][j][k];
				    model.addConstr (
					expr == 0,
					"flow_" + itos (i) + "_" + itos (j));
				  }

			      //arc flow upperbound
			      for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
				  {
				    GRBLinExpr expr = 0.0;
				    expr += theta[i][j] - Q * x[i][j];
				    model.addConstr (
					expr <= 0,
					"flowBound_" + itos (i) + "_"
					    + itos (j));
				  }

			      //==============set objective function in Gurobi================

			      GRBLinExpr obj = 0.0;
			      for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
				  {
				    obj += price * dis[i][j] * wt[i][j]
					* y[i][j];
				    obj -= cost * dis[i][j] * theta[i][j];
				    obj -= cost * vw * dis[i][j] * x[i][j];
				  }
			      model.setObjective (obj, GRB_MAXIMIZE);

			      //==============solve the model in Gurobi================
			      model.optimize ();

			      //write model to file
			      model.write (
				  "BPMP_" + itos (vehicleIndex) + ".lp");

			      //==============Extract solution from Gurobi================
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
				  printf ("Selected requests: \n");
				  for (i = 0; i < n; i++)
				    {
				      sol[i] = model.get (GRB_DoubleAttr_X,
							  y[i], n);
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
			      delete[] y[i];
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
