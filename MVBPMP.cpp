#include "gurobi_c++.h"
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
using namespace std::string_literals;
//to generate vehicle data file name

int numV = 3; //number of vehicles
int numNodes; //number of nodes in instance
bool ADD_CUTS = false;
bool PRINT_x_WHEN_INTEGER_SOL = false;
bool PRINT_y_WHEN_INTEGER_SOL = true;
bool SEARCH_K_BEST_SOL = true;

string
itos (int i)
{
  stringstream s;
  s << i;
  return s.str ();
}

class printIntSol : public GRBCallback
{
public:

  int n;
  int numV;
  //GRBVar xv[n][n][n];
  //GRBVar yv[n][n][n];
  GRBVar ***xv;
  GRBVar ***yv;

  printIntSol (GRBVar ***xvars, GRBVar ***yvars, int nvar, int numVvar)
  {
    xv = xvars;
    yv = yvars;
    n = nvar;
    numV = numVvar;
  }
protected:
  void
  callback ()
  {
    try
      {
	if (where == GRB_CB_MIPSOL)
	  {
	    // Found an integer feasible solution
	    printf ("\n======> an integer solution found.\n");
	    int i, j, q;
	    double ***x = NULL;
	    double ***y = NULL;
	    int *overlapVeh = NULL;
	    int sum = 0;

	    x = new double**[n];
	    y = new double**[n];
	    overlapVeh = new int[numV];
	    for (i = 0; i < numV; i++)
	      overlapVeh[i] = -1;
	    for (i = 0; i < n; i++)
	      {
		x[i] = new double*[n];
		y[i] = new double*[n];
		for (j = 0; j < n; j++)
		  {
		    x[i][j] = new double[numV];
		    y[i][j] = new double[numV];
		  }
	      }
	    for (i = 0; i < n; i++)
	      {
		for (j = 0; j < n; j++)
		  {
		    x[i][j] = getSolution (xv[i][j], numV);
		    y[i][j] = getSolution (yv[i][j], numV);
		    sum = 0;
		    for (q = 0; q < numV; q++)
		      {
			if (PRINT_x_WHEN_INTEGER_SOL)
			  if (x[i][j][q] > 0.5)
			    printf ("x: %3d ->%3d (v%d)\n", i + 1, j + 1,
				    q + 1);
			if (PRINT_y_WHEN_INTEGER_SOL)
			  if (y[i][j][q] > 0.5)
			    printf ("y: %3d ->%3d (v%d)\n", i + 1, j + 1,
				    q + 1);
		      }
		    if (ADD_CUTS)
		      {
			for (q = 0; q < numV; q++)
			  if (y[i][j][q] > 0.5)
			    {
			      overlapVeh[q] = 1;
			      sum += 1;
			    }
			GRBLinExpr expr = 0;
			if (sum > 1.5)
			  {
			    printf ("======> request conflict found \n");
			    for (q = 0; q < numV; q++)
			      if (overlapVeh[q] == 1)
				expr += yv[i][j][q];
			    addLazy (expr <= 1);
			  }
		      }
		  }
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
	    delete[] overlapVeh;
	  }
      }
    catch (GRBException e)
      {
	cout << "Error number: " << e.getErrorCode () << endl;
	cout << e.getMessage () << endl;
      }
    catch (...)
      {
	cout << "Error during callback" << endl;
      }
  }
};

int
main (int argc, char *argv[])
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
  auto it = filesystem::directory_iterator (argv[1]);
  for (const auto &entry : it)
    {
      filesystem::path path = entry.path ();
      if (entry.is_regular_file ())
	{

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

	  //==============reading vehicle info==============

	  string vehicleFileName = "vehicle_location_data/vehicle_data_"
	      + to_string (n) + "_" + to_string (numV) + "_1.txt";
	  printf ("reading vehicle file: %s ...\n", vehicleFileName.c_str ());
	  //exit (1);
	  route.readSingleVehicleFile (vehicleFileName, numV);

	  //vehicle and origin always starts from 0
	  int vehicle[numV];
	  int origin[numV];

	  printf ("vehicles (start from 0):\n");
	  for (int i = 0; i < numV; i++)
	    {
	      vehicle[i] = route.vehicle[i];
	      printf ("%d ", vehicle[i]);
	    }
	  printf ("\norigins (start from 0):\n");
	  for (int i = 0; i < numV; i++)
	    {
	      origin[i] = route.origin[i];
	      printf ("%d ", origin[i]);
	    }
	  printf ("\n");

	  //==============start calling Gurobi===============
	  int i, j, k, q;
	  int status, nSolutions;

	  GRBEnv *env = NULL;
	  //GRBVar x[n][n][numV];
	  //GRBVar y[n][n][numV];
	  GRBVar s[n][numV];
	  GRBVar u[n][n][n][numV];
	  GRBVar theta[n][n][numV];

	  GRBVar ***x = NULL;
	  GRBVar ***y = NULL;
	  x = new GRBVar**[n];
	  y = new GRBVar**[n];
	  for (i = 0; i < n; i++)
	    {
	      x[i] = new GRBVar*[n];
	      y[i] = new GRBVar*[n];
	      for (j = 0; j < n; j++)
		{
		  x[i][j] = new GRBVar[numV];
		  y[i][j] = new GRBVar[numV];
		}
	    }

	  try
	    {
	      env = new GRBEnv ();
	      GRBModel model = GRBModel (*env);

	      // Must set LazyConstraints parameter when using lazy constraints
	      model.set (GRB_IntParam_LazyConstraints, 1);

	      // Create binary decision variables
	      for (q = 0; q < numV; q++)
		{
		  for (i = 0; i < n; i++)
		    {
		      s[i][q] = model.addVar (0.0, n, 0.0, GRB_CONTINUOUS,
					      "s_" + itos (i) + "_" + itos (q));
		      for (j = 0; j < n; j++)
			{
			  x[i][j][q] = model.addVar (
			      0.0,
			      1.0,
			      0,
			      GRB_BINARY,
			      "x_" + itos (i) + "_" + itos (j) + "_"
				  + itos (q));
			  y[i][j][q] = model.addVar (
			      0.0,
			      1.0,
			      0,
			      GRB_BINARY,
			      "y_" + itos (i) + "_" + itos (j) + "_"
				  + itos (q));
			  theta[i][j][q] = model.addVar (
			      0.0,
			      GRB_INFINITY,
			      0.0,
			      GRB_CONTINUOUS,
			      "theta_" + itos (i) + "_" + itos (j) + "_"
				  + itos (q));

			  for (k = 0; k < n; k++)
			    {
			      string s = "u_" + itos (i) + "_" + itos (j) + "_"
				  + itos (k) + "_" + itos (q);
			      u[i][j][k][q] = model.addVar (0.0, GRB_INFINITY,
							    0.0, GRB_CONTINUOUS,
							    s);
			    }
			}
		    }
		  int ogn = origin[q];
		  for (i = 0; i < n; i++)
		    {
		      x[i][i][q].set (GRB_DoubleAttr_UB, 0);
		      y[i][i][q].set (GRB_DoubleAttr_UB, 0);
		      x[i][ogn][q].set (GRB_DoubleAttr_UB, 0);
		      y[i][ogn][q].set (GRB_DoubleAttr_UB, 0);
		      x[n - 1][i][q].set (GRB_DoubleAttr_UB, 0);
		      y[n - 1][i][q].set (GRB_DoubleAttr_UB, 0);
		      for (j = 0; j < n; j++)
			if (wt[i][j] == 0)
			  y[i][j][q].set (GRB_DoubleAttr_UB, 0);
		    }
		}

	      // set up constraints
	      for (q = 0; q < numV; q++)
		{
		  int ogn = origin[q];

		  //vehicle goes out of origins
		  GRBLinExpr expr1 = 0.0;
		  for (i = 0; i < n; i++)
		    expr1 += x[ogn][i][q];

		  model.addConstr (expr1 == 1, "origin_" + itos (q));

		  //vehicle goes back to node n
		  GRBLinExpr expr2 = 0.0;
		  for (i = 0; i < n - 1; i++)
		    expr2 += x[i][n - 1][q];
		  model.addConstr (expr2 == 1, "destination_" + itos (q));

		  //flow conservation
		  for (int k = 0; k < n - 1; k++)
		    {
		      if (k != ogn)
			{
			  GRBLinExpr expr = 0;
			  for (i = 0; i < n - 1; i++)
			    expr += x[i][k][q];
			  for (j = 1; j < n; j++)
			    expr -= x[k][j][q];
			  model.addConstr (
			      expr == 0,
			      "flow_conservation_" + itos (k) + "_" + itos (q));
			}
		    }

		  //distance
		  GRBLinExpr expr3 = 0.0;
		  for (i = 0; i < n - 1; i++)
		    for (j = 0; j < n; j++)
		      expr3 += dis[i][j] * x[i][j][q];
		  model.addConstr (expr3 <= route.DIS, "distance_" + itos (q));

		  //node degree less than 1
		  for (int j = 0; j < n - 1; j++)
		    {
		      GRBLinExpr expr = 0.0;
		      for (i = 0; i < n - 1; i++)
			expr += x[i][j][q];
		      model.addConstr (expr <= 1,
				       "indegree_" + itos (j) + "_" + itos (q));
		    }

		  //subtour elimination
		  for (i = 0; i < n - 1; i++)
		    for (j = 0; j < n; j++)
		      {
			GRBLinExpr expr = 0.0;
			expr += s[i][q] - s[j][q] + (n - 1) * x[i][j][q]
			    + (n - 3) * x[j][i][q];
			model.addConstr (
			    expr <= n - 2,
			    "subtour_" + itos (i) + "_" + itos (j) + "_"
				+ itos (q));
		      }

		  //arc flow
		  for (i = 0; i < n - 1; i++)
		    for (j = 0; j < n; j++)
		      {
			GRBLinExpr expr = 0.0;
			expr += wt[i][j] * y[i][j][q] - theta[i][j][q];
			for (k = 0; k < n - 1; k++)
			  {
			    if (k != ogn)
			      expr += u[i][k][j][q] + u[k][j][i][q]
				  - u[i][j][k][q];
			    else
			      expr += u[k][j][i][q];
			  }
			//when k==n-1
			if (j != n - 1)
			  expr += u[i][n - 1][j][q];
			model.addConstr (
			    expr == 0,
			    "flow_" + itos (i) + "_" + itos (j) + "_"
				+ itos (q));
		      }

		  //arc flow upperbound
		  for (i = 0; i < n - 1; i++)
		    for (j = 0; j < n; j++)
		      {
			GRBLinExpr expr = 0.0;
			expr += theta[i][j][q] - Q * x[i][j][q];
			model.addConstr (
			    expr <= 0,
			    "flowBound_" + itos (i) + "_" + itos (j) + "_"
				+ itos (q));
		      }

		  //add redundant constraint to test if same multi obj value
		  //caused by binary variable
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

	      //one vehicle for one cargo
	      if (!ADD_CUTS)
		{
		  for (i = 0; i < n - 1; i++)
		    for (j = 0; j < n; j++)
		      {
			GRBLinExpr expr = 0.0;
			for (q = 0; q < numV; q++)
			  expr += y[i][j][q];
			model.addConstr (expr <= 1,
					 "one-one" + itos (i) + "_" + itos (j));
		      }
		}

	      //set up objective
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

	      model.setObjective (obj, GRB_MAXIMIZE);
	      //model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
	      //printf ("price = %lf\n", price);
	      //printf ("cost = %lf\n", cost);
	      //printf ("vw = %lf\n", vw);

	      // Set callback function
	      //printIntSol cb = printIntSol (x, y, n, numV);
	      //model.setCallback (&cb);

	      if (SEARCH_K_BEST_SOL)
		{
		  // do a systematic search for the k-best solutions
		  model.set (GRB_IntParam_PoolSearchMode, 2);
		  // Limit how many solutions to collect
		  model.set (GRB_IntParam_PoolSolutions, 150);
		}

	      // Optimize model
	      model.optimize ();

	      //write model to file
	      model.write ("MVBPMP.lp");

	      // Status checking
	      status = model.get (GRB_IntAttr_Status);
	      if (status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE
		  || status == GRB_UNBOUNDED)
		{
		  cout << "The model cannot be solved "
		      << "because it is infeasible or unbounded" << endl;
		  return 1;
		}
	      if (status != GRB_OPTIMAL)
		{
		  cout << "Optimization was stopped with status " << status
		      << endl;
		  return 1;
		}

	      //retrieve K best objs
	      if (SEARCH_K_BEST_SOL)
		{
		  // Print number of solutions stored
		  nSolutions = model.get (GRB_IntAttr_SolCount);
		  cout << "Number of solutions found: " << nSolutions << endl;

		  // Print objective values of solutions
		  for (int e = 0; e < nSolutions; e++)
		    {
		      model.set (GRB_IntParam_SolutionNumber, e);
		      cout << model.get (GRB_DoubleAttr_PoolObjVal) << " ";
		      if (e % 10 == 9)
			cout << endl;
		    }
		  cout << endl;

		  // print
		  double oldobj = -1;
		  double newobj = 0;
		  //for (int e = 0; e < nSolutions; e++)
		  for (int e = 0; e < 10; e++)
		    {
		      model.set (GRB_IntParam_SolutionNumber, e);
		      newobj = model.get (GRB_DoubleAttr_PoolObjVal);
		      if (newobj == oldobj)
			continue;

		      oldobj = newobj;

		      printf ("=== %d (obj=%lf) ===\n", e + 1, newobj);
		      printf ("------ print x=1 ------\n");
		      for (q = 0; q < numV; q++)
			for (i = 0; i < n - 1; i++)
			  for (j = 0; j < n; j++)
			    {
			      double xval = x[i][j][q].get (GRB_DoubleAttr_Xn);
			      if (xval > 0.99 && xval < 1.01)
				printf ("%3d %3d  v%d\n", i + 1, j + 1, q + 1);
			    }
		      printf ("------ print y=1 ------\n");
		      for (q = 0; q < numV; q++)
			for (i = 0; i < n - 1; i++)
			  for (j = 0; j < n; j++)
			    {
			      double yval = y[i][j][q].get (GRB_DoubleAttr_Xn);
			      if (yval > 0.99 && yval < 1.01)
				printf ("%3d %3d  v%d (w=%lf)\n", i + 1, j + 1,
					q + 1, wt[i][j]);
			    }
		      printf ("------ print theta>0.000001 ------\n");
		      for (q = 0; q < numV; q++)
			for (i = 0; i < n - 1; i++)
			  for (j = 0; j < n; j++)
			    {
			      double thetaval = theta[i][j][q].get (
				  GRB_DoubleAttr_Xn);
			      if (thetaval > 0.000001)
				printf ("%3d %3d  v%d flow=%lf\n", i + 1, j + 1,
					q + 1, thetaval);
			    }
		      printf ("------ print u>0.000001 ------\n");
		      for (q = 0; q < numV; q++)
			{
			  int ogn = origin[q];
			  for (i = 0; i < n - 1; i++)
			    for (j = 0; j < n; j++)
			      for (k = 0; k < n - 1; k++)
				if (k != ogn)
				  {
				    double uval = u[i][j][k][q].get (
					GRB_DoubleAttr_Xn);
				    if (uval > 0.000001)
				      printf ("%3d %3d %3d  v%d\n", i + 1,
					      j + 1, k + 1, q + 1);
				  }
			}

		    }
		}

	      // Extract solution
	      if (model.get (GRB_IntAttr_SolCount) > 0)
		{
		  double ***solx = NULL;
		  double ***soly = NULL;
		  solx = new double**[n];
		  soly = new double**[n];
		  for (i = 0; i < n; i++)
		    {
		      solx[i] = new double*[n];
		      soly[i] = new double*[n];
		      for (j = 0; j < n; j++)
			{
			  solx[i][j] = model.get (GRB_DoubleAttr_X, x[i][j],
						  numV);
			  soly[i][j] = model.get (GRB_DoubleAttr_X, y[i][j],
						  numV);
			}
		    }

		  printf ("Selected arcs: \n");
		  for (q = 0; q < numV; q++)
		    for (i = 0; i < n; i++)
		      for (j = 0; j < n; j++)
			if (solx[i][j][q] > 0.9)
			  printf ("%d -- %d (%d)\n", i + 1, j + 1, q + 1);
		  printf ("Selected requests: \n");
		  for (q = 0; q < numV; q++)
		    for (i = 0; i < n; i++)
		      for (j = 0; j < n; j++)
			if (soly[i][j][q] > 0.9)
			  printf ("%d -- %d (%d)\n", i + 1, j + 1, q + 1);

		  for (i = 0; i < n; i++)
		    {
		      for (j = 0; j < n; j++)
			{
			  delete[] solx[i][j];
			  delete[] soly[i][j];
			}
		      delete[] solx[i];
		      delete[] soly[i];
		    }
		  delete[] solx;
		  delete[] soly;
		}
	    }
	  catch (GRBException e)
	    {
	      cout << "Error number: " << e.getErrorCode () << endl;
	      cout << e.getMessage () << endl;
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

	}
    }

  return 0;
}
