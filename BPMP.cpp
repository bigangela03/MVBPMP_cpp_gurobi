/* Copyright 2023, Gurobi Optimization, LLC */

/* Solve a traveling salesman problem on a randomly generated set of
 points using lazy constraints.   The base MIP model only includes
 'degree-2' constraints, requiring each node to have exactly
 two incident edges.  Solutions to this model may contain subtours -
 tours that don't visit every node.  The lazy constraint callback
 adds new constraints to cut them off. */

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
  GRBVar **xv;
  GRBVar **yv;
  int n;
  printIntSol (GRBVar **xvars, GRBVar **yvars, int nvar)
  {
    xv = xvars;
    yv = yvars;
    n = nvar;
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
	    int i, j;
	    double **x = new double*[n];
	    double **y = new double*[n];
	    printf ("x vars: ");
	    for (i = 0; i < n; i++)
	      {
		x[i] = getSolution (xv[i], n);
		for (j = 0; j < n; j++)
		  {
		    if (x[i][j] > 0.5)
		      printf ("%d->%d ", i, j);
		  }
	      }
	    printf ("\ny vars: ");
	    for (i = 0; i < n; i++)
	      {
		y[i] = getSolution (yv[i], n);
		for (j = 0; j < n; j++)
		  {
		    if (y[i][j] > 0.5)
		      printf ("%d->%d ", i, j);
		  }
	      }
	    printf ("\n");
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
	  int i;

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
	      int j, k;

	      env = new GRBEnv ();
	      GRBModel model = GRBModel (*env);

	      // Must set LazyConstraints parameter when using lazy constraints
	      model.set (GRB_IntParam_LazyConstraints, 1);

	      // Create binary decision variables
	      for (i = 0; i < n; i++)
		{
		  s[i] = model.addVar (0.0, n, 0.0, GRB_CONTINUOUS,
				       "s_" + itos (i));
		  for (j = 0; j < n; j++)
		    {
		      x[i][j] = model.addVar (0.0, 1.0, 0, GRB_BINARY,
					      "x_" + itos (i) + "_" + itos (j));
		      y[i][j] = model.addVar (0.0, 1.0, 0, GRB_BINARY,
					      "y_" + itos (i) + "_" + itos (j));
		      theta[i][j] = model.addVar (
			  0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
			  "theta_" + itos (i) + "_" + itos (j));

		      for (k = 0; k < n; k++)
			{
			  string s = "u_" + itos (i) + "_" + itos (j) + "_"
			      + itos (k);
			  u[i][j][k] = model.addVar (0.0, GRB_INFINITY, 0.0,
						     GRB_CONTINUOUS, s);
			}
		    }

		}

	      for (i = 0; i < n; i++)
		{
		  x[i][i].set (GRB_DoubleAttr_UB, 0);
		  y[i][i].set (GRB_DoubleAttr_UB, 0);
		  x[i][0].set (GRB_DoubleAttr_UB, 0);
		  y[i][0].set (GRB_DoubleAttr_UB, 0);
		  x[n - 1][i].set (GRB_DoubleAttr_UB, 0);
		  y[n - 1][i].set (GRB_DoubleAttr_UB, 0);
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
		  model.addConstr (expr == 0, "flow_conservation_" + itos (k));
		}

	      //distance
	      GRBLinExpr expr3 = 0.0;
	      for (i = 0; i < n - 1; i++)
		for (j = 1; j < n; j++)
		  expr3 += dis[i][j] * x[i][j];
	      model.addConstr (expr3 <= route.DIS, "distance");

	      //node degree less than 1
	      for (int k = 1; k < n - 1; k++)
		{
		  GRBLinExpr expr = 0.0;
		  for (i = 0; i < n - 1; i++)
		    expr += x[i][k];
		  model.addConstr (expr <= 1, "indegree_" + itos (k));
		}

	      //subtour elimination
	      for (i = 0; i < n - 1; i++)
		for (j = 1; j < n; j++)
		  {
		    GRBLinExpr expr = 0.0;
		    expr += s[i] - s[j] + (n - 1) * x[i][j] + (n - 3) * x[j][i];
		    model.addConstr (expr <= n - 2,
				     "s_" + itos (i) + "_" + itos (j));
		  }

	      //arc flow
	      for (i = 0; i < n - 1; i++)
		for (j = 1; j < n; j++)
		  {
		    GRBLinExpr expr = 0.0;
		    expr += wt[i][j] * y[i][j] - theta[i][j];
		    for (k = 1; k < n - 1; k++)
		      expr += u[i][k][j] + u[k][j][i] - u[i][j][k];
		    //when k==0
		    if (i != 0)
		      expr += u[0][j][i];
		    //when k==n-1
		    if (j != n - 1)
		      expr += u[i][n - 1][j];
		    model.addConstr (expr == 0,
				     "flow_" + itos (i) + "_" + itos (j));
		  }

	      //arc flow upperbound
	      for (i = 0; i < n - 1; i++)
		for (j = 1; j < n; j++)
		  {
		    GRBLinExpr expr = 0.0;
		    expr += theta[i][j] - Q * x[i][j];
		    model.addConstr (expr <= 0,
				     "flowBound_" + itos (i) + "_" + itos (j));
		  }

	      //set objective
	      GRBLinExpr obj = 0.0;
	      for (i = 0; i < n - 1; i++)
		for (j = 1; j < n; j++)
		  {
		    obj += price * dis[i][j] * wt[i][j] * y[i][j];
		    obj -= cost * dis[i][j] * theta[i][j];
		    obj -= cost * vw * dis[i][j] * x[i][j];
		  }
	      model.setObjective (obj, GRB_MAXIMIZE);
	      //model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
	      printf ("price = %lf\n", price);
	      printf ("cost = %lf\n", cost);
	      printf ("vw = %lf\n", vw);

	      // Set callback function
	      printIntSol cb = printIntSol (x, y, n);
	      model.setCallback (&cb);

	      // Optimize model
	      model.optimize ();

	      //write model to file
	      model.write ("BPMP.lp");

	      // Extract solution
	      if (model.get (GRB_IntAttr_SolCount) > 0)
		{
		  double **sol = new double*[n];
		  printf ("Selected arcs: \n");
		  for (i = 0; i < n; i++)
		    {
		      sol[i] = model.get (GRB_DoubleAttr_X, x[i], n);
		      for (int j = 0; j < n; j++)
			if (sol[i][j] > 0.9)
			  printf ("%d -- %d\n", i + 1, j + 1);
		    }
		  printf ("Selected requests: \n");
		  for (i = 0; i < n; i++)
		    {
		      sol[i] = model.get (GRB_DoubleAttr_X, y[i], n);
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
	      cout << "Error number: " << e.getErrorCode () << endl;
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

  return 0;
}
