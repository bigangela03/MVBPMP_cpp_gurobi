//Angela; Feb3, 2024; this code do not use callback function

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

double price;
double cost_flow;
double cost_vehicle;
double vw;
double Q;
double **wt;
double **dis;

double UB;
double LB;

#define epsilon 0.001

bool initialSolution = true;
bool printMore = false;
bool PRINT_UNBOUNDED_VECTOR = false;

bool DUAL_USE_GRB = true;

GRBEnv *env; //set gurobi environment as a global var
//and define the environment only once in main(),
//then use env.start() to start an empty env.
//when using env = new GRBEnv (), system shows the following message
//Set parameter TokenServer to value "sengr7lic2.smu.edu"

string
itos (int i)
{
  stringstream s;
  s << i;
  return s.str ();
}

bool
subProblemDual (double**, double**, int, double**, double**, double*);
void
deleteVar (double**, int);

/*
 class printIntSol : public GRBCallback
 {
 public:
 GRBVar **xv;
 GRBVar **yv;
 GRBVar *zv;
 int n;
 printIntSol (GRBVar **xvars, GRBVar **yvars, GRBVar *zvar, int nvar)
 {
 xv = xvars;
 yv = yvars;
 zv = zvar;
 n = nvar;
 }
 protected:
 void
 callback ()
 {
 try
 {
 int i, j;

 if (where == GRB_CB_MIPSOL)
 {
 // Found an integer feasible solution
 printf ("==========GRB_CB_MIPSOL=========\n");

 int nodecnt = (int) getDoubleInfo (GRB_CB_MIPSOL_NODCNT);
 double obj = getDoubleInfo (GRB_CB_MIPSOL_OBJ);
 double objbst = getDoubleInfo (GRB_CB_MIPSOL_OBJBST);
 double objbnd = getDoubleInfo (GRB_CB_MIPSOL_OBJBND);
 int solcnt = getIntInfo (GRB_CB_MIPSOL_SOLCNT);

 printf ("nodecnt=%d\n", nodecnt);
 printf ("obj=%lf\n", obj);
 printf ("objbst=%lf\n", objbst);
 printf ("objbnd=%lf\n", objbnd);
 printf ("solcnt=%d\n", solcnt);

 //call sub-problem
 double diffRate = (objbnd - objbst) / objbnd;
 if (diffRate < 0.1)
 {
 }

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
 }
 ;
 */

void
deleteVar (double **x, int n)
{
  for (int i = 0; i < n; i++)
    {
      delete[] x[i];
    }
  delete[] x;
}

bool
subProblemDual (double **xVal, double **yVal, int n, double **piVal,
		double **unbdPiVal, double *subProProfit)
{
  //printf ("********* start sub-problem *********\n");

  int i, j, k;

  GRBModel *SP = new GRBModel (*env);		    //sub-problem model

  (*SP).set (GRB_IntParam_OutputFlag, 0);
  //To obtain a more definitive conclusion, set the DualReductions parameter to 0 and reoptimize
  //from https://www.gurobi.com/documentation/9.1/refman/optimization_status_codes.html#sec:StatusCodes
  (*SP).set (GRB_IntParam_DualReductions, 0);
  //To obtain an unbound ray
  (*SP).set (GRB_IntParam_InfUnbdInfo, 1);

  GRBVar **piVar = new GRBVar*[n];
  for (i = 0; i < n; i++)
    piVar[i] = new GRBVar[n];

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	piVar[i][j] = SP->addVar (0.0, GRB_INFINITY, 0, GRB_CONTINUOUS,
				  "pi_" + itos (i) + "_" + itos (j));
      }

  for (i = 0; i < n; i++)
    {
      piVar[i][i].set (GRB_DoubleAttr_UB, 0);
      piVar[i][0].set (GRB_DoubleAttr_UB, 0);
      piVar[n - 1][i].set (GRB_DoubleAttr_UB, 0);
    }

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
	{
	  GRBLinExpr expr = 0.0;
	  expr += piVar[i][k] + piVar[k][j] - piVar[i][j];
	  SP->addConstr (
	      expr >= -cost_flow * (dis[i][k] + dis[k][j] - dis[i][j]),
	      "dual" + itos (i) + "_" + itos (j) + "_" + itos (k));
	}
//set objective
  GRBLinExpr obj = 0.0;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      obj += (Q * xVal[i][j] - wt[i][j] * yVal[i][j]) * piVar[i][j];
  SP->setObjective (obj, GRB_MINIMIZE);
//model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

// Optimize model
  SP->optimize ();

//write model to file
  SP->write ("BPMP_benders_SP.lp");

  //printf ("GRB_IntAttr_SolCount=%d\n", SP->get (GRB_IntAttr_SolCount));
  //printf ("GRB_IntAttr_Status=%d\n", SP->get (GRB_IntAttr_Status));

  bool findOptSol = false;
  int optimstatus = SP->get (GRB_IntAttr_Status);
  if (optimstatus == GRB_OPTIMAL)
    {
      cout << "Model is bounded." << endl;
      //printf ("positive pi vars:\n");
      for (i = 0; i < n; i++)
	{
	  // Extract solution
	  piVal[i] = SP->get (GRB_DoubleAttr_X, piVar[i], n);
	  /*
	   for (j = 0; j < n; j++)
	   if (pi[i][j] > 0)
	   printf ("pi[%d,%d] ", i, j, pi[i][j]);
	   */
	}
      findOptSol = true;
      //*subProProfit = SP->getObjective ();
      *subProProfit = SP->get (GRB_DoubleAttr_ObjVal);

    }
  /*
   else if (optimstatus == GRB_INF_OR_UNBD)
   {
   cout << "Model is GRB_INF_OR_UNBD" << endl;
   }*/
  else if (optimstatus == GRB_INFEASIBLE)
    {
      cout << "===> Subproblem is infeasible!" << endl;
      cout << "===> Terminating running." << endl;
      exit (1);
    }
  else if (optimstatus == GRB_UNBOUNDED)
    {
      cout << "Model is unbounded." << endl;
      /*not working
       GRBVar *varVector;
       varVector = SP->getVars ();
       printf ("%s: %f", varVector[0].get (GRB_StringAttr_VarName),
       varVector[0].get (GRB_DoubleAttr_UnbdRay));
       printf ("\n");
       /*
       for (i=0;i<n*n;i++)
       if (varVector[i].UnbdRay != 0)
       printf("%s: %f", varVector[i].VarName,varVector[i].UnbdRay);
       /*
       for (GRBVar v:SP->getVars())
       if (v.get(GRB_DoubleAttr_UnbdRay) != 0)
       printf("%s:%f  ", v.get(GRB_StringAttr_VarName),v.get(GRB_DoubleAttr_UnbdRay));
       */

      for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	  {
	    GRBVar v = SP->getVarByName ("pi_" + itos (i) + "_" + itos (j));
	    double temp = v.get (GRB_DoubleAttr_UnbdRay);
	    unbdPiVal[i][j] = temp;
	    if (PRINT_UNBOUNDED_VECTOR)
	      if (temp != 0)
		printf ("%d->%d %lf  ", i, j, temp);
	  }
      if (PRINT_UNBOUNDED_VECTOR)
	printf ("\n");
    }
  else
    {
      cout << "Optimization was stopped with status = " << optimstatus << endl;
      cout << "===> Terminating running." << endl;
      exit (1);
    }

  delete SP; //if model is defined without pointer, then no deed to delete model

  //printf ("********* end of sub-problem *********\n");
  return findOptSol;
}

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

  int i, j, k;
  int cutCountFea = 1;
  int cutCountUnb = 1;
  double profitInMaster;
  int status;

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

	  wt = new double*[n];
	  dis = new double*[n];
	  for (int i = 0; i < n; i++)
	    {
	      wt[i] = new double[n];
	      dis[i] = new double[n];
	    }
	  for (int i = 1; i <= n; i++)
	    for (int j = 1; j <= n; j++)
	      {
		wt[i - 1][j - 1] = route.w[i][j];
		dis[i - 1][j - 1] = route.d[i][j];
	      }
	  price = route.priceCharged;
	  cost_flow = route.travelCost;
	  cost_vehicle = route.travelCost;
	  vw = route.vehicleWeight;
	  Q = route.totalCapacity;

	  //start calling Gurobi

	  //GRBEnv *env = NULL;

	  GRBVar **x = NULL;
	  GRBVar **y = NULL;
	  GRBVar z;
	  GRBVar s[n];
	  //GRBVar u[n][n][n];
	  //GRBVar theta[n][n];

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

	      UB = 1000000;
	      LB = -1000000;

	      // Create decision variables
	      z = model.addVar (0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "z");
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

	      GRBLinExpr expr9 = z;
	      for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		  {

		    expr9 -= (price - cost_flow) * dis[i][j] * wt[i][j]
			* y[i][j];
		    expr9 += cost_vehicle * vw * dis[i][j] * x[i][j];

		  }
	      model.addConstr (expr9 <= 0, "cut_");

	      /*
	       //set objective
	       GRBLinExpr obj = 0.0;
	       for (i = 0; i < n - 1; i++)
	       for (j = 1; j < n; j++)
	       {
	       obj += (price - cost_flow) * dis[i][j] * wt[i][j] * y[i][j];
	       //obj -= cost * dis[i][j] * theta[i][j];
	       obj -= cost_vehicle * vw * dis[i][j] * x[i][j];
	       }
	       */

	      //constraints to restrict x and y variables
	      for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		  {
		    GRBLinExpr expr9 = 0;
		    GRBLinExpr expr10 = 0;
		    expr9 += y[i][j];
		    expr10 += y[i][j];
		    for (k = 0; k < n; k++)
		      {
			expr9 -= x[i][k];
			expr10 -= x[k][j];
		      }
		    model.addConstr (expr9 <= 0,
				     "x-y-link1_" + itos (i) + "_" + itos (j));
		    model.addConstr (expr10 <= 0,
				     "x-y-link2_" + itos (i) + "_" + itos (j));
		  }

	      GRBLinExpr obj = z;
	      model.setObjective (obj, GRB_MAXIMIZE);
	      //model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
	      printf ("price = %lf\n", price);
	      printf ("cost_flow = %lf\n", cost_flow);
	      printf ("cost_vehicle = %lf\n", cost_vehicle);
	      printf ("vw = %lf\n", vw);

	      // Must set LazyConstraints parameter when using lazy constraints
	      /*
	       model.set (GRB_IntParam_LazyConstraints, 1);
	       model.set (GRB_IntParam_Presolve, 0);
	       model.set (GRB_IntParam_Cuts, 0);
	       model.set (GRB_DoubleParam_Heuristics, 0);
	       model.set (GRB_IntParam_Threads, 1);

	       // Set callback function
	       printIntSol cb = printIntSol (x, y, &z, n);
	       model.setCallback (&cb);
	       */
	      model.set (GRB_IntParam_OutputFlag, 0); //Silent Mode

	      // Optimize model
	      model.optimize ();

	      double **xVal = new double*[n];
	      double **yVal = new double*[n];
	      double **piVal = new double*[n];
	      double **unbdPiVal = new double*[n];
	      for (i = 0; i < n; i++)
		{
		  piVal[i] = new double[n];
		  unbdPiVal[i] = new double[n];
		}

	      while (UB - LB > epsilon)
		{
		  //if (printMore)
		  printf ("========> UB = %lf, LB = %lf\n", UB, LB);
		  if (initialSolution)
		    {
		      for (i = 0; i < n; i++)
			{
			  xVal[i] = new double[n];
			  yVal[i] = new double[n];
			}

		      for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
			  {
			    xVal[i][j] = 0;
			    yVal[i][j] = 0;
			  }

		      xVal[0][n - 1] = 1;
		      yVal[0][n - 1] = 1;
		      /*
		       //the optimal solution for t20_01_data
		       xVal[0][6] = 1;
		       xVal[6][5] = 1;
		       xVal[5][4] = 1;
		       xVal[4][14] = 1;
		       xVal[14][19] = 1;
		       yVal[0][6] = 1;
		       yVal[6][5] = 1;
		       yVal[5][4] = 1;
		       yVal[4][14] = 1;
		       yVal[14][19] = 1;
		       */
		      //hereinafter, this bool var will be false all the time.
		      initialSolution = false;
		    }
		  else
		    {
		      if (model.get (GRB_IntAttr_SolCount) > 0)
			{
			  for (i = 0; i < n; i++)
			    {
			      xVal[i] = model.get (GRB_DoubleAttr_X, x[i], n);
			      yVal[i] = model.get (GRB_DoubleAttr_X, y[i], n);
			    }
			  if (printMore)
			    {
			      printf ("x vars: ");
			      for (i = 0; i < n; i++)
				{
				  xVal[i] = model.get (GRB_DoubleAttr_X, x[i],
						       n);
				  for (j = 0; j < n; j++)
				    if (xVal[i][j] > 0.5)
				      printf ("%d->%d ", i, j);
				}
			      printf ("\ny vars: ");
			      for (i = 0; i < n; i++)
				{
				  yVal[i] = model.get (GRB_DoubleAttr_X, y[i],
						       n);
				  for (j = 0; j < n; j++)
				    if (yVal[i][j] > 0.5)
				      printf ("%d->%d ", i, j);
				}
			      printf ("\n");
			    }

			}
		      // Status checking
		      status = model.get (GRB_IntAttr_Status);
		      if (status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE
			  || status == GRB_UNBOUNDED)
			{
			  cout << "The model cannot be solved "
			      << "because it is infeasible or unbounded"
			      << endl;
			  return 1;
			}
		      if (status != GRB_OPTIMAL)
			{
			  cout << "Optimization was stopped with status "
			      << status << endl;
			  return 1;
			}
		    }

		  //call sub-problem

		  if (DUAL_USE_GRB)
		    {
		      double subPlmObj;
		      bool findOptimalSol = subProblemDual (xVal, yVal, n,
							    piVal, unbdPiVal,
							    &subPlmObj);
		      if (findOptimalSol)
			{

			  for (i = 0; i < n; i++)
			    for (j = 0; j < n; j++)
			      {
				if (piVal[i][j] > 0)
				  printf ("pi[%d,%d]=%lf", i, j, piVal[i][j]);
			      }
			  printf ("\n");

			  GRBLinExpr expr = 0;
			  expr += z;
			  for (i = 0; i < n; i++)
			    for (j = 0; j < n; j++)
			      {
				expr -= (price - cost_flow) * dis[i][j]
				    * wt[i][j] * y[i][j];
				expr += cost_vehicle * vw * dis[i][j] * x[i][j];
				expr -= (Q * x[i][j] - wt[i][j] * y[i][j])
				    * piVal[i][j];
			      }
			  model.addConstr (
			      expr <= 0,
			      "cut_feasible_" + itos (cutCountFea++));
			  if (printMore)
			    printf ("*** fesible sol cut is added ***\n");

			  //update LB
			  profitInMaster = 0;
			  for (i = 0; i < n - 1; i++)
			    for (j = 1; j < n; j++)
			      {
				profitInMaster += (price - cost_flow)
				    * dis[i][j] * wt[i][j] * yVal[i][j];
				profitInMaster -= cost_vehicle * vw * dis[i][j]
				    * xVal[i][j];
			      }
			  printf ("profitInMaster=%lf\n", profitInMaster);
			  printf ("subPlmObj=%lf\n", subPlmObj);
			  if (profitInMaster + subPlmObj > LB)
			    LB = profitInMaster + subPlmObj;
			  //if (profitInMaster > LB)
			  //LB = profitInMaster;
			}
		      else
			{ //the sub-problem is unbounded
			  //find an extreme ray

			  if (printMore)
			    printf ("*** find an extreme ray ***\n");
			  for (i = 0; i < n; i++)
			    for (j = 0; j < n; j++)
			      {
				if (Q * xVal[i][j] - wt[i][j] * yVal[i][j] < 0)
				  {
				    piVal[i][j] = 1;
				  }
				else
				  piVal[i][j] = 0;
			      }
			  GRBLinExpr expr = 0;
			  for (i = 0; i < n; i++)
			    for (j = 0; j < n; j++)
			      expr += (Q * x[i][j] - wt[i][j] * y[i][j])
				  * piVal[i][j];
			  model.addConstr (
			      expr >= 0,
			      "cut_unbounded_" + itos (cutCountUnb++));
			  if (printMore)
			    printf ("*** unbounded sol cut is added ***\n");

			  /*
			   GRBLinExpr expr = 0;
			   for (i = 0; i < n; i++)
			   for (j = 0; j < n; j++)
			   expr += (Q * x[i][j] - wt[i][j] * y[i][j])
			   * unbdPiVal[i][j];
			   model.addConstr (expr >= 0,
			   "cut_unbounded_" + itos (cutCountUnb++));
			   printf ("*** unbounded sol cut is added ***\n");
			   */
			}

		      model.set ("MIPGapAbs", "0.9");
		      model.update ();
		      // Optimize model
		      model.optimize ();
		      //write model to file
		      model.write ("BPMP.lp");
		      //update UB
		      UB = model.get (GRB_DoubleAttr_ObjVal);
		    }
		  else{


		  }

		}

	      printf ("========> Solution found!\n");
	      printf ("========> UB = %lf, LB = %lf\n", UB, LB);

	      //write model to file
	      //model.write ("BPMP.lp");

	      // Extract solution
	      /*
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
	       */
	      //release memory
	      for (i = 0; i < n; i++)
		{
		  delete[] xVal[i];
		  delete[] yVal[i];
		  delete[] piVal[i];
		  delete[] unbdPiVal[i];
		}
	      delete[] xVal;
	      delete[] yVal;
	      delete[] piVal;
	      delete[] unbdPiVal;
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
	  //release memory
	  for (i = 0; i < n; i++)
	    {
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
