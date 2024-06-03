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
#include <time.h>

#include <iomanip>

#include "readData.h"

using namespace std;
using namespace std::chrono;
using namespace std::string_literals;
//to generate vehicle data file name

int numV = 3; //number of vehicles; will be overwritten when reading vehicle file
//number of nodes in instance
//in main(), we also use n to present number of nodes
//since it's easy to read variables with short subscript, like x[n][n][numV]
int NUM_NODES;

double INITIAL_U_COEFFICIENT;

//add cut only in callback function
//if callback function is not used in main()
//this bool has to be set to false, so the constraint can be added to model
//if this bool is set to true, the lazy constraint is only added in callback function
//also, PLEASE CHECK model.set (GRB_IntParam_LazyConstraints, 1)
//GRB_IntParam_LazyConstraints cann't be controled by using bool var
//so if not using lazy constraint, force ADD_UNIQUE_PICKUP_CONSTR_AS_LAZY to false to add
//constraint in model without callback function
//bool ADD_UNIQUE_PICKUP_CONSTR_AS_LAZY = false;
bool PRINT_x_WHEN_INTEGER_SOL = false;
bool PRINT_y_WHEN_INTEGER_SOL = false;
bool SEARCH_K_BEST_SOL = false;
bool USE_MULTISCENARIO = true;
//if USE_LAGRANGIAN_RELAX = true,
//USE_MULTISCENARIO will be forced to be true
//and SEARCH_K_BEST_SOL will be forced to be false
bool USE_LAGRANGIAN_RELAX = true;
//according to p174 on "Integer Programming" by Laurence A. Wolsey 1st edition
bool USE_LR_MULTIPLIER_TYPE_B = false;
bool USE_LR_MULTIPLIER_TYPE_C = true;

bool PRINT_VAR_VALUR = false;
bool PRINT_CONFLICT_PICKUP = true;
bool PRINT_LRmultiplierTypeB_update_process = false;

double bigM = 10000000;

string
itos (int i)
{
  stringstream s;
  s << i;
  return s.str ();
}
bool
isUniquePickup (double***); //this func is not used right now.
void
updateLRmultiplierTypeC (double***, double, double*, double**, double, double);
void
printVar (double***, double***, double****, int*);
void
updateLRmultiplierTypeB (double***, double**, double, int, int);
void
reportTime (clock_t, auto);

//====================================================================================
void
reportTime (clock_t begin, auto beginWallClock)
{
  //clock() gives cpu time on Linux, and wall time on Windows.
  //link: https://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows
  clock_t end = clock ();
  double computerTime = (double) (end - begin);
  double second = computerTime / CLOCKS_PER_SEC;
  printf ("CPU time (on Linux): %lf computer time,  %lf seconds\n",
	  computerTime, second);

  //wallclock time
  auto endWallClock = high_resolution_clock::now ();
  auto elapsedWallClock = duration_cast < std::chrono::nanoseconds
      > (endWallClock - beginWallClock);
  printf ("Wall clock time: %.3f seconds.\n", elapsedWallClock.count () * 1e-9);
}

void
updateLRmultiplierTypeB (double ***soly, double **LR_u, double LR_rou,
			 int countItr, int numArcs)
{
  int i, j, q;
  int n = NUM_NODES;
  double slack[n][n];
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	double sum = 0;
	for (q = 0; q < numV; q++)
	  sum += soly[i][j][q];
	slack[i][j] = 1 - sum;
      }
  double LR_miu = 1.0 / ((double) 10 * numArcs * (countItr + 1));

  //for (i = 0; i < countItr; i++)
  //LR_miu *= LR_rou;

  cout << "updated LR_miu = " << LR_miu << endl;

  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
	{
	  if (PRINT_LRmultiplierTypeB_update_process)
	    printf ("old LR_u_%d_%d = %lf  ", i, j, LR_u[i][j]);

	  LR_u[i][j] = LR_u[i][j] - LR_miu * slack[i][j];

	  if (PRINT_LRmultiplierTypeB_update_process)
	    {
	      printf ("slack_%d_%d = %lf  ", i, j, slack[i][j]);
	      printf ("after updating: LR_u_%d_%d = %lf  ", i, j, LR_u[i][j]);
	    }

	  if (LR_u[i][j] < 0)
	    LR_u[i][j] = 0;
	}
      if (PRINT_LRmultiplierTypeB_update_process)
	{
	  cout << endl;
	}
    }

}

void
printVar (double ***solx, double ***soly, double ****solu, int *origin)
{
  int i, j, k, q;
  int n = NUM_NODES;
  cout << "SOLUTION:" << endl;

  for (q = 0; q < numV; q++)
    {
      cout
	  << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
	  << endl;
      printf ("Vehicle %d starting at node %d:\n", q + 1, origin[q] + 1);

      printf ("x=1:\n");
      for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	  if (solx[i][j][q] > 0.99 && solx[i][j][q] < 1.01)
	    printf ("(%d,%d) ", i + 1, j + 1);
      printf ("\n");

      printf ("y=1:\n");
      for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	  if (soly[i][j][q] > 0.99 && soly[i][j][q] < 1.01)
	    printf ("(%d,%d) ", i + 1, j + 1);

      printf ("\n");
      cout << "u>0.000001:" << endl;
      for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	  for (k = 0; k < n; k++)
	    if (solu[i][j][k][q] > 0.000001)
	      printf ("(%d,%d,%d) ", i + 1, j + 1, k + 1);
      printf ("\n");
    }
  cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
      << endl;
}

bool
isUniquePickup (double ***soly)
{
  bool notFindConflict = true;
  int i, j, q1, q2;
  for (i = 0; i < NUM_NODES; i++)
    {
      for (j = 0; j < NUM_NODES; j++)
	{
	  for (q1 = 0; q1 < numV - 1; q1++)
	    {

	      for (q2 = q1 + 1; q2 < numV; q2++)
		{
		  //printf("%d %d (v%d) %lf",i,j,q1,soly[i][j][q1]);
		  //printf("%d %d (v%d) %lf\n",i,j,q1,soly[i][j][q2]);
		  if (soly[i][j][q1] + soly[i][j][q2] > 1.8)
		    {
		      printf ("find conflict for %d->%d by v%d and v%d", i + 1,
			      j + 1, q1 + 1, q2 + 1);
		      notFindConflict = false;
		      break;
		    }
		}
	      if (!notFindConflict)
		break;
	    }
	  if (!notFindConflict)
	    break;
	}
      if (!notFindConflict)
	break;
    }

  return notFindConflict;
}

void
updateLRmultiplierTypeC (double ***soly, double LR_lamda, double *LR_miu,
			 double **LR_u, double LB, double totalProfit)
{
  int i, j, q;
  //update LR multiplier LR_u
  double norm = 0;
  int n = NUM_NODES;
  double slack[n][n];
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	double sum = 0;
	for (q = 0; q < numV; q++)
	  sum += soly[i][j][q];
	slack[i][j] = 1 - sum;
	//printf("slack %d %d = %lf\n",i,j,slack[i][j]);
	norm += slack[i][j] * slack[i][j];
      }
  cout << "norm = " << norm << endl;

  *LR_miu = LR_lamda * (totalProfit - LB) / norm;
  cout << "updated LR_miu = " << *LR_miu << endl;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	LR_u[i][j] = LR_u[i][j] - (*LR_miu) * slack[i][j];
	if (LR_u[i][j] < 0)
	  LR_u[i][j] = 0;
      }
}

int
main (int argc, char *argv[])
{

  if (argc != 4)
    {
      //cout << "Usage: ./BPMP.x dataFolder" << endl;
      cout
	  << "Usage: ./mvbpmp_lr.x nodesDataNameAndPath numberOfVehicles vehicleDataNameAndPath"
	  << endl;
      return 1;
    }
  else
    {
      cout << "Nodes Data: " << argv[1] << endl;
      cout << "Number of Vehicles: " << argv[2] << endl;
      cout << "Vehicles Data: " << argv[3] << endl;
      numV = stoi(argv[2]);
    }

  if (USE_MULTISCENARIO && SEARCH_K_BEST_SOL)
    {
      printf (
	  "Can't do both USE_MULTISCENARIO and SEARCH_K_BEST_SOL right now!\n");
      printf ("Exit running.");
      exit (1);

      if (USE_LR_MULTIPLIER_TYPE_B && USE_LR_MULTIPLIER_TYPE_C)
	{
	  printf ("only one of the LR multiplier type can be selected.\n");
	  printf ("exit.\n");
	  exit (1);
	}

      if (!USE_LR_MULTIPLIER_TYPE_B && !USE_LR_MULTIPLIER_TYPE_C)
	{
	  printf ("one of the LR multiplier type must be selected.\n");
	  printf ("exit.\n");
	  exit (1);
	}
    }

  if (USE_LAGRANGIAN_RELAX)
    {
      USE_MULTISCENARIO = true;
      SEARCH_K_BEST_SOL = false;
      //when this is true, the constraint won't be added to model directly
      //when callback function is used, this will be added in callback
      //since we will relax it using Lagrangian Relaxation, set it as true
      //ADD_UNIQUE_PICKUP_CONSTR_AS_LAZY = true;
      cout
	  << "Since USE_LAGRANGIAN_RELAX = true, USE_MULTISCENARIO is forced to be true, "
	      "SEARCH_K_BEST_SOL is forced to be false\n" << endl;
    }

  clock_t beginTime, endTimeOfLastIteration;

  auto beginWallClock = high_resolution_clock::now ();
  auto endTimeOfLastIterationWallClock = high_resolution_clock::now ();
  readData route;
  //auto it = filesystem::directory_iterator (argv[1]);
  //for (const auto &entry : it)
  //{
  // filesystem::path path = entry.path ();

  //filesystem::path path = argv[1];

  //if (entry.is_regular_file ())
  // if (path.is_regular_file ())
  //{

  //string filename = path.string ();
  string filename = argv[1];
  printf ("reading file %s ...\n", filename.c_str ());
  route.readSingleFile (filename);

  //route.printStats();
  //route.printData();

  int n = route.numOfNode;
  NUM_NODES = route.numOfNode;

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

  /*
   string vehicleFileName = "vehicle_location_data/vehicle_data_"
   + to_string (n) + "_" + to_string (numV) + "_1.txt";
   */
  string vehicleFileName = argv[3];
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
      beginTime = clock ();

      env = new GRBEnv ();
      GRBModel model = GRBModel (*env);

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
		      0.0, 1.0, 0, GRB_BINARY,
		      "x_" + itos (i) + "_" + itos (j) + "_" + itos (q));
		  y[i][j][q] = model.addVar (
		      0.0, 1.0, 0, GRB_BINARY,
		      "y_" + itos (i) + "_" + itos (j) + "_" + itos (q));
		  theta[i][j][q] = model.addVar (
		      0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
		      "theta_" + itos (i) + "_" + itos (j) + "_" + itos (q));

		  for (k = 0; k < n; k++)
		    {
		      string s = "u_" + itos (i) + "_" + itos (j) + "_"
			  + itos (k) + "_" + itos (q);
		      u[i][j][k][q] = model.addVar (0.0, GRB_INFINITY, 0.0,
						    GRB_CONTINUOUS, s);
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
      GRBConstr *vehOriginConstr = 0;
      GRBConstr *vehDestConstr = 0;
      vehOriginConstr = new GRBConstr[numV];
      vehDestConstr = new GRBConstr[numV];

      GRBConstr **uniquePickup = new GRBConstr*[n];
      for (i = 0; i < n; i++)
	uniquePickup[i] = new GRBConstr[n];

      for (q = 0; q < numV; q++)
	{
	  int ogn = origin[q];

	  //vehicle goes out of origins
	  GRBLinExpr expr1 = 0.0;
	  for (i = 0; i < n; i++)
	    expr1 += x[ogn][i][q];

	  //model.addConstr (expr1 == 1, "origin_" + itos (q));
	  vehOriginConstr[q] = model.addConstr (expr1 == 1,
						"origin_" + itos (q));

	  //vehicle goes back to node n
	  GRBLinExpr expr2 = 0.0;
	  for (i = 0; i < n - 1; i++)
	    expr2 += x[i][n - 1][q];
	  //model.addConstr (expr2 == 1, "destination_" + itos (q));
	  vehDestConstr[q] = model.addConstr (expr2 == 1,
					      "destination_" + itos (q));

	  //flow conservation
	  for (int k = 0; k < n - 1; k++)
	    {
	      if (k != ogn)
		{
		  GRBLinExpr expr = 0;
		  for (i = 0; i < n - 1; i++)
		    expr += x[i][k][q];
		  //BE CAREFUL!
		  //I used j=1 to start which exclues node 0
		  //which cause that the optimal profit is lower!
		  for (j = 0; j < n; j++)
		    expr -= x[k][j][q];
		  model.addConstr (
		      expr == 0,
		      "flow_conservation_" + itos (k) + "_" + itos (q));
		}
	    }

	  //distance
	  GRBLinExpr expr3 = 0.0;
	  for (i = 0; i < n; i++)
	    for (j = 0; j < n; j++)
	      expr3 += dis[i][j] * x[i][j][q];
	  model.addConstr (expr3 <= route.DIS, "distance_" + itos (q));

	  //node degree less than 1
	  for (int j = 0; j < n; j++)
	    {
	      GRBLinExpr expr = 0.0;
	      for (i = 0; i < n; i++)
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
		    "subtour_" + itos (i) + "_" + itos (j) + "_" + itos (q));
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
		      expr += u[i][k][j][q] + u[k][j][i][q] - u[i][j][k][q];
		    else
		      expr += u[k][j][i][q];
		  }
		//when k==n-1
		if (j != n - 1)
		  expr += u[i][n - 1][j][q];
		model.addConstr (
		    expr == 0,
		    "flow_" + itos (i) + "_" + itos (j) + "_" + itos (q));
	      }

	  //arc flow upperbound
	  for (i = 0; i < n - 1; i++)
	    for (j = 0; j < n; j++)
	      {
		GRBLinExpr expr = 0.0;
		expr += theta[i][j][q] - Q * x[i][j][q];
		model.addConstr (
		    expr <= 0,
		    "flowBound_" + itos (i) + "_" + itos (j) + "_" + itos (q));
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
      //if (!ADD_UNIQUE_PICKUP_CONSTR_AS_LAZY)

      for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	  {
	    GRBLinExpr expr = 0.0;
	    for (q = 0; q < numV; q++)
	      expr += y[i][j][q];
	    //model.addConstr (expr <= 1,"one-one_" + itos (i) + "_" + itos (j));

	    uniquePickup[i][j] = model.addConstr (
		expr <= 1, "uniquePickup_" + itos (i) + "_" + itos (j));
	  }

      double yCoefficient[n][n];

      for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	  yCoefficient[i][j] = price * dis[i][j] * wt[i][j];

      //set up objective
      GRBLinExpr obj = 0.0;
      for (q = 0; q < numV; q++)
	{
	  for (i = 0; i < n - 1; i++)
	    for (j = 0; j < n; j++)
	      {
		obj += yCoefficient[i][j] * y[i][j][q];
		obj -= cost * dis[i][j] * theta[i][j][q];
		obj -= cost * vw * dis[i][j] * x[i][j][q];
	      }
	}

      model.setObjective (obj, GRB_MAXIMIZE);
      //model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
      //printf ("price = %lf\n", price);
      //printf ("cost = %lf\n", cost);
      //printf ("vw = %lf\n", vw);

      model.set (GRB_IntParam_OutputFlag, 0);

      /*
       //find dual value by solving relaxed model
       //most dual vaue are zeros

       for (q = 0; q < numV; q++)
       {
       for (i = 0; i < n; i++)
       for (j = 0; j < n; j++)
       {
       x[i][j][q].set (GRB_CharAttr_VType, GRB_CONTINUOUS);
       y[i][j][q].set (GRB_CharAttr_VType, GRB_CONTINUOUS);
       }
       }

       model.optimize ();
       int optimstatus = model.get (GRB_IntAttr_Status);
       if (optimstatus == GRB_OPTIMAL)
       {
       for (i = 1; i < n; i++)
       for (j = 1; j < n; j++)
       cout << uniquePickup[i][j].get (GRB_DoubleAttr_Pi)
       << endl;
       }
       exit (1);

       */
      /*
       GRBModel relax = model.relax ();
       relax.optimize ();
       */

      /*
       //========== solve the relaxed model to get dual value of uniquePickup constraint
       //but it's not working right now
       *
       GRBModel relax = model.relax ();
       // Optimize model
       relax.optimize ();

       double dualValue[n][n];
       for (i = 1; i < n; i++)
       for (j = 1; j < n; j++)
       cout << uniquePickup[i][j].get (GRB_DoubleAttr_Pi) << endl;

       int optimstatus = relax.get (GRB_IntAttr_Status);
       if (optimstatus == GRB_OPTIMAL)
       {
       cout << "Find optimal soltuion." << endl;
       for (i = 1; i < n; i++)
       for (j = 1; j < n; j++)
       cout << uniquePickup[i][j].get (GRB_DoubleAttr_Pi)
       << endl;

       //for (i = 0; i < n; i++)
       //for (j = 0; j < n; j++)
       //printf ("%d %d dual=%lf\n", i, j, dualValue[i][j]);
       }
       if (optimstatus == GRB_INFEASIBLE)
       {
       cout << "===> Relaxed problem is infeasible!" << endl;
       cout << "===> Terminating running." << endl;
       exit (1);
       }
       else if (optimstatus == GRB_UNBOUNDED)
       {
       cout << "Model is unbounded" << endl;
       cout << "===> Terminating running." << endl;
       exit (1);
       }
       */

      //**********the above is the original MVBPMP model**********
      //**********the lagrangian relaxed obj will be added later**********
      if (USE_LAGRANGIAN_RELAX)
	{
	  double UB = bigM;
	  double LB = -bigM;

	  double LR_gap_tolerance = 0.05;
	  double LR_complementarity_tolerance = 0.0001;
	  double LR_min_lamda = 0.0001;
	  double LR_miu_tolerance = 0.00000001;

	  double LR_rou = 0.4; //for LR multiplier type b

	  double LR_lamda = 2;
	  int LR_maxNumNoImprovement = 3;
	  int LR_numNoImprovement = 0;

	  //double LR_u[n][n];	      //the lagrangian multiplier
	  double **LR_u = new double*[n];
	  for (i = 0; i < n; i++)
	    LR_u[i] = new double[n];

	  //initialize LR_u
	  double numArcs = n * (n - 1);
	  for (i = 0; i < n; i++)
	    for (j = 0; j < n; j++)
	      LR_u[i][j] = 1 / numArcs;

	  double LR_lowestUpperBound = bigM;

	  //store the solution of MVBPMP_LR in each iteration
	  //double solx[n][n][numV];
	  //double soly[n][n][numV];
	  //double solu[n][n][n][numV];
	  double ***solx_d = new double**[n];	    //solution x for LR dual
	  double ***soly_d = new double**[n];
	  double ****solu_d = new double***[n];
	  double soltheta_d[n][n][numV];
	  for (i = 0; i < n; i++)
	    {
	      solx_d[i] = new double*[n];
	      soly_d[i] = new double*[n];
	      solu_d[i] = new double**[n];
	      for (j = 0; j < n; j++)
		{
		  solx_d[i][j] = new double[numV];
		  soly_d[i][j] = new double[numV];
		  solu_d[i][j] = new double*[n];
		  for (k = 0; k < n; k++)
		    solu_d[i][j][k] = new double[numV];
		}
	    }

	  int maxItr = 3;
	  int countItr = 0;

	  cout << "===> LR_gap_tolerance = " << LR_gap_tolerance << endl;

	  endTimeOfLastIteration = clock ();
	  endTimeOfLastIterationWallClock = high_resolution_clock::now ();
	  while ((UB - LB) / UB > LR_gap_tolerance)
	    {

	      countItr++;
	      //stop after a certain iterations
	      //if (countItr > maxItr)break;

	      cout << "============= Itr " << countItr << " =============="
		  << endl;
	      cout << "LB = " << LB << endl;
	      cout << "UB = " << UB << endl;

	      printf ("===> start multiscenarios\n");
	      //model.set (GRB_StringAttr_ModelName, "multiscenario");

	      //after constructing the base model, add scenarios
	      //Scenario 0: the base model
	      //scenario 1: the seperated BPMP for vehicle 1
	      //scenario 2: the seperated BPMP for vehicle 2
	      //scenario 3: the seperated BPMP for vehicle 3 (only for 3 vehicles instances)

	      model.set (GRB_IntAttr_NumScenarios, numV);

	      // Scenario 0: Base model, hence, nothing to do except giving the
	      //             scenario a name
	      //model.set (GRB_IntParam_ScenarioNumber, 0);
	      //model.set (GRB_StringAttr_ScenNName, "Base model");

	      for (int s = 0; s < numV; s++)
		{
		  //======== Scenario s ========

		  model.set (GRB_IntParam_ScenarioNumber, s);

		  string scenarioName = "BPMP_vehicle_" + itos (s + 1);
		  model.set (GRB_StringAttr_ScenNName, scenarioName.c_str ());

		  //======== preset all y vars LB and UB =============
		  for (q = 0; q < numV; q++)
		    for (i = 0; i < n; i++)
		      for (j = 0; j < n; j++)
			{
			  x[i][j][q].set (GRB_DoubleAttr_ScenNLB, 0.0);
			  x[i][j][q].set (GRB_DoubleAttr_ScenNUB, 1.0);
			  y[i][j][q].set (GRB_DoubleAttr_ScenNLB, 0.0);
			  y[i][j][q].set (GRB_DoubleAttr_ScenNUB, 1.0);

			  for (k = 0; k < n; k++)
			    {
			      u[i][j][k][q].set (GRB_DoubleAttr_ScenNLB, 0.0);
			      u[i][j][k][q].set (GRB_DoubleAttr_ScenNUB, 1.0);
			    }
			}
		  for (q = 0; q < numV; q++)
		    {
		      int ogn = origin[q];
		      for (i = 0; i < n; i++)
			{
			  x[i][i][q].set (GRB_DoubleAttr_ScenNUB, 0);
			  y[i][i][q].set (GRB_DoubleAttr_ScenNUB, 0);
			  x[i][ogn][q].set (GRB_DoubleAttr_ScenNUB, 0);
			  y[i][ogn][q].set (GRB_DoubleAttr_ScenNUB, 0);
			  x[n - 1][i][q].set (GRB_DoubleAttr_ScenNUB, 0);
			  y[n - 1][i][q].set (GRB_DoubleAttr_ScenNUB, 0);
			  for (j = 0; j < n; j++)
			    if (wt[i][j] == 0)
			      y[i][j][q].set (GRB_DoubleAttr_ScenNUB, 0);
			}
		    }

		  //============ relax unique pickup constraint ==============
		  for (i = 0; i < n; i++)
		    for (j = 0; j < n; j++)
		      uniquePickup[i][j].set (GRB_DoubleAttr_ScenNRHS,
					      GRB_INFINITY);

		  for (int q = 0; q < numV; q++)
		    {
		      //============ set up LR multiplier in obj ==============
		      for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
			  //update coefficient in obj
			  y[i][j][q].set (GRB_DoubleAttr_ScenNObj,
					  yCoefficient[i][j] - LR_u[i][j]);

		      //=== set up RHS for org and dest for current vehicle ===
		      if (q == s)
			{
			  vehOriginConstr[q].set (GRB_DoubleAttr_ScenNRHS, 1);
			  vehDestConstr[q].set (GRB_DoubleAttr_ScenNRHS, 1);
			}
		      else
			{
			  //===> make other two vehicles satisfy constraint
			  //change the RHS for other vehicles from 1 to 0
			  vehOriginConstr[q].set (GRB_DoubleAttr_ScenNRHS, 0);
			  vehDestConstr[q].set (GRB_DoubleAttr_ScenNRHS, 0);
			  //============ set up x, y, u var for other vehicles  ==============
			  for (i = 0; i < n; i++)
			    for (j = 0; j < n; j++)
			      {
				x[i][j][q].set (GRB_DoubleAttr_ScenNUB, 0.0);
				y[i][j][q].set (GRB_DoubleAttr_ScenNUB, 0.0);
				for (k = 0; k < n; k++)
				  u[i][j][k][q].set (GRB_DoubleAttr_ScenNUB,
						     0.0);
			      }
			}
		    }
		}

	      // Optimize model
	      model.optimize ();

	      model.write ("multiscenario.lp");

	      double scenObj[numV];

	      int nScenarios = model.get (GRB_IntAttr_NumScenarios);

	      if (numV != nScenarios)
		{
		  cout
		      << "The number of vehicles is not equal to the number of scenarios!"
		      << endl;
		  cout << "Exit running!" << endl;
		  exit (1);
		}

	      for (int s = 0; s < nScenarios; s++)
		{
		  int modelSense = GRB_MAXIMIZE;

		  // Set the scenario number to query the information for this scenario
		  model.set (GRB_IntParam_ScenarioNumber, s);

		  // collect result for the scenario
		  double scenNObjBound = model.get (
		      GRB_DoubleAttr_ScenNObjBound);
		  double scenNObjVal = model.get (GRB_DoubleAttr_ScenNObjVal);

		  cout << "------ Scenario " << s << " ("
		      << model.get (GRB_StringAttr_ScenNName) << ")" << endl;

		  // Check if we found a feasible solution for this scenario
		  if (modelSense * scenNObjVal >= GRB_INFINITY)
		    if (modelSense * scenNObjBound >= GRB_INFINITY)
		      {
			// Scenario was proven to be infeasible
			cout << endl << "INFEASIBLE" << endl;
			exit (1);
		      }

		    else
		      {
			// We did not find any feasible solution - should not happen in
			// this case, because we did not set any limit (like a time
			// limit) on the optimization process
			cout << endl << "NO SOLUTION" << endl;
			exit (1);
		      }
		  else
		    {
		      cout << "OBJ VALUE: " << scenNObjVal << endl;

		      scenObj[s] = scenNObjVal;

		      //x, y, u vars can only be positive when the vehicle
		      //index is equato to scene index according to
		      //the settting in scenarios
		      double scenNX;

		      for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
			  {
			    scenNX = x[i][j][s].get (GRB_DoubleAttr_ScenNX);
			    solx_d[i][j][s] = scenNX;

			    scenNX = y[i][j][s].get (GRB_DoubleAttr_ScenNX);
			    soly_d[i][j][s] = scenNX;

			    scenNX = theta[i][j][s].get (GRB_DoubleAttr_ScenNX);
			    soltheta_d[i][j][s] = scenNX;

			    for (k = 0; k < n; k++)
			      {
				scenNX = u[i][j][k][s].get (
				    GRB_DoubleAttr_ScenNX);
				solu_d[i][j][k][s] = scenNX;
			      }
			  }

		      if (PRINT_VAR_VALUR)
			printVar (solx_d, soly_d, solu_d, origin);
		    }
		}
	      //collect profit from all scenarios
	      //each scenario represent one vehicle's schedule
	      double totalProfit = 0;
	      for (i = 0; i < numV; i++)
		totalProfit += scenObj[i];
	      cout << "===> sum of OBJ VALUE: " << totalProfit << endl;

	      //add the lagrangian multiplier part sum_LR_u back to obj
	      //since this is constant and not considered in scenarios
	      //sum_LR_u_times_y is for printing here and future use when solution
	      //is feasible for the original problem
	      double sum_LR_u = 0;
	      double sum_LR_u_times_y = 0;
	      bool findConflictPickup = false;
	      for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		  {
		    sum_LR_u += LR_u[i][j];

		    double sum = 0;
		    for (q = 0; q < numV; q++)
		      sum += soly_d[i][j][q];
		    sum_LR_u_times_y += LR_u[i][j] * sum;

		    if (sum > 1.9)
		      {
			findConflictPickup = true;
			if (PRINT_CONFLICT_PICKUP)
			  printf ("===> conflict %d -> %d\n", i + 1, j + 1);
		      }
		  }

	      cout << "sum_LR_u = " << sum_LR_u << endl;
	      cout << "sum_LR_u_times_y = " << sum_LR_u_times_y << endl;

	      totalProfit += sum_LR_u;
	      cout << "totalProfit in obj including LR = " << totalProfit
		  << endl;

	      double profitInLRdual = totalProfit;

	      if (profitInLRdual < UB)
		UB = profitInLRdual;

	      cout << "===> UB = " << UB << endl;

	      if (UB >= LR_lowestUpperBound)
		LR_numNoImprovement++;
	      else
		{
		  LR_numNoImprovement = 0;
		  LR_lowestUpperBound = UB;
		  cout << "find lower UB" << endl;
		}

	      printf ("===> LR_numNoImprovement = %d\n", LR_numNoImprovement);

	      //===> update LR parameters for the next iteration
	      if (LR_numNoImprovement >= LR_maxNumNoImprovement)
		{
		  LR_lamda = 0.5 * LR_lamda;
		  LR_numNoImprovement = 0;
		}

	      cout << "===> LR_lamda = " << LR_lamda << endl;

	      if (LR_lamda <= LR_min_lamda)
		{
		  cout << "===> LR_lamda is less than LR_min_damda "
		      << LR_min_lamda << endl;
		  cout << "end loop";
		  break;
		}

	      //check if there is request conflict
	      //bool notFindConflict = isUniquePickup (soly_d);

	      //===> if there is no request conflict
	      //===> use the current solution without LR item as LB
	      if (!findConflictPickup)
		{
		  cout << "===> There is no request conflict in solution :)"
		      << endl;

		  //check if the LR_item is close to zero
		  //if yes, then this is the optimal solution
		  double LR_item = sum_LR_u - sum_LR_u_times_y;
		  cout << "===> the LR items added in obj = " << LR_item
		      << endl;

		  if ((LR_item >= 0 && LR_item <= LR_complementarity_tolerance)
		      || (LR_item < 0
			  && LR_item >= -LR_complementarity_tolerance))
		    {
		      LB = profitInLRdual - LR_item;
		      cout << "===> the abs(LR item) added in obj is less than "
			  << LR_complementarity_tolerance << endl;
		      cout << "===> the optimal solution is found :)" << endl;
		      cout << "===> end loop." << endl;

		      reportTime (endTimeOfLastIteration,
				  endTimeOfLastIterationWallClock);

		      printVar (solx_d, soly_d, solu_d, origin);

		      //A STORY:
		      //the following is the solution for t10_09_data
		      //at first my optimal obj is sometimes less than the correct optimum
		      //for example, in t10_09_data. After feeding the correct solution below
		      //into the model, the model is infeasible.
		      //After debug, I found that in flow_conservation constraint
		      //I set up (j=1; j<n; j++), so j=0 is missed, which cause the bug.
		      /*
		       //=========check profit=========

		       for (i = 0; i < n; i++)
		       for (j = 0; j < n; j++)
		       yCoefficient[i][j] = price * dis[i][j] * wt[i][j];

		       //set up objective
		       double objcheck = 0.0;
		       for (q = 0; q < numV; q++)
		       {
		       for (i = 0; i < n; i++)
		       for (j = 0; j < n; j++)
		       {
		       objcheck += yCoefficient[i][j]
		       * soly_d[i][j][q];
		       objcheck -= cost * dis[i][j]
		       * soltheta_d[i][j][q];
		       objcheck -= cost * vw * dis[i][j]
		       * solx_d[i][j][q];
		       }
		       }
		       cout << "recaculated profit = " << objcheck << endl;

		       //below is the solution from t10_09, 3 vehicles
		       for (q = 0; q < numV; q++)
		       for (i = 0; i < n; i++)
		       for (j = 0; j < n; j++)
		       {
		       solx_d[i][j][q] = 0;
		       soly_d[i][j][q] = 0;
		       soltheta_d[i][j][q] = 0;
		       }
		       solx_d[0][3][0] = 1;
		       solx_d[2][9][0] = 1;
		       solx_d[3][2][0] = 1;

		       soly_d[0][2][0] = 1;
		       soly_d[2][9][0] = 1;
		       soly_d[3][2][0] = 1;

		       solx_d[0][9][1] = 1;
		       solx_d[3][0][1] = 1;
		       solx_d[5][3][1] = 1;

		       soly_d[0][9][1] = 1;
		       soly_d[3][0][1] = 1;
		       soly_d[5][3][1] = 1;

		       solx_d[1][7][2] = 1;
		       solx_d[2][4][2] = 1;
		       solx_d[4][1][2] = 1;
		       solx_d[7][9][2] = 1;
		       soly_d[2][1][2] = 1;
		       soly_d[2][4][2] = 1;
		       soly_d[2][7][2] = 1;
		       soly_d[4][1][2] = 1;
		       soly_d[4][7][2] = 1;
		       soly_d[7][9][2] = 1;

		       soltheta_d[0][3][0] = 0.7;
		       soltheta_d[2][9][0] = 0.9;
		       soltheta_d[3][2][0] = 0.9;

		       soltheta_d[0][9][1] = 0.9;
		       soltheta_d[3][0][1] = 0.7;
		       soltheta_d[5][3][1] = 1;

		       soltheta_d[1][7][2] = 0.6;
		       soltheta_d[2][4][2] = 0.8;
		       soltheta_d[4][1][2] = 0.8;
		       soltheta_d[7][9][2] = 1;

		       objcheck = 0.0;
		       for (q = 0; q < numV; q++)
		       {
		       for (i = 0; i < n; i++)
		       for (j = 0; j < n; j++)
		       {
		       objcheck += yCoefficient[i][j]
		       * soly_d[i][j][q];
		       objcheck -= cost * dis[i][j]
		       * soltheta_d[i][j][q];
		       objcheck -= cost * vw * dis[i][j]
		       * solx_d[i][j][q];
		       }
		       }
		       cout << "recaculated profit = " << objcheck << endl;

		       //=========end=============
		       */
		      break;
		    }

		  //if LR_item is still big, then we need to keep iteration
		  cout << "===> Use this solution to calculate LB." << endl;

		  //since no request conflict
		  //so the solution for x, y, u are feasible for the original model
		  //but - u_(i,j)*( 1-sum_q y(i,j,q) ) was added to obj
		  //so we need to recalculte total profit
		  double profitForLB = profitInLRdual;
		  profitForLB += sum_LR_u_times_y;
		  profitForLB -= sum_LR_u;

		  if (profitForLB > LB)
		    {
		      LB = profitForLB;
		      cout << "===> find a better LB = " << LB << endl;
		    }

		  if (USE_LR_MULTIPLIER_TYPE_B)
		    updateLRmultiplierTypeB (soly_d, LR_u, LR_rou, countItr,
					     numArcs);

		  double LR_miu;
		  if (USE_LR_MULTIPLIER_TYPE_C)
		    updateLRmultiplierTypeC (soly_d, LR_lamda, &LR_miu, LR_u,
					     LB, profitInLRdual);

		  reportTime (endTimeOfLastIteration,
			      endTimeOfLastIterationWallClock);
		  endTimeOfLastIteration = clock ();
		  endTimeOfLastIterationWallClock =
		      high_resolution_clock::now ();

		  if ((UB - LB) / UB <= LR_gap_tolerance)
		    cout << "=== will stop loop because gap " << (UB - LB) / UB
			<< " < LR_gap_toleranc " << LR_gap_tolerance
			<< endl;

		  if (LR_miu <= LR_miu_tolerance)
		    {
		      cout
			  << "=== will stop loop because LR_miu < LR_miu_tolerance "
			  << LR_miu_tolerance << endl;
		      break;
		    }
		  else
		    continue; //jump to the next iteration of while loop
		}
	      //if there is request conflict, then call gurobi to find LB

	      //========================== calculate LB ========================//

	      //1. If we found requests picked up by multiple vehicles
	      //for example, r_(i,j) picked up by multiple vehicles
	      //then we let y_(i,j) undecided
	      //2. If r_(i,j) is picked up by only one vehicle, let y_(i,j)=1
	      //3. force other y_(i,j)=0
	      //in this way, we try to reallocate the request

	      //model.set (GRB_StringAttr_ModelName, "multiscenario_LB");
	      model.set (GRB_IntAttr_NumScenarios, 1);
	      model.set (GRB_IntParam_ScenarioNumber, 0);

	      model.set (GRB_StringAttr_ScenNName, "BPMP_LB");

	      //add unique pickup constraint back
	      for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		  uniquePickup[i][j].set (GRB_DoubleAttr_ScenNRHS, 1.0);

	      //preset all y vars LB and UB
	      for (q = 0; q < numV; q++)
		for (i = 0; i < n; i++)
		  for (j = 0; j < n; j++)
		    {
		      /*
		       if (soly_d[i][j][q] > 0.5)
		       {
		       y[i][j][q].set (GRB_DoubleAttr_ScenNLB, 1.0);
		       y[i][j][q].set (GRB_DoubleAttr_ScenNUB, 1.0);
		       }
		       else
		       */
		      if (soly_d[i][j][q] < 0.5)
			{
			  y[i][j][q].set (GRB_DoubleAttr_ScenNLB, 0.0);
			  y[i][j][q].set (GRB_DoubleAttr_ScenNUB, 0.0);
			}
		      else
			{
			  y[i][j][q].set (GRB_DoubleAttr_ScenNLB, 0.0);
			  y[i][j][q].set (GRB_DoubleAttr_ScenNUB, 1.0);
			}
		      x[i][j][q].set (GRB_DoubleAttr_ScenNLB, 0.0);
		      x[i][j][q].set (GRB_DoubleAttr_ScenNUB, 1.0);
		      for (k = 0; k < n; k++)
			{
			  u[i][j][k][q].set (GRB_DoubleAttr_ScenNLB, 0.0);
			  u[i][j][k][q].set (GRB_DoubleAttr_ScenNUB, 1.0);
			}
		    }

	      /*
	       //if there is conflicting requests (i,j) then let y_(i,j) undecided
	       for (i = 0; i < n; i++)
	       for (j = 0; j < n; j++)
	       {
	       for (int q1 = 0; q1 < numV - 1; q1++)
	       for (int q2 = q1 + 1; q2 < numV; q2++)
	       {
	       if (soly_d[i][j][q1] + soly_d[i][j][q2] > 1.9)
	       {
	       y[i][j][q1].set (GRB_DoubleAttr_ScenNLB,
	       0.0);
	       y[i][j][q1].set (GRB_DoubleAttr_ScenNUB,
	       1.0);
	       y[i][j][q2].set (GRB_DoubleAttr_ScenNLB,
	       0.0);
	       y[i][j][q2].set (GRB_DoubleAttr_ScenNUB,
	       1.0);
	       }
	       }
	       }
	       */
	      //set obj parameters, remove LR multiplier
	      for (q = 0; q < numV; q++)
		{
		  for (i = 0; i < n; i++)
		    for (j = 0; j < n; j++)
		      y[i][j][q].set (GRB_DoubleAttr_ScenNObj,
				      yCoefficient[i][j]);

		  vehOriginConstr[q].set (GRB_DoubleAttr_ScenNRHS, 1);
		  vehDestConstr[q].set (GRB_DoubleAttr_ScenNRHS, 1);
		}

	      model.optimize ();

	      nScenarios = model.get (GRB_IntAttr_NumScenarios);
	      for (int s = 0; s < nScenarios; s++)
		{
		  int modelSense = GRB_MAXIMIZE;
		  // Set the scenario number to query the information for this scenario
		  model.set (GRB_IntParam_ScenarioNumber, s);

		  // collect result for the scenario
		  double scenNObjBound = model.get (
		      GRB_DoubleAttr_ScenNObjBound);
		  double scenNObjVal = model.get (GRB_DoubleAttr_ScenNObjVal);

		  cout << "------ Scenario " << s << " ("
		      << model.get (GRB_StringAttr_ScenNName) << ")" << endl;

		  // Check if we found a feasible solution for this scenario
		  if (modelSense * scenNObjVal >= GRB_INFINITY)
		    if (modelSense * scenNObjBound >= GRB_INFINITY)
		      {
			cout << endl << "INFEASIBLE" << endl;
			exit (1);
		      }

		    else
		      {
			cout << endl << "NO SOLUTION" << endl;
			exit (1);
		      }
		  else
		    {
		      cout << "OBJ: " << scenNObjVal << endl;
		      //update LB
		      if (scenNObjVal > LB)
			LB = scenNObjVal;
		      cout << "===> LB = " << LB << endl;
		      double scenNX;

		      double ***solx = new double**[n];	//solution x for LR dual
		      double ***soly = new double**[n];
		      double ****solu = new double***[n];
		      for (i = 0; i < n; i++)
			{
			  solx[i] = new double*[n];
			  soly[i] = new double*[n];
			  solu[i] = new double**[n];
			  for (j = 0; j < n; j++)
			    {
			      solx[i][j] = new double[numV];
			      soly[i][j] = new double[numV];
			      solu[i][j] = new double*[n];
			      for (k = 0; k < n; k++)
				solu[i][j][k] = new double[numV];
			    }
			}

		      for (q = 0; q < numV; q++)
			for (i = 0; i < n; i++)
			  for (j = 0; j < n; j++)
			    {
			      scenNX = x[i][j][q].get (GRB_DoubleAttr_ScenNX);
			      solx[i][j][q] = scenNX;

			      scenNX = y[i][j][q].get (GRB_DoubleAttr_ScenNX);
			      soly[i][j][q] = scenNX;

			      for (k = 0; k < n; k++)
				{
				  scenNX = u[i][j][k][q].get (
				      GRB_DoubleAttr_ScenNX);
				  solu[i][j][k][q] = scenNX;
				}
			    }

		      if (PRINT_VAR_VALUR)
			printVar (solx, soly, solu, origin);

		      for (i = 0; i < n; i++)
			{
			  for (j = 0; j < n; j++)
			    {
			      delete[] solx[i][j];
			      delete[] soly[i][j];
			      for (k = 0; k < n; k++)
				delete[] solu[i][j][k];
			      delete[] solu[i][j];
			    }
			  delete[] solx[i];
			  delete[] soly[i];
			  delete[] solu[i];
			}
		    }
		}

	      if (USE_LR_MULTIPLIER_TYPE_B)
		updateLRmultiplierTypeB (soly_d, LR_u, LR_rou, countItr,
					 numArcs);

	      double LR_miu;
	      if (USE_LR_MULTIPLIER_TYPE_C)
		updateLRmultiplierTypeC (soly_d, LR_lamda, &LR_miu, LR_u, LB,
					 profitInLRdual);
	      if (LR_miu <= LR_miu_tolerance)
		{
		  cout
		      << "=== will stop loop because LR_miu < LR_miu_tolerance "
		      << LR_miu_tolerance << endl;
		  break;
		}

	      printf ("===> LB UB gap = %lf \n", (UB - LB) / UB);
	      if ((UB - LB) / UB <= LR_gap_tolerance)
		cout << "=== will stop loop because gap < LR_gap_toleranc "
		    << LR_gap_tolerance << endl;

	      reportTime (endTimeOfLastIteration,
			  endTimeOfLastIterationWallClock);
	      endTimeOfLastIteration = clock ();
	      endTimeOfLastIterationWallClock = high_resolution_clock::now ();

	    } //end of while loop

	  for (i = 0; i < n; i++)
	    {
	      for (j = 0; j < n; j++)
		{
		  delete[] solx_d[i][j];
		  delete[] soly_d[i][j];
		  for (k = 0; k < n; k++)
		    delete[] solu_d[i][j][k];
		  delete[] solu_d[i][j];
		}
	      delete[] solx_d[i];
	      delete[] soly_d[i];
	      delete[] solu_d[i];
	    }

	  cout << endl
	      << "===> The total time since setting up Gurobi environment: "
	      << endl;
	  reportTime (beginTime, beginWallClock);

	  clock_t end = clock ();
	  double second = (double) (end - beginTime) / CLOCKS_PER_SEC;

	  auto endWallClock = high_resolution_clock::now ();
	  auto elapsedWallClock = duration_cast < std::chrono::nanoseconds
	      > (endWallClock - beginWallClock);

	  cout << "===> SUMMARY:" << endl;
	  printf ("SolutionTimeCPU = %lf\n", second);
	  printf ("SolutionTimeWallClock = %.3f seconds\n",
		  elapsedWallClock.count () * 1e-9);
	  printf ("SolutionLB = %lf\n", LB);
	  printf ("SolutionUB = %lf\n", UB);
	  printf ("SolutionGap = %lf\n", (UB - LB) / UB);

	} //end of if using LR condition
    } //end of try opening Gurobi environment
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
  // }
  //}

  return 0;
}
