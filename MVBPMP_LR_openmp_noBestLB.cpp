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

#include <omp.h>

using namespace std;
using namespace std::chrono;
using namespace std::string_literals;
//to generate vehicle data file name

int numV = 3; //number of vehicles; will be overwritten when reading vehicle file
//number of nodes in instance
//in main(), we also use n to present number of nodes
//since it's easy to read variables with short subscript, like x[n][n][numV]
int NUM_NODES;

//when runing on HPC, please change the required cores in sbatch file accordingly
//for example, if requiredThreads = 9, and 3 vehicles in the instance, then the cores
//should be set up at least 9*3=27
int NUM_THREADS_VEH = 8; //number of threads for solving each vehicle's BPMP model

double INITIAL_U_COEFFICIENT;

bool PRINT_x_WHEN_INTEGER_SOL = false;
bool PRINT_y_WHEN_INTEGER_SOL = false;

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
	//printf ("slack=%lf %d %d\n", slack[i][j], i, j);
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
  //********* only read graph info and vehicle info by arguments
  //********* doesn't go through the graph info in the data folder
  if (argc != 4)
    {
      cout
	  << "Usage: ./mvbpmp_openmp.x nodesDataNameAndPath numberOfVehicles vehicleDataNameAndPath"
	  << endl;
      return 1;
    }

  if (argc == 4)
    {
      cout << "Nodes Data: " << argv[1] << endl;
      cout << "Number of Vehicles: " << argv[2] << endl;
      cout << "Vehicles Data: " << argv[3] << endl;
      numV = stoi (argv[2]);
    }

  clock_t beginTime, endTimeOfLastIteration;

  auto beginWallClock = high_resolution_clock::now ();
  auto endTimeOfLastIterationWallClock = high_resolution_clock::now ();

  readData route;

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

  //============== start Lagrangian Relaxation ==============

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
  for (int i = 0; i < n; i++)
    LR_u[i] = new double[n];

  //initialize LR_u
  double numArcs = n * (n - 1);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
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
  double profit_d[numV];
  for (int i = 0; i < n; i++)
    {
      solx_d[i] = new double*[n];
      soly_d[i] = new double*[n];
      solu_d[i] = new double**[n];
      for (int j = 0; j < n; j++)
	{
	  solx_d[i][j] = new double[numV];
	  soly_d[i][j] = new double[numV];
	  solu_d[i][j] = new double*[n];
	  for (int k = 0; k < n; k++)
	    solu_d[i][j][k] = new double[numV];
	}
    }

  int maxItr = 3;
  int countItr = 0;

  cout << "===> LR_gap_tolerance = " << LR_gap_tolerance << endl;

  beginTime = clock ();

  endTimeOfLastIterationWallClock = high_resolution_clock::now ();
  while ((UB - LB) / UB > LR_gap_tolerance)
    {

      countItr++;
      //stop after a certain iterations
      //if (countItr > maxItr)break;

      cout << "============= Itr " << countItr << " ==============" << endl;
      cout << "LB = " << LB << endl;
      cout << "UB = " << UB << endl;

//=======================> SOLVING BPMP FOR EACH VEHICLE IN PARALLEL HERE <======================

      // figure out how many can run at once based on available cores
      const int max_threads = omp_get_max_threads ();
      cout << "max_threads: " << max_threads << endl;
      const int numParallelVehicles = max_threads / NUM_THREADS_VEH;
      cout << "numParallelVehicles: " << numParallelVehicles << endl;

      if (numParallelVehicles < numV)
	{
	  cout << "Quit running: there are not enough threads available."
	      << endl;
	  exit (1);
	}

      // loop over vehicles
//#pragma omp parallel for num_threads(numParallelVehicles) default(shared)
#pragma omp parallel for num_threads(numV) default(shared)
      for (int vehicleIndex = 0; vehicleIndex < numV; vehicleIndex++)
	{

	  int nthreads = omp_get_num_threads ();
	  cout << "NUM_THREADS=" << nthreads << endl;

	  // start calling Gurobi
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

	      //set up logfile
	      model.set (
		  "LogFile",
		  "mvbpmp_" + itos (n) + "_" + itos (vehicleIndex) + ".log");
	      //send the log to a file only
	      model.set (GRB_IntParam_LogToConsole, 0);
	      //keep outputflag 1 to print out log
	      model.set (GRB_IntParam_OutputFlag, 1);

	      // set the number of threads to required number of Threads
	      model.set (GRB_IntParam_Threads, NUM_THREADS_VEH);

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

	      // force some x, y variables to be zeros
	      for (i = 0; i < n; i++)
		{
		  x[i][i].set (GRB_DoubleAttr_UB, 0);
		  y[i][i].set (GRB_DoubleAttr_UB, 0);
		  x[i][vehicleOrigin].set (GRB_DoubleAttr_UB, 0);
		  y[i][vehicleOrigin].set (GRB_DoubleAttr_UB, 0);
		  x[n - 1][i].set (GRB_DoubleAttr_UB, 0);
		  y[n - 1][i].set (GRB_DoubleAttr_UB, 0);
		}
	      // force some triples variables to be zeros
	      for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		  {
		    u[n - 1][i][j].set (GRB_DoubleAttr_UB, 0);
		    u[i][vehicleOrigin][j].set (GRB_DoubleAttr_UB, 0);
		    u[i][j][vehicleOrigin].set (GRB_DoubleAttr_UB, 0);
		    u[i][j][n - 1].set (GRB_DoubleAttr_UB, 0);
		  }

	      //==============generate constraints in Gurobi================
	      // vehicle goes out of vehicle's origin
	      GRBLinExpr expr1 = 0.0;
	      for (i = 0; i < n; i++)
		expr1 += x[vehicleOrigin][i];
	      model.addConstr (expr1 == 1, "origin");

	      // vehicle goes back to node n
	      GRBLinExpr expr2 = 0.0;
	      for (i = 0; i < n; i++)
		expr2 += x[i][n - 1];
	      model.addConstr (expr2 == 1, "destination");

	      //flow conservation (the last node is the destination)
	      for (k = 0; k < n - 1; k++)
		{
		  if (k != vehicleOrigin)
		    {
		      GRBLinExpr expr = 0;
		      for (i = 0; i < n; i++)
			expr += x[i][k];
		      for (j = 0; j < n; j++)
			expr -= x[k][j];
		      model.addConstr (expr == 0,
				       "flow_conservation_" + itos (k));
		    }
		}

	      // distance
	      GRBLinExpr expr3 = 0.0;
	      for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		  expr3 += dis[i][j] * x[i][j];
	      model.addConstr (expr3 <= route.DIS, "distance");

	      // node degree less than 1
	      for (k = 0; k < n; k++)
		{
		  GRBLinExpr expr = 0.0;
		  for (i = 0; i < n - 1; i++)
		    expr += x[i][k];
		  model.addConstr (expr <= 1, "indegree_" + itos (k));
		}

	      // subtour elimination
	      for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		  {
		    GRBLinExpr expr = 0.0;
		    expr += s[i] - s[j] + (n - 1) * x[i][j] + (n - 3) * x[j][i];
		    model.addConstr (expr <= n - 2,
				     "s_" + itos (i) + "_" + itos (j));
		  }

	      // arc flow
	      for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		  {
		    GRBLinExpr expr = 0.0;
		    expr += wt[i][j] * y[i][j] - theta[i][j];
		    for (k = 0; k < n; k++)
		      expr += u[i][k][j] + u[k][j][i] - u[i][j][k];
		    model.addConstr (expr == 0,
				     "flow_" + itos (i) + "_" + itos (j));
		  }

	      // arc flow upperbound
	      for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		  {
		    GRBLinExpr expr = 0.0;
		    expr += theta[i][j] - Q * x[i][j];
		    model.addConstr (expr <= 0,
				     "flowBound_" + itos (i) + "_" + itos (j));
		  }

	      //==============set objective function in Gurobi================

	      GRBLinExpr obj = 0.0;
	      for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		  {
		    obj += (price * dis[i][j] * wt[i][j] - LR_u[i][j])
			* y[i][j];
		    obj -= cost * dis[i][j] * theta[i][j];
		    obj -= cost * vw * dis[i][j] * x[i][j];
		  }
	      model.setObjective (obj, GRB_MAXIMIZE);

	      //==============solve the model in Gurobi================
	      model.optimize ();

	      // write model to file
	      model.write (
		  "BPMP_" + itos (n) + "_" + itos (vehicleIndex) + ".lp");

	      //==============Extract solution from Gurobi================

	      if (model.get (GRB_IntAttr_SolCount) > 0)
		{
		  profit_d[vehicleIndex] = model.get (GRB_DoubleAttr_ObjVal);

		  double **sol = new double*[n];

		  double ***solU = new double**[n]; //gurobi solution for triples var u
		  for (i = 0; i < n; i++)
		    solU[i] = new double*[n];

		  printf ("Selected arcs: \n");
		  for (i = 0; i < n; i++)
		    {
		      sol[i] = model.get (GRB_DoubleAttr_X, x[i], n);
		      for (j = 0; j < n; j++)
			{
			  if (sol[i][j] > 0.9)
			    printf ("%d -- %d\n", i + 1, j + 1);

			  solx_d[i][j][vehicleIndex] = sol[i][j];
			}
		    }

		  printf ("Selected requests: \n");
		  for (i = 0; i < n; i++)
		    {
		      sol[i] = model.get (GRB_DoubleAttr_X, y[i], n);
		      for (j = 0; j < n; j++)
			{
			  if (sol[i][j] > 0.9)
			    printf ("%d -- %d\n", i + 1, j + 1);

			  soly_d[i][j][vehicleIndex] = sol[i][j];
			}

		    }

		  for (i = 0; i < n; i++)
		    for (j = 0; j < n; j++)
		      {
			solU[i][j] = model.get (GRB_DoubleAttr_X, u[i][j], n);
			for (k = 0; k < n; k++)
			  solu_d[i][j][k][vehicleIndex] = solU[i][j][k];
		      }

		  for (i = 0; i < n; i++)
		    {
		      delete[] sol[i];
		      for (j = 0; j < n; j++)
			delete[] solU[i][j];
		      delete[] solU[i];
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
	      cout << "Error during optimization" << endl;
	    }

	  for (i = 0; i < n; i++)
	    {
	      delete[] x[i];
	      delete[] y[i];
	    }
	  delete env;
	}	      //end of openmp parallel running using Gurobi

      //collect profit from all scenarios
      //each scenario represent one vehicle's schedule
      double totalProfit = 0;
      for (int i = 0; i < numV; i++)
	{
	  totalProfit += profit_d[i];
	  cout << "profit for vehicle " << i << ": " << profit_d[i] << endl;
	}

      cout << "===> sum of OBJ VALUE: " << totalProfit << endl;

      //add the lagrangian multiplier part sum_LR_u back to obj
      //since this is constant and not considered in scenarios
      //sum_LR_u_times_y is for printing here and future use when solution
      //is feasible for the original problem
      double sum_LR_u = 0;
      double sum_LR_u_times_y = 0;
      bool findConflictPickup = false;
      for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
	  {
	    sum_LR_u += LR_u[i][j];

	    double sum = 0;
	    for (int q = 0; q < numV; q++)
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
      cout << "totalProfit in obj including LR = " << totalProfit << endl;

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
	  cout << "===> LR_lamda is less than LR_min_damda " << LR_min_lamda
	      << endl;
	  cout << "end loop";
	  break;
	}

      //===> if there is no request conflict
      //===> use the current solution without LR item as LB
      if (!findConflictPickup)
	{
	  cout << "===> There is no request conflict in solution :)" << endl;

	  //check if the LR_item is close to zero
	  //if yes, then this is the optimal solution
	  double LR_item = sum_LR_u - sum_LR_u_times_y;
	  cout << "===> the LR items added in obj = " << LR_item << endl;

	  if ((LR_item >= 0 && LR_item <= LR_complementarity_tolerance)
	      || (LR_item < 0 && LR_item >= -LR_complementarity_tolerance))
	    {
	      LB = profitInLRdual - LR_item;
	      cout << "===> the abs(LR item) added in obj is less than "
		  << LR_complementarity_tolerance << endl;
	      cout << "===> the optimal solution is found :)" << endl;
	      cout << "===> end loop." << endl;

	      reportTime (endTimeOfLastIteration,
			  endTimeOfLastIterationWallClock);

	      printVar (solx_d, soly_d, solu_d, origin);

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
	    updateLRmultiplierTypeB (soly_d, LR_u, LR_rou, countItr, numArcs);

	  double LR_miu;
	  if (USE_LR_MULTIPLIER_TYPE_C)
	    updateLRmultiplierTypeC (soly_d, LR_lamda, &LR_miu, LR_u, LB,
				     profitInLRdual);

	  reportTime (endTimeOfLastIteration, endTimeOfLastIterationWallClock);
	  endTimeOfLastIteration = clock ();
	  endTimeOfLastIterationWallClock = high_resolution_clock::now ();

	  if ((UB - LB) / UB <= LR_gap_tolerance)
	    cout << "=== will stop loop because gap " << (UB - LB) / UB
		<< " < LR_gap_toleranc " << LR_gap_tolerance << endl;

	  if (LR_miu <= LR_miu_tolerance)
	    {
	      cout << "=== will stop loop because LR_miu < LR_miu_tolerance "
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
		  if (soly_d[i][j][q] < 0.5) //if the cargo is not selected in LR dual, then do not consider them in MVBPMP
		    {
		      y[i][j][q].set (GRB_DoubleAttr_LB, 0.0);
		      y[i][j][q].set (GRB_DoubleAttr_UB, 0.0);
		    }
		}

	  // set up constraints
	  GRBConstr *vehOriginConstr = 0;
	  GRBConstr *vehDestConstr = 0;
	  vehOriginConstr = new GRBConstr[numV];
	  vehDestConstr = new GRBConstr[numV];

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

		    /*
		     expr += wt[i][j] * y[i][j][q] - theta[i][j][q];
		     for (k = 0; k < n; k++)
		     expr += u[i][k][j][q] + u[k][j][i][q] - u[i][j][k][q];
		     */
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
	  //if (!ADD_CUTS){}
	  for (i = 0; i < n - 1; i++)
	    for (j = 0; j < n; j++)
	      {
		GRBLinExpr expr = 0.0;
		for (q = 0; q < numV; q++)
		  expr += y[i][j][q];
		model.addConstr (expr <= 1,
				 "one-one" + itos (i) + "_" + itos (j));
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

	  // Optimize model
	  model.optimize ();

	  //write model to file
	  //model.write ("MVBPMP.lp");

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
	      cout << "Optimization was stopped with status " << status << endl;
	      return 1;
	    }

	  //=====================> READ SOLUTION FROM GUROBI <============================
	  // Extract solution
	  if (model.get (GRB_IntAttr_SolCount) > 0)
	    {

	      double objtemp = model.get (GRB_DoubleAttr_ObjVal);
	      cout << "OBJ: " << objtemp << endl;
	      //update LB
	      if (objtemp > LB)
		LB = objtemp;
	      cout << "===> LB = " << LB << endl;

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
		      solx[i][j] = model.get (GRB_DoubleAttr_X, x[i][j], numV);
		      soly[i][j] = model.get (GRB_DoubleAttr_X, y[i][j], numV);
		    }
		}

	      printf ("Selected arcs: \n");
	      for (q = 0; q < numV; q++)
		for (i = 0; i < n; i++)
		  for (j = 0; j < n; j++)
		    {
		      if (solx[i][j][q] > 0.9)
			printf ("%d -- %d (%d)\n", i + 1, j + 1, q + 1);
		    }

	      printf ("Selected requests: \n");
	      for (q = 0; q < numV; q++)
		for (i = 0; i < n; i++)
		  for (j = 0; j < n; j++)
		    {
		      if (soly[i][j][q] > 0.9)
			printf ("%d -- %d (%d)\n", i + 1, j + 1, q + 1);
		    }

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

      if (USE_LR_MULTIPLIER_TYPE_B)
	updateLRmultiplierTypeB (soly_d, LR_u, LR_rou, countItr, numArcs);

      double LR_miu;
      if (USE_LR_MULTIPLIER_TYPE_C)
	updateLRmultiplierTypeC (soly_d, LR_lamda, &LR_miu, LR_u, LB,
				 profitInLRdual);
      if (LR_miu <= LR_miu_tolerance)
	{
	  cout << "=== will stop loop because LR_miu < LR_miu_tolerance "
	      << LR_miu_tolerance << endl;
	  break;
	}

      printf ("===> LB UB gap = %lf \n", (UB - LB) / UB);
      if ((UB - LB) / UB <= LR_gap_tolerance)
	cout << "=== will stop loop because gap < LR_gap_toleranc "
	    << LR_gap_tolerance << endl;

      reportTime (endTimeOfLastIteration, endTimeOfLastIterationWallClock);
      endTimeOfLastIteration = clock ();
      endTimeOfLastIterationWallClock = high_resolution_clock::now ();

    } //end of while loop

  for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
	{
	  delete[] solx_d[i][j];
	  delete[] soly_d[i][j];
	  for (int k = 0; k < n; k++)
	    delete[] solu_d[i][j][k];
	  delete[] solu_d[i][j];
	}
      delete[] solx_d[i];
      delete[] soly_d[i];
      delete[] solu_d[i];
    }

  cout << endl << "===> The total time since setting up Gurobi environment: "
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

  return 0;
}
