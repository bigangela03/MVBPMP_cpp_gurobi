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
//for example, if NUM_THREADS_VEH = 9, and 3 vehicles in the instance, then the cores
//should be set up at least 9*3=27
int NUM_THREADS_VEH = 10; //number of threads for solving each vehicle's BPMP model

double INITIAL_U_COEFFICIENT;

bool USE_UB_FROM_LR = false;
//feed the LR solution to Gurobi and find the optimal solution
bool DO_STAGE_TWO = true;

bool PRINT_x_WHEN_INTEGER_SOL = false;
bool PRINT_y_WHEN_INTEGER_SOL = false;

//according to p174 on "Integer Programming" by Laurence A. Wolsey 1st edition
bool USE_LR_MULTIPLIER_TYPE_C = true;

bool PRINT_VAR_VALUR = false;
bool PRINT_CONFLICT_PICKUP = true;

double bigM = 10000000;
double epsilon = numeric_limits<double>::epsilon ();

string
itos (int i)
{
  stringstream s;
  s << i;
  return s.str ();
}

void
updateLRmultiplierTypeC_flow (double**, double, double*, double**, double,
			      double);
void
printVar (double**, double**, double***, int*);
void
reportTime (clock_t, auto);
void
storeBestLB (double**, double**, double***, double*, double**, double**,
	     double***, double*);

//====================================================================================
void
print2Darray (int n, double **array, bool reverse)
{
  int i, j;
  if (reverse)
    {
      for (i = 0; i < n; i++)
	{
	  for (j = 0; j < n; j++)
	    printf ("%f  ", -array[i][j]);
	  cout << endl;
	}
    }
  else
    {
      for (i = 0; i < n; i++)
	{
	  for (j = 0; j < n; j++)
	    printf ("%f  ", array[i][j]);
	  cout << endl;
	}
    }

}

void
findPath (int startNode, int endNode, int **FW_pred, vector<int> &path)
{
  int pred = FW_pred[startNode][endNode];
  path.insert (path.begin (), pred);
  cout << "pred = " << pred << endl;
  if (pred == 0)
    return;
  else
    findPath (0, pred, FW_pred, path);
}

void
storeBestLB (double **solx_d, double **soly_d, double ***solu_d, double *sols_d,
	     double **solx_best, double **soly_best, double ***solz_best,
	     double *sols_best, int n)
{
  int i, j, k, q;
  for (i = 0; i < n; i++)
    {
      sols_best[i] = sols_d[i];
      for (j = 0; j < n; j++)
	{
	  solx_best[i][j] = solx_d[i][j];
	  soly_best[i][j] = soly_d[i][j];

	  for (k = 0; k < n; k++)
	    solz_best[i][j][k] = solu_d[i][j][k];
	}
    }
}
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
printVar (double **solx, double **soly, double ***solu, int *origin)
{
  int i, j, k, q;
  int n = NUM_NODES;
  //cout << "SOLUTION:" << endl;

  cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
      << endl;
  printf ("Vehicle %d starting at node %d:\n", q + 1, origin[q] + 1);

  printf ("x=1:\n");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      if (solx[i][j] > 0.99 && solx[i][j] < 1.01)
	printf ("(%d,%d) ", i + 1, j + 1);
  printf ("\n");

  printf ("y=1:\n");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      if (soly[i][j] > 0.99 && soly[i][j] < 1.01)
	printf ("(%d,%d) ", i + 1, j + 1);

  printf ("\n");
  cout << "u>0.000001:" << endl;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
	if (solu[i][j][k] > 0.000001)
	  printf ("(%d,%d,%d) ", i + 1, j + 1, k + 1);
  printf ("\n");

  cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
      << endl;
}

void
updateLRmultiplierTypeC_flow (double **flow_slack, double LR_lamda,
			      double *LR_miu, double **LR_u, double LB,
			      double totalProfit)
{
  int i, j, q;
  //update LR multiplier LR_u
  double norm = 0;
  int n = NUM_NODES;
  double slack[n][n];
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      norm += flow_slack[i][j] * flow_slack[i][j];

  cout << "norm = " << norm << endl;

  *LR_miu = LR_lamda * (totalProfit - LB) / norm;
  cout << "updated LR_miu = " << *LR_miu << endl;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	LR_u[i][j] = LR_u[i][j] - (*LR_miu) * flow_slack[i][j];
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
	  << "Usage: ./bpmp_lr.x nodesDataNameAndPath numberOfVehicles vehicleDataNameAndPath"
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
  //route.printData ();

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

  bool findOptimalSolution = false;

  //double LR_u[n][n];	      //the lagrangian multiplier
  double **LR_u = new double*[n];
  for (int i = 0; i < n; i++)
    LR_u[i] = new double[n];

  //initialize LR_u
  double numArcs = n * (n - 1);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      LR_u[i][j] = 0.1 / numArcs;

  //the lagrangian multiplier for distance constraint
  double LR_dis_u = 0.1 / numArcs;

  double LR_lowestUpperBound = bigM;

  //store the solution of MVBPMP_LR in each iteration
  double **solx_d = new double*[n];	    //solution x for LR dual
  double **soly_d = new double*[n];
  double ***solz_d = new double**[n];
  double soltheta_d[n][n];
  //double sols_d[n][numV];
  double *sols_d = new double[n];
  double profit_d[numV];

  double **solx_best = new double*[n];
  double **soly_best = new double*[n];
  double ***solz_best = new double**[n];
  double *sols_best = new double[n];

  double **coeff_x_d = new double*[n]; //coefficient of x variables in LR dual
  double **coeff_y_d = new double*[n];
  double ***coeff_z_d = new double**[n];

  for (int i = 0; i < n; i++)
    {
      solx_d[i] = new double[n];
      soly_d[i] = new double[n];
      solz_d[i] = new double*[n];

      solx_best[i] = new double[n];
      soly_best[i] = new double[n];
      solz_best[i] = new double*[n];

      coeff_x_d[i] = new double[n];
      coeff_y_d[i] = new double[n];
      coeff_z_d[i] = new double*[n];

      for (int j = 0; j < n; j++)
	{
	  solz_d[i][j] = new double[n];
	  solz_best[i][j] = new double[n];
	  coeff_z_d[i][j] = new double[n];
	}
    }

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
	solx_d[i][j] = 0;
	soly_d[i][j] = 0;
	for (int k = 0; k < n; k++)
	  solz_d[i][j][k] = 0;
      }

  int currentVehOrigin = 0;

  int maxItr = 2;
  int countItr = 0;

  cout << "===> LR_gap_tolerance = " << LR_gap_tolerance << endl;

  beginTime = clock ();

  endTimeOfLastIterationWallClock = high_resolution_clock::now ();

  while ((UB - LB) / UB > LR_gap_tolerance)
    {

      countItr++;
      //stop after a certain iterations
      if (countItr > maxItr)
	break;

      cout << "============= Itr " << countItr << " ==============" << endl;
      cout << "LB = " << LB << endl;
      cout << "UB = " << UB << endl;

      int i, j, k;

      cout << "===> calculating vars' coefficients of LR dual " << endl;
      for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	  {

	    coeff_x_d[i][j] = -((cost * vw + LR_dis_u) * dis[i][j]
		- Q * LR_u[i][j]);

	    coeff_y_d[i][j] = ((price - cost) * dis[i][j] - LR_u[i][j])
		* wt[i][j];

	    for (k = 0; k < n; k++)
	      {
		coeff_z_d[i][j][k] = -(cost
		    * (dis[i][k] + dis[k][j] - dis[i][j]) + LR_u[i][k]
		    + LR_u[k][j] - LR_u[i][j]);
		if (coeff_z_d[i][j][k] > epsilon)
		  printf ("coeff_z_d[%d][%d][%d]=%.18lf\n", i, j, k,
			  coeff_z_d[i][j][k]);
	      }
	  }

      //some x, y and z variables are always zero in the model
      //so we have to force their coeffieients 0 in Floy-Warshall algorithm
      // force some x, y variables to be zeros
      for (i = 0; i < n; i++)
	{
	  coeff_x_d[i][i] = 0.0;
	  coeff_y_d[i][i] = 0.0;
	  //coeff_x_d[i][currentVehOrigin] = 0.0;
	  //coeff_y_d[i][currentVehOrigin] = 0.0;
	  coeff_x_d[n - 1][i] = 0.0;
	  coeff_y_d[n - 1][i] = 0.0;
	}
      // force some triples variables to be zeros
      for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	  {
	    coeff_z_d[n - 1][i][j] = 0.0;
	    coeff_z_d[i][currentVehOrigin][j] = 0.0;
	    coeff_z_d[i][j][currentVehOrigin] = 0.0;
	    coeff_z_d[i][j][n - 1] = 0.0;
	  }

      //=======================> SOLVING BPMP USING Floyd-Warshall algorithm <======================
      //use -coeff_x_d as distance in Floyd-Warshall algorithm for shortest distance
      //in case that some -coeff_x_d is negative, we use Floyd-Warshall to handle negative distance

      cout << "===> start Floyd-Warshall algorithm " << endl;
      //double FW_cost[n][n];
      //double FW_pred[n][n];
      double **FW_cost = new double*[n];
      int **FW_pred = new int*[n];
      for (i = 0; i < n; i++)
	{
	  FW_cost[i] = new double[n];
	  FW_pred[i] = new int[n];
	}

      for (i = 0; i < n; i++)
	{
	  FW_cost[i][i] = 0.0;
	  for (j = 0; j < n; j++)
	    {
	      FW_cost[i][j] = -coeff_x_d[i][j];
	      FW_pred[i][j] = i;
	    }
	}

      //vehicle can't go back to (n-1) or leave from currentVehOrigin
      for (i = 0; i < n; i++)
	{
	  FW_cost[i][currentVehOrigin] = bigM;
	  FW_cost[n - 1][i] = bigM;
	}

      for (k = 0; k < n; k++)
	for (i = 0; i < n; i++)
	  for (j = 0; j < n; j++)
	    {
	      //printf ("checking k=%d i=%d j=%d\n", k, i, j);
	      if (i != j && i != k && k != j)
		if (FW_cost[i][j] > FW_cost[i][k] + FW_cost[k][j])
		  {
		    FW_cost[i][j] = FW_cost[i][k] + FW_cost[k][j];
		    FW_pred[i][j] = FW_pred[k][j];
		    //cout << "find lower cost" << endl;
		    //printf ("FW_cost[%d][%d]=%lf\n", i, j, FW_cost[i][j]);
		  }
	    }

      bool reverse = true;      //print -array(i,j) when true
      cout << "cost array" << endl;
      print2Darray (n, coeff_x_d, reverse);

      cout << "min cost" << endl;
      print2Darray (n, FW_cost, false);

      cout << "===> The min cost: " << FW_cost[currentVehOrigin][n - 1] << endl;
      cout << "===> The max obj after Floyd-Warshall algorithm: "
	  << -FW_cost[currentVehOrigin][n - 1] << endl;

      for (i = 0; i < n; i++)
	{
	  if (FW_cost[i][i] < 0)
	    {
	      printf ("c(%d,%d)=%lf <0  ", i, j, FW_cost[i][j]);
	      cout << "There is negative cycle. Quit running." << endl;
	      return 0;
	    }
	}

      if (countItr >= maxItr)
	exit (1);

      //======find the x var value based on F-W algorithm result======
      //the nodes on path will be stored in path array

      cout << "===> find path based on Floyd-Warshall result" << endl;

      //int path[n];
      bool findConflict = false;
      vector<int> path;

      findPath (0, n - 1, FW_pred, path);
      path.push_back (n - 1);

      for (i = 0; i < path.size (); i++)
	cout << path[i] << " ";
      cout << endl;

      double totalDistnce = 0;
      for (i = 0; i < path.size () - 1; i++)
	{
	  int str = path[i];
	  int end = path[i + 1];
	  totalDistnce += dis[str][end];
	  solx_d[str][end] = 1;
	}

      cout << "Total distance = " << totalDistnce << endl;

      if (totalDistnce > route.DIS)
	findConflict = true;

      //======calcute the total profit for LR dual
      double totalProfit = -FW_cost[currentVehOrigin][n - 1];

      //======add the lagrangian multiplier part for y and z back to obj
      double sum_LR_y_z = 0;

      for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	  {
	    if (coeff_y_d[i][j] > 0)
	      {
		sum_LR_y_z += coeff_y_d[i][j];
		//======find y var value based on the y's LR dual coefficient
		soly_d[i][j] = 1;
		//printf ("coeff_y_d(%d,%d)=%lf\n", i, j, coeff_y_d[i][j]);
	      }
	    //coefficients for z_(i,j,k) in LR dual is <=Q
	    for (k = 0; k < n; k++)
	      if (coeff_z_d[i][j][k] > 0)
		{
		  sum_LR_y_z += Q * coeff_z_d[i][j][k];
		  solz_d[i][j][k] = Q;
		}
	  }

      sum_LR_y_z += LR_dis_u * route.DIS;

      cout << "sum_LR_y_z = " << sum_LR_y_z << endl;

      totalProfit += sum_LR_y_z;
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

      //======check the capacity constraint
      double **flow_slack = new double*[n];
      for (i = 0; i < n; i++)
	flow_slack[i] = new double[n];

      double LR_item = 0;
      for (i = 0; i < n - 1; i++)
	for (j = 1; j < n; j++)
	  {
	    double flow = 0.0;
	    flow += wt[i][j] * soly_d[i][j];
	    for (k = 0; k < n; k++)
	      flow += solz_d[i][k][j] + solz_d[k][j][i] - solz_d[i][j][k];

	    LR_item += LR_u[i][j] * (Q * solx_d[i][j] - flow);

	    flow_slack[i][j] = Q * solx_d[i][j] - flow;

	    if (flow > Q * solx_d[i][j])
	      findConflict = true;
	  }

      LR_item += LR_dis_u * (route.DIS - totalDistnce);

      //===> if there is no request conflict
      //===> use the current solution without LR item as LB
      if (!findConflict)
	{
	  cout << "===> There is no request conflict in solution :)" << endl;

	  //check if the LR_item is close to zero
	  //if yes, then this is the optimal solution

	  cout << "===> the LR items added in obj = " << LR_item << endl;

	  if ((LR_item >= 0 && LR_item <= LR_complementarity_tolerance)
	      || (LR_item < 0 && LR_item >= -LR_complementarity_tolerance))
	    {
	      LB = profitInLRdual - LR_item;

	      findOptimalSolution = true;

	      cout << "===> the abs(LR item) added in obj is less than "
		  << LR_complementarity_tolerance << endl;
	      cout << "===> the optimal solution is found :)" << endl;
	      cout << "===> end loop." << endl;

	      reportTime (endTimeOfLastIteration,
			  endTimeOfLastIterationWallClock);

	      cout << "**************** THE OPTIMAL SOLUTION ****************"
		  << endl;
	      printVar (solx_d, soly_d, solz_d, origin);

	      //no need to store the solution since the optimal solution is found
	      //and there is no need to feed solution to MVBPMP model to find optimal solution
	      break;
	    }

	  //if LR_item is still big, then we need to keep iteration
	  cout << "===> Use this solution to calculate LB." << endl;

	  //since no request conflict
	  //so the solution for x, y, u are feasible for the original model
	  //but - u_(i,j)*( 1-sum_q y(i,j,q) ) was added to obj
	  //so we need to recalculte total profit
	  double profitForLB = profitInLRdual;
	  profitForLB -= LR_item;

	  if (profitForLB > LB)
	    {
	      LB = profitForLB;
	      cout << "===> find a better LB = " << LB << endl;

	      //store the solution as the best LB
	      storeBestLB (solx_d, soly_d, solz_d, sols_d, solx_best, soly_best,
			   solz_best, sols_best, n);
	    }

	  double LR_miu;
	  if (USE_LR_MULTIPLIER_TYPE_C)
	    updateLRmultiplierTypeC_flow (flow_slack, LR_lamda, &LR_miu, LR_u,
					  LB, profitInLRdual);

	  reportTime (endTimeOfLastIteration, endTimeOfLastIterationWallClock);
	  endTimeOfLastIteration = clock ();
	  endTimeOfLastIterationWallClock = high_resolution_clock::now ();

	  if (UB - LB < 0)
	    {
	      cout << "UB = " << UB << endl;
	      cout << "LB = " << LB << endl;
	      cout
		  << "UB is less than LB, which is an error. Stop running and check code."
		  << endl;
	      exit (1);
	    }

	  if ((UB - LB) / UB <= LR_gap_tolerance)
	    {
	      cout << "=== will stop loop because gap " << (UB - LB) / UB
		  << " < LR_gap_toleranc " << LR_gap_tolerance << endl;

	      //here if(profitForLB > LB), it is already checked before
	      //and the new LB is stored
	      //if(profitForLB <= LB), then no need to store the sol?_d solution as best LB
	      //so to conclude, no need to storeBestLB here
	      //storeBestLB (solx_d, soly_d, solz_d, sols_d, solx_best, soly_best, solz_best, sols_best);
	    }
	  //this condition is added after adding OpenMP feature.
	  //if no added, the code will calculate MVBPMP for new LB.
	  //better LB might be found, but not necessary since gap is already within gap tolerance
	  else
	    continue;

	  if (LR_miu <= LR_miu_tolerance)
	    {
	      cout << "=== will stop loop because LR_miu < LR_miu_tolerance "
		  << LR_miu_tolerance << endl;

	      //here if(profitForLB > LB), it is already checked before
	      //and the new LB is stored
	      //if(profitForLB <= LB), then no need to store the sol?_d solution as best LB
	      //so to conclude, no need to storeBestLB here
	      //storeBestLB (solx_d, soly_d, solz_d, sols_d, solx_best, soly_best, solz_best, sols_best);

	      break;
	    }
	  else
	    continue; //jump to the next iteration of while loop
	}

      LB = 1;

      double LR_miu;
      if (USE_LR_MULTIPLIER_TYPE_C)
	updateLRmultiplierTypeC_flow (flow_slack, LR_lamda, &LR_miu, LR_u, LB,
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
	delete[] solz_d[i][j];

      delete[] solx_d[i];
      delete[] soly_d[i];
      delete[] solz_d[i];
    }
  delete[] solx_d;
  delete[] soly_d;
  delete[] solz_d;
  delete[] sols_d;

  if (!findOptimalSolution)
    {
      cout << "**************** THE BEST LB SOLUTION ****************" << endl;
      printVar (solx_best, soly_best, solz_best, origin);
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
