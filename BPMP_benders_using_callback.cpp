//Angela; Feb3, 2024; this code using callback still doesn't work.
//it finishes early with a very high objective value.
//Feb15, 2024; When testing on t10_01_data, it can show that
//adding the lazy constraint when where=MIPNODE, using the best solution S0
//stored in hashmap, doesn't cut S0 node from the tree, it is the best solution found
//Hence I have to gave up on benders decompositon using callback


#include "gurobi_c++.h"
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <string.h>

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
#include <unordered_map>

using namespace std;
using namespace std::chrono;
using namespace std::string_literals;

double price;
double cost_flow;
double cost_vehicle;
double vw;
double Q;
double **wt;
double **dis;

//bool initialSolution = true;
bool PRINT_VAR_VALUE = true;
bool PRINT_VAR_OF_OBJBST_HASHMAP = true;
bool PRINT_VAR_IN_DUAL_SUBPROBLEM = true;
bool PRINT_NODE_INFO = true;
bool PRINT_UB_LB = true;
bool PRINT_MIPSOL_MIPNODE_TITLE = true;
bool PRINT_UNBOUNDED_VECTOR = true;
double LB;
double UB;
bool MIP_SOL_FOUND_BUT_BIG_GAP = true;

GRBEnv *env; //set gurobi environment as a global var
//and define the environment only once in main(),
//then use env.start() to start an empty env.
//when using env = new GRBEnv (), system shows the following message
//Set parameter TokenServer to value "sengr7lic2.smu.edu"

//double BESTBOUND = 1000;
//BIG_M is also used as parameter in hashmap for the same obj
//like, obj, BIG_M+obj, 2*BIG_M + obj, and so on
//in this way we can store the same obj with different x and y var value in hashmap
#define BIG_M 100000
#define EPSILON 0.00001
#define UB_LB_TOLERANCE 0.01
#define OBJ_GAP_TOLERANCE 0.5

struct integerSol
{
  double **x;
  double **y;
};

//hashmap obj -> x, y
unordered_map<double, integerSol> hashmap;

string
itos (int i)
{
  stringstream s;
  s << i;
  return s.str ();
}

bool
subProblemDual (double**, double**, int, double**, double**);
void
deleteVar (double**, int);
void
insertToHashmap (double**, double**, double, unordered_map<double, integerSol>*,
		 int);
void
updateLB (double**, double**, double**, int);
double
mod (double);

class printIntSol : public GRBCallback
{
public:
  GRBVar **xVar;
  GRBVar **yVar;
  //GRBVar *zVar;
  GRBVar zVar;
  int n;
  //printIntSol (GRBVar **xvars, GRBVar **yvars, GRBVar *zvar, int nvar)
  printIntSol (GRBVar **xvars, GRBVar **yvars, GRBVar zvar, int nvar)
  {
    xVar = xvars;
    yVar = yvars;
    zVar = zvar;
    n = nvar;
  }
protected:
  void
  callback ()
  {
    try
      {
	//https://www.gurobi.com/documentation/9.1/refman/cb_codes.html
	//check the above website for the detailed info on callback code
	int i, j;

	double **xVal = new double*[n];
	double **yVal = new double*[n];
	double **piVal = new double*[n];
	double **unbdPiVal = new double*[n];
	for (i = 0; i < n; i++)
	  {
	    xVal[i] = new double[n];
	    yVal[i] = new double[n];
	    piVal[i] = new double[n];
	    unbdPiVal[i] = new double[n];
	  }

	//printf ("---------------------------------------------------------\n");

	if (where == GRB_CB_MIPSOL)
	  {
	    // Found an integer feasible solution
	    if (PRINT_MIPSOL_MIPNODE_TITLE)
	      printf ("==========GRB_CB_MIPSOL=========\n");

	    int nodecnt = (int) getDoubleInfo (GRB_CB_MIPSOL_NODCNT);
	    int solcnt = getIntInfo (GRB_CB_MIPSOL_SOLCNT);
	    double objbst = getDoubleInfo (GRB_CB_MIPSOL_OBJBST);
	    double objbnd = getDoubleInfo (GRB_CB_MIPSOL_OBJBND);
	    double obj = getDoubleInfo (GRB_CB_MIPSOL_OBJ);

	    //zVal = getSolution (*zVar);
	    double zVal = getSolution (zVar);
	    if (PRINT_NODE_INFO)
	      {
		printf ("nodecnt=%d\n", nodecnt);
		printf ("solcnt=%d\n", solcnt);
		printf ("objbst=%lf\n", objbst);
		printf ("objbnd=%lf\n", objbnd);
		printf ("obj=%lf\n", obj);
		printf ("zVal=%lf\n", zVal);
	      }

	    //read solution

	    for (i = 0; i < n; i++)
	      xVal[i] = getSolution (xVar[i], n);
	    for (i = 0; i < n; i++)
	      yVal[i] = getSolution (yVar[i], n);

	    if (PRINT_VAR_VALUE)
	      {
		printf ("x vars: ");
		for (i = 0; i < n; i++)
		  for (j = 0; j < n; j++)
		    if (xVal[i][j] > EPSILON)
		      printf ("%d->%d %lf  ", i, j, xVal[i][j]);
		printf ("\ny vars: ");
		for (i = 0; i < n; i++)
		  for (j = 0; j < n; j++)
		    if (yVal[i][j] > EPSILON)
		      printf ("%d->%d %lf  ", i, j, yVal[i][j]);
		printf ("\n");
	      }

	    double gap = (objbnd - obj) / objbnd;
	    printf ("gap = %lf\n", gap);
	    //if the gap is big, store solution and keep optimizing
	    if (gap > OBJ_GAP_TOLERANCE)
	      {
		//store the MIP solution
		insertToHashmap (xVal, yVal, obj, &hashmap, n);
		cout << "hashmap size after finishing inserting: "
		    << hashmap.size () << endl;
		MIP_SOL_FOUND_BUT_BIG_GAP = true;
		return;
	      }

	    //if the gap is small, update UB and call dual of SP
	    UB = obj;

	    if (PRINT_UB_LB)
	      {
		printf ("======> UB %lf\n", UB);
		printf ("======> LB %lf\n", LB);
	      }
	    if (UB - LB < UB_LB_TOLERANCE)
	      {
		printf ("======> UB(%lf) - LB(%lf) < %lf\n", UB, LB,
		UB_LB_TOLERANCE);
		printf ("terminate optimizing\n");
		abort ();
	      }

	    //======> call sub-problem

	    bool findOptimalSol = subProblemDual (xVal, yVal, n, piVal,
						  unbdPiVal);

	    if (findOptimalSol)
	      {
		/*
		 if (PRINT_VAR_VALUE)
		 {
		 for (i = 0; i < n; i++)
		 {
		 for (j = 0; j < n; j++)
		 if (piVal[i][j] > 0)
		 {
		 printf ("pi[%d,%d] %lf  ", i, j, piVal[i][j]);
		 //it's know that the best value of dual of SP is 0 when all pi[i][j]=0
		 //the unnecessary positive pi  is changed back to zero

		 //double coefficient = Q*xVal[i][j]-wt[i][j]*yVal[i][j];
		 // printf("coefficient of pi = %lf ", coefficient);
		 //if((coefficient>0 &&coefficient<EPSILON)||(coefficient<0 &&coefficient>EPSILON))
		 //piVal[i][j]=0;


		 }
		 }
		 printf ("\n");
		 }

		 GRBLinExpr expr = 0;
		 //expr += *zVar;
		 expr += zVar;
		 for (i = 0; i < n; i++)
		 for (j = 0; j < n; j++)
		 {
		 expr -= (price - cost_flow) * dis[i][j] * wt[i][j]
		 * yVar[i][j];
		 expr += cost_vehicle * vw * dis[i][j] * xVar[i][j];
		 expr -= (Q * xVar[i][j] - wt[i][j] * yVar[i][j])
		 * piVal[i][j];

		 }
		 addLazy (expr <= 0);
		 printf ("*** feasible sol cut is added ***\n");
		 */
		insertToHashmap (xVal, yVal, obj, &hashmap, n);
		cout << "hashmap size after finishing inserting: "
		    << hashmap.size () << endl;

		updateLB (xVal, yVal, piVal, n);
	      }
	    else
	      { //the dual of sub-problem is unbounded
		//the current x and y value will be cut after adding new cut
		//so no need to store x and y when dual of SP is unbounded
		//find an extreme ray
		//printf ("*** find an extreme ray ***\n");
		/*
		 for (i = 0; i < n; i++)
		 for (j = 0; j < n; j++)
		 {
		 if (Q * xVal[i][j] - wt[i][j] * yVal[i][j] < 0)
		 {
		 piVal[i][j] = 1;
		 //printf ("pi[%d,%d]=1\n", i, j, pi[i][j]);
		 }
		 else
		 piVal[i][j] = 0;
		 }
		 */

		GRBLinExpr expr = 0;
		for (i = 0; i < n; i++)
		  for (j = 0; j < n; j++)
		    expr += (Q * xVar[i][j] - wt[i][j] * yVar[i][j])
			* unbdPiVal[i][j];
		addLazy (expr >= 0);
		printf ("*** unbounded sol cut is added ***\n");

		double debugConV = 0;
		for (i = 0; i < n; i++)
		  for (j = 0; j < n; j++)
		    {
		      debugConV += (Q * xVal[i][j] - wt[i][j] * yVal[i][j])
			  * unbdPiVal[i][j];
		    }
		printf ("======> debugConV = %lf\n", debugConV);

		UB = BIG_M;
	      }

	    //initialSolution = false;

	    //for (auto it = hashmap.begin (); it != hashmap.end (); it++)
	    /*
	     for (const auto &entry : hashmap)
	     {
	     printf ("obj=%lf\n", entry.first);
	     for (i = 0; i < n; i++)
	     for (j = 0; j < n; j++)
	     if (entry.second.x[i][j] > 0)
	     printf ("%d->%d %.1lf  ", i, j, entry.second.x[i][j]);
	     printf ("\n");
	     }
	     */
	  }
	else if (where == GRB_CB_MIPNODE)
	  {
	    /***** DO NOT DELETE THIS PART. IT'S DIFFERENT WITH GRB_CB_MIPSOL *******/
	    // MIP node callback
	    int nodecnt = (int) getDoubleInfo (GRB_CB_MIPNODE_NODCNT);
	    double objbst = getDoubleInfo (GRB_CB_MIPNODE_OBJBST);
	    double objbnd = getDoubleInfo (GRB_CB_MIPNODE_OBJBND);
	    int solcnt = getIntInfo (GRB_CB_MIPNODE_SOLCNT);

	    if (PRINT_MIPSOL_MIPNODE_TITLE)
	      printf ("==========GRB_CB_GRB_CB_MIPNODE (new node)=========\n");
	    if (PRINT_NODE_INFO)
	      {
		printf ("nodecnt=%d\n", nodecnt);
		printf ("objbst=%lf\n", objbst);
		printf ("objbnd=%lf\n", objbnd);
		printf ("solcnt=%d\n", solcnt);
	      }
	    if (PRINT_UB_LB)
	      {
		printf ("======> UB %lf\n", UB);
		printf ("======> LB %lf\n", LB);
	      }

	    double gap = (objbnd - objbst) / objbnd;
	    printf ("gap=%lf\n", gap);

	    //======> call sub-problem

	    bool findObjInHashmap = false;
	    double objInHashmap;
	    for (const auto &entry : hashmap)
	      {
		if (mod (entry.first) == objbst)
		  {
		    findObjInHashmap = true;
		    objInHashmap = entry.first;
		  }
	      }

	    if (gap < OBJ_GAP_TOLERANCE && findObjInHashmap)
	      {
		cout << objbst << " found in hashmap!" << endl;

		UB = objbst;

		if (UB - LB < UB_LB_TOLERANCE)
		  {
		    printf ("======> UB(%lf) - LB(%lf) < %lf\n", UB, LB,
		    UB_LB_TOLERANCE);
		    printf ("terminate optimizing\n");
		    abort ();
		  }
		//read solution
		for (i = 0; i < n; i++)
		  for (j = 0; j < n; j++)
		    {
		      xVal[i][j] = hashmap[objbst].x[i][j];
		      yVal[i][j] = hashmap[objbst].y[i][j];
		    }

		if (PRINT_VAR_OF_OBJBST_HASHMAP)
		  {
		    printf ("x vars: ");
		    for (i = 0; i < n; i++)
		      for (j = 0; j < n; j++)
			if (xVal[i][j] > EPSILON)
			  printf ("%d->%d %lf  ", i, j, xVal[i][j]);
		    printf ("\ny vars: ");
		    for (i = 0; i < n; i++)
		      for (j = 0; j < n; j++)
			if (yVal[i][j] > EPSILON)
			  printf ("%d->%d %lf  ", i, j, yVal[i][j]);
		    printf ("\n");
		  }

		bool findOptimalSol = subProblemDual (xVal, yVal, n, piVal,
						      unbdPiVal);
		if (findOptimalSol)
		  {
		    /*
		     if (PRINT_VAR_VALUE)
		     {
		     for (i = 0; i < n; i++)
		     for (j = 0; j < n; j++)
		     {
		     if (piVal[i][j] > 0)
		     printf ("pi[%d,%d] ", i, j, piVal[i][j]);
		     printf ("\n");
		     }
		     }

		     GRBLinExpr expr = 0;
		     //expr += *zVar;
		     expr += zVar;
		     for (i = 0; i < n; i++)
		     for (j = 0; j < n; j++)
		     {
		     expr -= (price - cost_flow) * dis[i][j] * wt[i][j]
		     * yVar[i][j];
		     expr += cost_vehicle * vw * dis[i][j] * xVar[i][j];
		     expr -= (Q * xVar[i][j] - wt[i][j] * yVar[i][j])
		     * piVal[i][j];
		     }
		     addLazy (expr <= 0);
		     printf ("*** fesible sol cut is added ***\n");
		     */
		    updateLB (xVal, yVal, piVal, n);
		  }
		else
		  { //the sub-problem is unbounded
		    //find an extreme ray
		    //printf ("*** find an extreme ray ***\n");
		    /*
		     for (i = 0; i < n; i++)
		     for (j = 0; j < n; j++)
		     {
		     if (Q * xVal[i][j] - wt[i][j] * yVal[i][j] < 0)
		     {
		     piVal[i][j] = 1;
		     //printf ("pi[%d,%d]=1\n", i, j, pi[i][j]);
		     }
		     else
		     piVal[i][j] = 0;
		     }
		     */
		    GRBLinExpr expr = 0;
		    for (i = 0; i < n; i++)
		      for (j = 0; j < n; j++)
			expr += (Q * xVar[i][j] - wt[i][j] * yVar[i][j])
			    * unbdPiVal[i][j];
		    addLazy (expr >= 0);
		    printf ("*** unbounded sol cut is added ***\n");

		    double debugConV = 0;
		    for (i = 0; i < n; i++)
		      for (j = 0; j < n; j++)
			{
			  debugConV += (Q * xVal[i][j] - wt[i][j] * yVal[i][j])
			      * unbdPiVal[i][j];
			}
		    printf ("======> debugConV = %lf\n", debugConV);
		    //remove the objbst info from hashmap
		    //since the newly added cut will make this solution infeasible
		    hashmap.erase (objInHashmap);
		    printf ("======> objInHashmap = %lf is removed\n",
			    objInHashmap);

		    UB = BIG_M;
		  }

		/*
		 printf ("x vars: ");
		 for (i = 0; i < n; i++)
		 {
		 //getSolution in NOT GRB_OPTIMAL in MIPNODE
		 //might have x and y value over 1, like 19, 154.75 etc very strange value!!!
		 xVal[i] = getSolution (xVar[i], n);
		 //xVal[i] = getNodeRel (xVar[i], n);

		 for (j = 0; j < n; j++)
		 if (xVal[i][j] > 0.000001)
		 printf ("%d->%d %lf  ", i, j, xVal[i][j]);

		 }
		 printf ("\ny vars: ");
		 for (i = 0; i < n; i++)
		 {
		 yVal[i] = getSolution (yVar[i], n);
		 //yVal[i] = getNodeRel (yVar[i], n);
		 for (j = 0; j < n; j++)
		 if (yVal[i][j] > 0.000001)
		 printf ("%d->%d %lf  ", i, j, yVal[i][j]);
		 }
		 printf ("\n");
		 */

	      }
	  }
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
	cout << "Error during callback" << endl;
      }
  }
}
;

void
updateLB (double **xVal, double **yVal, double **piVal, int n)
{
  double newprofit = 0;
  int i, j;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	newprofit += (price - cost_flow) * dis[i][j] * wt[i][j] * yVal[i][j];
	newprofit -= cost_vehicle * vw * dis[i][j] * xVal[i][j];
	newprofit += (Q * xVal[i][j] - wt[i][j] * yVal[i][j]) * piVal[i][j];
      }
  if (newprofit > LB)
    LB = newprofit;
}

double
mod (double entry)
{
  int i;
  double reminder;
  for (i = 0; i < BIG_M; i++)
    {
      if (entry - i * BIG_M < BIG_M)
	{
	  reminder = entry - i * BIG_M;
	  break;
	}
    }
  if (reminder <= 0)
    {
      printf ("the mod of the entry %lf can't be negative (%lf).\n", entry,
	      reminder);
      exit (1);
    }
  return reminder;
}

void
insertToHashmap (double **xVal, double **yVal, double obj,
		 unordered_map<double, integerSol> *hashmap, int n)
{
  int i, j;
  struct integerSol sol;
  sol.x = new double*[n];
  sol.y = new double*[n];
  for (i = 0; i < n; i++)
    {
      sol.x[i] = new double[n];
      sol.y[i] = new double[n];
      for (j = 0; j < n; j++)
	{
	  sol.x[i][j] = xVal[i][j];
	  sol.y[i][j] = yVal[i][j];
	}
    }

  if (hashmap->find (obj) != hashmap->end ())
    {		    //if objbst alreay exists in hashmap
      int count = 0;
      for (const auto &entry : (*hashmap))
	{
	  double key = entry.first;
	  bool sameEntries = true;
	  //if obj is alreay in hashmap
	  //but x or y vars are different
	  //then new obj is added as BIG_M+obj
	  //if there are obj and BIG_M+obj exist in hashmap
	  //add new one as 2*BIG_M+obj, and so on
	  if (mod (key) == obj)
	    {

	      for (i = 0; i < n; i++)
		{
		  for (j = 0; j < n; j++)
		    if (entry.second.x[i][j] != sol.x[i][j]
			|| entry.second.y[i][j] != sol.y[i][j])
		      {
			sameEntries = false;
			break;
		      }
		  if (!sameEntries)
		    break;
		}
	      if (!sameEntries)
		count++;
	    }
	}

      hashmap->insert (
	{ (count + 1) * BIG_M + obj, sol });
    }
  else
    //if obj is not in hashmap now
    hashmap->insert (
      { obj, sol });

  cout << "the new integer solution is added to hashmap" << endl;
  cout << "hashmap size: " << hashmap->size () << endl;
}

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
		double **unbdPiVal)
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
      cout << "Find optimal soltuion." << endl;
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
      cout << "Model is unbounded" << endl;
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

	  LB = -BIG_M;
	  UB = BIG_M;

	  //======> start calling Gurobi <======

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

	      //extra cut
	      GRBLinExpr expr8 = z;
	      for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		  {

		    expr8 -= (price - cost_flow) * dis[i][j] * wt[i][j]
			* y[i][j];
		    expr8 += cost_vehicle * vw * dis[i][j] * x[i][j];

		  }
	      model.addConstr (expr8 <= 0, "cut_");

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

	      printf ("price = %lf\n", price);
	      printf ("cost_flow = %lf\n", cost_flow);
	      printf ("cost_vehicle = %lf\n", cost_vehicle);
	      printf ("vw = %lf\n", vw);

	      // Set callback function
	      //printIntSol cb = printIntSol (x, y, &z, n);
	      printIntSol cb = printIntSol (x, y, z, n);
	      model.setCallback (&cb);

	      // Must set LazyConstraints parameter when using lazy constraints
	      model.set (GRB_IntParam_LazyConstraints, 1);
	      //model.set (GRB_IntParam_Presolve, 0);
	      model.set (GRB_IntParam_Cuts, 0);
	      model.set (GRB_DoubleParam_Heuristics, 0);
	      model.set (GRB_IntParam_Threads, 1);
	      model.set (GRB_DoubleParam_MIPGap, 0);

	      //set initial solution
	      //the optimal solution for t20_01_data
	      /*
	       for (i = 0; i < n; i++)
	       for (j = 0; j < n; j++)
	       {
	       x[i][j].set (GRB_DoubleAttr_Start, 0);
	       y[i][j].set (GRB_DoubleAttr_Start, 0);
	       }
	       x[0][6].set (GRB_DoubleAttr_Start, 1);
	       x[6][5].set (GRB_DoubleAttr_Start, 1);
	       x[5][4].set (GRB_DoubleAttr_Start, 1);
	       x[4][14].set (GRB_DoubleAttr_Start, 1);
	       x[14][19].set (GRB_DoubleAttr_Start, 1);
	       y[0][6].set (GRB_DoubleAttr_Start, 1);
	       y[6][5].set (GRB_DoubleAttr_Start, 1);
	       y[5][4].set (GRB_DoubleAttr_Start, 1);
	       y[4][14].set (GRB_DoubleAttr_Start, 1);
	       y[14][19].set (GRB_DoubleAttr_Start, 1);
	       */
	      // Optimize model
	      model.optimize ();

	      //write model to file
	      model.write ("BPMP.lp");

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
