/* Copyright 2023, Gurobi Optimization, LLC */

/* Solve a traveling salesman problem on a randomly generated set of
 points using lazy constraints.   The base MIP model only includes
 'degree-2' constraints, requiring each node to have exactly
 two incident edges.  Solutions to this model may contain subtours -
 tours that don't visit every node.  The lazy constraint callback
 adds new constraints to cut them off. */

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

int numNodes; //number of nodes in instance

string
itos (int i)
{
  stringstream s;
  s << i;
  return s.str ();
}

void
//func (vector<vector<vector<vector<double>>>> &b)
//func (double (*b)[numNodes][numNodes][numNodes]) //can't use vars for dimension
func1 (double (*b)[3][3][3])
{
  printf ("in func1: i=%lf\n", b[0][1][0][1]);
}

void
func2 (double **b)
{
  printf ("in func2: %lf\n", b[0][1]);
}

void
func3 (double ****b)
{
  printf ("in function3: %lf\n", b[0][1][1][0]);
}

int
main (int argc, char *argv[])
{
  int i, j, k, q;
  int n;
  printf ("i am here\n");
  n = stoi (argv[1]);
  //n = 3;
  //double a[n][n][n][n];
  //vector < vector < vector<vector<double>> >> a;
  //double a[n][n][n][n];
  double a[3][3][3][3];
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
	for (q = 0; q < n; q++)
	  a[i][j][k][q] = i + j + k + q;

  func1 (a);

  //---------------------------------------
  double **dynamicArray = NULL;
  int ROWS = 2;
  int COLUMNS = 3;
  // memory allocated for elements of rows.
  dynamicArray = new double*[ROWS];
  // memory allocated for  elements of each column.
  for (int i = 0; i < ROWS; i++)
    {
      dynamicArray[i] = new double[COLUMNS];
    }
  dynamicArray[0][1] = 10;
  func2 (dynamicArray);
  // free the allocated memory
  for (int i = 0; i < ROWS; i++)
    {
      delete[] dynamicArray[i];
    }
  delete[] dynamicArray;

  //---------------------------------------
  double ****array = NULL;
  numNodes=stoi(argv[1]);
  int imax=numNodes;
  int jmax=numNodes;
  int kmax=numNodes;
  int qmax=numNodes;
  array = new double***[imax];
  for (i = 0; i < imax; i++)
    {
      array[i] = new double**[jmax];
      for (j = 0; j < jmax; j++)
	{
	  array[i][j] = new double*[kmax];
	  for (k = 0; k < kmax; k++)
	 	{
	 	  array[i][j][k] = new double[qmax];
	 	}
	}
    }

  for (i = 0; i < imax; i++)
     for (j = 0; j < jmax; j++)
       for (k = 0; k < kmax; k++)
 	for (q = 0; q < qmax; q++)
 	  array[i][j][k][q] = i + j + k + q;
  func3(array);

}
