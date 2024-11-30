//This is the example from openmp official website
//https://www.openmp.org/wp-content/uploads/Intro_To_OpenMP_Mattson.pdf
#include <omp.h>
static long num_steps = 10000000;
double step;
#define NUM_THREADS 3
#include <sstream>
#include <iostream>

#include <chrono> //to measure run time

using namespace std;
using namespace std::chrono;

int
main ()
{
  auto beginWallClock = high_resolution_clock::now ();
  int i, nthreads;
  double pi;
  //Promote scalar to an array dimensioned by number of threads to avoid race condition.
  double sum[NUM_THREADS];
  step = 1.0 / (double) num_steps;
  omp_set_num_threads (NUM_THREADS);
#pragma omp parallel
    {
      int i, id, nthrds;
      double x;
      id = omp_get_thread_num ();
      nthrds = omp_get_num_threads ();
      cout << "id=" << id << endl;
      cout << "nthrds=" << nthrds << endl;
      //Only one thread should copy the number of threads to the global value to make sure multiple threads writing to the same address don’t conflict.
      if (id == 0)
	nthreads = nthrds;
      //This is a common trick in SPMD programs to create a cyclic distribution of loop iterations
      //Angela: for id=0, the indices of i will be 0, 2, 4, 6, ... (nthrds=2)
      //for id=1; the indices are 1, 3, 5, 7, ...
      for (i = id, sum[id] = 0.0; i < num_steps; i = i + nthrds)
	{
	  x = (i + 0.5) * step;
	  sum[id] += 4.0 / (1.0 + x * x);
	}
    }
  for (i = 0, pi = 0.0; i < nthreads; i++)
    pi += sum[i] * step;

  cout << "pi = " << pi << endl;

  auto endWallClock = high_resolution_clock::now ();
  auto elapsedWallClock = duration_cast < std::chrono::nanoseconds
      > (endWallClock - beginWallClock);
  printf ("SolutionTimeWallClock = %.3f seconds\n",
	  elapsedWallClock.count () * 1e-9);

  return 1;
}
