#include <omp.h>
#include <sstream>
#include <iostream>

#include <chrono> //to measure run time

using namespace std;
using namespace std::chrono;

static long num_steps = 10000000;
double step;
#define PAD 8 // assume 64 byte L1 cache line size
#define NUM_THREADS 1
int
main ()
{
  auto beginWallClock = high_resolution_clock::now ();
  int i, nthreads;
  double pi;
  //Pad the array so each sum value is in a different cache line
  double sum[NUM_THREADS][PAD];
  step = 1.0 / (double) num_steps;
  omp_set_num_threads (NUM_THREADS);
#pragma omp parallel
    {
      int i, id, nthrds;
      double x;
      id = omp_get_thread_num ();
      nthrds = omp_get_num_threads ();
      if (id == 0)
	nthreads = nthrds;
      for (i = id, sum[id][0] = 0.0; i < num_steps; i = i + nthrds)
	{
	  x = (i + 0.5) * step;
	  sum[id][0] += 4.0 / (1.0 + x * x);
	}
    }
  for (i = 0, pi = 0.0; i < nthreads; i++)
    pi += sum[i][0] * step;

  cout << "pi = " << pi << endl;

  auto endWallClock = high_resolution_clock::now ();
  auto elapsedWallClock = duration_cast < std::chrono::nanoseconds
      > (endWallClock - beginWallClock);
  printf ("SolutionTimeWallClock = %.3f seconds\n",
	  elapsedWallClock.count () * 1e-9);

  return 1;
}
