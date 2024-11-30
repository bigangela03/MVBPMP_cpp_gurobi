#include <omp.h>
#include <sstream>
#include <iostream>

#include <chrono> //to measure run time

using namespace std;
using namespace std::chrono;

static long num_steps = 10000000;
double step;
#define NUM_THREADS 4
int
main ()
{
  auto beginWallClock = high_resolution_clock::now ();

  double pi;
  step = 1.0 / (double) num_steps;
  omp_set_num_threads (NUM_THREADS);
  int nthreads;
#pragma omp parallel
    {
      int i, id, nthrds;
      double x, sum; //Create a scalar local to each thread to accumulate partial sums.
      id = omp_get_thread_num ();
      nthrds = omp_get_num_threads ();
      if (id == 0)
	nthreads = nthrds;
      id = omp_get_thread_num ();
      nthrds = omp_get_num_threads ();
      for (i = id, sum = 0.0; i < num_steps; i = i + nthreads)
	{
	  x = (i + 0.5) * step;
	  sum += 4.0 / (1.0 + x * x);
	}
#pragma omp critical
      pi += sum * step;
      cout<<"id="<<id<<endl;
      cout<<"sum="<<sum<<endl;
    }

  cout << "pi = " << pi << endl;

  auto endWallClock = high_resolution_clock::now ();
  auto elapsedWallClock = duration_cast < std::chrono::nanoseconds
      > (endWallClock - beginWallClock);
  printf ("SolutionTimeWallClock = %.3f seconds\n",
	  elapsedWallClock.count () * 1e-9);

  return 1;
}
