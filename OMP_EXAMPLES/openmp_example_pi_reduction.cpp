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

  int i;
  double pi;
  double x, sum;
  step = 1.0 / (double) num_steps;
  omp_set_num_threads (NUM_THREADS);
#pragma omp parallel
    {
      int id, nthrds;

      double x;

#pragma omp for reduction(+:sum)
      for (i = 0; i < num_steps; i++)
	{
	  x = (i + 0.5) * step;
	  sum = sum + 4.0 / (1.0 + x * x);
	}
      pi = step * sum;
      id = omp_get_thread_num ();
      nthrds = omp_get_num_threads ();
      cout << "id=" << id << endl;
      cout << "nthrds=" << nthrds << endl;
    }

  cout << "pi = " << pi << endl;

  auto endWallClock = high_resolution_clock::now ();
  auto elapsedWallClock = duration_cast < std::chrono::nanoseconds
      > (endWallClock - beginWallClock);
  printf ("SolutionTimeWallClock = %.3f seconds\n",
	  elapsedWallClock.count () * 1e-9);

  return 1;
}
