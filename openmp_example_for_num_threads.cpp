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

 omp_set_dynamic (0);     // Explicitly disable dynamic teams
 //#pragma omp parallel for
#pragma omp parallel for num_threads(NUM_THREADS)
  for (i = 0; i < 3; i++)
    {

      int id = omp_get_thread_num ();
      int nthrds = omp_get_num_threads ();

#pragma omp critical
      cout << "--------------------" << endl;
      cout << "i=" << i << endl;
      cout << "id=" << id << endl;
      cout << "nthrds=" << nthrds << endl;
    }

  auto endWallClock = high_resolution_clock::now ();
  auto elapsedWallClock = duration_cast < std::chrono::nanoseconds
      > (endWallClock - beginWallClock);
  printf ("SolutionTimeWallClock = %.3f seconds\n",
	  elapsedWallClock.count () * 1e-9);

  return 1;
}
