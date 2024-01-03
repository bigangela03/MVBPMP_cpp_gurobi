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

int
main (int argc, char **argv)
{
  if (argc < 3)
    {
      printf ("Usage: ./main dataFolder output.csv. Terminating...\n");
      exit (1);
    }
  cout << argc << endl;
  for (int i = 1; i < argc; ++i)
    cout << argv[i] << endl;

  readData route;

  fstream csvFile (argv[2], ios::out); // Open the CSV file for writing

  // Check if the file was successfully opened
  if (!csvFile)
    {
      std::cerr << "Error opening a file." << std::endl;
      return 1;
    }

  csvFile
      << "Filename, Type of Algorithm, Profit, Profit after Local Search, CPU Time(ms), Run Time(ms), Route, Requests"
      << endl;

  auto it = filesystem::directory_iterator (argv[1]);

  for (const auto &entry : it)
    {
      filesystem::path path = entry.path ();
      if (entry.is_regular_file ())
	{
	  string filename = path.string ();
	  route.readSingleFile (filename);
	  printf ("===>file name: %s\n", filename.c_str ());

	  //route.printStats ();
	  //route.printData ();

	  int n = route.numOfNode;
	  double wt[n][n];
	  double dis[n][n];
	  for (int i = 1; i <= n; i++)
	    for (int j = 1; j <= n; j++)
	      {
		wt[i - 1][j - 1] = route.w[i][j];
		dis[i - 1][j - 1] = route.d[i][j];
	      }
	  for (int i = 0; i < n; i++)
	    {
	      for (int j = 0; j < n; j++)
		{
		  //printf("%.2f  ",dis[i][j]);
		}
	      //	printf("\n");
	    }



	  int numV = 3;

	  string vehicleFileName = "vehicle_location_data/vehicle_data_"
	      + to_string (n) + "_" + to_string (numV) + "_1.txt";
	  printf ("===>vehicle name: %s\n", vehicleFileName.c_str ());
	  //exit (1);
	  route.readSingleVehicleFile (vehicleFileName, numV);

	  //vehicle and origin always starts from 0
	  int vehicle[numV];
	  int origin[numV];

	  for (int i = 0; i < numV; i++)
	    {
	      vehicle[i] = route.vehicle[i];
	      //printf ("%d ", vehicle[i]);
	    }
	  //printf ("\n");
	  for (int i = 0; i < numV; i++)
	    {
	      origin[i] = route.origin[i];
	      //printf ("%d ", origin[i]);
	    }
	  //printf ("\n");
	}

    }
  return 0;
}
