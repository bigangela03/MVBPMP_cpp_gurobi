#ifndef READDATA_H
#define READDATA_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <filesystem>

using namespace std;
using matrix = vector<vector<double>>;

class readData
{

public:
  // the weight in tones of a load that a potential customer would like to ship
  // from node i to node j as w(ij)
  // the vehicle weight v tons when empty
  // carrying capacity of Q tons
  // incurs a travel cost of c dollars per mile per ton
  // empty vehicle incurs a cost of c v d(ij) when traversing arc(i,j)
  // where as fully loaded vehicle incurs a cost of c(v+Q) d(ij)
  // the customer pays p d w(ij) where p is the price charged (revenue received)
  // in dollar per mile per ton
  // DIS

  int numOfNode = 0;
  int travelCost = 0;
  int totalCapacity = 0;
  int DIS = 0;
  double vehicleWeight = 0;
  double priceCharged = 0;
  matrix w;
  matrix d;

  //int vehicle[];
  //int origin[][];
  int *vehicle = NULL;
  int *origin = NULL;

  void
  readFile (const string &path)
  {
    auto it = filesystem::directory_iterator (path);

    for (const auto &entry : it)
      {
	filesystem::path path = entry.path ();
	if (entry.is_regular_file ())
	  {
	    readSingleFile (path.string ());
	  }
      }
  }

  /**
   * This function read txt file
   * @param filename
   */
  void
  readSingleFile (const string &filename)
  {
    ifstream file (filename);
    if (!file.is_open ())
      {
	cerr << "Error in opening file: " << filename << endl;
	return;
      }

    w.clear (); // clear the matrix if the stores incorrect values
    d.clear (); // clear the matrix if the stores incorrect values

    string line;
    while (getline (file, line))
      {
	istringstream iss (line);
	string param, value, assignSymbol, n, p, c, Q, v, dis;
	string variable = "";
	int pos = 0;
	int check = 0;
	bool isWeight, isDemand;

	if (line.empty () || line[0] == '#')
	  {
	    // Skip empty lines or lines starting with '#'
	    continue;
	  }

	if (iss >> param >> value >> assignSymbol && assignSymbol == ":=")
	  {
	    if (param == "param" && value == "n")
	      {
		iss >> n;
		numOfNode = stoi (n);
	      }
	    else if (param == "param" && value == "p")
	      {
		iss >> p;
		priceCharged = stod (p);
	      }
	    else if (param == "param" && value == "c")
	      {
		iss >> c;
		travelCost = stoi (c);
	      }
	    else if (param == "param" && value == "Q")
	      {
		iss >> Q;
		totalCapacity = stoi (Q);
	      }
	    else if (param == "param" && value == "v")
	      {
		iss >> v;
		vehicleWeight = stod (v);
	      }
	    else if (param == "param" && value == "DIS")
	      {
		iss >> dis;
		DIS = stoi (dis);
	      }
	    else if (param == "param" && value == "w")
	      {
		isWeight = true;
		isDemand = false;
	      }
	    else if (param == "param" && value == "d")
	      {
		isWeight = false;
		isDemand = true;
	      }
	  }
	else if ((pos = line.find (":")))
	  {
	    if (pos != string::npos)
	      {
		param = line.substr (0, pos);
		variable = line.substr (pos + 2);
		variable.erase (0, value.find_first_not_of (";"));

		// Process the parameter and value as needed
		if (param == "param n")
		  {
		    numOfNode = stoi (variable);
		  }
		else if (param == "param p")
		  {
		    priceCharged = stod (variable);
		  }
		else if (param == "param c")
		  {
		    travelCost = stoi (variable);
		  }
		else if (param == "param Q")
		  {
		    totalCapacity = stoi (variable);
		  }
		else if (param == "param v")
		  {
		    vehicleWeight = stod (variable);
		  }
		else if (param == "param DIS")
		  {
		    DIS = stoi (variable);
		  }
		else if (param == "param w")
		  {
		    isWeight = true;
		    isDemand = false;
		    check = 1;
		  }
		else if (param == "param d")
		  {
		    isWeight = false;
		    isDemand = true;
		    check = 1;
		  }
	      }

	    vector < string > tokens = splitString (line, '\t');

	    if (tokens.size () == 3 && isWeight && !isDemand)
	      {
		int intI = stoi (tokens[0]);
		int intJ = stoi (tokens[1]);
		double doubleWeight = stod (tokens[2]);

		if (intI >= w.size ())
		  w.resize (intI + 1);
		if (intJ >= w[intI].size ())
		  w[intI].resize (intJ + 1);

		w[intI][intJ] = doubleWeight;
	      }
	    else if (tokens.size () > 3 && check != 1 && isWeight && !isDemand)
	      {
		vector<double> j;

		if (w.size () == 0)
		  {
		    j.push_back (0);
		    w.push_back (j);
		  }
		else
		  j.push_back (0);

		for (auto it = tokens.begin () + 1; it != tokens.end (); ++it)
		  {
		    if (*it == "")
		      {
			continue;
		      }
		    else
		      j.push_back (std::stod (*it));
		  }

		w.push_back (j);
	      }
	    else if (tokens.size () == 3 && isDemand && !isWeight)
	      {
		int intI = stoi (tokens[0]);
		int intJ = stoi (tokens[1]);
		double doubleDemand = stod (tokens[2]);

		if (intI >= d.size ())
		  d.resize (intI + 1);
		if (intJ >= d[intI].size ())
		  d[intI].resize (intJ + 1);

		d[intI][intJ] = doubleDemand;
	      }
	    else if (tokens.size () > 3 && check != 1 && !isWeight && isDemand)
	      {
		vector<double> j;
		if (d.size () == 0)
		  {
		    j.push_back (0);
		    d.push_back (j);
		  }
		else
		  j.push_back (0);

		for (auto it = tokens.begin () + 1; it != tokens.end (); ++it)
		  {
		    if (*it == "")
		      {
			continue;
		      }
		    else
		      j.push_back (std::stod (*it));
		  }
		d.push_back (j);
	      }
	  }
      }

    file.close ();
  }

  void
  readSingleVehicleFile (const string &filename, int numV)
  {

    ifstream file (filename);
    if (!file.is_open ())
      {
	cerr << "Error in opening file: " << filename << endl;
	return;
      }

    //int vehicle[];
    //int origin[][];

    //vehicles.clear ();
    //origin.clear ();
    vehicle = new int[numV];
    origin = new int[numV];

    string line;
    int countOrigins = 1;
    while (getline (file, line))
      {
	//printf ("------------ \n");
	istringstream iss (line);
	string param, value, assignSymbol;
	string variable = "";
	int pos = 0;
	int check = 0;
	bool isVehicle, isOrigin;

	//printf ("%s \n", line.c_str ());

	if (line.empty () || line[0] == '#')
	  {
	    // Skip empty lines or lines starting with '#'
	    continue;
	  }

	if (iss >> param >> value >> assignSymbol && assignSymbol == ":=")
	  {
	    if (param == "set" && value == "V")
	      {
		string s1, s2;
		string delimiter = ":=";
		if ((pos = line.find (delimiter)) != string::npos)
		  {
		    s1 = line.substr (pos + delimiter.length (),
				      line.length ());

		    delimiter = ";";
		    if ((pos = s1.find (delimiter)) != string::npos)
		      {
			s2 = s1.substr (0, pos);
			//printf ("%s \n", s2.c_str ());

			vector < string > tokens = splitString (s2, '\t');
			int tokensize = tokens.size ();
			int maxIndex = stoi (tokens[tokensize - 1]);
			if (maxIndex != numV)
			  {
			    printf ("the max (last one) vehicle index is %d\n",
				    maxIndex);
			    printf ("the number of vehicles is %d\n", numV);
			    printf (
				"These two numbers should be the same! Please check input data. Termininating...\n");
			    exit (1);
			  }
			for (int i = 0; i < tokensize; i++)
			  {
			    vehicle[i] = stoi (tokens[i]) - 1;
			    //printf ("%d\n", vehicle[i]);
			  }

		      }
		  }

	      }
	    else if (param == "param" && value == "origin")
	      {
		isOrigin = true;
		check = 1;
	      }
	  }
	else
	  {
	    vector < string > tokens = splitString (line, '\t');
	    //printf("--------- \n");

	    if (tokens.size () == 2)
	      {
		if (isOrigin)
		  {
		    //printf ("token0 %s \n", tokens[0].c_str ());
		    //printf ("token1 %s \n", tokens[1].c_str ());
		    int index = stoi (tokens[0]); //origin index
		    if (countOrigins == index)
		      {
			origin[index - 1] = stoi (tokens[1]) - 1;
			countOrigins++;
		      }
		    else
		      {
			printf (
			    "the index in origin data should be in format like 1 2 3 ...\n");
			printf ("Termininating...\n");
			exit (1);
		      }
		  }

	      }
	    else
	      {
		if (tokens[0] != ";")
		  {
		    printf (
			"origin data in vehicle file should have only two elements in each row\n");
		    printf ("the current has %lu elements\n", tokens.size ());
		    printf ("Termininating...\n");
		    exit (1);
		  }
	      }

	  }
      }

    file.close ();

  }

  /**
   * This function print out the stats if the file read in
   */
  void
  printStats ()
  {
    cout << "===============Stats===============" << endl;
    cout << "Number of Node: " << numOfNode << endl;
    cout << "Price charged: " << priceCharged << endl;
    cout << "Travel cost: " << travelCost << endl;
    cout << "Total Capacity: " << totalCapacity << endl;
    cout << "Vehicle Weight: " << vehicleWeight << endl;
    cout << "DIS Maximum distance that Vehicle travel: " << DIS << endl;
    cout << "Total W: " << (w.size () - 1) << endl;
    cout << "Total D: " << (d.size () - 1) << endl;
    cout << "===================================" << endl << endl;
  }

  void
  printData ()
  {

    cout << "Weight: " << endl;
    for (auto i : w)
      {
	for (auto j : i)
	  {
	    cout << j << " ";
	  }
	cout << endl;
      }

    cout << "Distance: " << endl;
    for (auto i : d)
      {
	for (auto j : i)
	  {
	    cout << j << " ";
	  }
	cout << endl;
      }
  }

  // string filename, algo, profit, cpu time, real time, route, requests

private:
  /**
   * This function split the string and return it as a vector of string seperate by delimiter
   * @param input string to split
   * @param delimiter which character to split at
   */
  vector<string>
  splitString (const string &input, char delimiter)
  {
    vector < string > tokens;
    istringstream iss (input);
    string token;

    // while (getline(iss, token, delimiter))
    // {
    //     token.erase(std::remove(token.begin(), token.end(), '\r'), token.end());
    //     tokens.push_back(token);
    // }

    while (iss >> token)
      {
	tokens.push_back (token);
      }

    return tokens;
  }

  static void
  printPathInfo (double &totalWeight, double &currDist, bool distRatio,
		 vector<int> path, int greedyNum)
  {
    cout << "Weight: " << totalWeight << endl;
    cout << "Distance: " << currDist << endl;
    cout << "Ratio: " << totalWeight / currDist << endl;
    cout << "Path for Greedy #" << greedyNum << " (distRatio = " << distRatio
	<< "): ";
    for (auto &c : path)
      cout << c << " ";
    cout << endl << endl;
  }
};

#endif
