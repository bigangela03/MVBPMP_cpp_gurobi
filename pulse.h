#ifndef PULSE_H
#define PULSE_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> //to get path to current directory
#include <math.h>
#include <time.h>
#include <malloc.h>
#include <string.h>

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
#include <unordered_map>
#include <unordered_set>

#include <limits> //to get the biggest value of double variables

vector<int> bestFoundRoutInPulse;
double bestFoundObjInPulse;
bool PRINT_FOR_DEBUG_PULSE = false;
double **B;
int numBoundingSections;
double *BsectionIndexLowerDist;

void pulse(int, vector<vector<int>> &, double, double, double **, double, vector<int>, vector<vector<double>> &, int);

void pulseAlgorithm(int n, double **dis, vector<vector<double>> xCoeff, double disLimit, vector<int> &nodesOnRoute,
                    int startNode, int endNode, vector<vector<int>> nodeNeighbors, double delta)
{
    int bigM = 1000000;
    double minDistForCheckingBounds = 8;

    vector<vector<int>> pathsToEndNode;
    vector<double> pathsObjToEndNode;
    vector<double> pathsDistToEndNode;

    vector<int> onePartialPath;
    double partialDist = 0;
    double partialObj = 0;

    // for example, delta=5, minDistForCheckingBounds=9, numBoundingSections=floor((20-9)/5)=2
    // B[1][2]=-2.3 means node distance section 1, for node 2 (node index starts from 0)
    //                    node 0, node 1, node 2,..., node n-1
    // 0 (dist>=10, <15)   -3.5    -1.8    ...
    // 1 (dist>=15, <20)     ..     ..     -2.3

    numBoundingSections = floor((disLimit - minDistForCheckingBounds) / delta);

    if (numBoundingSections < 1)
    {
        cout << "ERROR: in pulse, numBoundingSections is at lease one!";
        exit(1);
    }

    BsectionIndexLowerDist = new double[numBoundingSections];
    for (int i = 0; i < numBoundingSections; i++)
        BsectionIndexLowerDist[i] = disLimit - (numBoundingSections - i) * delta;

    // if (PRINT_FOR_DEBUG_PULSE)
    // {
    //     cout << "Pulse algorithm: numBoundingSections=" << numBoundingSections << endl;
    //     for (int i = 0; i < numBoundingSections; i++)
    //         cout << BsectionIndexLowerDist[i] << endl;
    // }

    // generate a numBoundingSections rows, n columns matrix B.
    // B[1][2]=-2.3 means node distance section 1, for node 2 (node index starts from 0)
    //                    node 0, node 1, node 2,..., node n-1
    // 0 (dist>=10, <15)   -3.5    -1.8    ...
    // 1 (dist>=15, <20)     ..     ..     -2.3

    // right now B matrix has no value assigned. And this procedure is to find value for B.
    // But pulse use B value, so we have to initialize B with some value. Use -bigM, so no route will be pruned because of Bound check

    B = new double *[numBoundingSections];
    for (int i = 0; i < numBoundingSections; i++)
        B[i] = new double[n];

    for (int i = 0; i < numBoundingSections; i++)
        for (int j = 0; j < n; j++)
            B[i][j] = -bigM;

    for (int i = 0; i < numBoundingSections; i++)
    {
        double currentDisforB = BsectionIndexLowerDist[i];
        for (int node = 0; node < n; node++)
        {
            // if (node == startNode || node == endNode)
            //     continue;

            partialObj = 0;
            partialDist = currentDisforB;
            if (node == startNode)
                partialDist = 0;

            bestFoundObjInPulse = bigM;
            bestFoundRoutInPulse.clear();

            onePartialPath.clear();
            onePartialPath.push_back(node);
            for (int neighbor : nodeNeighbors[node])
                pulse(neighbor, nodeNeighbors, partialObj, partialDist, dis, disLimit, onePartialPath, xCoeff, endNode);

            B[i][node] = bestFoundObjInPulse; // if no route found, B[i][node]=bigM as bestFoundObjInPulse's initial value
        }
    }

    // if (PRINT_FOR_DEBUG_PULSE)
    // {
    //     cout << "===> B matrix" << endl;
    //     for (int i = 0; i < numBoundingSections; i++)
    //     {
    //         cout << BsectionIndexLowerDist[i] << ": ";
    //         for (int j = 0; j < n; j++)
    //             cout << B[i][j] << " ";
    //         cout << endl;
    //     }
    // }

    // pulse(startNode, nodeNeighbors, partialObj, partialDist, dis, disLimit, onePartialPath, xCoeff, endNode);
    // cout << "===> start searching" << endl;
    partialObj = 0;
    partialDist = 0;
    onePartialPath.clear();
    bestFoundObjInPulse = bigM;
    bestFoundRoutInPulse.clear();
    onePartialPath.push_back(startNode);
    for (int neighbor : nodeNeighbors[startNode])
        pulse(neighbor, nodeNeighbors, partialObj, partialDist, dis, disLimit, onePartialPath, xCoeff, endNode);
}

void pulse(int currentNode, vector<vector<int>> &nodeNeighbors, double partialObj, double partialDist, double **dis,
           double disLimit, vector<int> onePartialPath, vector<vector<double>> &xCoeff, int endNode)
{
    // cout << "currentNode=" << currentNode << "; partialObj=" << partialObj << "; partialDist=" << partialDist << endl;

    // check Feasibility
    // if (onePartialPath.empty())
    // {
    //     onePartialPath.push_back(currentNode);
    //     for (int neighbor : nodeNeighbors[currentNode])
    //         pulse(neighbor, nodeNeighbors, partialObj, partialDist, dis, disLimit, onePartialPath, xCoeff, endNode);

    //     return;
    // }

    int numNodes = onePartialPath.size();
    int lastNode = onePartialPath[numNodes - 1];

    double distTemp = partialDist + dis[lastNode][currentNode];
    double objTemp = partialObj + xCoeff[lastNode][currentNode];

    if (currentNode == endNode)
    {
        if (objTemp < bestFoundObjInPulse)
        {
            bestFoundObjInPulse = objTemp;
            onePartialPath.push_back(currentNode);
            bestFoundRoutInPulse = onePartialPath;
            if (PRINT_FOR_DEBUG_PULSE)
            {
                cout << "===> find a better obj " << bestFoundObjInPulse << endl;
                for (int &e : bestFoundRoutInPulse)
                    cout << e << " ";
                cout << endl;
            }
        }
        return;
    }

    if (distTemp + dis[currentNode][endNode] > disLimit)
    {
        // if (PRINT_FOR_DEBUG_PULSE)
        //     cout << "===> distance is over limit. return." << endl;
        return;
    }

    // 2. check bound
    // numBoundingSections is at lease one
    if (distTemp >= BsectionIndexLowerDist[0])
    {
        if (distTemp >= BsectionIndexLowerDist[numBoundingSections - 1])
        {
            if (objTemp + B[numBoundingSections - 1][currentNode] >= bestFoundObjInPulse)
                return;
        }
        else
            for (int i = 1; i < numBoundingSections; i++)
            {
                if ((distTemp >= BsectionIndexLowerDist[i - 1]) && (distTemp < BsectionIndexLowerDist[i]))
                {
                    if (objTemp + B[i - 1][currentNode] >= bestFoundObjInPulse)
                        return;
                    else
                        break;
                }
            }
    }

    // 3. rollback
    if (numNodes >= 2)
    {
        int secondLastNode = onePartialPath[numNodes - 2];
        // if (PRINT_FOR_DEBUG_PULSE)
        //     cout << "secondLastNode=" << secondLastNode << endl;

        // if (find(nodeNeighbors[secondLastNode].begin(), nodeNeighbors[secondLastNode].end(), currentNode) != nodeNeighbors[secondLastNode].end())
        {
            // cost of route from node i(secondLastNode) directly to j (currentNode)
            double directObj = partialObj - xCoeff[secondLastNode][lastNode] + xCoeff[secondLastNode][currentNode];

            if (directObj < objTemp)
            {
                // if (PRINT_FOR_DEBUG_PULSE)
                //     cout << "===> rollback. return." << endl;
                return;
            }
        }
    }

    onePartialPath.push_back(currentNode);

    for (int neighbor : nodeNeighbors[currentNode])
    {
        // if neighbor is not in path
        if (find(onePartialPath.begin(), onePartialPath.end(), neighbor) == onePartialPath.end())
        {
            // if (PRINT_FOR_DEBUG_PULSE)
            //     cout << "start checking node " << currentNode << "'s neighbor " << neighbor << endl;
            pulse(neighbor, nodeNeighbors, objTemp, distTemp, dis, disLimit, onePartialPath, xCoeff, endNode);
        }
    }
}

#endif