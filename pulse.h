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

    BsectionIndexLowerDist = new double[numBoundingSections];
    for (int i = 0; i < numBoundingSections; i++)
        BsectionIndexLowerDist[i] = disLimit - (numBoundingSections - i) * delta;

    if (PRINT_FOR_DEBUG_PULSE)
    {
        cout << "numBoundingSections=" << numBoundingSections << endl;
        for (int i = 0; i < numBoundingSections; i++)
            cout << BsectionIndexLowerDist[i] << endl;
    }
    // generate a numBoundingSections rows, n columns matrix B.
    // B[1][2]=-2.3 means node distance section 1, for node 2 (node index starts from 0)
    //                    node 0, node 1, node 2,..., node n-1
    // 0 (dist>=10, <15)   -3.5    -1.8    ...
    // 1 (dist>=15, <20)     ..     ..     -2.3
    B = new double *[numBoundingSections];
    for (int i = 0; i < numBoundingSections; i++)
        B[i] = new double[n];

    // bound(B, BsectionIndexLowerDist);

    for (int i = 0; i < numBoundingSections; i++)
    {
        double currentDisforB = BsectionIndexLowerDist[i];
        for (int node = 0; node < n; node++)
        {
            if (node == startNode || node == endNode)
                continue;
            partialObj = 0;
            partialDist = currentDisforB;
            onePartialPath.clear();
            bestFoundObjInPulse = bigM;
            bestFoundRoutInPulse.clear();
            pulse(node, nodeNeighbors, partialObj, partialDist, dis, disLimit, onePartialPath, xCoeff, endNode);

            B[i][node] = bestFoundObjInPulse; // if no route found, B[i][node]=bigM as bestFoundObjInPulse's initial value
        }
    }
    partialObj = 0;
    partialDist = 0;
    onePartialPath.clear();
    bestFoundObjInPulse = bigM;
    bestFoundRoutInPulse.clear();

    if (PRINT_FOR_DEBUG_PULSE)
    {
        cout << "===> B matrix" << endl;
        for (int i = 0; i < numBoundingSections; i++)
        {
            cout << BsectionIndexLowerDist[i] << ": ";
            for (int j = 0; j < n; j++)
                cout << B[i][j] << " ";
            // printf("%12lf ", B[i][j]);
            cout << endl;
        }
    }
    pulse(startNode, nodeNeighbors, partialObj, partialDist, dis, disLimit, onePartialPath, xCoeff, endNode);
}

void pulse(int currentNode, vector<vector<int>> &nodeNeighbors, double partialObj, double partialDist, double **dis,
           double disLimit, vector<int> onePartialPath, vector<vector<double>> &xCoeff, int endNode)
{

    int numNodes = onePartialPath.size();
    int lastNode;

    if (PRINT_FOR_DEBUG_PULSE)
    {
        cout << "------pulse------" << endl;
        cout << "currentNode=" << currentNode << endl;
        cout << "numNodes=" << numNodes << endl;
    }

    // check Feasibility
    if (!onePartialPath.empty())
    {
        lastNode = onePartialPath[numNodes - 1];

        if (PRINT_FOR_DEBUG_PULSE)
            cout << "lastNode=" << lastNode << endl;

        if (partialDist + dis[lastNode][currentNode] > disLimit)
        {
            if (PRINT_FOR_DEBUG_PULSE)
                cout << "===> distance is over limit. return." << endl;
            return;
        }
    }
    else
    {
        onePartialPath.push_back(currentNode);
        for (int neighbor : nodeNeighbors[currentNode])
        {
            pulse(neighbor, nodeNeighbors, partialObj, partialDist, dis, disLimit, onePartialPath, xCoeff, endNode);
        }

        return;
    }

    // 2. check bound
    for (int i = 0; i < numBoundingSections; i++)
    {
        // cout << "===> check bound i=" << i << "; dist=" << BsectionIndexLowerDist[i] << endl;
        if (partialDist >= BsectionIndexLowerDist[i])
        {
            if (partialObj + B[i][lastNode] >= bestFoundObjInPulse)
                return;
            else
                break;
        }
    }

    // 3. rollback

    if (numNodes >= 2)
    {
        int secondLastNode = onePartialPath[numNodes - 2];
        if (PRINT_FOR_DEBUG_PULSE)
            cout << "secondLastNode=" << secondLastNode << endl;

        // if (find(nodeNeighbors[secondLastNode].begin(), nodeNeighbors[secondLastNode].end(), currentNode) != nodeNeighbors[secondLastNode].end())
        {
            // cost of route from node i(secondLastNode) directly to j (currentNode)
            double directObj = partialObj - xCoeff[secondLastNode][lastNode] + xCoeff[secondLastNode][currentNode];

            // cost of route from node i(secondLastNode) to k (lastNode) then j (currentNode)
            double detourObj = partialObj + xCoeff[lastNode][currentNode];

            if (directObj < detourObj)
            {
                if (PRINT_FOR_DEBUG_PULSE)
                    cout << "===> rollback. return." << endl;
                return;
            }
        }
    }

    onePartialPath.push_back(currentNode);
    partialObj += xCoeff[lastNode][currentNode];
    partialDist += dis[lastNode][currentNode];

    if (currentNode == endNode)
    {
        if (partialObj < bestFoundObjInPulse)
        {
            bestFoundObjInPulse = partialObj;
            bestFoundRoutInPulse = onePartialPath;
            if (PRINT_FOR_DEBUG_PULSE)
            {
                cout << "===> find a better obj " << bestFoundObjInPulse << endl;
                for (int &e : bestFoundRoutInPulse)
                    cout << e << " ";
                cout << endl;
            }
        }
    }

    for (int neighbor : nodeNeighbors[currentNode])
    {
        if (find(onePartialPath.begin(), onePartialPath.end(), neighbor) == onePartialPath.end())
        {
            if (PRINT_FOR_DEBUG_PULSE)
                cout << "start checking node " << currentNode << "'s neighbor " << neighbor << endl;
            pulse(neighbor, nodeNeighbors, partialObj, partialDist, dis, disLimit, onePartialPath, xCoeff, endNode);
        }
    }
}

#endif