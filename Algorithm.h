#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <fstream>
#include "readData.h"
#include <ctime>  //to measure CPU time
#include <chrono> //to measure run time
#include "Parameters.h"
#include "path.h"
#include <cstring>
#include <string>
#include <sstream>
#include <set>
#include <unordered_map>

using namespace std;
using namespace std::chrono;

using matrix = vector<vector<double>>;

class Algorithm
{
public:
    // Declare private attributes of Algorithm class here
    int numOfNode;
    int travelCost;
    int totalCapacity;
    int DIS;
    double vehicleWeight;
    double priceCharged;
    matrix w;
    matrix d;
    const readData &data;

public:
    // Constructor to copy values from readData class
    Algorithm(const readData &obj) : data(obj)
    {
        numOfNode = obj.numOfNode;
        travelCost = obj.travelCost;
        totalCapacity = obj.totalCapacity;
        DIS = obj.DIS;
        vehicleWeight = obj.vehicleWeight;
        priceCharged = obj.priceCharged;
        w = obj.w;
        d = obj.d;
    }

    Path GreedyAlgo1(fstream &csvFile, const string &file, bool distRatio_a = false, bool distRatio_b = false, bool locSrc = true)
    {
        clock_t CPUStart = clock();
        auto start = steady_clock::now();

        // add end node parameter <- num nodes
        int i = 1;               // current location
        int endNode = numOfNode; // find ending destination
        int numNodes = w.size();

        Path path(data);

        int nextPos = 0; // next node to travel to

        bool *nodesVisited = new bool[numNodes]{};
        nodesVisited[i] = 1;

        while (i != endNode)
        {                         // while the truck has not reached destination
            double maxWeight = 0; // max weight for current routes
            double maxRatio = 0;  // weights / distance ratio

            // iterate through each cargo request that start at the current location
            for (int j = 1; j < path.visited[i].size(); j++)
            {

                double cargoDist = d[i][j] + d[j][endNode]; // distance of cargo request + distance of request dropoff to end destination
                if (cargoDist == 0)
                {
                    double test = d[i][j];
                }
                // Tran added
                if (path.isFound(j))
                {
                    continue; // Skip this iteration if j is found in path
                }

                if (i == j)
                    continue;

                if ((!nodesVisited[j]) && (cargoDist + path.totalDist) <= DIS && w[i][j] <= totalCapacity)
                { // check if cargo request fits within capacity and distance restraints
                    if (!distRatio_a && !distRatio_b && maxWeight < w[i][j])
                    {                        // compare routes based on weight alone
                        maxWeight = w[i][j]; // weight of request
                        nextPos = j;         // destination of request
                    }
                    else if (distRatio_a && (w[i][j] / cargoDist) > maxRatio)
                    { // compare routes based on weight to distance ratio
                        maxRatio = w[i][j] / cargoDist;
                        maxWeight = w[i][j];
                        nextPos = j;
                    }
                    // Tran added
                    else if (distRatio_b && (w[i][j] / d[i][j]) > maxRatio)
                    { // compare routes based on weight to distance ratio
                        maxRatio = w[i][j] / d[i][j];
                        maxWeight = w[i][j];
                        nextPos = j;
                    }
                }
            }
            if (i != nextPos && nextPos != 0)
            {
                path.totalDist += d[i][nextPos];      // Add travel distance to cargo request desination
                path.visited[i][nextPos] = maxWeight; // add current capacity of truck (for profit calculation)
                nodesVisited[nextPos] = 1;
                // path.visited[nextPos][i] = 0;
                path.cargo[i][nextPos] = 1; // mark cargo request
            }

            if (maxWeight != 0) // go to destination  if no routes are found
                i = nextPos;
            else
            {
                path.totalDist += d[i][endNode];
                path.visited[i][endNode] = w[i][endNode];
                printf("\nvisited[%i][%i] -> %1.1f\n", i, endNode, w[i][endNode]);
                path.cargo[i][endNode] = 1;
                i = endNode;
            }
            path.push_back(i);
        }

        // Print results
        double greedyProfit = path.calculateProfit();

        string name = "";
        if (distRatio_a)
            name = "GreedyAlgo1_a_disRatio_a";
        else if (distRatio_b)
            name = "GreedyAlgo1_a_disRatio_b";
        else
            name = "GreedyAlgo1_a";

        // Do local search if toggled
        if (locSrc)
        {

            path.insertSearch_Tran_2();
            path.locSrcProfit = path.profit;
            path.replaceSearch();

            // localSearchProfit = path.fillCapacity(textFile);
            // localSearchProfit = path.fillCapacity_Daniel();
        }

        auto stop = steady_clock::now();
        clock_t CPUEnd = clock();
        double cpu_time_used = static_cast<double>(CPUEnd - CPUStart) * 1000.0 / CLOCKS_PER_SEC;
        double run_time = (duration<double>(stop - start).count()) * 1000.0;

        recordData(csvFile, file, name, greedyProfit, path, cpu_time_used, run_time);

        delete[] nodesVisited;
        return path;
    }

    /**
     * This is the copy of first greedy heuristic. The only difference is the function check the largest cargo excluding the last node.
     * @param distRatio if weight to distance ratio is used to compare routes
     * @return the path the truck will take as a vector
     */
    Path GreedyAlgo1_b(fstream &csvFile, const string &file, bool distRatio_a = false, bool distRatio_b = false, bool locSrc = true)
    {
        clock_t CPUStart = clock();
        auto start = steady_clock::now();

        int i = 1; // starting node / current location
        int endNode = numOfNode;
        int numNodes = w.size();

        Path path(data);

        // double currDist = 0;    // current dist traveled
        double totalWeight = 0; // total weight carried
        int nextPos = 0;        // next node

        while (i != endNode)
        {
            double maxWeight = 0; // max weight for current routes
            double maxRatio = 0;  // weights / distance

            for (int j = 1; j < path.visited[i].size() - 1; j++)
            {
                if (path.isFound(j))
                {
                    continue; // Skip this iteration if j is found in path
                }

                if (j != i)
                {
                    double cargoDist = d[i][j] + d[j][endNode];
                    if ((cargoDist + path.totalDist) <= DIS && w[i][j] <= totalCapacity)
                    {
                        // if the previous weight (@param tempWeight) is less than the current weight; update tempWeight is current weight.
                        if (!distRatio_a && !distRatio_b && maxWeight < w[i][j])
                        { // compare routes based on weight alone
                            maxWeight = w[i][j];
                            nextPos = j;
                        }
                        else if (distRatio_a && (w[i][j] / cargoDist) > maxRatio)
                        { // compare routes based on weight to distance ratio
                            maxRatio = w[i][j] / cargoDist;
                            maxWeight = w[i][j];
                            nextPos = j;
                        }
                        else if (distRatio_b && (w[i][j] / d[i][j]) > maxRatio)
                        { // compare routes based on weight to distance ratio
                            maxRatio = w[i][j] / d[i][j];
                            maxWeight = w[i][j];
                            nextPos = j;
                        }
                    }
                }
            }

            if (i != nextPos && nextPos != 0)
            {
                // Update cargo accepted to be 1
                path.totalDist += d[i][nextPos];      // Add travel distance
                path.visited[i][nextPos] = maxWeight; // @param tempWeight is the weight if the largest cargo weight found in the previous loop
                // nodesVisited[nextPos] = 1;
                path.visited[nextPos][i] = 0; // update the reverse cargo to be visited.
                path.cargo[i][nextPos] = 1;
            }

            if (maxWeight != 0) // If there is a route are found, update i = nextPos to start the new largest cargo weight search
                i = nextPos;
            else // go to end if no routes are found, then update i = endNode to finish the path
            {
                path.totalDist += d[i][endNode];
                path.visited[i][endNode] = w[i][endNode];
                printf("\nvisited[%i][%i] -> %1.1f\n", i, endNode, w[i][endNode]);
                path.cargo[i][endNode] = 1;
                maxWeight = w[i][endNode];
                i = endNode;
            }

            path.push_back(i);
            // totalWeight += maxWeight;
        }

        // Print results
        double greedyProfit = path.calculateProfit();

        string name = "";
        if (distRatio_a)
            name = "GreedyAlgo1_b_disRatio_a";
        else if (distRatio_b)
            name = "GreedyAlgo1_b_disRatio_b";
        else
            name = "GreedyAlgo1_b";

        if (locSrc)
        {

            path.insertSearch_Tran_2();
            path.locSrcProfit = path.profit;
            path.replaceSearch();
            // localSearchProfit = path.fillCapacity(textFile);
            // localSearchProfit = path.fillCapacity_Daniel();
        }

        auto stop = steady_clock::now();
        clock_t CPUEnd = clock();
        double cpu_time_used = static_cast<double>(CPUEnd - CPUStart) * 1000.0 / CLOCKS_PER_SEC;
        double run_time = (duration<double>(stop - start).count()) * 1000.0;

        recordData(csvFile, file, name, greedyProfit, path, cpu_time_used, run_time);

        return path;
    }

    /**
     * This is the second greedy heuristic that I mentioned. This algorithm looks through
     * all existing cargo requests and chooses the max load that is still within the limit of capacity and distance.
     * @param distRatio if weight to distance ratio is used to compare routes
     * @param locSrc if local search will be performed after greedy algorithm is found
     * @return the path the truck will take
     */
    Path GreedyAlgo2_b(fstream &csvFile, const string &file, bool distRatio = false, bool locSrc = true)
    {
        clock_t CPUStart = clock();
        auto start = steady_clock::now();

        int i = 1;               // current location
        int endNode = numOfNode; // destination
        int numNodes = w.size();

        Path path(data);

        vector<pair<int, double>> dropOffs; // list of saved cargo destinations and carresponding weights

        double capacity = (double)totalCapacity; // copy of capacity to be changed
        int *nodesVisited = new int[numNodes]{};

        double epsilon = std::numeric_limits<double>::epsilon();

        // double currDist = 0;    // current distance traveled
        double currWeight = 0; // total weight carried
        int nextPos = 0;       // next node to travel to
        double maxRatio = 0;

        printf("endNode=%d\n", endNode);
        printf("w size=%i\n", numNodes);
        printf("distRatio=%d\n", distRatio);
        w[1][1] = 0;

        while (i != endNode)
        {
            // Angela: add for testing
            printf("--------------- iteration with current location i=%d ---------------\n", i);
            // count++; if(count>2) exit(0);
            printf("==> available capacity= %.3f\n", capacity);

            double tempWeight = 0; // max weight for current routes
            // double maxRatio = 0;   // weight / distance
            int j_fin = 1, k_fin = 1;

            // iterate through every cargo request
            for (int j = 1; j < numNodes - 1; j++)
            {
                for (int k = 1; k < w[j].size() - 1; k++)
                {
                    if (k == j || k == i)
                        continue;

                    // ----calculate distance----
                    double totalDist = path.totalDist + d[i][j]; // distance to pickup
                    double totalWeight = currWeight + w[j][k];

                    if (!dropOffs.empty())
                    {
                        bool ifKInDropOffs = false; // if dropOff is already listed in dropOffs

                        totalDist += d[j][dropOffs[0].first]; // distance from pickup -> first dropoff
                        for (int c = 0; c < dropOffs.size() - 1; c++)
                        {
                            if (dropOffs[c].first == k)
                                ifKInDropOffs = true;
                            totalDist += d[dropOffs[c].first][dropOffs[c + 1].first]; // distance from dropoff -> dropoff
                        }

                        if (ifKInDropOffs)
                            totalDist += d[dropOffs[dropOffs.size() - 1].first][endNode]; // distance from last dropoff -> end
                        else
                            totalDist += d[dropOffs[dropOffs.size() - 1].first][k] + d[k][endNode]; // dist from last dropoff -> k -> endNode
                    }
                    else
                    {
                        totalDist += d[j][k] + d[k][endNode]; // dist from pickup -> dropoff -> endNode
                    }

                    if ((!nodesVisited[j] && nodesVisited[k] != 1 && path.cargo[j][k] == 0) && totalDist <= DIS && w[j][k] <= (capacity + epsilon))
                    { // if request fits within constraint of distance and capacity

                        if (!distRatio && w[j_fin][k_fin] < w[j][k])
                        { // compare routes using weight alone
                            j_fin = j;
                            k_fin = k;
                        }
                        else if (distRatio && (totalWeight / totalDist) > maxRatio)
                        { // compare routes using weight to distance ratio
                            cout << "max ratio: " << totalWeight << "/" << totalDist << endl;
                            maxRatio = totalWeight / totalDist;
                            j_fin = j;
                            k_fin = k;
                        }
                    }
                }
            }
            printf("==> after caculation, the selected cargo is w[%d][%d]= %.3f\n\n", j_fin, k_fin, w[j_fin][k_fin]);

            if (w[j_fin][k_fin] == 0)
            {
                if (!dropOffs.empty())
                { // no route found + cargo left to deliver = go to next dropoff location
                    path.visited[i][dropOffs[0].first] = (totalCapacity - capacity);
                    path.totalDist += d[i][dropOffs[0].first];

                    capacity += dropOffs[0].second;

                    nodesVisited[i] = 1;
                    i = dropOffs[0].first;
                    nodesVisited[i] = 0;
                    dropOffs.erase(dropOffs.begin());
                }
                else
                { // no route found + no cargo left to deliver = go to next delivery destination
                    path.visited[i][endNode] = (totalCapacity - capacity);
                    path.totalDist += d[i][endNode];
                    i = endNode; // go to destination
                }
            }
            else
            {                                                        // if a route is found
                path.visited[i][j_fin] = (totalCapacity - capacity); // add current capacity of truck for profit calculation]
                nodesVisited[i] = 1;
                nodesVisited[k_fin] = 2;
                path.totalDist += d[i][j_fin]; // distance to cargo pickup location
                capacity -= w[j_fin][k_fin];   // cargo load capacity

                // checks if the dropoff node is already listed
                auto nodeIter = find_if(dropOffs.begin(), dropOffs.end(), [k_fin](const pair<int, double> &pair)
                                        { return pair.first == k_fin; });
                // if it is already listed, just add to the weight to be dropped, otherwise add the new node
                if (nodeIter == dropOffs.end())
                    dropOffs.push_back(make_pair(k_fin, w[j_fin][k_fin])); // add delivery to list of deliveries
                else
                {
                    nodeIter->second += w[j_fin][k_fin];
                }

                currWeight += w[j_fin][k_fin]; // add new cargo to total weight
                path.cargo[j_fin][k_fin] = 1;  // update cargo list
                i = j_fin;                     // update current location to cargo request pickup location
            }

            path.push_back(i);
        }

        // Print results
        double greedyProfit = path.calculateProfit();

        string name = (distRatio) ? "GreedyAlgo2_b_disRatio" : "GreedyAlgo2_b";
        // recordText(textFile, file, name, path, 2);
        cout << "Final ratio: " << currWeight << "/" << path.totalDist << endl;
        cout << "   Cargo: ";
        for (int c = 1; c < path.cargo.size(); c++)
            for (int ca = 1; ca < path.cargo[c].size(); ca++)
                if (path.cargo[c][ca] > 0)
                    cout << c << ":" << ca << "=>" << w[c][ca] << " ";
        cout << endl;
        cout << "   Visited: ";
        for (int c = 1; c < path.visited.size(); c++)
            for (int ca = 1; ca < path.visited[c].size(); ca++)
                if (path.visited[c][ca] >= 0)
                    cout << c << ":" << ca << "=>(" << path.visited[c][ca] << ", " << d[c][ca] << ") ";
        cout << endl;

        if (locSrc)
        {

            path.insertSearch_Tran_2();
            path.locSrcProfit = path.profit;
            path.replaceSearch();
            // localSearchProfit = path.fillCapacity(textFile);
            // localSearchProfit = path.fillCapacity_Daniel();
        }

        auto stop = steady_clock::now();
        clock_t CPUEnd = clock();
        double cpu_time_used = static_cast<double>(CPUEnd - CPUStart) * 1000.0 / CLOCKS_PER_SEC;
        double run_time = (duration<double>(stop - start).count()) * 1000.0;

        recordData(csvFile, file, name, greedyProfit, path, cpu_time_used, run_time);

        return path;
    }

    /**
     * This is the second greedy heuristic that I mentioned. This algorithm looks through
     * all existing cargo requests and chooses the max load that is still within the limit of capacity and distance.
     * @param distRatio if weight to distance ratio is used to compare routes
     * @param locSrc if local search will be performed after greedy algorithm is found
     * @return the path the truck will take
     */
    Path GreedyAlgo2_c(fstream &csvFile, const string &file, bool locSrc = true)
    {
        clock_t CPUStart = clock();
        auto start = steady_clock::now();

        int i = 1;               // current location
        int endNode = numOfNode; // destination
        int numNodes = w.size();
        Path path(data);
        path.push_back(i);

        vector<pair<int, double>> dropOffs; // list of saved cargo destinations and carresponding weights

        double capacity = (double)totalCapacity; // copy of capacity to be changed
        int *nodesVisited = new int[numNodes]{};

        double epsilon = std::numeric_limits<double>::epsilon();

        double currDist = 0;    // current distance traveled
        double totalWeight = 0; // total weight carried
        int nextPos = 0;        // next node to travel to
        double currProfit = -1;

        double revenue = 0;
        double cost = 0;

        printf("endNode=%d\n", endNode);
        printf("w size=%i\n", numNodes);
        w[1][1] = 0;

        while (i != endNode)
        {
            // Angela: add for testing
            printf("--------------- iteration with current location i=%d ---------------\n", i);
            // count++; if(count>2) exit(0);
            printf("==> available capacity= %.3f\n", capacity);

            double tempWeight = 0; // max weight for current routes
            double maxRatio = 0;   // weight / distance
            int j_fin = 1, k_fin = 1;

            printf("==> currentProfit= %f\n", currProfit);

            // iterate through every cargo request
            for (int j = 1; j < numNodes - 1; j++)
            {
                for (int k = 1; k < w[j].size() - 1; k++)
                {
                    if (k == j || !(!nodesVisited[j] && nodesVisited[k] != 1) || !(w[j][k] <= (capacity + epsilon)))
                        continue;

                    // calculate profit if cargo is picked up
                    double routeDist = path.totalDist + d[i][j]; // distance to pickup
                    double totalRevenue = revenue + (w[j][k] * d[j][k]);
                    double tempCapacity = totalCapacity - capacity;
                    double totalCost = cost + (d[i][j] * (tempCapacity));
                    tempCapacity += w[j][k];

                    if (!dropOffs.empty())
                    {
                        bool ifKInDropOffs = false; // if dropOff is already listed in dropOffs

                        routeDist += d[j][dropOffs[0].first]; // distance from pickup -> first dropoff
                        totalCost += (d[j][dropOffs[0].first] * tempCapacity);
                        tempCapacity -= dropOffs[0].second;

                        for (int c = 0; c < dropOffs.size() - 1; c++)
                        {
                            if (dropOffs[c].first == k)
                                ifKInDropOffs = true;
                            routeDist += d[dropOffs[c].first][dropOffs[c + 1].first]; // distance from dropoff -> dropoff
                            totalCost += (d[dropOffs[c].first][dropOffs[c + 1].first] * tempCapacity);
                            tempCapacity -= dropOffs[c + 1].second;
                        }

                        if (ifKInDropOffs)
                        {
                            routeDist += d[dropOffs[dropOffs.size() - 1].first][endNode]; // distance from last dropoff -> end
                            totalCost += (d[dropOffs[dropOffs.size() - 1].first][endNode] * tempCapacity);
                        }
                        else
                        {
                            routeDist += d[dropOffs[dropOffs.size() - 1].first][k] + d[k][endNode]; // dist from last dropoff -> k -> endNode
                            totalCost += (d[dropOffs[dropOffs.size() - 1].first][k] * tempCapacity) + (d[k][endNode] * 0);
                        }
                    }
                    else
                    {
                        routeDist += d[j][k] + d[k][endNode]; // dist from pickup -> dropoff -> endNode
                        totalCost += (d[j][k] * tempCapacity) + (d[k][endNode] * 0);
                    }
                    double potentialProfit = (priceCharged * totalRevenue) - (travelCost * totalCost) - (travelCost * vehicleWeight * routeDist);

                    if (routeDist <= DIS && currProfit < potentialProfit)
                    { // compare routes using weight alone
                        j_fin = j;
                        k_fin = k;
                        currProfit = potentialProfit;
                        printf("====> potential profit for %i:%i= %f\n", j, k, potentialProfit);
                    }
                }
            }
            printf("==> after caculation, the selected cargo is w[%d][%d]= %.3f\n\n", j_fin, k_fin, w[j_fin][k_fin]);

            if (w[j_fin][k_fin] == 0)
            {

                if (!dropOffs.empty())
                { // no route found + cargo left to deliver = go to next dropoff location
                    path.visited[i][dropOffs[0].first] = (totalCapacity - capacity);
                    path.totalDist += d[i][dropOffs[0].first];
                    capacity += dropOffs[0].second;
                    cost += path.visited[i][dropOffs[0].first] * d[i][dropOffs[0].first];

                    nodesVisited[i] = 1;
                    i = dropOffs[0].first;
                    nodesVisited[i] = 0;
                    dropOffs.erase(dropOffs.begin());
                }
                else
                { // no route found + no cargo left to deliver = go to endNode
                    path.visited[i][endNode] = (totalCapacity - capacity);
                    path.totalDist += d[i][endNode];
                    i = endNode; // go to destination
                }
            }
            else
            {                                                        // if a route is found
                path.visited[i][j_fin] = (totalCapacity - capacity); // add current capacity of truck for profit calculation]
                nodesVisited[i] = 1;
                nodesVisited[k_fin] = 2;
                path.totalDist += d[i][j_fin]; // distance to cargo pickup location
                capacity -= w[j_fin][k_fin];   // cargo load capacity
                cost += path.visited[i][j_fin] * d[i][j_fin];

                // checks if the dropoff node is already listed
                auto nodeIter = find_if(dropOffs.begin(), dropOffs.end(), [k_fin](const pair<int, double> &pair)
                                        { return pair.first == k_fin; });
                // if it is already listed, just add to the weight to be dropped, otherwise add the new node
                if (nodeIter == dropOffs.end())
                    dropOffs.push_back(make_pair(k_fin, w[j_fin][k_fin])); // add delivery to list of deliveries
                else
                {
                    nodeIter->second += w[j_fin][k_fin];
                }

                totalWeight += w[j_fin][k_fin]; // add new cargo to total weight
                revenue += w[j_fin][k_fin] * d[j_fin][k_fin];
                path.cargo[j_fin][k_fin] = 1; // update cargo list
                i = j_fin;                    // update current location to cargo request pickup location
            }

            path.push_back(i);
        }

        // Print results
        double greedyProfit = path.calculateProfit();

        // path.locSrcProfit = profit;
        string name = "GreedyAlgo2_c";

        // perform local search if toggled
        if (locSrc)
        {
            path.insertSearch_Tran_2();
            path.locSrcProfit = path.profit;
            path.replaceSearch();
            // localSearchProfit = path.fillCapacity_Daniel();
        }

        auto stop = steady_clock::now();
        clock_t CPUEnd = clock();
        double cpu_time_used = static_cast<double>(CPUEnd - CPUStart) * 1000.0 / CLOCKS_PER_SEC;
        double run_time = (duration<double>(stop - start).count()) * 1000.0;

        recordData(csvFile, file, name, greedyProfit, path, cpu_time_used, run_time);
        return path;
    }
    /**
     * @brief Greedy Algorithm 3 is searching for the routes that generate largest profit excluding the last node.
     * In the GreadyAlgo3_a, the profit is calculated from the pickup to the nextNode itself
     * in the GreadyAlgo3_b, the profit is calculated from the pickup to the nextNode, nextNode to endNode.
     *
     * @param csvFile
     * @param textFile
     * @param file
     * @param a
     * @param b
     * @param locSrc
     * @return path object
     */
    Path GreedyAlgo3(fstream &csvFile, const string &file, bool a = false, bool b = false, bool locSrc = true)
    {

        clock_t CPUStart = clock();
        auto start = steady_clock::now();

        int i = 1;               // current location
        int endNode = numOfNode; // find ending destination

        Path path(data);

        int nextPos = 0; // next node to travel to

        while (i != endNode)
        {
            double maxProfit = 0;
            double curWeight = 0;

            for (int j = 1; j < numOfNode; j++)
            {

                if (path.isFound(j))
                {
                    continue; // Skip this iteration if j is found in path
                }

                if (j != i)
                {
                    double cargoDist = d[i][j] + d[j][endNode];
                    if (path.visited[i][j] == -1 && (cargoDist + path.totalDist) <= DIS && w[i][j] <= totalCapacity)
                    {

                        double tempProfit;

                        if (a)
                        {
                            tempProfit = (priceCharged * d[i][j] * w[i][j]) - (travelCost * d[i][j] * w[i][j]) - (travelCost * vehicleWeight * d[i][j]);
                        }
                        else if (b)
                        {
                            tempProfit = (priceCharged * (d[i][j] * w[i][j] + d[j][endNode] * w[j][endNode])) - (travelCost * (d[i][j] * w[i][j] + d[j][endNode] * w[j][endNode])) - (travelCost * vehicleWeight * (d[i][j] + d[j][endNode]));
                        }
                        if (tempProfit > maxProfit && tempProfit > 0)
                        {
                            maxProfit = tempProfit;
                            nextPos = j;
                            curWeight = w[i][j];
                        }
                    }
                }
            }

            if (i != nextPos && nextPos != 0)
            {
                // Update cargo accepted to be 1
                path.totalDist += d[i][nextPos];      // Add travel distance
                path.visited[i][nextPos] = curWeight; // @param tempWeight is the weight if the largest cargo weight found in the previous loop
                path.visited[nextPos][i] = 0;         // update the reverse cargo to be visited.
                path.cargo[i][nextPos] = 1;
            }

            if (maxProfit != 0) // If there is a route are found, update i = nextPos to start the new largest cargo weight search
                i = nextPos;
            else // go to end if no routes are found, then update i = endNode to finish the path
            {
                path.totalDist += d[i][endNode];
                path.visited[i][endNode] = w[i][endNode];
                path.visited[endNode][i] = 0; // update the reverse cargo to be visited.
                path.cargo[i][endNode] = 1;
                curWeight = w[i][endNode];
                i = endNode;
            }

            path.push_back(i);
        }

        // Get profit
        double greedyProfit = path.calculateProfit();

        string name = "";
        if (a)
            name = "GreedyAlgo3_a";
        else if (b)
            name = "GreedyAlgo3_b";

        if (locSrc)
        {

            path.insertSearch_Tran_2();
            path.locSrcProfit = path.profit;
            path.replaceSearch();
            // localSearchProfit = path.fillCapacity_Daniel();
        }

        auto stop = steady_clock::now();
        clock_t CPUEnd = clock();
        double cpu_time_used = static_cast<double>(CPUEnd - CPUStart) * 1000.0 / CLOCKS_PER_SEC;
        double run_time = (duration<double>(stop - start).count()) * 1000.0;

        recordData(csvFile, file, name, greedyProfit, path, cpu_time_used, run_time);

        return path;
    }

    Path GreedyAlgo4(fstream &csvFile, const string &file, bool locSrc = true)
    {

        clock_t CPUStart = clock();
        auto start = steady_clock::now();

        Path path(data);

        int endNode = numOfNode;
        int i = 1;

        int nextPos = 0;

        while (i != endNode)
        {

            bool noNode = true;
            double tempDist = DIS;
            for (int j = 1; j < numOfNode - 1; j++)
            {

                if (path.isFound(j))
                {
                    continue; // Skip this iteration if j is found in path
                }

                if (j != i)
                {

                    double cargoDist = d[i][j] + d[j][endNode];
                    if (path.visited[i][j] == -1 && (cargoDist + path.totalDist) <= DIS && w[i][j] <= totalCapacity)
                    {
                        if (d[i][j] < tempDist)
                        {
                            tempDist = d[i][j];
                            nextPos = j;
                            noNode = false;
                        }
                    }
                }
            }

            if (i != nextPos && nextPos != 0)
            {
                // Update cargo accepted to be 1
                path.totalDist += d[i][nextPos];          // Add travel distance
                path.visited[i][nextPos] = w[i][nextPos]; // @param tempWeight is the weight if the largest cargo weight found in the previous loop        // update the reverse cargo to be visited.
                path.cargo[i][nextPos] = 1;
            }

            if (!noNode) // If there is a route are found, update i = nextPos to start the new largest cargo weight search
                i = nextPos;
            else // go to end if no routes are found, then update i = endNode to finish the path
            {
                path.totalDist += d[i][endNode];
                path.visited[i][endNode] = w[i][endNode];
                path.cargo[i][endNode] = 1;
                i = endNode;
            }

            path.push_back(i);
        }

        double greedyProfit = path.calculateProfit();

        if (locSrc)
        {

            path.insertSearch_Tran_2();
            path.locSrcProfit = path.profit;
            path.replaceSearch();
            //        //localSearchProfit = path.insertSearch_Daniel(textFile);
            // localSearchProfit = path.fillCapacity_Daniel();
        }

        auto stop = steady_clock::now();
        clock_t CPUEnd = clock();
        double cpu_time_used = static_cast<double>(CPUEnd - CPUStart) * 1000.0 / CLOCKS_PER_SEC;
        double run_time = (duration<double>(stop - start).count()) * 1000.0;

        recordData(csvFile, file, "GreedyAlgo4", greedyProfit, path, cpu_time_used, run_time);

        return path;
    }

private:
    static void printPathInfo(double &totalWeight, double &currDist, bool distRatio, vector<int> path, vector<vector<int>> cargo, double greedyNum)
    {
        cout << "Weight: " << totalWeight << endl;
        cout << "Distance: " << currDist << endl;
        cout << "Ratio: " << totalWeight / currDist << endl;
        cout << "Path for Greedy #" << greedyNum << " (distRatio = " << distRatio << "): ";
        for (auto &c : path)
            cout << c << " ";
        cout << endl;
        cout << "Cargo Accepted: ";
        for (size_t i = 0; i < cargo.size(); ++i)
        {
            for (size_t j = 0; j < cargo[i].size(); ++j)
            {
                if (cargo[i][j] == 1)
                {
                    cout << i << ":" << j << " ";
                }
            }
        }
        cout << endl
             << endl;
    }

    /**
     * @brief this function write oupute to csv file
     *
     */
    void recordData(fstream &outputFile, const string &file, string algoType, double greedyProfit, const Path &path, double cpuTime, double realTime) //, vector<int> requested)
    {
        cout << file << "," << algoType << "," << greedyProfit << "," << path.locSrcProfit << "," << cpuTime << "," << realTime << ",";
        outputFile << file << "," << algoType << "," << greedyProfit << "," << path.locSrcProfit << "," << cpuTime << "," << realTime << ",";

        cout << path;
        outputFile << path;

        cout << ",";
        outputFile << ",";

        path.printCargo(outputFile);
        path.printCargo(cout);
        // output
    }

    void recordText(fstream &textFile, const string &file, string algoType, const Path &path, double greedyNum)
    {
        textFile << "Filename: " << file << endl;
        textFile << "Algorithm: " << algoType << endl;
        textFile << "Profit: " << path.profit << endl;
        textFile << "Weight: " << path.totalWeight() << endl;
        textFile << "Distance: " << path.totalDist << endl;
        textFile << "Path for Greedy #" << greedyNum << ": ";
        textFile << path;
        textFile << endl;
        textFile << "Cargo Accepted: ";
        path.printCargo(textFile);
        textFile << endl
                 << endl;
    }

    std::vector<std::string> split(std::string s, std::string delimiter)
    {
        size_t pos_start = 0, pos_end, delim_len = delimiter.length();
        std::string token;
        std::vector<std::string> res;

        while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos)
        {
            token = s.substr(pos_start, pos_end - pos_start);
            pos_start = pos_end + delim_len;
            res.push_back(token);
        }

        res.push_back(s.substr(pos_start));
        return res;
    }
};

#endif