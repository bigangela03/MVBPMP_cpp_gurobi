#ifndef DOMINANCE_BIT_HASHMAP_H
#define DOMINANCE_BIT_HASHMAP_H

//!!!!!!!!!!!!!!!!!!! NOT WORKING !!!!!!!!!!!!!!!!!!!!!!!!!!
// the data structre of hashmap is not working now!!!!!!!!!!!

// This code is only for complete graph right now, i.e., any two nodes are connected
// since I didn't consider the vertex incidence matrix, only iterate nodes from 0 to n-1

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

void printLabelSet(&unordered_set<long>, &unordered_map<long, vector<int>>, &unordered_map<long, double>, &unordered_map<long, double>, int);
void printOneLabel(long, &unordered_map<long, vector<int>>, &unordered_map<long, double>, &unordered_map<long, double>, int);
bool compareToLabel(vector<double>, vector<double>, int, double);
int checkIfLabelIsDominated(vector<double>, vector<double>, int);

void runDominance(int n, double **dis, vector<vector<double>> xCoeff, double disLimit, double *obj, vector<int> &nodesOnRoute)
{

    //===> set up parameters used in dominance.h
    int startingNode = 0;
    int endingNode = n - 1;
    bool PRINT_FOR_DEBUG = false;
    int bigM = 10000000; // bigM has to be big enough, so than when checking domincated labels inside F set, when comparing two equivalent labels, it helps to find label with more bigM
    double verySmallNum = 0.0000001;

    if (startingNode < 0 || startingNode > n - 1)
    {
        printf("ERROR: startingNode is supposed to be {0,...,n-1}\n");
        exit(1);
    }
    if (endingNode < 0 || endingNode > n - 1)
    {
        printf("ERROR: endingNode is supposed to be {0,...,n-1}\n");
        exit(1);
    }

    long *twoPow = new long[n];

    // for (int i = 0; i < n; i++)
    //     twoPow[i] = 1 << i;
    for (int i = 0; i < n; i++)
    {
        twoPow[i] = pow(2, i);
    }

    if (PRINT_FOR_DEBUG)
        cout << "===> start dominace()" << endl;

    if (PRINT_FOR_DEBUG)
    {
        cout << "cost of each arc" << endl;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                cout << xCoeff[i][j] << " ";
            cout << endl;
        }
    }

    vector<unordered_set<long>> bigLambda_allNodesLableSet;
    unordered_map<long, vector<int>> routeMap;
    unordered_map<long, double> distMap;
    unordered_map<long, double> costMap;
    unordered_map<long, double> numVisitedNodesMap;

    //===> initialize the lable for starting node
    // vector<double> oneLabel;
    // for (int i = 0; i < n + 3; i++)
    //     oneLabel.push_back(0);
    // oneLabel[startingNode] = 1;
    // oneLabel[n] = 1;
    // set<vector<double>> startingNodeLables;
    // startingNodeLables.insert(oneLabel);

    long B2Dtemp = twoPow[n - startingNode - 1]; // B2D, binary to decimal
    distMap[B2Dtemp] = 0;
    costMap[B2Dtemp] = 0;
    numVisitedNodesMap[B2Dtemp] = 1;
    vector<int> routeTemp;
    for (int i = 0; i < n; i++)
        routeTemp.push_back(0);
    routeTemp[startingNode] = 1;
    routeMap[B2Dtemp] = routeTemp;
    unordered_set<long> startingNodeLables;
    startingNodeLables.insert(B2Dtemp);

    //===> add empty lable set to all other nodes
    // for (int i = 0; i < n; i++)
    // {
    //     if (i == startingNode)
    //         bigLambda_allNodesLableSet.push_back(startingNodeLables);
    //     else
    //     {
    //         set<vector<double>> emptySet;
    //         bigLambda_allNodesLableSet.push_back(emptySet);
    //     }
    // }

    for (int i = 0; i < n; i++)
    {
        if (i == startingNode)
            bigLambda_allNodesLableSet.push_back(startingNodeLables);
        else
        {
            unordered_set<long> emptySet;
            bigLambda_allNodesLableSet.push_back(emptySet);
        }
    }

    //===> initialize the newly added labels set as well
    // vector<set<vector<double>>> newLabelSet(bigLambda_allNodesLableSet);
    vector<unordered_set<long>> newLabelSet(bigLambda_allNodesLableSet);

    unordered_set<int> E_toBeTreatedNodesSet;
    E_toBeTreatedNodesSet.insert(startingNode);

    int itrNum = 0; // iteration number

    while (E_toBeTreatedNodesSet.size() > 0) //===> Line 7 in ESPPRC(p)
    {
        itrNum++;

        if (PRINT_FOR_DEBUG)
        {
            cout << endl;
            cout << "======> dominace iteration " << itrNum << " <======" << endl;
            cout << "E_toBeTreatedNodesSet size = " << E_toBeTreatedNodesSet.size() << endl;
            for (auto &temp : E_toBeTreatedNodesSet)
                cout << temp << " ";
            cout << endl;
        }

        //===> only check the labels newly added in the last iteration
        vector<unordered_set<long>> labelsAddedInLastIteration(newLabelSet);
        for (int i = 0; i < newLabelSet.size(); i++)
            newLabelSet[i].clear();

        unordered_set<int> E_toBeTreatedNodesSetCopy;

        for (auto &node : E_toBeTreatedNodesSet)
            E_toBeTreatedNodesSetCopy.insert(node);

        // clear the original node set, and add the node with change later
        E_toBeTreatedNodesSet.clear();

        for (auto &node : E_toBeTreatedNodesSetCopy) //===> LINE 9 in ESPPRC(p)
        {

            if (labelsAddedInLastIteration[node].size() < 1)
            {
                printf("ERROR: There must be labels in E set!");
                exit(1);
            }

            // ===> ending node (n-1) has no successors
            if (node == n - 1)
                continue;

            if (PRINT_FOR_DEBUG)
                cout << "\n===> treaing node " << node << " in iteration " << itrNum << endl;

            for (int succ = 0; succ < n; succ++) //===> LINE 10 in ESPPRC(p)
            {
                // ===> there are arcs only out of startingNode
                if (succ == startingNode || succ == node)
                    continue;

                if (PRINT_FOR_DEBUG)
                {
                    cout << "\nsucc=" << succ << " (itr " << itrNum << ", node " << node << ")" << endl;
                    cout << "===> bigLambda_allNodesLableSet[" << node << "].size =" << bigLambda_allNodesLableSet[node].size() << endl;
                    printLabelSet(bigLambda_allNodesLableSet[node], distMap, costMap, bigM);
                    cout << "===> labelsAddedInLastIteration[" << node << "].size =" << labelsAddedInLastIteration[node].size() << endl;
                    printLabelSet(labelsAddedInLastIteration[node], distMap, costMap, bigM);
                    cout << "===> newLabelSet[" << node << "].size =" << newLabelSet[node].size() << endl;
                    printLabelSet(newLabelSet[node], distMap, costMap, bigM);
                }

                set<long> F_nodeSuccLabelSet; //===> LINE 11 in ESPPRC(p) PSEUDO-CODE
                                              // set<vector<double>> labelsOverDistLimit;

                for (auto &label : labelsAddedInLastIteration[node]) //===> LINE 12 in ESPPRC(p) PSEUDO-CODE
                {

                    // if (label[succ] == 0) //===> LINE 13 in ESPPRC(p) PSEUDO-CODE
                    if (label | twoPow[succ] == 0) //===> LINE 13 in ESPPRC(p) PSEUDO-CODE
                    {
                        if (PRINT_FOR_DEBUG)
                        {
                            cout << "find a label not visit succ yet" << endl;
                            printOneLabel(label, distMap, costMap, bigM);
                        }

                        //===> LINE 14 in ESPPRC(p) PSEUDO-CODE

                        //===> check resource (distance)
                        double currentDistance = distMap[label];

                        if (currentDistance + dis[node][succ] + dis[succ][endingNode] <= disLimit)
                        {
                            long labelTemp = label;

                            int numVisitedNodes = numVisitedNodesMap[label];
                            work here !!!!!!!!!!!!!!labelTemp[succ] = numVisitedNodes + 1;
                            labelTemp += twoPow[n - succ - 1];
                            labelTemp[n] = numVisitedNodes + 1;                   //!!! lable[n] is the number of visited nodes;
                            labelTemp[n + 1] = currentDistance + dis[node][succ]; //!!! label[n+1] is the consumed resouces;
                            labelTemp[n + 2] = label[n + 2] + xCoeff[node][succ]; //!!! lable[n+2] is the cost
                            F_nodeSuccLabelSet.insert(labelTemp);

                            if (PRINT_FOR_DEBUG)
                            {
                                cout << "visited number of nodes: labelTemp[n]=" << labelTemp[n] << endl;
                                cout << "traveled distance: labelTemp[n + 1]=" << labelTemp[n + 1] << endl;
                            }
                        }
                        else
                        {
                            if (PRINT_FOR_DEBUG)
                                cout << "distance is over limit" << endl;
                            // labelsOverDistLimit.insert(label);
                        }
                    }
                }

                // if (labelsOverDistLimit.size() > 0)
                // {
                //     cout << "the label before replacing to M" << endl;
                //     printLabelSet(bigLambda_allNodesLableSet[node], bigM);
                // }

                // for (auto &labelTemp : labelsOverDistLimit)
                // {
                //     auto it = bigLambda_allNodesLableSet[node].find(labelTemp);

                //     if (it != bigLambda_allNodesLableSet[node].end())
                //     {
                //         vector<double> temp = *it;
                //         // temp[succ] = 1;
                //         temp[succ] = bigM;
                //         bigLambda_allNodesLableSet[node].erase(it);
                //         bigLambda_allNodesLableSet[node].insert(temp);

                //         if (PRINT_FOR_DEBUG)
                //         {
                //             printf("lable replaced: ");
                //             printOneLabel(labelTemp, bigM);
                //             printf("new label     : ");
                //             printOneLabel(temp, bigM);
                //         }

                //         if (compareToLabel(temp, targetLabel, n, verySmallNum))
                //         {
                //             printf("==> target label is found in node %d for succ %d in iteration %d\n", node, succ, itrNum);
                //         }
                //     }
                //     else
                //     {
                //         cout << "Error: labels that exceed distance limit is not found." << endl;
                //         exit(1);
                //     }
                // }
                // if (labelsOverDistLimit.size() > 0)
                // {
                //     cout << "the label after replacing to M" << endl;
                //     printLabelSet(bigLambda_allNodesLableSet[node], bigM);
                // }

                //========> EFF function in LINE 15 and 16 in ESPPRC(p) PSEUDO-CODE <========
                {
                    //===> check if each label in F_nodeSuccLabelSet is dominated by bigLambda_allNodesLableSet[succ]
                    // if each label is dominated, then no need to add to bigLambda_allNodesLableSet[succ], so no change, and hasChange = false;
                    // if not, then the non-dominated label is added to bigLambda_allNodesLableSet[succ], so hasChange = true;

                    //===> check the labels of visiting nodes
                    // EFF() step (1): compare labels inside new-label-set
                    // EFF() step (2): remove dominated lables inside new-label-set
                    // EFF() step (3): compare labels of old and updated new-label-set without dominated ones
                    // EFF() step (4): remove dominated labels from old and updated new-label-set

                    //===> step (1)
                    if (PRINT_FOR_DEBUG)
                    {
                        cout << "# of labels in F set: " << F_nodeSuccLabelSet.size() << endl;
                        cout << "-> F set before removing dominated labels" << endl;
                        printLabelSet(F_nodeSuccLabelSet, bigM);
                    }

                    vector<vector<double>> Fset;
                    for (auto &labelTemp : F_nodeSuccLabelSet)
                    {
                        vector<double> labelToAdd;
                        for (auto &ele : labelTemp)
                            labelToAdd.push_back(ele);
                        Fset.push_back(labelToAdd);
                    }

                    set<vector<double>> dominatedLabelsAmongNewLabelSet;

                    if (Fset.size() >= 2)
                    { // since we are checking lables in a Set, so no labels are duplicate according to Set definition
                        for (int i = 0; i < Fset.size(); i++)
                            for (int j = i + 1; j < Fset.size(); j++)
                            {

                                vector<double> label1 = Fset[i];
                                vector<double> label2 = Fset[j];

                                // bool label1Dominated = true;
                                // bool label2Dominated = true;

                                // if (!(label1[n + 1] <= label2[n + 1] && label1[n + 2] <= label2[n + 2]))
                                // {
                                //     label2Dominated = false;
                                //     if (!(label1[n + 1] >= label2[n + 1] && label1[n + 2] >= label2[n + 2])) // new label is not dominated
                                //         continue;
                                // }
                                // else // old label is dominated
                                // {
                                //     // new label is not dominated
                                //     if (!(label1[n + 1] >= label2[n + 1] && label1[n + 2] >= label2[n + 2]))
                                //         label1Dominated = false;
                                // }
                                // // if both new label and old label are dominated, they are the same

                                // // keep checking other elements in the label
                                // if (label1Dominated)
                                //     if (label1[n] < label2[n]) // check the number of visted nodes first
                                //         continue;
                                //     else // then check each node
                                //         for (int i = 0; i < n; i++)
                                //             if (label1[i] < label2[i])
                                //             {
                                //                 label1Dominated = false;
                                //                 break;
                                //             }

                                // if (label2Dominated)
                                //     if (label2[n] < label1[n]) // check the number of visted nodes first
                                //         continue;
                                //     else // then check each node
                                //         for (int i = 0; i < n; i++)
                                //             if (label2[i] < label1[i])
                                //             {
                                //                 label2Dominated = false;
                                //                 break;
                                //             }

                                // // we are checking labels in a set, so this scenario doesn't happen
                                // if (label1Dominated && label2Dominated) // two lables are the same
                                //     dominatedLabelsAmongNewLabelSet.insert(label1);

                                // if (label1Dominated && !label2Dominated)
                                //     dominatedLabelsAmongNewLabelSet.insert(label1);

                                // if (!label1Dominated && label2Dominated)
                                //     dominatedLabelsAmongNewLabelSet.insert(label2);

                                int dominationResult = checkIfLabelIsDominated(label1, label2, n);
                                // 1. both are dominated, which means they are the same 2. only newLabel is dominated; 3. only oldLabel is dominated;
                                if (dominationResult == 1)
                                {
                                    printf("ERROR: two lables can't be the same in a set!");
                                    exit(1);
                                }
                                else if (dominationResult == 2)
                                    dominatedLabelsAmongNewLabelSet.insert(label1);
                                else if (dominationResult == 3)
                                    dominatedLabelsAmongNewLabelSet.insert(label2);
                            }
                    }

                    //===> step (2)
                    //===> remove the dominated labels from F_nodeSuccLabelSet
                    for (auto &labelTemp : dominatedLabelsAmongNewLabelSet)
                        F_nodeSuccLabelSet.erase(labelTemp);

                    if (PRINT_FOR_DEBUG)
                    {
                        if (dominatedLabelsAmongNewLabelSet.size() > 0)
                        {
                            cout << "-> dominatedLabelsAmongNewLabelSet" << endl;

                            printLabelSet(dominatedLabelsAmongNewLabelSet, bigM);

                            cout << "# of labels in F set after update: " << F_nodeSuccLabelSet.size() << endl;
                            cout << "-> F set after removing dominated labels" << endl;
                            printLabelSet(F_nodeSuccLabelSet, bigM);
                        }
                        else
                            cout << "no dominated lables found among F set" << endl;
                    }

                    //===> step (3)

                    set<vector<double>> dominatedNewLabelSet;
                    set<vector<double>> dominatedExistingLabelSet;

                    for (auto &newLabel : F_nodeSuccLabelSet)
                    {

                        // we can compare dist and cost first, if  not (d1>=d2, c1>=c2, or d1<=d2, c1<=c2), continue

                        // then, use some other methods to decide, like if any element violates dominated conditin, then break, so we don't need to check all the elements in a label.

                        for (auto &oldLabel : bigLambda_allNodesLableSet[succ])
                        {
                            int dominationResult = checkIfLabelIsDominated(newLabel, oldLabel, n);
                            // 1. both are dominated, which means they are the same 2. only newLabel is dominated; 3. only oldLabel is dominated;
                            if (dominationResult == 1)
                                dominatedNewLabelSet.insert(newLabel);
                            else if (dominationResult == 2)
                                dominatedNewLabelSet.insert(newLabel);
                            else if (dominationResult == 3)
                                dominatedExistingLabelSet.insert(oldLabel);

                            if (PRINT_FOR_DEBUG)
                            {
                                cout << "new label: ";
                                printOneLabel(newLabel, bigM);
                                cout << endl;
                                cout << "old label: ";
                                printOneLabel(oldLabel, bigM);
                                cout << endl;
                                if (dominationResult == 1)
                                    cout << "Two labels are the same." << endl;
                                else if (dominationResult == 2)
                                    cout << "New label is dominated by old label." << endl;
                                else if (dominationResult == 3)
                                    cout << "Old label is dominated." << endl;
                            }
                        }
                    }
                    if (PRINT_FOR_DEBUG)
                    {
                        cout << "# of labels at succ node: " << bigLambda_allNodesLableSet[succ].size() << endl;
                        cout << "-> labels before update succ=" << succ << " labels" << endl;
                        printLabelSet(bigLambda_allNodesLableSet[succ], bigM);

                        cout << "-> dominatedExistingLabelSet" << endl;
                        printLabelSet(dominatedExistingLabelSet, bigM);

                        cout << "-> dominatedNewLabelSet" << endl;
                        printLabelSet(dominatedNewLabelSet, bigM);

                        cout << "-> F set before update" << endl;
                        printLabelSet(F_nodeSuccLabelSet, bigM);
                    }

                    //===> step (4)
                    //===> remove the dominated labels for current and succsor nodes
                    for (auto &labelTemp : dominatedExistingLabelSet)
                        bigLambda_allNodesLableSet[succ].erase(labelTemp);
                    for (auto &labelTemp : dominatedNewLabelSet)
                        F_nodeSuccLabelSet.erase(labelTemp);

                    //===> add the non-dominated labels to bigLambda_allNodesLableSet[succ]
                    for (auto &labelTemp : F_nodeSuccLabelSet)
                    {
                        bigLambda_allNodesLableSet[succ].insert(labelTemp);
                        newLabelSet[succ].insert(labelTemp);
                    }

                    if (PRINT_FOR_DEBUG)
                    {
                        // cout << "======> ERROR: a conflict found!" << endl;
                        cout << "-> F_nodeSuccLabelSet" << endl;
                        printLabelSet(F_nodeSuccLabelSet, bigM);

                        cout << "# of labels at succ node: " << bigLambda_allNodesLableSet[succ].size() << endl;
                        cout << "-> labels after update succ=" << succ << " labels" << endl;
                        printLabelSet(bigLambda_allNodesLableSet[succ], bigM);
                    }

                    // {
                    //     for (auto &labTemp : F_nodeSuccLabelSet)
                    //     {
                    //         if (compareToLabel(labTemp, targetLabel, n, verySmallNum))
                    //         {
                    //             printf("==> target label is found in F_(%d,%d) in iteration %d\n", node, succ, itrNum);
                    //         }
                    //     }
                    // }

                    if (F_nodeSuccLabelSet.size() != 0) //===> LINE 16 in ESPPRC(p) PSEUDO-CODE//===> LINE 16 in ESPPRC(p) PSEUDO-CODE
                    {
                        E_toBeTreatedNodesSet.insert(succ); //===> LINE 17 in ESPPRC(p) PSEUDO-CODE
                        if (PRINT_FOR_DEBUG)
                            printf("succ=%d is added int set E\n", succ);
                    }
                }
            }
        }

        // since we are looking for the elementary shortest path, so for n nodes, at most n-1 arc can be added to path
        // since each iteration can add one arc, so n-1 iterations is enough

        if (itrNum >= n)
        {
            // break;
            printf("ERROR: in dominance.h, the while loop runs over n iterations, which is not supposed to happen.");
            exit(1);
        }
    }

    if (PRINT_FOR_DEBUG)
    {
        printf("\n-> the found visited nodes in label\n");
        printLabelSet(bigLambda_allNodesLableSet[n - 1], bigM);
    }

    double minCost = numeric_limits<double>::max();
    vector<double> bestLabel;

    for (auto &labelTemp : bigLambda_allNodesLableSet[n - 1])
    {
        if (labelTemp[n + 2] < minCost)
        {
            bestLabel.clear();
            minCost = labelTemp[n + 2];

            for (auto &ele : labelTemp)
                bestLabel.push_back(ele);
        }
    }

    if (PRINT_FOR_DEBUG)
    {
        cout << "best label:" << endl;
        for (auto &ele : bestLabel)
            cout << ele << " ";
        cout << endl;
    }

    *obj = minCost;

    for (int i = 0; i < n; i++)
        nodesOnRoute.push_back(0);

    //===> for example: label=[1, bigM, 2, 3, 3, 12.8, -1.28], n=4
    //===> it means, node 0 is visted first, then node 2, then node 3: 0->2->3

    double sum = 0;
    double numAllVistedNodes;
    for (int i = 0; i < n; i++)
    {
        if (bestLabel[i] < bigM - verySmallNum && bestLabel[i] > verySmallNum)
        {
            int index = (long)bestLabel[i];
            nodesOnRoute[index - 1] = i;
            numAllVistedNodes = bestLabel[n];
            sum += bestLabel[i];
        }
    }
    if (sum != numAllVistedNodes * (numAllVistedNodes + 1) / 2)
    {
        printf("ERROR: the number in first n elements of a label in dominace method should be in set{1, 2, .., numAllVistedNodes}\n");
        printf("the visit order jumped, for example, label=[1, bigM, 2, 4, #, dis, cost] (n=4)\n");
        printf("the right result is supposed to be label=[1,bigM,2,3, #,dis,cost]\n");
        exit(1);
    }

    if (PRINT_FOR_DEBUG)
    {
        cout << "selected route before removing 0s:" << endl;
        for (auto &ele : nodesOnRoute)
            cout << ele << " ";
        cout << endl;

        cout << "# of visited nodes: " << numAllVistedNodes << endl;
    }

    nodesOnRoute.erase(nodesOnRoute.begin() + (int)numAllVistedNodes, nodesOnRoute.begin() + n);

    if (PRINT_FOR_DEBUG)
    {
        cout << "selected route after removing 0s:" << endl;
        for (auto &ele : nodesOnRoute)
            cout << ele << " ";
        cout << endl;
    }
}

void printLabelSet(set<long> labelSet, &unordered_map<long, vector<int>> routeMap, &unordered_map<long, double> distMap, &unordered_map<long, double> costMap, double bigM)
{
    for (auto &labelTemp : labelSet)
    {
        for (auto &temp : routeMap[labelTemp])
            // cout << fixed << setprecision(2) << temp << " ";
            if (temp == bigM)
                cout << "M ";
            else
                cout << temp << " ";
        cout << cout << labelTemp << " ";
        cout << distMap[labelTemp] << " ";
        cout << costMap[labelTemp] << " ";
        cout << endl;
    }
}

work here !!!!!!!!!!!!!!!!!!void printOneLabel(long label, &unordered_map<long, double> distMap, &unordered_map<long, double> costMap, double bigM)
{
    // for (auto &temp : label)
    // cout << fixed << setprecision(2) << temp << " ";
    // if (temp == bigM)
    //     cout << "M ";
    // else
    cout << label << " ";
    cout << distMap[label] << " ";
    cout << costMap[label] << " ";
    cout << endl;
}

int checkIfLabelIsDominated(vector<double> newLabel, vector<double> oldLabel, int n)
{

    bool newLabelDominated = true;
    bool oldLabelDominated = true;

    int returnValue = 0;

    if (!(newLabel[n + 1] <= oldLabel[n + 1] && newLabel[n + 2] <= oldLabel[n + 2])) // old label is not dominated
    {
        oldLabelDominated = false;
        if (!(newLabel[n + 1] >= oldLabel[n + 1] && newLabel[n + 2] >= oldLabel[n + 2])) // new label is not dominated
            return returnValue;
    }
    else // old label is dominated
    {
        // new label is not dominated
        if (!(newLabel[n + 1] >= oldLabel[n + 1] && newLabel[n + 2] >= oldLabel[n + 2]))
            newLabelDominated = false;
    }

    // here, at least one is true:
    //  if both new label and old label are dominated, they are the same

    // keep checking other elements in the label
    if (newLabelDominated)
        // if (newLabel[n] < oldLabel[n]) // check the number of visted nodes first
        //     return returnValue;
        // else // then check each node
        //     for (int i = 0; i < n; i++)
        //         if (newLabel[i] < oldLabel[i])
        //         {
        //             newLabelDominated = false;
        //             break;
        //         }
        for (int i = n; i >= 0; i--) // check label[n]first, which is the number of visited nodes
            if (newLabel[i] < oldLabel[i])
            {
                newLabelDominated = false;
                break;
            }

    if (oldLabelDominated)
        if (oldLabel[n] < newLabel[n]) // check the number of visted nodes first
            return returnValue;
        else
            for (int i = 0; i < n; i++)
                if (oldLabel[i] < newLabel[i])
                {
                    oldLabelDominated = false;
                    break;
                }

    if (newLabelDominated && oldLabelDominated) // two lables are the same
        returnValue = 1;

    if (newLabelDominated && !oldLabelDominated)
        returnValue = 2;

    if (!newLabelDominated && oldLabelDominated)
        returnValue = 3;

    return returnValue;
}

#endif