#ifndef DOMINANCE_BIT_HASHMAP_H
#define DOMINANCE_BIT_HASHMAP_H

//!!!!!!!!!!!!!!!!!!! NOT WORKING !!!!!!!!!!!!!!!!!!!!!!!!!!
// the data structre of hashmap is not working now!!!!!!!!!!!

// I generate only one hashmap for all labels of all nodes; the hashmap maps key (binary representation of visited nodes to decimal) to visitited nodes, distance, and cost array.
//  But there is a problem.
// for a graph with node 0, 1, 2, ..., starting from node 0, at iteration 1, the label for node 1 is (1 1 0 ...0,d01,c01), label at node2 (1 0 1 0 ...0,d02,c02);
// in iteration 2, check node 1 first, its successor node2 has new label (1 1 1 0...0,d01+d12,c01+c12); then check node2, its successor node1 has new label (1 1 1 0...1,d02+d21,c02+c21)
// when generating keys, both (1 1 1 0...0,d01+d12,c01+c12) and (1 1 1 0...1,d02+d21,c02+c21) has the same key 2^(n-1)+2^(n-2)+2^(n-3)
// when (1 1 1 0...1,d02+d21,c02+c21) is added, its key will replace the old key, which cause bugs

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

void printLabelSet(unordered_set<long> &, unordered_map<long, vector<double>> &, double);
void printOneLabel(vector<double> &, double);

void runDominance(int n, double **dis, vector<vector<double>> xCoeff, double disLimit, double *obj, vector<int> &nodesOnRoute)
{

    //===> set up parameters used in dominance.h
    int startingNode = 0;
    int endingNode = n - 1;
    bool PRINT_FOR_DEBUG = true;
    double bigM = 10000000; // bigM has to be big enough, so than when checking domincated labels inside F set, when comparing two equivalent labels, it helps to find label with more bigM
    double verySmallNum = 0.0000001;

    if (PRINT_FOR_DEBUG)
        cout << "===> start dominace()" << endl;

    //===> !!! we reorder the label for node i to (Vi1, Vi2, ..., Vin, si, Ti, Ci) so that
    //!!! lable[0] to lable[n-1] represents visited nodes
    //!!! lable[n] is the number of visited nodes;
    //!!! label[n+1] is the consumed resouces;
    //!!! lable[n+2] is the cost

    //===> a label for node i is (Ti, si, Vi1, Vi2, ..., Vin, Ci) in Feillet etc paper
    //"an exact algorithm for the elementary shortest path problem with resource constraints: application to some vehicle routing problems" in 2004
    //===> Ti: the resource distance used
    //===> si: the number of nodes visited
    //===> Vij: 1, node j is visited; 0, node j is not visited
    //===> Ci: the cost
    //===> the length of a lable is 1+1+n+1=n+3

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

    long *twoPow = new long[n];
    for (int i = 0; i < n; i++)
        twoPow[i] = pow(2, i);

    vector<unordered_set<long>> bigLambda_allNodesLableSet;
    unordered_map<long, vector<double>> key2routeMap;

    //===> initialize the lable for starting node
    vector<double> oneLabel;
    for (int i = 0; i < n + 3; i++)
        oneLabel.push_back(0);
    oneLabel[startingNode] = 1;
    oneLabel[n] = 1;

    unordered_set<long> startingNodeKeys;

    long startingNodeKey = twoPow[n - startingNode - 1]; // B2D, binary to decimal
    startingNodeKeys.insert(startingNodeKey);
    key2routeMap[startingNodeKey] = oneLabel;

    //===> add empty lable set to all other nodes
    for (int i = 0; i < n; i++)
    {
        if (i == startingNode)
            bigLambda_allNodesLableSet.push_back(startingNodeKeys);
        else
        {
            unordered_set<long> emptySet;
            bigLambda_allNodesLableSet.push_back(emptySet);
        }
    }

    //===> initialize the newly added labels set as well
    vector<unordered_set<long>> newLabelSet(bigLambda_allNodesLableSet);

    // unordered_set<int> E_toBeTreatedNodesSet;
    set<int> E_toBeTreatedNodesSet;
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

        // unordered_set<int> E_toBeTreatedNodesSetCopy(E_toBeTreatedNodesSet);
        set<int> E_toBeTreatedNodesSetCopy(E_toBeTreatedNodesSet);

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
                    printLabelSet(bigLambda_allNodesLableSet[node], key2routeMap, bigM);
                    cout << "===> labelsAddedInLastIteration[" << node << "].size =" << labelsAddedInLastIteration[node].size() << endl;
                    printLabelSet(labelsAddedInLastIteration[node], key2routeMap, bigM);
                    cout << "===> newLabelSet[" << node << "].size =" << newLabelSet[node].size() << endl;
                    printLabelSet(newLabelSet[node], key2routeMap, bigM);
                }

                unordered_set<long> F_nodeSuccLabelSet; //===> LINE 11 in ESPPRC(p) PSEUDO-CODE
                                                        // unordered_set<long> labelsOverDistLimit;

                for (auto &key : labelsAddedInLastIteration[node]) //===> LINE 12 in ESPPRC(p) PSEUDO-CODE
                {

                    vector<double> &label = key2routeMap[key];

                    if (label[succ] == 0) //===> LINE 13 in ESPPRC(p) PSEUDO-CODE
                    {
                        if (PRINT_FOR_DEBUG)
                        {
                            cout << "find a label not visit succ yet" << endl;
                            printOneLabel(label, bigM);
                        }

                        //===> LINE 14 in ESPPRC(p) PSEUDO-CODE

                        //===> check resource (distance)
                        double currentDistance = label[n + 1];

                        if (currentDistance + dis[node][succ] + dis[succ][endingNode] <= disLimit)
                        {
                            vector<double> labelTemp;
                            for (int i = 0; i < n + 3; i++)
                                labelTemp.push_back(label[i]);

                            int numVisitedNodes = label[n];
                            // labelTemp[succ] = 1;                                  //!!! lable[0] to lable[n-1] represents visited nodes
                            labelTemp[succ] = numVisitedNodes + 1;
                            labelTemp[n] = numVisitedNodes + 1;                   //!!! lable[n] is the number of visited nodes;
                            labelTemp[n + 1] = currentDistance + dis[node][succ]; //!!! label[n+1] is the consumed resouces;
                            labelTemp[n + 2] = label[n + 2] + xCoeff[node][succ]; //!!! lable[n+2] is the cost
                            long keyTemp = key + twoPow[n - succ - 1];
                            F_nodeSuccLabelSet.insert(keyTemp);
                            key2routeMap[keyTemp] = labelTemp;

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
                        printLabelSet(F_nodeSuccLabelSet, key2routeMap, bigM);
                    }

                    // vector<vector<double>> Fset;
                    // for (auto &labelTemp : F_nodeSuccLabelSet)
                    // {
                    //     vector<double> labelToAdd;
                    //     for (auto &ele : labelTemp)
                    //         labelToAdd.push_back(ele);
                    //     Fset.push_back(labelToAdd);
                    // }

                    unordered_set<long> dominatedLabelsAmongNewLabelSet;

                    // if (Fset.size() >= 2)
                    if (F_nodeSuccLabelSet.size() >= 2)
                    { // since we are checking lables in a Set, so no labels are duplicate according to Set definition
                        // for (int i = 0; i < Fset.size(); i++)
                        //     for (int j = i + 1; j < Fset.size(); j++)
                        for (auto &key1 : F_nodeSuccLabelSet)
                            for (auto &key2 : F_nodeSuccLabelSet)
                                if (key1 != key2)
                                {

                                    vector<double> &oldLabel = key2routeMap[key1];
                                    vector<double> &newLabel = key2routeMap[key2];

                                    if (!(newLabel[n + 1] <= oldLabel[n + 1] && newLabel[n + 2] <= oldLabel[n + 2]))     // old label is not dominated
                                        if (!(newLabel[n + 1] >= oldLabel[n + 1] && newLabel[n + 2] >= oldLabel[n + 2])) // new label is not dominated
                                            continue;
                                    // now, at least one label has possibility to be dominated
                                    long result = (key1 | key2);
                                    if (result == key1)
                                        dominatedLabelsAmongNewLabelSet.insert(key1);
                                    else if (result == key2)
                                        dominatedLabelsAmongNewLabelSet.insert(key2);
                                }
                    }

                    //===> step (2)
                    //===> remove the dominated labels from F_nodeSuccLabelSet
                    for (auto &keyTemp : dominatedLabelsAmongNewLabelSet)
                    {
                        F_nodeSuccLabelSet.erase(keyTemp);
                        // key2routeMap.erase(keyTemp);
                    }

                    if (PRINT_FOR_DEBUG)
                    {
                        if (dominatedLabelsAmongNewLabelSet.size() > 0)
                        {
                            cout << "-> dominatedLabelsAmongNewLabelSet" << endl;

                            printLabelSet(dominatedLabelsAmongNewLabelSet, key2routeMap, bigM);

                            cout << "# of labels in F set after update: " << F_nodeSuccLabelSet.size() << endl;
                            cout << "-> F set after removing dominated labels" << endl;
                            printLabelSet(F_nodeSuccLabelSet, key2routeMap, bigM);
                        }
                        else
                            cout << "no dominated lables found among F set" << endl;
                    }

                    //===> step (3)

                    unordered_set<long> dominatedNewLabelSet;
                    unordered_set<long> dominatedExistingLabelSet;

                    for (auto &newKey : F_nodeSuccLabelSet)
                    {

                        // we can compare dist and cost first, if  not (d1>=d2, c1>=c2, or d1<=d2, c1<=c2), continue

                        // then, use some other methods to decide, like if any element violates dominated conditin, then break, so we don't need to check all the elements in a label.

                        for (auto &oldKey : bigLambda_allNodesLableSet[succ])
                        {
                            if (PRINT_FOR_DEBUG)
                            {
                                cout << "new label key " << newKey << ": ";
                                printOneLabel(key2routeMap[newKey], bigM);
                                cout << "old label key " << oldKey << ": ";
                                printOneLabel(key2routeMap[oldKey], bigM);
                            }
                            if (newKey != oldKey)
                            {

                                vector<double> &oldLabel = key2routeMap[oldKey];
                                vector<double> &newLabel = key2routeMap[newKey];

                                if (!(newLabel[n + 1] <= oldLabel[n + 1] && newLabel[n + 2] <= oldLabel[n + 2]))     // old label is not dominated
                                    if (!(newLabel[n + 1] >= oldLabel[n + 1] && newLabel[n + 2] >= oldLabel[n + 2])) // new label is not dominated
                                        continue;
                                // now, at least one label is dominated
                                long result = (newKey | oldKey);
                                if (result == newKey)
                                    dominatedNewLabelSet.insert(newKey);
                                else if (result == oldKey)
                                    dominatedExistingLabelSet.insert(oldKey);
                                if (PRINT_FOR_DEBUG)
                                {
                                    if (result == newKey)
                                        cout << "New label is dominated by old label." << endl;
                                    else if (result == oldKey)
                                        cout << "Old label is dominated." << endl;
                                }
                            }
                            else
                                // we need to add the same lable in dominated new labels to to remove it from F set later, so it won't be in the newly added labels in next iteration
                                dominatedNewLabelSet.insert(newKey);
                        }
                    }
                    if (PRINT_FOR_DEBUG)
                    {
                        cout << "# of labels at succ node: " << bigLambda_allNodesLableSet[succ].size() << endl;
                        cout << "-> labels before update succ=" << succ << " labels" << endl;
                        printLabelSet(bigLambda_allNodesLableSet[succ], key2routeMap, bigM);

                        cout << "-> dominatedExistingLabelSet" << endl;
                        printLabelSet(dominatedExistingLabelSet, key2routeMap, bigM);

                        cout << "-> dominatedNewLabelSet" << endl;
                        printLabelSet(dominatedNewLabelSet, key2routeMap, bigM);

                        cout << "-> F set before update" << endl;
                        printLabelSet(F_nodeSuccLabelSet, key2routeMap, bigM);
                    }

                    //===> step (4)
                    //===> remove the dominated labels for current and succsor nodes
                    for (auto &keyTemp : dominatedExistingLabelSet)
                    {
                        bigLambda_allNodesLableSet[succ].erase(keyTemp);
                        // key2routeMap.erase(keyTemp);
                    }
                    for (auto &keyTemp : dominatedNewLabelSet)
                    {
                        F_nodeSuccLabelSet.erase(keyTemp);
                        // key2routeMap.erase(keyTemp);
                    }

                    //===> add the non-dominated labels to bigLambda_allNodesLableSet[succ]
                    for (auto &keyTemp : F_nodeSuccLabelSet)
                    {
                        bigLambda_allNodesLableSet[succ].insert(keyTemp);
                        newLabelSet[succ].insert(keyTemp);
                    }

                    if (PRINT_FOR_DEBUG)
                    {
                        // cout << "======> ERROR: a conflict found!" << endl;
                        cout << "-> F_nodeSuccLabelSet" << endl;
                        printLabelSet(F_nodeSuccLabelSet, key2routeMap, bigM);

                        cout << "# of labels at succ node: " << bigLambda_allNodesLableSet[succ].size() << endl;
                        cout << "-> labels after update succ=" << succ << " labels" << endl;
                        printLabelSet(bigLambda_allNodesLableSet[succ], key2routeMap, bigM);
                    }

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
        printLabelSet(bigLambda_allNodesLableSet[n - 1], key2routeMap, bigM);
    }

    double minCost = numeric_limits<double>::max();
    vector<double> bestLabel;

    for (auto &keyTemp : bigLambda_allNodesLableSet[n - 1])
    {
        vector<double> &labelTemp = key2routeMap[keyTemp];
        if (labelTemp[n + 2] < minCost)
        {
            bestLabel.clear();
            minCost = labelTemp[n + 2];
            bestLabel = labelTemp;
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

void printLabelSet(unordered_set<long> &keySet, unordered_map<long, vector<double>> &key2routeMap, double bigM)
{
    for (auto &key : keySet)
    {
        for (auto &temp : key2routeMap[key])
            // cout << fixed << setprecision(2) << temp << " ";
            if (temp == bigM)
                cout << "M ";
            else
                cout << temp << " ";
        cout << endl;
    }
}

void printOneLabel(vector<double> &label, double bigM)
{
    for (auto &temp : label)
        // cout << fixed << setprecision(2) << temp << " ";
        if (temp == bigM)
            cout << "M ";
        else
            cout << temp << " ";
    cout << endl;
}

#endif