#ifndef DOMINANCE_INAB_SHORTLABEL_H
#define DOMINANCE_INAB_SHORTLABEL_H

// this version is based on dominance_inab.h. It test very short label (distance, cost, B2D, B2Dinab, predecessor), hopefully it can speed up running time
// but there are some problems after test. For example, the best route is 0-3-2-9, visited node label is (1 0 1 1 0 0 0 0 0 1), after checking other nodes'distance,
// it might be label like (1 M 1 1 M M M M M 1). If at node 2, there is a new label with the same cost and distnce, then this label will be removed, so the predecessor
// might be wrong and lead to a

// this version introduce bigM in labels

// this version dropped vector<vector<vector<double>>> newLabelSet, but only use vector<int>numNewlyAddedLabels (and vector<int> numNewlyAddedInLastIteration),
// which can refer to the index of labels added in last iteration and current iterations.It will make adding bigM to bigLambda_allNodesLableSet easier

// this version is basedon on dominance_doubleB2D_vector.h, with inaccessible nodes as 1 in labels, but since it's only changed in the current added label set, so it doesn't work right now.
// if we change corresponding labels in bigLambda_allNodesLableSet, it takes time to search. So I will drop the current added label set, and use only one set: bigLambda_allNodesLableSet

//  this version is based on dominance_doubleB2D.h, but we replace unordered_set with vector to try if it can speed up
// this version is based on dominance.h, but it has (n+4) elements in a label, including B2D value at the end of label
// dominance.h has only (n+3) elemens in a label.

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

void printLabelSet(vector<vector<double>>, double);
void printOneLabel(vector<double>, double);

void runDominance(int n, double **dis, vector<vector<double>> xCoeff, double disLimit, double *obj, vector<int> &nodesOnRoute)
{

    //===> set up parameters used in dominance.h
    int startingNode = 0;
    int endingNode = n - 1;
    bool PRINT_FOR_DEBUG = false;
    double bigM = 10000000; // bigM has to be big enough, so than when checking domincated labels inside F set, when comparing two equivalent labels, it helps to find label with more bigM
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
    long allOnes = 0;

    for (int i = 0; i < n; i++)
    {
        twoPow[i] = 1 << i;
        allOnes += twoPow[i];
    }

    if (PRINT_FOR_DEBUG)
        cout << "===> start dominace()" << endl;

    //===> !!! we reorder the label for node i to (si, Ti, Ci, B2D, B2Dinab, predecessor) so that
    //!!! label[0] is the consumed resouces;
    //!!! label[1] is the cost
    //!!! label[2] is the binary representation of the visited nodes B2D without inaccesible node bigM
    //!!! label[3] is the binary representation of the visited nodes B2D with bigM
    //!!! label[4] is the node before reaching current node in the label

    //===> a label for node i is (Ti, si, Vi1, Vi2, ..., Vin, Ci) in Feillet etc paper
    //"an exact algorithm for the elementary shortest path problem with resource constraints: application to some vehicle routing problems" in 2004
    //===> Ti: the resource distance used
    //===> si: the number of nodes visited
    //===> Vij: 1, node j is visited; 0, node j is not visited
    //===> Ci: the cost
    //===> the length of a lable is 1+1+n+1=n+3
    //===> add one more element B2D, the length is n+4 now
    //===> add B2Dinab and predecessor node, and remove Ti and number of visited nodes, the lenght is 5 now

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

    vector<vector<vector<double>>> bigLambda_allNodesLableSet;

    //===> initialize the lable for starting node
    vector<double> oneLabel;

    oneLabel.push_back(0);
    oneLabel.push_back(0);
    oneLabel.push_back(twoPow[n - startingNode - 1]);
    oneLabel.push_back(twoPow[n - startingNode - 1]);
    oneLabel.push_back(-1); //-1 means it has no predecessor

    vector<vector<double>> startingNodeLables;
    startingNodeLables.push_back(oneLabel);

    //===> add empty lable set to all other nodes
    for (int i = 0; i < n; i++)
    {
        if (i == startingNode)
            bigLambda_allNodesLableSet.push_back(startingNodeLables);
        else
        {
            vector<vector<double>> emptySet;
            bigLambda_allNodesLableSet.push_back(emptySet);
        }
    }

    // //===> initialize the newly added labels set as well
    // vector<vector<vector<double>>> newLabelSet(bigLambda_allNodesLableSet);

    unordered_set<int> E_toBeTreatedNodesSet;
    E_toBeTreatedNodesSet.insert(startingNode);

    int itrNum = 0; // iteration number

    // we are using vector<vector<vector<double>>> bigLambda_allNodesLableSet, so new labels added are alway in the end
    // therefore we can trace the number of newly added labels to tell which labels at the end are added in the last iteration
    vector<int> numNewlyAddedLabels(n, 0);
    numNewlyAddedLabels[startingNode] = 1;

    while (E_toBeTreatedNodesSet.size() > 0) //===> Line 7 in ESPPRC(p)
    {
        itrNum++;

        if (PRINT_FOR_DEBUG)
        {
            cout << endl;
            cout << "======> dominance iteration " << itrNum << " <======" << endl;
            cout << "E_toBeTreatedNodesSet size = " << E_toBeTreatedNodesSet.size() << endl;
            for (auto &temp : E_toBeTreatedNodesSet)
                cout << temp << " ";
            cout << endl;
        }

        vector<int> numNewlyAddedInLastIteration(numNewlyAddedLabels);
        for (int i = 0; i < n; i++)
            numNewlyAddedLabels[i] = 0;

        unordered_set<int> E_toBeTreatedNodesSetCopy(E_toBeTreatedNodesSet);

        // clear the original node set, and add the node with change later
        E_toBeTreatedNodesSet.clear();

        for (auto &node : E_toBeTreatedNodesSetCopy) //===> LINE 9 in ESPPRC(p)
        {

            // if (labelsAddedInLastIteration[node].size() < 1)
            // {
            //     printf("ERROR: There must be labels in E set!");
            //     exit(1);
            // }

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
                    printLabelSet(bigLambda_allNodesLableSet[node], bigM);
                    // cout << "===> labelsAddedInLastIteration[" << node << "].size =" << labelsAddedInLastIteration[node].size() << endl;
                    // printLabelSet(labelsAddedInLastIteration[node], bigM);
                    // cout << "===> newLabelSet[" << node << "].size =" << newLabelSet[node].size() << endl;
                    // printLabelSet(newLabelSet[node], bigM);
                }

                vector<vector<double>> F_nodeSuccLabelSet; //===> LINE 11 in ESPPRC(p) PSEUDO-CODE

                int nodeLabelSize = bigLambda_allNodesLableSet[node].size();

                for (int i = 0; i < numNewlyAddedInLastIteration[node]; i++) //===> LINE 12 in ESPPRC(p) PSEUDO-CODE
                {
                    vector<double> &label = bigLambda_allNodesLableSet[node][nodeLabelSize - 1 - i - numNewlyAddedLabels[node]];

                    long succB2D = twoPow[n - succ - 1];

                    // if (label[succ] == 0)
                    if (((long)label[3] & succB2D) == 0) //===> LINE 13 in ESPPRC(p) PSEUDO-CODE
                    {
                        if (PRINT_FOR_DEBUG)
                        {
                            cout << "number of labels added in the last iteration = " << numNewlyAddedInLastIteration[node] << endl;
                            cout << "find a label not visit succ yet" << endl;
                            printOneLabel(label, bigM);
                        }

                        //===> LINE 14 in ESPPRC(p) PSEUDO-CODE

                        //===> check resource (distance)
                        double currentDistance = label[0];

                        if (currentDistance + dis[node][succ] + dis[succ][endingNode] <= disLimit)
                        {
                            double distTemp = currentDistance + dis[node][succ]; //!!! label[n+1] is the consumed resouces;
                            double costTemp = label[1] + xCoeff[node][succ];     //!!! lable[n+2] is the cost
                            double B2Dtemp = label[2] + twoPow[n - succ - 1];
                            double B2Dinabtemp = label[3] + twoPow[n - succ - 1];

                            // check if it is dominated by any label in F set
                            bool dominatedTemp = false;
                            for (const auto &oldLabel : F_nodeSuccLabelSet)
                            {
                                if (distTemp >= oldLabel[0] && costTemp >= oldLabel[1])
                                    if (((long)B2Dinabtemp | (long)oldLabel[2]) == B2Dinabtemp)
                                    {
                                        dominatedTemp = true;
                                        break;
                                    }
                            }
                            if (!dominatedTemp)
                            {
                                vector<double> labelTemp;

                                labelTemp.push_back(distTemp); //!!! label[n+1] is the consumed resouces;
                                labelTemp.push_back(costTemp); //!!! lable[n+2] is the cost
                                labelTemp.push_back(B2Dtemp);
                                labelTemp.push_back(B2Dinabtemp);
                                labelTemp.push_back(node);

                                F_nodeSuccLabelSet.push_back(labelTemp);

                                if (PRINT_FOR_DEBUG)
                                {
                                    // cout << "visited number of nodes: labelTemp[n]=" << labelTemp[n] << endl;
                                    cout << "traveled distance: labelTemp[n + 1]=" << labelTemp[0] << endl;

                                    cout << "# of labels in F set: " << F_nodeSuccLabelSet.size() << endl;
                                    cout << "-> F set" << endl;
                                    printLabelSet(F_nodeSuccLabelSet, bigM);
                                }
                            }
                        }
                        else
                        {
                            label[3] += twoPow[n - succ - 1];

                            if (PRINT_FOR_DEBUG)
                                cout
                                    << "distance is over limit" << endl;
                            // labelsOverDistLimit.push_back(label);
                        }
                    }
                }

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
                    //===> step (1) (2) and (3) are combined as procedures below
                    for (auto newLabelIt = F_nodeSuccLabelSet.begin(); newLabelIt != F_nodeSuccLabelSet.end();)
                    {
                        bool getNewAddress = false;
                        int oldLabelVecIndex = 0;
                        for (auto oldLabelIt = bigLambda_allNodesLableSet[succ].begin(); oldLabelIt != bigLambda_allNodesLableSet[succ].end();)
                        {
                            int oldLabelSize = bigLambda_allNodesLableSet[succ].size();

                            if ((*newLabelIt)[0] >= (*oldLabelIt)[0] && (*newLabelIt)[1] >= (*oldLabelIt)[1])
                            {
                                if (PRINT_FOR_DEBUG)
                                {
                                    cout << "new label: ";
                                    printOneLabel(*newLabelIt, bigM);
                                    cout << "old label: ";
                                    printOneLabel(*oldLabelIt, bigM);
                                }

                                long result = ((long)(*newLabelIt)[3] | (long)(*oldLabelIt)[2]);
                                if (result == (long)(*newLabelIt)[3])
                                {
                                    newLabelIt = F_nodeSuccLabelSet.erase(newLabelIt);
                                    getNewAddress = true;
                                    break;

                                    if (PRINT_FOR_DEBUG)
                                        cout << "new label is dominated" << endl;
                                }
                                else
                                    ++oldLabelIt;
                            } // in case that old and new labels are the same, then only add in dominatedNewLabelSet
                            else if ((*newLabelIt)[0] <= (*oldLabelIt)[0] && (*newLabelIt)[1] <= (*oldLabelIt)[1]) // old label is possible to be dominated
                            {
                                long result = ((long)(*newLabelIt)[2] | (long)(*oldLabelIt)[3]);
                                if (result == (long)(*oldLabelIt)[3])
                                {
                                    // CHECK IF THE OLD LABEL IS ADDED IN THE CURRENT ITERATION WHEN CHECKING OTHER NODES' SUCCESSOR:
                                    // when oldLabelVecIndex< updated oldLabelSize, and oldLabelVecIndex>= updated oldLabelSize - updated number of newly added labels
                                    // then it means, the current old label is the one newly added. (newly added means added in the laster iteration) (oldLabelVecIndex starts from 0)
                                    // IF YES, THEN CHANGE THE NUMBER OF LABELS ADDED IN THIS ITERATION
                                    if (oldLabelVecIndex < oldLabelSize && oldLabelVecIndex >= (oldLabelSize - numNewlyAddedLabels[succ]))
                                        numNewlyAddedLabels[succ] -= 1;

                                    // CHECK IF THE OLD LABEL IS ADDED IN THE THE PREVIOUS  ITERATION WHEN CHECKING OTHER NODES' SUCCESSOR
                                    if (oldLabelVecIndex < (oldLabelSize - numNewlyAddedLabels[succ]) && oldLabelVecIndex >= (oldLabelSize - numNewlyAddedLabels[succ] - numNewlyAddedInLastIteration[succ]))
                                        numNewlyAddedInLastIteration[succ] -= 1;

                                    oldLabelVecIndex--;

                                    oldLabelIt = bigLambda_allNodesLableSet[succ].erase(oldLabelIt);

                                    if (PRINT_FOR_DEBUG)
                                        cout << "new label is dominated" << endl;
                                }
                                else
                                    ++oldLabelIt;
                            }
                            else
                                ++oldLabelIt;

                            oldLabelVecIndex++;
                        }
                        if (!getNewAddress)
                            ++newLabelIt;
                    }

                    //===> add the non-dominated labels to bigLambda_allNodesLableSet[succ]
                    for (auto &labelTemp : F_nodeSuccLabelSet)
                    {
                        bigLambda_allNodesLableSet[succ].push_back(labelTemp);
                        // newLabelSet[succ].push_back(labelTemp);
                        numNewlyAddedLabels[succ] += 1;
                    }

                    if (PRINT_FOR_DEBUG)
                    {
                        cout << "# of labels at succ node: " << bigLambda_allNodesLableSet[succ].size() << endl;
                        cout << "-> labels after update succ=" << succ << " labels" << endl;
                        printLabelSet(bigLambda_allNodesLableSet[succ], bigM);

                        // cout << "-> dominatedExistingLabelSet" << endl;
                        // printLabelSet(dominatedExistingLabelSet, bigM);

                        // cout << "-> dominatedNewLabelSet" << endl;
                        // printLabelSet(dominatedNewLabelSet, bigM);

                        cout << "-> F set before update" << endl;
                        printLabelSet(F_nodeSuccLabelSet, bigM);
                    }

                    // //===> step (4)
                    // //===> remove the dominated labels for current and succsor nodes

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
        if (labelTemp[1] < minCost)
        {
            bestLabel.clear();
            minCost = labelTemp[1];

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

    int numAllVistedNodes = 0;

    int currentNode = endingNode;
    int predecessor = bestLabel[4];
    double currentDist = bestLabel[0];
    double currentCost = bestLabel[1];
    vector<int> VNreverseOrder;
    int itrTemp = 0;

    while (predecessor >= 0) // starting node predecessor is -1; if use if(currentNode>=0), it won't work for starting node !=0, so that node 0 can be the intermediate node in route!
    {
        if (itrTemp > n)
        {
            exit(1);
        }
        itrTemp++;
        if (PRINT_FOR_DEBUG)
        {
            cout << "----------------" << endl;
            cout << "current node = " << currentNode << endl;
            cout << "predecessor = " << predecessor << endl;
        }

        if (numAllVistedNodes > n)
        {
            cout << "ERROR: two many nodes in route. " << endl;
            exit(1);
        }
        double distTemp = currentDist - dis[predecessor][currentNode];
        double costTemp = currentCost - xCoeff[predecessor][currentNode];
        for (auto &labelTemp : bigLambda_allNodesLableSet[predecessor])
        {
            // cout << "---a new label in set---" << endl;
            if (labelTemp[0] < distTemp + verySmallNum && labelTemp[0] > distTemp - verySmallNum)
                if (labelTemp[1] < costTemp + verySmallNum && labelTemp[1] > costTemp - verySmallNum)
                {
                    // bool temp = ((long)labelTemp[2] & twoPow[n - predecessor - 1]);
                    // cout << "labelTemp[2]=" << labelTemp[2] << endl;
                    // cout << "temp=" << temp << endl;
                    // if (temp == 1)
                    // {
                    VNreverseOrder.push_back(currentNode);
                    currentNode = predecessor;
                    predecessor = labelTemp[4];
                    currentDist = distTemp;
                    currentCost = costTemp;
                    numAllVistedNodes++;

                    break;
                    // }
                }
        }
    }
    VNreverseOrder.push_back(startingNode);

    for (int i = 0; i < numAllVistedNodes + 1; i++)
        nodesOnRoute.push_back(0);

    for (int i = 0; i < VNreverseOrder.size(); i++)
    {
        nodesOnRoute[numAllVistedNodes - i] = VNreverseOrder[i];
    }

    if (PRINT_FOR_DEBUG)
    {
        cout << "selected route:" << endl;
        for (auto &ele : nodesOnRoute)
            cout << ele << " ";
        cout << endl;
    }
}

void printLabelSet(vector<vector<double>> labelSet, double bigM)
{
    for (auto &labelTemp : labelSet)
    {
        for (auto &temp : labelTemp)
            // cout << fixed << setprecision(2) << temp << " ";
            if (temp == bigM)
                cout << "M ";
            else
                cout << temp << " ";
        cout << endl;
    }
}

void printOneLabel(vector<double> label, double bigM)
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