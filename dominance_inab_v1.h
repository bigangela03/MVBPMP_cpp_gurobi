#ifndef DOMINANCE_INAB_V1_H
#define DOMINANCE_INAB_V1_H

// this version is based on dominance_inab_v0.h. it dropped vector<vector<vector<double>>> newLabelSet, but only use vector<int>numNewlyAddedLabels (and vector<int> numNewlyAddedInLastIteration),
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
bool compareToLabel(vector<double>, vector<double>, int, double);

void runDominance(int n, double **dis, vector<vector<double>> xCoeff, double disLimit, double *obj, vector<int> &nodesOnRoute)
{

    //===> set up parameters used in dominance.h
    int startingNode = 0;
    int endingNode = n - 1;
    bool PRINT_FOR_DEBUG = false;
    bool CHECK_NODE_LABELS_END_OF_ITERATION = true;
    double bigM = 10000000; // bigM has to be big enough, so than when checking domincated labels inside F set, when comparing two equivalent labels, it helps to find label with more bigM
    double verySmallNum = 0.0000001;

    // vector<double> targetLabel;
    // targetLabel.push_back(1);
    // targetLabel.push_back(0);
    // targetLabel.push_back(0);
    // targetLabel.push_back(0);
    // targetLabel.push_back(3);
    // targetLabel.push_back(bigM);
    // targetLabel.push_back(2);
    // targetLabel.push_back(0);
    // targetLabel.push_back(0);
    // targetLabel.push_back(0);

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
        twoPow[i] = pow(2, i);
        allOnes += twoPow[i];
    }

    if (PRINT_FOR_DEBUG)
        cout << "===> start dominace()" << endl;

    //===> !!! we reorder the label for node i to (Vi1, Vi2, ..., Vin, si, Ti, Ci, B2D) so that
    //!!! label[0] to lable[n-1] represents visited nodes
    //!!! label[n] is the number of visited nodes;
    //!!! label[n+1] is the consumed resouces;
    //!!! label[n+2] is the cost
    //!!! label[n+3] is the binary representation of the visited nodes B2D

    //===> a label for node i is (Ti, si, Vi1, Vi2, ..., Vin, Ci) in Feillet etc paper
    //"an exact algorithm for the elementary shortest path problem with resource constraints: application to some vehicle routing problems" in 2004
    //===> Ti: the resource distance used
    //===> si: the number of nodes visited
    //===> Vij: 1, node j is visited; 0, node j is not visited
    //===> Ci: the cost
    //===> the length of a lable is 1+1+n+1=n+3
    //===> add one more element B2D, the length is n+4 now

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
    for (int i = 0; i < n + 4; i++)
        oneLabel.push_back(0);
    oneLabel[startingNode] = 1;
    oneLabel[n] = 1;
    oneLabel[n + 3] = twoPow[n - startingNode - 1];

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

        // //===> only check the labels newly added in the last iteration
        // vector<vector<vector<double>>> labelsAddedInLastIteration(newLabelSet);
        // for (int i = 0; i < newLabelSet.size(); i++)
        //     newLabelSet[i].clear();

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

                // for (auto &label : labelsAddedInLastIteration[node]) //===> LINE 12 in ESPPRC(p) PSEUDO-CODE
                for (int i = 0; i < numNewlyAddedInLastIteration[node]; i++)
                {
                    vector<double> &label = bigLambda_allNodesLableSet[node][nodeLabelSize - 1 - i - numNewlyAddedLabels[node]];

                    if (label[succ] == 0) //===> LINE 13 in ESPPRC(p) PSEUDO-CODE
                    {
                        if (PRINT_FOR_DEBUG)
                        {
                            cout << "number of labels added in the last iteration = " << numNewlyAddedInLastIteration[node] << endl;
                            cout << "find a label not visit succ yet" << endl;
                            printOneLabel(label, bigM);
                        }

                        //===> LINE 14 in ESPPRC(p) PSEUDO-CODE

                        //===> check resource (distance)
                        double currentDistance = label[n + 1];

                        if (currentDistance + dis[node][succ] + dis[succ][endingNode] <= disLimit)
                        {

                            double distTemp = currentDistance + dis[node][succ]; //!!! label[n+1] is the consumed resouces;
                            double costTemp = label[n + 2] + xCoeff[node][succ]; //!!! lable[n+2] is the cost
                            double B2Dtemp = label[n + 3] + twoPow[n - succ - 1];

                            // check if it is dominated by any label in F set
                            bool dominatedTemp = false;
                            for (const auto &oldLabel : F_nodeSuccLabelSet)
                            {
                                if (distTemp >= oldLabel[n + 1] && costTemp >= oldLabel[n + 2])
                                    if (((long)B2Dtemp | (long)oldLabel[n + 3]) == B2Dtemp)
                                    {
                                        dominatedTemp = true;
                                        break;
                                    }
                            }
                            if (!dominatedTemp)
                            {
                                vector<double> labelTemp;
                                for (int j = 0; j < n + 4; j++)
                                    labelTemp.push_back(label[j]);

                                int numVisitedNodes = label[n];
                                // labelTemp[succ] = 1;                                  //!!! lable[0] to lable[n-1] represents visited nodes
                                labelTemp[succ] = numVisitedNodes + 1;
                                labelTemp[n] = numVisitedNodes + 1; //!!! lable[n] is the number of visited nodes;
                                labelTemp[n + 1] = distTemp;        //!!! label[n+1] is the consumed resouces;
                                labelTemp[n + 2] = costTemp;        //!!! lable[n+2] is the cost
                                labelTemp[n + 3] = B2Dtemp;
                                F_nodeSuccLabelSet.push_back(labelTemp);

                                if (PRINT_FOR_DEBUG)
                                {
                                    cout << "visited number of nodes: labelTemp[n]=" << labelTemp[n] << endl;
                                    cout << "traveled distance: labelTemp[n + 1]=" << labelTemp[n + 1] << endl;

                                    cout << "# of labels in F set: " << F_nodeSuccLabelSet.size() << endl;
                                    cout << "-> F set" << endl;
                                    printLabelSet(F_nodeSuccLabelSet, bigM);
                                }
                            }
                        }
                        else
                        {
                            // label[succ] = bigM;
                            // label[n + 3] = label[n + 3] + twoPow[n - succ - 1];

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

                            if ((*newLabelIt)[n + 1] >= (*oldLabelIt)[n + 1] && (*newLabelIt)[n + 2] >= (*oldLabelIt)[n + 2])
                            {
                                if (PRINT_FOR_DEBUG)
                                {
                                    cout << "new label: ";
                                    printOneLabel(*newLabelIt, bigM);
                                    cout << "old label: ";
                                    printOneLabel(*oldLabelIt, bigM);
                                }

                                long result = ((long)(*newLabelIt)[n + 3] | (long)(*oldLabelIt)[n + 3]);
                                if (result == (*newLabelIt)[n + 3])
                                {
                                    // dominatedNewLabelSet.push_back(*newLabelIt);
                                    newLabelIt = F_nodeSuccLabelSet.erase(newLabelIt);
                                    getNewAddress = true;
                                    break;

                                    if (PRINT_FOR_DEBUG)
                                        cout << "new label is dominated" << endl;
                                }
                                else
                                    ++oldLabelIt;
                            } // in case that old and new labels are the same, then only add in dominatedNewLabelSet
                            else if ((*newLabelIt)[n + 1] <= (*oldLabelIt)[n + 1] && (*newLabelIt)[n + 2] <= (*oldLabelIt)[n + 2]) // old label is possible to be dominated
                            {
                                long result = ((long)(*newLabelIt)[n + 3] | (long)(*oldLabelIt)[n + 3]);
                                if (result == (*oldLabelIt)[n + 3])
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

                                    // dominatedExistingLabelSet.push_back(*oldLabelIt);
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
            if (CHECK_NODE_LABELS_END_OF_ITERATION)
            {
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

bool compareToLabel(vector<double> lab, vector<double> target, int n, double verySmallNum)
{
    bool result = true;
    for (int i = 0; i < n; i++)
    {
        if (lab[i] > target[i] + verySmallNum || lab[i] < target[i] - verySmallNum)
            result = false;
    }
    return result;
}

#endif