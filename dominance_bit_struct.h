#ifndef DOMINANCE_BIT_STRUCT_H
#define DOMINANCE_BIT_STRUCT_H

//!!!!!!!!!!!!!!!!!!! NOT WORKING !!!!!!!!!!!!!!!!!!!!!!!!!!
// the data structre of hashmap is not working now!!!!!!!!!!!

//!!!!!! READE ME BEFORE USING THE CODE !!!!!!
// This code has bug and can't find the correct solution right now.
// It is copied from dominance.h and then changed the code line by line by using struct Label instead of double array to represent a label
// The possible bug is that at the beginning of while loop, the newly added label vector of labels is copied to another vector, and then cleared for the new non-dominated labels in currect iteration
// But when copying, I am not sure if c++ does the deep or shallow copy. It needs time to figure out
// Then it occured to me that I actually don't need such complicated data structure like "struct", I can use the decimal value of binary representation of visited code (such as BVN=21 for (10101))
// I can use BVN alone in unordered_set, then do the hashmap from BVN to distance and a seperate hashmap from BVN to cost. OH Yeah!!!!:)

// One bug I found later is as follows:
//  bool operator<(const Label &other) const
//      {
//          return BVN < other.BVN;
//      }
// BVN actually is not unique , some migh thave the same BVN. Please find details in struct Label below!

// This code is only for complete graph right now, i.e., any two nodes are connected
// since I didn't consider the vertex incidence matrix, only iterate nodes from 0 to n-1

// #include <stdio.h>
// #include <stdlib.h>
// #include <unistd.h> //to get path to current directory
#include <cmath>
#include <time.h>
#include <malloc.h>
#include <string>

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
#include <unordered_set>
#include <string>
#include <ctime>  //to measure CPU time
#include <chrono> //to measure run time
#include <unordered_map>

#include <limits> //to get the biggest value of double variables
using namespace std;
using namespace std::chrono;
using namespace std::string_literals;

struct Label
{
    vector<int> visitedNodes;
    int numVisitedNodes;
    double distance;
    double cost;
    long BVN; // binary visited ndoes, for example,(10101) means visiting node 0, 2 and 4, the decimanl is 21, so BVN=21

    // bool operator==(const Label &other) const
    // {
    //     return BVN == other.BVN;
    // }

    // here is a problem when using struct in set
    // we need to define the key to keep uniqueness in a set.
    // but BVN might be the same for many labels, like (1023004) and (1032004) are both in (1011001) in BVN
    // so when the BVN is the same, the new label wont't be inserted because there is already one element
    // it is hard to figure out a good way to make each label unique in set.
    // The same problem happens to unordered_set of struct for operator==

    bool operator<(const Label &other) const
    {
        return BVN < other.BVN;
    }
};

// class for hash function
// class LabelHashFunction
// {
// public:
//     // BVN is returned as hash function
//     size_t operator()(const Label &lab) const
//     {
//         return lab.BVN;
//     }
// };

// Passing by reference (&) avoids unnecessary copying of the unordered_set.
// void printLabelSet(unordered_set<Label, LabelHashFunction> &labelSet, double bigM)
void printLabelSet(set<Label> &labelSet, double bigM)
{
    for (auto &labelTemp : labelSet)
    {
        for (auto &temp : labelTemp.visitedNodes)
            // cout << fixed << setprecision(2) << temp << " ";
            if (temp == bigM)
                cout << "M ";
            else
                cout << temp << " ";
        // cout << "| " << labelTemp.numVisitedNodes << " ";
        cout << labelTemp.numVisitedNodes << " ";
        cout << labelTemp.distance << " ";
        cout << labelTemp.cost << " " << endl;
        // cout << labelTemp.BVN << endl;
    }
}

void printOneLabel(Label label, double bigM)
{
    for (auto &temp : label.visitedNodes)
        // cout << fixed << setprecision(2) << temp << " ";
        if (temp == bigM)
            cout << "M ";
        else
            cout << temp << " ";
    cout << label.numVisitedNodes << " ";
    cout << label.distance << " ";
    cout << label.cost << " " << endl;
    // cout << label.BVN << endl;
}

void runDominance(int n, double **dis, vector<vector<double>> xCoeff, double disLimit, double *obj, vector<int> &nodesOnRoute)
{

    //===> set up parameters used in dominance.h
    int startingNode = 0;
    int endingNode = n - 1;
    bool PRINT_FOR_DEBUG = true;
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

    for (int i = 0; i < n; i++)
        twoPow[i] = pow(2, i);

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

    vector<set<Label>> bigLambda_allNodesLableSet;
    // vector<unordered_set<Label, LabelHashFunction>> bigLambda_allNodesLableSet;

    //===> initialize the lable for starting node

    struct Label startingNodeLabel;
    startingNodeLabel.cost = 0;
    startingNodeLabel.distance = 0;

    vector<int> visitedNodeVec;
    for (int i = 0; i < n; i++)
        visitedNodeVec.push_back(0);
    visitedNodeVec[startingNode] = 1;

    startingNodeLabel.numVisitedNodes = 1;
    startingNodeLabel.visitedNodes = visitedNodeVec;
    startingNodeLabel.BVN = twoPow[n - startingNode - 1];

    // unordered_set<Label, LabelHashFunction> startingNodeLables;
    set<Label> startingNodeLables;
    startingNodeLables.insert(startingNodeLabel);

    //===> add empty lable set to all other nodes
    for (int i = 0; i < n; i++)
    {
        if (i == startingNode)
            bigLambda_allNodesLableSet.push_back(startingNodeLables);
        else
        {
            // unordered_set<Label, LabelHashFunction> emptySet;
            set<Label> emptySet;
            bigLambda_allNodesLableSet.push_back(emptySet);
        }
    }

    //===> initialize the newly added labels set as well
    // vector<unordered_set<Label, LabelHashFunction>> newLabelSet(bigLambda_allNodesLableSet);
    vector<set<Label>> newLabelSet(bigLambda_allNodesLableSet);

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
        // vector<unordered_set<Label, LabelHashFunction>> labelsAddedInLastIteration(newLabelSet);
        vector<set<Label>> labelsAddedInLastIteration(newLabelSet);

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

            // ===> ending node  has no successors
            if (node == endingNode)
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
                    cout << "===> labelsAddedInLastIteration[" << node << "].size =" << labelsAddedInLastIteration[node].size() << endl;
                    printLabelSet(labelsAddedInLastIteration[node], bigM);
                    cout << "===> newLabelSet[" << node << "].size =" << newLabelSet[node].size() << endl;
                    printLabelSet(newLabelSet[node], bigM);
                }

                // unordered_set<Label, LabelHashFunction> F_nodeSuccLabelSet; //===> LINE 11 in ESPPRC(p) PSEUDO-CODE
                set<Label> F_nodeSuccLabelSet; // unordered_set<Label, LabelHashFunction> labelsOverDistLimit;

                for (auto &label : labelsAddedInLastIteration[node]) //===> LINE 12 in ESPPRC(p) PSEUDO-CODE
                {

                    if (label.visitedNodes[succ] == 0) //===> LINE 13 in ESPPRC(p) PSEUDO-CODE
                    {
                        if (PRINT_FOR_DEBUG)
                        {
                            cout << "find a label not visit succ yet" << endl;
                            printOneLabel(label, bigM);
                        }

                        //===> LINE 14 in ESPPRC(p) PSEUDO-CODE

                        //===> check resource (distance)
                        double currentDistance = label.distance;

                        if (currentDistance + dis[node][succ] + dis[succ][endingNode] <= disLimit)
                        {
                            Label labelTemp;
                            for (int i = 0; i < n; i++)
                                labelTemp.visitedNodes.push_back(label.visitedNodes[i]);

                            int numTemp = label.numVisitedNodes;
                            labelTemp.visitedNodes[succ] = numTemp + 1;
                            labelTemp.numVisitedNodes = numTemp + 1;
                            labelTemp.distance = currentDistance + dis[node][succ];
                            labelTemp.cost = label.cost + xCoeff[node][succ];
                            labelTemp.BVN = label.BVN + twoPow[n - succ - 1];

                            // check if it is dominated
                            bool dominatedTemp = false;
                            for (const auto &oldLabel : F_nodeSuccLabelSet)
                            {
                                if (labelTemp.distance >= oldLabel.distance && labelTemp.cost >= oldLabel.cost)
                                    if ((labelTemp.BVN | oldLabel.BVN) == labelTemp.BVN)
                                    {
                                        dominatedTemp = true;
                                        break;
                                    }
                            }
                            if (!dominatedTemp)
                                F_nodeSuccLabelSet.insert(labelTemp);

                            if (PRINT_FOR_DEBUG)
                            {
                                cout << "visited number of nodes: labelTemp[n]=" << labelTemp.numVisitedNodes << endl;
                                cout << "traveled distance: labelTemp[n + 1]=" << labelTemp.distance << endl;

                                cout << "# of labels in F set: " << F_nodeSuccLabelSet.size() << endl;
                                cout << "-> F set" << endl;
                                printLabelSet(F_nodeSuccLabelSet, bigM);
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

                    // the process in step (1) and (2) are moved to the procedure when creating new labels at successors
                    //  //===> step (1)
                    //  if (PRINT_FOR_DEBUG)
                    //  {
                    //      cout << "# of labels in F set: " << F_nodeSuccLabelSet.size() << endl;
                    //      cout << "-> F set before removing dominated labels" << endl;
                    //      printLabelSet(F_nodeSuccLabelSet, bigM);
                    //  }

                    // // unordered_set<Label, LabelHashFunction> dominatedLabelsAmongNewLabelSet;
                    // set<Label> dominatedLabelsAmongNewLabelSet;

                    // if (F_nodeSuccLabelSet.size() >= 2)
                    // {
                    //     for (const auto &newLabel : F_nodeSuccLabelSet)
                    //         for (const auto &oldLabel : F_nodeSuccLabelSet)
                    //         {
                    //             if (newLabel.distance >= oldLabel.distance && newLabel.cost >= oldLabel.cost)
                    //             {
                    //                 long result = (newLabel.BVN | oldLabel.BVN);
                    //                 if (result == newLabel.BVN)
                    //                     dominatedLabelsAmongNewLabelSet.insert(newLabel);
                    //             } // in case that old and new labels are the same, then only add in dominatedNewLabelSet
                    //             else if (newLabel.distance <= oldLabel.distance && newLabel.cost <= oldLabel.cost) // old label is possible to be dominated
                    //             {
                    //                 long result = (newLabel.BVN | oldLabel.BVN);
                    //                 if (result == oldLabel.BVN)
                    //                     dominatedLabelsAmongNewLabelSet.insert(oldLabel);
                    //             }
                    //         }
                    // }

                    // //===> step (2)
                    // //===> remove the dominated labels from F_nodeSuccLabelSet
                    // for (auto &labelTemp : dominatedLabelsAmongNewLabelSet)
                    //     F_nodeSuccLabelSet.erase(labelTemp);

                    // if (PRINT_FOR_DEBUG)
                    // {
                    //     if (dominatedLabelsAmongNewLabelSet.size() > 0)
                    //     {
                    //         cout << "-> dominatedLabelsAmongNewLabelSet" << endl;

                    //         printLabelSet(dominatedLabelsAmongNewLabelSet, bigM);

                    //         cout << "# of labels in F set after update: " << F_nodeSuccLabelSet.size() << endl;
                    //         cout << "-> F set after removing dominated labels" << endl;
                    //         printLabelSet(F_nodeSuccLabelSet, bigM);
                    //     }
                    //     else
                    //         cout << "no dominated lables found among F set" << endl;
                    // }

                    //===> step (3)

                    // unordered_set<Label, LabelHashFunction> dominatedNewLabelSet;
                    // unordered_set<Label, LabelHashFunction> dominatedExistingLabelSet;
                    set<Label> dominatedNewLabelSet;
                    set<Label> dominatedExistingLabelSet;

                    for (auto &newLabel : F_nodeSuccLabelSet)
                    {

                        // we can compare dist and cost first, if  not (d1>=d2, c1>=c2, or d1<=d2, c1<=c2), continue

                        // then, use some other methods to decide, like if any element violates dominated conditin, then break, so we don't need to check all the elements in a label.

                        for (auto &oldLabel : bigLambda_allNodesLableSet[succ])
                        {
                            if (PRINT_FOR_DEBUG)
                            {
                                cout << "new label: ";
                                printOneLabel(newLabel, bigM);
                                cout << "old label: ";
                                printOneLabel(oldLabel, bigM);
                            }

                            // if (newLabel.BVN != oldLabel.BVN)
                            // {
                            //     if (!(newLabel.distance <= oldLabel.distance && newLabel.cost <= oldLabel.cost))     // old label is not dominated
                            //         if (!(newLabel.distance >= oldLabel.distance && newLabel.cost >= oldLabel.cost)) // new label is not dominated
                            //             continue;
                            //     // now, at least one label is dominated
                            //     long result = (newLabel.BVN | oldLabel.BVN);
                            //     if (result == newLabel.BVN)
                            //         dominatedNewLabelSet.insert(newLabel);
                            //     else if (result == oldLabel.BVN)
                            //         dominatedExistingLabelSet.insert(oldLabel);
                            //     if (PRINT_FOR_DEBUG)
                            //     {
                            //         if (result == newLabel.BVN)
                            //             cout << "New label is dominated by old label." << endl;
                            //         else if (result == oldLabel.BVN)
                            //             cout << "Old label is dominated." << endl;
                            //     }
                            // }
                            // else
                            //     dominatedNewLabelSet.insert(newLabel); // we need to add the same lable in dominated new labels to to remove it from F set later, so it won't be in the newly added labels in next iteration

                            if (newLabel.distance >= oldLabel.distance && newLabel.cost >= oldLabel.cost)
                            {
                                long result = (newLabel.BVN | oldLabel.BVN);
                                if (result == newLabel.BVN)
                                {
                                    dominatedNewLabelSet.insert(newLabel);

                                    if (PRINT_FOR_DEBUG)
                                        cout << "new label is dominated" << endl;
                                }
                            } // in case that old and new labels are the same, then only add in dominatedNewLabelSet
                            else if (newLabel.distance <= oldLabel.distance && newLabel.cost <= oldLabel.cost) // old label is possible to be dominated
                            {
                                long result = (newLabel.BVN | oldLabel.BVN);
                                if (result == oldLabel.BVN)
                                {
                                    dominatedExistingLabelSet.insert(oldLabel);

                                    if (PRINT_FOR_DEBUG)
                                        cout << "new label is dominated" << endl;
                                }
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
                        cout << "-> F_nodeSuccLabelSet" << endl;
                        printLabelSet(F_nodeSuccLabelSet, bigM);

                        cout << "# of labels at succ node: " << bigLambda_allNodesLableSet[succ].size() << endl;
                        cout << "-> labels after update succ=" << succ << " labels" << endl;
                        printLabelSet(bigLambda_allNodesLableSet[succ], bigM);
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
        printf("\nthe found visited nodes in label\n");
        printLabelSet(bigLambda_allNodesLableSet[n - 1], bigM);
    }

    double minCost = numeric_limits<double>::max();
    vector<int> bestLabel;
    int numAllVistedNodes;

    for (auto &labelTemp : bigLambda_allNodesLableSet[n - 1])
    {
        if (labelTemp.cost < minCost)
        {
            bestLabel.clear();
            minCost = labelTemp.cost;
            numAllVistedNodes = labelTemp.numVisitedNodes;

            for (auto &ele : labelTemp.visitedNodes)
                bestLabel.push_back(ele);
        }
    }

    if (PRINT_FOR_DEBUG)
    {
        cout << "best label:" << endl;
        for (auto &ele : bestLabel)
            cout << ele << " ";
        cout << endl;
        cout << "min cost: " << minCost;
    }

    *obj = minCost;

    for (int i = 0; i < n; i++)
        nodesOnRoute.push_back(0);

    //===> for example: label=[1, bigM, 2, 3, 3, 12.8, -1.28], n=4
    //===> it means, node 0 is visted first, then node 2, then node 3: 0->2->3

    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        if (bestLabel[i] < bigM - verySmallNum && bestLabel[i] > verySmallNum)
        {
            int index = (long)bestLabel[i];
            nodesOnRoute[index - 1] = i;
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

#endif