#ifndef DOMINACE_H
#define DOMINACE_H

// This code is only for complete graph right now, i.e., any two nodes are connected
// since I didn't consider the vertex incidence matrix, only iterate nodes from 0 to n-1

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> //to get path to current directory
#include <math.h>
#include <time.h>
#include <malloc.h>
#include <string.h>

#include "gurobi_c++.h"
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

#include <limits> //to get the biggest value of double variables

vector<int>
runDominace(int n, double **dis, double **xCoeff, double disLimit, double *obj)
{
    vector<int> selectedRoute;

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

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cout << xCoeff[i][j] << " ";
        cout << endl;
    }

    int startingNode = 0;

    vector<set<vector<double>>> bigLambda_allNodesLableSet;

    //===> initialize the lable for starting node
    vector<double> oneLabel;
    for (int i = 0; i < n + 3; i++)
        oneLabel.push_back(0);
    oneLabel[startingNode] = 1;
    oneLabel[n] = 1;

    set<vector<double>> startingNodeLables;
    startingNodeLables.insert(oneLabel);

    //===> add empty lable set to all other nodes
    for (int i = 0; i < n; i++)
    {
        if (i == startingNode)
            bigLambda_allNodesLableSet.push_back(startingNodeLables);
        else
        {
            set<vector<double>> emptySet;
            bigLambda_allNodesLableSet.push_back(emptySet);
        }
    }

    // for (auto &vec : bigLambda_allNodesLableSet[startingNode])
    // {
    //     for (auto &num : vec)
    //         cout << num << " " << endl;
    //     cout << endl;
    // }

    set<int> E_toBeTreatedNodesSet;
    E_toBeTreatedNodesSet.insert(startingNode);

    int itrNum = 0; // iteration number

    while (E_toBeTreatedNodesSet.size() > 0) //===> Line 7 in ESPPRC(p)
    {
        itrNum++;

        cout << endl;
        cout << "======> dominace iteration " << itrNum << " <======" << endl;
        cout << "E_toBeTreatedNodesSet size = " << E_toBeTreatedNodesSet.size() << endl;
        for (auto &temp : E_toBeTreatedNodesSet)
            cout << temp << " ";
        cout << endl;

        set<int> E_toBeTreatedNodesSetCopy;

        for (auto &node : E_toBeTreatedNodesSet)
            E_toBeTreatedNodesSetCopy.insert(node);

        for (auto &node : E_toBeTreatedNodesSetCopy) //===> LINE 9 in ESPPRC(p)
        {
            cout << "===> treaing node " << node << endl;
            for (int succ = 0; succ < n; succ++) //===> LINE 10 in ESPPRC(p)
            {
                if (succ == startingNode || succ == node)
                    continue;
                cout << endl;
                cout << "succ=" << succ << endl;
                set<vector<double>> F_nodeSuccLabelSet; //===> LINE 11 in ESPPRC(p) PSEUDO-CODE

                set<vector<double>> labelsOverDistLimit;

                cout << "bigLambda_allNodesLableSet[" << node << "].size =" << bigLambda_allNodesLableSet[node].size() << endl;
                for (auto &label : bigLambda_allNodesLableSet[node]) //===> LINE 12 in ESPPRC(p) PSEUDO-CODE
                {
                    cout << "-> ";
                    for (auto &ele : label)
                        cout << ele << " ";
                    cout << endl;
                    if (label[succ] == 0) //===> LINE 13 in ESPPRC(p) PSEUDO-CODE
                    {
                        cout << "find a label not visit succ yet" << endl;
                        //===> LINE 14 in ESPPRC(p) PSEUDO-CODE

                        //===> check resource (distance)
                        double currentDistance = label[n + 1];

                        vector<double> labelTemp;
                        for (int i = 0; i < n + 3; i++)
                            labelTemp.push_back(label[i]);

                        if (currentDistance + dis[node][succ] <= disLimit)
                        {
                            labelTemp[succ] = 1;                                  //!!! lable[0] to lable[n-1] represents visited nodes
                            labelTemp[n] = label[n] + 1;                          //!!! lable[n] is the number of visited nodes;
                            labelTemp[n + 1] = currentDistance + dis[node][succ]; //!!! label[n+1] is the consumed resouces;
                            labelTemp[n + 2] = label[n + 2] + xCoeff[node][succ]; //!!! lable[n+2] is the cost
                            F_nodeSuccLabelSet.insert(labelTemp);

                            cout << " visited number of nodes: labelTemp[n]=" << labelTemp[n] << endl;
                            cout << "traveled distance: labelTemp[n + 1]=" << labelTemp[n + 1] << endl;
                        }
                        else
                        {
                            cout << "distance is over limit" << endl;
                            labelsOverDistLimit.insert(label);
                        }
                    }
                }

                for (auto &labelTemp : labelsOverDistLimit)
                {
                    auto it = bigLambda_allNodesLableSet[node].find(labelTemp);

                    if (it != bigLambda_allNodesLableSet[node].end())
                    {
                        vector<double> temp = *it;
                        temp[succ] = 1;
                        bigLambda_allNodesLableSet[node].erase(it);
                        bigLambda_allNodesLableSet[node].insert(temp);
                    }
                    else
                    {
                        cout << "Error: labels that exceed distance limit is not found." << endl;
                        exit(1);
                    }
                }

                //========> EFF function in LINE 15 and 16 in ESPPRC(p) PSEUDO-CODE <========
                {
                    //===> check if each label in F_nodeSuccLabelSet is dominated by bigLambda_allNodesLableSet[succ]
                    // if each label is dominated, then no need to add to bigLambda_allNodesLableSet[succ], so no change, and hasChange = false;
                    // if not, then the non-dominated label is added to bigLambda_allNodesLableSet[succ], so hasChange = true;

                    set<vector<double>> dominatedNewLabelSet;
                    set<vector<double>> dominatedExistingLabelSet;

                    //===> find the dominated new labels
                    for (auto &newLabel : F_nodeSuccLabelSet)
                    {
                        // bool newLabelIsDominated = true;
                        for (auto &oldLabel : bigLambda_allNodesLableSet[succ])
                        {
                            int numBiggerEleInOldLabel = 0;
                            int numBiggerEleInNewLabel = 0;
                            int numEqual = 0;
                            for (int i = 0; i < n + 3; i++)
                            {
                                if (newLabel[i] <= oldLabel[i])
                                    numBiggerEleInOldLabel++;

                                if (newLabel[i] == oldLabel[i])
                                    numEqual++;

                                if (newLabel[i] >= oldLabel[i])
                                    numBiggerEleInNewLabel++;
                            }

                            if (numBiggerEleInNewLabel == (n + 3) && numEqual != (n + 3))
                                dominatedNewLabelSet.insert(newLabel);

                            if (numBiggerEleInOldLabel == (n + 3) && numEqual != (n + 3))
                                dominatedExistingLabelSet.insert(oldLabel);

                            // if(numEqual==(n+3)), it means oldLabel is the same as newLabel
                            // when we add newLabel into old label set, the set only keep one
                            // so we don't need to add this label to either toBeRemoved* label set.
                            // if we don't consider numEqual != (n+3) in the above 2 condistions, then
                            // this equal label will be added to both toBeRemoved set and finally removed from label set.
                        }
                    }

                    cout << "-> succ node labels" << endl;
                    for (auto &labelTemp : bigLambda_allNodesLableSet[succ])
                        for (auto &temp : labelTemp)
                            cout << temp << " ";
                    cout << endl;

                    cout << "-> dominatedExistingLabelSet" << endl;
                    for (auto &labelTemp : dominatedExistingLabelSet)
                        for (auto &temp : labelTemp)
                            cout << temp << " ";
                    cout << endl;

                    cout << "-> dominatedNewLabelSet" << endl;
                    for (auto &labelTemp : dominatedNewLabelSet)
                        for (auto &temp : labelTemp)
                            cout << temp << " ";
                    cout << endl;

                    //===> remove the dominated labels for current and succsor nodes
                    for (auto &labelTemp : dominatedExistingLabelSet)
                        bigLambda_allNodesLableSet[succ].erase(labelTemp);
                    for (auto &labelTemp : dominatedNewLabelSet)
                        F_nodeSuccLabelSet.erase(labelTemp);

                    //===> add the non-dominated labels to bigLambda_allNodesLableSet[succ]
                    for (auto &labelTemp : F_nodeSuccLabelSet)
                        bigLambda_allNodesLableSet[succ].insert(labelTemp);

                    if (dominatedExistingLabelSet.size() != 0 || F_nodeSuccLabelSet.size() != 0) //===> LINE 16 in ESPPRC(p) PSEUDO-CODE
                        E_toBeTreatedNodesSet.insert(succ);                                      //===> LINE 17 in ESPPRC(p) PSEUDO-CODE
                }
            }

            E_toBeTreatedNodesSet.erase(node); //===> LINE 19 in ESPPRC(p) PSEUDO-CODE
        }
    }

    cout << endl;

    cout << "the found visited nodes in label" << endl;
    double minCost = numeric_limits<double>::max();
    for (auto &labelTemp : bigLambda_allNodesLableSet[n - 1])
    {
        for (auto &ele : labelTemp)
        {
            cout << ele << " ";
            if (labelTemp[n + 2] < minCost)
                minCost = labelTemp[n + 2];
        }
        cout << endl;
    }

    *obj = minCost;

    // work here, get selected route from labels!!!!!!!!!!!!!!!!!!!
    return selectedRoute;
}

#endif