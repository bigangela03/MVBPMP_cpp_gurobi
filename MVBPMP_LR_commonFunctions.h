#ifndef MVBPMP_LR_COMMONFUNCTIONS_H
#define MVBPMP_LR_COMMONFUNCTIONS_H

double updateUB(double, double, double);
void printBestLB(double, double, double);
void printBestUB(double, double, double);

string
itos(int i)
{
  stringstream s;
  s << i;
  return s.str();
}

double updateUB(double UB, double UBinLR, double UBinGRB)
{
  double bestValue;
  if (UB > UBinGRB && UBinLR > UBinGRB)
  {
    bestValue = UBinGRB;
    printf("===> UBinGRB = %lf is a better upper bound\n", UBinGRB);
  }
  else if (UB > UBinLR && UBinGRB > UBinLR)
  {
    bestValue = UBinLR;
    printf("===> UBinLR = %lf is a better upper bound\n", UBinLR);
  }
  else if (UB > UBinLR && UBinGRB > UBinLR - 0.00000001 && UBinGRB < UBinLR + 0.00000001)
  {
    bestValue = UBinGRB;
    printf("===> both UBinLR = %lf and UBinGRB = %lf are better upper bounds\n", UBinLR, UBinGRB);
  }
  else
  {
    cout << "Did not update UB." << endl;
    printf("UB=%lf, UBinLR=%lf, UBinGRB=%lf\n", UB, UBinLR, UBinGRB);
    bestValue = UB;
  }

  return bestValue;
}

void printBestLB(double LB, double LBinLR, double LBinGRB)
{
  if (LBinLR == LB && LBinGRB < LB)
    // cout << "Best LB is found by LR" << endl;
    cout << "BestLowerBoundSolver = LR" << endl;
  else if (LBinLR < LB && LBinGRB == LB)
    // cout << "Best LB is found by GRB in parallel" << endl;
    cout << "BestLowerBoundSolver = GRB" << endl;
  else if (LBinLR == LB && LBinGRB == LB)
    // cout << "Best LB is found by both LR and GRB in parallel" << endl;
    cout << "BestLowerBoundSolver = LR_GRB" << endl;
  else if (LBinLR > LB - 0.00000001 && LBinLR < LB + 0.00000001 && LBinGRB > LB - 0.00000001 && LBinGRB < LB + 0.00000001)
    cout << "BestLowerBoundSolver = LR_GRB" << endl;
  else
    printf("ERROR: LB=%lf, LBinLR=%lf, LBinGRB=%lf\n", LB, LBinLR, LBinGRB);
}

void printBestUB(double UB, double UBinLR, double UBinGRB)
{
  if (UBinLR == UB && UBinGRB > UB)
    // cout << "Best UB is found by LR" << endl;
    cout << "BestUpperBoundSolver = LR" << endl;
  else if (UBinLR > UB && UBinGRB == UB)
    // cout << "Best UB is found by GRB in parallel" << endl;
    cout << "BestUpperBoundSolver = GRB" << endl;
  else if (UBinLR == UB && UBinGRB == UB)
    // cout << "Best UB is found by both LR and GRB in parallel" << endl;
    cout << "BestUpperBoundSolver = LR_GRB" << endl;
  else if (UBinLR > UB - 0.00000001 && UBinLR < UB + 0.00000001 && UBinGRB > UB - 0.00000001 && UBinGRB < UB + 0.00000001)
    cout << "BestUpperBoundSolver = LR_GRB" << endl;
  else
    printf("ERROR: UB=%lf, UBinLR=%lf, UBinGRB=%lf\n", UB, UBinLR, UBinGRB);
}

// I tried moving solving MVBPMP here so that same start node and diff start node can share the code
// but then I found too many variabls have to be passed to here as arguments
// so I gave up.

// class printIntSol : public GRBCallback
// {
// public:
//   int n;
//   int numV;
//   GRBVar ***xv;
//   GRBVar ***yv;
//   GRBVar ****uv;
//   GRBVar **sv;
//   double *UBinGRB;
//   double *LBinGRB;
//   double ***solx_GRB;
//   double ***soly_GRB;

//   printIntSol(GRBVar ***xvars, GRBVar ***yvars, GRBVar ****uvars, GRBVar **svars, int nvar, int numVvar,
//               double *UBinGRBarg, double *LBinGRBarg, double ***solx_GRBarg, double ***soly_GRBarg)
//   {
//     xv = xvars;
//     yv = yvars;
//     uv = uvars;
//     sv = svars;
//     n = nvar;
//     numV = numVvar;
//     UBinGRB = UBinGRBarg;
//     LBinGRB = LBinGRBarg;
//     solx_GRB = solx_GRBarg;
//     soly_GRB = soly_GRBarg;
//   }

// protected:
//   void
//   callback()
//   {
//     try
//     {
//       if (where == GRB_CB_MIPSOL)
//       {
//         // Found an integer feasible solution
//         double objTemp = getDoubleInfo(GRB_CB_MIPSOL_OBJ);

//         if (objTemp > *LBinGRB)
//           *LBinGRB = objTemp;

//         printf("--------> GRB: OBJ = %lf.\n", objTemp);

//         *UBinGRB = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);
//         printf("--------> GRB:  UB = %lf.\n", *UBinGRB);

//         int i, j, q;

//         for (i = 0; i < n; i++)
//         {
//           // sols_GRB[i] = getSolution(sv[i], numV);
//           for (j = 0; j < n; j++)
//           {
//             solx_GRB[i][j] = getSolution(xv[i][j], numV);
//             soly_GRB[i][j] = getSolution(yv[i][j], numV);
//             // for (int k = 0; k < n; k++)
//             //   solu_GRB[i][j][k] = getSolution(uv[i][j][k], numV);

//             for (q = 0; q < numV; q++)
//             {
//               // if (PRINT_x_WHEN_INTEGER_SOL)
//               if (solx_GRB[i][j][q] > 0.5)
//                 printf("x: %3d ->%3d (v%d)\n", i + 1, j + 1,
//                        q + 1);
//               // if (PRINT_y_WHEN_INTEGER_SOL)
//               if (soly_GRB[i][j][q] > 0.5)
//                 printf("y: %3d ->%3d (v%d)\n", i + 1, j + 1,
//                        q + 1);
//             }
//           }
//         }
//       }
//     }
//     catch (GRBException e)
//     {
//       cout << "Error number: " << e.getErrorCode() << endl;
//       cout << e.getMessage() << endl;
//     }
//     catch (...)
//     {
//       cout << "Error during callback" << endl;
//     }
//   }
// };

// void solveMVBPMPinParallel(double UB, double UBinLR, double UBinGRB, double LB, double LBinLR, double LBinGRB,
//                            double ***solx_GRB, double ***soly_GRB,
//                            vector<vector<vector<int>>> &allVehiclesInaccNeighbors,
//                            vector<vector<vector<int>>> &allVehiclesNodeNeighbors)
// {

//   printf("===> solving MVBPMP with GRB <===\n");

//   int i, j, k, q;
//   int status, nSolutions;

//   GRBEnv *env = NULL;
//   // GRBVar x[n][n][numV];
//   // GRBVar y[n][n][numV];
//   // GRBVar s[n][numV];
//   // GRBVar u[n][n][n][numV];
//   GRBVar theta[n][n][numV];

//   GRBVar ***x = NULL;
//   GRBVar ***y = NULL;
//   x = new GRBVar **[n];
//   y = new GRBVar **[n];
//   for (i = 0; i < n; i++)
//   {
//     x[i] = new GRBVar *[n];
//     y[i] = new GRBVar *[n];
//     for (j = 0; j < n; j++)
//     {
//       x[i][j] = new GRBVar[numV];
//       y[i][j] = new GRBVar[numV];
//     }
//   }

//   GRBVar ****u = new GRBVar ***[n];
//   GRBVar **s = new GRBVar *[n];

//   for (int i = 0; i < n; i++)
//   {
//     u[i] = new GRBVar **[n];
//     s[i] = new GRBVar[numV];
//     for (int j = 0; j < n; j++)
//     {
//       u[i][j] = new GRBVar *[n];
//       for (int k = 0; k < n; k++)
//         u[i][j][k] = new GRBVar[numV];
//     }
//   }

//   try
//   {

//     env = new GRBEnv();
//     GRBModel model = GRBModel(*env);

//     // Create binary decision variables
//     for (q = 0; q < numV; q++)
//     {
//       for (i = 0; i < n; i++)
//       {
//         s[i][q] = model.addVar(0.0, n, 0.0, GRB_CONTINUOUS,
//                                "s_" + itos(i) + "_" + itos(q));
//         for (j = 0; j < n; j++)
//         {
//           x[i][j][q] = model.addVar(
//               0.0, 1.0, 0, GRB_BINARY,
//               "x_" + itos(i) + "_" + itos(j) + "_" + itos(q));
//           y[i][j][q] = model.addVar(
//               0.0, 1.0, 0, GRB_BINARY,
//               "y_" + itos(i) + "_" + itos(j) + "_" + itos(q));
//           theta[i][j][q] = model.addVar(
//               0.0,
//               GRB_INFINITY,
//               0.0,
//               GRB_CONTINUOUS,
//               "theta_" + itos(i) + "_" + itos(j) + "_" + itos(q));

//           for (k = 0; k < n; k++)
//           {
//             string s = "u_" + itos(i) + "_" + itos(j) + "_" + itos(k) + "_" + itos(q);
//             u[i][j][k][q] = model.addVar(0.0, GRB_INFINITY, 0.0,
//                                          GRB_CONTINUOUS, s);
//           }
//         }
//       }

//       int ogn = origin[q];
//       for (i = 0; i < n; i++)
//       {
//         x[i][i][q].set(GRB_DoubleAttr_UB, 0);
//         y[i][i][q].set(GRB_DoubleAttr_UB, 0);
//         x[i][ogn][q].set(GRB_DoubleAttr_UB, 0);
//         y[i][ogn][q].set(GRB_DoubleAttr_UB, 0);
//         x[n - 1][i][q].set(GRB_DoubleAttr_UB, 0);
//         y[n - 1][i][q].set(GRB_DoubleAttr_UB, 0);
//         for (j = 0; j < n; j++)
//           if (wt[i][j] == 0)
//             y[i][j][q].set(GRB_DoubleAttr_UB, 0);
//       }

//       // if (ADD_PREPROCESS)
//       {
//         for (int i = 0; i < n; i++)
//         {
//           vector<int> nbsTemp = allVehiclesInaccNeighbors[q][i];
//           for (j = 0; j < nbsTemp.size(); j++)
//           {
//             x[i][nbsTemp[j]][q].set(GRB_DoubleAttr_UB, 0);
//             // y[i][nbsTemp[j]].set(GRB_DoubleAttr_UB, 0);
//             for (k = 0; k < n; k++)
//             {
//               u[i][nbsTemp[j]][k][q].set(GRB_DoubleAttr_UB, 0);
//               u[i][k][nbsTemp[j]][q].set(GRB_DoubleAttr_UB, 0);
//               u[k][nbsTemp[j]][i][q].set(GRB_DoubleAttr_UB, 0);
//             }
//           }
//         }

//         // int numViolatedTriples = 0;
//         for (int i = 0; i < n; i++)
//           if (i != ogn && i != endNode)
//             for (auto &k : allVehiclesNodeNeighbors[q][i])
//               for (auto &j : allVehiclesNodeNeighbors[q][k])
//                 if (i != j)
//                   if (dis[ogn][i] + dis[i][k] + dis[k][j] + dis[j][endNode] > disLimit)
//                   {
//                     u[i][j][k][q].set(GRB_DoubleAttr_UB, 0);
//                     // numViolatedTriples++;
//                   }
//         // cout << "numViolatedTriples=" << numViolatedTriples << endl;
//       }
//     }

//     // set up constraints
//     GRBConstr *vehOriginConstr = 0;
//     GRBConstr *vehDestConstr = 0;
//     vehOriginConstr = new GRBConstr[numV];
//     vehDestConstr = new GRBConstr[numV];

//     for (q = 0; q < numV; q++)
//     {
//       int ogn = origin[q];

//       // vehicle goes out of origins
//       GRBLinExpr expr1 = 0.0;
//       for (i = 0; i < n; i++)
//         expr1 += x[ogn][i][q];

//       // model.addConstr (expr1 == 1, "origin_" + itos (q));
//       vehOriginConstr[q] = model.addConstr(expr1 == 1,
//                                            "origin_" + itos(q));

//       // vehicle goes back to node n
//       GRBLinExpr expr2 = 0.0;
//       for (i = 0; i < n - 1; i++)
//         expr2 += x[i][n - 1][q];
//       // model.addConstr (expr2 == 1, "destination_" + itos (q));
//       vehDestConstr[q] = model.addConstr(expr2 == 1,
//                                          "destination_" + itos(q));

//       // flow conservation
//       for (int k = 0; k < n - 1; k++)
//       {
//         if (k != ogn)
//         {
//           GRBLinExpr expr = 0;
//           for (i = 0; i < n - 1; i++)
//             expr += x[i][k][q];

//           // BE CAREFUL!
//           // I used j=1 to start which exclues node 0
//           // which cause that the optimal profit is lower!
//           for (j = 0; j < n; j++)
//             expr -= x[k][j][q];
//           model.addConstr(
//               expr == 0,
//               "flow_conservation_" + itos(k) + "_" + itos(q));
//         }
//       }

//       // distance
//       GRBLinExpr expr3 = 0.0;
//       for (i = 0; i < n - 1; i++)
//         for (j = 0; j < n; j++)
//           expr3 += dis[i][j] * x[i][j][q];
//       model.addConstr(expr3 <= route.DIS, "distance_" + itos(q));

//       // node degree less than 1
//       for (int j = 0; j < n - 1; j++)
//       {
//         GRBLinExpr expr = 0.0;
//         for (i = 0; i < n - 1; i++)
//           expr += x[i][j][q];
//         model.addConstr(expr <= 1,
//                         "indegree_" + itos(j) + "_" + itos(q));
//       }

//       // subtour elimination
//       for (i = 0; i < n - 1; i++)
//         for (j = 0; j < n; j++)
//         {
//           GRBLinExpr expr = 0.0;
//           expr += s[i][q] - s[j][q] + (n - 1) * x[i][j][q] + (n - 3) * x[j][i][q];
//           model.addConstr(
//               expr <= n - 2,
//               "subtour_" + itos(i) + "_" + itos(j) + "_" + itos(q));
//         }

//       // arc flow
//       for (i = 0; i < n - 1; i++)
//         for (j = 0; j < n; j++)
//         {
//           GRBLinExpr expr = 0.0;

//           expr += wt[i][j] * y[i][j][q] - theta[i][j][q];
//           for (k = 0; k < n - 1; k++)
//           {
//             if (k != ogn)
//               expr += u[i][k][j][q] + u[k][j][i][q] - u[i][j][k][q];
//             else
//               expr += u[k][j][i][q];
//           }
//           // when k==n-1
//           if (j != n - 1)
//             expr += u[i][n - 1][j][q];
//           model.addConstr(
//               expr == 0,
//               "flow_" + itos(i) + "_" + itos(j) + "_" + itos(q));

//           /*
//            expr += wt[i][j] * y[i][j][q] - theta[i][j][q];
//            for (k = 0; k < n; k++)
//            expr += u[i][k][j][q] + u[k][j][i][q] - u[i][j][k][q];
//            */
//         }

//       // arc flow upperbound
//       for (i = 0; i < n - 1; i++)
//         for (j = 0; j < n; j++)
//         {
//           GRBLinExpr expr = 0.0;
//           expr += theta[i][j][q] - Q * x[i][j][q];
//           model.addConstr(
//               expr <= 0,
//               "flowBound_" + itos(i) + "_" + itos(j) + "_" + itos(q));
//         }

//       // add redundant constraint to test if same multi obj value
//       // caused by binary variable
//       /*
//        for (i = 0; i < n - 1; i++)
//        for (j = 0; j < n; j++)
//        for (k = 0; k < n - 1; k++)
//        {
//        GRBLinExpr expr = 0.0;
//        if (k != ogn)
//        expr += u[i][j][k][q] - x[i][k][q];
//        model.addConstr (
//        expr <= 0,
//        "unique_triples_" + itos (i) + "_" + itos (j)
//        + "_" + itos (q));
//        }
//        */
//     }

//     // one vehicle for one cargo
//     // if (!ADD_CUTS){}
//     for (i = 0; i < n - 1; i++)
//       for (j = 0; j < n; j++)
//       {
//         GRBLinExpr expr = 0.0;
//         for (q = 0; q < numV; q++)
//           expr += y[i][j][q];
//         model.addConstr(expr <= 1,
//                         "one-one" + itos(i) + "_" + itos(j));
//       }

//     // set up objective
//     GRBLinExpr obj = 0.0;
//     for (q = 0; q < numV; q++)
//     {
//       for (i = 0; i < n - 1; i++)
//         for (j = 0; j < n; j++)
//         {
//           obj += price * dis[i][j] * wt[i][j] * y[i][j][q];
//           obj -= cost * dis[i][j] * theta[i][j][q];
//           obj -= cost * vw * dis[i][j] * x[i][j][q];
//         }
//     }

//     model.setObjective(obj, GRB_MAXIMIZE);

//     // Set callback function
//     printIntSol cb = printIntSol(x, y, u, s, n, numV, &UBinGRB, &LBinGRB, solx_GRB, soly_GRB);
//     model.setCallback(&cb);

//     model.set(GRB_IntParam_OutputFlag, 0);

//     model.set(GRB_IntParam_Threads, NUM_THREADS_VEH);

//     // Optimize model
//     model.optimize();

//     // write model to file
//     // model.write ("MVBPMP.lp");

//     // Status checking
//     status = model.get(GRB_IntAttr_Status);
//     if (status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE || status == GRB_UNBOUNDED)
//     {
//       cout << "The model cannot be solved "
//            << "because it is infeasible or unbounded" << endl;
//       exit(1);
//     }
//     if (status != GRB_OPTIMAL)
//     {
//       cout << "Optimization was stopped with status " << status << endl;
//       exit(1);
//     }
//     if (status == GRB_OPTIMAL)
//     {
//       // Gurobi found optimal solution before LR converges
//       // print out optimal solution and stop
//       cout << "Optimization was completed by Gurobi. Optimal solution found. " << endl;

//       if (model.get(GRB_IntAttr_SolCount) > 0)
//       {

//         double runtime = model.get(GRB_DoubleAttr_Runtime);
//         cout << "Gurobi in parallel Runtime: " << runtime << " seconds" << endl;

//         double objtemp = model.get(GRB_DoubleAttr_ObjVal);
//         UB = objtemp;
//         LB = objtemp;

//         //====== start defining solutions for x, y, u, and s ======
//         double ***solx = new double **[n]; // solution x for LR dual
//         double ***soly = new double **[n];
//         double ****solu = new double ***[n];
//         double **sols = new double *[n];

//         for (int i = 0; i < n; i++)
//         {
//           solx[i] = new double *[n];
//           soly[i] = new double *[n];
//           solu[i] = new double **[n];
//           sols[i] = new double[numV];

//           for (int j = 0; j < n; j++)
//           {
//             solx[i][j] = new double[numV];
//             soly[i][j] = new double[numV];
//             solu[i][j] = new double *[n];

//             for (int k = 0; k < n; k++)
//               solu[i][j][k] = new double[numV];
//           }
//         }

//         //====== start reading solutions for x, y, u, and s ======
//         // solx[i][j] = model.get (GRB_DoubleAttr_X, x[i][j],numV);
//         // soly[i][j] = model.get (GRB_DoubleAttr_X, y[i][j],numV);
//         for (q = 0; q < numV; q++)
//           for (i = 0; i < n; i++)
//           {
//             sols[i][q] = s[i][q].get(GRB_DoubleAttr_X);
//             for (j = 0; j < n; j++)
//             {
//               solx[i][j][q] = x[i][j][q].get(GRB_DoubleAttr_X);
//               soly[i][j][q] = y[i][j][q].get(GRB_DoubleAttr_X);
//               for (k = 0; k < n; k++)
//                 solu[i][j][k][q] = u[i][j][k][q].get(
//                     GRB_DoubleAttr_X);
//             }
//           }

//         // store the solution as the best LB
//         // storeBestLB(solx, soly, solu, sols, solx_best, soly_best,solu_best, sols_best, n);
//         storeBestLB(solx, soly, solx_best, soly_best, n);

//         for (i = 0; i < n; i++)
//         {
//           for (j = 0; j < n; j++)
//           {
//             for (k = 0; k < n; k++)
//               delete[] solu[i][j][k];

//             delete[] solx[i][j];
//             delete[] soly[i][j];
//             delete[] solu[i][j];
//           }
//           delete[] solx[i];
//           delete[] soly[i];
//           delete[] solu[i];
//           delete[] sols[i];
//         }
//         delete[] solx;
//         delete[] soly;
//         delete[] solu;
//         delete[] sols;
//       }
//     }
//   }
//   catch (GRBException e)
//   {
//     cout << "Error number: " << e.getErrorCode() << endl;
//     cout << e.getMessage() << endl;
//   }
//   catch (...)
//   {
//     cout << "Error during optimization" << endl;
//   }

//   delete env;

//   cout << "**************** THE OPTIMAL SOLUTION ****************" << endl;
//   printVar(solx_best, soly_best, solu_best, origin);
//   cout << endl
//        << "===> The total time: " << endl;
//   reportTime(beginTime, beginWallClock);

//   clock_t end = clock();
//   double second = (double)(end - beginTime) / CLOCKS_PER_SEC;

//   auto endWallClock = high_resolution_clock::now();
//   auto elapsedWallClock = duration_cast<std::chrono::nanoseconds>(endWallClock - beginWallClock);

//   cout << "===> SUMMARY:" << endl;
//   printf("SolutionTimeCPU = %lf\n", second);
//   printf("SolutionTimeWallClock = %.3f seconds\n",
//          elapsedWallClock.count() * 1e-9);
//   printf("SolutionLB = %lf\n", LB);
//   printf("SolutionUB = %lf\n", UB);
//   printf("SolutionGap = %lf\n", (UB - LB) / LB);
//   // if program finish in this thread (thread 1), it means Gurobi finishes and find optimal solution before LR finishes
//   // cout << "Best LB is found by GRB in parallel (find optimal solution)" << endl;
//   cout << "BestLowerBoundSolver = GRB" << endl;
//   cout << "BestUpperBoundSolver = GRB" << endl;
//   exit(1);
//   // return 0;
// }

#endif