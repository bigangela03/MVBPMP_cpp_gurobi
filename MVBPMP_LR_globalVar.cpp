#include "MVBPMP_LR_globalVar.h"

// we assume than if the next best profit is within range below, it is the same as the best one
// for example, if p1=1.999999, p2=1.999998, then by scaling up by a factor of 2500
// p1=1.999999*2500=$4999.9975, p2=999998*2500=$4999.995, the diff is less than 1 cent
// so we assume these two profits are the same
double pftTolerance = 0.000001;

double LR_gap_tolerance = 0.05;
// double LR_complementarity_tolerance = 0.0001;
double LR_complementarity_tolerance = pftTolerance;

// double LR_min_lamda = 0.0001; //remove this stopping criteria, only keep miu tolerance
double LR_miu_tolerance = 0.00000001;

double LR_rou = 0.4; // for LR multiplier type b

double LR_lamda = 2;
int LR_maxNumNoImprovement = 3;

double bigM = 10000000;

double LBinGRB = -bigM;
double UBinGRB = bigM;
double LBzero = 0.0000000001;
double LBinLR = -bigM;
double UBinLR = bigM;

int TIME_LIMIT = 3600;
