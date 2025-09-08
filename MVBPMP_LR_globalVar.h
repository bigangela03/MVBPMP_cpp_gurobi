#ifndef MVBPMP_LR_GLOBALVAR_H
#define MVBPMP_LR_GLOBALVAR_H

// we assume than if the next best profit is within range below, it is the same as the best one
// for example, if p1=1.999999, p2=1.999998, then by scaling up by a factor of 2500
// p1=1.999999*2500=$4999.9975, p2=999998*2500=$4999.995, the diff is less than 1 cent
// so we assume these two profits are the same
extern double pftTolerance;

extern double LR_gap_tolerance;
extern double LR_complementarity_tolerance;
// extern double LR_min_lamda;//remove this stopping criteria, only keep miu tolerance
extern double LR_miu_tolerance;
extern double LR_rou;

extern double LR_lamda;
extern int LR_maxNumNoImprovement;

extern double bigM;
extern int TIME_LIMIT;
extern double LBinGRB;
extern double UBinGRB;
extern double LBzero;
extern double LBinLR;
extern double UBinLR;

#endif