#ifndef GA_OBJFUNC_GREEDY_H
#define GA_OBJFUNC_GREEDY_H

// GA2 only use greedy algorithm when calculating profit
/*the distance should be less than or equal to maxdis, but I wrote it as distanc<maxdis which lead to bad solution. May25, 2011*/
/*delete A in this program. if distance>limt, set fitness=0*/
/*no extra test printing compared with GA30A.c*/
/*in function(), change the condition to if profit>=0, give value to critter->cust[i][j]*/
/*in function mutation() there is something wrong. it has corrected in this file*/
/*change all 8 and related number to NC*/
/*last line of objfunc() to calculate profitness if distance is more than MAXDIS should be pro not dis, the previous files are wrong*/
/*use DATA1*/
/*It's based on G10 and only change the data to the new sequence according to paper description*/
/*Line 351-363 is modified and now the processing speed is very quick*/
/*It's based on G9.c and delete some memory control sentences*/
/*sizeof(int)is four under Linux*/
/*It's based on G8.c and just keep the final result output to increase speed(it's very very slow with so many printf)*/
/*add clock() to calculate processing time*/
/*it's based on G7.c and change all the random() functions to rand()*/
/*if x=1.999, i=(int)x; then i==1; which means compulsory change always round down number value*/
/*this edition is for output every step to the file on disk*/
/*I don't know why but it always give the same solution route*/
/*it is based on G5test, just delete all the test points and change position of value of min in objfunc()*/
/*add the selected customer for every route*/

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

#include <cstdio> //to use printf

#include <iostream> //input and output
#include <fstream>	//read file
#include <sstream>
#include <algorithm>
#include <vector>
#include <climits>
#include <set>
#include <string>
#include <ctime>  //to measure CPU time
#include <chrono> //to measure run time
#include <unordered_map>

#include "readData.h"

using namespace std;
using namespace std::chrono;
using namespace std::string_literals;

#define NUM_RESTART 10
// if there is no improvement after NUM_NO_IMPRO, restart the procedure
#define pcross 0.9
#define pmutation 0.1

#define MAXRUN 10
/*5 times to tune parameter, 10 times to get average and compare with other algorithms*/
// #define MAXINSTANCE 10

#define TAU 0.7 // the probability of choose a gene from parent 1

bool UseGreedy = true;
bool UseGurobi = false;

bool UseNewCrossover = true;

bool printForDebug = false;

/*the number of all cities, including starting and ending points*/
int NumAllCities;

double P;
double C;
double Q;
/*#define A 0(1/(10000*NumAllCities)) the parameter to decide how much percent of profit kept if distance exceeds maxdis*/
/*global variable*/
double maxdis;
double wv;

// #define popsize 200
// #define maxgen 2000
int popsize; // it has to be an even number
int maxgen;
int NUM_NO_IMPRO = 50;

struct individual
{
	// int single[NumAllCities - 2];
	int *single;
	double fitness;
	double distance;
	int xsite;
	int parent[2];
	int cus[100][2]; /*the selected customer in the route*/
	int cusnum;		 /*the number of customers in the route*/
};
struct bestever
{
	// int single[NumAllCities - 2];
	int *single;
	double fitness;
	int generation;
	double distance;
	int cus[100][2]; /*the selected customer in the route*/
	int cusnum;		 /*the number of customers in the route*/
};

unordered_map<string, individual> hashmap;

struct individual *oldpop;
struct individual *newpop;
struct bestever bestfit;
double sumfitness;
double maxFitness;
double avgFitnessGlobal;
double minFitness;
int gen;
int run;
long nmutation;
long ncross;
int jcross;
// double demand[NumAllCities][NumAllCities];
// double dis[NumAllCities][NumAllCities];
double **demand;
double **dis;

int status = 0; /*the parameter for gurobi*/

// FILE *outfp;
// FILE *fp1;
// FILE *fp2;
/*definition of function*/

void copyIndividual(struct individual *, struct individual);
void moveMinusOneToEnd(struct individual *);
int flip(double);
void initialize(),
	initpop();
void initreport(),
	generation(),
	generationNewCrossover();
void report();
void preselect();
void statistics(struct individual *);
int select1();
void objfunc(struct individual *);
int crossover(int *, int *, int *, int *);
void NewCrossover(int *, int *, int *);
void mutation(int *);
void delet(int *, int *);
void initmalloc();
void freeall();

string
itos(int i)
{
	stringstream s;
	s << i;
	return s.str();
}

void copyIndividual(struct individual *indiv, struct individual chroma)
{
	int i, j;
	(*indiv).single = new int[NumAllCities - 2];
	for (i = 0; i < NumAllCities - 2; i++)
		(*indiv).single[i] = chroma.single[i];
	(*indiv).fitness = chroma.fitness;
	(*indiv).distance = chroma.distance;
	(*indiv).xsite = chroma.xsite;
	(*indiv).parent[0] = chroma.parent[0];
	(*indiv).parent[1] = chroma.parent[1];
	(*indiv).cusnum = chroma.cusnum;
	for (i = 0; i < 100; i++)
		for (j = 0; j < 2; j++)
			(*indiv).cus[i][j] = chroma.cus[i][j];
}

void moveMinusOneToEnd(struct individual *critter)
{

	int i, j, count1 = 0, count2 = 0;

	for (i = 0; i < NumAllCities - 2; i++)
	{
		if (critter->single[i] < 0)
			count1++;
	}

	for (i = 0; i < NumAllCities - 2; i++)
	{
		if (critter->single[i] < 0)
		{
			for (j = i; j < NumAllCities - 3; j++)
				critter->single[j] = critter->single[j + 1];
			i = i - 1;
			critter->single[NumAllCities - 3] = -1;
			count2++;
			if (count1 == count2)
				break;
		}
	}
}

void initialize()
{
	initmalloc();
	nmutation = 0;
	ncross = 0;
	bestfit.single = new int[NumAllCities - 2];
	bestfit.fitness = 0;
	bestfit.generation = 0;
	bestfit.distance = 0;
	bestfit.cusnum = 0;
	initpop();
	statistics(oldpop);
	/*initreport();*/
}

void initpop()
{
	int i, j, k, m, n = 1, l, count = 0;

	for (j = 0; j < popsize; j++)
	{
		oldpop[j].single = new int[NumAllCities - 2];
		// although we don't need new pop right now
		// but we have to define the size of new pop
		// so we do it here along with oldpop
		newpop[j].single = new int[NumAllCities - 2];
		for (k = 0; k < NumAllCities - 2; k++)
		{
			if (n == 1)
			{
				if (m = flip(0.8))
				{
					oldpop[j].single[k] = 1 + (int)((double)(NumAllCities - 2 - k) * rand() / (RAND_MAX + 1.0));

					/*random(8-k)+1;*/
				}
				else
				{
					i = k;
					oldpop[j].single[k] = -1;
					n = 0;
				}
			}
			if (n == 0)
			{
				for (m = i; m < NumAllCities - 2; m++)
				{
					if (m != NumAllCities - 2)
						oldpop[j].single[m] = -1;
				}
			}
		}
		n = 1;
		oldpop[j].parent[0] = 0;
		oldpop[j].parent[1] = 0;
		oldpop[j].xsite = 0;

		if (printForDebug)
		{
			for (k = 0; k < NumAllCities - 2; k++)
				printf("%d ", oldpop[j].single[k]);
			printf("\n");
		}

		objfunc(&(oldpop[j]));

		l = flip(0.05);

		if ((l == 1) || (count != 0))
		{
			if (oldpop[j].distance > maxdis)
			{
				j = j - 1;
				count++;
			} /*make sure at least 10% are feasible solutions*/
			else
				count = 0;
		}
	}
}

void initreport()
{

	// fprintf(outfp, "the parameters of GA:\n");
	// fprintf(outfp, "-------------------------------\n");
	// fprintf(outfp, "popsize=  %d\n", popsize);
	// fprintf(outfp, "maxgen=   %d\n", maxgen);
	// fprintf(outfp, "pcross=   %f\n", pcross);
	// fprintf(outfp, "pmutation=%f\n", pmutation);
	// fprintf(outfp, "-------------------------------\n");

	printf("the parameters of GA:\n");
	printf("-------------------------------\n");
	printf("popsize=  %d\n", popsize);
	printf("maxgen=   %d\n", maxgen);
	printf("pcross=   %f\n", pcross);
	printf("pmutation=%f\n", pmutation);
	printf("-------------------------------\n");
}

void generationNewCrossover()
{

	int i, mate1, mate2, curInd = 0, k; // curInd: current individual

	preselect();
	do
	{
		mate1 = select1();
		mate2 = select1();

		/*fprintf(outfp,"mate1=%d,mate2=%d\n",mate1,mate2);*/

		if (printForDebug)
		{
			printf("mate1=%d,mate2=%d\n", mate1, mate2);

			for (i = 0; i < NumAllCities - 2; i++)
			{
				printf("%d ", oldpop[mate1].single[i]);
			}
			printf("\n");
			for (i = 0; i < NumAllCities - 2; i++)
			{
				printf("%d ", oldpop[mate2].single[i]);
			}
			printf("\n");
		}

		// check if result if valid
		for (i = 0; i < NumAllCities - 2; i++)
		{
			if (oldpop[mate1].single[i] > NumAllCities - 2)
			{
				printf("ERROR in generation()!\n");
				printf("oldpop[%d].single[%d] = %d\n", mate1, i,
					   oldpop[mate1].single[i]);
				printf("it's not supposed to be greater than NumAllCities=%d\n",
					   NumAllCities);
				printf("terminatering...\n");
				exit(1);
			}
		}

		NewCrossover(oldpop[mate1].single, oldpop[mate2].single,
					 newpop[curInd].single);

		if (printForDebug)
		{
			printf("after crossover\n");
			for (i = 0; i < NumAllCities - 2; i++)
			{
				printf("%d ", newpop[curInd].single[i]);
			}
			printf("\n");
		}

		mutation(newpop[curInd].single);

		if (printForDebug)
		{
			printf("after mutaion\n");
			for (i = 0; i < NumAllCities - 2; i++)
			{
				printf("%d ", newpop[curInd].single[i]);
			}
			printf("\n");
		}

		moveMinusOneToEnd(&newpop[curInd]);

		string chrom1;

		for (i = 0; i < NumAllCities - 2; i++)
		{
			chrom1 += itos(newpop[curInd].single[i]) + " ";
		}
		// printf ("chrom1=%s; chrom2=%s\n", chrom1.c_str (),chrom2.c_str ());
		// cout<<chrom1<<"  |  "<<chrom2<<endl;

		// if chrom1 is in hashmap, copy the info from hashmap
		if (hashmap.find(chrom1) != hashmap.end())
			copyIndividual(&newpop[curInd], hashmap[chrom1]);
		else // if chrom1 is not in hashmap, call objfunc(), and add it to hashmap
		{
			objfunc(&(newpop[curInd]));
			newpop[curInd].parent[0] = mate1 + 1;
			newpop[curInd].parent[1] = mate2 + 1;

			individual ind;
			copyIndividual(&ind, newpop[curInd]);
			hashmap.insert(
				{chrom1, ind});
		}

		if (printForDebug)
		{
			cout << "================hashmap size: " << hashmap.size() << endl;
			for (const auto &entry : hashmap)
			{
				string key = entry.first;
				individual ind = entry.second;
				cout << "chro: " << key << endl;
				cout << "profit: " << ind.fitness << endl;
			}

			printf("after obj()\n");
			for (i = 0; i < NumAllCities - 2; i++)
			{
				printf("%d ", newpop[curInd].single[i]);
			}
			printf("\n");
		}

		curInd = curInd + 1;
	} while (curInd < popsize);
}

void generation()
{
	int i, mate1, mate2, curInd = 0, k; // curInd: current individual

	preselect();
	do
	{
		mate1 = select1();
		mate2 = select1();

		/*fprintf(outfp,"mate1=%d,mate2=%d\n",mate1,mate2);*/

		if (printForDebug)
		{
			printf("mate1=%d,mate2=%d\n", mate1, mate2);

			for (i = 0; i < NumAllCities - 2; i++)
			{
				printf("%d ", oldpop[mate1].single[i]);
			}
			printf("\n");
			for (i = 0; i < NumAllCities - 2; i++)
			{
				printf("%d ", oldpop[mate2].single[i]);
			}
			printf("\n");
		}

		for (i = 0; i < NumAllCities - 2; i++)
		{
			if (oldpop[mate1].single[i] > NumAllCities - 2)
			{
				printf("ERROR in generation()!\n");
				printf("oldpop[%d].single[%d] = %d\n", mate1, i,
					   oldpop[mate1].single[i]);
				printf("it's not supposed to be greater than NumAllCities=%d\n",
					   NumAllCities);
				printf("terminatering...\n");
				exit(1);
			}
		}

		jcross = crossover(oldpop[mate1].single, oldpop[mate2].single,
						   newpop[curInd].single, newpop[curInd + 1].single);

		if (printForDebug)
		{
			printf("after crossover\n");
			for (i = 0; i < NumAllCities - 2; i++)
			{
				printf("%d ", newpop[curInd].single[i]);
			}
			printf("\n");
			for (i = 0; i < NumAllCities - 2; i++)
			{
				printf("%d ", newpop[curInd + 1].single[i]);
			}
			printf("\n");
		}

		mutation(newpop[curInd].single);
		mutation(newpop[curInd + 1].single);

		if (printForDebug)
		{
			printf("after mutaion\n");
			for (i = 0; i < NumAllCities - 2; i++)
			{
				printf("%d ", newpop[curInd].single[i]);
			}
			printf("\n");
			for (i = 0; i < NumAllCities - 2; i++)
			{
				printf("%d ", newpop[curInd + 1].single[i]);
			}
			printf("\n");
		}

		moveMinusOneToEnd(&newpop[curInd]);
		moveMinusOneToEnd(&newpop[curInd + 1]);

		string chrom1;
		string chrom2;
		for (i = 0; i < NumAllCities - 2; i++)
		{
			chrom1 += itos(newpop[curInd].single[i]) + " ";
			chrom2 += itos(newpop[curInd + 1].single[i]) + " ";
		}
		// printf ("chrom1=%s; chrom2=%s\n", chrom1.c_str (),chrom2.c_str ());
		// cout<<chrom1<<"  |  "<<chrom2<<endl;

		// if chrom1 is in hashmap, copy the info from hashmap
		if (hashmap.find(chrom1) != hashmap.end())
			copyIndividual(&newpop[curInd], hashmap[chrom1]);
		else // if chrom1 is not in hashmap, call objfunc(), and add it to hashmap
		{
			objfunc(&(newpop[curInd]));
			newpop[curInd].parent[0] = mate1 + 1;
			newpop[curInd].parent[1] = mate2 + 1;
			newpop[curInd].xsite = jcross;

			individual ind;
			copyIndividual(&ind, newpop[curInd]);
			hashmap.insert(
				{chrom1, ind});
		}

		// check chrom2
		if (hashmap.find(chrom2) != hashmap.end())
			copyIndividual(&newpop[curInd + 1], hashmap[chrom2]);
		else // if chrom1 is not in hashmap, call objfunc(), and add it to hashmap
		{
			objfunc(&(newpop[curInd + 1]));
			newpop[curInd + 1].parent[0] = mate1 + 1;
			newpop[curInd + 1].parent[1] = mate2 + 1;
			newpop[curInd + 1].xsite = jcross;

			individual ind;
			copyIndividual(&ind, newpop[curInd + 1]);
			hashmap.insert(
				{chrom2, ind});
		}

		if (printForDebug)
		{
			cout << "================hashmap size: " << hashmap.size() << endl;
			for (const auto &entry : hashmap)
			{
				string key = entry.first;
				individual ind = entry.second;
				cout << "chro: " << key << endl;
				cout << "profit: " << ind.fitness << endl;
			}

			printf("after obj()\n");
			for (i = 0; i < NumAllCities - 2; i++)
			{
				printf("%d ", newpop[curInd].single[i]);
			}
			printf("\n");
			for (i = 0; i < NumAllCities - 2; i++)
			{
				printf("%d ", newpop[curInd + 1].single[i]);
			}
			printf("\n");
		}

		curInd = curInd + 2;
	} while (curInd < (popsize - 1));
}

void initmalloc()
{
	unsigned nbytes;
	int j;
	// struct individual.single= new double[NumAllCities-2];
	nbytes = popsize * sizeof(struct individual);

	// if ((oldpop = (struct individual *)malloc(nbytes)) == NULL)
	// 	// fprintf(outfp, "malloc:out of memory making oldpop!!\n");
	// 	printf("malloc:out of memory making oldpop!!\n");
	// if ((newpop = (struct individual *)malloc(nbytes)) == NULL)
	// 	// fprintf(outfp, "malloc:out of memory making newpop!!\n");
	// 	printf("malloc:out of memory making newpop!!\n");

	oldpop = (struct individual *)malloc(nbytes);
	if (oldpop == NULL)
		printf("malloc:out of memory making oldpop!!\n");

	newpop = (struct individual *)malloc(nbytes);
	if (newpop == NULL)
		printf("malloc:out of memory making newpop!!\n");
}

void freeall()
{
	int i;
	free(oldpop);
	free(newpop);
}

void report()
{
	int i, j, num = 0;
	/*fprintf(outfp,"\nreport:\n");
	 fprintf(outfp,"statics of %d generation:\n",gen);
	 fprintf(outfp,"total crossover=%ld,total mutation=%ld\n",ncross,nmutation);
	 fprintf(outfp,"min fitness:%f,max fitness:%f,avg fitness:%.2f\n",minFitness,maxFitness,avgFitnessGlobal);
	 if(minFitness==0)  fprintf(outfp,"there is no route fit for constraint in this generation.\n");*/
	if (bestfit.cusnum == 0)
	{
		// fprintf(outfp, "(%d)\n", run + 1);
		// fprintf(outfp, "No feasible solution.\n");

		printf("(%d)\n", run + 1);
		printf("No feasible solution.\n");
	}
	else
	{
		// fprintf(outfp, "(%d)\n", run + 1);
		// /*fprintf(outfp,"best individual ever found:\ngeneration:%d   ",bestfit.generation);*/
		// fprintf(outfp, "BEST INDIVIDUAL: generation: %d;", bestfit.generation);
		// fprintf(outfp, " fitness:%f;", bestfit.fitness);
		// fprintf(outfp, " distance:%f.\n", bestfit.distance);
		// fprintf(outfp, "best route :");

		printf("(%d)\n", run + 1);
		printf("BEST INDIVIDUAL: generation: %d;", bestfit.generation);
		printf(" fitness:%f;", bestfit.fitness);
		printf(" distance:%f.\n", bestfit.distance);
		printf("best route :");

		for (i = 0; i < NumAllCities - 2; i++)
			// fprintf(outfp, "%d,", bestfit.single[i]);
			printf("%d,", bestfit.single[i]);
		// fprintf(outfp, "\n");
		// fprintf(outfp, "ROUTE:");

		printf("\n");
		printf("ROUTE:");

		for (j = 0; j < bestfit.cusnum; j++)
		{
			// fprintf(outfp, " customer%d: %d-->%d ", j + 1, bestfit.cus[j][0],
			// 		bestfit.cus[j][1]);
			printf(" customer%d: %d-->%d ", j + 1, bestfit.cus[j][0],
				   bestfit.cus[j][1]);
		}
		// fprintf(outfp, "\n");
		printf("\n");
	}
}

void preselect()
{
	int j;
	sumfitness = 0;
	for (j = 0; j < popsize; j++)
	{
		sumfitness = sumfitness + oldpop[j].fitness;
	}
}

int select1()
{
	double k, pick;
	double sum;
	int i;

	k = rand();

	pick = k / RAND_MAX;

	sum = 0;
	if (sumfitness != 0)
	{
		for (i = 0; (sum < pick) && (i < popsize); i++)
		{
			sum += (double)oldpop[i].fitness / sumfitness;
		}
	}
	else
		i = 1 + (int)((double)popsize * rand() / (RAND_MAX + 1.0)); /*random(popsize)+1;*/
	return (i - 1);
}

void statistics(struct individual *pop)
{
	int i, j, h, m = 0, k = 0;
	sumfitness = 0;

	for (j = 0; j < popsize; j++)
	{
		if (pop[j].distance <= maxdis)
		{
			minFitness = pop[j].fitness;
			maxFitness = pop[j].fitness;
			m = 1;
			break;
		}
	}

	if (m == 0)
	{
		minFitness = 0;
		maxFitness = 0;
	}

	if (m > 0)
	{
		for (j = 0; j < popsize; j++)
		{
			if (pop[j].distance <= maxdis)
			{
				k++;
				sumfitness = sumfitness + pop[j].fitness;
				if (pop[j].fitness > maxFitness)
					maxFitness = pop[j].fitness;
				if (pop[j].fitness < minFitness)
					minFitness = pop[j].fitness;
				if (pop[j].fitness > bestfit.fitness)
				{
					for (i = 0; i < NumAllCities - 2; i++)
					{
						bestfit.single[i] = pop[j].single[i];
					}
					bestfit.fitness = pop[j].fitness;
					bestfit.generation = gen;
					bestfit.distance = pop[j].distance;
					bestfit.cusnum = pop[j].cusnum;

					for (h = 0; h < pop[j].cusnum; h++)
					{
						bestfit.cus[h][0] = pop[j].cus[h][0];
						bestfit.cus[h][1] = pop[j].cus[h][1];
					}
				}
			}
		}
	}
	if (k == 0)
		avgFitnessGlobal = 0;
	else
		avgFitnessGlobal = (double)sumfitness / k; /*this avg is the avg of pop that satisfy the distance constrain*/
}

void objfunc(struct individual *critter)
{
	/*int status = 0;*/
	double objval;
	double tp;
	// printf ("==> beginning of objfunc()\n");
	/*cus: the number of all potential customers on the route; such as 0->2->9,
	 * there are 3 custoemrs 0->2, 0->9, 2->9 in total*/
	int loc, cus;
	/*seg: number of segmetations, such as 0 -> 2 -> 9 (route[0]->route[1]->route[2]),
	 * there are 2 segmentations and 3 locations and 3 customers*/
	int seg;

	double pro = 0.0, dist = 0.0, dmn;
	int count1 = 0, count2 = 0, num = 0, i, j, k, m, n, h, s, l, pprocus = 0;
	/*pprocus: positive profit customer numbers*/
	int sta, end;
	int order[NumAllCities - 2], sequence[NumAllCities - 2],
		fatesingle[NumAllCities - 2], route[NumAllCities];
	// route[] is the real route
	int ncities = 0;
	// the number of cities between starting and ending points(not include starting and ending points)
	double objmax;
	int numwht, numwht1;
	// pprofit: profit of potential customer
	// sprofit: profit of selected customer
	// prorate: rate of profit and demand
	double pprofit[NumAllCities][NumAllCities],
		sprofit[NumAllCities][NumAllCities], prorate[NumAllCities][NumAllCities],
		weight[NumAllCities - 1];

	for (i = 0; i < NumAllCities - 1; i++) /*define weigh*/
	{
		weight[i] = 0;
	}

	for (i = 0; i < NumAllCities; i++) /*define profit and rate*/
		for (j = 0; j < NumAllCities; j++)
		{
			pprofit[i][j] = -1;
			sprofit[i][j] = -1;
			prorate[i][j] = -1;
		}

	for (i = 0; i < NumAllCities - 2; i++) /*change the code to the real route*/
	{
		order[i] = -1;
		sequence[i] = i + 1;
	}

	moveMinusOneToEnd(critter);
	/*
	 for (i = 0; i < NumAllCities - 2; i++)
	 {
	 if (critter->single[i] < 0)
	 count1++;
	 }

	 for (i = 0; i < NumAllCities - 2; i++)
	 {
	 if (critter->single[i] < 0)
	 {
	 for (j = i; j < NumAllCities - 3; j++)
	 critter->single[j] = critter->single[j + 1];
	 i = i - 1;
	 critter->single[NumAllCities - 3] = -1;
	 count2++;
	 if (count1 == count2)
	 break;
	 }
	 }
	 */
	if (printForDebug)
	{
		printf("===> inside objfunc() after reasign\n");
		for (i = 0; i < NumAllCities - 2; i++)
		{
			printf("%d ", critter->single[i]);
		}
		printf("\n");
	}

	for (i = 0; i < NumAllCities - 2; i++)
		fatesingle[i] = (*critter).single[i];

	for (i = 0; i < NumAllCities - 2; i++)
	{
		j = fatesingle[i];
		if (j < 0)
			break;
		else
		{
			order[i] = sequence[j - 1];
			m = j - 1;
			delet(sequence, &m);
		}
	}

	for (i = 0; i < NumAllCities - 2; i++) /*find the number of cities visited except starting point and destination*/
	{
		if (order[i] > 0)
			ncities = ncities + 1;
	}

	/*create route[10] including starting and ending locations*/
	route[0] = 0;
	route[ncities + 1] = NumAllCities - 1;

	for (i = 1; i < NumAllCities; i++)
	{
		if (i < ncities + 1)
			route[i] = order[i - 1];
		if (i > ncities + 1)
			route[i] = -1;
	}

	for (i = 0; i < ncities + 1; i++)
		dist = dist + dis[route[i]][route[i + 1]];
	critter->distance = dist;

	if (dist > maxdis)
	{
		critter->fitness = 0;
		// printf ("= middle of objfunc(). not feasible solution. return.\n");
		return;
	}

	/*
	 printf ("\n=========> route <=========\n");
	 for (i = 0; i < ncities + 2; i++)
	 {
	 printf ("%d  ", route[i] + 1);
	 }
	 printf ("\n");
	 printf ("dist=%lf\n", dist);
	 */

	if (UseGreedy)
	{
		for (i = 0; i < ncities + 1; i++)
			for (j = 0; j < ncities + 1 - i; j++)
			{
				// i is the number of cities between two points
				//(j<ncities+1-i) means that when there are i cities between w points, the last starting point is from route[ncities-i]
				pprofit[route[j]][route[j + i + 1]] = demand[route[j]][route[j + i + 1]] * dis[route[j]][route[j + i + 1]] * P;

				for (k = 0; k < i + 1; k++)
				{
					pprofit[route[j]][route[j + i + 1]] = pprofit[route[j]][route[j + i + 1]] - demand[route[j]][route[j + i + 1]] * dis[route[j + k]][route[j + k + 1]] * C;
				}

				pprofit[route[j]][route[j + i + 1]] = pprofit[route[j]][route[j + i + 1]] - wv * C * dist;

				if (pprofit[route[j]][route[j + i + 1]] >= 0)
				{
					pprocus = pprocus + 1;
					if (pprofit[route[j]][route[j + i + 1]] == 0)
						prorate[route[j]][route[j + i + 1]] = 0;
					else
						prorate[route[j]][route[j + i + 1]] =
							pprofit[route[j]][route[j + i + 1]] / demand[route[j]][route[j + i + 1]];
				}
			}

		for (s = 0; s < pprocus; s++)
		{
			numwht = 0;
			objmax = 0;
			// find the customer of smallest profit rate prorate[][]
			for (i = 0; i < ncities + 1; i++)
			{
				for (j = i + 1; j < ncities + 2; j++)
				{
					if (prorate[route[i]][route[j]] >= objmax)
					{
						objmax = prorate[route[i]][route[j]];
						sta = i;
						end = j; /*h and h are the i numbers of route[i]*/
						dmn = demand[route[i]][route[j]];
					}
				}
			}

			prorate[route[sta]][route[end]] = -1;

			for (i = sta; i < end; i++)
			{
				weight[i] = weight[i] + dmn;
				if (weight[i] > Q)
				{
					numwht = numwht + 1;
				}
			}

			if (numwht > 0)
			{
				for (i = sta; i < end; i++)
					weight[i] = weight[i] - dmn;
			}

			if (numwht == 0)
			{
				sprofit[route[sta]][route[end]] = pprofit[route[sta]][route[end]];
				pprofit[route[sta]][route[end]] = -1;
			}
		}

		for (i = 0; i < NumAllCities; i++)
			for (j = 0; j < NumAllCities; j++)
			{
				if (sprofit[i][j] >= 0)
				{
					pro = pro + sprofit[i][j];
					critter->cus[num][0] = i;
					critter->cus[num][1] = j;
					num = num + 1;
				}
			}
		critter->fitness = pro;
		critter->cusnum = num;
		// printf("cusnum=%d\n",num);
		// if (dist>maxdis) critter->fitness=A*(critter->fitness); */
		/*Here should be pro not dis, the previous files are wrong*/
	}

	if (UseGurobi)
	{
	}

	if (printForDebug)
	{
		printf("===> after gurobi and at the end of objfunc()\n");
		for (i = 0; i < NumAllCities - 2; i++)
		{
			printf("%d ", critter->single[i]);
		}
		printf("\n");
		printf("==> end of objfunc()\n");
		// exit (1);
	}
}

void delet(int *se, int *r)
{

	int i, j;
	i = *r;
	for (j = i; j < NumAllCities - 3; j++)
		se[j] = se[j + 1];
}

void mutation(int *child)
{
	int i, j, count;
	for (i = 0; i < NumAllCities - 2; i++)
	{
		count = 0;
		if (flip(pmutation))
		{
			if (child[i] == -1)
			{
				if (i == 0)
					child[0] = 1 + (int)((double)(NumAllCities - 2) * rand() / (RAND_MAX + 1.0)); /*random(8)+1;*/
				else
				{
					for (j = 0; j < i; j++)
					{
						if (child[j] < 0)
							count++;
					}
					child[i] = 1 + (int)((double)(NumAllCities - 2 - i + count) * rand() / (RAND_MAX + 1.0)); /*random(8-i+count)+1;count=0;*/
				}
			}
			else
			{
				child[i] = -1;
			}

			nmutation++;
		}
	}
}

void NewCrossover(int *parent1, int *parent2, int *child)
{
	int i;
	for (i = 0; i < NumAllCities - 2; i++)
	{
		if (flip(TAU))
			child[i] = parent1[i];
		else
			child[i] = parent2[i];
	}
}

int
// crossover (int parent1[NumAllCities - 2], int parent2[NumAllCities - 2], int child1[NumAllCities - 2],
//	   int child2[NumAllCities - 2])
crossover(int *parent1, int *parent2, int *child1, int *child2)
{

	int i, j, m = 0, n = 0, crossmax;
	if (flip(pcross))
	{
		for (i = 0; i < NumAllCities - 2; i++)
		{
			if (parent1[i] > 0)
				n = n + 1;
			if (parent2[i] > 0)
				m = m + 1;
		}

		if (n > m)
			crossmax = n;
		else
			crossmax = m;

		j = 1 + (int)((double)crossmax * rand() / (RAND_MAX + 1.0));	  /*random(crossmax);*/
		int j2 = 1 + (int)((double)crossmax * rand() / (RAND_MAX + 1.0)); /*random(crossmax);*/

		// swap j and j2 if j>j2
		if (j > j2)
		{
			int tmp = j;
			j = j2;
			j2 = tmp;
		}

		// cout<<"crossmax="<<crossmax<<endl;
		// cout<<"j="<<j<<" j2="<<j2<<endl;

		ncross++;
		for (i = 0; i < NumAllCities - 2; i++)
		{
			child1[i] = parent1[i];
			child2[i] = parent2[i];
		}
		/*
		 for (i = j; i < NumAllCities - 2; i++)
		 {
		 child1[i] = parent2[i];
		 child2[i] = parent1[i];
		 }
		 */
		//=== new feature: only swap genes between j and j2
		for (i = j; i <= j2; i++)
		{
			child1[i] = parent2[i];
			child2[i] = parent1[i];
		}
	}
	else
	{
		for (i = 0; i < NumAllCities - 2; i++)
		{
			child1[i] = parent1[i];
			child2[i] = parent2[i];
		}
		j = -1;
	}

	return (j);
}

int flip(double prob)
{
	double i;
	int j;
	j = 1 + (int)((double)100.0 * rand() / (RAND_MAX + 1.0)); /*random(100)+1;*/

	i = ((double)j) / 100;
	if (i < prob)
		return (1);
	else
		return (0);
}

//================> genetic algorithm main process <================

vector<int>
runga(int NumAllCities_val, double P_val, double C_val, double Q_val, double maxdis_val,
	  double wv_val, double **demand_val, double **dis_val)
{

	vector<int> bestChromosome;

	NumAllCities = NumAllCities_val;
	P = P_val;
	C = C_val;
	Q = Q_val;
	maxdis = maxdis_val;
	wv = wv_val;
	demand = demand_val;
	dis = dis_val;

	// cout << "#cities" << NumAllCities << endl;
	// cout << "P=" << P << endl;
	// cout << "C=" << C << endl;
	// cout << "Q=" << Q << endl;
	// cout << "maxdis=" << maxdis << endl;
	// for (int i = 0; i < NumAllCities; i++)
	// 	for (int j = 0; j < NumAllCities; j++)
	// 	{
	// 		printf("demand[%d][%d]=%lf ", i, j, demand[i][j]);
	// 		printf("dis[%d][%d]=%lf ", i, j, dis[i][j]);
	// 	}

	if ((UseGreedy == true && UseGurobi == true) || (UseGreedy == false && UseGurobi == false))
	{
		cout << "One of the methods must be true: greedy or gurobi." << endl;
		cout << "Please modify it in source code." << endl;
		return bestChromosome;
	}
	if (UseGreedy)
		printf("Obj calculation method: Greedy\n");
	else
		printf("Obj calculation method: Gurobi\n");

	int bsindex = (int)NumAllCities / 10 - 1;

	if (UseGreedy == false && UseGurobi == true)
	{
	}

	int i, j, z = 0;
	struct individual *temp;
	// double finalfitness[100];
	// double avgfitness, s1 = 0, s2 = 0, s3 = 0;
	// double avgs, jfc;
	// double sd;					 /*sd: standard deviation of every soltuion from best known solution*/
	double mainmax = 0, mainmin; /*the maximum and minimu value of one data set, and one data set could be run for a number of times "MAXRUN"*/
	time_t t;
	clock_t begin, end;
	double minute;
	double second;
	double tm[MAXRUN];								  /*record the run time of different runs*/
	double mintime, maxtime = 0, avgtime, totime = 0; /*totime=total time of MAXRUN times runnings*/
	// ds is now passed to runga() by argument
	// int ds;											  /*data set*/
	double dsfitness[10]; /*best fitness of every  data set. If running for more than 1 times then it is the maximum value of MAXRUN runnings*/
	double dstime[10];	  /*average computaion time of every data set. If running for more than 1 time then it is the average time of MAXRUN times*/
	double avgds;		  /*average value of 10 best fitnesses of data set*/

	// BESTSOL[i][] is the optimal obj value or best obj value found
	// for (i+1)*10 cities instances
	//[1][10] has to be the best obj value found for 20 ndoes
	//[2][10] has to be the best obj value found for 30 ndoes
	// double BESTSOL[7][10] =
	// 	{
	// 		{1.76214, 1.31778, 1.23058, 0.79174, 1.42076, 1.63532, 1.50894, 1.32678,
	// 		 1.37084, 0.91972},
	// 		{1.65874, 1.61320, 1.61530, 1.81528, 1.77702, 1.77340, 1.68558, 1.94288,
	// 		 1.68794, 1.66976},
	// 		{1.75576, 1.76158, 1.72974, 1.89554, 1.72964, 1.76224, 1.77678, 1.65434,
	// 		 1.64152, 1.67608},
	// 		{1.55192, 1.36283, 1.59519, 1.69155, 1.69395, 1.66352, 1.47887, 1.49914,
	// 		 1.45319, 1.63109},
	// 		{1.86315, 1.82456, 1.87558, 1.97310, 1.77176, 1.92162, 1.86842, 1.79938,
	// 		 1.83202, 1.83838},
	// 		{1.9411222, 1.89862, 1.8608388, 1.928323, 1.8470248, 1.9389642, 1.772803,
	// 		 1.8856282, 1.8543142, 1.9684246},
	// 		{1.882262, 1.944, 1.8882498, 1.90455, 1.8954476, 1.9201386, 1.89167,
	// 		 1.9656782, 1.9230848, 1.92505}};

	if (NumAllCities == 10 || NumAllCities == 20)
	{
		// popsize is a EVEN number!
		// If it's a odd number, like 25, then in generation(), it will refer to oldpop[i]
		// which doesn't exist!
		if (UseGurobi)
		{
			popsize = 26;
			maxgen = 50;
		}
		if (UseGreedy)
		{
			popsize = 50;
			maxgen = 100;
		}
		NUM_NO_IMPRO = 50;
	}
	else
	{
		// bsindex is the data set instance index
		// it's supposed to be dsindex, will fix it later
		// when NumAllCities=30, bsindex=2, so popsize=100, maxgen=400
		// when NumAllCities=40, bsindex=3, so popsize=200, maxgen=800......
		// when NumAllCities=50, bsindex=4, so popsize=400, maxgen=1600
		popsize = 25 * pow(2, bsindex);
		maxgen = 100 * pow(2, bsindex);
		// when = 30,150
		// when=40, 200
		// when=50, 250
		NUM_NO_IMPRO = 50 + (bsindex) * 50;

		//==============
		popsize = 500;
		maxgen = 1000;

		NUM_NO_IMPRO = 800;
	}

	if (popsize % 2 != 0)
	{
		printf("popsize has to be an even number!\n");
		printf("else, error will occur when running!\n");
		exit(1);
	}

	// string outputfile = "ga_" + to_string(NumAllCities) + "_" + to_string(ds) + "_output.txt";
	// outfp = fopen(outputfile.c_str(), "w");

	// if (outfp == NULL)
	// {
	// 	printf("File Open Error:\n");
	// 	exit(1);
	// }

	// fprintf(outfp, "pcross    = %lf\n", pcross);
	// fprintf(outfp, "pmutation = %lf\n", pmutation);
	// fprintf(outfp, "popsize   = %lf\n", popsize);
	// fprintf(outfp, "maxgen    = %lf\n", maxgen);
	printf("pcross    = %lf\n", pcross);
	printf("pmutation = %lf\n", pmutation);
	printf("popsize   = %lf\n", popsize);
	printf("maxgen    = %lf\n", maxgen);

	unsigned seed;
	// seed = (unsigned)time(&t);
	// srand (9);
	// seed=10011;
	// seed = 1704994534;
	seed = 300;
	srand(seed);
	printf("Random seed is %u\n", seed);
	// for (ds = 0; ds < MAXINSTANCE; ds++)
	//  for (ds = 5; ds < 6; ds++)
	{
		// fprintf(outfp, "========== instance %d ========\n", ds);
		// printf("========== instance %d ========\n", ds);
		// printf("========== i am here ========\n");
		// cout<<"I am here"<<endl;

		mainmax = 0;
		maxtime = 0;
		totime = 0;
		z = 0;

		hashmap.clear();
		printf("NumAllCities=%d, P=%lf, C=%lf, Q=%lf, maxdis=%lf, wv=%lf\n",
			   NumAllCities, P, C, Q, maxdis, wv);
		double oldProfit = 0;
		int numNoImpro = 0;

		for (run = 0; run < NUM_RESTART; run++)
		{
			// printf("The %dth run\n", run + 1);
			// cout<<"the run "<<run<<endl;
			begin = clock();
			initialize();

			oldProfit = 0;
			for (gen = 0; gen < maxgen; gen++)
			{

				if (UseNewCrossover)
					generationNewCrossover();
				else
					generation();

				statistics(newpop);
				if (gen == maxgen - 1)
					report();
				temp = newpop;
				newpop = oldpop;
				oldpop = temp;

				if (bestfit.fitness > oldProfit)
				{
					numNoImpro = 0;
					oldProfit = bestfit.fitness;
				}
				else
				{
					numNoImpro++;
					if (numNoImpro >= NUM_NO_IMPRO)
					{
						// cout << "Break before maxgen at gen " << gen
						// 	 << " since no improvement for long time" << endl;
						break;
					}
				}
			}
			// finalfitness[run] = bestfit.fitness;
			end = clock();
			tm[run] = (double)(end - begin);
			cout << "mainmax=" << mainmax << endl;
			if (bestfit.fitness > mainmax)
			{
				mainmax = bestfit.fitness;

				for (int i = 0; i < NumAllCities - 2; i++)
				{
					bestChromosome.push_back(bestfit.single[i]);
				}
			}
		}
		for (i = 0; i < NUM_RESTART; i++)
		{
			second = tm[i] / CLOCKS_PER_SEC;
			// fprintf(outfp, "Processing time: %lf computer time,  %lf seconds\n", tm[i], second);
			totime = totime + tm[i];
		}
		// fprintf(outfp, "Total Processing time: %lf seconds\n", totime / CLOCKS_PER_SEC);
		printf("Total Processing time: %lf seconds\n", totime / CLOCKS_PER_SEC);
	}

	// delete env;

	// fclose(outfp);
	freeall();

	//========> decode chromosome to real route <========
	vector<int> finalRoute;
	finalRoute.push_back(0);
	if (mainmax > 0)
	{
		/*create route[10] including starting and ending locations*/
		int count1 = 0, count2 = 0;

		int order[NumAllCities - 2], sequence[NumAllCities - 2],
			fatesingle[NumAllCities - 2];

		for (i = 0; i < NumAllCities - 2; i++) /*change the code to the real route*/
		{
			order[i] = -1;
			sequence[i] = i + 1;
		}
		// move all -1 to end
		for (i = 0; i < NumAllCities - 2; i++)
		{
			if (bestChromosome[i] < 0)
				count1++;
		}
		cout << "count1=" << count1 << end;
		for (i = 0; i < NumAllCities - 2; i++)
		{
			if (bestChromosome[i] < 0)
			{
				for (j = i; j < NumAllCities - 3; j++)
					bestChromosome[j] = bestChromosome[j + 1];
				i = i - 1;
				bestChromosome[NumAllCities - 3] = -1;
				count2++;
				if (count1 == count2)
					break;
			}
		}
		// end of moving -1 to end
		// if (printForDebug)
		if (true)
		{
			printf("===> inside objfunc() after reasign\n");
			for (i = 0; i < NumAllCities - 2; i++)
			{
				printf("%d ", bestChromosome[i]);
				cout << bestChromosome[i] << " " << endl;
			}
			printf("\n");
		}

		for (i = 0; i < NumAllCities - 2; i++)
			fatesingle[i] = bestChromosome[i];

		for (i = 0; i < NumAllCities - 2; i++)
		{
			j = fatesingle[i];
			if (j < 0)
				break;
			else
			{
				order[i] = sequence[j - 1];
				int m = j - 1;
				delet(sequence, &m);
			}
		}
		int ncities = 0;
		for (i = 0; i < NumAllCities - 2; i++) /*find the number of cities visited except starting point and destination*/
		{
			if (order[i] > 0)
				ncities = ncities + 1;
		}

		/*create route[10] including starting and ending locations*/

		for (i = 1; i < NumAllCities; i++)
		{
			if (i < ncities + 1)
				finalRoute.push_back(order[i - 1]);
			if (i == ncities + 1)
				finalRoute.push_back(NumAllCities - 1);
		}
	}
	else //if no route with positive profit found, then use 0 -> depot as the route
		finalRoute.push_back(NumAllCities - 1);

	return finalRoute;
}

#endif