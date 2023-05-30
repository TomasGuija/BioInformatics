
#ifndef DistancesF
#define DistancesF

#include <stdlib.h>
#include <stdio.h>
#include "Utils.h"
typedef float TypeFloat;

typedef struct DISTANCES {
	char **name;
	TypeNumber number;
	TypeFloat *table;
} TypeDistance;



TypeDistance diffCorr(TypeDistance dist0, TypeDistance dist1);
TypeDistance diffRank(TypeDistance dist0, TypeDistance dist1);
void printWorst(FILE *fo, TypeDistance dist0, TypeDistance dist1);
TypeFloat computeKendall(TypeDistance dist0, TypeDistance dist1);
TypeFloat computeSpearman(TypeDistance dist0, TypeDistance dist1);
TypeDistance readDistances(FILE *f);
TypeDistance readDistancesNexus(FILE *f);
TypeDistance readDistancesNexusBis(FILE *f);
TypeFloat computeVarianceRow(TypeDistance dist);	
TypeFloat computeVariance(TypeDistance dist);
TypeFloat computeMean(TypeDistance dist);
TypeCoeffReg computeCorrelation(TypeDistance dist0, TypeDistance dist1);
TypeCoeffReg computeCorrelationLog(TypeDistance dist0, TypeDistance dist1);
TypeFloat distance(TypeDistance d, TypeNumber i, TypeNumber j);
TypeFloat distanceClassesMin(TypeDistance d, TypeNumber *tab1, TypeNumber n1, TypeNumber *tab2, TypeNumber n2);
TypeFloat distanceClassesMoy(TypeDistance d, TypeNumber *tab1, TypeNumber n1, TypeNumber *tab2, TypeNumber n2);
void printDistanceRaw(FILE *f, TypeDistance dist);
void printDistanceNexus(FILE *f, TypeDistance dist);
void printDistancePhylip(FILE *f, TypeDistance dist);
void printDistanceTable(FILE *f, TypeDistance dist);
void printPoint(FILE *f, TypeDistance dist0, TypeDistance dist1);

#endif
