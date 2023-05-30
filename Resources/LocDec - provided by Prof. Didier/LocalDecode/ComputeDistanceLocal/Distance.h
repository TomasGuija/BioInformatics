
#ifndef DistancesF
#define DistancesF

#include <stdlib.h>
#include <stdio.h>
#include "Sequences.h"
typedef float TypeFloat;

typedef struct DISTANCES {
	char **name;
	TypeNumber number;
	TypeFloat *table;
} TypeDistance;



TypeDistance readDistances(FILE *f);
TypeDistance readDistancesNexus(FILE *f);
TypeFloat computeVarianceRow(TypeDistance dist);	
TypeFloat computeVariance(TypeDistance dist);
TypeFloat computeMean(TypeDistance dist);
TypeFloat computeCorrelation(TypeDistance dist0, TypeDistance dist1);
TypeFloat computeCorrelationLog(TypeDistance dist0, TypeDistance dist1);
TypeFloat distance(TypeDistance d, TypeNumber i, TypeNumber j);
TypeFloat distanceClassesMin(TypeDistance d, TypeNumber *tab1, TypeNumber n1, TypeNumber *tab2, TypeNumber n2);
TypeFloat distanceClassesMoy(TypeDistance d, TypeNumber *tab1, TypeNumber n1, TypeNumber *tab2, TypeNumber n2);
void printDistanceRaw(FILE *f, TypeDistance dist);
void printDistanceNexus(FILE *f, TypeDistance dist);
void printDistancePhylip(FILE *f, TypeDistance dist);
void printDistanceTable(FILE *f, TypeDistance dist);

#endif
