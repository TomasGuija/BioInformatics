#ifndef DistAlnF
#define DistAlnF

#include <stdlib.h>
#include <stdio.h>
#include "Alignement.h"
#include "Distance.h"

typedef TypeFloat TypeDistFunction(TypeAlignment aln);

TypeFloat computeMatches(TypeAlignment aln);
TypeFloat computeNorm1(TypeAlignment aln);
TypeFloat computeNorm2(TypeAlignment aln);
void printMatches(TypeAlignment aln);

TypeDistance computeWholeDistancePair(TypeAlignment aln, TypeDistFunction *distfunc);

#endif
