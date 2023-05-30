#ifndef DistSeqF
#define DistSeqF

#include <stdlib.h>
#include <stdio.h>
#include "Sequences.h"
#include "Distance.h"

typedef TypeFloat TypeDistFunction(TypeSetOfSequences);

TypeFloat computePham(TypeSetOfSequences set);
TypeFloat computeTrivial(TypeSetOfSequences set);

TypeDistance computeWholeDistancePair(TypeSetOfSequences set,  TypePosition order, TypeDistFunction *distfunc, int local);

#endif
