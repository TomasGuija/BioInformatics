#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "Utils.h"
#include "DistAln.h"

#define SPECIAL -1

void printMatches(TypeAlignment aln) {
	TypePosition p = 0, p0=0, count = 0;
	while(p<aln.size) {
		while(p<aln.size && aln.sequence[0][p] == aln.empty)
			p++;
		if((aln.sequence[0][p]<aln.cardinal) && (aln.sequence[0][p] == aln.sequence[1][p])) {
			TypePosition start = p0;
			printf("%ld start %ld ", ++count, p0);
			while(p<aln.size && aln.sequence[0][p]<aln.cardinal && aln.sequence[0][p] == aln.sequence[1][p]) {
				p++;
				p0++;
				while(p<aln.size && aln.sequence[0][p] == aln.empty)
					p++;
			}
			printf("end %ld size %ld\n", p0-1, p0-start);
		} else {
			p++;
			p0++;
		}
	}
}

TypeFloat computeMatches(TypeAlignment aln) {
	TypePosition p;
	double dist = 0;
	
	for(p=0; p<aln.size; p++)
		if(aln.sequence[0][p] != aln.empty && aln.sequence[0][p]<aln.cardinal)
			if(aln.sequence[0][p] == aln.sequence[1][p])
				dist++;
	return dist;
}

TypeFloat computeMatchesBis(TypeAlignment aln) {
	TypePosition p, l =0, dist = 0;
	
	for(p=0; p<aln.size; p++)
		if(aln.sequence[0][p] != aln.empty && aln.sequence[0][p]<aln.cardinal && aln.sequence[1][p] != aln.empty && aln.sequence[1][p]<aln.cardinal) {
			l++;
			if(aln.sequence[0][p] == aln.sequence[1][p])
				dist++;
		}
	return ((double)dist)/((double)l);
}

TypeFloat computeNorm1(TypeAlignment aln) {
	TypePosition p, l =0, dist = 0;
	
	for(p=0; p<aln.size; p++)
		if(aln.sequence[0][p] != aln.empty && aln.sequence[0][p]<aln.cardinal && aln.sequence[1][p] != aln.empty && aln.sequence[1][p]<aln.cardinal) {
			l++;
			if(aln.sequence[0][p] == aln.sequence[1][p])
				dist++;
		}
	return 1.0-((double)dist)/((double)l);
}

TypeFloat computeNorm2(TypeAlignment aln) {
	TypePosition p, l0=0, l1=0, dist = 0;
	
	for(p=0; p<aln.size; p++) {
		if(aln.sequence[0][p] != aln.empty && aln.sequence[0][p]<aln.cardinal)
			l0++;
		if(aln.sequence[1][p] != aln.empty && aln.sequence[1][p]<aln.cardinal)
			l1++;
		if(aln.sequence[0][p] != aln.empty && aln.sequence[0][p]<aln.cardinal && aln.sequence[1][p] != aln.empty && aln.sequence[1][p]<aln.cardinal) {
			if(aln.sequence[0][p] == aln.sequence[1][p])
				dist++;
		}
	}
	return 1.0-((double)dist)/((double)MIN(l0, l1));
}


TypeDistance computeWholeDistancePair(TypeAlignment aln, TypeDistFunction *distfunc) {
	TypeDistance dist;
	TypeNumber i, j;
	TypeAlignment atmp;

	dist.number = aln.number;
	dist.name = aln.name;	
	dist.table = (TypeFloat*) monmalloc(((dist.number*(dist.number-1))/2)*sizeof(TypeFloat));

	atmp.number = 2;
	atmp.empty = aln.empty;
	atmp.size = aln.size;
	atmp.sequence = (TypeSymbol **) monmalloc(2*sizeof(TypeSymbol*));
	atmp.cardinal = aln.cardinal;
	
	for(i=1; i<dist.number; i++) {
		for(j=0; j<i; j++) {
#ifdef DEBUG
printf("\n\n%s (%ld) vs %s (%ld)\n", dist.name[i], aln.size, dist.name[j], aln.size);
#endif
			atmp.sequence[0] = aln.sequence[j];
			atmp.sequence[1] = aln.sequence[i];
#ifdef DEBUG
printMatches(atmp);
printf("\n");
#endif
			dist.table[(i*(i-1))/2+j] = distfunc(atmp);
#ifdef DEBUG
printf("\ndist %f\n", dist.table[(i*(i-1))/2+j]);
#endif
		}
	}
	monfree((void*)atmp.sequence);
	return dist;
}
