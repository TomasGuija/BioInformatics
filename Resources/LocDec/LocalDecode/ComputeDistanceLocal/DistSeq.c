#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "Utils.h"
#include "Decode.h"
#include "DistSeq.h"


#define SPECIAL -1

typedef struct TRANS {
	struct TRANS *next;
	TypeSymbol to;
	TypePosition n;
} TypeTrans;

static TypeTrans *newTrans(TypeSymbol to, TypeTrans *next);
static void incTrans(TypeTrans **t, TypeSymbol c);
static int incTransBis(TypeTrans **t, TypeSymbol c);
static TypePosition getTrans(TypeTrans *t, TypeSymbol c);
static double getMarkovLogProb(TypeSymbol *s, TypeSymbol card, TypePosition size, TypePosition *occ, TypeTrans **trans, TypePosition sizer);

TypeTrans *newTrans(TypeSymbol to, TypeTrans *next) {
	TypeTrans *t;
	t = (TypeTrans*) monmalloc(sizeof(TypeTrans));
	t->to = to;
	t->next = next;
	t->n = 1;
}

void incTrans(TypeTrans **t, TypeSymbol c) {
	while(*t != NULL && (*t)->to < c)
		t = &((*t)->next);
	if(*t != NULL && (*t)->to == c)
		(*t)->n++;
	else
		*t = newTrans(c, *t);
}

int incTransBis(TypeTrans **t, TypeSymbol c) {
	while(*t != NULL && (*t)->to < c)
		t = &((*t)->next);
	if(*t != NULL && (*t)->to == c) {
		(*t)->n++;
		return 0;
	} else {
		*t = newTrans(c, *t);
		return 1;
	}
}

TypePosition getTrans(TypeTrans *t, TypeSymbol c) {
	while(t != NULL && t->to < c)
		t = t->next;
	if(t != NULL && t->to == c)
		return t->n;
	else
		return 1;
}

double getMarkovLogProb(TypeSymbol *s, TypeSymbol card, TypePosition size, TypePosition *occ, TypeTrans **trans, TypePosition sizer) {
	TypePosition p, b;
	double logprob;
	for(b=0; b<size && s[b]>=card; b++)
		;
	if(b<size)
		logprob = log(occ[s[b]])-log(sizer);
	else
		logprob = 0;
	for(p=b+1; p<size; p++) {
		if(s[p]>=card) {
			while(p<size && s[p]>=card)
				p++;
			if(p<size)
				logprob += log(occ[s[p]])-log(sizer);
		} else {
			logprob += log(getTrans(trans[s[p-1]], s[p]))-log(occ[s[p-1]]);
		}
	}
	return logprob;
}


TypeFloat computeTrivial(TypeSetOfSequences set) {
	TypePosition *occ0, *occ1, p, l0=0, l1=0, dist = 0, count = 0,
	*nocc, minocc, maxocc, tmp;
	TypeSymbol c;
	
	occ0 =(TypePosition*) monmalloc(set.cardinal*sizeof(TypePosition));
	occ1 =(TypePosition*) monmalloc(set.cardinal*sizeof(TypePosition));
	for(c=0; c<set.cardinal; c++) {
		occ0[c] = 0;
		occ1[c] = 0;
	}
	for(p=0; p<set.size[0]; p++)
		if(set.sequence[0][p]<set.cardinal) {
			l0++;
			occ0[set.sequence[0][p]]++;
		}
	for(p=0; p<set.size[1]; p++)
		if(set.sequence[1][p]<set.cardinal) {
			l1++;
			occ1[set.sequence[1][p]]++;
		}
	for(c=0; c<set.cardinal; c++) {
		dist += MIN(occ0[c], occ1[c]);
	}
	monfree((void*)occ0);
	monfree((void*)occ1);
	return 1-(((double) dist)/((double) MIN(l0, l1)));
}



TypeFloat computePham(TypeSetOfSequences set) {
	TypePosition *occ0, *occ1, l0=0, l1=0, n, p, norm = 0, max;
	TypeTrans **state0, **state1;
	TypeSymbol c;
	double LogProb00, LogProb01, LogProb10, LogProb11,dist = 0;
	
	state0 = (TypeTrans**) monmalloc(set.cardinal*sizeof(TypeTrans*));
	state1 = (TypeTrans**) monmalloc(set.cardinal*sizeof(TypeTrans*));
	occ0 =(TypePosition*) monmalloc(set.cardinal*sizeof(TypePosition));
	occ1 =(TypePosition*) monmalloc(set.cardinal*sizeof(TypePosition));
	for(c=0; c<set.cardinal; c++) {
		state0[c] = NULL;
		state1[c] = NULL;
		occ0[c] = 0;
		occ1[c] = 0;
	}
	for(p=0; p<set.size[0]; p++)
		if(set.sequence[0][p]<set.cardinal) {
			l0++;
			occ0[set.sequence[0][p]]++;
		}
	for(p=0; p<set.size[1]; p++)
		if(set.sequence[1][p]<set.cardinal) {
			l1++;
			occ1[set.sequence[1][p]]++;
		}
	for(p=1; p<set.size[0]; p++)
		if(set.sequence[0][p-1]<set.cardinal && set.sequence[0][p]<set.cardinal)
			incTrans(&(state0[set.sequence[0][p-1]]), set.sequence[0][p]);
	for(p=1; p<set.size[1]; p++)
		if(set.sequence[1][p-1]<set.cardinal && set.sequence[1][p]<set.cardinal)
			incTrans(&(state1[set.sequence[1][p-1]]), set.sequence[1][p]);
	LogProb00 = getMarkovLogProb(set.sequence[0], set.cardinal,  set.size[0], occ0, state0, l0);
	LogProb01 = getMarkovLogProb(set.sequence[0], set.cardinal, set.size[0], occ1, state1, l1);
	LogProb10 = getMarkovLogProb(set.sequence[1], set.cardinal, set.size[1], occ0, state0, l0);
	LogProb11 = getMarkovLogProb(set.sequence[1], set.cardinal, set.size[1], occ1, state1, l1);
	for(c=0; c<set.cardinal; c++) {
		while(state0[c] != NULL) {
			TypeTrans *tmp;
			tmp = state0[c]->next;
			monfree((void*) state0[c]);
			state0[c] = tmp;
		}
		while(state1[c] != NULL) {
			TypeTrans *tmp;
			tmp = state1[c]->next;
			monfree((void*) state1[c]);
			state1[c] = tmp;
		}
	}
	monfree((void*)state0);
	monfree((void*)state1);
	monfree((void*)occ0);
	monfree((void*)occ1);
	return 1-exp(((LogProb10-LogProb11)/((double)set.size[1])+(LogProb01-LogProb00)/((double)set.size[0]))/2);
}



TypeDistance computeWholeDistancePair(TypeSetOfSequences set,  TypePosition order, TypeDistFunction *distfunc, int local) {
	TypeDistance dist;
	TypeNumber i, j;
	TypeSetOfSequences stmp, dec;


	dist.number = set.number;
	dist.name = set.name;	
	dist.table = (TypeFloat*) monmalloc(((dist.number*(dist.number-1))/2)*sizeof(TypeFloat));

	stmp.number = 2;
	stmp.size = (TypePosition*) monmalloc(2*sizeof(TypePosition));;
	stmp.sequence = (TypeSymbol **) monmalloc(2*sizeof(TypeSymbol*));
	stmp.cardinal = set.cardinal;
	
	dec.number = 2;
	dec.size = (TypePosition*) monmalloc(2*sizeof(TypePosition));
	dec.sequence = (TypeSymbol**) monmalloc(2*sizeof(TypeSymbol*));

	if(order<=1) {
		for(i=1; i<dist.number; i++)
			for(j=0; j<i; j++) {
				dec.sequence[0] = set.sequence[j]; dec.size[0] = set.size[j];
				dec.sequence[1] = set.sequence[i]; dec.size[1] = set.size[i];
				dist.table[(i*(i-1))/2+j] = distfunc(dec);
			}
	} else {
		for(i=1; i<dist.number; i++) {
			for(j=0; j<i; j++) {
				TypeSuffixTree *suffixTree;
				TypeOneSequence one, resOne;
				TypeLocalTree *localTree;
				TypePosition posOne, posRes;
				TypeSymbol *table, symb;

				TypeNumber n, ind = (i*(i-1))/2+j;
				stmp.sequence[0] = set.sequence[j]; stmp.size[0] = set.size[j];
				stmp.sequence[1] = set.sequence[i]; stmp.size[1] = set.size[i];
				dec.size[0] = stmp.size[0];
				dec.size[1] = stmp.size[1];
				dec.sequence = (TypeSymbol**) monmalloc(dec.number*sizeof(TypeSymbol*));
				dec.sequence[0] = (TypeSymbol*) monmalloc(dec.size[0]*sizeof(TypeSymbol));
				dec.sequence[1] = (TypeSymbol*) monmalloc(dec.size[1]*sizeof(TypeSymbol));
				dist.table[ind] = POS_INFTY;
				one = toOne(stmp);
				suffixTree = computeSuffixTree(one);
				if(local)
					localTree = computeLocalTree(suffixTree, one, order);
				else
					localTree = initLocalTree(suffixTree, one, order);
				resOne = computeSequence(localTree);
				freeLocalTree(localTree);
				table = (TypeSymbol*) monmalloc(resOne.cardinal*sizeof(TypeSymbol));
				for(symb=0; symb<resOne.cardinal; symb++)
					table[symb] = SPECIAL;
				dec.cardinal = 0;
				posOne = 0;
				if(local) {
					for(n=0; n<2; n++) {
						for(posRes=0; posRes<dec.size[n]; posRes++) {
							if(stmp.sequence[n][posRes] < set.cardinal) {
								if(table[resOne.sequence[posOne]] == SPECIAL)
									table[resOne.sequence[posOne]] = dec.cardinal++;
								dec.sequence[n][posRes] = table[resOne.sequence[posOne]];
							}
							posOne++;
						}
						posOne++;
					}
					for(n=0; n<2; n++)
						for(posRes=0; posRes<dec.size[n]; posRes++)
							if(stmp.sequence[n][posRes] >= set.cardinal)
								dec.sequence[n][posRes] = stmp.sequence[n][posRes]-set.cardinal+dec.cardinal;
				} else {
					for(n=0; n<2; n++) {
						dec.size[n] = stmp.size[n]-order+1;
						for(posRes=0; posRes<dec.size[n]; posRes++) {
							if(stmp.sequence[n][posRes] < set.cardinal) {
								if(table[resOne.sequence[posOne]] == SPECIAL)
									table[resOne.sequence[posOne]] = dec.cardinal++;
								dec.sequence[n][posRes] = table[resOne.sequence[posOne]];
							}
							posOne++;
						}
						posOne+=order;
					}
					for(n=0; n<2; n++)
						for(posRes=0; posRes<stmp.size[n]; posRes++)
							if(stmp.sequence[n][posRes] >= set.cardinal){
								TypePosition p = posRes-order+1, fin = (posRes<dec.size[n])?posRes:dec.size[n]-1;
								for(p=(p>0)?p:0; p<=fin; p++)
									dec.sequence[n][p] = stmp.sequence[n][posRes]-set.cardinal+dec.cardinal;
							}
				}
				monfree((void*)resOne.sequence);
				monfree((void*)table);
				dist.table[ind] = distfunc(dec);
				monfree((void*)one.sequence);
				monfree((void*)dec.sequence[0]);
				monfree((void*)dec.sequence[1]);
				freeSuffixTree(suffixTree);
			}
		}
	}
	monfree((void*)dec.size);
	monfree((void*)dec.sequence);
	monfree((void*)stmp.size);
	monfree((void*)stmp.sequence);
	return dist;
}





















































TypeDistance computeWholeDistancePairBis(TypeSetOfSequences set,  TypePosition orderStart,  TypePosition orderEnd, TypeDistFunction *distfunc, int local) {
	TypeDistance dist;
	TypeNumber i, j;
	TypeSetOfSequences stmp, dec;


	dist.number = set.number;
	dist.name = set.name;	
	dist.table = (TypeFloat*) monmalloc(((dist.number*(dist.number-1))/2)*sizeof(TypeFloat));

	stmp.number = 2;
	stmp.size = (TypePosition*) monmalloc(2*sizeof(TypePosition));;
	stmp.sequence = (TypeSymbol **) monmalloc(2*sizeof(TypeSymbol*));
	stmp.cardinal = set.cardinal;
	
	dec.number = 2;
	dec.size = (TypePosition*) monmalloc(2*sizeof(TypePosition));
	dec.sequence = (TypeSymbol**) monmalloc(2*sizeof(TypeSymbol*));

	for(i=1; i<dist.number; i++) {
		for(j=0; j<i; j++) {
			TypeSuffixTree *suffixTree;
			TypeOneSequence one;
			TypePosition o;
			TypeNumber ind = (i*(i-1))/2+j;
			stmp.sequence[0] = set.sequence[j]; stmp.size[0] = set.size[j];
			stmp.sequence[1] = set.sequence[i]; stmp.size[1] = set.size[i];
			dec.size[0] = stmp.size[0];
			dec.size[1] = stmp.size[1];
			dec.sequence = (TypeSymbol**) monmalloc(dec.number*sizeof(TypeSymbol*));
			dec.sequence[0] = (TypeSymbol*) monmalloc(dec.size[0]*sizeof(TypeSymbol));
			dec.sequence[1] = (TypeSymbol*) monmalloc(dec.size[1]*sizeof(TypeSymbol));
			dist.table[ind] = POS_INFTY;
			one = toOne(stmp);
			suffixTree = computeSuffixTree(one);
			for(o=orderStart; o<=orderEnd; o++) {
				TypeLocalTree *localTree;
				TypeNumber n;
				TypePosition posOne, posRes;
				TypeSymbol *table, symb;
				TypeOneSequence resOne;
				TypeFloat disto;

				if(local)
					localTree = computeLocalTree(suffixTree, one, o);
				else
					localTree = initLocalTree(suffixTree, one, o);
				resOne = computeSequence(localTree);
				freeLocalTree(localTree);
				table = (TypeSymbol*) monmalloc(resOne.cardinal*sizeof(TypeSymbol));
				for(symb=0; symb<resOne.cardinal; symb++)
					table[symb] = SPECIAL;
				dec.cardinal = 0;
				posOne = 0;
				if(local) {
					for(n=0; n<2; n++) {
						for(posRes=0; posRes<dec.size[n]; posRes++) {
							if(stmp.sequence[n][posRes] < set.cardinal) {
								if(table[resOne.sequence[posOne]] == SPECIAL)
									table[resOne.sequence[posOne]] = dec.cardinal++;
								dec.sequence[n][posRes] = table[resOne.sequence[posOne]];
							}
							posOne++;
						}
						posOne++;
					}
					for(n=0; n<2; n++)
						for(posRes=0; posRes<dec.size[n]; posRes++)
							if(stmp.sequence[n][posRes] >= set.cardinal)
								dec.sequence[n][posRes] = stmp.sequence[n][posRes]-set.cardinal+dec.cardinal;
				} else {
					for(n=0; n<2; n++) {
						dec.size[n] = stmp.size[n]-o+1;
						for(posRes=0; posRes<dec.size[n]; posRes++) {
							if(stmp.sequence[n][posRes] < set.cardinal) {
								if(table[resOne.sequence[posOne]] == SPECIAL)
									table[resOne.sequence[posOne]] = dec.cardinal++;
								dec.sequence[n][posRes] = table[resOne.sequence[posOne]];
							}
							posOne++;
						}
						posOne+=o;
					}
					for(n=0; n<2; n++)
						for(posRes=0; posRes<stmp.size[n]; posRes++)
							if(stmp.sequence[n][posRes] >= set.cardinal){
								TypePosition p = posRes-o+1, fin = (posRes<dec.size[n])?posRes:dec.size[n]-1;
								for(p=(p>0)?p:0; p<=fin; p++)
									dec.sequence[n][p] = stmp.sequence[n][posRes]-set.cardinal+dec.cardinal;
							}
				}
				monfree((void*)resOne.sequence);
				monfree((void*)table);
				disto = distfunc(dec);
				if(disto<dist.table[ind])
					dist.table[ind] = disto;
			}
			monfree((void*)one.sequence);
			monfree((void*)dec.sequence[0]);
			monfree((void*)dec.sequence[1]);
			freeSuffixTree(suffixTree);
		}

	}
	monfree((void*)dec.size);
	monfree((void*)dec.sequence);
	monfree((void*)stmp.size);
	monfree((void*)stmp.sequence);
	return dist;
}
