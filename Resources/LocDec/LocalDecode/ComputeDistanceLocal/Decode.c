#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "Sequences.h"
#include "Utils.h"
#include "Decode.h"

#define INFTY LONG_MAX
#define SPECIAL -1

#define LEFT 1
#define RIGHT 1 << 1
#define DONE 1 << 2
#define NEW 1 << 3

static TypePosition ident;

typedef struct TMPLOCALTREE {
	TypeLocalNode **buffer;
	TypePosition size, bufferSize, incSize;
} TypeTmpLocalTree;

static int isLeft(TypeLocalNode *node);
static void leftOn(TypeLocalNode *node);
static void leftOff(TypeLocalNode *node);
static int isRight(TypeLocalNode *node);
static void rightOn(TypeLocalNode *node);
static void rightOff(TypeLocalNode *node);
static int isDone(TypeLocalNode *node);
static void doneOn(TypeLocalNode *node);
static void doneOff(TypeLocalNode *node);
static int isNew(TypeLocalNode *node);
static void newOn(TypeLocalNode *node);
static void newOff(TypeLocalNode *node);
static TypeLink *newLink(TypePosition position, TypeLink *next);
static void freeLink(TypeLink *link);
static TypeSuffixNode *newSuffixNode(TypePosition start, TypePosition end);
static TypeLocalNode *newLocalNode(TypePosition pos);
static void freeNode(TypeSuffixNode *node);
static void printSuffixNode(FILE *f, TypeSuffixNode* node, TypePosition depth);
static TypeLocalNode* getAncestor(TypeLocalNode* node);
static void printLocalNode(FILE *f, TypeLocalNode* node, int *done, TypePosition depth);
static TypeSuffixNode *findTransition(TypeSuffixNode *node, TypeSymbol symbol, TypeSuffixTree *suffixTree);
static void addTransition(TypeSuffixNode *node, TypeSuffixNode *toadd, TypeSuffixTree *suffixTree);
static void freeLocalNode(TypeLocalNode *node, TypeLocalTree *localTree);

static void canonize(TypeSuffixNode *node, TypePosition k, TypePosition p, TypeSuffixNode **sres, TypePosition *kres, TypeSuffixTree *suffixTree);
static void update(TypeSuffixNode *s, TypePosition k, TypePosition i, TypeSuffixNode **sres, TypePosition *kres, TypeSuffixTree *suffixTree);
static int test_and_split(TypeSuffixNode *s, TypePosition k, TypePosition p, TypeSymbol t, TypeSuffixTree *suffixTree, TypeSuffixNode **sres);

static void writeSymbol(TypeLocalNode *node, TypeSymbol *sequence, TypeSymbol symb);
static void parse_and_check(TypeSuffixNode *suffixNode, TypeLocalTree *local, TypeTmpLocalTree *tmpLocal, TypePosition depth, TypePosition order);
static void parse_and_set(TypeSuffixNode *suffixNode, TypeLocalTree *local, TypeLocalNode *node, TypePosition depth);
static TypeLocalTree *iterateLocalTree(TypeTmpLocalTree *cur, TypeTmpLocalTree *new);

TypeDistance computeDistance(TypeSetOfSequences set) {
	TypeDistance dist;
	TypeNumber i, j;
	TypePosition **occ, n, p;
	
	occ = (TypePosition**) monmalloc(set.number*sizeof(TypePosition*));
	for(n=0; n<set.number; n++) {
		occ[n] = (TypePosition*) monmalloc(set.cardinal*sizeof(TypePosition));
		for(p=0; p<set.cardinal; p++)
			occ[n][p] = 0;
		for(p=0; p<set.size[n]; p++)
			occ[n][set.sequence[n][p]]++;
	}
	dist.number = set.number;
	dist.name = set.name;
	dist.table = (TypeFloat*) monmalloc(((dist.number*(dist.number-1))/2)*sizeof(TypeFloat));
	for(i=1; i<dist.number; i++)
		for(j=0; j<i; j++) {
			TypeNumber ind = (i*(i-1))/2+j;
			TypeFloat sizeMin = (TypeFloat) MIN(set.size[i], set.size[j]);
			dist.table[ind] = 0.0;
			for(p=0; p<set.cardinal; p++)
				dist.table[ind] += ((TypeFloat) MIN(occ[i][p], occ[j][p]));
			dist.table[ind] = 1-dist.table[ind]/sizeMin;
		}
	for(n=0; n<set.number; n++)
		monfree((void*) occ[n]);
	monfree((void*) occ);
	return dist;
}



/*print a Fasta format file including all the sequences of the set - each decoded symbol is designed by 'X<number>' 
where X is the symbol in the same position in the non-decoded sequence and <number> an identifiant number*/
void printDecodedSequencesFasta(FILE *f, TypeSetOfSequences set, TypeSetOfSequences dec, int sizeLine) {
	TypeNumber n;
	TypeSymbol *nIdent, *sIdent, s;
	
	nIdent = (TypeSymbol*) monmalloc(dec.cardinal*sizeof(TypeSymbol));
	for(s=0; s<dec.cardinal; s++)
		nIdent[s] = SPECIAL;
	sIdent = (TypeSymbol*) monmalloc(set.cardinal*sizeof(TypeSymbol));
	for(s=0; s<set.cardinal; s++)
		sIdent[s] = 0;
	for(n=0; n<set.number; n++) {
		TypePosition position;
		if(set.name != NULL)
			fprintf(f, ">%s\n", set.name[n]);
		else
			fprintf(f, ">Sequence %d\n", n+1);
		for(position=0; position<set.size[n]; position += sizeLine) {
			TypePosition j, maxLine;
			maxLine = position+sizeLine;
			maxLine = maxLine>set.size[n]?set.size[n]:maxLine;
			if(set.table != NULL) {
				for(j=position; j<maxLine; j++) {
					if(nIdent[dec.sequence[n][j]] == SPECIAL)
						nIdent[dec.sequence[n][j]] = sIdent[set.sequence[n][j]]++;
					fprintf(f, "%c%ld ", set.table[set.sequence[n][j]], nIdent[dec.sequence[n][j]]);
				}
			} else {
				for(j=position; j<maxLine; j++) {
					if(nIdent[dec.sequence[n][j]] == SPECIAL)
						nIdent[dec.sequence[n][j]] = sIdent[set.sequence[n][j]]++;
					fprintf(f, "%ld-%ld ", set.sequence[n][j], nIdent[set.sequence[n][j]]);
				}
			}
			fprintf(f, "\n");
		}
		fprintf(f, "\n\n");
	}
	monfree((void*) nIdent);
	monfree((void*) sIdent);
}

/*print a Fasta format file including all the sequences of the set - each decoded symbol is designed by 'X<number>' 
where X is the symbol in the same position in the non-decoded sequence and <number> an identifiant number*/
void printDecodedSequencesFastaNumber(FILE *f, TypeSetOfSequences dec, int sizeLine) {
	TypeNumber n;
	for(n=0; n<dec.number; n++) {
		TypePosition position;
		if(dec.name != NULL)
			fprintf(f, ">%s\n", dec.name[n]);
		else
			fprintf(f, ">Sequence %d\n", n+1);
		for(position=0; position<dec.size[n]; position += sizeLine) {
			TypePosition j, maxLine;
			maxLine = position+sizeLine;
			maxLine = maxLine>dec.size[n]?dec.size[n]:maxLine;
			for(j=position; j<maxLine; j++) {
					fprintf(f, "%ld ", dec.sequence[n][j]);
			}
			fprintf(f, "\n");
		}
		fprintf(f, "\n\n");
	}
}

/*order-Local-decode the set of sequences set*/
TypeSetOfSequences slideDecode(TypeSetOfSequences set,  TypePosition order) {
	TypeOneSequence resOne, one;
	TypeSetOfSequences res;
	TypeNumber n;
	TypePosition posOne, posRes;
	TypeSymbol *table, symb;
	
	one = toOne(set);
	resOne = slideDecodeOneSequence(one, order);
	monfree((void*)one.sequence);
	
	table = (TypeSymbol*) monmalloc(resOne.cardinal*sizeof(TypeSymbol));
	for(symb=0; symb<resOne.cardinal; symb++)
		table[symb] = SPECIAL;
	res.name = set.name;
	res.table = set.table;
	res.number = set.number;
	res.size = (TypePosition*) monmalloc(res.number*sizeof(TypePosition));
	res.sequence = (TypeSymbol**) monmalloc(res.number*sizeof(TypeSymbol*));
	res.cardinal = 0;
	posOne = 0;
	for(n=0; n<res.number; n++) {
		res.size[n] = set.size[n]-order+1;
		res.sequence[n] = (TypeSymbol*) monmalloc(res.size[n]*sizeof(TypeSymbol));
		for(posRes=0; posRes<res.size[n]; posRes++) {
			if(set.sequence[n][posRes] < set.cardinal) {
				if(table[resOne.sequence[posOne]] == SPECIAL)
					table[resOne.sequence[posOne]] = res.cardinal++;
				res.sequence[n][posRes] = table[resOne.sequence[posOne]];
			}
			posOne++;
		}
		posOne+=order;
	}
	for(n=0; n<res.number; n++)
		for(posRes=0; posRes<set.size[n]; posRes++)
			if(set.sequence[n][posRes] >= set.cardinal){
				TypePosition p = posRes-order+1, fin = (posRes<res.size[n])?posRes:res.size[n]-1;
				for(p=(p>0)?p:0; p<=fin; p++)
					res.sequence[n][p] = set.sequence[n][posRes]-set.cardinal+res.cardinal;
			}
	monfree((void*)resOne.sequence);
	monfree((void*)table);
	return res;
}

/*order-Local-decode the set of sequences set*/
TypeSetOfSequences localDecode(TypeSetOfSequences set,  TypePosition order) {
	TypeOneSequence resOne, one;
	TypeSetOfSequences res;
	TypeNumber n;
	TypePosition posOne, posRes;
	TypeSymbol *table, symb;
	
	one = toOne(set);
	resOne = localDecodeOneSequence(one, order);
	monfree((void*)one.sequence);
	
	table = (TypeSymbol*) monmalloc(resOne.cardinal*sizeof(TypeSymbol));
	for(symb=0; symb<resOne.cardinal; symb++)
		table[symb] = SPECIAL;
	res.name = set.name;
	res.table = set.table;
	res.number = set.number;
	res.size = (TypePosition*) monmalloc(res.number*sizeof(TypePosition));
	res.sequence = (TypeSymbol**) monmalloc(res.number*sizeof(TypeSymbol*));
	res.cardinal = 0;
	posOne = 0;
	for(n=0; n<res.number; n++) {
		res.size[n] = set.size[n];
		res.sequence[n] = (TypeSymbol*) monmalloc(res.size[n]*sizeof(TypeSymbol));
		for(posRes=0; posRes<res.size[n]; posRes++) {
			if(set.sequence[n][posRes] < set.cardinal) {
				if(table[resOne.sequence[posOne]] == SPECIAL)
					table[resOne.sequence[posOne]] = res.cardinal++;
				res.sequence[n][posRes] = table[resOne.sequence[posOne]];
			}
			posOne++;
		}
		posOne++;
	}
	for(n=0; n<res.number; n++)
		for(posRes=0; posRes<res.size[n]; posRes++)
			if(set.sequence[n][posRes] >= set.cardinal)
				res.sequence[n][posRes] = set.sequence[n][posRes]-set.cardinal+res.cardinal;	
	monfree((void*)resOne.sequence);
	monfree((void*)table);
	return res;
}

/*order-Local-decode the sequence seq*/
TypeOneSequence slideDecodeOneSequence(TypeOneSequence seq, TypePosition order) {
	TypeLink **link;
	TypeSuffixTree *suffixTree;
	TypeLocalTree *localTree;
	TypeOneSequence res;
	TypePosition pos;
	suffixTree = computeSuffixTree(seq);
	localTree = initLocalTree(suffixTree, seq, order);
	freeSuffixTree(suffixTree);
	res = computeSequence(localTree);
	freeLocalTree(localTree);
	return res;
}

/*order-Local-decode the sequence seq*/
TypeOneSequence localDecodeOneSequence(TypeOneSequence seq, TypePosition order) {
	TypeSuffixTree *suffixTree;
	TypeLocalTree *localTree;
	TypeOneSequence res;
	TypePosition pos;
	suffixTree = computeSuffixTree(seq);
	localTree = computeLocalTree(suffixTree, seq, order);
	freeSuffixTree(suffixTree);
	res = computeSequence(localTree);
	freeLocalTree(localTree);
	return res;
}

/*compute the sequence from the position graph link in the following way:
each position in a conex component has the same symbol*/
TypeOneSequence computeSequence(TypeLocalTree *local) {
	TypeOneSequence seq;
	TypePosition pos;
	seq.size = local->size;
	seq.sequence = (TypeSymbol*) monmalloc(seq.size*sizeof(TypeSymbol));
	seq.cardinal = 0;
	for(pos = 0; pos<seq.size; pos++)
		seq.sequence[pos] = SPECIAL;
	for(pos = 0; pos<seq.size; pos++) {
		if(seq.sequence[pos] == SPECIAL)
			writeSymbol(getAncestor(local->leaf[pos]), seq.sequence, seq.cardinal++);
	}	
	return seq;
}

void writeSymbol(TypeLocalNode *node, TypeSymbol *sequence, TypeSymbol symb) {
	if(node->trans == NULL) {
		sequence[node->position] = symb;
	} else {
		for(node=node->trans; node!=NULL; node = node->next)
			writeSymbol(node, sequence, symb);
	}
}	

/*recursively parse all the position connected to pos in the graph link
used by computeSequence*/
void parse_and_check(TypeSuffixNode *suffixNode, TypeLocalTree *local, TypeTmpLocalTree *tmpLocal, TypePosition depth, TypePosition order) {
	if(suffixNode->end != INFTY) {
		TypeSuffixNode *tmp;
		depth += suffixNode->end-suffixNode->start+1;
		if(depth>=order) {
			TypeLocalNode *node;
			if(tmpLocal->size >= tmpLocal->bufferSize) {
				tmpLocal->bufferSize += tmpLocal->incSize;
				tmpLocal->buffer = (TypeLocalNode**) monrealloc(tmpLocal->buffer, tmpLocal->bufferSize*sizeof(TypeLocalNode*));
			}
			node = newLocalNode(ident++);
			tmpLocal->buffer[tmpLocal->size++] = node;
			if(depth == order) {
				newOn(node);
			}
			for(tmp=suffixNode->trans; tmp!=NULL; tmp=tmp->next)
				parse_and_set(tmp, local, node, depth);
		} else {
			for(tmp=suffixNode->trans; tmp!=NULL; tmp=tmp->next)
				parse_and_check(tmp, local, tmpLocal, depth, order);
		}
	}
}

/*recursively parse all the position connected to pos in the graph link
used by computeSequence*/
void parse_and_set(TypeSuffixNode *suffixNode, TypeLocalTree *local, TypeLocalNode *node, TypePosition depth) {
	if(suffixNode->end == INFTY) {
		TypeLocalNode *tmp;
		local->leaf[suffixNode->start-depth]->ancestor = node;
		local->leaf[suffixNode->start-depth]->next = node->trans;
		node->trans = local->leaf[suffixNode->start-depth];
		if(suffixNode->start > depth) {
			local->leaf[suffixNode->start-depth-1]->nlink = node;
			local->leaf[suffixNode->start-depth-1]->pnext = node->ptrans;
			node->ptrans = local->leaf[suffixNode->start-depth-1];
		}
	} else {
		TypeSuffixNode *tmp;
		depth += suffixNode->end-suffixNode->start+1;
		for(tmp=suffixNode->trans; tmp!=NULL; tmp=tmp->next)
			parse_and_set(tmp, local, node, depth);
	}
}

TypeLocalTree *initLocalTree(TypeSuffixTree *suffixTree, TypeOneSequence seq, TypePosition order) {
	TypePosition pos, size, k, maxSize = 0;
	TypeLocalTree *local;
	TypeTmpLocalTree *cur, *next;
	
	local = (TypeLocalTree*) monmalloc(sizeof(TypeLocalTree));
	local->size = seq.size;
	local->leaf = (TypeLocalNode**) monmalloc(local->size*sizeof(TypeLocalNode*));
	local->leaf[0] = newLocalNode(0);
	for(pos=1; pos<seq.size; pos++) {
		local->leaf[pos] = newLocalNode(pos);
		local->leaf[pos]->ptrans = local->leaf[pos-1];
		local->leaf[pos-1]->nlink = local->leaf[pos];
	}
	ident = seq.size;
	cur = (TypeTmpLocalTree *) monmalloc(sizeof(TypeTmpLocalTree));
	cur->size = 0;
	cur->incSize = (seq.size/10>10)?seq.size/10:10;
	cur->bufferSize = cur->incSize;
	cur->buffer = (TypeLocalNode**) monmalloc(cur->bufferSize*sizeof(TypeLocalNode*));
	parse_and_check(suffixTree->root, local, cur, -1, order);
	monfree((void*)cur->buffer);
	monfree((void*)cur);
	return local;
}

TypeLocalTree *computeLocalTree(TypeSuffixTree *suffixTree, TypeOneSequence seq, TypePosition order) {
	TypePosition pos, size, k, maxSize = 0;
	TypeLocalTree *local;
	TypeTmpLocalTree *cur, *next;
	
	local = (TypeLocalTree*) monmalloc(sizeof(TypeLocalTree));
	local->size = seq.size;
	local->leaf = (TypeLocalNode**) monmalloc(local->size*sizeof(TypeLocalNode*));
	local->leaf[0] = newLocalNode(0);
	for(pos=1; pos<seq.size; pos++) {
		local->leaf[pos] = newLocalNode(pos);
		local->leaf[pos]->ptrans = local->leaf[pos-1];
		local->leaf[pos-1]->nlink = local->leaf[pos];
	}
	ident = seq.size;
	cur = (TypeTmpLocalTree *) monmalloc(sizeof(TypeTmpLocalTree));
	cur->size = 0;
	cur->incSize = (seq.size/10>10)?seq.size/10:10;
	cur->bufferSize = cur->incSize;
	cur->buffer = (TypeLocalNode**) monmalloc(cur->bufferSize*sizeof(TypeLocalNode*));
	parse_and_check(suffixTree->root, local, cur, -1, order);
	for(pos=0; pos<cur->size; pos++) {
		if(!isNew(cur->buffer[pos])) {
			TypeLocalNode *tmp, **prev;
			cur->buffer[pos]->nlink = cur->buffer[pos]->trans->nlink;
			for(tmp=cur->buffer[pos]->trans; tmp != NULL; tmp = tmp->next)
				doneOn(tmp);
			prev = &cur->buffer[pos]->nlink->ptrans;
			for(tmp = cur->buffer[pos]->nlink->ptrans; tmp != NULL; tmp = tmp->pnext) {
				if(!isDone(tmp)) {
					*prev = tmp;
					prev = &tmp->pnext;
				}
			}
			*prev = tmp;
			cur->buffer[pos]->pnext = cur->buffer[pos]->nlink->ptrans;
			cur->buffer[pos]->nlink->ptrans = cur->buffer[pos];
			for(tmp=cur->buffer[pos]->trans; tmp!=NULL; tmp = tmp->next)
				doneOff(tmp);
		}
	}
	size = 0;
	for(pos=0; pos<cur->size; pos++) {
		if(isNew(cur->buffer[pos])) {
			cur->buffer[size++] = cur->buffer[pos];
			newOff(cur->buffer[pos]);
		}
	}
	cur->size = size;
	next = (TypeTmpLocalTree *) monmalloc(sizeof(TypeTmpLocalTree));
	next->bufferSize = cur->size;
	next->buffer = (TypeLocalNode**) monmalloc(next->bufferSize*sizeof(TypeLocalNode*));
	cur->incSize = seq.size;
	next->incSize = seq.size;
	for(k=1; k<order; k++) {
		TypeTmpLocalTree *tmp;
		TypePosition i;	
		next->size = 0;
		iterateLocalTree(cur, next);
		tmp = cur; cur = next; next = tmp;
/*printf("order %ld end %ld nodes\n", k, ident);*/
	}
	monfree((void*)cur->buffer);
	monfree((void*)next->buffer);
	monfree((void*)cur);
	monfree((void*)next);
	return local;
}

TypeLocalTree *iterateLocalTree(TypeTmpLocalTree *cur, TypeTmpLocalTree *new) {
	TypePosition ind, sizer, sizel;
	TypeLocalNode **savel, **saver;
	savel = (TypeLocalNode **) monmalloc(cur->incSize*sizeof(TypeLocalNode *));
	saver = (TypeLocalNode **) monmalloc(cur->incSize*sizeof(TypeLocalNode *));
	for(ind=0; ind<cur->size; ind++) {
		if(!isDone(cur->buffer[ind])) {
			TypePosition i, indl;
			sizer = 0; sizel = 1;
			savel[0] = cur->buffer[ind]; leftOn(cur->buffer[ind]);
			for(indl=0; indl<sizel; indl++) {
				if(savel[indl]->nlink == NULL) {
					TypeLocalNode *tmp;
					doneOn(savel[indl]);
					for(tmp=savel[indl]->trans; tmp!=NULL; tmp = tmp->next) {		
						if(!isRight(tmp->nlink)/* && !isNew(tmp->nlink)*/) {
							TypeLocalNode *bis;
							saver[sizer++] = tmp->nlink; rightOn(tmp->nlink);
							for(bis=tmp->nlink->ptrans; bis!=NULL; bis=bis->pnext) {
								TypeLocalNode *tostack;
								tostack = (bis->ancestor != NULL && !isNew(bis->ancestor))?bis->ancestor:bis;
								if(!isLeft(tostack)/* && !isNew(tostack)*/) {
									savel[sizel++] = tostack; leftOn(tostack);
								}
							}
						}
					}
				}
			}
if(sizer<2)
printf("Number %ld nodes\n", sizer);
			new->buffer[new->size] = newLocalNode(ident++); newOn(new->buffer[new->size]);
			new->buffer[new->size]->trans = saver[0];
			saver[0]->ancestor = new->buffer[new->size];
			for(i=1; i<sizer; i++) {
				saver[i-1]->next = saver[i];
				saver[i]->ancestor = new->buffer[new->size];
			}
			saver[sizer-1]->next = NULL;
			new->buffer[new->size]->ptrans = savel[0];
			savel[0]->nlink = new->buffer[new->size];
			for(i=1; i<sizel; i++) {
				savel[i-1]->pnext = savel[i];
				savel[i]->nlink = new->buffer[new->size];
			}
			savel[sizel-1]->pnext = NULL;
			new->size++;
			for(i=0; i<sizer; i++) {
				rightOff(saver[i]);
			}
			for(i=0; i<sizel; i++) {
				leftOff(savel[i]);
			}
		}
	}	
	for(ind=0; ind<new->size; ind++) {
		newOff(new->buffer[ind]);
	}
	for(ind=0; ind<cur->size; ind++) {
		doneOff(cur->buffer[ind]);
	}
	monfree((void*)savel);
	monfree((void*)saver);
}

/*Compute local graph of the given order of positions of sequence seq
by computing the suffix tree of seq with Ukkonen's algorithm*/
TypeSuffixTree *computeSuffixTree(TypeOneSequence seq) {
	TypeSuffixNode *s;
	TypePosition i, k, depth, offset;
	TypeSymbol t;
	TypeSuffixTree *suffixTree;
	
	suffixTree = (TypeSuffixTree*) monmalloc(sizeof(TypeSuffixTree));	
	suffixTree->bottom = newSuffixNode(0, -1);
	suffixTree->root = newSuffixNode(-1, -1);
	suffixTree->root->suffix = suffixTree->bottom;
	suffixTree->sequence = seq.sequence;
	s = suffixTree->root;
	k = 0;
	for(i=0; i<seq.size; i++) {
		update(s, k, i, &s, &k, suffixTree);
		canonize(s, k, i, &s, &k, suffixTree);
	}
	return suffixTree;
}

/*test whether or not a state with canonical reference pair (s(k,p)) is the end-point (has a t transition i positions later)
unchanged from Ukkonen's algorithm*/
int test_and_split(TypeSuffixNode *s, TypePosition k, TypePosition p, TypeSymbol t, TypeSuffixTree *suffixTree, TypeSuffixNode **sres) {
	if(k<=p) {
		TypeSuffixNode **prev, *s2;
		for(prev = &(s->trans); *prev != NULL && suffixTree->sequence[(*prev)->start] < suffixTree->sequence[k]; prev = &((*prev)->next))
		;
		s2 = *prev;
		if(suffixTree->sequence[s2->start+p-k+1] == t) {
			*sres = s;
			return 1;
		} else {
			*prev = newSuffixNode(s2->start, s2->start+p-k);
			(*prev)->next = s2->next;
			s2->next = NULL;
			s2->start += p-k+1;
			(*prev)->trans = s2;
			*sres = *prev;
			return 0;
		}
	} else {
		*sres = s;
		return (findTransition(s, t, suffixTree) != NULL);
	}		
}

/*return the canonical reference (sres,(kres,p)) of the state represented by the reference pair (s,(k,p))
unchanged from Ukkonen's algorithm*/
void canonize(TypeSuffixNode *s, TypePosition k, TypePosition p, TypeSuffixNode **sres, TypePosition *kres, TypeSuffixTree *suffixTree) {
	*sres = s;
	*kres = k;
	if(p>=k) {
		TypeSuffixNode *s2;
		s2 = findTransition(*sres, suffixTree->sequence[*kres], suffixTree);
		while(s2->end-s2->start <= p-*kres) {
			*kres += s2->end-s2->start+1;
			*sres = s2;
			if(*kres<=p)
				s2 = findTransition(*sres, suffixTree->sequence[*kres], suffixTree);
		}
	}
}

/*update the suffix suffixTree by adding the ith caracter of the sequence.
(s,(k,i-1)) is the canonical reference pair of the active point,
(sres, (kres,i)) return a reference pair of the next active point
suffixTree stores basic information relative to the sequence and the suffix suffixTree
depth is the depth of the active point (the length of the string spelled out from the root to this one)
order is the given order of decoding
link stores the local graph of this order
only changes from Ukkonen's algorithm are lines starting with /'star' 'star'/ 
concerning update of depth and local graph by adding links*/	
void update(TypeSuffixNode *s, TypePosition k, TypePosition i, TypeSuffixNode **sres, TypePosition *kres, TypeSuffixTree *suffixTree) {
	TypeSuffixNode *oldr, *r;
	int end_point;
	oldr = suffixTree->root;
	end_point = test_and_split(s, k, i-1, suffixTree->sequence[i], suffixTree, &r);
	while(!end_point) {
		addTransition(r, newSuffixNode(i, INFTY), suffixTree);
		if(oldr != suffixTree->root)
			oldr->suffix = r;
		oldr = r;
		canonize(s->suffix, k, i-1, &s, &k, suffixTree);
		end_point = test_and_split(s, k, i-1, suffixTree->sequence[i], suffixTree, &r);
	}
	if(oldr != suffixTree->root)
		oldr->suffix = s;
	*sres = s; *kres = k;
}


/*Basic utilities*/

/*return the son of node with label starting by symbol if it exists, NULL otherwise*/
TypeSuffixNode *findTransition(TypeSuffixNode *node, TypeSymbol symbol, TypeSuffixTree *suffixTree) {
	TypeSuffixNode *tmp;
	if(node == suffixTree->bottom)
		return suffixTree->root;
	for(tmp=node->trans; tmp != NULL && suffixTree->sequence[tmp->start] < symbol; tmp=tmp->next)
	;
	if(tmp != NULL && suffixTree->sequence[tmp->start] == symbol)
		return tmp;
	return NULL;		
}

/*insert toadd as son of node in lexicographic order - assume node and toadd aren't NULL*/
void addTransition(TypeSuffixNode *node, TypeSuffixNode *toadd, TypeSuffixTree *suffixTree) {
	TypeSuffixNode *cur;
	if(node->trans == NULL || suffixTree->sequence[node->trans->start]>suffixTree->sequence[toadd->start]) {
		toadd->next = node->trans;
		node->trans = toadd;
		return;
	}
	for(cur = node->trans; cur->next != NULL && suffixTree->sequence[cur->next->start]<suffixTree->sequence[toadd->start]; cur=cur->next)
	;
	toadd->next = cur->next;
	cur->next = toadd;
}	

/*allocated and fill a new link*/
TypeLink *newLink(TypePosition position, TypeLink *next) {
	TypeLink *link;
	link = (TypeLink*) monmalloc(sizeof(TypeLink));
	link->position = position;
	link->next = next;
	return link;
}

/*free a link*/
void freeLink(TypeLink *link) {
	if(link != NULL) {
		freeLink(link->next);
		monfree((void*)link);
	}
}

/*allocated and fill a new node*/
TypeSuffixNode *newSuffixNode(TypePosition start, TypePosition end) {
	TypeSuffixNode *node;
	node = (TypeSuffixNode*) monmalloc(sizeof(TypeSuffixNode));
	node->trans = NULL;
	node->next = NULL;
	node->suffix = NULL;
	node->start = start;
	node->end = end;
	return node;
}

/*free a node*/
void freeNode(TypeSuffixNode *node) {
	if(node != NULL) {
		TypeSuffixNode *tmp;
		tmp = node->trans;
		while(tmp != NULL) {
			TypeSuffixNode *retmp;
			retmp= tmp;
			tmp = tmp->next;
			freeNode(retmp);
		}
		monfree((void*)node);
	}
}
void freeSuffixTree(TypeSuffixTree *suffixTree) {
	freeNode(suffixTree->root);
	freeNode(suffixTree->bottom);
	monfree(suffixTree);
}
/*allocated and fill a new node*/
TypeLocalNode *newLocalNode(TypePosition pos) {
	TypeLocalNode *node;
	node = (TypeLocalNode*) monmalloc(sizeof(TypeLocalNode));
	node->position = pos;
	node->trans = NULL;
	node->next = NULL;
	node->ancestor = NULL;
	node->ptrans = NULL;
	node->pnext = NULL;
	node->nlink = NULL;
	node->state = 0;
	return node;
}

void freeLocalTree(TypeLocalTree *localTree) {
	TypePosition pos;
	
	for(pos=0; pos<localTree->size; pos++) {
		if(localTree->leaf[pos] != NULL)
			freeLocalNode(getAncestor(localTree->leaf[pos]), localTree);
	}
	monfree((void*)localTree->leaf);
	monfree((void*)localTree);
}
/*free a node*/
void freeLocalNode(TypeLocalNode *node, TypeLocalTree *localTree) {
	if(node != NULL) {
		TypeLocalNode *tmp;
		tmp = node->trans;
		if(tmp == NULL && localTree != NULL)
			localTree->leaf[node->position] = NULL;
		while(tmp != NULL) {
			TypeLocalNode *retmp;
			retmp= tmp;
			tmp = tmp->next;
			freeLocalNode(retmp, localTree);
		}
		monfree((void*)node);
	}
}

TypeLocalNode* getAncestor(TypeLocalNode* node) {
	if(node->ancestor != NULL)
		return getAncestor(node->ancestor);
	else
		return node;
}

void printLocalTree(FILE *f, TypeLocalTree *localTree) {
	TypePosition pos;
	int *done;
	
	done = (int*) monmalloc(localTree->size*sizeof(int));
	for(pos=0; pos<localTree->size; pos++) {
		done[pos] = 0;
	}
	for(pos=0; pos<localTree->size; pos++) {
		if(!done[pos]) {
			done[pos] = 1;
			printLocalNode(f, getAncestor(localTree->leaf[pos]), done, 0);
		}
	}
	monfree((void*) done);
}

/*print recursively the local node*/
void printLocalNode(FILE *f, TypeLocalNode* node, int *done, TypePosition depth) {
	TypePosition d = depth;
	TypeSymbol t;
	TypeLocalNode *tmp;
	
	if(node == NULL)
		return;
	if(depth>=0) {
		TypeLocalNode *bof;
		/*Print the branches coming from higher nodes.*/
		while(d>0) {
			fprintf(f, "|");
			d--;
		}
		fprintf(f, "+%ld ", node->position);
		if(node->ancestor == NULL) {
			fprintf(f, "prev ");
			for(tmp=node->ptrans; tmp!=NULL; tmp = tmp->pnext)
				fprintf(f, " %ld", tmp->position);
		}
			if(node->nlink != NULL)
				fprintf(f, " - next %ld", node->nlink->position);
			else
			 fprintf(f, " - next none");
		fprintf(f, "\n");
	}
	if(done != NULL && node->trans == NULL)
		done[node->position] = 1;
	for(tmp=node->trans; tmp!=NULL; tmp = tmp->next) {
		printLocalNode(f, tmp, done, depth+1);
	}
}

/*print the suffix tree*/
void printSuffixTree(FILE *f, TypeSuffixTree *suffixTree) {
	printSuffixNode(f, suffixTree->root, 0);
}

/*print recursively the suffix node*/
void printSuffixNode(FILE *f, TypeSuffixNode* node, TypePosition depth) {
	TypePosition d = depth;
	TypeSymbol t;
	TypeSuffixNode *tmp;
	
	if(node == NULL)
		return;
	if(depth>=0) {
		/*Print the branches coming from higher nodes.*/
		while(d>0) {
			fprintf(f, "|");
			d--;
		}
		fprintf(f, "+");
		fprintf(f, "[%ld, ", node->start);
		if(node->end != INFTY)
			fprintf(f, "%ld", node->end);
		else
			fprintf(f, "infty");
		fprintf(f, "]\n");
	}
	for(tmp=node->trans; tmp!=NULL; tmp = tmp->next) {
		printSuffixNode(f, tmp, depth+1);
	}
}

/*print all the sequences in a fasta file f*/
void printSequenceDebug(FILE *f, TypeOneSequence s, int sizeLine) {
		TypePosition position;
		for(position=0; position<s.size; position += sizeLine) {
			TypePosition j, maxLine;
			maxLine = position+sizeLine;
			maxLine = maxLine>s.size?s.size:maxLine;
			for(j=position; j<maxLine; j++)
				fprintf(f, " %ld", j);
			fprintf(f, "\n");
			for(j=position; j<maxLine; j++) {
				if(j>9)
					fprintf(f, " ");
				fprintf(f, " %ld", s.sequence[j]);
			}
			fprintf(f, "\n\n");
		}
		fprintf(f, "\n\n");
}



int isLeft(TypeLocalNode *node) {
	return (node->state & LEFT);
}

void leftOn(TypeLocalNode *node) {
	node->state = node->state | LEFT;
}

void leftOff(TypeLocalNode *node) {
	node->state = node->state & ~LEFT;
}

int isRight(TypeLocalNode *node) {
	return (node->state & RIGHT);
}

void rightOn(TypeLocalNode *node) {
	node->state = node->state | RIGHT;
}

void rightOff(TypeLocalNode *node) {
	node->state = node->state & ~RIGHT;
}
int isDone(TypeLocalNode *node) {
	return (node->state & DONE);
}

void doneOn(TypeLocalNode *node) {
	node->state = node->state | DONE;
}

void doneOff(TypeLocalNode *node) {
	node->state = node->state & ~DONE;
}
int isNew(TypeLocalNode *node) {
	return (node->state & NEW);
}

void newOn(TypeLocalNode *node) {
	node->state = node->state | NEW;
}

void newOff(TypeLocalNode *node) {
	node->state = node->state & ~NEW;
}

