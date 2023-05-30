
#ifndef SequencesF
#define SequencesF

#include <stdlib.h>
#include <stdio.h>

typedef /*unsigned*/ long TypePosition;
typedef unsigned int TypeSymbol;
typedef unsigned int TypeNumber;
/*Structure storing a set of sequence*/
typedef struct AMBIGUITY {
	char *car;
	TypeSymbol **set, *size, number;
} TypeAmbiguity;

/*Structure storing one sequence*/
typedef struct ONESEQUENCE {
/*cardinal of the alphabet*/
	TypePosition size;
	TypeSymbol cardinal;
	TypeSymbol *sequence;
	TypeAmbiguity ambiguity;
} TypeOneSequence;

/*Structure storing a set of sequence*/
typedef struct SEQUENCES {
	char **name, *table, *title;
/*number of TypeSetOfSequences include in the set*/
	TypeNumber number;
/*cardinal of the alphabet*/
	TypeSymbol cardinal;
/*size[i] size of the sequence i*/
	TypePosition *size;
/*table of all the sequences*/
	TypeSymbol **sequence;
	TypeAmbiguity ambiguity;
} TypeSetOfSequences;

TypeOneSequence toOne(TypeSetOfSequences set);
void printSequencesFasta(FILE *f, TypeSetOfSequences s, int sizeLine);
TypeSetOfSequences readSequencesFasta(FILE *f, char *table, int canInc);
TypeSetOfSequences readSequencesRaw(FILE *f);
void writeSequences(FILE *f, TypeSetOfSequences seq);
void freeSequences(TypeSetOfSequences seq);
TypeAmbiguity getXNAAmbiguity();
TypeAmbiguity getProteinAmbiguity();
TypeSetOfSequences clone(TypeSetOfSequences s);

#define SEP -1
#define FIN -1
#define DNA "ACGTYRMKWSBDHVN" /*"ACGT" + "YRMKWSBDHVN"*/
#define RNA "ACGUYRMKWSBDHVN" /*"ACGU" + "YRMKWSBDHVN"*/
#define PRO "DEGNQCSTYAVLIPFMWKRHBZX" /*"DEGNQCSTYAVLIPFMWKRH" + "X"*/

#endif
