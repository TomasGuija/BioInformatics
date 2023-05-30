
#ifndef AlignementF
#define AlignementF

#include <stdlib.h>
#include <stdio.h>
#include "Utils.h"

#define EMPTY 255

typedef enum TAF {
	fasta=0,
	clustal,
	msf,
	markx,
	srs,
	unknown
} TypeAlignmentFile;

typedef struct COMMENT {
	int number;
	char **text;
} TypeComment;

typedef struct ALIGNMENT {
	TypeNumber number;
	TypePosition size;
	TypeSymbol **sequence, empty, cardinal;
	TypeComment comment;
	char **name, *table, *title;
	TypeAmbiguity ambiguity;
} TypeAlignment;


TypeAlignment readAlignement(FILE *f, char *table, int canInc);
TypeAlignment readAlignementFasta(FILE *f, char *table, int canInc);
TypeAlignment readAlignementMsf(FILE *f, char *table, int canInc);
TypeAlignment readAlignementClustal(FILE *f, char *table, int canInc);
void printAlignmentFasta(FILE *f, TypeAlignment al, int sizeLine);
void printAlignmentMsf(FILE *f, TypeAlignment al, int sizeLine);
void printAlignmentSrs(FILE *f, TypeAlignment al, int sizeLine);
void printAlignmentMarkX(FILE *f, TypeAlignment al, int sizeLine);
void printHeadPair(FILE *f, TypeAlignment al);
void printHeadMulti(FILE *f, TypeAlignment al);
#endif
