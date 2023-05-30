#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "Sequences.h"
#include "Utils.h"


#define MARK '>'
#define INC_BUFFER_SEQ 10
#define INC_BUFFER_SYMBOL 500
#define INC_BUFFER_ANC 10
#define INC_BUFFER_TMP 30
#define SIZE_MIN 1
static void fixNames(char **name, TypeNumber N);

/*return the single sequence wich is the concatenation of all the sequences in
set, each followed by a special caracter (all differents) Fix all the ambiguity character 
as unknown and all differents*/

TypeOneSequence toOne(TypeSetOfSequences set) {
	TypeOneSequence one;
	TypePosition totalSize, ind;
	TypeNumber n;
	
	one.cardinal = set.cardinal;
	one.size = set.number;
	for(n=0; n<set.number; n++)
		one.size += set.size[n];
	one.sequence = (TypeSymbol*) monmalloc(one.size*sizeof(TypeSymbol));
	ind = 0;
	for(n=0; n<set.number; n++) {
		TypePosition p;
		for(p=0; p<set.size[n]; p++)
			if(set.sequence[n][p]<set.cardinal)
				one.sequence[ind++] = set.sequence[n][p];
			else
				one.sequence[ind++] = one.cardinal++;
		one.sequence[ind++] = one.cardinal++;
	}
	return one;
}

/*read all the sequences from a fasta file f*/
TypeSetOfSequences readSequencesFasta(FILE *f, char *table, int canInc) {
	char c;
	int i, code[256];
	TypeNumber sizeBufferSeq;
	TypeSetOfSequences seq;

	seq.cardinal = 0;
	for(i=0; i<256; i++)
		code[i] = -1;
	if(table != NULL) {		
		for(i=0; table[i] != '\0'; i++) {
			code[table[i]] = i;
			seq.cardinal++;
		}
	}
	seq.number = 0;
	sizeBufferSeq = INC_BUFFER_SEQ;
	seq.sequence = (TypeSymbol **) monmalloc(sizeBufferSeq*sizeof(TypeSymbol*));
	seq.size = (TypePosition *) monmalloc(sizeBufferSeq*sizeof(TypePosition));
	seq.name = (char **) monmalloc(sizeBufferSeq*sizeof(char *));
	do {
		c = getc(f);
	} while(c != EOF && c != MARK);
	while(c != EOF) {
		TypeSymbol *bufferSymbol;
		char *bufferName;
		TypePosition sizebufferName = INC_BUFFER_TMP, sizeTmp = 0, sizeBufferSymbol, position;
		bufferName = (char *) monmalloc(sizebufferName*sizeof(char));
		c = getc(f);
		while(c != EOF && !IsLineSeparator(c)) {
			if(sizeTmp>=sizebufferName) {
				sizebufferName += INC_BUFFER_TMP;
				bufferName = (char *) monrealloc((void *) bufferName, sizebufferName*sizeof(char));
			}
			bufferName[sizeTmp++] = c;
			c = getc(f);
		}
		if(sizeTmp>=sizebufferName) {
			sizebufferName += INC_BUFFER_TMP;
			bufferName = (char *) monrealloc((void *) bufferName, sizebufferName*sizeof(char));
		}
		bufferName[sizeTmp++] ='\0';
		while(IsLineSeparator(c)) {
			c = getc(f);
		}
		sizeBufferSymbol = INC_BUFFER_SYMBOL;
		bufferSymbol = (TypeSymbol*) monmalloc(sizeBufferSymbol*sizeof(TypeSymbol));
		position = 0;
		while(c != EOF && c != MARK) {
			TypeSymbol symbol;
			if(isalpha(c)) {
				if(code[(unsigned char) toupper(c)] == -1) {
					if(table == NULL || canInc) {
						code[(unsigned char) toupper(c)] = seq.cardinal++;
					} else {
						char *msg = (char*) monmalloc(100*sizeof(char));
						sprintf(msg, "Model and TypeSetOfSequences aren't compatible (symbol %c not in the model)", toupper(c));
						exitProg(ErrorInit, msg);
					}
				}
				symbol = code[(unsigned char) toupper(c)];
				if(position>=sizeBufferSymbol) {
					sizeBufferSymbol += INC_BUFFER_SYMBOL;
					bufferSymbol = (TypeSymbol*) monrealloc(bufferSymbol, sizeBufferSymbol*sizeof(TypeSymbol));
				}
				bufferSymbol[position++] = symbol;
			}
			c = getc(f);
		}
		if(position>=SIZE_MIN) {
			if(seq.number>=sizeBufferSeq) {
				sizeBufferSeq += INC_BUFFER_SEQ;
				seq.sequence = (TypeSymbol **) monrealloc(seq.sequence, sizeBufferSeq*sizeof(TypeSymbol*));
				seq.size = (TypePosition *) monrealloc(seq.size, sizeBufferSeq*sizeof(TypePosition));
				seq.name = (char **) monrealloc(seq.name, sizeBufferSeq*sizeof(char *));
			}
			seq.sequence[seq.number] = (TypeSymbol*) monrealloc(bufferSymbol, position*sizeof(TypeSymbol));
			seq.name[seq.number] = (char *) monrealloc((void *) bufferName, sizeTmp*sizeof(char));
			seq.size[seq.number] = position;
			seq.number++;
		}
	}
	seq.sequence = (TypeSymbol **) monrealloc(seq.sequence, seq.number*sizeof(TypeSymbol*));
	seq.size = (TypePosition *) monrealloc(seq.size, seq.number*sizeof(TypePosition));
	seq.name = (char **) monrealloc(seq.name, seq.number*sizeof(char *));
	if(table == NULL || canInc) {
		if(table == NULL) {
			seq.table = (char *) monmalloc((seq.cardinal+1)*sizeof(char *));
		} else {
			seq.table = (char *) monrealloc(table, (seq.cardinal+1)*sizeof(char *));
		}
		for(i=0; i<256; i++)
			if(code[i] != -1)
				seq.table[code[i]] = (char) i;
		seq.table[seq.cardinal] = '\0';
	} else {
		seq.table = table;
	}
	fixNames(seq.name, seq.number);
	return seq;
}

void fixNames(char **name, TypeNumber N) {
	TypeNumber n;
	int i;
	
	for(n=0; n<N; n++)
		for(i=0; name[n][i]!='\0'; i++)
			if(isspace(name[n][i]))
				name[n][i] = '_';
}
/*print all the sequences in a fasta file f*/
void printSequencesFasta(FILE *f, TypeSetOfSequences s, int sizeLine) {
	TypeNumber n;
	for(n=0; n<s.number; n++) {
		TypePosition position;
		if(s.name != NULL)
			fprintf(f, ">%s\n", s.name[n]);
		else
			fprintf(f, ">Sequence %d\n", n+1);
		for(position=0; position<s.size[n]; position += sizeLine) {
			TypePosition j, maxLine;
			maxLine = position+sizeLine;
			maxLine = maxLine>s.size[n]?s.size[n]:maxLine;
			if(s.table != NULL) {
				for(j=position; j<maxLine; j++)
					fprintf(f, "%c", s.table[s.sequence[n][j]]);
			} else {
				for(j=position; j<maxLine; j++)
					fprintf(f, "%d", s.sequence[n][j]);
			}
			fprintf(f, "\n");
		}
		fprintf(f, "\n\n");
	}
}

/*read all the sequences from a raw file f*/
TypeSetOfSequences readSequencesRaw(FILE *f) {
	char c;
	int i;
	TypeNumber sizeBufferSeq;

	TypeSetOfSequences seq;

	seq.cardinal = 0; seq.number = 0;
	sizeBufferSeq = INC_BUFFER_SEQ;
	seq.sequence = (TypeSymbol **) monmalloc(sizeBufferSeq*sizeof(TypeSymbol*));
	seq.size = (TypePosition *) monmalloc(sizeBufferSeq*sizeof(TypePosition));
	c = getc(f);
	while(c != EOF) {
		TypeSymbol *bufferSymbol;
		char *bufferName;
		TypePosition sizeBufferSymbol, position;
		sizeBufferSymbol = INC_BUFFER_SYMBOL;
		bufferSymbol = (TypeSymbol*) monmalloc(sizeBufferSymbol*sizeof(TypeSymbol));
		position = 0;
		while(c != EOF && !((c >= '0') && (c <= '9'))) {
			c = getc(f);
		}
		while(c != EOF && !IsLineSeparator(c)) {
			TypeSymbol symbol = c - '0';
			while(((c=getc(f)) != EOF) && (c >= '0') && (c <= '9'))
				symbol = symbol*10+c-'0';
			if(symbol>seq.cardinal)
				seq.cardinal = symbol;
			if(position>=sizeBufferSymbol) {
				sizeBufferSymbol += INC_BUFFER_SYMBOL;
				bufferSymbol = (TypeSymbol*) monrealloc(bufferSymbol, sizeBufferSymbol*sizeof(TypeSymbol));
			}
			bufferSymbol[position++] = symbol;
			while(c != EOF && !IsLineSeparator(c) && !((c >= '0') && (c <= '9'))) {
				c = getc(f);
			}
		}
		if(position>=SIZE_MIN) {
			if(seq.number>=sizeBufferSeq) {
				sizeBufferSeq += INC_BUFFER_SEQ;
				seq.sequence = (TypeSymbol **) monrealloc(seq.sequence, sizeBufferSeq*sizeof(TypeSymbol*));
				seq.size = (TypePosition *) monrealloc(seq.size, sizeBufferSeq*sizeof(TypePosition));
			}
			seq.sequence[seq.number] = (TypeSymbol*) monrealloc(bufferSymbol, position*sizeof(TypeSymbol));
			seq.size[seq.number] = position;
			seq.number++;
		}
	}
	seq.sequence = (TypeSymbol **) monrealloc(seq.sequence, seq.number*sizeof(TypeSymbol*));
	seq.size = (TypePosition *) monrealloc(seq.size, seq.number*sizeof(TypePosition));
	seq.cardinal++;
	seq.name = NULL;
	seq.table = NULL;
	return seq;
}


/*free "seq"*/
void freeSequences(TypeSetOfSequences seq) {
	TypeNumber n;
	monfree((void *) seq.size);
	for(n=0; n<seq.number; n++)
		monfree((void *) seq.sequence[n]);
	monfree((void *) seq.sequence);
	if(seq.name != NULL) {
		for(n=0; n<seq.number; n++)
			monfree((void *) seq.name[n]);
		monfree((void *) seq.name);
	}
	if(seq.table != NULL)
		monfree((void *) seq.table);
}

/*write the sequence "seq" in "f"*/
void writeSequences(FILE *f, TypeSetOfSequences seq) {
	TypeNumber n;
	for(n=0; n<seq.number; n++) {
		TypePosition position;
		for(position=0; position<seq.size[n]; position++)
			fprintf(f, "%ld ", seq.sequence[n][position]);
		fprintf(f, "\n");
	}
}

/*return a "hard-copy" of "s"*/
TypeSetOfSequences clone(TypeSetOfSequences seq) {
	TypePosition i;
	TypeNumber n;
	TypeSetOfSequences res;
	res.name = seq.name;
	res.table = seq.table;
	res.number = seq.number;
	res.cardinal = seq.cardinal;
	res.size = (TypePosition*) monmalloc(res.number*sizeof(TypePosition));
	res.sequence = (TypeSymbol**) monmalloc(res.number*sizeof(TypeSymbol*));
	for(n=0; n<res.number; n++) {
		res.size[n] = seq.size[n];
		res.sequence[n] = (TypeSymbol*) monmalloc(res.size[n]*sizeof(TypeSymbol));
		for(i=0; i<res.size[n]; i++)
			res.sequence[i] = seq.sequence[i];
	}
	return res;
}

TypeAmbiguity getXNAAmbiguity() {
	TypeAmbiguity ambi;
	static TypeSymbol 
		tmp0[2] = {3, 1},
		tmp1[2] = {0, 2},
		tmp2[2] = {0, 1},
		tmp3[2] = {3, 2},
		tmp4[2] = {3, 0},
		tmp5[2] = {1, 2},
		tmp6[3] = {3, 1, 2},
		tmp7[3] = {3, 0, 2},
		tmp8[3] = {3, 1, 0},
		tmp9[3] = {1, 0, 2},
		tmp10[4] = {3, 1, 0, 2};
	ambi.number = 11;
	ambi.size = (TypeSymbol*) monmalloc(ambi.number*sizeof(TypeSymbol));
	ambi.set = (TypeSymbol**) monmalloc(ambi.number*sizeof(TypeSymbol*));
	ambi.size[0] = 2;
	ambi.size[1] = 2;
	ambi.size[2] = 2;
	ambi.size[3] = 2;
	ambi.size[4] = 2;
	ambi.size[5] = 2;
	ambi.size[6] = 3;
	ambi.size[7] = 3;
	ambi.size[8] = 3;
	ambi.size[9] = 3;
	ambi.size[10] = 4;
	ambi.set[0] = tmp0;
	ambi.set[1] = tmp1;
	ambi.set[2] = tmp2;
	ambi.set[3] = tmp3;
	ambi.set[4] = tmp4;
	ambi.set[5] = tmp5;
	ambi.set[6] = tmp6;
	ambi.set[7] = tmp7;
	ambi.set[8] = tmp8;
	ambi.set[9] = tmp9;
	ambi.set[10] = tmp10;
	return ambi;
}

TypeAmbiguity getProteinAmbiguity() {
	TypeAmbiguity ambi;
	static TypeSymbol 
			tmp0[2] = {0, 3},
			tmp1[2] = {1, 4},
			tmp2[20] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
	ambi.number = 3;
	ambi.size = (TypeSymbol*) monmalloc(ambi.number*sizeof(TypeSymbol));
	ambi.set = (TypeSymbol**) monmalloc(ambi.number*sizeof(TypeSymbol*));
	ambi.size[0] = 2;
	ambi.size[1] = 2;
	ambi.size[2] = 20;
	ambi.set[0] = tmp0;
	ambi.set[1] = tmp1;
	ambi.set[2] = tmp2;
	return ambi;
}
//             01234567890123456789
//#define PRO "DEGNQCSTYAVLIPFMWKRHBZX" /*"DEGNQCSTYAVLIPFMWKRH" + "X"*/
