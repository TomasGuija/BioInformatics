#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "Utils.h"
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


void printPreambleNexus(FILE *f, char **name, int number) {
	int i;
	fprintf(f, "#NEXUS\n\n");
	fprintf(f, "BEGIN taxa;\n");
	fprintf(f, "    DIMENSIONS ntax=%ld;\n", number);
	fprintf(f, "TAXLABELS\n");
	for(i=0; i<number; i++)
		fprintf(f, "'%s'\n", name[i]);
	fprintf(f, ";\n");
	fprintf(f, "END;\n");
}

double log2(double x) {
	return log(x)/log(2);
}

void exitProg(TypeExit code, char *message) {
	switch(code)
	{
		case ErrorArgument:
			printf("Error when reading arguments\n");
			break;
		case ErrorInit:
			printf("Problem of initialization\n");
			break;
		case ErrorReading:
			printf("Problem when reading a file\n");
			break;
		case ErrorWriting:
			printf("Problem when writing a file\n");
			break;
		case ErrorMemory:
			printf("Not enougth memory\n");
			break;
		case ErrorExec:
			printf("Problem during execution...\n");
			break;
	}
	if(message != NULL)
		printf("%s\n", message);
	if(code == ExitOk)
		exit(EXIT_SUCCESS);
	else
		exit(EXIT_FAILURE);
}

/****************************************************/
/************** Allocations ********************/
/****************************************************/
void *monmalloc(long size) {
	void *point;
	if(size<0 || !(point = malloc(size))) {
		char tmp[200];
		sprintf(tmp, "Try to allocate %ld bytes\n", size);
		exitProg(ErrorMemory, tmp);
	}
	return point;
} 

void *monrealloc(void *in, long size) {
	void *point;
	if(size<0 || !(point = realloc(in, size))) {
		char tmp[200];
		sprintf(tmp, "Try to reallocate %ld bytes\n", size);
		exitProg(ErrorMemory, tmp);
	}
	return point;
}

void monfree(void *in) {
	if(in != 0L)
		free(in);
}
int IsSeparator(char c) {
	return (c == ' ' || c == '\t' || c == '\n' || c == '\r');
}
int IsItemSeparator(char c) {
	return (c == ' ' || c == '\t' || c == ';');
}
int IsLineSeparator(char c) {
	return (c == '\n' || c == '\r');
}

int readLine(FILE *f, char *buffer) {
	char c;
	int tot = 0;
	while((c=getc(f))!=EOF && IsLineSeparator(c));
	while(c!=EOF && !IsLineSeparator(c)) {
		buffer[tot++] = c;
		c=getc(f);
	}
	if(c!=EOF)
		buffer[tot++] = '\n';
	buffer[tot] = '\0';
	return tot;
}

int readItem(FILE *f, char *buffer) {
	char c;
	int tot = 0;
	while((c=getc(f))!=EOF && IsSeparator(c));
	while(c!=EOF && !IsSeparator(c)) {
		buffer[tot++] = c;
		c=getc(f);
	}
	buffer[tot++] = '\n';
	buffer[tot] = '\0';
	return tot;
}

char *truncFileName(char* name) {
	int i;
	for(i=strlen(name)-1; i>=0 && name[i] != '.'; i--);
	if(i>=0)
	 name[i] = '\0';
	return name;
}

char *getExtension(char* name) {
	int i, length;
	
	length = strlen(name);
	for(i=length-1; i>0 && name[i] != '.'; i--);
	if(i>0)
		return &(name[i+1]);
	else
		return &(name[length]);
}

int tokenize(char *src, char **dest) {
	int indice = 0, pos = 0;
	do {
		while(src[pos] != '\0' && IsSeparator(src[pos]))
			 pos++;
		if(src[pos] != '\0') {
			dest[indice++] = &(src[pos]);
		}
		while(src[pos] != '\0' && !IsSeparator(src[pos]))
			 pos++;
		if(src[pos] != '\0') {
			src[pos++] = '\0';
		}
	} while(src[pos] != '\0');
	dest[indice] = NULL;
	return indice;
}

int find(char *src, char **dest, int size) {
	int i;
	for(i=0; i<size && strcmp(src, dest[i]) != 0; i++);
	return i;
}

void fixSpace(char *src) {
	int i;
	for(i=0; src[i] != '\0'; i++)
		if(IsSeparator(src[i]))
			src[i] = '_';
}
