#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "Alignement.h"
#include "Distance.h"
#include "Utils.h"
#include "DistAln.h"


#define SIZE_BUFFER_CHAR 300
#define HELPMESSAGE "usage: align [options] <input file> <output file>\n\nThe input file has to be in Clustal format.\n\nOptions:\n\n-h                  output help\n\n-f <f, s, m, x>     select the format of output for distances matrix\n                    (r: raw, t: text table, p: phylip, n: nexus)\n\n-s <d, r, p>        select the type of alphabet (d: DNA, r: RNA, p: Protein)\n                    maybe useful to handle ambiguity characters\n"
#define AMBIMESSAGE "The ambiguity codes are for DNA and RNA:\n                                         * Y: T/U or C\n                                         * R: A or G\n                                         * M: A or C\n                                         * K: T/U or G\n                                         * W: T/U or A\n                                         * S: C or G\n                                         * B: T/U, C, or G\n                                         * D: T/U, A, or G\n                                         * H: T/U, C, or A\n                                         * V: C, A, or G\n                                         * N: T/U, C, A, or G\n\nand for protein:\n                                         * B: D or N\n                                         * Q or E\n                                         * X: Unidentified\n"

int main(int argc, char **argv) {		
	char option[256], inputFileName[SIZE_BUFFER_CHAR], outputFileName[SIZE_BUFFER_CHAR], bufferOutput[SIZE_BUFFER_CHAR], *table, 
	outputFormat = 'r', typeAlphabet = '?';
	TypeAlignment aln;
		
	FILE *fi, *fo;
	int i = 1, typeDist = 0;

	for(i=0; i<256; i++)
		option[i] = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		int j;
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['f']) {
			option['f'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &outputFormat) == 1)
				i++;
		}
		if(option['s']) {
			option['s'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &typeAlphabet) == 1)
				i++;
		}
		if(option['h']) {
			printf("%s\n", HELPMESSAGE);
			exitProg(ExitOk, NULL);
		}
	}
	if (i>=argc || sscanf(argv[i++], "%s", inputFileName) != 1) exitProg(ErrorArgument, HELPMESSAGE);
	if (i>=argc || sscanf(argv[i++], "%s", outputFileName) != 1) exitProg(ErrorArgument, HELPMESSAGE);

	switch(typeAlphabet) {
		case 'd':
			table = (char*) monmalloc((strlen(DNA)+1)*sizeof(char));
			strcpy(table, DNA);
			break;
		case 'r':
			table = (char*) monmalloc((strlen(RNA)+1)*sizeof(char));
			strcpy(table, RNA);
			break;
		case 'p':
			table = (char*) monmalloc((strlen(PRO)+1)*sizeof(char));
			strcpy(table, PRO);
			break;
		case '?':
		default:
			table = (char*) monmalloc(sizeof(char));
			table[0] = '\0';
	}
	if(fi = fopen(inputFileName, "r")) {
		aln = readAlignementClustal(fi, table, typeAlphabet == '?');
		switch(typeAlphabet) {
			case 'd':
			case 'r':
				aln.ambiguity = getXNAAmbiguity();
				break;
			case 'p':
				aln.ambiguity = getProteinAmbiguity();
				break;
			case '?':
			default:
				aln.ambiguity.number = 0;
		}
		aln.cardinal -= aln.ambiguity.number;
		fclose(fi);
	} else {
		exitProg(ErrorReading, inputFileName);
	}
	if(fo = fopen(outputFileName, "w")) {
		TypeDistance dist;
		dist = computeWholeDistancePair(aln, computeNorm1);
		switch(outputFormat) {
			case 't':
				printDistanceTable(fo, dist);
				break;
			case 'r':
				printDistanceRaw(fo, dist);
				break;
			case 'p':
				printDistancePhylip(fo, dist);
				break;
			case 'n':
				printDistanceNexus(fo, dist);
		}
		fclose(fo);
	} else {
		exitProg(ErrorWriting, outputFileName);
	}
	exitProg(ExitOk,"Bye");
	return 0;
}
