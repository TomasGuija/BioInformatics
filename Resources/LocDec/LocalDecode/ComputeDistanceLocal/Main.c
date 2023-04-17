#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "Sequences.h"
#include "Distance.h"
#include "Decode.h"
#include "Utils.h"
#include "DistSeq.h"


#define SIZE_BUFFER_CHAR 300
#define MAX_FUNC 2

#define HELPMESSAGE "\nusage: local [options] <input file> <output file>\n\nThe input file has to be in Fasta format.\n\nOptions:\n\n-h                  output help\n\n-f <f, s, m, x>     select the format of output for distances matrix\n                    (r: raw, t: text table, p: phylip, n: nexus)\n\n-s <d, r, p>        select the type of alphabet (d: DNA, r: RNA, p: Protein)\n                    maybe useful to handle ambiguity characters\n\n-d <l, s>           select the type of decoding:\n                     * l: use local decoding\n                     * s: use sequence of sliding blocks\n\n-m <number>           select the type of dissimilarity:\n                     * 0: trivial dist\n                     * 1: Pham\n-o <number>         set the order of local decoding/length of the blocks\n"
#define AMBIMESSAGE "The ambiguity codes are for DNA and RNA:\n                                         * Y: T/U or C\n                                         * R: A or G\n                                         * M: A or C\n                                         * K: T/U or G\n                                         * W: T/U or A\n                                         * S: C or G\n                                         * B: T/U, C, or G\n                                         * D: T/U, A, or G\n                                         * H: T/U, C, or A\n                                         * V: C, A, or G\n                                         * N: T/U, C, A, or G\n\nand for protein:\n                                         * B: D or N\n                                         * Q or E\n                                         * X: Unidentified\n"
static double meanLength(TypeSetOfSequences set);

double meanLength(TypeSetOfSequences set) {
	TypePosition sum = 0;
	TypeNumber n;
	for(n=0; n<set.number; n++)
		sum += set.size[n];
	return ((double) sum)/((double)set.number);
}

int main(int argc, char **argv) {		
	TypePosition order=10;
	char option[256], inputFileName[SIZE_BUFFER_CHAR], outputFileName[SIZE_BUFFER_CHAR], bufferOutput[SIZE_BUFFER_CHAR], *table, 
	outputFormat = 'r', typeDec = 'l', typeAlphabet = '?';
	TypeSetOfSequences set;
	TypeDistFunction *distfunc[MAX_FUNC]={computeTrivial, computePham};
	FILE *fi, *fo;
	int i, typeDist = 0;

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
		if(option['d']) {
			option['d'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &typeDec) == 1)
				i++;
		}
		if(option['m']) {
			option['m'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &typeDist) == 1)
				i++;
			if(typeDist >= MAX_FUNC)
				typeDist = 0;
		}
		if(option['h']) {
			printf("%s\n", HELPMESSAGE);
			exitProg(ExitOk, NULL);
		}
		if(option['o']) {
			option['o'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%ld", &order) == 1)
				i++;
		}
	}
	if (i>=argc || sscanf(argv[i++], "%s", inputFileName) != 1) exitProg(ErrorArgument, HELPMESSAGE);
	if (i>=argc || sscanf(argv[i++], "%s", outputFileName) != 1) exitProg(ErrorArgument, HELPMESSAGE);
	
	if(typeDist == 1)
		order = 1;

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
		set = readSequencesFasta(fi, table, typeAlphabet == '?');
		switch(typeAlphabet) {
			case 'd':
			case 'r':
				set.ambiguity = getXNAAmbiguity();
				break;
			case 'p':
				set.ambiguity = getProteinAmbiguity();
				break;
			case '?':
			default:
				set.ambiguity.number = 0;
		}
		set.cardinal -= set.ambiguity.number;
		fclose(fi);
	} else {
		exitProg(ErrorReading, inputFileName);
	}
	
	if(fo = fopen(outputFileName, "w")) {
		TypeSetOfSequences dec;
		TypeDistance dist;
		dist = computeWholeDistancePair(set,  order, distfunc[typeDist], typeDec=='l');
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
			default:
				printDistanceNexus(fo, dist);
		}
		fclose(fo);
	} else {
		exitProg(ErrorWriting, outputFileName);
	}
	exitProg(ExitOk,NULL);
	return 0;
}
