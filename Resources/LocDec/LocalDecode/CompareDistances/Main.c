#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "Distance.h"
#include "Utils.h"


#define SIZE_BUFFER_CHAR 300
#define MAX_FUNC 9

#define HELPMESSAGE "\nusage: compare [options] <input file1> <input file2>\n\nreturns the name of input file name followed by the Pearson coefficient\n\nThe input files have to be in Nexus format.\n\nOptions:\n\n-h                  output help\n\n-o <output file>    write a table with in each line values for dist1 and dist2\n\n-f <output file>    write distances matrix dist1-dist2\n\n-f <f, s, m, x>     select the format of output for distances matrix\n                    (r: raw, t: text table, p: phylip, n: nexus)\n"

int main(int argc, char **argv) {		
	char option[256], *name, inputFileName1[SIZE_BUFFER_CHAR], inputFileName2[SIZE_BUFFER_CHAR], 
	outputFileName[SIZE_BUFFER_CHAR],outputFileCorrName[SIZE_BUFFER_CHAR], bufferOutput[SIZE_BUFFER_CHAR], *table, 
	outputFormat = 'r', typeOut = 'c';
	int output = 0, outputCorr = 0, i;
	FILE *fi1, *fi2, *fo;

	for(i=0; i<256; i++)
		option[i] = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		int j;
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['o']) {
			option['o'] = 0;
			output = 1;
			if((i+1)<argc && sscanf(argv[i+1], "%s", outputFileName) == 1)
				i++;
		}
		if(option['c']) {
			option['c'] = 0;
			outputCorr = 1;
			if((i+1)<argc && sscanf(argv[i+1], "%s", outputFileCorrName) == 1)
				i++;
		}
		if(option['f']) {
			option['f'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &outputFormat) == 1)
				i++;
		}
		if(option['h']) {
			printf("%s\n", HELPMESSAGE);
			exitProg(ExitOk, NULL);
		}
	}
	if (i>=argc || sscanf(argv[i++], "%s", inputFileName1) != 1) exitProg(ErrorArgument, HELPMESSAGE);
	if (i>=argc || sscanf(argv[i++], "%s", inputFileName2) != 1) exitProg(ErrorArgument, HELPMESSAGE);

	if((fi1 = fopen(inputFileName1, "r")) && (fi2 = fopen(inputFileName2, "r"))) {
		TypeDistance dist1, dist2;
		TypeCoeffReg corr;
		char *name, *tmp;
		dist1 = readDistancesNexus(fi1);
		dist2 = readDistancesNexus(fi2);
		corr = computeCorrelation(dist1, dist2);
		name = strrchr(inputFileName1, '/');
		if(name == NULL)
			name = inputFileName1;
		else
			name++;
		tmp = strrchr(name, '-');
		if(tmp != NULL)
			*tmp = '\0';
		printf("%f\n", corr.r);
		if(output && (fo = fopen(outputFileName, "w"))) {
			printPoint(fo, dist1, dist2);
			fclose(fo);
		}
		if(outputCorr && (fo = fopen(outputFileCorrName, "w"))) {
			TypeDistance dist = diffCorr(dist1, dist2);
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
		}
	} else {
		exitProg(ErrorReading, inputFileName1);
	}

	exitProg(ExitOk,NULL);
}
