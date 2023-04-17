#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "Utils.h"
#include "Alignement.h"

#define MARK '>'
#define INC_BUFFER_SEQ 10
#define INC_BUFFER_SYMBOL 10000
#define INC_BUFFER_ANC 10
#define INC_BUFFER_TMP 30
#define SIZE_MIN 1

#define MARGIN 8
#define LEFT 15
#define RIGHT 6
#define INC 10
#define EMPTY_CAR '-'

typedef TypeAlignment TypeReaderAlign(FILE*, char*, int);

static int isempty(char c) {
	return c == '-' || c =='_' || c == '.';
}



/*Read Alignement Clustal*/
TypeAlignmentFile getAlignmentFileType(FILE *f) {
	TypeReaderAlign reader;
	int ntoken, sizeline;
	char *buffer, **token;
	TypeAlignmentFile type = fasta;

	buffer = (char*) monmalloc(INC_BUFFER_SYMBOL*sizeof(char));
	token = (char**) monmalloc(INC_BUFFER_SYMBOL*sizeof(char*));

	while((sizeline = readLine(f, buffer)) != 0 && buffer[0] == '#');
	if(sizeline == 0 || (ntoken=tokenize(buffer, token)) == 0)
		return unknown;
	if(token[0][0] == '>')
		type = fasta;
	if(strcmp(token[0], "CLUSTAL") == 0)
		type = clustal;
	if(strcmp(token[0], "PileUp") == 0 ||
		strcmp(token[0], "!!AA_MULTIPLE_ALIGNMENT") == 0 ||
		strcmp(token[0], "!!NA_MULTIPLE_ALIGNMENT") == 0 ||
		strcmp(token[0], "MSF") == 0)
		type = msf;
	rewind(f);
	monfree((void*) buffer);
	monfree((void*) token);
	return type;
}

TypeAlignment readAlignement(FILE *f, char *table, int canInc) {
	switch(getAlignmentFileType(f)) {
		case clustal:
			return readAlignementClustal(f, table, canInc);
			break;
		case msf:
			return readAlignementMsf(f, table, canInc);
			break;
		case fasta:
		default:
			return readAlignementFasta(f, table, canInc);
			break;
	}
}

/*Read Alignement FASTA*/
TypeAlignment readAlignementFasta(FILE *f, char *table, int canInc) {
	char c;
	int i, code[256];
	TypeNumber sizeBufferSeq, n;
	TypeAlignment align;
	TypePosition *bufferSize, max, position;

	align.cardinal = 0;
	for(i=0; i<256; i++)
		code[i] = -1;
	if(table != NULL) {		
		for(i=0; table[i] != '\0'; i++) {
			code[table[i]] = i;
			align.cardinal++;
		}
	}	
	align.number = 0; align.empty = EMPTY;
	sizeBufferSeq = INC_BUFFER_SEQ;
	align.sequence = (TypeSymbol **) monmalloc(sizeBufferSeq*sizeof(TypeSymbol*));
	bufferSize = (TypePosition *) monmalloc(sizeBufferSeq*sizeof(TypePosition));
	align.name = (char **) monmalloc(sizeBufferSeq*sizeof(char *));
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
						code[(unsigned char) toupper(c)] = align.cardinal++;
					} else {
						char *msg;
						msg = (char*) monmalloc(100*sizeof(char));
						sprintf(msg, "Model and sequences aren't compatible (symbol %c not in the model)", toupper(c));
						exitProg(ErrorInit, msg);
					}
				}
				symbol = code[(unsigned char) toupper(c)];
				if(position>=sizeBufferSymbol) {
					sizeBufferSymbol += INC_BUFFER_SYMBOL;
					bufferSymbol = (TypeSymbol*) monrealloc(bufferSymbol, sizeBufferSymbol*sizeof(TypeSymbol));
				}
				bufferSymbol[position++] = symbol;
			} else {
				if(isempty(c)) {
					if(position>=sizeBufferSymbol) {
						sizeBufferSymbol += INC_BUFFER_SYMBOL;
						bufferSymbol = (TypeSymbol*) monrealloc(bufferSymbol, sizeBufferSymbol*sizeof(TypeSymbol));
					}
					bufferSymbol[position++] = align.empty;
				}
			}
			c = getc(f);
		}
		if(position>=SIZE_MIN) {
			if(align.number>=sizeBufferSeq) {
				sizeBufferSeq += INC_BUFFER_SEQ;
				align.sequence = (TypeSymbol **) monrealloc(align.sequence, sizeBufferSeq*sizeof(TypeSymbol*));
				bufferSize = (TypePosition *) monrealloc(bufferSize, sizeBufferSeq*sizeof(TypePosition));
				align.name = (char **) monrealloc(align.name, sizeBufferSeq*sizeof(char *));
			}
			align.sequence[align.number] = (TypeSymbol*) monrealloc(bufferSymbol, position*sizeof(TypeSymbol));
			align.name[align.number] = (char *) monrealloc((void *) bufferName, sizeTmp*sizeof(char));
			bufferSize[align.number] = position;
			align.number++;
		} else {
			monfree((void*)bufferSymbol);
		}
	}
	align.sequence = (TypeSymbol **) monrealloc(align.sequence, align.number*sizeof(TypeSymbol*));
	bufferSize = (TypePosition *) monrealloc(bufferSize, align.number*sizeof(TypePosition));
	max = 0;
	for(n=0; n<align.number; n++) {
		if(bufferSize[n]>max) {
			max = bufferSize[n];
		}
	}
	align.size = max;
	for(n=0; n<align.number; n++) {
		align.sequence[n] = (TypeSymbol*) monrealloc(align.sequence[n], align.size*sizeof(TypeSymbol));
		for(position=bufferSize[n]; position<align.size; position++) {
			align.sequence[n][position] = align.empty;
		}
	}
	monfree((void*)bufferSize);
	align.name = (char **) monrealloc(align.name, align.number*sizeof(char *));
	if(table == NULL || canInc) {
		if(table == NULL) {
			align.table = (char *) monmalloc((align.cardinal+1)*sizeof(char *));
		} else {
			align.table = (char *) monrealloc(table, (align.cardinal+1)*sizeof(char *));
		}
		for(i=0; i<256; i++)
			if(code[i] != -1)
				align.table[code[i]] = (char) i;
		align.table[align.cardinal] = '\0';
	} else {
		align.table = table;
	}
	align.comment.number = 0;
	align.comment.text = NULL;
	return align;
}

/*Read Alignement MSF*/
TypeAlignment readAlignementMsf(FILE *f, char *table, int canInc) {
	char c;
	int i, code[256], ntoken, sizeline;
	TypeNumber sizeBufferSeq, n;
	TypeAlignment align;
	TypePosition *bufferSize, max, position, sizeBufferSymbol;
	char *buffer, **token;

	buffer = (char*) monmalloc(INC_BUFFER_SYMBOL*sizeof(char));
	token = (char**) monmalloc(INC_BUFFER_SYMBOL*sizeof(char*));
	align.cardinal = 0;
	for(i=0; i<256; i++)
		code[i] = -1;
	if(table != NULL) {		
		for(i=0; table[i] != '\0'; i++) {
			code[table[i]] = i;
			align.cardinal++;
		}
	}	
	align.number = 0; align.empty = EMPTY;
	sizeBufferSeq = INC_BUFFER_SEQ;
	sizeBufferSymbol = INC_BUFFER_SYMBOL;
	align.sequence = (TypeSymbol **) monmalloc(sizeBufferSeq*sizeof(TypeSymbol*));
	bufferSize = (TypePosition *) monmalloc(sizeBufferSeq*sizeof(TypePosition));
	align.name = (char **) monmalloc(sizeBufferSeq*sizeof(char *));
	sizeBufferSymbol = INC_BUFFER_SYMBOL;

	do {
		while((sizeline = readLine(f, buffer)) != 0 && 
			(ntoken=tokenize(buffer, token)) != 0 && 
			strcmp(token[0], "Name:") != 0 && 
			strcmp(token[0], "//") != 0)
			;
		if(sizeline>0 && ntoken > 0 && strcmp(token[0], "Name:") == 0) {
			if(align.number>=sizeBufferSeq) {
				sizeBufferSeq += INC_BUFFER_SEQ;
				align.sequence = (TypeSymbol **) monrealloc(align.sequence, sizeBufferSeq*sizeof(TypeSymbol*));
				bufferSize = (TypePosition *) monrealloc(bufferSize, sizeBufferSeq*sizeof(TypePosition));
				align.name = (char **) monrealloc(align.name, sizeBufferSeq*sizeof(char *));
			}
			if(ntoken<2)
				exitProg(ErrorReading, "This is not a valid MSF file");
			align.name[align.number] = (char *) monmalloc((strlen(token[1])+1)*sizeof(char));
			strcpy(align.name[align.number], token[1]);
			align.sequence[align.number] = (TypeSymbol*) monmalloc(sizeBufferSymbol*sizeof(TypeSymbol));
			bufferSize[align.number] = 0;
			align.number++;
		}
	} while(sizeline != 0 && ntoken>0 && strcmp(token[0], "//") != 0);

	align.sequence = (TypeSymbol **) monrealloc(align.sequence, align.number*sizeof(TypeSymbol*));
	bufferSize = (TypePosition *) monrealloc(bufferSize, align.number*sizeof(TypePosition));
	align.name = (char **) monrealloc(align.name, align.number*sizeof(char *));
	
	do {
		TypeNumber current;
		while((sizeline = readLine(f, buffer)) != 0 && 
			(ntoken=tokenize(buffer, token)) != 0 && 
			(current = find(token[0], align.name, align.number))>=align.number)
			;
		if(sizeline > 0 && ntoken > 0 && current<align.number) {
		for(i=1; i<ntoken; i++) {
			int j;
			for(j=0; token[i][j] != '\0'; j++) {
				TypeSymbol symbol;
				if(isalpha(token[i][j])) {
					if(code[(unsigned char) toupper(token[i][j])] == -1) {
						if(table == NULL || canInc) {
							code[(unsigned char) toupper(token[i][j])] = align.cardinal++;
						} else{
							char *msg = (char*) monmalloc(100*sizeof(char));
							sprintf(msg, "Model and sequences aren't compatible (symbol %c not in the model)", toupper(token[i][j]));
							exitProg(ErrorInit, msg);
						}
					}
					symbol = code[(unsigned char) toupper(token[i][j])];
					if(bufferSize[current]>=sizeBufferSymbol) {
						TypeNumber n;
						sizeBufferSymbol += INC_BUFFER_SYMBOL;
						for(n=0; n<align.number; n++) {
							align.sequence[n] = (TypeSymbol*) monrealloc(align.sequence[n], sizeBufferSymbol*sizeof(TypeSymbol));
						}
					}
					align.sequence[current][bufferSize[current]++] = symbol;
				} else {
					if(isempty(token[i][j])) {
						if(bufferSize[current]>=sizeBufferSymbol) {
							TypeNumber n;
							sizeBufferSymbol += INC_BUFFER_SYMBOL;
							for(n=0; n<align.number; n++) {
								align.sequence[n] = (TypeSymbol*) monrealloc(align.sequence[n], sizeBufferSymbol*sizeof(TypeSymbol));
							}
						}
						align.sequence[current][bufferSize[current]++] = align.empty;
					}
				}
			}
		}
		}	
	} while(sizeline != 0);
	max = 0;
	for(n=0; n<align.number; n++) {
		if(bufferSize[n]>max) {
			max = bufferSize[n];
		}
	}
	align.size = max;
	for(n=0; n<align.number; n++) {
		align.sequence[n] = (TypeSymbol*) monrealloc(align.sequence[n], align.size*sizeof(TypeSymbol));
		for(position=bufferSize[n]; position<align.size; position++) {
			align.sequence[n][position] = align.empty;
		}
	}
	monfree((void*)bufferSize);
	if(table == NULL || canInc) {
		if(table == NULL) {
			align.table = (char *) monmalloc((align.cardinal+1)*sizeof(char *));
		} else {
			align.table = (char *) monrealloc(table, (align.cardinal+1)*sizeof(char *));
		}
		for(i=0; i<256; i++)
			if(code[i] != -1)
				align.table[code[i]] = (char) i;
		align.table[align.cardinal] = '\0';
	} else {
		align.table = table;
	}
	align.comment.number = 0;
	align.comment.text = NULL;
	monfree((void*) buffer);
	monfree((void*) token);
	return align;
}

/*Read Alignement Clustal*/
TypeAlignment readAlignementClustal(FILE *f, char *table, int canInc) {
	char c;
	int i, code[256], ntoken, sizeline;
	TypeNumber sizeBufferSeq, n;
	TypeAlignment align;
	TypePosition *bufferSize, max, position, sizeBufferSymbol;
	char *buffer, **token;

	buffer = (char*) monmalloc(INC_BUFFER_SYMBOL*sizeof(char));
	token = (char**) monmalloc(INC_BUFFER_SYMBOL*sizeof(char*));
	align.cardinal = 0;
	for(i=0; i<256; i++)
		code[i] = -1;
	if(table != NULL) {		
		for(i=0; table[i] != '\0'; i++) {
			code[table[i]] = i;
			align.cardinal++;
		}
	}	
	align.number = 0; align.empty = EMPTY;
	sizeBufferSeq = INC_BUFFER_SEQ;
	sizeBufferSymbol = INC_BUFFER_SYMBOL;
	align.sequence = (TypeSymbol **) monmalloc(sizeBufferSeq*sizeof(TypeSymbol*));
	bufferSize = (TypePosition *) monmalloc(sizeBufferSeq*sizeof(TypePosition));
	align.name = (char **) monmalloc(sizeBufferSeq*sizeof(char *));
	sizeBufferSymbol = INC_BUFFER_SYMBOL;

	while((sizeline = readLine(f, buffer)) != 0 && 
		(ntoken=tokenize(buffer, token)) != 0 && 
		strcmp(token[0], "CLUSTAL") != 0)
		;
	do {
		TypeNumber current;
		while((sizeline = readLine(f, buffer)) != 0 && IsSeparator(buffer[0]));
#ifdef DEBUG
printf("line %d %s\n",sizeline, buffer);
#endif
		if(sizeline != 0) {
			ntoken = tokenize(buffer, token);
			current = find(token[0], align.name, align.number);
#ifdef DEBUG
printf("token %s %d %d\n", token[0], current, align.number);
#endif
			if(current>=align.number) {
				if(align.number>=sizeBufferSeq) {
					sizeBufferSeq += INC_BUFFER_SEQ;
					align.sequence = (TypeSymbol **) monrealloc(align.sequence, sizeBufferSeq*sizeof(TypeSymbol*));
					bufferSize = (TypePosition *) monrealloc(bufferSize, sizeBufferSeq*sizeof(TypePosition));
					align.name = (char **) monrealloc(align.name, sizeBufferSeq*sizeof(char *));
				}
				align.name[align.number] = (char *) monmalloc((strlen(token[0])+1)*sizeof(char));
				strcpy(align.name[align.number], token[0]);
				align.sequence[align.number] = (TypeSymbol*) monmalloc(sizeBufferSymbol*sizeof(TypeSymbol));
				bufferSize[align.number] = 0;
				align.number++;
			}
			for(i=1; i<ntoken; i++) {
				int j;
				if(token[i][0] != '[') {
					for(j=0; token[i][j] != '\0'; j++) {
						TypeSymbol symbol;
						if(isalpha(token[i][j])) {
							if(code[(unsigned char) toupper(token[i][j])] == -1) {
							if(table == NULL || canInc) {
								code[(unsigned char) toupper(token[i][j])] = align.cardinal++;
							} else {
								char *msg = (char*) monmalloc(100*sizeof(char));
								sprintf(msg, "Model and sequences aren't compatible (symbol %c not in the model)", toupper(token[i][j]));
								exitProg(ErrorInit, msg);
							}
							}
							symbol = code[(unsigned char) toupper(token[i][j])];
							if(bufferSize[current]>=sizeBufferSymbol) {
								TypeNumber n;
								sizeBufferSymbol += INC_BUFFER_SYMBOL;
								for(n=0; n<align.number; n++) {
									align.sequence[n] = (TypeSymbol*) monrealloc(align.sequence[n], sizeBufferSymbol*sizeof(TypeSymbol));
								}
							}
							align.sequence[current][bufferSize[current]++] = symbol;
						} else {
							if(isempty(token[i][j])) {
								if(bufferSize[current]>=sizeBufferSymbol) {
									TypeNumber n;
									sizeBufferSymbol += INC_BUFFER_SYMBOL;
									for(n=0; n<align.number; n++) {
										align.sequence[n] = (TypeSymbol*) monrealloc(align.sequence[n], sizeBufferSymbol*sizeof(TypeSymbol));
									}
								}
								align.sequence[current][bufferSize[current]++] = align.empty;
							}
						}
					}
				}
			}
		}	
	} while(sizeline != 0);
	monfree((void*) buffer);
	monfree((void*) token);
#ifdef DEBUG
printf("ENd %d \n", align.number);
#endif

	align.sequence = (TypeSymbol **) monrealloc(align.sequence, align.number*sizeof(TypeSymbol*));
	bufferSize = (TypePosition *) monrealloc(bufferSize, align.number*sizeof(TypePosition));
	align.name = (char **) monrealloc(align.name, align.number*sizeof(char *));
	max = 0;
	for(n=0; n<align.number; n++) {
		if(bufferSize[n]>max) {
			max = bufferSize[n];
		}
	}
	align.size = max;
	for(n=0; n<align.number; n++) {
		align.sequence[n] = (TypeSymbol*) monrealloc(align.sequence[n], align.size*sizeof(TypeSymbol));
		for(position=bufferSize[n]; position<align.size; position++) {
			align.sequence[n][position] = align.empty;
		}
	}
	monfree((void*)bufferSize);
	if(table == NULL || canInc) {
		if(table == NULL) {
			align.table = (char *) monmalloc((align.cardinal+1)*sizeof(char *));
		} else {
			align.table = (char *) monrealloc(table, (align.cardinal+1)*sizeof(char *));
		}
		for(i=0; i<256; i++)
			if(code[i] != -1)
				align.table[code[i]] = (char) i;
		align.table[align.cardinal] = '\0';
	} else {
		align.table = table;
	}
	align.comment.number = 0;
	align.comment.text = NULL;
	return align;
}

void printAlignmentFasta(FILE *f, TypeAlignment al, int sizeLine) {
	TypePosition position;
	TypeNumber n;
	for(n=0; n<al.number; n++) {
		if(al.name != NULL)
			fprintf(f, ">%s\n", al.name[n]);
		else
			fprintf(f, ">Sequence %d\n", n+1);
		for(position=0; position<al.size; position += sizeLine) {
			TypePosition j, maxLine;
			maxLine = position+sizeLine;
			maxLine = maxLine>al.size?al.size:maxLine;
			if(al.table != NULL) {
				for(j=position; j<maxLine; j++)
					if(al.sequence[n][j] != al.empty)
						fprintf(f, "%c", al.table[al.sequence[n][j]]);
					else
						fprintf(f, "%c", EMPTY_CAR);
			} else {
				for(j=position; j<maxLine; j++)
					if(al.sequence[n][j] != al.empty)
						fprintf(f, "%d", al.sequence[n][j]);
					else
						fprintf(f, "%c", EMPTY_CAR);
			}
			fprintf(f, "\n");
		}
		fprintf(f, "\n");
	}
}

void printAlignmentMsf(FILE *f, TypeAlignment al, int sizeLine) {
	TypePosition position;
	TypeNumber n;
	int i, margin=2;
	
	fprintf(f, "PileUp\n");
	fprintf(f, "MSF: %ld Type:? ..\n\n", al.size);
	for(n=0; n<al.number; n++) {
		fixSpace(al.name[n]);
		if(strlen(al.name[n])>margin)
			margin = strlen(al.name[n]);
		fprintf(f, "Name: %s\tLen: %ld\tCheck: %ld\tWeight: %f\n", al.name[n], al.size, al.size, 1.0);
	}
	margin += 2;
	fprintf(f, "\n//\n\n");
	for(position=0; position<al.size; position += sizeLine) {
		TypePosition j, maxLine;
		TypeNumber n;
		int i;
		maxLine = position+sizeLine;
		maxLine = maxLine>al.size?al.size:maxLine;
		for(n=0; n<al.number; n++) {
			fprintf(f, "%s", al.name[n]);
			for(i=strlen(al.name[n]); i<margin; i++)
				fprintf(f, " ");
			if(al.table != NULL) {
				for(j=position; j<maxLine; j+=INC) {
					TypePosition maxGroup, k;
					maxGroup = j+INC;
					maxGroup = maxGroup>maxLine?maxLine:maxGroup;
					for(k=j; k<maxGroup; k++) {
						if(al.sequence[n][k] != al.empty) {
							fprintf(f, "%c", al.table[al.sequence[n][k]]);
						} else {
							fprintf(f, "%c", EMPTY_CAR);
						}
					}
					fprintf(f, " ");
				}
			} else {
				for(j=position; j<maxLine; j+=INC) {
					TypePosition maxGroup, k;
					maxGroup = j+INC;
					maxGroup = maxGroup>al.size?al.size:maxGroup;
					for(k=j; k<maxGroup; k++) {
						if(al.sequence[n][k] != al.empty) {
							fprintf(f, "%d", al.sequence[n][k]);
						} else {
							fprintf(f, "%c", EMPTY_CAR);
						}
					}
					fprintf(f, " ");
				}
			}
			fprintf(f, "\n");
		}
		fprintf(f, "\n");
	}
}

void printAlignmentSrs(FILE *f, TypeAlignment al, int sizeLine) {
	TypePosition position, *pos, nident=0, ngap=0;
	TypeNumber n;
	int i, margin = 0;

	pos = (TypePosition*) monmalloc(al.number*sizeof(TypePosition));
	for(n=0; n<al.number; n++)
		pos[n] = 1;
	if(al.number == 2)
		printHeadPair(f, al);
	else	
		printHeadMulti(f, al);
	for(n=0; n<al.number; n++) {
		fixSpace(al.name[n]);
		if(strlen(al.name[n])>margin)
			margin = strlen(al.name[n]);
	}
	margin += 2;
	for(position=0; position<al.size; position += sizeLine) {
		TypePosition j, maxLine;
		TypeNumber n;
		int i;
		maxLine = position+sizeLine;
		maxLine = maxLine>al.size?al.size:maxLine;
		for(n=0; n<al.number; n++) {
			fprintf(f, "\n%s:", al.name[n]);
			for(i=strlen(al.name[n]); i<margin-(1+(int)floor(log(pos[n])/log(10))); i++)
				fprintf(f, " ");
			fprintf(f, "%ld  ", pos[n]);
			if(al.table != NULL) {
				for(j=position; j<maxLine; j++) {
					if(al.sequence[n][j] != al.empty) {
						fprintf(f, "%c", al.table[al.sequence[n][j]]);
						pos[n]++;
					} else
						fprintf(f, "%c", EMPTY_CAR);
				}
			} else {
				for(j=position; j<maxLine; j++)
					if(al.sequence[n][j] != al.empty) {
						fprintf(f, "%d", al.sequence[n][j]);
						pos[n]++;
					} else
						fprintf(f, "%c", EMPTY_CAR);
			}
			for(i=0; i<RIGHT-(1+(int)floor(log(pos[n]-1)/log(10))); i++)
				fprintf(f, " ");
			fprintf(f, "%ld\n", pos[n]-1);
		}
		fprintf(f, "\n\n");
	}
	monfree((void*) pos);
}

void printAlignmentMarkX(FILE *f, TypeAlignment al, int sizeLine) {
	TypePosition position, pos1 = 1, pos2 = 1, nident=0, ngap=0;
	int i, margin = 0, n;

	printHeadPair(f, al);
	for(n=0; n<2; n++) {
		fixSpace(al.name[n]);
		if(strlen(al.name[n])>margin)
			margin = strlen(al.name[n]);
	}
	margin += 2;
	for(position=0; position<al.size; position += sizeLine) {
		TypePosition j, maxLine;
		TypeNumber n;
		int i;
		maxLine = position+sizeLine;
		maxLine = maxLine>al.size?al.size:maxLine;
		fprintf(f, "\n\n");
		for(i=0; i<margin; i++)
				fprintf(f, " ");
		for(j=position; j<maxLine;) {
			if(al.sequence[0][j]!=al.empty && (pos1)%INC == 0) {
				TypePosition inc, k;
				fprintf(f, "%ld", pos1);
				inc = 1+(TypePosition) floor(log(pos1)/log(10));
				inc = (inc<(maxLine-j))?inc:maxLine-j;
				for(k=0; k<inc; k++)
					if(al.sequence[0][j+k] != al.empty)
						pos1++;
				j += inc;
			} else {
				fprintf(f, " ");
				if(al.sequence[0][j] != al.empty)
					pos1++;
 				j++;
 			}
		}
		fprintf(f, "\n%s:", al.name[0]);
		for(i=strlen(al.name[0])+1; i<margin; i++)
			fprintf(f, " ");
		if(al.table != NULL) {
			for(j=position; j<maxLine; j++)
				if(al.sequence[0][j] != al.empty)
					fprintf(f, "%c", al.table[al.sequence[0][j]]);
				else
					fprintf(f, "%c", EMPTY_CAR);
		} else {
			for(j=position; j<maxLine; j++)
				if(al.sequence[0][j] != al.empty)
					fprintf(f, "%d", al.sequence[1][j]);
				else
					fprintf(f, "%c", EMPTY_CAR);
		}
		fprintf(f, "\n");
		for(i=0; i<margin; i++)
			fprintf(f, " ");
		for(j=position; j<maxLine; j++) {
			if(al.sequence[0][j] == al.empty || al.sequence[1][j] == al.empty) {
				fprintf(f, " ");
			} else {
				if(al.sequence[0][j] == al.sequence[1][j])
					fprintf(f, "|");
				else
					fprintf(f, "x");
			}
		}
		fprintf(f, "\n%s:", al.name[1]);
		for(i=strlen(al.name[1])+1; i<margin; i++)
			fprintf(f, " ");
		if(al.table != NULL) {
			for(j=position; j<maxLine; j++)
				if(al.sequence[1][j] != al.empty)
					fprintf(f, "%c", al.table[al.sequence[1][j]]);
				else
					fprintf(f, "%c", EMPTY_CAR);
		} else {
			for(j=position; j<maxLine; j++)
				if(al.sequence[1][j] != al.empty)
					fprintf(f, "%d", al.sequence[1][j]);
				else
					fprintf(f, "%c", EMPTY_CAR);
		}
		fprintf(f, "\n");
		for(i=0; i<margin; i++)
			fprintf(f, " ");
		for(j=position; j<maxLine;) {
			if(al.sequence[1][j]!=al.empty && (pos2)%INC == 0) {
				TypePosition inc, k;
				fprintf(f, "%ld", pos2);
				inc = 1+(TypePosition) floor(log(pos2)/log(10));
				inc = (inc<(maxLine-j))?inc:maxLine-j;
				for(k=0; k<inc; k++)
					if(al.sequence[1][j+k]!=al.empty)
						pos2++;
					j += inc;
			} else {
				fprintf(f, " ");
				if(al.sequence[1][j]!=al.empty)
					pos2++;
 				j++;
 			}
		}
	}
	fprintf(f, "\n");
}

void printHeadPair(FILE *f, TypeAlignment al) {
	time_t curTime;
	TypePosition position, pos1 = 1, pos2 = 1, nident=0, ngap=0, *length;
	int i, n;
	
	length = (TypePosition *) monmalloc(2*sizeof(TypePosition));
	for(n=0; n<2; n++)
		length[n] = 0;
	for(position=0; position<al.size; position ++)
		for(n=0; n<2; n++)
			if(al.sequence[n][position] != al.empty)
				length[n]++;
	for(position=0; position<al.size; position ++) {
		if(al.sequence[0][position] == al.sequence[1][position])
			nident++;
		if(al.sequence[0][position] == al.empty || al.sequence[1][position] == al.empty)
			ngap++;
	}

	fprintf(f, "######################################################################\n");
	fprintf(f, "# Best Contextual Ancestor Align\n");
	fprintf(f, "# \n");
	fprintf(f, "######################################################################\n");
	curTime = time(&curTime);
	fprintf(f, "# %s", asctime(localtime(&curTime)));
	fprintf(f, "# Aligned Sequences: %ld\n", 2);
	if(al.name != NULL) {
		TypeNumber n;
		for(n=0; n<2; n++)
			fprintf(f, "# %d: %s (length %ld)\n", n+1, al.name[n], length[n]);
	} else {
		TypeNumber n;
		for(n=0; n<2; n++)
			fprintf(f, "# %d: sequence %d (length %ld)\n", n+1, n+1, length[n]);
	}
	monfree((void*) length);
	fprintf(f, "# Length: %ld\n", al.size);
	fprintf(f, "# Identity: %ld/%ld (%.2f\%)\n", nident, al.size, ((float) 100*nident)/((float)al.size));
	fprintf(f, "# Gaps: %ld/%ld (%.2f\%)\n", ngap, al.size, ((float) 100*ngap)/((float)al.size));
	for(i=0; i<al.comment.number; i++)
		fprintf(f, "# %s\n", al.comment.text[i]);
	fprintf(f, "######################################################################\n");
}

void printHeadMulti(FILE *f, TypeAlignment al) {
	time_t curTime;
	TypePosition position, *pos, nident=0, ngap=0, *length;
	TypeNumber n;
	int i;
	
	length = (TypePosition *) monmalloc(al.number*sizeof(TypePosition));
	for(n=0; n<al.number; n++)
		length[n] = 0;
	for(position=0; position<al.size; position ++) {
		for(n=0; n<al.number; n++)
			if(al.sequence[n][position] == al.empty)
				ngap++;
			else
				length[n]++;
	}
	fprintf(f, "######################################################################\n");
	fprintf(f, "# Best Contextual Ancestor Align\n");
	fprintf(f, "# \n");
	fprintf(f, "######################################################################\n");
	curTime = time(&curTime);
	fprintf(f, "# %s", asctime(localtime(&curTime)));
	fprintf(f, "# Aligned Sequences: %ld\n", al.number);
	if(al.name != NULL) {
		TypeNumber n;
		for(n=0; n<al.number; n++)
			fprintf(f, "# %d: %s (length %ld)\n", n+1, al.name[n], length[n]);
	} else {
		TypeNumber n;
		for(n=0; n<al.number; n++)
			fprintf(f, "# %d: sequence %d (length %ld)\n", n+1, n+1, length[n]);
	}
	monfree((void*) length);
	fprintf(f, "# Length: %ld\n", al.size);
	fprintf(f, "# Gaps: %ld/%ld (%.2f\%)\n", ngap, al.size, ((float) 100*ngap)/((float)al.size));
	for(i=0; i<al.comment.number; i++)
		fprintf(f, "# %s\n", al.comment.text[i]);
	fprintf(f, "######################################################################\n");
}


