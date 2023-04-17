#include <math.h>
#include "Distance.h"
#include "Utils.h"


#define INC_BUFFER_DIST 1000
#define SIZE_BUFFER_TypeFloat 100
#define INC_BUFFER_ANC 10
#define INC_BUFFER 50


TypeFloat computeCorrelation(TypeDistance dist0, TypeDistance dist1) {
	TypeNumber i, j;
	double s0 = 0.0, s02 = 0.0, s1 = 0.0, s12 = 0.0, sp = 0.0;
	if(dist0.number != dist1.number)
		return -2;
	for(i=1; i<dist0.number; i++)
		for(j=0; j<i; j++) {
			s0 += dist0.table[(i*(i-1))/2+j];
			s1 += dist1.table[(i*(i-1))/2+j];
			s02 += pow(dist0.table[(i*(i-1))/2+j], 2);
			s12 += pow(dist1.table[(i*(i-1))/2+j], 2);
			sp += dist0.table[(i*(i-1))/2+j]*dist1.table[(i*(i-1))/2+j];
		}
	s0 *= 2.0; s1 *= 2; s02 *= 2.0; s12 *= 2; sp*=2;
	return (pow((double) dist0.number, 2)*sp-s0*s1)/(sqrt(pow((double) dist0.number, 2)*s02-pow(s0,2))*sqrt(pow((double) dist0.number, 2)*s12-pow(s1,2)));
}

TypeFloat computeCorrelationLog(TypeDistance dist0, TypeDistance dist1) {
	TypeNumber i, j;
	double s0 = 0.0, s02 = 0.0, s1 = 0.0, s12 = 0.0, sp = 0.0;
	if(dist0.number != dist1.number)
		return -2;
	for(i=1; i<dist0.number; i++)
		for(j=0; j<i; j++) {
			s0 += log(dist0.table[(i*(i-1))/2+j]+1);
			s1 += log(dist1.table[(i*(i-1))/2+j]+1);
			s02 += pow(log(dist0.table[(i*(i-1))/2+j]+1), 2);
			s12 += pow(log(dist1.table[(i*(i-1))/2+j]+1), 2);
			sp += log(dist0.table[(i*(i-1))/2+j]+1)*log(dist1.table[(i*(i-1))/2+j]+1);
		}
	s0 *= 2.0; s1 *= 2; s02 *= 2.0; s12 *= 2; sp*=2;
	return (pow((double) dist0.number, 2)*sp-s0*s1)/(sqrt(pow((double) dist0.number, 2)*s02-pow(s0,2))*sqrt(pow((double) dist0.number, 2)*s12-pow(s1,2)));
}

TypeFloat computeMean(TypeDistance dist) {	
	TypeNumber i, j;
	double mean= 0.0;
	for(i=1; i<dist.number; i++)
		for(j=0; j<i; j++)
			mean += dist.table[(i*(i-1))/2+j];
	return (2.0*mean)/pow((double) dist.number, 2);
}

TypeFloat computeVariance(TypeDistance dist) {	
	TypeNumber i, j;
	double var = 0.0,  mean = computeMean(dist);
	for(i=1; i<dist.number; i++)
		for(j=0; j<i; j++)
			var += pow(2, dist.table[(i*(i-1))/2+j]-mean);
	return (2.0*var)/(pow((double) dist.number, 2)-1.0);
}

TypeFloat computeVarianceRow(TypeDistance dist) {	
	TypeNumber i, j;
	double row = 0.0;
	for(i=0; i<dist.number; i++) {
		double var = 0.0,  mean = 0;
		for(j=0; j<dist.number; j++)
			mean += dist.table[(i*(i-1))/2+j];
		mean /= (double) dist.number;
		for(j=0; j<dist.number; j++)
			var += pow(2, dist.table[(i*(i-1))/2+j]-mean);
		var /= (double) dist.number-1.0;
		row += var;
	}
	return row;
}


void printDistanceRaw(FILE *f, TypeDistance dist) {
	TypeNumber i, j;
	for(i=1; i<dist.number; i++) {
		fprintf(f, "'%s'", dist.name[i]);
		for(j=0; j<i; j++)
			fprintf(f, "\t%f",dist.table[(i*(i-1))/2+j]);
		fprintf(f, "\n");
	}
}

void printDistanceTable(FILE *f, TypeDistance dist) {
	TypeNumber i, j;
	fprintf(f, "dist");
	for(i=0; i<dist.number; i++)
		fprintf(f, "\t%s", dist.name[i]);
	fprintf(f, "\n");
	for(i=0; i<dist.number; i++) {
		fprintf(f, "%s", dist.name[i]);
		for(j=0; j<i; j++)
			fprintf(f, "\t%f",dist.table[(i*(i-1))/2+j]);
		fprintf(f, "\t0.00");
		for(j=i+1; j<dist.number; j++)
			fprintf(f, "\t%f",dist.table[(j*(j-1))/2+i]);
		fprintf(f, "\n");
	}
}

void printDistancePhylip(FILE *f, TypeDistance dist) {
	TypeNumber i, j;
	fprintf(f, "%d\n", dist.number);
	for(i=0; i<dist.number; i++) {
		int k, l;
		for(k=0; k<10 && dist.name[i][k] != '\0'; k++)
			fprintf(f, "%c", dist.name[i][k]);
		for(; k<10; k++)
			fprintf(f, " ");
		for(j=0; j<i; j++)
			fprintf(f, "  %f",dist.table[(i*(i-1))/2+j]);
		fprintf(f, "  0.0\n");
		for(j=i+1; j<dist.number; j++)
			fprintf(f, "  %f",dist.table[(j*(j-1))/2+i]);
		fprintf(f, "\n");
	}
}


void printDistanceNexus(FILE *f, TypeDistance dist) {
	TypeNumber i, j;
	printPreambleNexus(f, dist.name, dist.number);
	fprintf(f, "BEGIN distances;\n");
	fprintf(f, "    DIMENSIONS ntax=%d;\n", dist.number);
	fprintf(f, "    FORMAT\ntriangle=LOWER\ndiagonal\nlabels\nmissing=?\n;\n");
	fprintf(f, "MATRIX\n");
	for(i=0; i<dist.number; i++) {
		fprintf(f, "'%s'", dist.name[i]);
		for(j=0; j<i; j++)
			fprintf(f, "\t%f",dist.table[(i*(i-1))/2+j]);
		fprintf(f, "\t0.0\n");
	}
	fprintf(f, ";\n");
	fprintf(f, "END;\n");
}


/*
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
*/

TypeDistance readDistancesNexus(FILE *f) {
	char c, *buffer;
	int sizeBuffer;
	TypeNumber sizeName;
	long sizeBufferDist, indTable = 0;
	TypeDistance dist;
	
	dist.number = 0;
	sizeName = INC_BUFFER;
	dist.name = (char **) monmalloc(sizeName*sizeof(char*));
	sizeBufferDist = INC_BUFFER_DIST;
	dist.table = (TypeFloat *) monmalloc(sizeBufferDist*sizeof(TypeFloat));
	sizeBuffer = INC_BUFFER;
	buffer = (char*) monmalloc(sizeBuffer*sizeof(char));
	while((c=getc(f))!=EOF && IsLineSeparator(c))
		;
	if(c==EOF)
		exitProg(ErrorReading, "Empty file!");
	do {
		int tot = 0;
		while(c!=EOF && IsLineSeparator(c))
			c=getc(f);
		for(tot=0; tot<7 && c!=EOF && !IsLineSeparator(c);) {
			buffer[tot++] = c;
			c=getc(f);
		}
		while(c!=EOF && !IsLineSeparator(c))
			c=getc(f);
	} while(strncmp(buffer, "MATRIX", 6));
	while(c!=EOF && !IsLineSeparator(c))
		c=getc(f);
	while(c!=EOF && IsLineSeparator(c))
		c=getc(f);
	do {
		TypeNumber i;
		do {
			int tot = 0;
			while(c!=EOF && isspace(c))
				c=getc(f);
			if(c == '\'') {
				int ind, sizeTmp = INC_BUFFER;
				c=getc(f);
				if(dist.number>=sizeName) {
					sizeName += INC_BUFFER;
					dist.name = (char **) monrealloc((void *) dist.name, sizeName*sizeof(char*));
				}
				dist.name[dist.number] = (char *) monmalloc(sizeTmp*sizeof(char));
				for(ind=0; c!=EOF && !IsLineSeparator(c) && c != '\''; ind++) {
					if(ind >= sizeTmp) {
						sizeTmp += INC_BUFFER;
						dist.name[dist.number] = (char *) monrealloc((void *) dist.name[dist.number], sizeTmp*sizeof(char));
					}
					dist.name[dist.number][ind] = c;
					c=getc(f);
				}
				if(c != '\'')
					exitProg(ErrorReading, "Problem name");
				dist.name[dist.number] = (char *) monrealloc((void *) dist.name[dist.number], (ind+1)*sizeof(char));
				dist.name[dist.number][ind] = '\0';
			} else {
				while(c!=EOF && !isspace(c))
					c=getc(f);
			}
		} while(c != '\'');
		c=getc(f);
		for(i=0; i<dist.number; i++) {
			int tot = 0;
			while(c != EOF && isspace(c)) {
				c = getc(f);
			}
			for(tot=0; tot<sizeBuffer && c != EOF && !isspace(c);) {
				if(c != ',')
					buffer[tot++] = c;
				else
					buffer[tot++] = '.';
				c = getc(f);
			}
			if(indTable >= sizeBufferDist) {
				sizeBufferDist += INC_BUFFER_DIST;
				dist.table = (TypeFloat *) monrealloc((void *) dist.table, sizeBufferDist*sizeof(TypeFloat));
			}
			if(sscanf(buffer, "%f", &(dist.table[indTable++])) != 1)
				exitProg(ErrorReading, "Problem number");
		}
		dist.number++;
		while(c!=EOF && !IsLineSeparator(c))
			c=getc(f);
		while(c!=EOF && IsLineSeparator(c))
			c=getc(f);
	} while(c != EOF && c != ';' );
	free((void*) buffer);
	dist.name = (char**) monrealloc((void *) dist.name, dist.number*sizeof(char*));
	dist.table = (TypeFloat*) monrealloc((void *) dist.table, ((dist.number*(dist.number-1))/2)*sizeof(TypeFloat));
	return dist;
}


TypeDistance readDistances(FILE *f) {
	char c, *bufferTypeFloat;
	TypeNumber sizeTmp = 0, sizeBufferDist;
	TypeFloat *bufferDist;
	TypeDistance dist;

	sizeBufferDist = INC_BUFFER_DIST;
	bufferDist = (TypeFloat *) monmalloc(sizeBufferDist*sizeof(TypeFloat));
	bufferTypeFloat = (char *) monmalloc(SIZE_BUFFER_TypeFloat*sizeof(char));
	c = getc(f);
	if(c==EOF)
		exitProg(ErrorReading, "Empty file!");
	do {
		int indTmp = 0;
		TypeFloat distTmp;
		while(c != EOF && IsSeparator(c)) {
			c = getc(f);
		}
		while(c != EOF && indTmp<SIZE_BUFFER_TypeFloat  && !IsSeparator(c)) {
			if(c != ',')
				bufferTypeFloat[indTmp++] = c;
			else
				bufferTypeFloat[indTmp++] = '.';
			c = getc(f);
		}
		if(indTmp<SIZE_BUFFER_TypeFloat)
			bufferTypeFloat[indTmp++] = ' ';
		if(indTmp>1 && sscanf(bufferTypeFloat, "%f", &distTmp) == 1) {
			if(sizeTmp >= sizeBufferDist) {
				sizeBufferDist += INC_BUFFER_DIST;
				bufferDist = (TypeFloat *) monrealloc((void *) bufferDist, sizeBufferDist*sizeof(TypeFloat));
			}
			bufferDist[sizeTmp++] = distTmp;
		}
	} while(c != EOF );
	free((void*) bufferTypeFloat);
	dist.number = (TypeNumber) floor((-1+sqrt(8*sizeTmp+1))/2);
	dist.table = (TypeFloat*) monrealloc((void *) bufferDist, ((dist.number*(dist.number-1))/2)*sizeof(TypeFloat));
	return dist;
}


TypeFloat distance(TypeDistance d, TypeNumber i, TypeNumber j) {
	if(i == j || i>=d.number || j>=d.number || i<0 || j<0)
		return 0;
	if(i<j)
		return d.table[(j*(j-1))/2+i];
	else
		return d.table[(i*(i-1))/2+j];
}

TypeFloat distanceClassesMin(TypeDistance d, TypeNumber *tab1, TypeNumber n1, TypeNumber *tab2, TypeNumber n2) {
	TypeNumber i, j;
	TypeFloat min;
	if(n1 == 0 || n2 == 0)
		return 0.0;
	min = distance(d, tab1[0], tab2[0]);
	for(i=0; i<n1; i++)
		for(j=0; j<n2; j++)
			if(min>distance(d, tab1[i], tab2[j]))
				min = distance(d, tab1[i], tab2[j]);
	return min;
}

TypeFloat distanceClassesMoy(TypeDistance d, TypeNumber *tab1, TypeNumber n1, TypeNumber *tab2, TypeNumber n2) {
	TypeNumber i, j;
	TypeFloat moy;
	if(n1 == 0 || n2 == 0)
		return 0.0;
	moy = 0.0;
	for(i=0; i<n1; i++)
		for(j=0; j<n2; j++)
				moy += distance(d, tab1[i], tab2[j]);
	return moy/(n1*n2);
}






