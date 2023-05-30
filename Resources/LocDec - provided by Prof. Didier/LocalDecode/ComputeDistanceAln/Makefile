
CC = gcc

LDFLAGS = -L/usr/local/lib -lm 

all: align


align: Main.c
	$(CC)  Utils.c Alignement.c Distance.c DistAln.c Main.c $(LDFLAGS) -o $@

