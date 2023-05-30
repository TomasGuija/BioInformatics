
CC = gcc

LDFLAGS = -L/usr/local/lib -lm 


ALL = compare 

all: compare


compare: Main.c
	$(CC)  Utils.c  Distance.c Main.c $(LDFLAGS) -o $@
