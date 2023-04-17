
CC = gcc

LDFLAGS = -L/usr/local/lib -lm 


all: local

local: Main.c
	$(CC)  Utils.c Sequences.c Distance.c Decode.c DistSeq.c Main.c $(LDFLAGS) -o $@	

