CC=gcc
CXX=g++
CFLAGS=-g -Wall -O2 -I..
CXXFLAGS=$(CFLAGS)
PROGS= minimizer

all:$(PROGS)

clean:
		rm -fr $(PROGS) a.out

minimizer:minimizer.c kseq.h
		$(CC) $(CFLAGS) -o $@ minimizer.c -lz

