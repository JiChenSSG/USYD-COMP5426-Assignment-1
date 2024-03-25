CC=gcc
DEBUGFLAGS=-g -DDEBUG
CFLAGS=-g -Wall -O3
OMPFLAGS =-fopenmp

all: gepp

gepp: gepp.c
	$(CC) $(CFLAGS) -o gepp gepp.c $(OMPFLAGS)

debug: gepp.c
	$(CC) $(CFLAGS) $(DEBUGFLAGS) -o gepp gepp.c $(OMPFLAGS)

clean:
	rm -f gepp

