CC=gcc
DEBUGFLAGS=-DDEBUG
CFLAGS=-g -Wall -O2
OMPFLAGS =-fopenmp

all: gepp omp

gepp: gepp.c
	$(CC) $(CFLAGS) -o gepp gepp.c $(OMPFLAGS)

omp: gepp.c
	$(CC) $(CFLAGS) -o gepp_omp gepp_omp.c $(OMPFLAGS)

example: example.c
	$(CC) $(CFLAGS) -o example example.c $(OMPFLAGS)

debug: gepp.c
	$(CC) $(CFLAGS) $(DEBUGFLAGS) -o gepp_d gepp.c $(OMPFLAGS)

clean:
	rm -f gepp example gepp_d

