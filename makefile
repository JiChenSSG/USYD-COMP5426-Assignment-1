CC=gcc
DEBUGFLAGS=-DDEBUG
CFLAGS=-g -Wall -O3
OMPFLAGS =-fopenmp

all: gepp debug

gepp: gepp.c
	$(CC) $(CFLAGS) -o gepp gepp.c $(OMPFLAGS)

example: example.c
	$(CC) $(CFLAGS) -o example example.c $(OMPFLAGS)

debug: gepp.c
	$(CC) $(CFLAGS) $(DEBUGFLAGS) -o gepp_d gepp.c $(OMPFLAGS)

clean:
	rm -f gepp example gepp_d

