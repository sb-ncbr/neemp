# NEEMP - Makefile
# by Tomas Racek (tom@krab1k.net)
# 2013, 2014

CC=icc
CFLAGS=-Wall -Wcheck -Wremarks -std=c99 -g -Ofast -ipo -Xhost -wd981 -fopenmp
# -wd981 disables warning about operands evaluated in an unspecified order
EXTRA_DEFINE=-DUSE_MKL

sources=$(wildcard *.c)
headers=$(wildcard *.h)
objects=$(sources:.c=.o)
libraries=-mkl
binaries=neemp
manpage=neemp.1

all: $(sources) $(headers) neemp

neemp-gnu: EXTRA_DEFINE=
neemp-gnu: CC=gcc
neemp-gnu: CFLAGS=-Wall -Wextra -std=c99 -pedantic -O3 -march=native -g -fopenmp
neemp-gnu: libraries=-lm -fopenmp
neemp-gnu: $(objects) neemp

neemp: $(objects)
	$(CC) $(objects) $(libraries) -o neemp

.c.o: $(headers) $(sources)
	$(CC) $(CFLAGS) $(EXTRA_DEFINE) -c $<
man: neemp
	help2man -N ./neemp -o $(manpage)
clean:
	rm -f $(objects) $(manpage) $(binaries)
