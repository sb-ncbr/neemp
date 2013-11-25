# NEEMP - Makefile
# by Tomas Racek (tom@krab1k.net)
# 2013

CC=icc
CFLAGS=-Wall -Wcheck -Wremarks -std=c99 -g -Ofast -ipo -Xhost -wd981 -DDEBUG
# -wd981 disables warning about operands evaluated in an unspecified order
EXTRA_DEFINE=-DUSE_MKL

sources=$(wildcard *.c)
headers=$(wildcard *.h)
objects=$(sources:.c=.o)
binaries=neemp
manpage=neemp.1

all: $(sources) $(headers) neemp man

neemp: $(objects)
	$(CC) $(objects) -mkl -o neemp

neemp-nomkl: EXTRA_DEFINE=
neemp-nomkl: $(objects) man
	$(CC) $(objects) $(EXTRA_DEFINE) -o neemp
.c.o: $(headers) $(sources)
	$(CC) $(CFLAGS) $(EXTRA_DEFINE) -c $<
man: neemp
	help2man -N ./neemp -o $(manpage)
clean:
	rm -f $(objects) $(manpage) $(binaries)
