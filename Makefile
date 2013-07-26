# NEEMP - Makefile
# by Tomas Racek (tom@krab1k.net)
# 2013

CC=icc
CFLAGS=-Wall -Wcheck -Wremarks -std=c99 -g -ipo -Xhost #-Ofast

sources=$(wildcard *.c)
headers=$(wildcard *.h)
objects=$(sources:.c=.o)

all: $(sources) $(headers) neemp

neemp: $(objects)
	$(CC) $(objects) -o neemp
.c.o: $(headers) $(sources)
	$(CC) $(CFLAGS) -c $<
clean:
	rm -f $(objects) neemp
