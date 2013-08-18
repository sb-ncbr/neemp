# NEEMP - Makefile
# by Tomas Racek (tom@krab1k.net)
# 2013

CC=icc
CFLAGS=-Wall -Wcheck -Wremarks -std=c99 -g -ipo -Xhost -wd981 #-Ofast

# -wd981 disables warning about operands evaluated in an unspecified order

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
