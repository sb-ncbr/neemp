#  please change the variables FC and FFLAGS to match your fortran compiler

CC = icc
CFLAGS=-Wall -Wcheck -Wremarks -std=c99 -g -Ofast -ipo -Xhost -wd981 -qopenmp
ICC_DISABLE_WARNS=-wd11021
sources=$(wildcard *.c)
headers=$(wildcard *.h)
objects=$(sources:.c=.o)
libraries=-mkl -lm -lz

all: $(objects)

gnu: EXTRA_DEFINE=
gnu: CC=gcc
gnu: ICC_DISABLE_WARNS=
gnu: CFLAGS=-Wall -Wextra -std=c99 -pedantic -O3 -march=native -g -gdwarf-3 -fopenmp
gnu: libraries=-lm -lz
gnu: $(objects)


%.o: %.c
	$(CC) $(CFLAGS) $(EXTRA_DEFINE) -c $<

clean:
	rm -f *.o
