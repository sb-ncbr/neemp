# NEEMP - Makefile
# by Tomas Racek (tom@krab1k.net)
# 2013, 2014

all:
	$(MAKE) -C src all

neemp:
	$(MAKE) -C src neemp

neemp-gnu:
	$(MAKE) -C src neemp-gnu

man: neemp
	help2man -N ./neemp > neemp.1

clean:
	$(MAKE) -C src clean
	rm -f ./neemp ./neemp.1
