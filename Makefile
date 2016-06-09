# NEEMP - Makefile
# by Tomas Racek (tom@krab1k.net)
# 2013, 2014

all:
	$(MAKE) -C externals/newuoa
	$(MAKE) -C externals/lhs
	$(MAKE) -C src all

neemp:
	$(MAKE) -C externals/newuoa
	$(MAKE) -C externals/lhs
	$(MAKE) -C src neemp

neemp-gnu:
	$(MAKE) -C externals/newuoa gnu
	$(MAKE) -C externals/lhs gnu
	$(MAKE) -C src neemp-gnu

man: neemp
	help2man -N ./neemp > neemp.1

clean:
	$(MAKE) -C externals/newuoa clean
	$(MAKE) -C externals/lhs clean
	$(MAKE) -C src clean
	rm -f ./neemp ./neemp.1
