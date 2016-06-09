#  please change the variables FC and FFLAGS to match your fortran compiler

FC = ifort #gfortran
FFLAGS = -g

objects = bigden.o biglag.o newuoa.o newuob.o trsapp.o update.o
all: $(objects)
gnu: FC = gfortran
gnu: $(objects)

bigden.o: bigden.f
	$(FC) $(FFLAGS) -c bigden.f
biglag.o: biglag.f
	$(FC) $(FFLAGS) -c biglag.f
trsapp.o: trsapp.f
	$(FC) $(FFLAGS) -c trsapp.f
newuoa.o: newuoa.f
	$(FC) $(FFLAGS) -c newuoa.f
newuob.o: newuob.f
	$(FC) $(FFLAGS) -c newuob.f
update.o: update.f
	$(FC) $(FFLAGS) -c update.f

clean:
	rm -f *.o main
