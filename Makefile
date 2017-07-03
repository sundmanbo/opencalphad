OBJS=metlib3.o oclablas.o ocnum.o gtp3.o matsmin.o lmdif1lib.o smp2.o pmon6.o  
VER=oc4A
TQPATH=./TQ4lib/Cpp
TQPATHC=$(TQPATH)
TQPATHCPP=$(TQPATH)/Matthias
CFLAGS=-fbounds-check -finit-local-zero
#UNCOMMENT THE NEXT LINE IF YOU WANT TO COMPILE OPENCALPHAD IN PARALLEL
#CFLAGS+=-O2 -fopenmp

#==============================================================================#

#Available compilation flags: all, tq, examples, clean

.PHONY : all tq tqexamples clean

#Compiles OpenCalphad to use as standalone Thermodynamic Equilibrium Calculation
#software.

all:
	gfortran -o linkoc linkocdate.F90
	./linkoc
	make $(OBJS)
	ar sq lib$(VER)eq.a $(OBJS)
	gfortran -o $(VER) $(CFLAGS) pmain1.F90 lib$(VER)eq.a

#Compiles the OCASI interfaces of OpenCalphad, so third party software can
#interact with OpenCalphad. Interfaces are provided in C++ and Fortran, and
#additional interfaces are available in C and Python.

tq: 
	make $(OBJS) liboctq.o liboctqisoc.o
	ar sq liboctq-f90.a liboctq.o $(OBJS) liboctq.o
	ar sq liboctq-isoc.a liboctqisoc.o liboctq.o $(OBJS)

#Compiles the TQ interface and various examples. 

tqexamples:
	make tq
	make -C $(TQPATH)

#Removes all binary files that were created in the compiling step.

clean:
	rm -f *.a *.o *.mod $(VER) linkoc
	make -C $(TQPATH) clean

#==============================================================================#

metlib3.o:	utilities/metlib3.F90
	gfortran -c $(CFLAGS) utilities/metlib3.F90

oclablas.o:	numlib/oclablas.F90
	gfortran -c $(CFLAGS) numlib/oclablas.F90

ocnum.o:	numlib/ocnum.F90
	gfortran -c $(CFLAGS) numlib/ocnum.F90

lmdif1lib.o:      numlib/lmdif1lib.F90
	gfortran -c $(CFLAGS) numlib/lmdif1lib.F90  

gtp3.o:	models/gtp3.F90
	gfortran -c $(CFLAGS) models/gtp3.F90

matsmin.o:	minimizer/matsmin.F90
	gfortran -c $(CFLAGS) minimizer/matsmin.F90

smp2.o:		stepmapplot/smp2.F90
	gfortran -c $(CFLAGS) stepmapplot/smp2.F90

pmon6.o:	userif/pmon6.F90
	gfortran -c $(CFLAGS) userif/pmon6.F90
	
liboctq.o:	$(TQPATH)/liboctq.F90
	gfortran -c -g $(CFLAGS) $(TQPATH)/liboctq.F90

liboctqisoc.o:	$(TQPATHCPP)/liboctqisoc.F90
	gfortran -c -g $(CFLAGS) $(TQPATHCPP)/liboctqisoc.F90
