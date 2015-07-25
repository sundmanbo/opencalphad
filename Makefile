OBJS=metlib3.o tpfun4.o lukasnum.o gtp3.o matsmin.o smp1.o pmon6.o
VER=oc3A
TQPATH=./TQ3lib-clean
TQPATHC=$(TQPATH)/C/cexample1
TQPATHCPP=$(TQPATH)/isoC-matthias
CFLAGS=-fbounds-check -finit-local-zero

#==============================================================================#

#Available compilation flags: all, tq, examples, clean

#Compiles OpenCalphad to use as standalone Thermodynamic Equilibrium Calculation
#software.

all:
	gfortran -o linkoc linkocdate.F90
	./linkoc
	make $(OBJS) 
	ar sq lib$(VER)eq.a $(OBJS)
	gfortran -o $(VER) $(CFLAGS) pmain1.F90 lib$(VER)eq.a

#Compiles the TQ interfaces of OpenCalphad, so third party software can interact
#with OpenCalphad. Interfaces are provided in C++ and Fortran, and additional
#interfaces are available in C and Python.

tq: 
	make $(OBJS) liboctq.o liboctqc.o liboctqisoc.o
	ar sq liboctq-f90.a liboctq.o $(OBJS) liboctq.o
	ar sq liboctq-c.a liboctqc.o liboctq.o $(OBJS)
	ar sq liboctq-isoc.a liboctqisoc.o liboctq.o $(OBJS)

#Compiles the TQ interface and various examples. 

examples:
	make -C $(TQPATH)

#Removes all binary files that were created in the compiling step.

clean:
	rm -f *.a *.o *.mod $(VER)
	make -C $(TQPATH) clean

#==============================================================================#

metlib3.o:	utilities/metlib3.F90
	gfortran -c $(CFLAGS) utilities/metlib3.F90

tpfun4.o:	utilities/tpfun4.F90
	gfortran -c $(CFLAGS) utilities/tpfun4.F90

lukasnum.o:	numlib/lukasnum.F90
	gfortran -c $(CFLAGS) numlib/lukasnum.F90

gtp3.o:	models/gtp3.F90
	gfortran -c $(CFLAGS) models/gtp3.F90

matsmin.o:	minimizer/matsmin.F90
	gfortran -c $(CFLAGS) minimizer/matsmin.F90

smp1.o:		stepmapplot/smp1.F90
	gfortran -c $(CFLAGS) stepmapplot/smp1.F90

pmon6.o:	userif/pmon6.F90
	gfortran -c $(CFLAGS) userif/pmon6.F90
	
liboctq.o:	$(TQPATH)/liboctq.F90
	gfortran -c -g $(CFLAGS) $(TQPATH)/liboctq.F90

liboctqc.o:	$(TQPATHC)/liboctqc.F90
	gfortran -c $(CFLAGS) $(TQPATHC)/liboctqc.F90

liboctqisoc.o:	$(TQPATHCPP)/liboctqisoc.F90
	gfortran -c -g $(CFLAGS) $(TQPATHCPP)/liboctqisoc.F90
	
