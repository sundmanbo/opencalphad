OBJS=metlib3.o tpfun4.o lukasnum.o gtp3.o matsmin.o smp1.o pmon6.o
EXE=oc3A
FLAGS=-fbounds-check -finit-local-zero

all:
	gfortran -o linkoc linkocdate.F90
	./linkoc
	make $(EXE)
	make -C TQ2lib

clean:
	rm -f *.o *.a *.mod $(EXE) linkoc
	make -C TQ2lib clean

metlib3.o:	utilities/metlib3.F90
	gfortran -c $(FLAGS) utilities/metlib3.F90

tpfun4.o:	utilities/tpfun4.F90
	gfortran -c $(FLAGS) utilities/tpfun4.F90

lukasnum.o:	numlib/lukasnum.F90
	gfortran -c $(FLAGS) numlib/lukasnum.F90

gtp3.o:	models/gtp3.F90
	gfortran -c $(FLAGS) models/gtp3.F90

matsmin.o:	minimizer/matsmin.F90
	gfortran -c $(FLAGS) minimizer/matsmin.F90

smp1.o:		stepmapplot/smp1.F90
	gfortran -c $(FLAGS) stepmapplot/smp1.F90

pmon6.o:	userif/pmon6.F90
	gfortran -c $(FLAGS) userif/pmon6.F90

$(EXE): 
	make $(OBJS)
	ar sq liboceq.a metlib3.o tpfun4.o lukasnum.o gtp3.o matsmin.o
	gfortran -o $(EXE) -fbounds-check  pmain1.F90 pmon6.o smp1.o liboceq.a
