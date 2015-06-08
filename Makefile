OBJS=metlib3.o tpfun4.o lukasnum.o gtp3.o matsmin.o smp1.o pmon6.o
EXE=oc3A

all:
	gfortran -o linkoc linkocdate.F90
	./linkoc
	make $(EXE)

examples:
	make -C ./TQ3lib-clean

clean:
	rm -f *.o *.mod $(EXE)

examplesclean:
	make -C ./TQ3lib-clean clean

metlib3.o:	utilities/metlib3.F90
	gfortran -c -fbounds-check  -finit-local-zero utilities/metlib3.F90

tpfun4.o:	utilities/tpfun4.F90
	gfortran -c -fbounds-check  -finit-local-zero utilities/tpfun4.F90

lukasnum.o:	numlib/lukasnum.F90
	gfortran -c -fbounds-check  -finit-local-zero numlib/lukasnum.F90

gtp3.o:	models/gtp3.F90
	gfortran -c -fbounds-check -finit-local-zero models/gtp3.F90

matsmin.o:	minimizer/matsmin.F90
	gfortran -c -fbounds-check  -finit-local-zero minimizer/matsmin.F90

smp1.o:		stepmapplot/smp1.F90
	gfortran -c -fbounds-check  -finit-local-zero stepmapplot/smp1.F90

pmon6.o:	userif/pmon6.F90
	gfortran -c -fbounds-check  -finit-local-zero userif/pmon6.F90

$(EXE): 
	make $(OBJS)
# liboceq.a
	ar sq liboceq.a metlib3.o tpfun4.o lukasnum.o gtp3.o matsmin.o
# oc2A
	gfortran -o $(EXE) -fbounds-check  pmain1.F90 pmon6.o smp1.o liboceq.a
