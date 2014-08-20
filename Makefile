OBJS=metlib3.o tpfun4.o lukasnum.o pmod25.o matsmin.o smp1.o pmon6.o
EXE=oc1A

all:
	gfortran -o linkoc linkocdate.F90
	./linkoc
	make $(EXE)


metlib3.o:	utilities/metlib3.F90
	gfortran -c -fbounds-check utilities/metlib3.F90

tpfun4.o:	utilities/tpfun4.F90
	gfortran -c -fbounds-check utilities/tpfun4.F90

lukasnum.o:	numlib/lukasnum.F90
	gfortran -c -fbounds-check numlib/lukasnum.F90

pmod25.o:	models/pmod25.F90
	gfortran -c -fbounds-check models/pmod25.F90

matsmin.o:	minimizer/matsmin.F90
	gfortran -c -fbounds-check minimizer/matsmin.F90

smp1.o:		stepmapplot/smp1.F90
	gfortran -c -fbounds-check stepmapplot/smp1.F90

pmon6.o:	userif/pmon6.F90
	gfortran -c -fbounds-check userif/pmon6.F90

$(EXE): $(OBJS) linkocdate.F90
	gfortran -o $(EXE) -fbounds-check  pmain1.F90 $(OBJS)
