OBJS=metlib3.o tpfun4.o lukasnum.o pmod25.o matsmin.o smp3.o pmon6.o
EXE=oc2A
GF=gfortran
OPTS= -g -fbounds-check -finit-local-zero
SUBDIRS= utilities numlib models minimizer stepmapplot userif
.PHONY: clean
all:
	$(GF) -o linkoc linkocdate.F90
	./linkoc
	make $(EXE)
clean:
	rm *.a *.o *.mod oc2A	
metlib3.o:	utilities/metlib3.F90
	$(GF) -c $(OPTS) utilities/metlib3.F90

tpfun4.o:	utilities/tpfun4.F90
	$(GF) -c $(OPTS) utilities/tpfun4.F90

lukasnum.o:	numlib/lukasnum.F90
	$(GF) -c $(OPTS) numlib/lukasnum.F90

pmod25.o:	models/pmod25.F90
	$(GF) -c $(OPTS) models/pmod25.F90

matsmin.o:	minimizer/matsmin.F90
	$(GF) -c $(OPTS) minimizer/matsmin.F90

smp3.o:   stepmapplot/smp3.F90
	$(GF) -c $(OPTS) stepmapplot/smp3.F90

pmon6.o:	userif/pmon6.F90
	$(GF) -c $(OPTS) userif/pmon6.F90

$(EXE): 
	make $(OBJS)
# liboceq.a
	ar sq liboceq.a metlib3.o tpfun4.o lukasnum.o pmod25.o matsmin.o
# oc2A
	$(GF) -o $(EXE) -fbounds-check  pmain1.F90 pmon6.o smp3.o liboceq.a
