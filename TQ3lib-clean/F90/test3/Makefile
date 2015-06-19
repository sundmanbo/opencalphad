EXE=tqex3
FLAGS=-fbounds-check -finit-local-zero
LIBS=../../..

.PHONY : all clean

all:
	make -C $(LIBS) tq
	make $(EXE)

clean:
	rm -f *.o *.mod $(EXE)

$(EXE):
	gfortran -o $(EXE) TQ3-OU.F90 -I$(LIBS) $(LIBS)/liboctq-f90.a
