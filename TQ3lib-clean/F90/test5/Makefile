EXE=tqex5
FLAGS=-fbounds-check -finit-local-zero
LIBS=../../..

.PHONY : all clean

all:
	make -C $(LIBS) tq
	make $(EXE)

clean:
	rm -f *.o *.mod $(EXE)

$(EXE):
	gfortran -o $(EXE) TQ5-crfemob.F90 -I$(LIBS) $(LIBS)/liboctq-f90.a
