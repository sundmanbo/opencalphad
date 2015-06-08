OBJS=liboctq.o
EXE=tqex5
FLAGS=-fbounds-check -finit-local-zero

.PHONY : all clean

all:
	make $(EXE)

clean:
	rm -f *.o *.mod $(EXE)

liboctq.o:	../../liboctq.F90
	gfortran -c $(FLAGS) -I../../.. ../../liboctq.F90

$(EXE):
	make $(OBJS)
	gfortran -o $(EXE) TQ5-crfemob.F90 liboctq.o ../../../liboceq.a
