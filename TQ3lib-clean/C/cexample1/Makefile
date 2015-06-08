OBJS=liboctq.o liboctqc.o 
EXE=tqexc1
FLAGS=-fbounds-check -finit-local-zero

.PHONY : all clean

all:
	make $(EXE)

clean:
	rm -f *.o *.mod $(EXE)

liboctq.o:	../../liboctq.F90
	gfortran -c $(FLAGS) -I../../.. ../../liboctq.F90

liboctqc.o:	liboctqc.F90
	gfortran -c $(FLAGS) liboctqc.F90

$(EXE):
	make $(OBJS)
	gcc -c tqexc1.c
	gfortran -o $(EXE) tqexc1.o liboctqc.o liboctq.o ../../../liboceq.a
