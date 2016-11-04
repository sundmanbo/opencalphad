OBJS=liboctq.o liboctqisoc.o
EXE=tqintf
LIBS=../..

.PHONY : all clean

all:
	make $(EXE)

clean:
	rm -f *.o *.mod $(EXE)

liboctq.o:	../liboctq.F90
	gfortran -c -g -I$(LIBS) ../liboctq.F90

liboctqisoc.o:	liboctqisoc.F90
	gfortran -c -g -I$(LIBS) liboctqisoc.F90

$(EXE): $(OBJS)
	gcc -o $(EXE) -lstdc++ $(EXE).cpp liboctqisoc.o liboctq.o ../../liboceq.a -lgfortran -lm
	gcc -c -I../.. $(EXE).cpp
