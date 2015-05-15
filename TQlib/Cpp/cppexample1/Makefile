OBJS=liboctq.o liboctqisoc.o
EXE=tqintf
LIBS=../../..
ISOB=../..

.PHONY : all clean

all:
	make $(EXE)

clean:
	rm -f *.o *.mod $(EXE)

liboctq.o:	$(ISOB)/liboctq.F90
	gfortran -c -g -fbounds-check -finit-local-zero -I$(LIBS) $(ISOB)/liboctq.F90

liboctqisoc.o:	$(ISOB)/C/liboctqisoc.F90
	gfortran -c -g -fbounds-check -finit-local-zero -I$(LIBS) $(ISOB)/C/liboctqisoc.F90

$(EXE): $(OBJS)
	gcc -o $(EXE) -lstdc++ $(EXE).cpp liboctqisoc.o liboctq.o $(LIBS)/liboceq.a -lgfortran -lm
	gcc -c -I$(LIBS) $(EXE).cpp
