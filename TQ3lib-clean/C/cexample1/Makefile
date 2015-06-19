EXE=tqexc1
LIBS=../../..

.PHONY : all clean

all:
	make -C $(LIBS) tq
	make $(EXE)

clean:
	rm -f *.o *.mod $(EXE)

$(EXE):
	gcc -c tqexc1.c
	gfortran -o $(EXE) tqexc1.o -L$(LIBS) -loctq-c
