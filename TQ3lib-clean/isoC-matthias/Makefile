EXE=tqintf
LIBS=../..

.PHONY : all clean

all:
	make -C $(LIBS) tq
	make $(EXE)

clean:
	rm -f *.o *.mod $(EXE)

$(EXE):
	g++ -o $(EXE) -lstdc++ $(EXE).cpp -L$(LIBS) -loctq-isoc -lgfortran -lm
