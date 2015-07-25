EXE=feni-tqiso
LIBS=../..

.PHONY : all clean

all:
	make -C $(LIBS) tq
	make $(EXE)

clean:
	rm -f *.o *.mod $(EXE)

$(EXE):
	gcc -o $(EXE) feni-tqiso.c -L$(LIBS) -loctq-isoc -lgfortran -lm./
