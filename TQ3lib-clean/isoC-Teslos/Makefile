OBJS=liboctq.o liboctqisoc.o
EXE=feni-tqiso

.PHONY : all clean

all:
	make $(EXE)

clean:
	rm -f *.o *.mod $(EXE)

liboctq.o:	../liboctq.F90
	gfortran -c -g -fbounds-check -finit-local-zero -I../.. ../liboctq.F90

liboctqisoc.o:	liboctqisoc.F90
	gfortran -c -g -fbounds-check -finit-local-zero -I../.. liboctqisoc.F90

$(EXE): $(OBJS)
	gcc -o $(EXE) feni-tqiso.c liboctqisoc.o liboctq.o ../../liboceq.a -lgfortran -lm
	gcc -c -I../.. feni-tqiso.c


