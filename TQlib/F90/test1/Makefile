OBJS=TQ1-crfe.o
INCLUDE=-I../../ -I../../../

.PHONY:
	clean

all: $(OBJS)
	gfortran $(INCLUDE) -o tq1 TQ1-crfe.o ../../liboctq.o ../../../liboceq.a

TQ1-crfe.o: TQ1-crfe.F90
	gfortran $(INCLUDE) -c -fbounds-check -finit-local-zero TQ1-crfe.F90
clean:
	rm -rf *.o tq1
