OBJS=TQ2-feni.o
INCLUDE=-I../../ -I../../../

all: $(OBJS)
	gfortran $(INCLUDE) -o tq2 TQ2-feni.o ../../liboctq.o ../../../liboceq.a

TQ2-feni.o: TQ2-feni.F90
	gfortran $(INCLUDE) -c -fbounds-check -finit-local-zero TQ2-feni.F90
	
