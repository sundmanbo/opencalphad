OBJS=fenitq2.o
INCLUDE=-I../../ -I../../../

all: $(OBJS)
	gfortran $(INCLUDE) -o tq2 fenitq2.o ../../liboctq.o ../../../liboceq.a

fenitq2.o: fenitq2.F90
	gfortran $(INCLUDE) -c -fbounds-check -finit-local-zero fenitq2.F90
	
