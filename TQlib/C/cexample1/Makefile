OBJ = fenitq2c.o liboctqisoc.o liboctq.o
INCLUDE = -I../../ -I../../../
FFLAGS = -g
GF=gfortran

all: $(OBJ)
	gcc -o test fenitq2c.o liboctqisoc.o liboctq.o ../../../liboceq.a -lgfortran  
    
clean:
	rm *.o 
liboctq.o:  ../../liboctq.F90
	$(GF) $(FFLAGS) -I../../../ -c ../../liboctq.F90
#liboctqc.o: ../liboctqc.F90
#	$(GF) $(FFLAGS) $(INCLUDE) -c ../liboctqc.F90
    
fenitq2c.o: fenitq2c.c
	gcc $(FFLAGS) $(INCLUDE) -c fenitq2c.c

liboctqisoc.o: ../liboctqisoc.F90
	$(GF) $(FFLAGS) $(INCLUDE) -c ../liboctqisoc.F90


