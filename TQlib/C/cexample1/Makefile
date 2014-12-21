OBJ = liboctqc.o fenitq2c.o test.o
INCLUDE = -I../../ -I../../../

all: $(OBJ)
	gcc -o test fenitq2c.o test.o ../../liboctq.o ../../../liboceq.a -lgfortran  
    
liboctqc.o: ../liboctqc.F90
	gfortran $(INCLUDE) -c ../liboctqc.F90
    
fenitq2c.o: fenitq2c.c
	gcc $(INCLUDE) -c fenitq2c.c

test.o: ../test.F90
	gfortran $(INCLUDE) -c ../test.F90


