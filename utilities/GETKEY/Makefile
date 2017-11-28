OBJS=getkey.o bintxt.o 
EXE=tin

all: $(EXE)


getkey.o:
	cc -c -DBSD getkey.c
	gfortran -c M_getkey.F90

bintxt.o:
	gfortran -c bintxt.F90

$(EXE): $(OBJS) 
	#	gfortran -o $(EXE) testbintxt2.F90 bintxt.o getkey-BSD.o

