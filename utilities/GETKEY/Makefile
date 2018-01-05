OBJS=getkey.o bintxt.o 
EXE=tin

all: $(EXE)


getkey.o:
# uncomment the appropriate line below
# for MAC
#	cc -c -DBSD getkey.c
# for Linux
#	cc -c -DLinux getkey.c
# for other UNIX systems
#	cc -c -DG77 getkey.c
# for CYGWIN on Windows
#	cc -c -DCYGWIN getkey.c
	gfortran -c M_getkey.F90

bintxt.o:
	gfortran -c bintxt.F90

$(EXE): $(OBJS) 
	#	gfortran -o $(EXE) testbintxt2.F90 bintxt.o getkey.o

