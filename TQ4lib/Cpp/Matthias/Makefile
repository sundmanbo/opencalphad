EXE=cpptest
LIBS=../../..
CFLAGS=-lstdc++
#UNCOMMENT THE NEXT LINE IF YOU WANT TO COMPILE OPENCALPHAD IN PARALLEL
#CFLAGS+=-fopenmp

#==============================================================================#

#Available compilation flags: all, clean

.PHONY : all clean

#Compiles OpenCalphad's OCASI examples as standalone binary executable
#software.

all:
	make $(EXE)

#Removes all binary files that were created in the compiling step.

clean:
	rm -f *.o *.mod $(EXE)

#==============================================================================#

$(EXE): $(OBJS)
	g++ -o $(EXE) $(CFLAGS) $(EXE).cpp -L$(LIBS) ../../../liboctq-isoc.a -lgfortran -lm
	g++ -c -I../.. $(EXE).cpp
