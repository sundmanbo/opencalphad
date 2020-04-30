This is a simple C++ example using TQ with isoC binding to Fortran

You must first nave installed OC so you have the libraries liboceq.a
and liboceqplus.mod.  

On Windows use linkmake, first rename it to linkmake.cmd and then execute it

On UNIX systems use Makefile

These command files will first copy the liboceq.a and liboceq.mod to
this directory.  It assumes they are two directories above as in the
zip file you downloaded.

The the Fortran TQ library liboctq.F90 is copied here from one
directory above.

Then it will compile gfortran -c liboctq.F90

Then it will compile gfortran -c liboctqisoc.F90

Then it will link and compile the main program tqintf.cpp

g++ -o tqintf -lstdc++ tqintf.cpp liboctqiso.o liboctq.o liboceq.a -lgfortran -lm




