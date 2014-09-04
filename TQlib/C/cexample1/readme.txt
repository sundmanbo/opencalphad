
This is an example of a C program calling the OC TQ library

The main program is testc1.c

As I know almost nothing about C it is very primitive.  I did not even
manage to loop the calculation asking for several values of the
conditions for calculations.

IThe other files are

crfe.TDB     the database file for the Cr-Fe system
liboceq.a    a compiled fortran library with the OC routines
liboceq.mod  the module information about this library
liboctq.F90  the Fortran TQ library (preliminary version)
liboctqc.F90 a Fortran subroutine callable from C to access the TQ library
linkc1.txt   a Windows command file to compile and link the test program
readme.txt   this file
testc1.c     the main program in C

The linkc1.txt must be changed to a linkc1.cmd to be executed in a
Windows environment.  On linux or other OS one can create makefiles or
other scripts to do the same.  The liboctq.F90 is compiled and this
need access to the liboceq.mod file.  Then the liboctqc.F90 is compiled
and it need access to the liboctq.mod file.  Finally the testc1.c is
compiled and linked with liboctqc.o, liboctq.o and liboceq.a

The files liboceq.a and liboceq.mod are created when compiling the
main OC program.  They should normally be copied from that compilation
but are provided here to simplify.

The C program will initiate the OC workspace, read the crfe.TDB file and
write some information about the system.  It will then set the conditions
P=1E5 (1 bar) and N=1 (system size 1 mole).  Then it will ask for a
temperature and 800 (Kelvin) is a good value.  Then it will ask for
a mole fraction of Cr and 0.3 is a good value.

Then there will be some output from OC which will eventually be removed
but can be useful while testing.

If the calculation is successful the total Gibbs energy and the
chemical potentials of the components will be listed.

If a C expert helps one can the loop to calculate again with another
temperature and composition but with my limited understanding of C
I have failed to implement this.

The LIBOCTQ and the LIBOCTQC libraries are currently very limited,
much more can be implemented.  See also the Fortran examples which are
more elaborated.  For example they allow arrays of data to be
extracted from the OC library, I do not know how to do this in C.

Good luck and do not hesitate to ask questions and suggest changes.

bo.sundman@gmail.com
