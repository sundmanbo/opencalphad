
This is an example of a C program calling the OC TQ library

The main program is testc1.c

As I know almost nothing about C it is very primitive.  I did not even
manage to loop the calculation asking for several values of the
conditions for calculations.

You must copy the following files from the OC or TQ folders

liboceq.a    compiled Fortran library with the basic OC routines
liboctq.o    compiled Fortran subroutines for the TQ interface
liboctq.mod  module information about the TQ and EQ packages

The other files are

crfe.TDB       the database file for the Cr-Fe system
liboctqc.F90   a Fortran subroutine callable from C to access the TQ library
linktqc1.txt   a Windows command file to compile and link the test program
readme.txt     this file
tqexc1.c       a main program in C for testing the C interface

The linktqc1.txt must be changed to a linktqc1.cmd to be executed in a
Windows environment.  On linux or other OS one can create makefiles or
other scripts to do the same.  

The liboctq.F90 should be compiled with access to the liboceq.mod file
compiled and this generates a liboctq.o and liboctq.mod files.  When
compiling the liboctqc.F90 it need access to liboctq.mod file.

When compiling and linking tqexc1.c it need access to liboctqc.o,
liboctq.o and liboceq.a

The files liboceq.a and liboceq.mod are created when compiling the
main OC program.  They should normally be copied from that compilation
but are provided here to simplify.

The C program will initiate the OC workspace, read the crfe.TDB file
and write some information about the system.  It will then set the
conditions P=1E5 (1 bar) and N=1 (system size 1 mole) and then it will
ask for a temperature. I suggest 800 (Kelvin) as a first value.  Then
it will ask for a mole fraction of Cr and 0.3 is a good value.

Then there will be some output from OC which will eventually be removed
but can be useful while testing.

If the calculation is successful the total Gibbs energy and the
chemical potentials of the components will be listed.

If a C expert helps one can the loop to calculate again with another
temperature and composition but with my limited understanding of C
I have failed to implement this.

The liboctq.F90 and the liboctqc.F90 libraries are currently very
limited, much more can be implemented.  See also the Fortran examples
which are more elaborated.  For example they allow arrays of data to
be extracted from the OC library, I do not know how to do this in C.

Good luck and do not hesitate to ask questions and suggest changes.

bo.sundman@gmail.com

