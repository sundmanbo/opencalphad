
This is an example of a C program calling the OC TQ library

The main program is testc1.c

As I know almost nothing about C it is very primitive.  I did not even
manage to loop the calculation asking for several values of the
conditions for calculations.

You must copy or compile the following files from the OC or TQ
folders.  It is done by the linktqc1 command file.

liboceq.a        compiled Fortran library with the basic OC routines
liboctqplus.mod  module information about the EQ packages
liboctq.F90      the Fortran library with the basic TQ routines

The files are

crfe.TDB       the database file for the Cr-Fe system
liboctqc.F90   a Fortran subroutine callable from C to access the TQ library
linktqc1.txt   a Windows command file to compile and link the test program
readme.txt     this file
tqexc1.c       a main program in C for testing the C interface

The linktqc1.txt must be changed to a linktqc1.cmd to be executed in a
Windows environment.  On linux or other OS one can create a Makefile
or other scripts to do the same.

The liboctq.F90 must be compiled with access to the liboceqplus.mod
file compiled and this generates a liboctq.o and a liboctq.mod file.
When compiling the liboctqc.F90 it need access to liboctq.mod file.

When compiling and linking tqexc1.c it need access to liboctqc.o,
liboctq.o and liboceq.a

The files liboceq.a and liboceqplus.mod are created when compiling the
main OC program.  They are copied from that compilation.

The C program will initiate the OC workspace, read the crfe.TDB file
and write some information about the system.  It will then set the
conditions P=1E5 (1 bar) and N=1 (system size 1 mole) and then it will
ask for a temperature. I suggest 800 (Kelvin) as a first value.  Then
it will ask for a mole fraction of Cr and 0.3 is a good value.

Then there will be some output from OC which will eventually be removed
but can be useful while testing.

If the calculation is successful the total Gibbs energy and the
chemical potentials of the components will be listed.

If a C expert helps you to read input the program can loop to
calculate again with another temperature and composition but with my
limited understanding of C I have failed to implement this.

The liboctq.F90 and the liboctqc.F90 libraries are currently very
limited, much more can be implemented.  See also the Fortran examples
and isoC examples for C++ which are more elaborated.  For example they
allow arrays of data to be extracted from the OC library, I do not
know how to do this in C.

Good luck
