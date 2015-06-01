The TQ library is to access Open Calphad from application software

Bo Sundman 20 May 2015

There is a Fortran version and a tentative iso-C version

The main library is liboctq.F90, this should be compiled using
liboceq.mov created when compiling the mian OC library liboceq.F90 For
linking application one also need liboceq.a.  Both of these are
generated when compiling and linking the main OC program.

The initial iso-C version of the library is provided by Teslos and it
has been extended to handle more calls.

Content:

Directoties:

- C has a very primitive C connection to OC

- F90 has 4 Fortran examples.  The 4th is new showing the use of some
additional features.

- isoC-Teslos example with liboctqisoc.F90 library with connection to C++

- isoC-Matthias example with liboctqisoc.F90 library with connection to C++

Files:

liboctq.F90 is the Fortran TQ interface

readme.txt is this file

To use the examples you must copy liboceq.a and liboceq.mod from the
main OC directory to this directory.  Then compile liboctq.F90 with
gfortran -c liboctq.F90.  You must then move the liboceq.a, liboctq.o
and liboctq.mod to the directory with the example programs and use the
linkexj.txt files found there (on Windows rename the file to .cmd to
execute them as batch files).

For the isoC examples you must also compile the liboctqisoc.F90 module.






