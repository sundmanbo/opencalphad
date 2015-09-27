The TQ library is to access Open Calphad from application software

Bo Sundman 27 September 2015

There is a Fortran version and a tentative iso-C version

If you are not familiar with compiling and linking software and do not
understand the intructions here please ask someone close to you for
help.  The instructions here are very brief but I am too busy to
answer questions about handling such things.

To link any of the examples you must first compile the OC main
program.  This compilation generates two files you need, liboceq.a and
liboceqplus.mod.  Both of these files are needed for these
applications.

The initial iso-C version of the library is provided by Teslos and it
has been extended by Matthias Stratmann at RUB, Germany and Christophe
Sigli at Constellium, France to handle more calls to different
subroutines.  There are different versions of this file on the isoC
directories as we are still in the development statge.

Directories:

- C has a very primitive C connection to OC

- F90 has 4 Fortran examples, all very simple cases.

- isoC-Teslos example with the firrst liboctqisoc.F90 library using
  the isoC binding for F90 and C++.  Teslos is the user name of a
  programmer at github who is unknown to me.

- isoC-Matthias example with an extended liboctqisoc.F90 library with
  connection to C++.

- isoC-Sigli is a more extensive example also using parallelization
  with an extended version of liboctq.F90 and liboctqisoc.F90.  As
  this is a very volatile development version some of these are not
  yet integrated in the basic liboctq.F90.

Files on this directory:

liboctq.F90 is the basic Fortran TQ interface.  Extended versions are
available in the subdirectories.

readme.txt is this file

To use the examples you must copy liboceq.a and liboceqplus.mod from
the main OC directory to this directory.  Then compile the liboctq.F90
with gfortran -c liboctq.F90 or compile those version available on the
subdirectores, they may contain extensions.  

You must then move the liboceq.a, liboceqplus.mod, liboctq.o and
liboctq.mod to the directory with the example programs and use the
linkexj.txt files found there.  Althoug sometimes this copying is done
automatically with this link file.  on Windows you can rename the
linkexj.txt file to .cmd and execute it as batch files.  On Linux you
can use it to create a Makefile.

For the isoC examples you must also compile the liboctqisoc.F90 file
on that directory but that is usually included in the linkexj.txt
file.  As always with C++ code there are many additional files and
libraries needed but they are normally provided when you install the
C++ compiler.


