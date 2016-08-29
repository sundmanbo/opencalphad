The Open Calphad Application Software Interface (OCASI)

Based on the TQ standard for interfacing thermodynamic software

Bo Sundman 29 August 2016

There is a Fortran version and a tentative iso-C version for C++.  
In the future it may be possible to merge these.

If you are not familiar with compiling and linking software and do not
understand the intructions here please ask someone close to you for
help.  The instructions here are very brief but I am too busy to
answer questions about handling such things.

To link any of the examples you must first compile the OC main
program.  This compilation generates two files you need:
LIBOCEQ.A and LIBOCEQPLUS.MOD.  
Both of these files are needed for these applications.

The initial iso-C version of the library was provided by Teslos in
2014 and it has been extended by Matthias Stratmann at RUB, Germany
and Christophe Sigli at Constellium, France to handle more calls to
different subroutines.  As things are still under development there
may be slightly different versions on various subdirectories.

Files on this directory:

LIBOCTQ.F90 is the basic Fortran TQ interface.  It calls many
subroutines and functions in the LIBOCEQ.A library.

README.TXT is this file.  There are specific readme.txt files on the
subdirectories.

Subdirectories:

- F90 has one Fortran examples, more will be added shortly.

- Cpp has one C++ example provided by Matthias Stratmann at RUB,
Germany.  More will be added shortly.



