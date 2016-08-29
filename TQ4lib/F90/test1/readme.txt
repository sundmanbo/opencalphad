The Open Calphad Application Software Interface (OCASI)

Based on the TQ standard for interfacing thermodynamic software

Bo Sundman 29 August 2016

This is the test1 case for the Fortran OC/TQ interface.

If you are not familiar with compiling and linking software and do not
understand the intructions here please ask someone close to you for
help.  The instructions here are very brief but I am too busy to
answer questions about handling such things.

To link this example you must first compile the OC main program.  This
compilation generates two files you need: LIBOCEQ.A and
LIBOCEQPLUS.MOD.  Both of these files are needed for these
applications.

Files on this directory:

- crfe.TDB is a small database in the TDB format.

- linktqtest1 is a text file without extention which you can use as a
base for a Makefile or on Windows simply add the extention .cmd and
execute as a batch file.  In this file there are some additional
comments and instructions.  If you do not understand these instruction
plasee ask some local guru for help.

output.pdf is a PDF file with an example of the output you should
expect to have when you successfully compile, link and run this test
program.

readmel.txt is this file.  There are specific readme.txt files on the
subdirectories.

TQ1-crfe.F90 is the source program written in Fortran95/08.


