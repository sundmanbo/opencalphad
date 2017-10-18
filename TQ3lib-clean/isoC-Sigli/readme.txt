This is an example using OC and the OC application software interface

Provided by Christophe Sigli November 2015

The example calculates liquidus/solidus/transformation temperatures
and a Scheil solidification simulation using the COST507 database for
aluminium.

NOTE THAT THE COST507 DATABASE IS NOT VALIDATED FOR MULTICOMPONENT ALLOYS.
The results must not be taken as accurate.

To compile the program (on Windows) change the linkpm.txt file to a
cmd file.  On UNIX change it to a Makefile.  It uses the compiled
versions of liboceq.a from the main OC directory.

The main program sigli.exe requires an input file, here the input.txt
is provied.  If you want to change alloy or what the program should
calculate you must change this input file.

The calculated results are listed on the screen and in the output.txt file

