Fortran 90 test program 6

This is a simple use of a routine that can calculate an equilibrium
for a given T, P and amounts of moles of the elemenets, N(el).  The
application program does not have to handle any of the OC data
structures, all is done inside the routne.

The routine accepts an input command string which can set various
things specified in a character variable for example open a database.
The number of elements and their amouts must be specified and the
values of T and P.  It will then calculate an equilibrium and in
another character variable different results can be extracted as
specified by state variables like N(LIQUID,*) meaning the amount of
the elements in the liquid phase, or MU(*) meaning the chemical
potentials of the elements.

The routine initiates OC when a database is specified in the call and
in subsequent calls, at different T, P or amounts and where the is no
database specified, it will use the system already defined.

With some modifucations this a routine could be used in more complex
main program to simulate a reactor process where phases are separated
to different local equilibria.  It could also handle heat controlled
processes.

