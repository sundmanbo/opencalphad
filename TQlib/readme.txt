Using the TQ library

Bo Sundman 22 December 2014

From F90 one can use the liboctc.F90 library directly.  Only a few
subroutines are implemented at present.  There are 3 test cases in the
F90 folder.

From C or C++ I do not know how to handle the F90 modules and
structures.  So I have created a single F90 function

liboctqc1

which calls the liboctq.F90 module.  Inside this subroutine the F90
structures and various tq subroutines can be accessed.  

In the call to liboctqc1 an integer variable selects with F90
subroutine to call and values are returned as arguments.

The function value returned is 0 if no error, otherwise the error
code.
