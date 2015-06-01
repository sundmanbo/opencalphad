This is an advanced C++ example using TQ with isoC binding to Fortran

First copy the liboceq.a and liboceq.mov here

Then copy liboctq.F90 here

Then compile gfortran -c liboctq.F90

Then compile gfortran -c liboctqisoc.F90

Then link and compile

g++ -o tqintf -lstdc++ tqintf.cpp liboctqiso.o liboctq.o liboceq.a -lgfortran -lm

All this on linkmake.txt which should be renamed to linkmake.cmd



