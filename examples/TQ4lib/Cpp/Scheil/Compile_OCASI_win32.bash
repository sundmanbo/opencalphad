#!/bin/bash
#
# You must adjust this command file to your local environment
#
COMPILER=i686-w64-mingw32
GFORTRAN=${COMPILER}-gfortran
GCPP=${COMPILER}-g++

set -x

rm *.o

# Copy the libraries from OC (compiled with parallel)
cp ../F90/{liboceq.a,liboceq.mod,/liboceqplus.mod} .

${GFORTRAN} -c -fopenmp liboctq.F90
${GFORTRAN} -c -fopenmp liboctqisoc.F90
${GCPP} -o scheil.exe -static -fopenmp -lstdc++ example_OCASI.cpp liboctqisoc.o liboctq.o liboceq.a -lgfortran -lm -lquadmath 
