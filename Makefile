OBJS=getkey.o M_getkey.o ftinyopen.o tinyopen.o tinyfiledialogs.o metlib3.o oclablas.o ocnum.o gtp3.o matsmin.o lmdif1lib.o smp2.o pmon6.o
LIBS=liboctq.o liboctqisoc.o liboctqcpp.o
EXE=oc5P

FC=gfortran
C=gcc
CPP=g++
FCOPT= -O2 -fopenmp -fPIC
# for debugging
#FCOPT= -fbounds-check -finit-local-zero
# no parallel
#FCOPT= -O2
#FC=ifort
#FCOPT= -check bounds -zero

EX1PATH = ./TQ4lib/F90/crfe
EX2PATH = ./TQ4lib/F90/feni
EX3PATH = ./TQ4lib/Cpp/Matthias/crfe
EX4PATH = ./TQ4lib/Cpp/Matthias/feni
EX5PATH = ./TQ4lib/Cpp/Matthias

#==============================================================================#

#Available compilation flags: all, OCASI, OCASIEXAMPLES, clean

.PHONY : all OCASI OCASIEXAMPLES clean

#Compiles OpenCalphad to use as standalone Thermodynamic Equilibrium Calculation
#software.
# ************************************
# OC now requires GNUPLOT 5.0 or later
# ************************************

# To have the command line editing and history feature on your OS 
# you must uncomment the appropriate line after the header getkey.o:

all: $(OBJS) $(EXE)

#Compiles the OCASI interfaces of OpenCalphad, so third party software can
#interact with OpenCalphad. Interfaces are provided in C++ and Fortran, and
#additional interfaces are available in C and Python.

OCASI:
	make $(OBJS) $(LIBS)
	ar sq liboctq-f90.a liboctq.o $(OBJS)
	ar sq liboctq-isoc.a liboctqisoc.o liboctq.o $(OBJS)
	ar sq liboctqcpp.a liboctqcpp.o liboctqisoc.o liboctq.o $(OBJS)
	$(CPP) -shared liboctqcpp.o liboctqisoc.o liboctq.o $(OBJS) -o $(EXE)_OCASI.so -lgfortran

#Compiles the OCASI interface and various examples.

OCASIEXAMPLES:
	make OCASI
	$(FC)  -o $(EX1PATH)/tqex1 $(EX1PATH)/TQ1-crfe.F90 liboctq-f90.a -fopenmp
	$(FC)  -o $(EX2PATH)/tqex2 $(EX2PATH)/TQ2-feni.F90 liboctq-f90.a -fopenmp
	$(CPP) -o $(EX3PATH)/tqex1 $(EX3PATH)/tqex1.cpp liboctqcpp.a -lgfortran -fopenmp
	$(CPP) -o $(EX4PATH)/tqex2 $(EX4PATH)/tqex2.cpp liboctqcpp.a -lgfortran -fopenmp
	$(CPP) -o $(EX5PATH)/tqex3 $(EX5PATH)/tqex3.cpp liboctqcpp.a -lgfortran -fopenmp

clean:
	rm -r *.mod *.a $(LIBS) $(OBJS) linkoc $(EXE)_OCASI.so $(EXE) $(EX1PATH)/tqex1 $(EX2PATH)/tqex2 $(EX3PATH)/tqex1 $(EX4PATH)/tqex2 $(EX5PATH)/tqex3

#==============================================================================#

# To have the command line editing and history feature on your OS
# you must uncomment the appropriate line after the header getkey.o:

getkey.o:
	echo "Do not forget to uncomment a line below for your OS"
	# compile utilities/GETKEY for command line editing
	# uncomment the line for the kind of Linux system you have
	# Mac >>
	#$(C) -c $(FCOPT)-DBSD utilities/GETKEY/getkey.c
	# Linux >>
	$(C) -c $(FCOPT) -DLinux utilities/GETKEY/getkey.c
	# other UNIX systems >>
	#$(C) -c $(FCOPT) -DG77 utilities/GETKEY/getkey.c
	# CYGWIN >> 
	#$(C) -c $(FCOPT) -DCYGWIN utilities/GETKEY/getkey.c

M_getkey.o:
	$(FC) -c $(FCOPT) utilities/GETKEY/M_getkey.F90

tinyfiledialogs.o:
	$(C)  -c $(FCOPT) utilities/TINYFILEDIALOGS/tinyfiledialogs.c

tinyopen.o:
	$(C)  -c $(FCOPT) utilities/TINYFILEDIALOGS/tinyopen.c

ftinyopen.o:
	$(FC) -c $(FCOPT) utilities/TINYFILEDIALOGS/ftinyopen.F90

metlib3.o:	utilities/metlib3.F90
	# lixed for command line editing, tinyfd for open files
	# lixhlp for browser help
	$(FC) -c $(FCOPT) -Dlixed -Dtinyfd -Dlixhlp utilities/metlib3.F90

oclablas.o:	numlib/oclablas.F90
	$(FC) -c $(FCOPT) numlib/oclablas.F90

ocnum.o:	numlib/ocnum.F90
	$(FC) -c $(FCOPT) numlib/ocnum.F90

lmdif1lib.o:      numlib/lmdif1lib.F90
	$(FC) -c $(FCOPT) numlib/lmdif1lib.F90

gtp3.o:	models/gtp3.F90
	$(FC) -c $(FCOPT) models/gtp3.F90

matsmin.o:	minimizer/matsmin.F90
	$(FC) -c $(FCOPT) minimizer/matsmin.F90

smp2.o:		stepmapplot/smp2.F90
	# Remove -Dnotwin if compiled on windows
	$(FC) -c $(FCOPT) -Dnotwin stepmapplot/smp2.F90

pmon6.o:	userif/pmon6.F90
	# default wxt graphical driver
	# use -Dqtplt fot Qt or -Daqplt for aqua plot drivers, nothing for wxt
	# use -Dlixhlp for online help in Linux
	$(FC) -c $(FCOPT) -Dqtplt -Dlixhlp userif/pmon6.F90
	
liboctq.o:	./TQ4lib/Cpp/liboctq.F90
	$(FC) -c -g $(FCOPT) ./TQ4lib/Cpp/liboctq.F90

liboctqisoc.o:	./TQ4lib/Cpp/Matthias/liboctqisoc.F90
	$(FC) -c -g $(FCOPT) ./TQ4lib/Cpp/Matthias/liboctqisoc.F90

liboctqcpp.o:	./TQ4lib/Cpp/Matthias/liboctqcpp.cpp
	$(CPP) -c -g $(FCOPT) ./TQ4lib/Cpp/Matthias/liboctqcpp.cpp

$(EXE): 
	$(FC) -o linkoc linkocdate.F90
	./linkoc
	# create library liboceq.a
	ar sq liboceq.a metlib3.o oclablas.o ocnum.o gtp3.o matsmin.o lmdif1lib.o
	# If getkey.o is undefined below uncomment a line above at getkey.o
	$(FC) -o $(EXE) $(FCOPT) pmain1.F90 $(OBJS) liboceq.a
