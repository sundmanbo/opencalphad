#
# Modified 2020.04.25 for the new direcory structure.  Bo Sundman
#
OBJS=getkey.o M_getkey.o ftinyopen.o tinyopen.o tinyfiledialogs.o metlib4.o oclablas.o ocnum.o minpack1.o  gtp3.o matsmin.o smp2.o pmon6.o
LIBS=liboctq.o liboctqisoc.o liboctqcpp.o
EXE=oc6P

#
# IMPORTANT 1. select getkey version 
#           2. select GNUPLOT terminal
#
#	echo "Do not forget to uncomment a line for your OS"
#
#=============================================================================#
# original provided by Matthias Strathmann including OC examples
# NOTE getkey.o : you have to select which getkey routine to compile
#=============================================================================#

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

EX1PATH = ./examples/TQ4lib/F90/crfe
EX2PATH = ./examples/TQ4lib/F90/feni
EX3PATH = ./examples/TQ4lib/Cpp/Matthias/crfe
EX4PATH = ./examples/TQ4lib/Cpp/Matthias/feni
EX5PATH = ./examples/TQ4lib/Cpp/Matthias

#=============================================================================#

#Available compilation flags: all, OCASI, OCASIEXAMPLES, clean

.PHONY : all OCASI OCASIEXAMPLES clean

#Compiles OpenCalphad to use as standalone Thermodynamic Equilibrium Calculation
#software.
# ************************************
# OC now requires GNUPLOT 5.2 or later
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

# IMPORTRANT 1:
# To have the command line editing and history feature on your OS
# you must uncomment the appropriate line after the header getkey.o:

getkey.o:
	echo "Do not forget to uncomment a line below for your OS"
	# compile utilities/GETKEY for command line editing
	# uncomment the line for the kind of Linux system you have
	# Mac >>
	#$(C) -c $(FCOPT) -DBSD src/utilities/GETKEY/getkey.c
	# Linux >>
	#$(C) -c $(FCOPT) -DLinux src/utilities/GETKEY/getkey.c
	# other UNIX systems >>
	#$(C) -c $(FCOPT) -DG77 src/utilities/GETKEY/getkey.c
	# CYGWIN >> 
	#$(C) -c $(FCOPT) -DCYGWIN src/utilities/GETKEY/getkey.c

# If you have not uncommented any getkey.c line above COMMENT next line
# and also remove the -Dlixed option for the metlib4.F90
M_getkey.o:
	$(FC) -c $(FCOPT) src/utilities/GETKEY/M_getkey.F90

tinyfiledialogs.o:
	$(C)  -c $(FCOPT) src/utilities/TINYFILEDIALOGS/tinyfiledialogs.c

tinyopen.o:
	$(C)  -c $(FCOPT) src/utilities/TINYFILEDIALOGS/tinyopen.c

ftinyopen.o:
	$(FC) -c $(FCOPT) src/utilities/TINYFILEDIALOGS/ftinyopen.F90

metlib4.o:	src/utilities/metlib4.F90
	# lixed for command line editing,
	# tinyfd for open files
	# lixhlp for browser help on Linux and MacOS
	$(FC) -c $(FCOPT) -Dlixed -Dtinyfd -Dlixhlp src/utilities/metlib4.F90

oclablas.o:	src/numlib/oclablas.F90
	$(FC) -c $(FCOPT) src/numlib/oclablas.F90

ocnum.o:	src/numlib/ocnum.F90
	$(FC) -c $(FCOPT) src/numlib/ocnum.F90

minpack1.o:      src/numlib/minpack1.F90
	$(FC) -c $(FCOPT) src/numlib/minpack1.F90

gtp3.o:	src/models/gtp3.F90
	$(FC) -c $(FCOPT) src/models/gtp3.F90

matsmin.o:	src/minimizer/matsmin.F90
	$(FC) -c $(FCOPT) src/minimizer/matsmin.F90

smp2.o:		src/stepmapplot/smp2.F90
	# Remove -Dnotwin if compiled on Windows (for spawning)
	$(FC) -c $(FCOPT) -Dnotwin src/stepmapplot/smp2.F90

# IMPORTANT 2: select GNUPLOT terminal
pmon6.o:	src/userif/pmon6.F90
	# default wxt graphical driver
	# use -Dqtplt for Qt or (also smaller window)
	# use -Daqplt for aqua plot drivers (also smaller window)
	# use -Dx11 for X11 plot drivers DEFAULT
	# use -Dlixhlp for online help in Linux
	# use -Dmachlp for online help in MacOS (browser)
	$(FC) -c $(FCOPT) -Dx11 -Dlixhlp src/userif/pmon6.F90

liboctq.o:	./examples//TQ4lib/Cpp/liboctq.F90
	$(FC) -c -g $(FCOPT) ./examples/TQ4lib/Cpp/liboctq.F90

liboctqisoc.o:	./examples/TQ4lib/Cpp/Matthias/liboctqisoc.F90
	$(FC) -c -g $(FCOPT) ./examples/TQ4lib/Cpp/Matthias/liboctqisoc.F90

liboctqcpp.o:	./examples/TQ4lib/Cpp/Matthias/liboctqcpp.cpp
	$(CPP) -c -g $(FCOPT) ./examples/TQ4lib/Cpp/Matthias/liboctqcpp.cpp

$(EXE): 
	# Add date of linking to main program
	$(FC) -o linkoc src/linkocdate.F90
	./linkoc
	# create library liboceq.a
	mkdir -p libs
	ar sq libs/liboceq.a metlib4.o oclablas.o ocnum.o gtp3.o matsmin.o minpack1.o

	# If getkey.o is undefined below
	# you have forgotten to uncomment a line above at getkey.o !!
	# static: $(FC) -o $(EXE) $(FCOPT) -static-libgfortran pmain1.F90 $(OBJS) liboceq.a
	#$(FC) -o $(EXE) $(FCOPT) src/pmain1.F90 $(OBJS) libs/liboceq.a
	$(FC) -o $(EXE) $(FCOPT) src/pmain1.F90 pmon6.o smp2.o ftinyopen.o tinyopen.o tinyfiledialogs.o getkey.o libs/liboceq.a

