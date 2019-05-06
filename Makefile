OBJS=M_getkey.o tinyopen.o metlib3.o oclablas.o ocnum.o gtp3.o matsmin.o lmdif1lib.o smp2.o pmon6.o
EXE=oc5P

FC=gfortran
C=gcc
CPP=g++
#FCOPT= -O2 -fopenmp
# for debugging
#FCOPT= -fbounds-check -finit-local-zero
# no parallel
#FCOPT= -O2
# for OCASI
FCOPT= -O2 -fPIC 
#-fopenmp 
#FC=ifort
#FCOPT= -check bounds -zero

# To have the command line editing and history feature on your OS 
# you must uncomment the appropriate line after the header getkey.o:

all: $(OBJS) $(EXE)

M_getkey.o:
	echo "Do not forget to uncomment a line below for your OS"
	# compile utilities/GETKEY for command line editing
	# uncomment the line for the kind of Linux system you have
	# Mac >>
	#$(C) -c $(FCOPT)-DBSD utilities/GETKEY/getkey.c
	# Linux >>
	#$(C) -c $(FCOPT) -DLinux utilities/GETKEY/getkey.c
	# other UNIX systems >>
	#$(C) -c $(FCOPT)-DG77 utilities/GETKEY/getkey.c
	# CYGWIN >>
	#gcc -c $(FCOPT)-DCYGWIN utilities/GETKEY/getkey.c
	$(FC) -c $(FCOPT)utilities/GETKEY/M_getkey.F90

tinyopen.o:
	# compile utilities/TINYFILEDIALOGS for popup openfile windows 
	$(C) -c $(FCOPT) utilities/TINYFILEDIALOGS/tinyfiledialogs.c
	$(C) -c $(FCOPT) utilities/TINYFILEDIALOGS/tinyopen.c	
	$(FC) -c $(FCOPT) utilities/TINYFILEDIALOGS/ftinyopen.F90

metlib3.o:	utilities/metlib3.F90
	# lixed for command line editing, tinyfd for open files
	# lixhlp for browser help
	$(FC) -c $(FCOPT) -Dlixhlp utilities/metlib3.F90

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
	#
	# If getkey.o is undefined below uncomment a line above at getkey.o
	#
	$(FC) -o $(EXE) $(FCOPT) pmain1.F90 pmon6.o smp2.o getkey.o ftinyopen.o tinyopen.o tinyfiledialogs.o liboceq.a

OCASI:
	make $(OBJS) liboctq.o liboctqisoc.o liboctqcpp.o
	ar sq liboctq-f90.a liboctq.o $(OBJS) liboctq.o
	ar sq liboctq-isoc.a liboctqisoc.o liboctq.o $(OBJS)
	ar sq liboctqcpp.a liboctqcpp.o liboctqisoc.o liboctq.o $(OBJS)
	ld -shared $(OBJS) -o $(EXE)_OCASI.so

OCASIEXAMPLES:
	make OCASI
	$(FC) -o ./TQ4lib/F90/crfe/tqex1 ./TQ4lib/F90/crfe/TQ1-crfe.F90 liboctq-f90.a
	$(FC) -o ./TQ4lib/F90/feni/tqex2 ./TQ4lib/F90/feni/TQ2-feni.F90 liboctq-f90.a
	$(CPP) -o ./TQ4lib/Cpp/Matthias/crfe/tqex1 ./TQ4lib/Cpp/Matthias/crfe/tqex1.cpp liboctqcpp.a -lgfortran
	$(CPP) -o ./TQ4lib/Cpp/Matthias/feni/tqex2 ./TQ4lib/Cpp/Matthias/feni/tqex2.cpp liboctqcpp.a -lgfortran
	$(CPP) -o ./TQ4lib/Cpp/Matthias/tqex3 ./TQ4lib/Cpp/Matthias/tqex3.cpp liboctqcpp.a -lgfortran
	
clean:
	rm -r *.o *.mod linkoc liboceq.a $(EXE)
