OBJS=getkey.o tinyopen.o metlib3.o oclablas.o ocnum.o gtp3.o matsmin.o lmdif1lib.o smp2.o pmon6.o
EXE=oc5P

FC=gfortran
FCOPT= -O2 -fopenmp
# for debugging
#FCOPT= -fbounds-check -finit-local-zero
# no parallel
#FCOPT= -O2

#FC=ifort
#FCOPT= -check bounds -zero

# To have the command line editing and history feature on your OS 
# you must uncomment the appropriate line after the header getkey.o:

all: $(OBJS) $(EXE)

getkey.o:
	echo "Do not forget to uncomment a line below for your OS"
	# compile utilities/GETKEY for command line editing
	# uncomment the line for the kind of Linux system you have
	# Mac >>
	#gcc -c -DBSD utilities/GETKEY/getkey.c
	# Linux >>
	#gcc -c -DLinux utilities/GETKEY/getkey.c
	# other UNIX systems >>
	#gcc -c -DG77 utilities/GETKEY/getkey.c
	# CYGWIN >>
	#gcc -c -DCYGWIN utilities/GETKEY/getkey.c
	gfortran -c utilities/GETKEY/M_getkey.F90
	cp utilities/GETKEY/M_getkey.mod .

tinyopen.o:
	# compile utilities/TINYFILEDIALOGS for popup openfile windows 
	gcc -c utilities/TINYFILEDIALOGS/tinyfiledialogs.c
	gcc -c utilities/TINYFILEDIALOGS/tinyopen.c	
	gfortran -c utilities/TINYFILEDIALOGS/ftinyopen.F90
	cp utilities/TINYFILEDIALOGS/ftinyopen.mod .

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
	$(FC) -c $(FCOPT) stepmapplot/smp2.F90

pmon6.o:	userif/pmon6.F90
	# default wxt graphical driver
	# use -Dqtplt fot Qt or -Daqplt for aqua plot drivers, nothing for wxt
	# use -Dlixhlp for online helo in Linux
	$(FC) -c $(FCOPT) -Dqtplt -Dlixhlp userif/pmon6.F90

$(EXE): 
	$(FC) -o linkoc linkocdate.F90
	./linkoc

# create library liboceq.a
	ar sq liboceq.a metlib3.o oclablas.o ocnum.o gtp3.o matsmin.o lmdif1lib.o
	#
	# If getkey.o is undefined below uncomment a line above at getkey.o
	#
	$(FC) -o $(EXE) $(FCOPT) pmain1.F90 pmon6.o smp2.o getkey.o ftinyopen.o tinyopen.o tinyfiledialogs.o liboceq.a

clean:
	rm -r *.o *.mod linkoc liboceq.a $(EXE)
