OBJS=tinyfiledialogs.o tinyopen.o ftinyopen.o
EXE=gopen

all: $(EXE)

tinyfiledialogs.o:
	gcc -c tinyfiledialogs.c

tinyopen.o:
	gcc -c tinyopen.c

ftinyopen.o:
	gfortran -c ftinyopen.F90

$(EXE): $(OBJS)
	gfortran -o $(EXE) main.F90 ftinyopen.o tinyopen.o tinyfiledialogs.o

clean: 
	rm *.o *.exe
