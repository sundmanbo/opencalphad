#include "octqc.h"
#include <stdio.h>
#include <stdlib.h>
#define MAXEL 10
#define MAXPH 10

extern void c_tqini(int, void** );
extern void c_tqrfil(char *, void** );
extern void c_tqgcom(int*, char* comp[24], void**);
extern void c_tqrpfil(char *, int, char**, void**);
extern void c_tqgnp(int*,void**);
extern void c_tqgpn(int*,char **, void**);
extern void examine_gtp_equilibrium_data(void*);
//extern void c_tqgnp(int, gtp_equilibrium_data** );


int main(int argc, char **argv) 
{
  int i, nel=2, nph, ics, ip;
  void *ceq = 0;
  char cmpname[MAXEL][24];
  char phnames[MAXPH][24];
  char target[60];

  double value, tp[2], dummy, mel[MAXEL], mf[MAXEL], volume;
  double xf[MAXEL], pxf[10*MAXPH], npf[MAXPH], mu[MAXEL];
  gtp_equilibrium_data *dceq;
  // set some defaults
	int n=0;
	char* filename="FENI.TDB";
	char* selel[2] = {"FE","NI"};

	//initilize tq 
	c_tqini(n, &ceq);
  	c_tqrfil(filename, &ceq);
	dceq = (gtp_equilibrium_data *) ceq;
	printf("Equilibrium : %s\n", dceq->eqname);
  //	c_tqrpfil(filename, nel, selel, &ceq);
  	c_tqgcom(&n, cmpname, &ceq);
  // read database file
  //examine_gtp_equilibrium_data(ceq);
	
  // find out about the system
  // number of elements and their names
  	
  // number of phases and their names
  	c_tqgnp(&n,&ceq);
	printf("and %i phases: ",n);
	for(i = 1; i <= n; i++) {
		c_tqgpn(&i, phnames[i], &ceq);
		printf(" %s,",phnames[i]);
	}
  // set default values
  	tp[0] = 1.0e3;
	tp[1] = 1.0e5;

	for (i = 1; i <= nel; i++)
		xf[i] = 1.0/(double)nel;

  // ask for conditions
  	printf("Give conditions: ");
	dummy = tp[0];
	
	if(tp[0] < 1.0) {
		printf("Temperature must be larger than 1K\n");
		tp[0] = 1.0;
	}
	dummy = tp[1];
	if(tp[1] < 1.0) {
		printf("Pressure must be larger than 1Pa\n");
		tp[1] = 1.0;
	}

	for (i=1; i <= nel-1; i++) {
		sprintf(quest, "Mole fraction of %s:",cmpname[i]);
		dummy = xf[i];
		// call gparrd
		if (xf[i] < 1.0e-6) {
			printf("Fraction must be larger than 1.0E-6\n");
			xf[i] = 1.0e-6;
		}
	}
	// calculate the equilibria
	target = " ";
	n1 = 0;
	n2 = 0;

	c_tqce(target, n1, n2, value, ceq);
	// list some results
	// amount of all phases
	statevar = "NP";
	n1 = -1;
	n2 = 0;
	n3 = sizeof(npf);
	c_tqgetv(statevar, n1,n2,n3, npf, ceq);
	printf("Amount of %i phases: ", n3,);
	mm=0;
	for (i=1; i <= nph; i++) {
		for( ics=1; ics <= noofcs(i); ics++)
		       mm++;
		       if(npf[mm] > 0.0) { 
// the phase is stable if it has a positive amount ... it can be stable with 0
 				phstable=10*i+ics;
		 		if (ics > 1) 
					printf("%s %i %lf", phnames[i], ics, npf[mm]);
				else
                    printf("%s %lf", phnames[i], npf[mm]);
                //composition of stable phase n2=-1 means all fractions
                statevar = 'X';
                n2=-1;
                //Use extended phase index: 10*phase number + compset number
                n4=sizeof(pxf);
                c_tqgetv(statevar, 10*n+ics, n2, n4, pxf, ceq);
                

							
  //c_tqgnp(nph, &ceq);
	return 0;
}

