#include "octqc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXEL 10
#define MAXPH 10

extern void	c_tqini(int, void **);
extern void	c_tqrfil(char *, void **);
extern void	c_tqgcom(int *, char[MAXEL][24], void **);
extern void	c_tqrpfil(char *, int, char **, void **);
extern void	c_tqgnp(int *, void **);
extern void	c_tqgpn(int, char *, void **);
extern void	c_tqgetv(char *, int, int, int *, double *, void **);
extern void	examine_gtp_equilibrium_data(void *);
//extern void	c_tqgnp(int, gtp_equilibrium_data **);

int
main(int argc, char **argv)
{
	int		i         , k, n1, n2, n3, n4, nel = 2, nph, mm, ics, ip, phstable,
			cnum          [MAXEL + 3];
	void           *ceq = 0;
	//char		**cmpname   ;
	char            cmpname   [MAXPH][24];
	char		phnames   [MAXPH][24];
	char		target    [60];
	char		statevar  [60];
	char		quest     [60];
	char		condition [60];
	double		one = 1.0;
	double		value  , tp[2], dummy, mel[MAXEL], mf[MAXEL], volume;
	double		xf      [MAXEL], pxf[10 * MAXPH], npf[MAXPH], mu[MAXEL];
	gtp_equilibrium_data *dceq;
	//set some defaults
	int		n = 0;
	char           *filename = "FENI.TDB";
	char           *selel[2] = {"FE", "NI"};
	
	//initilize tq
	c_tqini(n, &ceq);
	c_tqrfil(filename, &ceq);
	//c_tqrpfil(filename, nel, selel, &ceq);
	//read database file
	// examine_gtp_equilibrium_data(ceq);

	//find out about the system
	// number of elements and their names
	c_tqgcom(&nel, cmpname, &ceq);
	printf("System with %i elements: ", nel);
	for (i = 0; i < nel; i++)
		printf("%s, ",cmpname[i]);
	
	// number of phases and their names
	c_tqgnp(&nph, &ceq);

	printf("and %i phases: ", nph);
	for (i = 1; i <= nph; i++) {
		c_tqgpn(i, phnames[i], &ceq);
		printf(" %s,", phnames[i]);
	}
	//set default values
	tp[0] = 1.0e3;
	tp[1] = 1.0e5;

	for (i = 1; i <= nel; i++)
		xf[i] = 1.0 / (double)nel;

	//ask for conditions
	printf("\nGive conditions: \n");
	dummy = tp[0];

	if (tp[0] < 1.0) {
		printf("Temperature must be larger than 1K\n");
		tp[0] = 1.0;
	}
	dummy = tp[1];
	if (tp[1] < 1.0) {
		printf("Pressure must be larger than 1Pa\n");
		tp[1] = 1.0;
	}
	for (i = 1; i <= nel - 1; i++) {
		sprintf(quest, "Mole fraction of %s:", cmpname[i]);
		dummy = xf[i];
		//call gparrd
		if (xf[i] < 1.0e-6) {
			printf("Fraction must be larger than 1.0E-6\n");
			xf[i] = 1.0e-6;
		}
	}
	//set conditions
	n1 = 0;
	n2 = 0;
	strcpy(condition, "T");
	c_tqsetc(condition, n1, n2, tp[0], &(cnum[1]), &ceq);
	strcpy(condition, "P");
	c_tqsetc(condition, n1, n2, tp[1], &(cnum[2]), &ceq);
	strcpy(condition, "N");
	c_tqsetc(condition, n1, n2, one, &(cnum[3]), &ceq);
	for (i = 1; i <= nel - 1; i++) {
		strcpy(condition, "X");
		c_tqsetc(condition, i, n2, xf[i], &(cnum[3 + i]), &ceq);
	}
	//calculate the equilibria
	strcpy(target, " ");
	n1 = 0;
	n2 = 0;

	c_tqce(&target, n1, n2, &value, &ceq);
	//list some results
	// amount of all phases
	strcpy(statevar, "NP");
	n1 = -1;
	n2 = 0;
	n3 = sizeof(npf) / sizeof(npf[0]);
	c_tqgetv(statevar, n1, n2, &n3, npf, &ceq);
	printf("Amount of %i phases: ", n3);
	for (i = 0; i < n3; i++)
		printf("%lf ", npf[i]);

	mm = 0;
	for (i = 1; i <= nph; i++) {
		for (ics = 1; ics <= c_noofcs(i); ics++) {
			if (npf[mm] > 0.0) {
				//the phase is stable if it
				//has a positive amount...it can be stable with 0
				phstable = 10 * i + ics;
				if (ics > 1)
					printf("\nStable phase: %s#%i, amount: %lf", phnames[i], ics, npf[mm]);
				else
					printf("\nStable phase: %s, amount: %lf", phnames[i], npf[mm]);
				//composition of stable phase n2 = -1 means all fractions
				strcpy(statevar, "X");
				n2 = -1;
				//Use extended phase index:10 * phase number + compset number
				n4 = sizeof(pxf)/sizeof(pxf[0]);
				c_tqgetv(statevar, 10 * i + ics, n2, &n4, pxf, &ceq);
				printf(" mole fractions:\n"); 
				for (k = 0; k < n4; k++)
					printf(" %s : %lf , ", cmpname[k], pxf[k]);
			}
			mm++;
		}
	}

	//c_tqgnp(nph, &ceq);
	return 0;
}
