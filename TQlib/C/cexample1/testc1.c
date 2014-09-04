/* sample program calling the OC Fortran library from C */

#include <stdio.h>

extern int liboctqc1(int *select, int *k, char *fn, double *value, int *ceq);

main()
{ 
  int i,j,k,p,ierr,select,nel,nph,ic,loop,ceq;
  char *fn;
  // define an array for 20 component names, each 24 character long ...
  char elem[20][24];
  char phnames[20][24];
  char dummy[24];
  double value,val,mucr,mufe;
  // inititate
  printf("starting\n" );
  ceq=1;
  //------------- initiate OC, all arguments except select dummy
  select=1;
  ierr=liboctqc1_(&select,&k,fn,value,&ceq);
  //  printf("Error code returned for call 1: %i \n",ierr);
  //-------------- read database, fn specifies file name
  // read the TDB file
  select=2;
  fn="crfe.TDB ";
  printf("File to be opened: %s \n",fn);
  ierr=liboctqc1_(&select,&nel,fn,value,&ceq);
  //  printf("Error code returned for call 2: %i \n",ierr);
  //  printf("Equilibrium identifier: %i \n",ceq);
  //-------------- get component nsmaes, k is number, elem is names
  // get the number and names of components
  select=3;
  ierr=liboctqc1_(&select,&nel,&elem,value,&ceq);
  //  printf("Error code returned for call 3: %i \n",ierr);
  printf("Number of components %i \n",nel);
  for (i=0; i<nel; i=i+1)
    {
      printf("Element %i: %s \n",i+1,elem[i]);
    }
  //-------------- get phase nsmaes, nph is number, elem is names
  // get the number and names of components
  select=4;
  ierr=liboctqc1_(&select,&nph,fn,value,&ceq);
  //  printf("Error code returned for call 4: %i \n",ierr);
  printf("Number of phases %i \n",nph);
  select=5;
  for (i=0; i<nph; i=i+1)
    {
      ierr=liboctqc1_(&select,&i,&phnames[i][0],value,&ceq);
      //      printf("Error code returned for call 5: %i \n",ierr);
      printf("Phase name %i %s \n",i+1,phnames[i]);
      //      printf("Phase %i: %s \n",i+1,phnames[i][1]);
    }
  //  printf("Error code returned for call 5: %i \n",ierr);
  //-----------------------------------------
  printf("Setting conditions P=1E5 \n");
  select=6;
  value=1.0E5;
  fn="P";
  ic=0;
  ierr=liboctqc1_(&select,&ic,fn,&value,&ceq);
  //  printf("Error code returned for call 6: %i \n",ierr);
  //-----------------------------------------
  printf("Setting conditions N=1 \n");
  value=1.0;
  fn="N";
  ierr=liboctqc1_(&select,&ic,fn,&value,&ceq);
  //  printf("Error code returned for call 6: %i \n",ierr);
  loop=0;
  while(loop=0);
    {
      // goto here
      //-----------------------------------------
      printf("Condition for T? \n");
      scanf("%lf",&value);
      fn="T";
      ierr=liboctqc1_(&select,&ic,fn,&value,&ceq);
      //      printf("Error code returned for call 6: %i \n",ierr);
      //-----------------------------------------
      printf("Setting conditions x(CR) \n");
      scanf("%lf",&value);
      fn="X";
      ic=1;
      ierr=liboctqc1_(&select,&ic,fn,&value,&ceq);
      //      printf("Error code returned for call 6: %i \n",ierr);
      //-----------------------------------------
      printf("Calculating equilibrium \n");
      select=7;
      fn=" ";
      ic=0;
      ierr=liboctqc1_(&select,&ic,fn,&value,&ceq);
      //      printf("Error code returned for call 7: %i \n",ierr);
      //-----------------------------------------
      printf("\nCalculation successful \n");
      select=8;
      fn="GM";
      ierr=liboctqc1_(&select,&ic,fn,&mucr,&ceq);
      //      printf("Error code returned for call 7: %i \n",ierr);
      printf("   Total Gibbs energy is       %e \n",mucr);
      //-----------------------------------------
      //      printf("Extraction chemical potentials \n");
      select=8;
      fn="MU";
      ic=1;
      ierr=liboctqc1_(&select,&ic,fn,&mucr,&ceq);
      //      printf("Error code returned for call 7: %i \n",ierr);
      printf("   Chemical potential of %s is %e \n",elem[ic-1],mucr);
      ic=2;
      ierr=liboctqc1_(&select,&ic,fn,&mufe,&ceq);
      //      printf("Error code returned for call 7: %i \n",ierr);
      printf("   Chemical potential of %s is %e \n",elem[ic-1],mufe);
      printf("Any more calculations?  0 = YES,  anything else = NO \n");
      scanf("%i",&loop);
      printf("input value %i \n",loop);
    }
  // the while loop does not work??
    printf("Program terminated %i \n",loop);
  // if the user answer YES jump back to "goto here" else quit
  //-----------------------------------------
}
