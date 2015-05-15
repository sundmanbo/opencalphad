#include "../../octqc.h"
#include "tqintf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXEL 10
#define MAXPH 10

int main(int argc, char **argv)
{
    // Initialize and read TDB data
    Initialize(&ceq);                                                           //
    ReadElements(fname, &ceq);                                                  //
    ReadPhases(&ceq);                                                           //

    // Set Conditions
    double T = 1.0e3;
    double P = 1.0e5;
    double N = 1.0;
    double xf[c_nel] = {0.4,0.6};

    SetTemperature(T, &ceq);                                                    //
    SetPressure(P, &ceq);                                                       //
    SetMoles(N, &ceq);                                                          //
    SetComposition(xf, &ceq);                                                   //

    // Calculate Equilibrium
    CalculateEquilibrium(&ceq);                                                 //

    // List Equilibrium
    ListPhaseFractions(&ceq);                                                   //
    ListConstituentFractions(&ceq);                                             //

    return 0;
}
