#include "tqintf.h"

void example_1(string fname, double t, double p, double n, vector<double> x);
void example_2(string fname, int phidx, double t, double p, double n,
                                            vector<double> x, vector<double> y);
void example_3(string fname, int phidx, double t, double p, double n,
                                                              vector<double> y);


int main(int argc, char **argv)
{

    /******************* EXAMPLE 1 *******************

      Similar to /C/cexample1:
    calls the TQ-Interface for any given thermodynamic
    database file {FILENAME}, sets conditions specified
    in this function {T}, {P}, {N} and {X[*]} and prints
    the output of the thermodynamic equilibrium calculation.

    string FILENAME = "FENI.TDB";                                               // Name of the thermodynamic database file (*.TDB, *.tdb)
    int I = 1;                                                                  // Number of Phase
    double T = 778.0;                                                           // Temperature in K
    double P = 1.0e5;                                                           // Pressure in Pa
    double N = 1.0;                                                             // Number of moles
    vector<double> X;                                                           // Concentration array for phase I
        X.push_back(0.6);   // manual override of X[0]
        X.push_back(0.4);   // manual override of X[1]
      //X.push_back(0.3);   //  .. and so on ..   X[2]

    example_1(FILENAME, T, P, N, X);                                            // Call to Example 1

    /**************************************************/




    /******************* EXAMPLE 2 *******************

      Similar to /F90/test4::
    calls the TQ-Interface for any given thermodynamic
    database file {FILENAME}, suspends all phases except
    phase {I}, sets conditions for this single phase as
    specified in this function {T}, {P}, {N} and {X[*]}
    and prints the output of the thermodynamic values
    like the Gibbs energy, the partial derivative of the
    Gibbs energy with respect to every single site fraction
    without doing a thermodyamic equilibrium calculation.*/

    string FILENAME = "steel1.TDB";                                             // Name of the thermodynamic database file (*.TDB, *.tdb)
    int I = 2;                                                                  // Number of Phase
    double T = 8.0e2;                                                           // Temperature in K
    double P = 1.0e5;                                                           // Pressure in Pa
    double N = 1.0;                                                             // Number of moles
    vector<double> X;                                                           // Concentration array for the system
        X.push_back(0.3);        // manual override of X[0]
    vector<double> Y;                                                           // Constituents array for phase I
        Y.push_back(0.197577);   // manual override of Y[0]
        Y.push_back(0.802423);   // manual override of Y[1]
        Y.push_back(1);          // manual override of Y[2]

    example_2(FILENAME, I, T, P, N, X, Y);                                      // Call to Example 2

    /**************************************************/





    /******************* EXAMPLE 3 *******************

    experimental - work in progress

    *************************************************
    
    string FILENAME = "steel1.TDB";                                             // Name of the thermodynamic database file (*.TDB, *.tdb)
    int I = 2;                                                                  // Number of Phase
    double T = 8.0e2;                                                           // Temperature in K
    double P = 1.0e5;                                                           // Pressure in Pa
    double N = 1.0;                                                             // Number of moles
    vector<double> X;                                                           // Concentration array for the system
        X.push_back(0.3);        // manual override of X[0]

    example_3(FILENAME, I, T, P, N, X);                                         // Call to Example 3

    /**************************************************/

    return 0;
}

/********************************** EXAMPLE 1 *********************************/

void example_1(string fname, double t, double p, double n, vector<double> x)
{
    void *ceq = 0;                                                              // Pointer to the OpenCalphad storage
    vector<vector<double> > elfract;                                            // Array including all equilibrium compositions
    vector<string> phnames;                                                     // Array including all phase names
    vector<double> phfract;                                                     // Array including all phase fractions

    //-----------------------Initialize and read TDB data-----------------------

    Initialize(&ceq);                                                           // Initialize OpenCalphad and allocate memory
    ReadDatabase(fname, &ceq);                                                  // Define TDB-file and read elements
    ReadPhases(phnames, &ceq);                                                  // Read Phases data
    SetTemperature(t, &ceq);                                                    // Set Temperature
    SetPressure(p, &ceq);                                                       // Set Pressure
    SetMoles(n, &ceq);                                                          // Set Number of moles
    SetComposition(x, &ceq);                                                    // Set Composition of the system

    //---------------------------Calculate Equilibrium--------------------------
 
    CalculateEquilibrium(&ceq);                                                 // Calculate a phase equilibrium

    //-------------------------------List Results-------------------------------

    ListPhaseFractions(phnames, phfract, &ceq);                                 // Write output of the amount of stable phases
    ListConstituentFractions(phnames, phfract, elfract, &ceq);                  // Write output of the composition of each stable phase

}

/********************************** EXAMPLE 2 *********************************/

void example_2(string fname, int phidx, double t, double p, double n,
                                             vector<double> x, vector<double> y)
{
    void *ceq = 0;                                                              // Pointer to the OpenCalphad storage
    vector<string> elnames;                                                     // Array including selected elements
        elnames.push_back("CR");
        elnames.push_back("FE");
    vector<string> phnames;                                                     // Array including all phase names
    vector<double> phfract;                                                     // Array including all phase fractions
    vector<vector<double> > elfract;                                            // Array including all equilibrium compositions

    //-----------------------Initialize and read TDB data-----------------------

    Initialize(&ceq);                                                           // Initialize OpenCalphad and allocate memory
    ReadDatabaseLimited(fname, elnames, &ceq);                                  // Define TDB-file and read only selected elements
    ReadPhases(phnames, &ceq);                                                  // Read Phases data
    SetTemperature(t, &ceq);                                                    // Set Temperature
    SetPressure(p, &ceq);                                                       // Set Pressure
    SetMoles(n, &ceq);                                                          // Set Number of moles
    SetComposition(x, &ceq);                                                    // Set Composition of the system

    //---------------------------Calculate Equilibrium--------------------------
 
    CalculateEquilibrium(&ceq);   

    //-------------------------------List Results-------------------------------

    ListPhaseFractions(phnames, phfract, &ceq);                                 // Write output of the amount of stable phases
    ListConstituentFractions(phnames, phfract, elfract, &ceq);                  // Write output of the composition of each stable phase
    ListExtConstituentFractions(phidx, phnames, &ceq);                          // Write output of the constituents of a given phase

    //----------------------------Change Parameters-----------------------------

    SetConstituents(phidx, y, &ceq);                                            // Set Constituents of the phase

    //-------------------------------List Results-------------------------------

    GetGibbsData(phidx, &ceq);                                                  // Write output of the thermodynamic values of the given parameters
};

/********************************** EXAMPLE 3 *********************************/

void example_3(string fname, int phidx, double t, double p, double n,
                                                               vector<double> x)
{
    void *ceq = 0;                                                              // Pointer to the OpenCalphad storage
    vector<string> elnames;                                                     // Array including selected elements
        elnames.push_back("CR");
        elnames.push_back("FE");
    vector<string> phnames;                                                     // Array including all phase names
    vector<double> phfract;                                                     // Array including all phase fractions
    vector<vector<double> > elfract;                                            // Array including all equilibrium compositions

    //-----------------------Initialize and read TDB data-----------------------

    Initialize(&ceq);                                                           // Initialize OpenCalphad and allocate memory
    ReadDatabaseLimited(fname, elnames, &ceq);                                  // Define TDB-file and read only selected elements
    ReadPhases(phnames, &ceq);                                                  // Read Phases data
    SetTemperature(t, &ceq);                                                    // Set Temperature
    SetPressure(p, &ceq);                                                       // Set Pressure
    SetMoles(n, &ceq);                                                          // Set Number of moles
    SetComposition(x, &ceq);                                                    // Set Composition of the system

    //---------------------------Calculate Equilibrium--------------------------
 
    CalculateEquilibrium(&ceq);   

    //-------------------------------List Results-------------------------------

    ListPhaseFractions(phnames, phfract, &ceq);                                 // Write output of the amount of stable phases
    ListConstituentFractions(phnames, phfract, elfract, &ceq);                  // Write output of the composition of each stable phase
    ListExtConstituentFractions(phidx, phnames, &ceq);                          // Write output of the constituents of a given phase

    for(int i = 0; i < 10; i++)
    {
        cout << "========== " << i << " / 10 ==========" << endl;
        double constit = i/10.0;
        vector<double> y;                                                       // Constituents array for phase I
            y.push_back(constit);     // manual override of Y[0]
            y.push_back(1-constit);   // manual override of Y[1]
            y.push_back(1);           // manual override of Y[2]

        //--------------------------Change Parameters---------------------------

        SetConstituents(phidx, y, &ceq);                                        // Set Constituents of the phase

        //-----------------------------List Results-----------------------------

        GetGibbsData(phidx, &ceq);                                              // Write output of the thermodynamic values of the given parameters
    }
};
