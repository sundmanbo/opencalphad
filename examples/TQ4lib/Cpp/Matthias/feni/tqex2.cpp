#include "../liboctqcpp.h"
#include <cmath>
using namespace std;

int main(int argc, char *argv[])
{
    /* Bugtracker:
    1) No equilph1d(phtup,ceq%tpval,xknown,mu,.TRUE.,nend,mugrad,mobilities,ceq)
       function in liboctq.F90
    */

    liboctqcpp OCASI;
    void * ceq = 0;
    string filename = "TQ4lib/Cpp/Matthias/feni/FENI.TDB";

    cout << endl;
    cout << "Calculation of equilibria and mobility data in Fe-Ni system" << endl;
    cout << endl;

    cout << "Fictitious ln(mobility data) in the TDB file:" << endl;
    cout << "PARAMETER MQ&FE(LIQUID,FE;0) 298.15 -10000*IQRT-18; 6000 N BOS !" << endl;
    cout << "PARAMETER MQ&FE(LIQUID,NI;0) 298.15 -12000*IQRT-19; 6000 N BOS !" << endl;
    cout << "PARAMETER MQ&NI(LIQUID,NI;0) 298.15 -11000*IQRT-18; 6000 N BOS !" << endl;
    cout << "PARAMETER MQ&NI(LIQUID,FE;0) 298.15 -13000*IQRT-19; 6000 N BOS !" << endl;
    cout << "PARAMETER MQ&FE(FCC_A1,FE:VA;0) 298.15 -30000*IQRT-30; 6000 N BOS !" << endl;
    cout << "PARAMETER MQ&FE(FCC_A1,NI:VA;0) 298.15 -32000*IQRT-32; 6000 N BOS !" << endl;
    cout << "PARAMETER MQ&NI(FCC_A1,NI:VA;0) 298.15 -33000*IQRT-32; 6000 N BOS !" << endl;
    cout << "PARAMETER MQ&NI(FCC_A1,FE:VA;0) 298.15 -25000*IQRT-31; 6000 N BOS !" << endl;
    cout << endl;

    OCASI.tqini(0,&ceq);

    vector<string> Elements;
    Elements.push_back("FE");
    Elements.push_back("NI");
    vector<string> elnames =

    OCASI.tqrpfil("TQ4lib/Cpp/Matthias/feni/FENI.TDB", Elements, &ceq);

    int nel = OCASI.tqgcom(&ceq).size();
    cout << "System with " << nel << " elements: ";
    for(int i = 0; i < nel; i++)
    {
        cout << OCASI.tqgcom(&ceq)[i];
        if(i < nel-1)
        {
            cout << ", ";
        }
    }
    cout << endl;


    int phasetuples = OCASI.tqgnp(&ceq);

    cout << "and " << phasetuples << " phases: ";
    for(int i = 0; i < phasetuples; i++)
    {
        cout << OCASI.tqgpn(i+1, &ceq);
        if(i < phasetuples-1)
        {
            cout << ", ";
        }
    }
    cout << endl;
    cout << endl;
    cout << "Give conditions:" << endl;

    cout << "Temperature: /1000/:" << endl;
    OCASI.tqsetc("T", 0, 0, 1000, &ceq);
    cout << "Pressure: /100000/:" << endl;
    OCASI.tqsetc("P", 0, 0, 100000, &ceq);
    cout << "Mole fractions for FE: /0.5/:" << endl;
    OCASI.tqsetc("X", 1, 0, 0.5, &ceq);
    OCASI.tqsetc("N", 0, 0, 1.0, &ceq);

    OCASI.tqce(&ceq);

    cout << endl;
    cout << "Successfull calculation" << endl;
    cout << endl;


    int equphasetuples = OCASI.tqgnp(&ceq);
    cout << "Amount of " << equphasetuples << " phases: ";
    for(int i = 0; i < equphasetuples; i++)
    {
        cout << OCASI.PhaseFractions(&ceq)[i];
        if(i < equphasetuples-1)
        {
            cout << ", ";
        }
    }
    cout << endl;

    for(unsigned int i = 0; i < equphasetuples; i++)
    if(OCASI.PhaseFractions(&ceq)[i] > 0.0)
    {
        cout << "Stable phase: " << OCASI.tqgpn(i+1, &ceq) << ", amount: "
             << OCASI.PhaseFractions(&ceq)[i] << ", mole fractions:" << endl;

        for(int j = 0; j < nel; j++)
        {
            cout << OCASI.tqgcom(&ceq)[j] << ": "
                 << OCASI.tqgetv("X", i+1, -1, 4, &ceq)[j];
            if(j < nel-1)
            {
                cout << ", ";
            }
        }
        cout << endl;
        cout << endl;
    }

    cout << "System volume: " << OCASI.tqgetv("V", 0, 0, &ceq) << endl;
    cout << endl;

    cout << "Component, mole fraction, chemical potentials, lnac = mu/RT"
         << endl;
    for(int i = 0; i < elnames.size(); i++)
    {
        cout << OCASI.tqgcom(&ceq)[i] << "    "
             << OCASI.tqgetv("X", i+1, 0, &ceq) << "    "
             << OCASI.tqgetv("MU", i+1, 0, &ceq) << "    "
             << OCASI.tqgetv("MU", i+1, 0, &ceq)/8.31451
                /OCASI.tqgetv("T", 0, 0, &ceq) << endl;
    }
    cout << endl;

    cout << "LN(mobility of component in phase) and exp(..):" << endl;
    for(int i = 0; i < equphasetuples; i++)
    for(int j = 0; j < nel; j++)
    {
        string temp = "MQ&" + OCASI.tqgcom(&ceq)[j];
        cout << temp << "(" << OCASI.tqgpn(i+1, &ceq) << ") = "
             << OCASI.tqgetv(temp, i+1, j+1, &ceq) << " "
             << exp(OCASI.tqgetv(temp, i+1, j+1, &ceq)) << endl;
    }

    cout << endl;
    cout << "Calculating Darken stability matrix, dG_A/dN_B for phase  2:" << endl;

    /*Calculating Darken stability matrix, dG_A/dN_B for phase  2:
    Calculation required    6 its

    Chemical potential derivative matrix, dG_I/dN_J for   2 endmembers
                  1           2
      1  1.2100E+04 -1.2100E+04
      2 -1.2100E+04  1.2100E+04

    LN(mobility) values for  2 components
      1 -3.4728E+01 -3.4988E+01*/

    cout << endl;
    cout << "Any more calculations? /N/:" << endl;
    cout << endl;
    cout << "Auf wiedersehen" << endl;

    return 0;
}
