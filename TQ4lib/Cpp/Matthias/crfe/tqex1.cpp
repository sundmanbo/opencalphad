#include "../liboctqcpp.h"
using namespace std;

int main(int argc, char *argv[])
{
    /* Bugtracker:
    1) OCASI.tqgetv("MUS", -1, 0, nel, &ceq)[i] returns error:
         Unknown state variable: MUS                 >:<MUS
     3F Error entering get_many_svar         8888   0.0000000000000000
     3F Failed decode statevar in get_many_svar
    */

    liboctqcpp OCASI;
    void * ceq = 0;
    string filename = "TQ4lib/Cpp/Matthias/crfe/crfe.TDB";

    OCASI.tqini(0,&ceq);

    cout << "Reading all elements from the database file: " << filename << endl;

    OCASI.tqrfil(filename, ceq);

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

    cout << "Temperature: /800/:" << endl;
    OCASI.tqsetc("T", 0, 0, 800, &ceq);
    cout << "Pressure: /100000/:" << endl;
    OCASI.tqsetc("P", 0, 0, 100000, &ceq);
    cout << "Mole fractions for CR: /0.25/:" << endl;
    OCASI.tqsetc("X", 1, 0, 0.25, &ceq);
    OCASI.tqsetc("N", 0, 0, 1.0, &ceq);

    OCASI.tqce(&ceq);

    int equphasetuples = OCASI.tqgnp(&ceq);

    cout << endl;
    cout << "Successfull calculation" << endl;

    cout << "Tuple index  Phase name                 Amount" << endl;
    for(int i = 0; i < equphasetuples; i++)
    {
        cout << i+1 << "    " << OCASI.tqgpn(i+1, &ceq) << "    "
             << OCASI.PhaseFractions(&ceq)[i] << endl;
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

    cout << "Component, mole fraction,  chemical potential (SER)   BCC" << endl;
    for(unsigned int i = 0; i < nel; i++)
    {
        cout << OCASI.tqgcom(&ceq)[i] << "   "
             << OCASI.tqgetv("X", -1, 0, nel, &ceq)[i] << "    "
             << OCASI.tqgetv("MU", -1, 0, nel, &ceq)[i] << "    "
             << OCASI.tqgetv("MU", -1, 0, nel, &ceq)[i] << endl;
    }
    cout << endl;

    cout << "Mole fractions of all components in stable phases:" << endl;
    cout << "X(*,*):";
    for(unsigned int i = 0; i < 4; i++)
    {
        cout << " " << OCASI.tqgetv("X(*,*)", -1, -1, 4, &ceq)[i];
    }
    cout << endl;
    cout << "Mole fraction of a component in all phases, also those unstable:" << endl;
    cout << "in phase tuple order!" << endl;
    cout << "X(*,CR):";
    for(unsigned int i = 0; i < 4; i++)
    {
        cout << " " << OCASI.tqgetv("X(*,*)", -1, 1, 4, &ceq)[i];
    }
    cout << endl;
    OCASI.tqlr(0, &ceq);
    cout << endl;
    cout << "Any more calculations? /N/:" << endl;
    cout << endl;
    cout << "Auf wiedersehen" << endl;
    return 0;
}
