#include "../liboctqcpp.h"
using namespace std;

int main(int argc, char *argv[])
{
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

    cout << endl;
    cout << "Successfull calculation" << endl;
    cout << endl;

    return 0;
}
