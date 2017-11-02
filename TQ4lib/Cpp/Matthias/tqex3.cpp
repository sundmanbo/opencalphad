#include "liboctqcpp.h"
using namespace std;

/*******************************************************************************
 * The following routine is for testing the interface functionality and produces
 * debug output. However it also demonstrates a usecase and can be used as a
 * starting point for a new implementation of OpenCalphad C++ interface.
*******************************************************************************/

int main(int argc, char *argv[])
{
    /* Bugtracker:
    1) OCASI.tqgetv("N", 0, 0, &ceq) returns 0 even though it was set to 1.0
    2) OCASI.tqgetv("X", 1, 0, &ceq) returns 0 even though it was set to 0.3
    3) OCASI.tqgpci(4, "CR", &ceq) breaks, when phase 4 is BCC_A2_AUTO#2
    4) OCASI.tqgccf only returns "tqgccf not implemented yet"
    5) OCASI.tqgnpc only returns "tqgnpc not implemented yet"
    6) OCASI.tqgpci only returns "tqgpci not implemented yet"
    7) OCASI.tqgpcs only returns "tqgpcs not implemented yet"
    8) OCASI.reset_conditions only resets a single condition, even though its
       name implies different. Additionally it uses the console!!! to ask for a
       new value for the condition.
    9) OCASI.tqcph2 breaks with an error
   10) OCASI.tqcceq("test", &newceq, &ceq) seems to copy ceq, but not all
       information! OCASI.tqgpn(1, &newceq) seems to have no information, and
       OCASI.tqlc(0,&newceq) crashes.
   11) tqdceq, tqcceq and tqselceq require a CEQ-name! Where is it set?
    */

    liboctqcpp OCASI;
    void * ceq = 0;

    //=================
    OCASI.tqini(0,&ceq);
    //=================

    cout << "-> Adress of ceq-Storage: [" << &ceq << "]" << endl;

/*  uncomment this for reading full database with all elemens

    //=====================================================
    vector<string> elnames
    = OCASI.tqrfil("TQ4lib/Cpp/Matthias/FECRMNC.TDB", ceq);
    //=====================================================

    cout << "-> Element Data: [";
    for(int i = 0; i < elnames.size(); i++)
    {
        cout << elnames[i];
        if(i < elnames.size()-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << " [" << &ceq << "]" << endl;*/

    vector<string> Elements;
    Elements.push_back("CR");
    Elements.push_back("FE");
    vector<string> elnames2 =

    //===============================================================
    OCASI.tqrpfil("TQ4lib/Cpp/Matthias/FECRMNC.TDB", Elements, &ceq);
    //===============================================================

    cout << "-> Element Data: [";
    for(int i = 0; i < elnames2.size(); i++)
    {
        cout << elnames2[i];
        if(i < elnames2.size()-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << " [" << &ceq << "]" << endl;
    vector<string> elnames3 =

    //=================
    OCASI.tqgcom(&ceq);
    //=================

    cout << "-> Element Data: [";
    for(int i = 0; i < elnames3.size(); i++)
    {
        cout << elnames3[i];
        if(i < elnames3.size()-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << " [" << &ceq << "]" << endl;
    int phasetuples =

    //================
    OCASI.tqgnp(&ceq);
    //================

    cout << "-> Number of phasetuples: [" << phasetuples << "]" << endl;
    vector<string> PhNames(phasetuples);
    cout << "-> Phase Data: [";
    for(int i = 0; i < phasetuples; i++)
    {
        PhNames[i] =

        //=====================
        OCASI.tqgpn(i+1, &ceq);
        //=====================

        cout << PhNames[i] << "[" <<

        //==================
        OCASI.tqgpi(PhNames[i], &ceq)
        //==================

        << "]";
        if(i < phasetuples-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << " [" << &ceq << "]" <<
    endl;

    double T = 800;
    double P = 100000;
    double N = 1.0;
    double XCR = 0.3;

    //===================================
    OCASI.tqsetc("T", 0, 0, T, &ceq);
    OCASI.tqsetc("P", 0, 0, P, &ceq);
    OCASI.tqsetc("N", 0, 0, N, &ceq);
    OCASI.tqsetc("X", 1, 0, XCR, &ceq);
    //===================================

    cout << "-> Set Temperature to: [" << T << "]" << " [" << &ceq << "]" << endl;
    cout << "-> Set Ambient Pressure to: [" << P << "]" << " [" << &ceq << "]" << endl;
    cout << "-> Set Moles to: [" << N << "]" << " [" << &ceq << "]" << endl;
    cout << "-> Set X(1) to: [" << XCR << "]" << " [" << &ceq << "]" << endl;

    //=========================================
    T = OCASI.tqgetv("T", 0, 0, &ceq);
    P = OCASI.tqgetv("P", 0, 0, &ceq);
    N = OCASI.tqgetv("N", 0, 0, &ceq);
    XCR = OCASI.tqgetv("X", 1, 0, &ceq);
    //=========================================

    cout << "-> Temperature set to: [" << T << "]" << " [" << &ceq << "]" << endl;
    cout << "-> Ambient Pressure set to: [" << P << "]" << " [" << &ceq << "]" << endl;
    cout << "-> Moles set to: [" << N << "]" << " [" << &ceq << "]" << endl;
    cout << "-> X(1) set to: [" << XCR << "]" << " [" << &ceq << "]" << endl;

    //===============
    OCASI.tqce(&ceq);
    //===============

    cout << "-> Calculated Equilibrium [" << ceq << "]" << endl;
    vector<double> EquPhFr =

    //===========================
    OCASI.PhaseFractions(&ceq);
    //===========================

    phasetuples =

    //================
    OCASI.tqgnp(&ceq);
    //================

    cout << "-> Number of phasetuples: [" << phasetuples << "]" << endl;
    PhNames.resize(phasetuples);
    cout << "-> Phase Data: [";
    for(int i = 0; i < phasetuples; i++)
    {
        PhNames[i] =

        //=====================
        OCASI.tqgpn(i+1, &ceq);
        //=====================

        cout << PhNames[i] << "[" <<

        //===========================
        OCASI.tqgpi(PhNames[i], &ceq)
        //===========================

        << "]";
        if(i < phasetuples-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << " [" << &ceq << "]" <<
    endl;

    cout << "-> Phase Fractions: [";
    for (unsigned int i = 0; i < EquPhFr.size(); i++)
    {
        cout << PhNames[i] << ": " << EquPhFr[i];
        if(i < EquPhFr.size()-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << " [" << &ceq << "]" <<
    endl;

    for(unsigned int phase = 0; phase < EquPhFr.size(); phase++)
    if(EquPhFr[phase] > 0.0)
    {
        cout << "-> Constituent Fractions for " << PhNames[phase] << " [";
        vector<double> ConFr =

        //==========================================
        OCASI.ConstituentFractions(phase+1, &ceq);
        //==========================================

        for(unsigned int i = 0; i < ConFr.size(); i++)
        {
            cout << elnames3[i] << ": " << ConFr[i];
            if(i < elnames3.size()-1)
            {
                cout << ", ";
            }
        }
        cout << "]" << " [" << &ceq << "]" <<
        endl;
    }

    vector<int> ncons;
    vector<int> sites;
    double moles;

    for(int k = 0; k < phasetuples; k++)
    {
    if(k == phasetuples-1)
    {
        cout << "TQGPCN CANNOT BE CALLED FOR PHASE BCC_A2_AUTO#2, BECAUSE IT WILL BREAK!" << endl;
        break;
    }
    vector<double> y = OCASI.tqgphc1(k+1, ncons, sites, moles, &ceq);

    cout << "-> Extended Constituent Fractions for " << PhNames[k]
         << " [" << moles << " moles of atoms/formula unit]";
    int consti = 0;
    for(unsigned int i = 0; i < ncons.size(); i++)
    {
        cout << " [";
        for(int j = 0; j < ncons[i]; j++)
        {
            string cname =
            //==============================
            OCASI.tqgpcn(k+1, consti+1, &ceq);
            //==============================

            cout << cname << "[" <<

            //============================
            //OCASI.tqgpci(k+1, cname, &ceq) // TODO: not yet implemented
            j+1
            //============================

            << "]: " << y[consti];
            if(j < ncons[i]-1)
            {
                cout << ", ";
            }
            consti += 1;
        }
        cout << "]_(" << sites[i] << ")";
    }
    cout << endl;
    }

    cout << "-> For Phase " << PhNames[1] << ":" << endl;

    vector<double> Y2;
    Y2.push_back(0.197577);
    Y2.push_back(0.802423);
    Y2.push_back(1);

    //=========================
    OCASI.tqsphc1(2, Y2, &ceq);
    //=========================

    cout << "-> Set Constituents to: [";
    for(int i = 0; i < Y2.size(); i++)
    {
        cout << i << ": " << Y2[i];
        if(i < Y2.size()-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << endl;

    vector<double> G_TP;
    vector<double> G_Y;
    vector<double> G_YT;
    vector<double> G_YP;
    vector<double> G_YY;

    double G =
    //=================================================
    OCASI.tqcph1(2, G_TP, G_Y, G_YT, G_YP, G_YY, &ceq);
    //=================================================

    cout << "-> Read Gibbs Energy G: [" << G << "]" << endl;
    cout << "-> Read Gibbs Data G: [";
    for(int i = 0; i < 5; i++)
    {
        cout << G_TP[i];
        if(i < 4)
        {
            cout << ", ";
        }
    }
    cout << "]" << endl;

    cout << "-> Read Gibbs Data dGdY: [";
    for(unsigned int i = 0; i < G_Y.size(); i++)
    {
        cout << G_Y[i];
        if(i < G_Y.size()-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << endl;

    cout << "-> Read Gibbs Data d2GdYdT: [";
    for(unsigned int i = 0; i < G_YT.size(); i++)
    {
        cout << G_YT[i];
        if(i < G_YT.size()-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << endl;

    cout << "-> Read Gibbs Data d2GdYdP: [";
    for(unsigned int i = 0; i < G_YP.size(); i++)
    {
        cout << G_YP[i];
        if(i < G_YP.size()-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << endl;

    int kk=G_Y.size()*(G_Y.size()+1)/2;

    cout << "-> Read Gibbs Data d2GdY2: [";
    for(int i = 0; i < kk; i++)
    {
        cout << G_YY[i];
        if(i < kk-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << endl;

    OCASI.tqgccf(1, &ceq);
    OCASI.tqgnpc(1, &ceq);
    OCASI.tqgpci(1, "CR", &ceq);
    OCASI.tqphtupsts(1, -2, 0.0, &ceq);
    double mass;
    OCASI.tqgpcs(1, 1, mass, & ceq);
    //OCASI.tqcph2(1, 0, &ceq);
    OCASI.tqlr(0, &ceq);
    OCASI.tqlc(0, &ceq);
    OCASI.reset_conditions("T", 1000.0, &ceq);
    OCASI.tqlc(0, &ceq);



    int * newceq = 0;
    //==================================
    OCASI.tqcceq("test", &newceq, &ceq);
    //==================================

    cout << "-> Copy CEQ@" << &ceq << " to NEWCEQ@" << &newceq << endl;
    vector<string> elnames4 =

    //=================
    OCASI.tqgcom(&newceq);
    //=================

    cout << "-> Element Data: [";
    for(int i = 0; i < elnames4.size(); i++)
    {
        cout << elnames4[i];
        if(i < elnames4.size()-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << " [" << &newceq << "]" << endl;
    phasetuples =

    //================
    OCASI.tqgnp(&newceq);
    //================

    cout << "-> Number of phasetuples: [" << phasetuples << "]" << endl;
    PhNames.resize(phasetuples);
    cout << "-> Phase Data: [";
    for(int i = 0; i < phasetuples; i++)
    {
        PhNames[i] =

        //=====================
        OCASI.tqgpn(i+1, &newceq);
        //=====================

        cout << PhNames[i] << "[" <<

        //===========================
        OCASI.tqgpi(PhNames[i], &newceq)
        //===========================

        << "]";
        if(i < phasetuples-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << " [" << &newceq << "]" <<
    endl;

    return 0;
}
