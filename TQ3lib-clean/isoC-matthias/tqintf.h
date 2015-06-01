#define MAXEL 10
#define MAXPH 20
#include "octqc.h"
#include <string>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include <vector>

extern"C"
{
    void c_tqini(int, void *);                                                  // initiates the OC package
    void c_tqrfil(char *, void *);                                              // read all elements from a TDB file
    //void c_tqgcom(int *, char[MAXEL][24], void **);                           // get system component names. At present the elements
    void c_tqrpfil(char *, int, char **, void *);                               // read TDB file with selection of elements
    //void c_tqgnp(int *, void **);                                             // get total number of phases and composition sets
    void c_tqgpn(int, char *, void *);                                          // get name of phase+compset tuple with index phcsx 
    void c_tqgetv(char *, int, int, int *, double *, void *);                   // get equilibrium results using state variables
    void c_tqsetc(char *, int, int, double, int *, void *);                     // set condition
    void c_tqce(char *, int, int, double *, void *);                            // calculate quilibrium with possible target
    //void c_tqgnp(int, gtp_equilibrium_data **);                               // get total number of phases and composition sets
    void examine_gtp_equilibrium_data(void *);                                  //
    //void c_getG(int, void *);
    //void c_calcg(int, int, int, int, void *);
    void c_tqgphc1(int, int * , int *, int *, double *, double *, double *,
                                                                        void *);
    void c_tqsphc1(int, double *, double *, void *);
    void c_tqcph1(int, int, int *, double *, double *, double *, double *, double *, void *);
}

extern"C" int  c_ntup;                                                          //
extern"C" int  c_nel;                                                           // number of elements
extern"C" int  c_maxc;                                                          //
extern"C" char *c_cnam[24];                                                     // character array with all element names
extern"C" double c_gval[24];
extern"C" int c_noofcs(int);

using namespace std;

void Initialize(void *ceq)
{
   int n = 0;

    //===============
    c_tqini(n, ceq);
    //===============

   cout << "-> Adress of ceq-Storage: [" << &ceq << "]" <<
   endl;
};

void ReadDatabase(string fname, void *ceq)
{
    char *filename = strcpy((char*)malloc(fname.length()+1), fname.c_str());

    //======================
    c_tqrfil(filename, ceq);
    //======================

    cout << "-> Element Data: [";
    for(int i = 0; i < c_nel; i++)
    {
        cout << c_cnam[i];
        if(i < c_nel-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << " [" << &ceq << "]" <<
    endl;
};

void ReadDatabaseLimited(string fname, vector<string> elnames, void *ceq)
{
    char *filename = strcpy((char*)malloc(fname.length()+1), fname.c_str());
    char *selel[elnames.size()];
    for(int i = 0; i < elnames.size(); i++)
    {
        char *tempchar
             = strcpy((char*)malloc(elnames[i].length()+1), elnames[i].c_str());       
        selel[i] = tempchar;   
    }

    //==============================================
    c_tqrpfil(filename, elnames.size(), selel, ceq);
    //==============================================

    cout << "-> Element Data: [";
    for(int i = 0; i < c_nel; i++)
    {
        cout << c_cnam[i];
        if(i < c_nel-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << " [" << &ceq << "]" <<
    endl;


};

void ReadPhases(vector<string> &phnames, void *ceq)
{
    phnames.clear();

    for(int i = 1; i < c_ntup+1; i++)
    {
        char phn[24];

        //==========================
        c_tqgpn(i, phn, ceq);
        //==========================

        phnames.push_back(phn);
    }

    cout << "-> Phase Data: [";
    for(int i = 0; i < phnames.size(); i++)
    {
        cout << phnames[i];
        if(i < phnames.size()-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << " [" << &ceq << "]" <<
    endl;
};

void SetTemperature(double T, void *ceq)
{
    int cnum;
    int n1 = 0;
    int n2 = 0;
    char par[60] = "T";
    if (T < 1.0) T = 1.0;

    //=========================================
    c_tqsetc(par, n1, n2, T, &cnum, ceq);
    //=========================================

    cout << "-> Set Temperature to: [" << T << "]" << " [" << &ceq << "]" <<
    endl;
};

void SetPressure(double P, void *ceq)
{
    int cnum;
    int n1 = 0;
    int n2 = 0;
    char par[60] = "P";
    if (P < 1.0) P = 1.0;

    //=========================================
    c_tqsetc(par, n1, n2, P, &cnum, ceq);
    //=========================================

    cout << "-> Set Pressure to: [" << P << "]" << " [" << &ceq << "]" <<
    endl;
};

void SetMoles(double N, void *ceq)
{
    int cnum;
    int n1 = 0;
    int n2 = 0;
    char par[60] = "N";

    //=========================================
    c_tqsetc(par, n1, n2, N, &cnum, ceq);
    //=========================================

    cout << "-> Set Moles to: [" << N << "]" << " [" << &ceq << "]" <<
    endl;
};

void SetComposition(vector<double> X, void *ceq)
{
    int cnum;
    int n1 = 0;
    int n2 = 0;
    char par[60] = "X";

    for (int i = 0; i < c_nel; i++)
    {
        if (X[i] < 1.0e-6) X[i] = 1.0e-6;                                       // Check and fix, if composition is below treshold

        if(i < c_nel - 1)
        {                                                                       // Set and print composition, if element 'i' is not the reference/(last) element
            //==================================================
            c_tqsetc(par, i+1, n2, X[i], &cnum, ceq);
            //==================================================

            cout << "-> Set Composition of " << c_cnam[i] << " to: [" <<
                         X[i] << "]" << " [" << &ceq << "]" <<
            endl;
        }
        else
        {                                                                       // Print composition, if element 'i' is the reference/(last) element
           double X_ref = 1;
            for(int j = 0; j < i; j++)
            {
                X_ref -= X[j];
            }

            cout << "-> Set Composition of " << c_cnam[i] << " to: [" <<
                         X_ref << "]" << " [" << &ceq << "]" <<
            endl;
        }
    }
};

void SetConstituents(int phidx, vector<double> y, void *ceq)
{
    int stable1 = phidx;
    double extra[MAXPH];
    double yfr[y.size()];
    for(int i = 0; i < y.size(); i++)
    {
        yfr[i] = y[i];
    }

    //===============================
    c_tqsphc1(stable1,yfr,extra,ceq);
    //===============================

    cout << "-> Set Constituents to: [";
    for(int i = 0; i < y.size(); i++)
    {
        cout << i << ": " << yfr[i];
        if(i < y.size()-1)
        {
            cout << ", ";        
        }
    }
    cout << "]" << endl;
};


void SelectSinglePhase(int PhIdx, void *ceq)
{
    //
};

void CalculateEquilibrium(void *ceq)
{
    char target[60] = " ";
    int null1 = 0;
    int null2 = 0;
    double val;

    //======================================
    c_tqce(target, null1, null2, &val, ceq);
    //======================================

    cout << "-> Calculated Equilibrium [" << ceq << "]"
         << endl;
};

void GetGibbsData(int phidx, void *ceq)
{
    int n2 = 2;
    int n3;
    double gtp[6];
    double dgdy[100];
    double d2gdydt[100];
    double d2gdydp[100];
    double d2gdy2[100];

    //=================================================================
    c_tqcph1(phidx, n2, &n3, gtp, dgdy, d2gdydt, d2gdydp, d2gdy2, ceq);
    //=================================================================

    cout << "-> Read Gibbs Data G: [";
    for(int i = 0; i < 6; i++)
    {
        cout << gtp[i];
        if(i < 5)
        {
            cout << ", ";
        }
    }
    cout << "]" << endl;

    cout << "-> Read Gibbs Data dGdY: [";
    for(int i = 0; i < n3; i++)
    {
        cout << dgdy[i];
        if(i < n3-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << endl;

    cout << "-> Read Gibbs Data d2GdYdT: [";
    for(int i = 0; i < n3; i++)
    {
        cout << d2gdydt[i];
        if(i < n3-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << endl;

    cout << "-> Read Gibbs Data d2GdYdP: [";
    for(int i = 0; i < n3; i++)
    {
        cout << d2gdydp[i];
        if(i < n3-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << endl;

    int kk=n2*(n2+1)/2;

    cout << "-> Read Gibbs Data d2GdY2: [";
    for(int i = 0; i < kk; i++)
    {
        cout << d2gdy2[i];
        if(i < kk-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << endl;
};

void ListPhaseFractions(vector<string> phnames, vector<double>& phfract,
                                                                      void *ceq)
{
    double npf[MAXPH];
    char statevar[60] = "NP";
    int n1 = -1;
    int n2 =  0;
    int n3 = MAXPH;//sizeof(npf) / sizeof(npf[0]);

    //========================================
    c_tqgetv(statevar, n1, n2, &n3, npf, ceq);
    //========================================

    for(int i = 0; i < n3; i++)
    phfract.push_back(npf[i]);

    cout << "-> Phase Fractions: [";
    for (int i = 0; i < n3; i++)
    {
        cout << phnames[i] << ": " << phfract[i]; 
        if(i < n3-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << " [" << &ceq << "]" <<
    endl;
};

void ListConstituentFractions(vector<string> phnames, vector<double> phfract,
                                     vector<vector<double> > elfract, void *ceq)
{
    elfract.clear();
    elfract.resize(phnames.size());
    double pxf[10*MAXPH];
    for (int i = 1; i < c_ntup+1; i++)
    {
        if (phfract[i-1] > 0.0)
        {
            char* statevar = "X";
            int n1 =  0;
            int n2 = -1;                                                        //composition of stable phase n2 = -1 means all fractions
            int n4 = sizeof(pxf)/sizeof(pxf[0]);

            //=======================================
            c_tqgetv(statevar, i, n2, &n4, pxf, ceq);
            //=======================================

            for (int k = 0; k < n4; k++)
            {
                elfract[i-1].push_back(pxf[k]);
            }
            cout << "-> Constituent Fractions for " << phnames[i-1] <<
                         " [";

            for (int k = 0; k < n4; k++)
            {
                cout << c_cnam[k] << ": " << elfract[i-1][k];
                if(k < n4-1)
                {
                    cout << ", ";
                }
            }
            cout << "]" << " [" << &ceq << "]" <<
            endl;
        }
    }
};

void ListExtConstituentFractions(int phidx, vector<string> phnames, void *ceq)
{
    int stable1 = phidx;
    int nlat;
    int nlatc[MAXPH];
    int conlista[MAXPH];
    double yfr[MAXPH];
    double sites[MAXPH];
    double extra[MAXPH];

    //======================================================================
    c_tqgphc1(stable1, &nlat, nlatc, conlista, yfr, sites, extra, ceq);
    //======================================================================

    cout << "-> Extended Constituent Fractions for " << phnames[stable1-1]
         << " [" << extra[0] << " moles of atoms/formula unit]";
    int consti = 0;
    for(int i = 0; i < nlat; i++)
    {
        cout << " [";
        for(int j = 0; j < nlatc[i]; j++)
        {
            cout << "Const. " << consti << ": " << yfr[consti];
            if(j < nlatc[i]-1)
            {
                cout << ", ";
            }
            consti += 1;
        }
        cout << "]_(" << sites[i] << ")";
    }
    cout << endl;
};
