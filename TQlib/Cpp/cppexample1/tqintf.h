#define MAXEL 10
#define MAXPH 10
#include <string>
#include <stdlib.h>
#include <iostream>
#include <cstring>

extern"C"
{
    void c_tqini(int, void *);                                                  // initiates the OC package
    void c_tqrfil(char *, void *);                                              // read all elements from a TDB file
    //void c_tqgcom(int *, char[MAXEL][24], void **);                           // get system component names. At present the elements
    //void c_tqrpfil(char *, int, char **, void **);                            // read TDB file with selection of elements
    //void c_tqgnp(int *, void **);                                             // get total number of phases and composition sets
    void c_tqgpn(int, char *, void *);                                          // get name of phase+compset tuple with index phcsx 
    void c_tqgetv(char *, int, int, int *, double *, void *);                   // get equilibrium results using state variables
    void c_tqsetc(char *, int, int, double, int *, void *);                     // set condition
    void c_tqce(char *, int, int, double *, void *);                            // calculate quilibrium with possible target
    //void c_tqgnp(int, gtp_equilibrium_data **);                               // get total number of phases and composition sets
    void examine_gtp_equilibrium_data(void *);                                  //	
}

extern"C" int  c_ntup;                                                          //
extern"C" int  c_nel;                                                           // number of elements
extern"C" int  c_maxc;                                                          //
extern"C" char *c_cnam[24];                                                     // character array with all element names

std::string fname = "FENI.TDB";                                                 //
void *ceq = 0;                                                                  //

// Internal Parameters

char phnames[MAXPH][24];
int cnum[MAXEL+3];
double pxf[10*MAXPH];
double npf[MAXPH];
gtp_equilibrium_data *dceq;

void Initialize(void *ceq)
{
   int n = 0;

    /*
    std::cout << "tqini(" << n << ", &ceq)" <<
    std::endl;
    */

    //===============
    c_tqini(n, ceq);
    //===============

   std::cout << "-> Adress of ceq-Storage: [" << &ceq << "]" <<
   std::endl;
};

void ReadElements(std::string fname, void *ceq)
{
    char *filename = strcpy((char*)malloc(fname.length()+1), fname.c_str());

    /*
    std::cout << "tqrfil(" << filename << ", &ceq)" <<
    std::endl;
    */

    //======================
    c_tqrfil(filename, ceq);
    //======================

    std::cout << "-> Element Data: [";
    for(int i = 0; i < c_nel; i++)
    {
        std::cout << c_cnam[i];
        if(i < c_nel-1)
        {
            std::cout << ", ";
        }
    }
    std::cout << "] [" << &ceq << "]" <<
    std::endl;
};

void ReadPhases(void *ceq)
{      
    /*
    std::cout << "tqgpn(&ceq)" <<
    std::endl;
    */

    for(int i = 1; i < c_ntup+1; i++)
    {
        //==========================
        c_tqgpn(i, phnames[i], ceq);
        //==========================
    }

    std::cout << "-> Phase Data: [";
    for(int i = 1; i < c_ntup+1; i++)
    {
        std::cout << phnames[i];
        if(i < c_ntup)
        {
            std::cout << ", ";
        }
    }
    std::cout << "] [" << &ceq << "]" <<
    std::endl;
};

void SetTemperature(double T, void *ceq)
{
    int n1 = 0;
    int n2 = 0;
    char par[60] = "T";
    if (T < 1.0) T = 1.0;

    /*
    std::cout << "tqsetc(" << par << ", " << n1 << ", " << n2 << ", " <<
                     T << ", &(cnum[1]), &ceq)" <<
    std::endl;
    */

    //=========================================
    c_tqsetc(par, n1, n2, T, &(cnum[1]), ceq);
    //=========================================

    std::cout << "-> Set Temperature to: [" << T << "] [" << &ceq << "]" <<
    std::endl;
};

void SetPressure(double P, void *ceq)
{
    int n1 = 0;
    int n2 = 0;
    char par[60] = "P";
    if (P < 1.0) P = 1.0;

    std::cout << "tqsetc(" << par << ", " << n1 << ", " << n2 << ", " <<
                 P << ", &(cnum[2]), &ceq)" <<
    std::endl;

    //=========================================
    c_tqsetc(par, n1, n2, P, &(cnum[1]), ceq);
    //=========================================

    std::cout << "-> Set Pressure to: [" << P << "] [" << &ceq << "]" <<
    std::endl;
};

void SetMoles(double N, void *ceq)
{
    int n1 = 0;
    int n2 = 0;
    char par[60] = "N";

    /*
    std::cout << "tqsetc(" << par << ", " << n1 << ", " << n2 << ", " <<
                 N << ", &(cnum[3]), &ceq)" <<
    std::endl;
    */

    //=========================================
    c_tqsetc(par, n1, n2, N, &(cnum[1]), ceq);
    //=========================================

    std::cout << "-> Set Moles to: [" << N << "] [" << &ceq << "]" <<
    std::endl;
};

void SetComposition(double xf[MAXEL], void *ceq)
{
    int n1 = 0;
    int n2 = 0;
    char par[60] = "X";

    for (int i = 0; i < c_nel - 1; i++)
    {
        if (xf[i] < 1.0e-6) xf[i] = 1.0e-6;

        /*
        std::cout << "tqsetc(" << par << ", " << i+1 << ", " << n2 <<
                     ", " << xf[i] << ", &(cnum[" << 4+i << "]), &ceq)" <<
        std::endl;
        */

        //==================================================
        c_tqsetc(par, i+1, n2, xf[i], &(cnum[4 + i]), ceq);
        //==================================================

        std::cout << "-> Set Composition of " << c_cnam[i] << " to: [" <<
                     xf[i] << "] [" << &ceq << "]" <<
        std::endl;
    }
};

void CalculateEquilibrium(void *ceq)
{
    char target[60] = " ";
    int null1 = 0;
    int null2 = 0;
    double val;

    /*
    std::cout << "tqce(" << target << ", " << null1 << ", " << null2 <<
                 ", &val, &ceq)" <<
    std::endl;
    */

    //======================================
    c_tqce(target, null1, null2, &val, ceq);
    //======================================
};

void ListPhaseFractions(void *ceq)
{
    char statevar[60] = "NP";
    int n1 = -1;
    int n2 =  0;
    int n3 = sizeof(npf) / sizeof(npf[0]);

    /*
    std::cout << "tqgetv(" << statevar << ", " << n1 << ", " << n2 <<
                 ", &n3, npf, " << ", ceq)" <<
    std::endl;
    */

    //========================================
    c_tqgetv(statevar, n1, n2, &n3, npf, ceq);
    //========================================

    std::cout << "-> Phase Fractions: [";
    for (int i = 0; i < n3; i++)
    {
        std::cout << phnames[i+1] << ": " << npf[i]; 
        if(i < n3-1)
        {
            std::cout << ", ";
        }
    }
    std::cout << "] [" << &ceq << "]" <<
    std::endl;
};

void ListConstituentFractions(void *ceq)
{
    for (int i = 1; i <= c_ntup; i++)
    {
        if (npf[i] > 0.0)
        {
            char* statevar = "X";
            int n1 =  0;
            int n2 = -1;                                                        //composition of stable phase n2 = -1 means all fractions
            int n4 = sizeof(pxf)/sizeof(pxf[0]);

            /*
            std::cout << "tqgetv(" << statevar << ", " << i << ", " << n2 <<
                         ", &n4, pxf, ceq);" <<
            std::endl;
            */

            //=======================================
            c_tqgetv(statevar, i, n2, &n4, pxf, ceq);
            //=======================================

            std::cout << "-> Constituent Fractions for " << phnames[i] <<
                         " [";

            for (int k = 0; k < n4; k++)
            {
                std::cout << c_cnam[k] << ": " << pxf[k];
                if(k < n4-1)
                {
                    std::cout << ", ";
                }
            }
            std::cout << "] [" << &ceq << "]" <<
            std::endl;
        }
    }
};
