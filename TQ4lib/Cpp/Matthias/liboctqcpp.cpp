#include "liboctqcpp.h"
using namespace std;

#define MAXEL 24;
#define MAXPH 20

void liboctqcpp::tqini(int n, void * ceq)
{
    //==============
    c_tqini(n, ceq);
    //==============
};

std::vector<std::string> liboctqcpp::tqrfil(std::string fname, void * ceq)
{
    char *filename = strcpy((char*)malloc(fname.length()+1),fname.c_str());
    //======================
    c_tqrfil(filename, ceq);
    //======================

    ntup = c_ntup;
    nel = c_nel;
    cnam.resize(nel);
    for(int i = 0; i < nel; i++)
    cnam[i] = c_cnam[i];

    //=================
    return tqgcom(ceq);
    //=================
};

std::vector<std::string> liboctqcpp::tqrpfil(std::string fname, std::vector<std::string> elnames, void * ceq)
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

    ntup = c_ntup;
    nel = c_nel;
    cnam.resize(nel);
    for(int i = 0; i < nel; i++)
    cnam[i] = c_cnam[i];

    //=================
    return tqgcom(ceq);
    //=================
};

int liboctqcpp::tqgcn(void * ceq)
{
    int n = MAXEL;
    char elnames[24];
    //=========================
    c_tqgcom(&n, elnames, ceq);
    //=========================
    return n;
};

std::vector<std::string> liboctqcpp::tqgcom(void * ceq)
{
    int n = MAXEL;
    char elnames[24];
    std::vector<std::string> result;
    //=========================
    c_tqgcom(&n, elnames, ceq);
    //=========================
    result.resize(n);
    for(int i = 0; i < n; i++)
    {
        char temp[3];
        for(int j = 0; j < 2; j++)
        {
            temp[j] = elnames[j+i*2];
            if(temp[j] == ' ')
            temp[j] = 0;
        }
        temp[2] = 0;
        string temp2(temp);
        result[i] = temp2;
    }
    return result;
};

int liboctqcpp::tqgnp(void * ceq)
{
    int n;
    //===============
    c_tqgnp(&n, ceq);
    //===============
    return n;
};

std::string liboctqcpp::tqgpn(int i, void * ceq)
{
    char phname[24];
    std::string result;
    //======================
    c_tqgpn(i, phname, ceq);
    //======================
    result = phname;
    return result;
};

int liboctqcpp::tqgpi(std::string pname, void * ceq)
{
    char *phasename = strcpy((char*)malloc(pname.length()+1), pname.c_str());
    int i;
    //=========================
    c_tqgpi(&i, phasename, ceq);
    //=========================
    return i;
};

std::string liboctqcpp::tqgpcn2(int phidx, int i, void * ceq)
{
    //---------------------------
    return tqgpcn(phidx, i, ceq);
    //---------------------------
};

std::string liboctqcpp::tqgpcn(int phidx, int i, void * ceq)
{
    char constituentname[24];
    std::string result;
    //=======================================
    c_tqgpcn(phidx, i, constituentname, ceq);
    //=======================================
    result = constituentname;
    return result;
};

int liboctqcpp::tqgpci(int phidx, std::string cname, void * ceq)
{
    int c;
    char *constituent = strcpy((char*)malloc(cname.length()+1), cname.c_str());
    //====================================
    c_tqgpci(phidx, &c, constituent, ceq); //TODO: c_tqgpci is not implemented in liboctq.F90!!
    //====================================
    return c;
};

void liboctqcpp::tqgpcs(int, int, double *, double *, void *)
{

};

void liboctqcpp::tqgccf(int, int *, char *, double *, double *, void *)
{

};

void liboctqcpp::tqgnpc(int, int *, void *)
{

};

void liboctqcpp::tqphtupsts(int, int, double, void *)
{

};

void liboctqcpp::tqsetc(string par, int n1, int n2, double val, void * ceq)
{
    int cnum;
    char *name = strcpy((char*)malloc(par.length()+1), par.c_str());
    //======================================
    c_tqsetc(name, n1, n2, val, &cnum, ceq);
    //======================================
};

void liboctqcpp::tqce(void * ceq)
{
    char target[60] = " ";
    double val;
    //==============================
    c_tqce(target, 0, 0, &val, ceq);
    //==============================
};

double liboctqcpp::tqgetv(string par, int n1, int n2, void * ceq)
{
    char *name = strcpy((char*)malloc(par.length()+1), par.c_str());
    double val;
    int cnum;
    //=======================================
    c_tqgetv(name, n1, n2, &cnum, &val, ceq);
    //=======================================
    return val;
};

std::vector<double> liboctqcpp::tqgetv(string par, int n1, int n2, int n3, void * ceq)
{
    std::vector<double> results(n3);
    char *name = strcpy((char*)malloc(par.length()+1), par.c_str());
    double val[n3];
    //====================================
    c_tqgetv(name, n1, n2, &n3, val, ceq);
    //====================================
    for(int i = 0; i < n3; i++)
    results[i] = val[i];
    return results;
};

std::vector<double> liboctqcpp::tqgphc1(int phIdx, std::vector<int>& ncons,
                                        std::vector<int>& sites, double& moles,
                                        void * ceq)
{
    int nlat;
    int nlatc[MAXPH];//TODO: MAXPH is misleading
    int conlista[MAXPH];
    double yfr[MAXPH];
    double site[MAXPH];
    double extra[MAXPH];
    //==============================================================
    c_tqgphc1(phIdx, &nlat, nlatc, conlista, yfr, site, extra, ceq);
    //==============================================================
    ncons.resize(nlat);
    sites.resize(nlat);
    int nc = 0;
    for(unsigned int i = 0; i < ncons.size(); i++)
    {
        ncons[i] = nlatc[i];
        sites[i] = site[i];
        nc += nlatc[i];
    }
    vector<double> y(nc, 0);
    for(unsigned int i = 0; i < nc; i++)
    y[i] = yfr[i];
    moles = extra[0];
    return y;
};

void liboctqcpp::tqsphc1(int phidx, std::vector<double> y, void * ceq)
{
    double extra[MAXPH];
    double yfr[y.size()];
    for(int i = 0; i < y.size(); i++)
    yfr[i] = y[i];
    //================================
    c_tqsphc1(phidx, yfr, extra, ceq);
    //================================
};

double liboctqcpp::tqcph1(int phidx, std::vector<double>& G_TP,
                          std::vector<double>& G_Y, std::vector<double>& G_YT,
                          std::vector<double>& G_YP, std::vector<double>& G_YY,
                          void * ceq)
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
    double G = gtp[0];
    G_TP.resize(5);
    G_Y.resize(n3);
    G_YT.resize(n3);
    G_YP.resize(n3);
    int yy = n3*(n3+1.)/2.;
    G_YY.resize(yy);
    for(int i = 0; i < 5; i++)
    G_TP[i] = gtp[i+1];
    //GibbsEnergy_TP[0] G.T
    //GibbsEnergy_TP[1] G.P
    //GibbsEnergy_TP[2] G.T.T
    //GibbsEnergy_TP[3] G.T.P
    //GibbsEnergy_TP[4] G.P.P

    for(int i = 0; i < n3; i++)
    {
        G_Y[i] = dgdy[i];
        //GibbsEnergy_Y[0] G.Y0
        //GibbsEnergy_Y[1] G.Y1

        G_YT[i] = d2gdydt[i];
        //GibbsEnergy_YT[0] G.Y0.T
        //GibbsEnergy_YT[1] G.Y1.T

        G_YP[i] = d2gdydp[i];
        //GibbsEnergy_YP[0] G.Y0.P
        //GibbsEnergy_YP[1] G.Y1.P
    }

    for(int i = 0; i < yy; i++)
    G_YY[i] = d2gdy2[i];
    //GibbsEnergy_YY[0] G.Y0.Y0
    //GibbsEnergy_YY[0] G.Y0.Y1
    //GibbsEnergy_YY[0] G.Y1.Y1
    //GibbsEnergy_YY[0] G.Y1.Y2
    //GibbsEnergy_YY[0] G.Y2.Y2
    //GibbsEnergy_YY[0] G.Y2.Y3

    return G;
};

void liboctqcpp::tqcph2(int, int, int *, int *, void *)
{

};

void liboctqcpp::tqdceq(char *)
{

};

void liboctqcpp::tqcceq(char *, int *, void *, void *)
{

};

void liboctqcpp::tqselceq(char *, void *)
{

};

void liboctqcpp::reset_conditions(char *, void *)
{

};

void liboctqcpp::Change_Status_Phase(char *, int, double, void *)
{

};

void liboctqcpp::tqlr(int, void *)
{

};

void liboctqcpp::tqlc(int, void *)
{

};

std::vector<double> liboctqcpp::PhaseFractions(void *ceq)
{
    std::vector<double> results =
    //----------------------------
    tqgetv("NP", -1, 0, tqgnp(ceq), ceq);
    //----------------------------
    return results;
};

std::vector<double> liboctqcpp::ConstituentFractions(int phase, void *ceq)
{
    std::vector<double> results =
    //-------------------------------
    tqgetv("X", phase, -1, tqgcn(ceq), ceq);
    //-------------------------------
    return results;
};

/*******************************************************************************
 * The following routine is for testing the interface functionality and produces
 * debug output. However it also demonstrates a usecase and can be used as a
 * starting point for a new implementation of OpenCalphad C++ interface.
*******************************************************************************/

int main(int argc, char *argv[])
{
    liboctqcpp OCASI;
    int n = 0;
    void * ceq = 0;

    //=================
    OCASI.tqini(n,&ceq);
    //=================

    cout << "-> Adress of ceq-Storage: [" << &ceq << "]" << endl;
/*
    //=====================================================
    std::vector<std::string> elnames
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
    cout << "]" << " [" << &ceq << "]" << endl;

*/
    vector<string> Elements;
    Elements.push_back("CR");
    Elements.push_back("FE");

    //================================================================
    std::vector<std::string> elnames2
    = OCASI.tqrpfil("TQ4lib/Cpp/Matthias/FECRMNC.TDB", Elements, &ceq);
    //= OCASI.tqrpfil("TQ4lib/Cpp/Matthias/steel1.TDB", Elements, &ceq);
    //================================================================

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

    //===============================
    std::vector<std::string> elnames3
    = OCASI.tqgcom(&ceq);
    //===============================

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

    //=================================
    int phasetuples = OCASI.tqgnp(&ceq);
    //=================================
    cout << "-> Number of phasetuples: [" << phasetuples << "]" << endl;

    std::vector<std::string> PhNames(phasetuples);

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
    std::vector<double> EquPhFr =

    //===========================
    OCASI.PhaseFractions(&ceq);
    //===========================

    //=================================
    phasetuples = OCASI.tqgnp(&ceq);
    //=================================
    cout << "-> Number of phasetuples: [" << phasetuples << "]" << endl;

    PhNames.resize(phasetuples);

    cout << "-> Phase Data: [";
    for(int i = 0; i < phasetuples; i++)
    {
        //==========================
        PhNames[i]=OCASI.tqgpn(i+1, &ceq);
        cout << PhNames[i];
        //==========================
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
        std::vector<double> ConFr =

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

    std::vector<int> ncons;
    std::vector<int> sites;
    double moles;

    for(int k = 0; k < phasetuples; k++)
    {
    if(k == phasetuples-1)
    {
        cout << "TQGPCN CANNOT BE CALLED FOR PHASE BCC_A2_AUTO#2, BECAUSE IT WILL BREAK!" << endl;
        break;
    }
    std::vector<double> y = OCASI.tqgphc1(k+1, ncons, sites, moles, &ceq);

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
            //OCASI.tqgpci(k+1, cname, &ceq)
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

    std::vector<double> G_TP;
    std::vector<double> G_Y;
    std::vector<double> G_YT;
    std::vector<double> G_YP;
    std::vector<double> G_YY;

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

    return 0;
}
