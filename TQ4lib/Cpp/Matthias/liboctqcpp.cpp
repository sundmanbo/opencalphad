#include "liboctqcpp.h"
using namespace std;

#define MAXEL 24;

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

std::vector<std::string> liboctqcpp::tqgcom(void * ceq)
{
    int n = MAXEL;
    char elnames[24];
    std::vector<std::string> result;
    //=======================
    c_tqgcom(&n, elnames, ceq);
    //=======================
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
    //==============
    c_tqgnp(&n, ceq);
    //==============
    return n;
};

std::string liboctqcpp::tqgpn(int i, void * ceq)
{
    char phname[24];
    std::string result;
    //=====================
    c_tqgpn(i, phname, ceq);
    //=====================
    result = phname;
    return result;
};

int liboctqcpp::tqgpi(std::string pname, void * ceq)
{
    char *phasename = strcpy((char*)malloc(pname.length()+1), pname.c_str());
    int *i;
    //=========================
    c_tqgpi(i, phasename, ceq);
    //=========================
};

void liboctqcpp::tqgpcn2(int, int, char *, void *)
{

};

void liboctqcpp::tqgpcn(int, int, char *, void *)
{

};

void liboctqcpp::tqgpci(int, int *, char *, void *)
{

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
    //=========================================
    c_tqsetc(name, n1, n2, val, &cnum, ceq);
    //=========================================
};

void liboctqcpp::tqce(void * ceq)
{
    char target[60] = " ";
    int null1 = 0;
    int null2 = 0;
    double val;
    //======================================
    c_tqce(target, null1, null2, &val, ceq);
    //======================================
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
    //=======================================
    c_tqgetv(name, n1, n2, &n3, val, ceq);
    //=======================================
    for(int i = 0; i < n3; i++)
    results[i] = val[i];
    return results;
};

void liboctqcpp::tqgphc1(int, int * , int *, int *, double *, double *, double *, void *)
{

};

void liboctqcpp::tqsphc1(int, double *, double *, void *)
{

};

void liboctqcpp::tqcph1(int, int, int *, double *, double *, double *, double *, double *, void *)
{

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
    int nph = c_ntup; //TODO: Remove c_ntup
    std::vector<double> results(nph);
    results = tqgetv("NP", -1, 0,nph, ceq);
    return results;
};

std::vector<double> liboctqcpp::ConstituentFractions(int phase, void *ceq)
{
    int nel = c_nel; //TODO: Remove c_ntup
    std::vector<double> results(nel);
    results = tqgetv("X", phase, -1,nel, ceq);
    return results;
};

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
    Elements.push_back("C");
    Elements.push_back("CR");
    Elements.push_back("FE");
    Elements.push_back("MN");

    //================================================================
    std::vector<std::string> elnames2
    = OCASI.tqrpfil("TQ4lib/Cpp/Matthias/FECRMNC.TDB", Elements, &ceq);
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

    std::vector<std::string> PhNames(c_ntup);

    cout << "-> Phase Data: [";
    for(int i = 0; i < c_ntup; i++)
    {
        //==========================
        PhNames[i]=OCASI.tqgpn(i+1, &ceq);
        cout << PhNames[i];
        //==========================
        if(i < c_ntup-1)
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
    cout << "-> Set X(CR) to: [" << XCR << "]" << " [" << &ceq << "]" << endl;

    //=========================================
    T = OCASI.tqgetv("T", 0, 0, &ceq);
    P = OCASI.tqgetv("P", 0, 0, &ceq);
    N = OCASI.tqgetv("N", 0, 0, &ceq);
    XCR = OCASI.tqgetv("X", 1, 0, &ceq);
    //=========================================

    cout << "-> Temperature set to: [" << T << "]" << " [" << &ceq << "]" << endl;
    cout << "-> Ambient Pressure set to: [" << P << "]" << " [" << &ceq << "]" << endl;
    cout << "-> Moles set to: [" << N << "]" << " [" << &ceq << "]" << endl;
    cout << "-> X(CR) set to: [" << XCR << "]" << " [" << &ceq << "]" << endl;

    //===============
    OCASI.tqce(&ceq);
    //===============

    cout << "-> Calculated Equilibrium [" << ceq << "]" << endl;
    std::vector<double> EquPhFr =

    //===========================
    OCASI.PhaseFractions(&ceq);
    //===========================

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

    return 0;
}
