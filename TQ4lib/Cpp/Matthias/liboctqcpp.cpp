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

vector<string> liboctqcpp::tqrfil(string fname, void * ceq)
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
    free(filename);

    //=================
    return tqgcom(ceq);
    //=================
};

vector<string> liboctqcpp::tqrpfil(string fname, vector<string> elnames, void * ceq)
{
    char *filename = strcpy((char*)malloc(fname.length()+1), fname.c_str());
    char *selel[elnames.size()];
    char *tempchar;
    for(int i = 0; i < elnames.size(); i++)
    {
        tempchar
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
    free(filename);
    free(tempchar);
    //free(selel);
    vector<string> asdf = tqgcom(ceq);
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

vector<string> liboctqcpp::tqgcom(void * ceq)
{
    int n = MAXEL;
    char elnames[24];
    vector<string> result;
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

string liboctqcpp::tqgpn(int i, void * ceq)
{
    char phname[24];
    string result;
    //======================
    c_tqgpn(i, phname, ceq);
    //======================
    result = phname;
    return result;
};

int liboctqcpp::tqgpi(string pname, void * ceq)
{
    char *phasename = strcpy((char*)malloc(pname.length()+1), pname.c_str());
    int i;
    //=========================
    c_tqgpi(&i, phasename, ceq);
    //=========================
    free(phasename);
    return i;
};

string liboctqcpp::tqgpcn2(int phidx, int i, void * ceq)
{
    //---------------------------
    return tqgpcn(phidx, i, ceq);
    //---------------------------
};

string liboctqcpp::tqgpcn(int phidx, int i, void * ceq)
{
    char constituentname[24];
    string result;
    //=======================================
    c_tqgpcn(phidx, i, constituentname, ceq); //TODO: c_tqgpcn is not implemented in liboctq.F90!!
    //=======================================
    result = constituentname;
    return result;
};

int liboctqcpp::tqgpci(int phidx, string cname, void * ceq)
{
    int c;
    char *constituent = strcpy((char*)malloc(cname.length()+1), cname.c_str());
    //====================================
    c_tqgpci(phidx, &c, constituent, ceq); //TODO: c_tqgpci is not implemented in liboctq.F90!!
    //====================================
    free(constituent);
    return c;
};

vector<double> liboctqcpp::tqgpcs(int phidx, int con, double& mass, void * ceq)
{
    vector<double> result;
    double * stoi;
    c_tqgpcs(phidx, con, stoi, &mass, ceq);
    free(stoi);
    return result;
};

void liboctqcpp::tqgccf(int comp, void * ceq)
{
    int nel;
    char * elnames;
    double stoi;
    double mass;
    //===============================================
    c_tqgccf(comp, &nel, elnames, &stoi, &mass, ceq); //TODO: c_tqgccf is not implemented in liboctq.F90
    //===============================================
    free(elnames);
};

int liboctqcpp::tqgnpc(int phidx, void * ceq)
{
    int nc;
    //========================
    c_tqgnpc(phidx, &nc, ceq); //TODO: c_tqgnpc is not implemented in liboctq.F90
    //========================
    return nc;
};

void liboctqcpp::tqphtupsts(int phidx, int newstatus, double val, void * ceq)
{
    // phidx < 0 means "all phases"
    // newstatus -4 hidden
    // newstatus -3 suspended
    // newstatus -2 dormant
    // newstatus -1 TODO: entered?
    // newstatus  0 TODO: entered?
    // newstatus  1 entered
    // newstatus  2 fix

    //=====================================
    c_tqphtupsts(phidx, newstatus, val, ceq);
    //=====================================
};

void liboctqcpp::tqsetc(string par, int n1, int n2, double val, void * ceq)
{
    int cnum;
    char *name = strcpy((char*)malloc(par.length()+1), par.c_str());
    //======================================
    c_tqsetc(name, n1, n2, val, &cnum, ceq);
    //======================================
    free(name);
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
    free(name);
    return val;
};

vector<double> liboctqcpp::tqgetv(string par, int n1, int n2, int n3, void * ceq)
{
    vector<double> results(n3);
    char *name = strcpy((char*)malloc(par.length()+1), par.c_str());
    double val[n3];
    //====================================
    c_tqgetv(name, n1, n2, &n3, val, ceq);
    //====================================
    for(int i = 0; i < n3; i++)
    results[i] = val[i];
    free(name);
    return results;
};

vector<double> liboctqcpp::tqgphc1(int phIdx, vector<int>& ncons,
                                        vector<int>& sites, double& moles,
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

void liboctqcpp::tqsphc1(int phidx, vector<double> y, void * ceq)
{
    double extra[MAXPH];
    double yfr[y.size()];
    for(int i = 0; i < y.size(); i++)
    yfr[i] = y[i];
    //================================
    c_tqsphc1(phidx, yfr, extra, ceq);
    //================================
};

double liboctqcpp::tqcph1(int phidx, vector<double>& G_TP,
                          vector<double>& G_Y, vector<double>& G_YT,
                          vector<double>& G_YP, vector<double>& G_YY,
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

int liboctqcpp::tqcph2(int phidx, int type, void * ceq)
{
    int lokres; //ceq%phase_varres(lokres)% with all results
    int n3;
    //=======================================
    c_tqcph2(phidx, type, &n3, &lokres, ceq);
    //=======================================
    return lokres;
};

void liboctqcpp::tqdceq(string name)
{
    char *ceqname = strcpy((char*)malloc(name.length()+1), name.c_str());
    //================
    c_tqdceq(ceqname);
    //================
    free(ceqname);
};

int liboctqcpp::tqcceq(string name, void * newceq, void * ceq)
{
    char *ceqname = strcpy((char*)malloc(name.length()+1), name.c_str());
    int n1;
    //==================================
    c_tqcceq(ceqname, &n1, newceq, ceq);
    //==================================
    return n1;
};

void liboctqcpp::tqselceq(string name, void * ceq)
{
    char *ceqname = strcpy((char*)malloc(name.length()+1), name.c_str());
    //=======================
    c_tqselceq(ceqname, ceq);
    //=======================
    free(ceqname);
};

void liboctqcpp::reset_conditions(string condition, double newval, void * ceq)
{
    char *cond = strcpy((char*)malloc(condition.length()+1), condition.c_str());
    //============================
    c_reset_conditions(cond, ceq); //TODO: send newval to liboctq.F90
    //============================
    free(cond);
};

void liboctqcpp::Change_Status_Phase(string phname, int newstatus, double val, void * ceq)
{
    //--------------------------------------------------
    tqphtupsts(tqgpi(phname, ceq), newstatus, val, ceq);
    //--------------------------------------------------
};

void liboctqcpp::tqlr(int lut, void * ceq)
{
    //===============
    c_tqlr(lut, ceq);
    //===============
};

void liboctqcpp::tqlc(int lut, void * ceq)
{
    //===============
    c_tqlc(lut, ceq);
    //===============
};

vector<double> liboctqcpp::PhaseFractions(void *ceq)
{
    vector<double> results =
    //----------------------------
    tqgetv("NP", -1, 0, tqgnp(ceq), ceq);
    //----------------------------
    return results;
};

vector<double> liboctqcpp::ConstituentFractions(int phase, void *ceq)
{
    vector<double> results =
    //-------------------------------
    tqgetv("X", phase, -1, tqgcn(ceq), ceq);
    //-------------------------------
    return results;
};
