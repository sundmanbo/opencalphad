#include <string>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <vector>
#include <cstdio>

extern"C" int  c_nel;
extern"C" int  c_maxc;
extern"C" int  c_maxp;
extern"C" char * c_cnam[24];
extern"C" char * cnames[25];
extern"C" int c_ntup;
extern"C" int c_noofcs(int);
extern"C" int c_ierr(void);
extern"C"
{
    void c_tqini(int, void *);
    void c_tqrfil(char *, void *);
    void c_tqrpfil(char *, int, char **, void *);
    void c_tqgcom(int *, char *, void *);
    void c_tqgnp(int *, void *);
    void c_tqgpn(int, char *, void *);
    void c_tqgpi(int *, char *, void *);
    void c_tqgpcn2(int, int, char *, void *);
    void c_tqgpcn(int, int, char *, void *);
    void c_tqgpci(int, int *, char *, void *);
    void c_tqgpcs(int, int, double *, double *, void *);
    void c_tqgccf(int, int *, char *, double *, double *, void *);
    void c_tqgnpc(int, int *, void *);
    void c_tqphtupsts(int, int, double, void *);
    void c_tqsetc(char *, int, int, double, int *, void *);
    void c_tqce(char *, int, int, double *, void *);
    void c_tqgetv(char *, int, int, int *, double *, void *);
    void c_tqgphc1(int, int * , int *, int *, double *, double *, double *, void *);
    void c_tqsphc1(int, double *, double *, void *);
    void c_tqcph1(int, int, int *, double *, double *, double *, double *, double *, void *);
    void c_tqcph2(int, int, int *, int *, void *);
    void c_tqdceq(char *);
    void c_tqcceq(char *, int *, void *, void *);
    void c_tqselceq(char *, void *);
    void c_reset_conditions(char *, void *);
    void c_Change_Status_Phase(char *, int, double, void *);
    void c_tqlr(int, void *);
    void c_tqlc(int, void *);
}

class liboctqcpp
{
    public:
    void * ceq2;
    int ntup;
    int nel;
    std::vector<std::string> cnames;
    std::vector<std::string> cnam;
    void tqini(int n, void * ceq);
    std::vector<std::string> tqrfil(std::string fname, void * ceq);
    std::vector<std::string> tqrpfil(std::string fname, std::vector<std::string> elnames, void * ceq);
    int tqgcn(void * ceq);
    std::vector<std::string> tqgcom(void * ceq);
    int tqgnp(void * ceq);
    std::string tqgpn(int i, void * ceq);
    int tqgpi(std::string pname, void * ceq);
    std::string tqgpcn(int phidx, int i, void * ceq);
    std::string tqgpcn2(int phidx, int i, void * ceq);
    int tqgpci(int phidx, std::string cname, void * ceq);
    void tqgpcs(int, int, double *, double *, void *);
    std::vector<double> tqgpcs(int phidx, int con, double& mass, void * ceq);
    void tqgccf(int comp, void * ceq);
    int tqgnpc(int phidx, void * ceq);
    void tqphtupsts(int, int, double, void *);
    void tqsetc(std::string, int, int, double, void *);
    void tqce(void *);
    double tqgetv(std::string, int, int, void *);
    std::vector<double> tqgetv(std::string, int, int, int, void *);
    std::vector<double> tqgphc1(int phIdx, std::vector<int>& ncons, std::vector<int>& sites, double& moles, void * ceq);
    void tqsphc1(int phidx, std::vector<double> y, void * ceq);
    double tqcph1(int phidx, std::vector<double>& G_TP,
                  std::vector<double>& G_Y, std::vector<double>& G_YT,
                  std::vector<double>& G_YP, std::vector<double>& G_YY,
                  void * ceq);
    int tqcph2(int phidx, int type, void * ceq);
    void tqdceq(std::string);
    int tqcceq(std::string name, void * newceq, void * ceq);
    void tqselceq(std::string, void *);
    void reset_conditions(std::string condition, double newval, void * ceq);
    void Change_Status_Phase(std::string phname, int newstatus, double val, void * ceq);
    void tqlr(int, void *);
    void tqlc(int, void *);
    std::vector<double> PhaseFractions(void *);
    std::vector<double> ConstituentFractions(int phase, void *ceq);
};
