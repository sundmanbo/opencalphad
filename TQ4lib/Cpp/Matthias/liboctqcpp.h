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
extern"C" int  c_ntup;
extern"C" int c_noofcs(int);
extern"C" int c_ierr(void);
extern"C"
{
    void c_tqini(int, void *);
    void c_tqrfil(char *, void *);
    void c_tqrpfil(char *, int, char **, void *);
    void c_tqgcom(int, char **, void *);
    void c_tqgnp(int, void *);
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
    int ntup;
    int nel;
    vector<string> cnames;
    vector<string> cnam;
    void tqini(int n, void * ceq)
    {
        //==============
        c_tqini(n, ceq);
        //==============
    };

    void tqrfil(string fname, void * ceq)
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
    };

    void tqrpfil(string fname, vector<string> elnames, void * ceq)
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
        for(int i = 0; i < nel; i++)
        cnam[i] = c_cnam[i];
    };

    void tqgcom(int, char **, void *)
    {

    };

    void tqgnp(int, void *)
    {

    };

    void tqgpn(int, char *, void *)
    {

    };

    void tqgpi(int *, char *, void *)
    {

    };

    void tqgpcn2(int, int, char *, void *)
    {

    };

    void tqgpcn(int, int, char *, void *)
    {

    };

    void tqgpci(int, int *, char *, void *)
    {

    };

    void tqgpcs(int, int, double *, double *, void *)
    {

    };

    void tqgccf(int, int *, char *, double *, double *, void *)
    {

    };

    void tqgnpc(int, int *, void *)
    {

    };

    void tqphtupsts(int, int, double, void *)
    {

    };

    void tqsetc(char *, int, int, double, int *, void *)
    {

    };

    void tqce(char *, int, int, double *, void *)
    {

    };

    void tqgetv(char *, int, int, int *, double *, void *)
    {

    };

    void tqgphc1(int, int * , int *, int *, double *, double *, double *, void *)
    {

    };

    void tqsphc1(int, double *, double *, void *)
    {

    };

    void tqcph1(int, int, int *, double *, double *, double *, double *, double *, void *)
    {

    };

    void tqcph2(int, int, int *, int *, void *)
    {

    };

    void tqdceq(char *)
    {

    };

    void tqcceq(char *, int *, void *, void *)
    {

    };

    void tqselceq(char *, void *)
    {

    };

    void reset_conditions(char *, void *)
    {

    };

    void Change_Status_Phase(char *, int, double, void *)
    {

    };

    void tqlr(int, void *)
    {

    };

    void tqlc(int, void *)
    {

    };

};
