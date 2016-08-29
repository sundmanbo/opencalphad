#if !defined __OCASI__
#define __OCASI__

/* Modification history
160829 Bo Sundman Update
2015-2016 Matthias Stratmann and Cristophe Sigli Modifications
2014 Teslos? First version

This contains the structure of TYPE variables in OC needed for the OC/TQ OCASI interface 

NOTE there is also a c_gtp_equilibrium_data structure defined in liboctqisoc.F90 */

typedef struct {
  int forcenewcalc;
  double tpused[2];
  double results[6];
} tpfun_parres;

typedef struct {
  int splink, phlink, status;
  char refstate[16];
  int *endmember;
  double tpref[2];
  double chempot[2];
  double mass, molat;
} gtp_components;

typedef struct {
  int lokph, compset, ixphase, lokvares, nextcs;
} gtp_phasetuple;

typedef struct {
  int statevarid, norm, unit, phref, argtyp;
  int phase, compset, component, constituent;
  double coeff;
  int oldstv;
} gtp_state_variable;

typedef struct {
  int latd, ndd, tnoofxfr, tnoofyfr, varreslink, totdis;
  char id;
  double *dsites;
  int *nooffr;
  int *splink;
  int *y2x;
  double *dxidyj;
  double fsites;
} gtp_fraction_set;

//struct gtp_fraction_set;

typedef struct {
  int nextfree, phlink, status2, phstate,phtupx;
  double abnorm[3];
  char prefix[4], suffix[4];
  int *constat;
  double *yfr;
  double *mmyfr;
  double *sites;
  double *dpqdy;
  double *d2pqdvay;
  //struct gtp_fraction_set disfra;
  double amfu, netcharge, dgm;
  int nprop;
  int *listprop;
  double **gval;
  double ***dgval;
  double **d2gval;
  double curlat[3][3];
  double **cinvy;
  double *cxmol;
  double **cdxmol;
} gtp_phase_varres;

typedef struct gtp_condition {
  int noofterms, statev, active, iunit, nid, iref, seqz, experimenttype;
  int symlink1, symlink2;
  int **indices;
  double *condcoeff;
  double *prescribed, current, uncertainity;
  // should this be a struct ??
  gtp_state_variable *statvar;
  struct gtp_condition *next, *previous;
} gtp_condition;

typedef struct {
  int status, multiuse, eqno, next;
  char eqname[24], comment[72];
  double tpval[2], rtn;
  double weight;
  double *svfunres;
  gtp_condition *lastcondition, *lastexperiment;
  gtp_components *complist;
  double **compstoi, **invcompstoi;
  gtp_phase_varres *phase_varres;
  tpfun_parres *eq_tpres;
  double *cmuval;
  double xconv;
  double gmindif;
  int maxiter;
  char eqextra[80];
  int sysmatdim, nfixmu, nfixph;
  int *fixmu;
  int *fixph;
  double **savesysmat;
} gtp_equilibrium_data; 
 
#endif

