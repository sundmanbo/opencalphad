#if !defined __OCTQC__
#define __OCTQC__

typedef struct {
	int statevarid, norm, unit, phref, argtyp;
	int phase, compset, component, constituent;
	double coeff;
	int oldstv;
} gtp_state_variable;

struct gtp_fraction_set;

typedef struct {
	int nextfree, phlink, status2, phstate;
	double abnorm[2];
	char prefix[4], suffix[4];
	int *constat;
	double *yfr;
	double *mmyfr;
	double *sites;
	double *dpqdy;
	double *d2pqdvay;
	//struct gtp_fraction_set disfra;
	double amfu, netcharge, dgm, amcom, damount;
	int nprop, ncc;
	int *listprop;
	double **gval;
	double ***dgval;
	double **d2gval;
	double curlat[3][3];
	double **cinvy;
	double *cxmol;
	double **cdxmol;
} gtp_phase_varres;

typedef struct gtp_fraction_set {
    int latd, ndd, tnoofxfr, tnoofyfr, varreslink, totdis;
    char id;
    double *dsites;
    int *nooffr;
    int *splink;
    int *y2x;
    double *dxidyj;
    double fsites;
    gtp_phase_varres *phdapointer;
} gtp_fraction_set;

typedef struct {
	int splink, phlink, status;
	char refstate[16];
	int endmember;
	double tpref[2];
	double chempot[2];
	double mass, molat;
} gtp_components;

typedef struct gtp_condition {
	int noofterms, statev, active, iunit, nid, iref, seqz;
	int symlink1, symlink2;
	int **indices;
	double *condcoeff;
	double *prescribed, current, uncertainity;
	struct gtp_condition *next, *previous;
	gtp_state_variable *statvar;
} gtp_condition;

typedef struct {
	double tpused[2];
	double results[6];
} tpfun_parres;

typedef struct {
	int status, multiuse, eqno, next;
	char eqname[24];
	double tpval[2], rtn;
	double *svfunres;
	gtp_condition *lastcondition, *lastexperiment;
	gtp_components *complist, **compstoi, **invcompstoi;
	gtp_phase_varres *phase_varres;
	tpfun_parres *eq_tpres;
	double *cmuval;
	double xconv;
	double gmindif;
	int maxiter;
	int sysmatdim, nfixmu, nfixph;
	int *fixmu;
	int *fixph;
	double **savesysmat;
} gtp_equilibrium_data; 
 
#endif

