#include "../octqc.h"
extern void	c_tqini(int, void **);
extern void	c_tqrfil(char *, void **);
extern void	c_tqgcom(int *, char[][24], void **);
extern void	c_tqrpfil(char *, int, char **, void **);
extern void	c_tqgnp(int *, void **);
extern void	c_tqgpn(int, char *, void **);
extern void	c_tqgetv(char *, int, int, int *, double *, void **);
extern void	examine_gtp_equilibrium_data(void *);
extern int      c_ntup;
extern int      c_nel;
extern int     	c_maxc;
extern char     *c_cnam[];
//extern void	c_tqgnp(int, gtp_equilibrium_data **);

