#ifndef __PYX_HAVE__pyoctq
#define __PYX_HAVE__pyoctq


#ifndef __PYX_HAVE_API__pyoctq

#ifndef __PYX_EXTERN_C
  #ifdef __cplusplus
    #define __PYX_EXTERN_C extern "C"
  #else
    #define __PYX_EXTERN_C extern
  #endif
#endif

#ifndef DL_IMPORT
  #define DL_IMPORT(_T) _T
#endif

__PYX_EXTERN_C DL_IMPORT(struct gtp_equilibrium_data) *ceq;

#endif /* !__PYX_HAVE_API__pyoctq */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initpyoctq(void);
#else
PyMODINIT_FUNC PyInit_pyoctq(void);
#endif

#endif /* !__PYX_HAVE__pyoctq */
