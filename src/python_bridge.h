/* Generated by Cython 0.28.2 */

#ifndef __PYX_HAVE__python_bridge
#define __PYX_HAVE__python_bridge


#ifndef __PYX_HAVE_API__python_bridge

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

__PYX_EXTERN_C void setup_params(void);

#endif /* !__PYX_HAVE_API__python_bridge */

/* WARNING: the interface of the module init function changed in CPython 3.5. */
/* It now returns a PyModuleDef instance instead of a PyModule instance. */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initpython_bridge(void);
#else
PyMODINIT_FUNC PyInit_python_bridge(void);
#endif

#endif /* !__PYX_HAVE__python_bridge */
