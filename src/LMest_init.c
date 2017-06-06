#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(back)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(bwforback)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(for_mult)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(nr_multilogit)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(prob_multilogif)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(sum_y)(void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"back",            (DL_FUNC) &F77_NAME(back),            13},
    {"bwforback",       (DL_FUNC) &F77_NAME(bwforback),        8},
    {"for_mult",        (DL_FUNC) &F77_NAME(for_mult),        12},
    {"nr_multilogit",   (DL_FUNC) &F77_NAME(nr_multilogit),   10},
    {"prob_multilogif", (DL_FUNC) &F77_NAME(prob_multilogif),  9},
    {"sum_y",           (DL_FUNC) &F77_NAME(sum_y),            6},
    {NULL, NULL, 0}
};

void R_init_LMest(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
