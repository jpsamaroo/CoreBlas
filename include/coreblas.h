/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 **/
#ifndef COREBLAS_CORE_BLAS_H
#define COREBLAS_CORE_BLAS_H

#include <stdio.h>
#include "coreblas_workspace.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef OPENBLAS_USE64BITINT
typedef BLASLONG blasint;
#else
typedef int blasint;
#endif
/***************************************************************************//**
 * This is just for translating enums into appropriate single characters; we
 * will only return the first character of the string; for compatibility with
 * earlier Fortran code.
 ******************************************************************************/
static const char *lapack_constants[] = {
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",

    "", "", "", "", "", "", "", "", "", "",
    "",
    "NoTrans",                              ///< 111: CoreBlasNoTrans
    "Trans",                                ///< 112: CoreBlasTrans
    "ConjTrans",                            ///< 113: CoreBlasConjTrans

    "", "", "", "", "", "", "",
    "Upper",                                ///< 121: CoreBlasUpper
    "Lower",                                ///< 122: CoreBlasLower
    "General",                              ///< 123: CoreBlasGeneral

    "", "", "", "", "", "", "",
    "NonUnit",                              ///< 131: CoreBlasNonUnit
    "Unit",                                 ///< 132: CoreBlasUnit

    "", "", "", "", "", "", "", "",
    "Left",                                 ///< 141: CoreBlasLeft
    "Right",                                ///< 142: CoreBlasRight

    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "",
    "One",                                  ///< 171: CoreBlasOneNorm
    "",                                     ///< 172: CoreBlasRealOneNorm
    "Two",                                  ///< 173: CoreBlasTwoNorm
    "Frobenius",                            ///< 174: CoreBlasFrobeniusNorm
    "Infinity",                             ///< 175: CoreBlasInfNorm
    "",                                     ///< 176: CoreBlasRealInfNorm
    "Maximum",                              ///< 177: CoreBlasMaxNorm
    "",                                     ///< 178: CoreBlasRealMaxNorm

    "", "", "", "",                         // 182. 
    "", "", "", "", "", "", "", "", "", "", // 192.
    "", "", "", "", "", "", "", "", "", "", // 202.
    "", "", "", "", "", "", "", "", "", "", // 212.
    "", "", "", "", "", "", "", "", "", "", // 222.
    "", "", "", "", "", "", "", "", "", "", // 232.
    "", "", "", "", "", "", "", "", "", "", // 242.
    "", "", "", "", "", "", "", "", "", "", // 252.
    "", "", "", "", "", "", "", "", "", "", // 262.
    "", "", "", "", "", "", "", "", "", "", // 272.
    "", "", "", "", "", "", "", "", "", "", // 282.
    "", "", "", "", "", "", "", "", "", "", // 292.
    "", "", "", "", "", "", "", "", "", "", // 302.
    "", "", "", "", "", "", "", "", "", "", // 312.
    "", "", "", "", "", "", "", "", "", "", // 322.
    "", "", "", "", "", "", "", "", "", "", // 332.
    "", "", "", "", "", "", "", "", "", "", // 342.
    "", "", "", "", "", "", "", "", "", "", // 352.
    "", "", "", "", "", "", "", "", "", "", // 362.
    "", "", "", "", "", "", "", "", "", "", // 372.
    "", "", "", "", "", "", "", "", "", "", // 382.
    "", "", "", "", "", "", "", "",         // 390.
    "Forward",                              ///< 391: CoreBlasForward
    "Backward",                             ///< 392: CoreBlasBackward
    "", "", "", "", "", "", "", "",         // 400.
    "Columnwise",                           ///< 401: CoreBlasColumnwise
    "Rowwise"                               ///< 402: CoreBlasRowwise
    "", "", "", "", "", "", "", "", "", "", // 412.
    "", "", "", "", "", "", "", "", "", "", // 422.
    "", "", "", "", "", "", "", "", "", "", // 432.
    "", "", "", "", "", "", "", "", "", "", // 442.
    "", "", "", "", "", "", "", "", "", "", // 452.
    "", "", "", "", "", "", "", "", "", "", // 462.
    "", "", "", "", "", "", "", "", "", "", // 472.
    "", "", "", "", "", "", "", "", "", "", // 482.
    "", "", "", "", "", "", "", "", "", "", // 492.
    "", "", "", "", "", "", "", "",         // 500.
    "W",                                    ///< 501: CoreBlasW.
    "A2",                                   ///< 502: CoreBlasA2.
    ""                                      ///< ...: Invalid const.
};

/***************************************************************************//**
 * @retval LAPACK character constant corresponding to COREBLAS constant
 * @ingroup coreblas_const
 ******************************************************************************/
static inline char lapack_const(int coreblas_const) {
    int entries = sizeof(lapack_constants)/sizeof(char*);
    if (coreblas_const < 0 || coreblas_const >= entries)
        return lapack_constants[entries-1][0];
    return lapack_constants[coreblas_const][0];
}

#define coreblas_error(msg) \
        coreblas_error_func_line_file(__func__, __LINE__, __FILE__, msg)

static inline void coreblas_error_func_line_file(
    char const *func, int line, const char *file, const char *msg)
{
    fprintf(stderr,
            "COREBLAS ERROR at %d of %s() in %s: %s\n",
            line, func, file, msg);
}

#ifdef __cplusplus
}  // extern "C"
#endif

#include "coreblas_s.h"
#include "coreblas_d.h"
#include "coreblas_ds.h"
#include "coreblas_c.h"
#include "coreblas_z.h"
#include "coreblas_zc.h"

#endif // COREBLAS_CORE_BLAS_H
