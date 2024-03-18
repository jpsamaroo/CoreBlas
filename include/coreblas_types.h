/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 **/
#ifndef COREBLAS_TYPES_H
#define COREBLAS_TYPES_H

#include <complex.h>

/*
 * RELEASE is a, b, c
 */
#define COREBLAS_VERSION_MAJOR   23
#define COREBLAS_VERSION_MINOR   8
#define COREBLAS_VERSION_PATCH   2

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************/
#if defined(COREBLAS_HAVE_MKL) || defined(COREBLAS_WITH_MKL)
#define lapack_complex_float coreblas_complex32_t
#define lapack_complex_double coreblas_complex64_t
#endif

/***************************************************************************//**
 *
 *  Some CBLAS routines take scalars by value in real arithmetic
 *  and by pointer in complex arithmetic.
 *  In precision generation, CBLAS_SADDR is removed from real arithmetic files.
 *
 **/
#ifndef CBLAS_SADDR
#if defined(COREBLAS_WITH_OPENBLAS)
#define CBLAS_SADDR(var) ((void*)&(var))
#else
#define CBLAS_SADDR(var) &(var)
#endif
#endif

/******************************************************************************/
enum {
    CoreBlasByte          = 0,
    CoreBlasInteger       = 1,
    CoreBlasRealFloat     = 2,
    CoreBlasRealDouble    = 3,
    CoreBlasComplexFloat  = 4,
    CoreBlasComplexDouble = 5
};

/***************************************************************************//**
 *
 *  COREBLAS constants - CBLAS & LAPACK.
 *  The naming and numbering is consistent with:
 *
 *    - CBLAS - http://www.netlib.org/blas/blast-forum/cblas.tgz,
 *    - LAPACKE - http://www.netlib.org/lapack/lapwrapc/.
 *
 *  During precision generation, CoreBlas_ConjTrans is conveted to CoreBlasTrans,
 *  while CoreBlasConjTrans is preserved.
 *
 **/
enum {
    CoreBlasInvalid       = -1,

    CoreBlasNoTrans       = 111,
    CoreBlasTrans         = 112,
    CoreBlasConjTrans     = 113,
    CoreBlas_ConjTrans    = CoreBlasConjTrans,

    CoreBlasUpper         = 121,
    CoreBlasLower         = 122,
    CoreBlasGeneral       = 123,
    CoreBlasGeneralBand   = 124,

    CoreBlasNonUnit       = 131,
    CoreBlasUnit          = 132,

    CoreBlasLeft          = 141,
    CoreBlasRight         = 142,

    CoreBlasOneNorm       = 171,
    CoreBlasRealOneNorm   = 172,
    CoreBlasTwoNorm       = 173,
    CoreBlasFrobeniusNorm = 174,
    CoreBlasInfNorm       = 175,
    CoreBlasRealInfNorm   = 176,
    CoreBlasMaxNorm       = 177,
    CoreBlasRealMaxNorm   = 178,

    CoreBlasNoVec         = 301,
    CoreBlasVec           = 302,
    CoreBlasCount         = 303,
    CoreBlasIVec          = 304,
    CoreBlasAllVec        = 305,
    CoreBlasSomeVec       = 306,

    CoreBlasRangeAll      = 351,
    CoreBlasRangeV        = 352,
    CoreBlasRangeI        = 353,

    CoreBlasForward       = 391,
    CoreBlasBackward      = 392,

    CoreBlasColumnwise    = 401,
    CoreBlasRowwise       = 402,

    CoreBlasW             = 501,
    CoreBlasA2            = 502,
    CoreBlas_Const_Limit  // Ensure always last.
};

enum {
    CoreBlasSuccess = 0,
    CoreBlasErrorNotInitialized,
    CoreBlasErrorNotSupported,
    CoreBlasErrorIllegalValue,
    CoreBlasErrorOutOfMemory,
    CoreBlasErrorNullParameter,
    CoreBlasErrorInternal,
    CoreBlasErrorSequence,
    CoreBlasErrorComponent,
    CoreBlasErrorEnvironment
};

enum {
    CoreBlasInplace,
    CoreBlasOutplace
};

enum {
    CoreBlasFlatHouseholder,
    CoreBlasTreeHouseholder
};

enum {
    CoreBlasDisabled = 0,
    CoreBlasEnabled = 1
};

enum {
    CoreBlasTuning,
    CoreBlasNb,
    CoreBlasIb,
    CoreBlasInplaceOutplace,
    CoreBlasNumPanelThreads,
    CoreBlasHouseholderMode
};

/******************************************************************************/
typedef int coreblas_enum_t;

typedef float  _Complex coreblas_complex32_t;
typedef double _Complex coreblas_complex64_t;

/******************************************************************************/
coreblas_enum_t coreblas_eigt_const(char lapack_char);
coreblas_enum_t coreblas_job_const(char lapack_char);
coreblas_enum_t coreblas_range_const(char lapack_char);
coreblas_enum_t coreblas_diag_const(char lapack_char);
coreblas_enum_t coreblas_direct_const(char lapack_char);
coreblas_enum_t coreblas_norm_const(char lapack_char);
coreblas_enum_t coreblas_side_const(char lapack_char);
coreblas_enum_t coreblas_storev_const(char lapack_char);
coreblas_enum_t coreblas_trans_const(char lapack_char);
coreblas_enum_t coreblas_uplo_const(char lapack_char);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // COREBLAS_TYPES_H
