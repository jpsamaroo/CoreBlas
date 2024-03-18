/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zlascl.c, normal z -> s, Mon Mar 18 06:35:34 2024
 *
 **/

#include <coreblas.h>
#include "coreblas_types.h"
#include "core_lapack.h"

/******************************************************************************/
__attribute__((weak))
void coreblas_slascl(coreblas_enum_t uplo,
                 float cfrom, float cto,
                 int m, int n,
                 float *A, int lda)
{
    // LAPACKE_slascl is not available in LAPACKE < 3.6.0
    int kl;
    int ku;
    int info;
    char type = lapack_const(uplo);
    LAPACK_slascl(&type,
                  &kl, &ku,
                  &cfrom, &cto,
                  &m, &n,
                  A, &lda, &info);
}

/******************************************************************************/
void coreblas_kernel_slascl(coreblas_enum_t uplo,
                     float cfrom, float cto,
                     int m, int n,
                     float *A, int lda)
{

    coreblas_slascl(uplo,
                cfrom, cto,
                m, n,
                A, lda);
}
