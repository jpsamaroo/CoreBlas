/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zlascl.c, normal z -> c, Mon Mar 18 06:35:34 2024
 *
 **/

#include <coreblas.h>
#include "coreblas_types.h"
#include "core_lapack.h"

/******************************************************************************/
__attribute__((weak))
void coreblas_clascl(coreblas_enum_t uplo,
                 float cfrom, float cto,
                 int m, int n,
                 coreblas_complex32_t *A, int lda)
{
    // LAPACKE_clascl is not available in LAPACKE < 3.6.0
    int kl;
    int ku;
    int info;
    char type = lapack_const(uplo);
    LAPACK_clascl(&type,
                  &kl, &ku,
                  &cfrom, &cto,
                  &m, &n,
                  A, &lda, &info);
}

/******************************************************************************/
void coreblas_kernel_clascl(coreblas_enum_t uplo,
                     float cfrom, float cto,
                     int m, int n,
                     coreblas_complex32_t *A, int lda)
{

    coreblas_clascl(uplo,
                cfrom, cto,
                m, n,
                A, lda);
}
