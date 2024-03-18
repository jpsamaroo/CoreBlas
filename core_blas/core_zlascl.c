/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @precisions normal z -> c d s
 *
 **/

#include <coreblas.h>
#include "coreblas_types.h"
#include "core_lapack.h"

/******************************************************************************/
__attribute__((weak))
void coreblas_zlascl(coreblas_enum_t uplo,
                 double cfrom, double cto,
                 int m, int n,
                 coreblas_complex64_t *A, int lda)
{
    // LAPACKE_zlascl is not available in LAPACKE < 3.6.0
    int kl;
    int ku;
    int info;
    char type = lapack_const(uplo);
    LAPACK_zlascl(&type,
                  &kl, &ku,
                  &cfrom, &cto,
                  &m, &n,
                  A, &lda, &info);
}

/******************************************************************************/
void coreblas_kernel_zlascl(coreblas_enum_t uplo,
                     double cfrom, double cto,
                     int m, int n,
                     coreblas_complex64_t *A, int lda)
{

    coreblas_zlascl(uplo,
                cfrom, cto,
                m, n,
                A, lda);
}
