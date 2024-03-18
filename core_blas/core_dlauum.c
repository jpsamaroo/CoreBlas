/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zlauum.c, normal z -> d, Mon Mar 18 06:35:35 2024
 *
 **/

#include <coreblas.h>
#include "coreblas_types.h"
#include "core_lapack.h"

/***************************************************************************//**
 *
 * @ingroup core_lauum
 *
 *  Computes the product U * U^T or L^T * L, where the triangular
 *  factor U or L is stored in the upper or lower triangular part of
 *  the array A.
 *
 *  If uplo = 'U' or 'u' then the upper triangle of the result is stored,
 *  overwriting the factor U in A.
 *  If uplo = 'L' or 'l' then the lower triangle of the result is stored,
 *  overwriting the factor L in A.

 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = CoreBlasUpper: Upper triangle of A is stored;
 *          = CoreBlasLower: Lower triangle of A is stored.
 *
 *
 * @param[in] n
 *          The order of the matrix A. n >= 0.
 *
 * @param[in,out] A
 *          On entry, the triangular factor U or L.
 *          On exit, if uplo = 'U', the upper triangle of A is
 *          overwritten with the upper triangle of the product U * U^T;
 *          if uplo = 'L', the lower triangle of A is overwritten with
 *          the lower triangle of the product L^T * L.

 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,n).
 *
 * @param[out] info
 *          - 0 on successful exit
 *          - < 0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
__attribute__((weak))
int coreblas_dlauum(coreblas_enum_t uplo,
                int n,
                double *A, int lda)
{
    return LAPACKE_dlauum_work(LAPACK_COL_MAJOR,
                        lapack_const(uplo), n, A, lda);
}

/******************************************************************************/
void coreblas_kernel_dlauum(coreblas_enum_t uplo,
                     int n,
                     double *A, int lda)
{

    int info = coreblas_dlauum(uplo, n, A, lda);
    if (info != CoreBlasSuccess) {
        coreblas_error("core_dlauum() failed");
    }
}
