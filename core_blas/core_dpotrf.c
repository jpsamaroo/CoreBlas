/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zpotrf.c, normal z -> d, Mon Mar 18 06:35:36 2024
 *
 **/

#include <coreblas.h>
#include "coreblas_types.h"
#include "core_lapack.h"

/***************************************************************************//**
 *
 * @ingroup core_potrf
 *
 *  Performs the Cholesky factorization of a symmetric positive definite
 *  matrix A. The factorization has the form
 *
 *    \f[ A = L \times L^T, \f]
 *    or
 *    \f[ A = U^T \times U, \f]
 *
 *  where U is an upper triangular matrix and L is a lower triangular matrix.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - CoreBlasUpper: Upper triangle of A is stored;
 *          - CoreBlasLower: Lower triangle of A is stored.
 *
 * @param[in] n
 *          The order of the matrix A. n >= 0.
 *
 * @param[in,out] A
 *          On entry, the symmetric positive definite matrix A.
 *          If uplo = CoreBlasUpper, the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the strictly
 *          lower triangular part of A is not referenced.
 *          If uplo = CoreBlasLower, the leading N-by-N lower triangular part of A
 *          contains the lower triangular part of the matrix A, and the strictly
 *          upper triangular part of A is not referenced.
 *          On exit, if return value = 0, the factor U or L from the Cholesky
 *          factorization A = U^T*U or A = L*L^T.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,n).
 *
 ******************************************************************************/
__attribute__((weak))
int coreblas_dpotrf(coreblas_enum_t uplo,
                 int n,
                 double *A, int lda)
{
    return LAPACKE_dpotrf_work(LAPACK_COL_MAJOR,
                               lapack_const(uplo),
                               n,
                               A, lda);
}

/******************************************************************************/
void coreblas_kernel_dpotrf(coreblas_enum_t uplo,
                     int n,
                     double *A, int lda,
                     int iinfo)
{

    int info = coreblas_dpotrf(uplo,
                           n,
                           A, lda);
    if (info != 0)
        coreblas_error("core_dpotrf() failed");
}
