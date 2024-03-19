/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zhegst.c, normal z -> d, Mon Mar 18 06:35:36 2024
 *
 **/

#include <coreblas.h>
#include "coreblas_types.h"
#include "core_lapack.h"

/***************************************************************************//**
 *
 * @ingroup core_hegst
 *
 *  Reduces a complex symmetric-definite generalized eigenproblem to standard
 *  form.
 *
 *  If ITYPE = 1, the problem is A*x = lambda*B*x,
 *  and A is overwritten by inv(U^T)*A*inv(U) or inv(L)*A*inv(L^T)
 *
 *  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
 *  B*A*x = lambda*x, and A is overwritten by U*A*U^T or L^T*A*L.
 *
 *******************************************************************************
 *
 * @param[in] itype
 *          = 1: compute inv(U^T)*A*inv(U) or inv(L)*A*inv(L^T);
 *          = 2 or 3: compute U*A*U^T or L^T*A*L.
 *
 * @param[in] uplo
 *          If CoreBlasUpper, upper triangle of A is stored and B is factored as
 *          U^T*U;
 *          If CoreBlasLower, lower triangle of A is stored and B is factored as
 *          L*L^T.
 *
 * @param[in] n
 *          The order of the matrices A and B.  N >= 0.
 *
 * @param[in,out] A
 *          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
 *          N-by-N upper triangular part of A contains the upper
 *          triangular part of the matrix A, and the strictly lower
 *          triangular part of A is not referenced.  If UPLO = 'L', the
 *          leading N-by-N lower triangular part of A contains the lower
 *          triangular part of the matrix A, and the strictly upper
 *          triangular part of A is not referenced.
 *
 *          On exit, if INFO = 0, the transformed matrix, stored in the
 *          same format as A.
 *
 * @param[in] lda
 *          The leading dimension of the array A.  LDA >= max(1,N).
 *
 * @param[in,out] B
 *          The triangular factor from the Cholesky factorization of B,
 *          as returned by DPOTRF.
 *
 * @param[in] ldb
 *          The leading dimension of the array B.  LDB >= max(1,N).
 *
 ******************************************************************************/
__attribute__((weak))
int coreblas_dsygst(int itype, coreblas_enum_t uplo,
                int n,
                double *A, int lda,
                double *B, int ldb)
{
    int info = LAPACKE_dsygst_work(
        LAPACK_COL_MAJOR,
        itype,
        lapack_const(uplo),
        n, A, lda, B, ldb );
    return info;
}

/******************************************************************************/
void coreblas_kernel_dsygst(int itype, coreblas_enum_t uplo,
                     int n,
                     double *A, int lda,
                     double *B, int ldb)
{

        coreblas_dsygst(itype, uplo,
                    n,
                    A, lda,
                    B, ldb);
}