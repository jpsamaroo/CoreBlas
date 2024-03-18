/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zgemm.c, normal z -> d, Mon Mar 18 06:35:30 2024
 *
 **/

#include <coreblas.h>
#include "coreblas_types.h"
#include "core_lapack.h"

/***************************************************************************//**
 *
 * @ingroup core_gemm
 *
 *  Performs one of the matrix-matrix operations
 *
 *    \f[ C = \alpha [op( A )\times op( B )] + \beta C, \f]
 *
 *  where op( X ) is one of:
 *    \f[ op( X ) = X,   \f]
 *    \f[ op( X ) = X^T, \f]
 *    \f[ op( X ) = X^T, \f]
 *
 *  alpha and beta are scalars, and A, B and C  are matrices, with op( A )
 *  an m-by-k matrix, op( B ) a k-by-n matrix and C an m-by-n matrix.
 *
 *******************************************************************************
 *
 * @param[in] transa
 *          - CoreBlasNoTrans:   A is not transposed,
 *          - CoreBlasTrans:     A is transposed,
 *          - CoreBlasConjTrans: A is conjugate transposed.
 *
 * @param[in] transb
 *          - CoreBlasNoTrans:   B is not transposed,
 *          - CoreBlasTrans:     B is transposed,
 *          - CoreBlasConjTrans: B is conjugate transposed.
 *
 * @param[in] m
 *          The number of rows of the matrix op( A ) and of the matrix C.
 *          m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix op( B ) and of the matrix C.
 *          n >= 0.
 *
 * @param[in] k
 *          The number of columns of the matrix op( A ) and the number of rows
 *          of the matrix op( B ). k >= 0.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          An lda-by-ka matrix, where ka is k when transa = CoreBlasNoTrans,
 *          and is m otherwise.
 *
 * @param[in] lda
 *          The leading dimension of the array A.
 *          When transa = CoreBlasNoTrans, lda >= max(1,m),
 *          otherwise, lda >= max(1,k).
 *
 * @param[in] B
 *          An ldb-by-kb matrix, where kb is n when transb = CoreBlasNoTrans,
 *          and is k otherwise.
 *
 * @param[in] ldb
 *          The leading dimension of the array B.
 *          When transb = CoreBlasNoTrans, ldb >= max(1,k),
 *          otherwise, ldb >= max(1,n).
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in,out] C
 *          An ldc-by-n matrix. On exit, the array is overwritten by the m-by-n
 *          matrix ( alpha*op( A )*op( B ) + beta*C ).
 *
 * @param[in] ldc
 *          The leading dimension of the array C. ldc >= max(1,m).
 *
 ******************************************************************************/
__attribute__((weak))
void coreblas_dgemm(coreblas_enum_t transa, coreblas_enum_t transb,
                int m, int n, int k,
                double alpha, const double *A, int lda,
                                          const double *B, int ldb,
                double beta,        double *C, int ldc)
{
    cblas_dgemm(CblasColMajor,
                (CBLAS_TRANSPOSE)transa, (CBLAS_TRANSPOSE)transb,
                m, n, k,
                (alpha), A, lda,
                                    B, ldb,
                (beta),  C, ldc);
}

/******************************************************************************/
void coreblas_kernel_dgemm(
    coreblas_enum_t transa, coreblas_enum_t transb,
    int m, int n, int k,
    double alpha, const double *A, int lda,
                              const double *B, int ldb,
    double beta,        double *C, int ldc)
{
    int ak;
    if (transa == CoreBlasNoTrans)
        ak = k;
    else
        ak = m;

    int bk;
    if (transb == CoreBlasNoTrans)
        bk = n;
    else
        bk = k;


    coreblas_dgemm(transa, transb,
               m, n, k,
               alpha, A, lda,
                      B, ldb,
               beta,  C, ldc);
}
