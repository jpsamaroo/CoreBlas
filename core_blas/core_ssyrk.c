/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zsyrk.c, normal z -> s, Mon Mar 18 06:35:37 2024
 *
 **/

#include <coreblas.h>
#include "coreblas_types.h"
#include "core_lapack.h"

/***************************************************************************//**
 *
 * @ingroup core_syrk
 *
 *  Performs one of the symmetric rank k operations
 *
 *    \f[ C = \alpha A \times A^T + \beta C, \f]
 *    or
 *    \f[ C = \alpha A^T \times A + \beta C, \f]
 *
 *  where alpha and beta are scalars, C is an n-by-n symmetric
 *  matrix, and A is an n-by-k matrix in the first case and a k-by-n
 *  matrix in the second case.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - CoreBlasUpper: Upper triangle of C is stored;
 *          - CoreBlasLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          - CoreBlasNoTrans: \f[ C = \alpha A \times A^T + \beta C; \f]
 *          - CoreBlasTrans:   \f[ C = \alpha A^T \times A + \beta C. \f]
 *
 * @param[in] n
 *          The order of the matrix C. n >= 0.
 *
 * @param[in] k
 *          If trans = CoreBlasNoTrans, number of columns of the A matrix;
 *          if trans = CoreBlasTrans, number of rows of the A matrix.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          A is an lda-by-ka matrix.
 *          If trans = CoreBlasNoTrans, ka = k;
 *          if trans = CoreBlasTrans,   ka = n.
 *
 * @param[in] lda
 *          The leading dimension of the array A.
 *          If trans = CoreBlasNoTrans, lda >= max(1, n);
 *          if trans = CoreBlasTrans,   lda >= max(1, k).
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in,out] C
 *          C is an ldc-by-n matrix.
 *          On exit, the uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 * @param[in] ldc
 *          The leading dimension of the array C. ldc >= max(1, n).
 *
 ******************************************************************************/
__attribute__((weak))
void coreblas_ssyrk(coreblas_enum_t uplo, coreblas_enum_t trans,
                int n, int k,
                float alpha, const float *A, int lda,
                float beta,        float *C, int ldc)
{
    cblas_ssyrk(CblasColMajor,
                (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                n, k,
                (alpha), A, lda,
                (beta),  C, ldc);
}

/******************************************************************************/
void coreblas_kernel_ssyrk(
    coreblas_enum_t uplo, coreblas_enum_t trans,
    int n, int k,
    float alpha, const float *A, int lda,
    float beta,        float *C, int ldc)
{
    int ak;
    if (trans == CoreBlasNoTrans)
        ak = k;
    else
        ak = n;


    coreblas_ssyrk(uplo, trans,
               n, k,
               alpha, A, lda,
               beta,  C, ldc);

}
