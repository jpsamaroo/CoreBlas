/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zher2k.c, normal z -> c, Mon Mar 18 06:35:33 2024
 *
 **/

#include <coreblas.h>
#include "coreblas_types.h"
#include "core_lapack.h"

#undef REAL
#define COMPLEX
/***************************************************************************//**
 *
 * @ingroup core_her2k
 *
 *  Performs one of the Hermitian rank 2k operations
 *
 *    \f[ C = \alpha A \times B^H + conjg( \alpha ) B \times A^H + \beta C, \f]
 *    or
 *    \f[ C = \alpha A^H \times B + conjg( \alpha ) B^H \times A + \beta C, \f]
 *
 *  where alpha is a complex scalar, beta is a real scalar,
 *  C is an n-by-n Hermitian matrix, and A and B are n-by-k matrices
 *  in the first case and k-by-n matrices in the second case.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - CoreBlasUpper: Upper triangle of C is stored;
 *          - CoreBlasLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          - CoreBlasNoTrans:
 *            \f[ C = \alpha A \times B^H
 *                  + conjg( \alpha ) B \times A^H + \beta C; \f]
 *          - CoreBlasConjTrans:
 *            \f[ C = \alpha A^H \times B
 *                  + conjg( \alpha ) B^H \times A + \beta C. \f]
 *
 * @param[in] n
 *          The order of the matrix C. n >= zero.
 *
 * @param[in] k
 *          If trans = CoreBlasNoTrans, number of columns of the A and B matrices;
 *          if trans = CoreBlasConjTrans, number of rows of the A and B matrices.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          An lda-by-ka matrix.
 *          If trans = CoreBlasNoTrans,   ka = k;
 *          if trans = CoreBlasConjTrans, ka = n.
 *
 * @param[in] lda
 *          The leading dimension of the array A.
 *          If trans = CoreBlasNoTrans,   lda >= max(1, n);
 *          if trans = CoreBlasConjTrans, lda >= max(1, k).
 *
 * @param[in] B
 *          An ldb-by-kb matrix.
 *          If trans = CoreBlasNoTrans,   kb = k;
 *          if trans = CoreBlasConjTrans, kb = n.
 *
 * @param[in] ldb
 *          The leading dimension of the array B.
 *          If trans = CoreBlasNoTrans,   ldb >= max(1, n);
 *          if trans = CoreBlasConjTrans, ldb >= max(1, k).
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in,out] C
 *          An ldc-by-n matrix.
 *          On exit, the uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 * @param[in] ldc
 *          The leading dimension of the array C. ldc >= max(1, n).
 *
 ******************************************************************************/
__attribute__((weak))
void coreblas_cher2k(coreblas_enum_t uplo, coreblas_enum_t trans,
                 int n, int k,
                 coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                                           const coreblas_complex32_t *B, int ldb,
                  float beta,                   coreblas_complex32_t *C, int ldc)
{
    cblas_cher2k(CblasColMajor,
                 (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                 n, k,
                 CBLAS_SADDR(alpha), A, lda,
                                     B, ldb,
                 beta,               C, ldc);
}

/******************************************************************************/
void coreblas_kernel_cher2k(
    coreblas_enum_t uplo, coreblas_enum_t trans,
    int n, int k,
    coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                              const coreblas_complex32_t *B, int ldb,
    float beta,                    coreblas_complex32_t *C, int ldc)
{
    int ak;
    int bk;
    if (trans == CoreBlasNoTrans) {
        ak = k;
        bk = k;
    }
    else {
        ak = n;
        bk = n;
    }

        coreblas_cher2k(uplo, trans,
                    n, k,
                    alpha, A, lda,
                           B, ldb,
                    beta,  C, ldc);
}
