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

/***************************************************************************//**
 *
 * @ingroup core_syr2k
 *
 *  Performs one of the symmetric rank 2k operations
 *
 *    \f[ C = \alpha A \times B^T + \alpha B \times A^T + \beta C, \f]
 *    or
 *    \f[ C = \alpha A^T \times B + \alpha B^T \times A + \beta C, \f]
 *
 *  where alpha and beta are scalars,
 *  C is an n-by-n symmetric matrix, and A and B are n-by-k matrices
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
 *            \f[ C = \alpha A \times B^T + \alpha B \times A^T + \beta C; \f]
 *          - CoreBlasTrans:
 *            \f[ C = \alpha A^T \times B + \alpha B^T \times A + \beta C. \f]
 *
 * @param[in] n
 *          The order of the matrix C. n >= zero.
 *
 * @param[in] k
 *          If trans = CoreBlasNoTrans, number of columns of the A and B matrices;
 *          if trans = CoreBlasTrans, number of rows of the A and B matrices.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          An lda-by-ka matrix.
 *          If trans = CoreBlasNoTrans, ka = k;
 *          if trans = CoreBlasTrans,   ka = n.
 *
 * @param[in] lda
 *          The leading dimension of the array A.
 *          If trans = CoreBlasNoTrans, lda >= max(1, n);
 *          if trans = CoreBlasTrans,   lda >= max(1, k).
 *
 * @param[in] B
 *          An ldb-by-kb matrix.
 *          If trans = CoreBlasNoTrans, kb = k;
 *          if trans = CoreBlasTrans,   kb = n.
 *
 * @param[in] ldb
 *          The leading dimension of the array B.
 *          If trans = CoreBlasNoTrans, ldb >= max(1, n);
 *          if trans = CoreBlasTrans,   ldb >= max(1, k).
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
void coreblas_zsyr2k(coreblas_enum_t uplo, coreblas_enum_t trans,
                 int n, int k,
                 coreblas_complex64_t alpha, const coreblas_complex64_t *A, int lda,
                                           const coreblas_complex64_t *B, int ldb,
                 coreblas_complex64_t beta,        coreblas_complex64_t *C, int ldc)
{
    #ifdef COREBLAS_USE_64BIT_BLAS
        cblas_zsyr2k64_(CblasColMajor,
                 (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                 n, k,
                 CBLAS_SADDR(alpha), A, lda,
                                     B, ldb,
                 CBLAS_SADDR(beta),  C, ldc);
    #else
        cblas_zsyr2k(CblasColMajor,
                 (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                 n, k,
                 CBLAS_SADDR(alpha), A, lda,
                                     B, ldb,
                 CBLAS_SADDR(beta),  C, ldc);
    #endif

}
