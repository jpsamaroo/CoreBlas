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
 * @ingroup core_trsm
 *
 *  Solves one of the matrix equations
 *
 *    \f[ op( A )\times X  = \alpha B, \f] or
 *    \f[ X \times op( A ) = \alpha B, \f]
 *
 *  where op( A ) is one of:
 *    \f[ op( A ) = A,   \f]
 *    \f[ op( A ) = A^T, \f]
 *    \f[ op( A ) = A^H, \f]
 *
 *  alpha is a scalar, X and B are m-by-n matrices, and
 *  A is a unit or non-unit, upper or lower triangular matrix.
 *  The matrix X overwrites B.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          - CoreBlasLeft:  op(A)*X = B,
 *          - CoreBlasRight: X*op(A) = B.
 *
 * @param[in] uplo
 *          - CoreBlasUpper: A is upper triangular,
 *          - CoreBlasLower: A is lower triangular.
 *
 * @param[in] transa
 *          - CoreBlasNoTrans:   A is not transposed,
 *          - CoreBlasTrans:     A is transposed,
 *          - CoreBlasConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          - CoreBlasNonUnit: A has non-unit diagonal,
 *          - CoreBlasUnit:    A has unit diagonal.
 *
 * @param[in] m
 *          The number of rows of the matrix B. m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix B. n >= 0.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          The lda-by-ka triangular matrix,
 *          where ka = m if side = CoreBlasLeft,
 *            and ka = n if side = CoreBlasRight.
 *          If uplo = CoreBlasUpper, the leading k-by-k upper triangular part
 *          of the array A contains the upper triangular matrix, and the
 *          strictly lower triangular part of A is not referenced.
 *          If uplo = CoreBlasLower, the leading k-by-k lower triangular part
 *          of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced.
 *          If diag = CoreBlasUnit, the diagonal elements of A are also not
 *          referenced and are assumed to be 1.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,k).
 *
 * @param[in,out] B
 *          On entry, the ldb-by-n right hand side matrix B.
 *          On exit, if return value = 0, the ldb-by-n solution matrix X.
 *
 * @param[in] ldb
 *          The leading dimension of the array B. ldb >= max(1,m).
 *
 ******************************************************************************/
__attribute__((weak))
void coreblas_ztrsm(coreblas_enum_t side, coreblas_enum_t uplo,
                coreblas_enum_t transa, coreblas_enum_t diag,
                int m, int n,
                coreblas_complex64_t alpha, const coreblas_complex64_t *A, int lda,
                                                coreblas_complex64_t *B, int ldb)
{
    cblas_ztrsm(CblasColMajor,
                (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
                (CBLAS_TRANSPOSE)transa, (CBLAS_DIAG)diag,
                m, n,
                CBLAS_SADDR(alpha), A, lda,
                                    B, ldb);
}

/******************************************************************************/
void coreblas_kernel_ztrsm(
    coreblas_enum_t side, coreblas_enum_t uplo,
    coreblas_enum_t transa, coreblas_enum_t diag,
    int m, int n,
    coreblas_complex64_t alpha, const coreblas_complex64_t *A, int lda,
                                    coreblas_complex64_t *B, int ldb)
{
    int ak;
    if (side == CoreBlasLeft)
        ak = m;
    else
        ak = n;

    coreblas_ztrsm(side, uplo,
               transa, diag,
               m, n,
               alpha, A, lda,
                B, ldb);
}
