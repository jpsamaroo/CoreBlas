/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_ztrtri.c, normal z -> s, Mon Mar 18 06:35:38 2024
 *
 **/

#include <coreblas.h>
#include "coreblas_types.h"
#include "core_lapack.h"

/***************************************************************************//**
 *
 * @ingroup core_trtri
 *
 *  Computes the inverse of an upper or lower
 *  triangular matrix A.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = CoreBlasUpper: Upper triangle of A is stored;
 *          = CoreBlasLower: Lower triangle of A is stored.
 *
 * @param[in] diag
 *          = CoreBlasNonUnit: A is non-unit triangular;
 *          = CoreBlasUnit:    A is unit triangular.
 *
 * @param[in] n
 *          The order of the matrix A. n >= 0.
 *
 * @param[in,out] A
 *          On entry, the triangular matrix A.  If uplo = 'U', the
 *          leading n-by-n upper triangular part of the array A
 *          contains the upper triangular matrix, and the strictly
 *          lower triangular part of A is not referenced.  If uplo =
 *          'L', the leading n-by-n lower triangular part of the array
 *          A contains the lower triangular matrix, and the strictly
 *          upper triangular part of A is not referenced.  If diag =
 *          'U', the diagonal elements of A are also not referenced and
 *          are assumed to be 1.  On exit, the (triangular) inverse of
 *          the original matrix.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,n).
 *
 * @retval CoreBlasSuccess on successful exit
 * @retval < 0 if -i, the i-th argument had an illegal value
 * @retval > 0 if i, A(i,i) is exactly zero.  The triangular
 *          matrix is singular and its inverse can not be computed.
 *
 ******************************************************************************/
__attribute__((weak))
int coreblas_strtri(coreblas_enum_t uplo, coreblas_enum_t diag,
                 int n,
                 float *A, int lda)
{
    return LAPACKE_strtri_work(LAPACK_COL_MAJOR,
                        lapack_const(uplo), lapack_const(diag),
                        n, A, lda);
}

/******************************************************************************/
void coreblas_kernel_strtri(coreblas_enum_t uplo, coreblas_enum_t diag,
                     int n,
                     float *A, int lda,
                     int iinfo)
{

    int info = coreblas_strtri(uplo, diag,
                           n, A, lda);
    if (info != 0)
        coreblas_error("core_strtri() failed");

}
