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

#include <math.h>

/***************************************************************************//**
 *
 * @ingroup core_lange
 *
 *  Calculates max, one, infinity or Frobenius norm of a given matrix.
 *
 *******************************************************************************
 *
 * @param[in] norm
 *          - CoreBlasMaxNorm: Max norm
 *          - CoreBlasOneNorm: One norm
 *          - CoreBlasInfNorm: Infinity norm
 *          - CoreBlasFrobeniusNorm: Frobenius norm
 *
 * @param[in] m
 *          The number of rows of the matrix A. m >= 0. When m = 0,
 *          the returned value is set to zero.
 *
 * @param[in] n
 *          The number of columns of the matrix A. n >= 0. When n = 0,
 *          the returned value is set to zero.
 *
 * @param[in] A
 *          The m-by-n matrix A.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,m).
 *
 * @param[in] work
 *          The auxiliary work array.
 *
 * @param[out] value
 *          The specified norm of the given matrix A
 *
 ******************************************************************************/
__attribute__((weak))
void coreblas_zlange(coreblas_enum_t norm, int m, int n,
                 const coreblas_complex64_t *A, int lda,
                 double *work, double *value)
{
    #ifdef COREBLAS_USE_64BIT_BLAS
        *value = LAPACKE_zlange_work64_(LAPACK_COL_MAJOR,
                                 lapack_const(norm),
                                 m, n, A, lda, work);
    #else
        *value = LAPACKE_zlange_work(LAPACK_COL_MAJOR,
                                 lapack_const(norm),
                                 m, n, A, lda, work);
    #endif

}
/******************************************************************************/
void coreblas_zlange_aux(coreblas_enum_t norm, int m, int n,
                         const coreblas_complex64_t *A, int lda,
                         double *value)
{
    switch (norm) {
    case CoreBlasOneNorm:
        for (int j = 0; j < n; j++) {
            value[j] = cabs(A[lda*j]);
            for (int i = 1; i < m; i++) {
                value[j] += cabs(A[lda*j+i]);
            }
        }
        break;
    case CoreBlasInfNorm:
        for (int i = 0; i < m; i++)
            value[i] = 0.0;
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                value[i] += cabs(A[lda*j+i]);
            }
        }

        break;
    }
}