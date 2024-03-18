/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zlange.c, normal z -> d, Mon Mar 18 06:35:33 2024
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
void coreblas_dlange(coreblas_enum_t norm, int m, int n,
                 const double *A, int lda,
                 double *work, double *value)
{
    *value = LAPACKE_dlange_work(LAPACK_COL_MAJOR,
                                 lapack_const(norm),
                                 m, n, A, lda, work);
}

/******************************************************************************/
void coreblas_kernel_dlange(coreblas_enum_t norm, int m, int n,
                     const double *A, int lda,
                     double *work, double *value)
{

    coreblas_dlange(norm, m, n, A, lda, work, value);
}

/******************************************************************************/
void coreblas_kernel_dlange_aux(coreblas_enum_t norm, int m, int n,
                         const double *A, int lda,
                         double *value)
{
    switch (norm) {
    case CoreBlasOneNorm:
        for (int j = 0; j < n; j++) {
            value[j] = fabs(A[lda*j]);
            for (int i = 1; i < m; i++) {
                value[j] += fabs(A[lda*j+i]);
            }
        }
        break;
    case CoreBlasInfNorm:
        for (int i = 0; i < m; i++)
            value[i] = 0.0;
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                value[i] += fabs(A[lda*j+i]);
            }
        }

        break;
    }
}
