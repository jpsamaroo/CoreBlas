/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zlag2c.c, mixed zc -> ds, Mon Mar 18 06:35:42 2024
 *
 **/

#include <coreblas.h>
#include "core_lapack.h"
#include "coreblas_types.h"

/***************************************************************************//**
 *
 * @ingroup core_lag2
 *
 *  Converts m-by-n matrix A from double complex to single complex precision.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows of the matrix A.
 *          m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix A.
 *          n >= 0.
 *
 * @param[in] A
 *          The lda-by-n matrix in double complex precision to convert.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A.
 *          lda >= max(1,m).
 *
 * @param[out] As
 *          On exit, the converted ldas-by-n matrix in single complex precision.
 *
 * @param[in] ldas
 *          The leading dimension of the matrix As.
 *          ldas >= max(1,m).
 *
 ******************************************************************************/
__attribute__((weak))
int coreblas_dlag2s(int m, int n,
                 double *A,  int lda,
                 float *As, int ldas)
{
    int info;
    info = LAPACKE_dlag2s_work(LAPACK_COL_MAJOR, m, n, A, lda, As, ldas);
    return info;
}

/******************************************************************************/
void coreblas_kernel_dlag2s(int m, int n,
                     double *A,  int lda,
                     float *As, int ldas)
{

    int info;
    info = coreblas_dlag2s(m, n, A, lda, As, ldas);
    if (info != 0) {
        // Value will be 1, so it doesn't matter which tile sets status.
        coreblas_error("core_dgeadd() failed");                    
    }
}
