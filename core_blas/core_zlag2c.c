/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @precisions mixed zc -> ds
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
int coreblas_zlag2c(int m, int n,
                 coreblas_complex64_t *A,  int lda,
                 coreblas_complex32_t *As, int ldas)
{
    int info;
    #ifdef COREBLAS_USE_64BIT_BLAS
        info = LAPACKE_zlag2c_work64_(LAPACK_COL_MAJOR, m, n, A, lda, As, ldas);
    #else
        info = LAPACKE_zlag2c_work(LAPACK_COL_MAJOR, m, n, A, lda, As, ldas);
    #endif

    return info;
}