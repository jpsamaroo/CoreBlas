/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zlansy.c, normal z -> d, Mon Mar 18 06:35:34 2024
 *
 **/

#include <coreblas.h>
#include "coreblas_types.h"
#include "core_lapack.h"

#include <math.h>

/******************************************************************************/
__attribute__((weak))
void coreblas_dlansy(coreblas_enum_t norm, coreblas_enum_t uplo,
                 int n,
                 const double *A, int lda,
                 double *work, double *value)
{
    *value = LAPACKE_dlansy_work(LAPACK_COL_MAJOR,
                                 lapack_const(norm),
                                 lapack_const(uplo),
                                 n, A, lda, work);
}

/******************************************************************************/
void coreblas_kernel_dlansy(coreblas_enum_t norm, coreblas_enum_t uplo,
                     int n,
                     const double *A, int lda,
                     double *work, double *value)
{

    coreblas_dlansy(norm, uplo, n, A, lda, work, value);
}

/******************************************************************************/
void coreblas_kernel_dlansy_aux(coreblas_enum_t norm, coreblas_enum_t uplo,
                         int n,
                         const double *A, int lda,
                         double *value)
{
    switch (norm) {
    case CoreBlasOneNorm:
    case CoreBlasInfNorm:
        if (uplo == CoreBlasUpper) {
            for (int i = 0; i < n; i++)
                value[i] = 0.0;
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < j; i++) {
                    value[i] += fabs(A[lda*j+i]);
                    value[j] += fabs(A[lda*j+i]);
                }
                value[j] += fabs(A[lda*j+j]);
            }
        }
        else { // CoreBlasLower
            for (int i = 0; i < n; i++)
                value[i] = 0.0;
            for (int j = 0; j < n; j++) {
                value[j] += fabs(A[lda*j+j]);
                for (int i = j+1; i < n; i++) {
                    value[i] += fabs(A[lda*j+i]);
                    value[j] += fabs(A[lda*j+i]);
                }
            }
        }
        break;
    }
}
