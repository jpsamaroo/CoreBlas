/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zlantr.c, normal z -> d, Mon Mar 18 06:35:34 2024
 *
 **/

#include <coreblas.h>
#include "coreblas_types.h"
#include "coreblas_internal.h"
#include "core_lapack.h"

#include <math.h>

/******************************************************************************/
__attribute__((weak))
void coreblas_dlantr(coreblas_enum_t norm, coreblas_enum_t uplo, coreblas_enum_t diag,
                 int m, int n,
                 const double *A, int lda,
                 double *work, double *value)
{
    // Due to a bug in LAPACKE < 3.6.1, this function always returns zero.
    // *value = LAPACKE_dlantr_work(LAPACK_COL_MAJOR,
    //                              lapack_const(norm), lapack_const(uplo),
    //                              lapack_const(diag),
    //                              m, n, A, lda, work);

    // Calling LAPACK directly instead.
    char nrm = lapack_const(norm);
    char upl = lapack_const(uplo);
    char dia = lapack_const(diag);
    *value = LAPACK_dlantr(&nrm, &upl, &dia, &m, &n, A, &lda, work);
}

/******************************************************************************/
void coreblas_kernel_dlantr(coreblas_enum_t norm, coreblas_enum_t uplo, coreblas_enum_t diag,
                     int m, int n,
                     const double *A, int lda,
                     double *work, double *value)
{

    coreblas_dlantr(norm, uplo, diag, m, n, A, lda, work, value);
}

/******************************************************************************/
void coreblas_kernel_dlantr_aux(coreblas_enum_t norm, coreblas_enum_t uplo,
                         coreblas_enum_t diag,
                         int m, int n,
                         const double *A, int lda,
                         double *value)
{
    switch (norm) {
    case CoreBlasOneNorm:
        if (uplo == CoreBlasUpper) {
            if (diag == CoreBlasNonUnit) {
                for (int j = 0; j < n; j++) {
                    value[j] = fabs(A[lda*j]);
                    for (int i = 1; i < imin(j+1, m); i++) {
                        value[j] += fabs(A[lda*j+i]);
                    }
                }
            }
            else { // CoreBlasUnit
                int j;
                for (j = 0; j < imin(n, m); j++) {
                    value[j] = 1.0;
                    for (int i = 0; i < j; i++) {
                        value[j] += fabs(A[lda*j+i]);
                    }
                }
                for (; j < n; j++) {
                    value[j] = fabs(A[lda*j]);
                    for (int i = 1; i < m; i++) {
                        value[j] += fabs(A[lda*j+i]);
                    }
                }
            }
        }
        else { // CoreBlasLower
            if (diag == CoreBlasNonUnit) {
                int j;
                for (j = 0; j < imin(n, m); j++) {
                    value[j] = fabs(A[lda*j+j]);
                    for (int i = j+1; i < m; i++) {
                        value[j] += fabs(A[lda*j+i]);
                    }
                }
                for (; j < n; j++)
                    value[j] = 0.0;
            }
            else { // CoreBlasUnit
                int j;
                for (j = 0; j < imin(n, m); j++) {
                    value[j] = 1.0;
                    for (int i = j+1; i < m; i++) {
                        value[j] += fabs(A[lda*j+i]);
                    }
                }
                for (; j < n; j++)
                    value[j] = 0.0;
            }
        }
        break;
    case CoreBlasInfNorm:
        if (uplo == CoreBlasUpper) {
            if (diag == CoreBlasNonUnit) {
                for (int i = 0; i < m; i++)
                    value[i] = 0.0;
                for (int j = 0; j < n; j++) {
                    for (int i = 0; i < imin(j+1, m); i++) {
                        value[i] += fabs(A[lda*j+i]);
                    }
                }
            }
            else { // CoreBlasUnit
                int i;
                for (i = 0; i < imin(m, n); i++)
                    value[i] = 1.0;
                for (; i < m; i++)
                    value[i] = 0.0;
                int j;
                for (j = 0; j < imin(n, m); j++) {
                    for (i = 0; i < j; i++) {
                        value[i] += fabs(A[lda*j+i]);
                    }
                }
                for (; j < n; j++) {
                    for (i = 0; i < m; i++) {
                        value[i] += fabs(A[lda*j+i]);
                    }
                }
            }
        }
        else { // CoreBlasLower
            if (diag == CoreBlasNonUnit) {
                for (int i = 0; i < m; i++)
                    value[i] = 0.0;
                for (int j = 0; j < imin(n, m); j++) {
                    for (int i = j; i < m; i++) {
                        value[i] += fabs(A[lda*j+i]);
                    }
                }
            }
            else { // CoreBlasUnit
                int i;
                for (i = 0; i < imin(m, n); i++)
                    value[i] = 1.0;
                for (; i < m; i++)
                    value[i] = 0.0;
                for (int j = 0; j < imin(n, m); j++) {
                    for (i = j+1; i < m; i++) {
                        value[i] += fabs(A[lda*j+i]);
                    }
                }
            }
        }
        break;
    }
}
