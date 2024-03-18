/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_ztrssq.c, normal z -> s, Mon Mar 18 06:35:38 2024
 *
 **/

#include <coreblas.h>
#include "coreblas_types.h"
#include "coreblas_internal.h"
#include "core_lapack.h"

#include <math.h>

/******************************************************************************/
// This computation also shows up in coreblas_ssyssq() and can be factored out.
// LAPACK does real and imag components separately in slassq.
static inline void ssq(float value, float *scale, float *sumsq)
{
    float absa = fabsf(value);
    if (absa != 0.0) { // != propagates nan
        if (*scale < absa) {
            *sumsq = 1.0 + *sumsq*((*scale/absa)*(*scale/absa));
            *scale = absa;
        }
        else {
            *sumsq = *sumsq + ((absa/(*scale))*(absa/(*scale)));
        }
    }
}

/******************************************************************************/
__attribute__((weak))
void coreblas_strssq(coreblas_enum_t uplo, coreblas_enum_t diag,
                 int m, int n,
                 const float *A, int lda,
                 float *scale, float *sumsq)
{
    if (uplo == CoreBlasUpper) {
        if (diag == CoreBlasNonUnit) {
            for (int j = 0; j < n; j++) {
                ssq(A[lda*j], scale, sumsq);
                for (int i = 1; i < imin(j+1, m); i++) {
                    ssq(A[lda*j+i], scale, sumsq);
                }
            }
        }
        else { // CoreBlasUnit
            int j;
            for (j = 0; j < imin(n, m); j++) {
                ssq(1.0, scale, sumsq);
                for (int i = 0; i < j; i++) {
                    ssq(A[lda*j+i], scale, sumsq);
                }
            }
            for (; j < n; j++) {
                ssq(A[lda*j], scale, sumsq);
                for (int i = 1; i < m; i++) {
                    ssq(A[lda*j+i], scale, sumsq);
                }
            }
        }
    }
    else { // CoreBlasLower
        if (diag == CoreBlasNonUnit) {
            for (int j = 0; j < imin(n, m); j++) {
                ssq(A[lda*j+j], scale, sumsq);
                for (int i = j+1; i < m; i++) {
                    ssq(A[lda*j+i], scale, sumsq);
                }
            }
        }
        else { // CoreBlasUnit
            for (int j = 0; j < imin(n, m); j++) {
                ssq(1.0, scale, sumsq);
                for (int i = j+1; i < m; i++) {
                    ssq(A[lda*j+i], scale, sumsq);
                }
            }
        }
    }
}

/******************************************************************************/
void coreblas_kernel_strssq(coreblas_enum_t uplo, coreblas_enum_t diag,
                     int m, int n,
                     const float *A, int lda,
                     float *scale, float *sumsq)
{

    *scale = 0.0;
    *sumsq = 1.0;
    coreblas_strssq(uplo, diag, m, n, A, lda, scale, sumsq);

}
