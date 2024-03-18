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
#include "coreblas_internal.h"
#include "core_lapack.h"

#include <math.h>

/******************************************************************************/
// This computation also shows up in coreblas_zsyssq() and can be factored out.
// LAPACK does real and imag components separately in zlassq.
static inline void ssq(coreblas_complex64_t value, double *scale, double *sumsq)
{
    double absa = cabs(value);
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
void coreblas_ztrssq(coreblas_enum_t uplo, coreblas_enum_t diag,
                 int m, int n,
                 const coreblas_complex64_t *A, int lda,
                 double *scale, double *sumsq)
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
void coreblas_kernel_ztrssq(coreblas_enum_t uplo, coreblas_enum_t diag,
                     int m, int n,
                     const coreblas_complex64_t *A, int lda,
                     double *scale, double *sumsq)
{

    *scale = 0.0;
    *sumsq = 1.0;
    coreblas_ztrssq(uplo, diag, m, n, A, lda, scale, sumsq);

}
