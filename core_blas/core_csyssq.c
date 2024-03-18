/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zsyssq.c, normal z -> c, Mon Mar 18 06:35:37 2024
 *
 **/

#include <coreblas.h>
#include "coreblas_types.h"
#include "core_lapack.h"

#include <math.h>

/******************************************************************************/
__attribute__((weak))
void coreblas_csyssq(coreblas_enum_t uplo,
                 int n,
                 const coreblas_complex32_t *A, int lda,
                 float *scale, float *sumsq)
{
    int ione = 1;
    if (uplo == CoreBlasUpper) {
        for (int j = 1; j < n; j++)
            // TODO: Inline this operation.
            LAPACK_classq(&j, &A[lda*j], &ione, scale, sumsq);
    }
    else { // CoreBlasLower
        for (int j = 0; j < n-1; j++) {
            int len = n-j-1;
            // TODO: Inline this operation.
            LAPACK_classq(&len, &A[lda*j+j+1], &ione, scale, sumsq);
        }
    }
    *sumsq *= 2.0;
    for (int i = 0; i < n; i++) {
        // diagonal is complex, don't ignore complex part
        float absa = cabsf(A[lda*i+i]);
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
}

/******************************************************************************/
void coreblas_kernel_csyssq(coreblas_enum_t uplo,
                     int n,
                     const coreblas_complex32_t *A, int lda,
                     float *scale, float *sumsq)
{

    *scale = 0.0;
    *sumsq = 1.0;
    coreblas_csyssq(uplo, n, A, lda, scale, sumsq);

}

/******************************************************************************/
void coreblas_kernel_csyssq_aux(int m, int n,
                         const float *scale, const float *sumsq,
                         float *value)
{

    float scl = 0.0;
    float sum = 1.0;
    for (int j = 0; j < n; j++) {
        for (int i = j+1; i < n; i++) {
            int idx = m*j+i;
            if (scl < scale[idx]) {
                sum = sumsq[idx] +
                    sum*((scl/scale[idx])*(scl/scale[idx]));
                scl = scale[idx];
            }
            else if (scl > 0.) {
                sum = sum +
                    sumsq[idx]*((scale[idx]/scl)*(scale[idx]/scl));
            }
        }
    }
    sum = 2.0*sum;
    for (int j = 0; j < n; j++) {
        int idx = m*j+j;
        if (scl < scale[idx]) {
            sum = sumsq[idx] + sum*((scl/scale[idx])*(scl/scale[idx]));
            scl = scale[idx];
        }
        else if (scl > 0.) {
            sum = sum + sumsq[idx]*((scale[idx]/scl)*(scale[idx]/scl));
        }
    }
    *value = scl*sqrtf(sum);
}
