/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zgessq.c, normal z -> s, Mon Mar 18 06:35:32 2024
 *
 **/

#include <math.h>

#include <coreblas.h>
#include "coreblas_types.h"
#include "core_lapack.h"

/******************************************************************************/
__attribute__((weak))
void coreblas_sgessq(int m, int n,
                 const float *A, int lda,
                 float *scale, float *sumsq)
{
    int ione = 1;
    for (int j = 0; j < n; j++) {
        // TODO: Inline this operation.
        LAPACK_slassq(&m, &A[j*lda], &ione, scale, sumsq);
    }
}

/******************************************************************************/
void coreblas_kernel_sgessq(int m, int n,
                     const float *A, int lda,
                     float *scale, float *sumsq)
{
    *scale = 0.0;
    *sumsq = 1.0;
    coreblas_sgessq(m, n, A, lda, scale, sumsq);
}

/******************************************************************************/
void coreblas_kernel_sgessq_aux(int n,
                         const float *scale, const float *sumsq,
                         float *value)
{

    float scl = 0.0;
    float sum = 1.0;
    for (int i = 0; i < n; i++) {
        if (scl < scale[i]) {
            sum = sumsq[i] + sum*((scl/scale[i])*(scl/scale[i]));
            scl = scale[i];
        }
        else if (scl > 0.) {
            sum = sum + sumsq[i]*(scale[i]/scl)*(scale[i]/scl);
        }
    }
    *value = scl*sqrtf(sum);
}
