/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zgessq.c, normal z -> d, Mon Mar 18 06:35:32 2024
 *
 **/

#include <math.h>

#include <coreblas.h>
#include "coreblas_types.h"
#include "core_lapack.h"

/******************************************************************************/
__attribute__((weak))
void coreblas_dgessq(int m, int n,
                 const double *A, int lda,
                 double *scale, double *sumsq)
{
    int ione = 1;
    for (int j = 0; j < n; j++) {
        // TODO: Inline this operation.
        LAPACK_dlassq(&m, &A[j*lda], &ione, scale, sumsq);
    }
}

/******************************************************************************/
void coreblas_kernel_dgessq(int m, int n,
                     const double *A, int lda,
                     double *scale, double *sumsq)
{
    *scale = 0.0;
    *sumsq = 1.0;
    coreblas_dgessq(m, n, A, lda, scale, sumsq);
}

/******************************************************************************/
void coreblas_kernel_dgessq_aux(int n,
                         const double *scale, const double *sumsq,
                         double *value)
{

    double scl = 0.0;
    double sum = 1.0;
    for (int i = 0; i < n; i++) {
        if (scl < scale[i]) {
            sum = sumsq[i] + sum*((scl/scale[i])*(scl/scale[i]));
            scl = scale[i];
        }
        else if (scl > 0.) {
            sum = sum + sumsq[i]*(scale[i]/scl)*(scale[i]/scl);
        }
    }
    *value = scl*sqrt(sum);
}
