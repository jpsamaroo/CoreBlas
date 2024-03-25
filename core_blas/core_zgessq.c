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

#include <math.h>

#include <coreblas.h>
#include "coreblas_types.h"
#include "core_lapack.h"

/******************************************************************************/
__attribute__((weak))
void coreblas_zgessq(int m, int n,
                 const coreblas_complex64_t *A, int lda,
                 double *scale, double *sumsq)
{
    *scale = 0.0;
    *sumsq = 1.0;
    int ione = 1;
    for (int j = 0; j < n; j++) {
        // TODO: Inline this operation.
        #ifdef COREBLAS_USE_64BIT_BLAS
            LAPACK_zlassq64_(&m, &A[j*lda], &ione, scale, sumsq);
        #else
            LAPACK_zlassq(&m, &A[j*lda], &ione, scale, sumsq);
        #endif 
    }
}

/******************************************************************************/
void coreblas_zgessq_aux(int n,
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