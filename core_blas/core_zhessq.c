/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @precisions normal z -> c
 *
 **/

#include <coreblas.h>
#include "coreblas_types.h"
#include "core_lapack.h"

#include <math.h>

/******************************************************************************/
__attribute__((weak))
void coreblas_zhessq(coreblas_enum_t uplo,
                 int n,
                 const coreblas_complex64_t *A, int lda,
                 double *scale, double *sumsq)
{
    *scale = 0.0;
    *sumsq = 1.0;
    int ione = 1;
    if (uplo == CoreBlasUpper) {
        for (int j = 1; j < n; j++)
            // TODO: Inline this operation.
            #ifdef COREBLAS_USE_64BIT_BLAS
                LAPACK_zlassq64_(&j, &A[lda*j], &ione, scale, sumsq);
            #else
                LAPACK_zlassq(&j, &A[lda*j], &ione, scale, sumsq);
            #endif
            
    }
    else { // CoreBlasLower
        for (int j = 0; j < n-1; j++) {
            int len = n-j-1;
            // TODO: Inline this operation.
            #ifdef COREBLAS_USE_64BIT_BLAS
                LAPACK_zlassq64_(&len, &A[lda*j+j+1], &ione, scale, sumsq);
            #else
                LAPACK_zlassq(&len, &A[lda*j+j+1], &ione, scale, sumsq);
            #endif
        }
    }
    *sumsq *= 2.0;
    for (int i = 0; i < n; i++) {
        // diagonal is real, ignore imaginary part
        if (creal(A[lda*i+i]) != 0.0) { // != propagates nan
            double absa = fabs(creal(A[lda*i+i]));
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