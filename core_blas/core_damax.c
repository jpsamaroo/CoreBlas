/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_dzamax.c, normal z -> d, Mon Mar 18 06:35:32 2024
 *
 **/

#include <coreblas.h>
#include "coreblas_types.h"
#include "core_lapack.h"

#include <math.h>

/******************************************************************************/
void coreblas_kernel_damax(int colrow, int m, int n,
                     const double *A, int lda,
                     double *values)
{
    switch (colrow) {
    case CoreBlasColumnwise:

        for (int j = 0; j < n; j++) {
            values[j] = fabs(A[lda*j]);
            for (int i = 1; i < m; i++) {
                double tmp = fabs(A[lda*j+i]);
                if (tmp > values[j])
                    values[j] = tmp;
            }
        }
        break;
    case CoreBlasRowwise:
 
        for (int i = 0; i < m; i++)
            values[i] = fabs(A[i]);
        for (int j = 1; j < n; j++) {
            for (int i = 0; i < m; i++) {
                double tmp = fabs(A[lda*j+i]);
                if (tmp > values[i])
                    values[i] = tmp;
            }
        }
        break;
    }
}
