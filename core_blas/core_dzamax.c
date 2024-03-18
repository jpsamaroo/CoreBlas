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
#include "core_lapack.h"

#include <math.h>

/******************************************************************************/
void coreblas_kernel_dzamax(int colrow, int m, int n,
                     const coreblas_complex64_t *A, int lda,
                     double *values)
{
    switch (colrow) {
    case CoreBlasColumnwise:

        for (int j = 0; j < n; j++) {
            values[j] = coreblas_dcabs1(A[lda*j]);
            for (int i = 1; i < m; i++) {
                double tmp = coreblas_dcabs1(A[lda*j+i]);
                if (tmp > values[j])
                    values[j] = tmp;
            }
        }
        break;
    case CoreBlasRowwise:
 
        for (int i = 0; i < m; i++)
            values[i] = coreblas_dcabs1(A[i]);
        for (int j = 1; j < n; j++) {
            for (int i = 0; i < m; i++) {
                double tmp = coreblas_dcabs1(A[lda*j+i]);
                if (tmp > values[i])
                    values[i] = tmp;
            }
        }
        break;
    }
}
