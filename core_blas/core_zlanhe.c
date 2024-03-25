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
void coreblas_zlanhe(coreblas_enum_t norm, coreblas_enum_t uplo,
                 int n,
                 const coreblas_complex64_t *A, int lda,
                 double *work, double *value)
{
    #ifdef COREBLAS_USE_64BIT_BLAS
        *value = LAPACKE_zlanhe_work64_(LAPACK_COL_MAJOR,
                                 lapack_const(norm),
                                 lapack_const(uplo),
                                 n, A, lda, work);
    #else
        *value = LAPACKE_zlanhe_work(LAPACK_COL_MAJOR,
                                 lapack_const(norm),
                                 lapack_const(uplo),
                                 n, A, lda, work);
    #endif

}

/******************************************************************************/
void coreblas_zlanhe_aux(coreblas_enum_t norm, coreblas_enum_t uplo,
                         int n,
                         const coreblas_complex64_t *A, int lda,
                         double *value)
{
    switch (norm) {
    case CoreBlasOneNorm:
    case CoreBlasInfNorm:

        if (uplo == CoreBlasUpper) {
            for (int i = 0; i < n; i++)
                value[i] = 0.0;
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < j; i++) {
                    value[i] += cabs(A[lda*j+i]);
                    value[j] += cabs(A[lda*j+i]);
                }
                value[j] += fabs(creal(A[lda*j+j]));
            }
        }
        else { // CoreBlasLower
            for (int i = 0; i < n; i++)
                value[i] = 0.0;
            for (int j = 0; j < n; j++) {
                value[j] += fabs(creal(A[lda*j+j]));
                for (int i = j+1; i < n; i++) {
                    value[i] += cabs(A[lda*j+i]);
                    value[j] += cabs(A[lda*j+i]);
                }
            }
        }
        break;
    }
}