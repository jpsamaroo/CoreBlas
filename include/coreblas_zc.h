/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @precisions mixed zc -> ds
 *
 **/
#ifndef COREBLAS_CORE_BLAS_ZC_H
#define COREBLAS_CORE_BLAS_ZC_H

#include "coreblas_types.h"
#include "coreblas_workspace.h"

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************/
int coreblas_zlag2c(int m, int n,
                 coreblas_complex64_t *A,  int lda,
                 coreblas_complex32_t *As, int ldas);

void coreblas_clag2z(int m, int n,
                 coreblas_complex32_t *As, int ldas,
                 coreblas_complex64_t *A,  int lda);

/******************************************************************************/
void coreblas_kernel_zlag2c(int m, int n,
                     coreblas_complex64_t *A,  int lda,
                     coreblas_complex32_t *As, int ldas);

void coreblas_kernel_clag2z(int m, int n,
                     coreblas_complex32_t *As, int ldas,
                     coreblas_complex64_t *A,  int lda);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // COREBLAS_CORE_BLAS_ZC_H
