/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from include/coreblas_zc.h, mixed zc -> ds, Mon Mar 18 06:35:29 2024
 *
 **/
#ifndef COREBLAS_CORE_BLAS_DS_H
#define COREBLAS_CORE_BLAS_DS_H

#include "coreblas_types.h"
#include "coreblas_workspace.h"

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************/
int coreblas_dlag2s(int m, int n,
                 double *A,  int lda,
                 float *As, int ldas);

void coreblas_slag2d(int m, int n,
                 float *As, int ldas,
                 double *A,  int lda);

/******************************************************************************/
void coreblas_kernel_dlag2s(int m, int n,
                     double *A,  int lda,
                     float *As, int ldas);

void coreblas_kernel_slag2d(int m, int n,
                     float *As, int ldas,
                     double *A,  int lda);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // COREBLAS_CORE_BLAS_DS_H
