/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @precisions normal z -> s d c
 *
 **/
#ifndef COREBLAS_CORE_BLAS_Z_H
#define COREBLAS_CORE_BLAS_Z_H

#include "coreblas_types.h"
#include "coreblas_workspace.h"

#ifdef __cplusplus
extern "C" {
#endif

#define COMPLEX

/******************************************************************************/
#ifdef COMPLEX
double coreblas_dcabs1(coreblas_complex64_t alpha);
#endif

void coreblas_zgbtype1cb(coreblas_enum_t uplo, int n, int nb,
                      coreblas_complex64_t *A, int lda,
                      coreblas_complex64_t *VQ, coreblas_complex64_t *TAUQ,
                      coreblas_complex64_t *VP, coreblas_complex64_t *TAUP,
                      int st, int ed, int sweep, int Vblksiz, int WANTZ,
                      coreblas_complex64_t *work);
    
void coreblas_zgbtype2cb(coreblas_enum_t uplo, int n, int nb,
                      coreblas_complex64_t *A, int lda,
                      coreblas_complex64_t *VQ, coreblas_complex64_t *TAUQ,
                      coreblas_complex64_t *VP, coreblas_complex64_t *TAUP,
                      int st, int ed, int sweep, int Vblksiz, int WANTZ,
                      coreblas_complex64_t *work);
    
void coreblas_zgbtype3cb(coreblas_enum_t uplo, int n, int nb,
                      coreblas_complex64_t *A, int lda,
                      coreblas_complex64_t *VQ, coreblas_complex64_t *TAUQ,
                      coreblas_complex64_t *VP, coreblas_complex64_t *TAUP,
                      int st, int ed, int sweep, int Vblksiz, int WANTZ,
                      coreblas_complex64_t *work);
    
int coreblas_zgeadd(coreblas_enum_t transa,
                int m, int n,
                coreblas_complex64_t alpha, const coreblas_complex64_t *A, int lda,
                coreblas_complex64_t beta,        coreblas_complex64_t *B, int ldb);

int coreblas_zgelqt(int m, int n, int ib,
                coreblas_complex64_t *A, int lda,
                coreblas_complex64_t *T, int ldt,
                coreblas_complex64_t *tau,
                coreblas_complex64_t *work);

void coreblas_zgemm(coreblas_enum_t transa, coreblas_enum_t transb,
                int m, int n, int k,
                coreblas_complex64_t alpha, const coreblas_complex64_t *A, int lda,
                                          const coreblas_complex64_t *B, int ldb,
                coreblas_complex64_t beta,        coreblas_complex64_t *C, int ldc);

int coreblas_zgeqrt(int m, int n, int ib,
                coreblas_complex64_t *A, int lda,
                coreblas_complex64_t *T, int ldt,
                coreblas_complex64_t *tau,
                coreblas_complex64_t *work);

void coreblas_zgessq(int m, int n,
                 const coreblas_complex64_t *A, int lda,
                 double *scale, double *sumsq);

//void coreblas_zgetrf(coreblas_desc_t A, int *ipiv, int ib, int rank, int size,
//                 volatile int *max_idx, volatile coreblas_complex64_t *max_val,
//                 volatile int *info);

int coreblas_zhegst(int itype, coreblas_enum_t uplo,
                int n,
                coreblas_complex64_t *A, int lda,
                coreblas_complex64_t *B, int ldb);

void coreblas_zhemm(coreblas_enum_t side, coreblas_enum_t uplo,
                int m, int n,
                coreblas_complex64_t alpha, const coreblas_complex64_t *A, int lda,
                                          const coreblas_complex64_t *B, int ldb,
                coreblas_complex64_t beta,        coreblas_complex64_t *C, int ldc);

void coreblas_zher2k(coreblas_enum_t uplo, coreblas_enum_t trans,
                 int n, int k,
                 coreblas_complex64_t alpha, const coreblas_complex64_t *A, int lda,
                                           const coreblas_complex64_t *B, int ldb,
                 double beta,                    coreblas_complex64_t *C, int ldc);

void coreblas_zherk(coreblas_enum_t uplo, coreblas_enum_t trans,
                int n, int k,
                double alpha, const coreblas_complex64_t *A, int lda,
                double beta,        coreblas_complex64_t *C, int ldc);

void coreblas_zhessq(coreblas_enum_t uplo,
                 int n,
                 const coreblas_complex64_t *A, int lda,
                 double *scale, double *sumsq);

void coreblas_zsyssq(coreblas_enum_t uplo,
                 int n,
                 const coreblas_complex64_t *A, int lda,
                 double *scale, double *sumsq);

void coreblas_zlacpy(coreblas_enum_t uplo, coreblas_enum_t transa,
                 int m, int n,
                 const coreblas_complex64_t *A, int lda,
                       coreblas_complex64_t *B, int ldb);

void coreblas_zlacpy_lapack2tile_band(coreblas_enum_t uplo,
                                  int it, int jt,
                                  int m, int n, int nb, int kl, int ku,
                                  const coreblas_complex64_t *A, int lda,
                                        coreblas_complex64_t *B, int ldb);

void coreblas_zlacpy_tile2lapack_band(coreblas_enum_t uplo,
                                  int it, int jt,
                                  int m, int n, int nb, int kl, int ku,
                                  const coreblas_complex64_t *B, int ldb,
                                        coreblas_complex64_t *A, int lda);

void coreblas_zlange(coreblas_enum_t norm,
                 int m, int n,
                 const coreblas_complex64_t *A, int lda,
                 double *work, double *result);

void coreblas_zlanhe(coreblas_enum_t norm, coreblas_enum_t uplo,
                 int n,
                 const coreblas_complex64_t *A, int lda,
                 double *work, double *value);

void coreblas_zlansy(coreblas_enum_t norm, coreblas_enum_t uplo,
                 int n,
                 const coreblas_complex64_t *A, int lda,
                 double *work, double *value);

void coreblas_zlantr(coreblas_enum_t norm, coreblas_enum_t uplo, coreblas_enum_t diag,
                 int m, int n,
                 const coreblas_complex64_t *A, int lda,
                 double *work, double *value);

int coreblas_zlarfb_gemm(coreblas_enum_t side, coreblas_enum_t trans, int direct, int storev,
                     int M, int N, int K,
                     const coreblas_complex64_t *V, int LDV,
                     const coreblas_complex64_t *T, int LDT,
                     coreblas_complex64_t *C, int LDC,
                     coreblas_complex64_t *WORK, int LDWORK);

void coreblas_zlascl(coreblas_enum_t uplo,
                 double cfrom, double cto,
                 int m, int n,
                 coreblas_complex64_t *A, int lda);

void coreblas_zlaset(coreblas_enum_t uplo,
                 int m, int n,
                 coreblas_complex64_t alpha, coreblas_complex64_t beta,
                 coreblas_complex64_t *A, int lda);
/*
void coreblas_zgeswp(coreblas_enum_t colrow,
                 coreblas_desc_t A, int k1, int k2, const int *ipiv, int incx);

void coreblas_zheswp(int rank, int num_threads,
                 int uplo, coreblas_desc_t A, int k1, int k2, const int *ipiv,
                 int incx);
*/
int coreblas_zlauum(coreblas_enum_t uplo,
                int n,
                coreblas_complex64_t *A, int lda);

int coreblas_zpamm(coreblas_enum_t op, coreblas_enum_t side, coreblas_enum_t storev,
               int m, int n, int k, int l,
               const coreblas_complex64_t *A1, int lda1,
                     coreblas_complex64_t *A2, int lda2,
               const coreblas_complex64_t *V,  int ldv,
                     coreblas_complex64_t *W,  int ldw);

int coreblas_zparfb(coreblas_enum_t side, coreblas_enum_t trans, coreblas_enum_t direct,
                coreblas_enum_t storev,
                int m1, int n1, int m2, int n2, int k, int l,
                      coreblas_complex64_t *A1,   int lda1,
                      coreblas_complex64_t *A2,   int lda2,
                const coreblas_complex64_t *V,    int ldv,
                const coreblas_complex64_t *T,    int ldt,
                      coreblas_complex64_t *work, int ldwork);

int coreblas_zpemv(coreblas_enum_t trans, int storev,
               int m, int n, int l,
               coreblas_complex64_t alpha,
               const coreblas_complex64_t *A, int lda,
               const coreblas_complex64_t *X, int incx,
               coreblas_complex64_t beta,
               coreblas_complex64_t *Y, int incy,
               coreblas_complex64_t *work);

int coreblas_zpotrf(coreblas_enum_t uplo,
                int n,
                coreblas_complex64_t *A, int lda);

void coreblas_zsymm(coreblas_enum_t side, coreblas_enum_t uplo,
                int m, int n,
                coreblas_complex64_t alpha, const coreblas_complex64_t *A, int lda,
                                          const coreblas_complex64_t *B, int ldb,
                coreblas_complex64_t beta,        coreblas_complex64_t *C, int ldc);

void coreblas_zsyr2k(
    coreblas_enum_t uplo, coreblas_enum_t trans,
    int n, int k,
    coreblas_complex64_t alpha, const coreblas_complex64_t *A, int lda,
                              const coreblas_complex64_t *B, int ldb,
    coreblas_complex64_t beta,        coreblas_complex64_t *C, int ldc);

void coreblas_zsyrk(coreblas_enum_t uplo, coreblas_enum_t trans,
                int n, int k,
                coreblas_complex64_t alpha, const coreblas_complex64_t *A, int lda,
                coreblas_complex64_t beta,        coreblas_complex64_t *C, int ldc);

int coreblas_ztradd(coreblas_enum_t uplo, coreblas_enum_t transa,
                int m, int n,
                coreblas_complex64_t alpha, const coreblas_complex64_t *A, int lda,
                coreblas_complex64_t beta,        coreblas_complex64_t *B, int ldb);

void coreblas_ztrmm(coreblas_enum_t side, coreblas_enum_t uplo,
                coreblas_enum_t transa, coreblas_enum_t diag,
                int m, int n,
                coreblas_complex64_t alpha, const coreblas_complex64_t *A, int lda,
                                                coreblas_complex64_t *B, int ldb);

void coreblas_ztrsm(coreblas_enum_t side, coreblas_enum_t uplo,
                coreblas_enum_t transa, coreblas_enum_t diag,
                int m, int n,
                coreblas_complex64_t alpha, const coreblas_complex64_t *A, int lda,
                                                coreblas_complex64_t *B, int ldb);

void coreblas_ztrssq(coreblas_enum_t uplo, coreblas_enum_t diag,
                 int m, int n,
                 const coreblas_complex64_t *A, int lda,
                 double *scale, double *sumsq);

int coreblas_ztrtri(coreblas_enum_t uplo, coreblas_enum_t diag,
                int n,
                coreblas_complex64_t *A, int lda);

int coreblas_ztslqt(int m, int n, int ib,
                coreblas_complex64_t *A1, int lda1,
                coreblas_complex64_t *A2, int lda2,
                coreblas_complex64_t *T,  int ldt,
                coreblas_complex64_t *tau,
                coreblas_complex64_t *work);

int coreblas_ztsmlq(coreblas_enum_t side, coreblas_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      coreblas_complex64_t *A1,   int lda1,
                      coreblas_complex64_t *A2,   int lda2,
                const coreblas_complex64_t *V,    int ldv,
                const coreblas_complex64_t *T,    int ldt,
                      coreblas_complex64_t *work, int ldwork);

int coreblas_ztsmqr(coreblas_enum_t side, coreblas_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      coreblas_complex64_t *A1,   int lda1,
                      coreblas_complex64_t *A2,   int lda2,
                const coreblas_complex64_t *V,    int ldv,
                const coreblas_complex64_t *T,    int ldt,
                      coreblas_complex64_t *work, int ldwork);

int coreblas_ztsqrt(int m, int n, int ib,
                coreblas_complex64_t *A1, int lda1,
                coreblas_complex64_t *A2, int lda2,
                coreblas_complex64_t *T,  int ldt,
                coreblas_complex64_t *tau,
                coreblas_complex64_t *work);

int coreblas_zttlqt(int m, int n, int ib,
                coreblas_complex64_t *A1, int lda1,
                coreblas_complex64_t *A2, int lda2,
                coreblas_complex64_t *T,  int ldt,
                coreblas_complex64_t *tau,
                coreblas_complex64_t *work);

int coreblas_zttmlq(coreblas_enum_t side, coreblas_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      coreblas_complex64_t *A1,   int lda1,
                      coreblas_complex64_t *A2,   int lda2,
                const coreblas_complex64_t *V,    int ldv,
                const coreblas_complex64_t *T,    int ldt,
                      coreblas_complex64_t *work, int ldwork);

int coreblas_zttmqr(coreblas_enum_t side, coreblas_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      coreblas_complex64_t *A1,   int lda1,
                      coreblas_complex64_t *A2,   int lda2,
                const coreblas_complex64_t *V,    int ldv,
                const coreblas_complex64_t *T,    int ldt,
                      coreblas_complex64_t *work, int ldwork);

int coreblas_zttqrt(int m, int n, int ib,
                coreblas_complex64_t *A1, int lda1,
                coreblas_complex64_t *A2, int lda2,
                coreblas_complex64_t *T,  int ldt,
                coreblas_complex64_t *tau,
                coreblas_complex64_t *work);

int coreblas_zunmlq(coreblas_enum_t side, coreblas_enum_t trans,
                int m, int n, int k, int ib,
                const coreblas_complex64_t *A,    int lda,
                const coreblas_complex64_t *T,    int ldt,
                      coreblas_complex64_t *C,    int ldc,
                      coreblas_complex64_t *work, int ldwork);

int coreblas_zunmqr(coreblas_enum_t side, coreblas_enum_t trans,
                int m, int n, int k, int ib,
                const coreblas_complex64_t *A,    int lda,
                const coreblas_complex64_t *T,    int ldt,
                      coreblas_complex64_t *C,    int ldc,
                      coreblas_complex64_t *work, int ldwork);

#undef COMPLEX

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // COREBLAS_CORE_BLAS_Z_H
