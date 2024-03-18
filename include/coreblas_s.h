/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from include/coreblas_z.h, normal z -> s, Mon Mar 18 06:35:29 2024
 *
 **/
#ifndef COREBLAS_CORE_BLAS_S_H
#define COREBLAS_CORE_BLAS_S_H

#include "coreblas_types.h"
#include "coreblas_workspace.h"

#ifdef __cplusplus
extern "C" {
#endif

#define REAL

/******************************************************************************/
#ifdef COMPLEX
float fabsf(float alpha);
#endif

void coreblas_sgbtype1cb(coreblas_enum_t uplo, int n, int nb,
                      float *A, int lda,
                      float *VQ, float *TAUQ,
                      float *VP, float *TAUP,
                      int st, int ed, int sweep, int Vblksiz, int WANTZ,
                      float *work);
    
void coreblas_sgbtype2cb(coreblas_enum_t uplo, int n, int nb,
                      float *A, int lda,
                      float *VQ, float *TAUQ,
                      float *VP, float *TAUP,
                      int st, int ed, int sweep, int Vblksiz, int WANTZ,
                      float *work);
    
void coreblas_sgbtype3cb(coreblas_enum_t uplo, int n, int nb,
                      float *A, int lda,
                      float *VQ, float *TAUQ,
                      float *VP, float *TAUP,
                      int st, int ed, int sweep, int Vblksiz, int WANTZ,
                      float *work);
    
int coreblas_sgeadd(coreblas_enum_t transa,
                int m, int n,
                float alpha, const float *A, int lda,
                float beta,        float *B, int ldb);

int coreblas_sgelqt(int m, int n, int ib,
                float *A, int lda,
                float *T, int ldt,
                float *tau,
                float *work);

void coreblas_sgemm(coreblas_enum_t transa, coreblas_enum_t transb,
                int m, int n, int k,
                float alpha, const float *A, int lda,
                                          const float *B, int ldb,
                float beta,        float *C, int ldc);

int coreblas_sgeqrt(int m, int n, int ib,
                float *A, int lda,
                float *T, int ldt,
                float *tau,
                float *work);

void coreblas_sgessq(int m, int n,
                 const float *A, int lda,
                 float *scale, float *sumsq);

//void coreblas_sgetrf(coreblas_desc_t A, int *ipiv, int ib, int rank, int size,
//                 volatile int *max_idx, volatile float *max_val,
//                 volatile int *info);

int coreblas_ssygst(int itype, coreblas_enum_t uplo,
                int n,
                float *A, int lda,
                float *B, int ldb);

void coreblas_ssymm(coreblas_enum_t side, coreblas_enum_t uplo,
                int m, int n,
                float alpha, const float *A, int lda,
                                          const float *B, int ldb,
                float beta,        float *C, int ldc);

void coreblas_ssyr2k(coreblas_enum_t uplo, coreblas_enum_t trans,
                 int n, int k,
                 float alpha, const float *A, int lda,
                                           const float *B, int ldb,
                 float beta,                    float *C, int ldc);

void coreblas_ssyrk(coreblas_enum_t uplo, coreblas_enum_t trans,
                int n, int k,
                float alpha, const float *A, int lda,
                float beta,        float *C, int ldc);

void coreblas_ssyssq(coreblas_enum_t uplo,
                 int n,
                 const float *A, int lda,
                 float *scale, float *sumsq);

void coreblas_ssyssq(coreblas_enum_t uplo,
                 int n,
                 const float *A, int lda,
                 float *scale, float *sumsq);

void coreblas_slacpy(coreblas_enum_t uplo, coreblas_enum_t transa,
                 int m, int n,
                 const float *A, int lda,
                       float *B, int ldb);

void coreblas_slacpy_lapack2tile_band(coreblas_enum_t uplo,
                                  int it, int jt,
                                  int m, int n, int nb, int kl, int ku,
                                  const float *A, int lda,
                                        float *B, int ldb);

void coreblas_slacpy_tile2lapack_band(coreblas_enum_t uplo,
                                  int it, int jt,
                                  int m, int n, int nb, int kl, int ku,
                                  const float *B, int ldb,
                                        float *A, int lda);

void coreblas_slange(coreblas_enum_t norm,
                 int m, int n,
                 const float *A, int lda,
                 float *work, float *result);

void coreblas_slansy(coreblas_enum_t norm, coreblas_enum_t uplo,
                 int n,
                 const float *A, int lda,
                 float *work, float *value);

void coreblas_slansy(coreblas_enum_t norm, coreblas_enum_t uplo,
                 int n,
                 const float *A, int lda,
                 float *work, float *value);

void coreblas_slantr(coreblas_enum_t norm, coreblas_enum_t uplo, coreblas_enum_t diag,
                 int m, int n,
                 const float *A, int lda,
                 float *work, float *value);

int coreblas_slarfb_gemm(coreblas_enum_t side, coreblas_enum_t trans, int direct, int storev,
                     int M, int N, int K,
                     const float *V, int LDV,
                     const float *T, int LDT,
                     float *C, int LDC,
                     float *WORK, int LDWORK);

void coreblas_slascl(coreblas_enum_t uplo,
                 float cfrom, float cto,
                 int m, int n,
                 float *A, int lda);

void coreblas_slaset(coreblas_enum_t uplo,
                 int m, int n,
                 float alpha, float beta,
                 float *A, int lda);
/*
void coreblas_sgeswp(coreblas_enum_t colrow,
                 coreblas_desc_t A, int k1, int k2, const int *ipiv, int incx);

void coreblas_ssyswp(int rank, int num_threads,
                 int uplo, coreblas_desc_t A, int k1, int k2, const int *ipiv,
                 int incx);
*/
int coreblas_slauum(coreblas_enum_t uplo,
                int n,
                float *A, int lda);

int coreblas_spamm(coreblas_enum_t op, coreblas_enum_t side, coreblas_enum_t storev,
               int m, int n, int k, int l,
               const float *A1, int lda1,
                     float *A2, int lda2,
               const float *V,  int ldv,
                     float *W,  int ldw);

int coreblas_sparfb(coreblas_enum_t side, coreblas_enum_t trans, coreblas_enum_t direct,
                coreblas_enum_t storev,
                int m1, int n1, int m2, int n2, int k, int l,
                      float *A1,   int lda1,
                      float *A2,   int lda2,
                const float *V,    int ldv,
                const float *T,    int ldt,
                      float *work, int ldwork);

int coreblas_spemv(coreblas_enum_t trans, int storev,
               int m, int n, int l,
               float alpha,
               const float *A, int lda,
               const float *X, int incx,
               float beta,
               float *Y, int incy,
               float *work);

int coreblas_spotrf(coreblas_enum_t uplo,
                int n,
                float *A, int lda);

void coreblas_ssymm(coreblas_enum_t side, coreblas_enum_t uplo,
                int m, int n,
                float alpha, const float *A, int lda,
                                          const float *B, int ldb,
                float beta,        float *C, int ldc);

void coreblas_ssyr2k(
    coreblas_enum_t uplo, coreblas_enum_t trans,
    int n, int k,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,        float *C, int ldc);

void coreblas_ssyrk(coreblas_enum_t uplo, coreblas_enum_t trans,
                int n, int k,
                float alpha, const float *A, int lda,
                float beta,        float *C, int ldc);

int coreblas_stradd(coreblas_enum_t uplo, coreblas_enum_t transa,
                int m, int n,
                float alpha, const float *A, int lda,
                float beta,        float *B, int ldb);

void coreblas_strmm(coreblas_enum_t side, coreblas_enum_t uplo,
                coreblas_enum_t transa, coreblas_enum_t diag,
                int m, int n,
                float alpha, const float *A, int lda,
                                                float *B, int ldb);

void coreblas_strsm(coreblas_enum_t side, coreblas_enum_t uplo,
                coreblas_enum_t transa, coreblas_enum_t diag,
                int m, int n,
                float alpha, const float *A, int lda,
                                                float *B, int ldb);

void coreblas_strssq(coreblas_enum_t uplo, coreblas_enum_t diag,
                 int m, int n,
                 const float *A, int lda,
                 float *scale, float *sumsq);

int coreblas_strtri(coreblas_enum_t uplo, coreblas_enum_t diag,
                int n,
                float *A, int lda);

int coreblas_stslqt(int m, int n, int ib,
                float *A1, int lda1,
                float *A2, int lda2,
                float *T,  int ldt,
                float *tau,
                float *work);

int coreblas_stsmlq(coreblas_enum_t side, coreblas_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      float *A1,   int lda1,
                      float *A2,   int lda2,
                const float *V,    int ldv,
                const float *T,    int ldt,
                      float *work, int ldwork);

int coreblas_stsmqr(coreblas_enum_t side, coreblas_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      float *A1,   int lda1,
                      float *A2,   int lda2,
                const float *V,    int ldv,
                const float *T,    int ldt,
                      float *work, int ldwork);

int coreblas_stsqrt(int m, int n, int ib,
                float *A1, int lda1,
                float *A2, int lda2,
                float *T,  int ldt,
                float *tau,
                float *work);

int coreblas_sttlqt(int m, int n, int ib,
                float *A1, int lda1,
                float *A2, int lda2,
                float *T,  int ldt,
                float *tau,
                float *work);

int coreblas_sttmlq(coreblas_enum_t side, coreblas_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      float *A1,   int lda1,
                      float *A2,   int lda2,
                const float *V,    int ldv,
                const float *T,    int ldt,
                      float *work, int ldwork);

int coreblas_sttmqr(coreblas_enum_t side, coreblas_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      float *A1,   int lda1,
                      float *A2,   int lda2,
                const float *V,    int ldv,
                const float *T,    int ldt,
                      float *work, int ldwork);

int coreblas_sttqrt(int m, int n, int ib,
                float *A1, int lda1,
                float *A2, int lda2,
                float *T,  int ldt,
                float *tau,
                float *work);

int coreblas_sormlq(coreblas_enum_t side, coreblas_enum_t trans,
                int m, int n, int k, int ib,
                const float *A,    int lda,
                const float *T,    int ldt,
                      float *C,    int ldc,
                      float *work, int ldwork);

int coreblas_sormqr(coreblas_enum_t side, coreblas_enum_t trans,
                int m, int n, int k, int ib,
                const float *A,    int lda,
                const float *T,    int ldt,
                      float *C,    int ldc,
                      float *work, int ldwork);

/******************************************************************************/
void coreblas_kernel_samax(int colrow, int m, int n,
                     const float *A, int lda,
                     float *values);

void coreblas_kernel_sgeadd(
    coreblas_enum_t transa, int m, int n,
    float alpha, const float *A, int lda,
    float beta,        float *B, int ldb);

void coreblas_kernel_sgelqt(int m, int n, int ib,
                     float *A, int lda,
                     float *T, int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_sgemm(
    coreblas_enum_t transa, coreblas_enum_t transb,
    int m, int n, int k,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,        float *C, int ldc);

void coreblas_kernel_sgeqrt(int m, int n, int ib,
                     float *A, int lda,
                     float *T, int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_sgessq(int m, int n,
                     const float *A, int lda,
                     float *scale, float *sumsq);

void coreblas_kernel_sgessq_aux(int n,
                         const float *scale, const float *sumsq,
                         float *value);

void coreblas_kernel_ssygst(int itype, coreblas_enum_t uplo,
                     int n,
                     float *A, int lda,
                     float *B, int ldb);

void coreblas_kernel_ssymm(
    coreblas_enum_t side, coreblas_enum_t uplo,
    int m, int n,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,        float *C, int ldc);

void coreblas_kernel_ssyr2k(
    coreblas_enum_t uplo, coreblas_enum_t trans,
    int n, int k,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,                    float *C, int ldc);

void coreblas_kernel_ssyrk(coreblas_enum_t uplo, coreblas_enum_t trans,
                    int n, int k,
                    float alpha, const float *A, int lda,
                    float beta,        float *C, int ldc);

void coreblas_kernel_ssyssq(coreblas_enum_t uplo,
                     int n,
                     const float *A, int lda,
                     float *scale, float *sumsq);

void coreblas_kernel_ssyssq(coreblas_enum_t uplo,
                     int n,
                     const float *A, int lda,
                     float *scale, float *sumsq);

void coreblas_kernel_ssyssq_aux(int m, int n,
                         const float *scale, const float *sumsq,
                         float *value);

void coreblas_kernel_slacpy(coreblas_enum_t uplo, coreblas_enum_t transa,
                     int m, int n,
                     const float *A, int lda,
                           float *B, int ldb);

void coreblas_kernel_slacpy_lapack2tile_band(coreblas_enum_t uplo,
                                      int it, int jt,
                                      int m, int n, int nb, int kl, int ku,
                                      const float *A, int lda,
                                            float *B, int ldb);

void coreblas_kernel_slacpy_tile2lapack_band(coreblas_enum_t uplo,
                                      int it, int jt,
                                      int m, int n, int nb, int kl, int ku,
                                      const float *B, int ldb,
                                            float *A, int lda);

void coreblas_kernel_slange(coreblas_enum_t norm,
                     int m, int n,
                     const float *A, int lda,
                     float *work, float *result);

void coreblas_kernel_slange_aux(coreblas_enum_t norm,
                         int m, int n,
                         const float *A, int lda,
                         float *value);

void coreblas_kernel_slansy(coreblas_enum_t norm, coreblas_enum_t uplo,
                     int n,
                     const float *A, int lda,
                     float *work, float *value);

void coreblas_kernel_slansy_aux(coreblas_enum_t norm, coreblas_enum_t uplo,
                         int n,
                         const float *A, int lda,
                         float *value);

void coreblas_kernel_slansy(coreblas_enum_t norm, coreblas_enum_t uplo,
                     int n,
                     const float *A, int lda,
                     float *work, float *value);

void coreblas_kernel_slansy_aux(coreblas_enum_t norm, coreblas_enum_t uplo,
                         int n,
                         const float *A, int lda,
                         float *value);

void coreblas_kernel_slantr(coreblas_enum_t norm, coreblas_enum_t uplo, coreblas_enum_t diag,
                     int m, int n,
                     const float *A, int lda,
                     float *work, float *value);

void coreblas_kernel_slantr_aux(coreblas_enum_t norm, coreblas_enum_t uplo,
                         coreblas_enum_t diag,
                         int m, int n,
                         const float *A, int lda,
                         float *value);

void coreblas_kernel_slascl(coreblas_enum_t uplo,
                     float cfrom, float cto,
                     int m, int n,
                     float *A, int lda);

void coreblas_kernel_slaset(coreblas_enum_t uplo,
                     int mb, int nb,
                     int i, int j,
                     int m, int n,
                     float alpha, float beta,
                     float *A);

void coreblas_kernel_slauum(coreblas_enum_t uplo,
                     int n,
                     float *A, int lda);

void coreblas_kernel_spotrf(coreblas_enum_t uplo,
                     int n,
                     float *A, int lda,
                     int iinfo);

void coreblas_kernel_ssymm(
    coreblas_enum_t side, coreblas_enum_t uplo,
    int m, int n,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,        float *C, int ldc);

void coreblas_kernel_ssyr2k(
    coreblas_enum_t uplo, coreblas_enum_t trans,
    int n, int k,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,        float *C, int ldc);

void coreblas_kernel_ssyrk(
    coreblas_enum_t uplo, coreblas_enum_t trans,
    int n, int k,
    float alpha, const float *A, int lda,
    float beta,        float *C, int ldc);

void coreblas_kernel_stradd(
    coreblas_enum_t uplo, coreblas_enum_t transa,
    int m, int n,
    float alpha, const float *A, int lda,
    float beta,        float *B, int ldb);

void coreblas_kernel_strmm(
    coreblas_enum_t side, coreblas_enum_t uplo,
    coreblas_enum_t transa, coreblas_enum_t diag,
    int m, int n,
    float alpha, const float *A, int lda,
                                    float *B, int ldb);

void coreblas_kernel_strsm(
    coreblas_enum_t side, coreblas_enum_t uplo,
    coreblas_enum_t transa, coreblas_enum_t diag,
    int m, int n,
    float alpha, const float *A, int lda,
                                    float *B, int ldb);

void coreblas_kernel_strssq(coreblas_enum_t uplo, coreblas_enum_t diag,
                     int m, int n,
                     const float *A, int lda,
                     float *scale, float *sumsq);

void coreblas_kernel_strtri(coreblas_enum_t uplo, coreblas_enum_t diag,
                     int n,
                     float *A, int lda,
                     int iinfo);

void coreblas_kernel_stslqt(int m, int n, int ib,
                     float *A1, int lda1,
                     float *A2, int lda2,
                     float *T,  int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_stsmlq(coreblas_enum_t side, coreblas_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           float *A1, int lda1,
                           float *A2, int lda2,
                     const float *V,  int ldv,
                     const float *T,  int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_stsmqr(coreblas_enum_t side, coreblas_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           float *A1, int lda1,
                           float *A2, int lda2,
                     const float *V, int ldv,
                     const float *T, int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_stsqrt(int m, int n, int ib,
                     float *A1, int lda1,
                     float *A2, int lda2,
                     float *T,  int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_sttlqt(int m, int n, int ib,
                     float *A1, int lda1,
                     float *A2, int lda2,
                     float *T,  int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_sttmlq(coreblas_enum_t side, coreblas_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           float *A1, int lda1,
                           float *A2, int lda2,
                     const float *V,  int ldv,
                     const float *T,  int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_sttmqr(coreblas_enum_t side, coreblas_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           float *A1, int lda1,
                           float *A2, int lda2,
                     const float *V, int ldv,
                     const float *T, int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_sttqrt(int m, int n, int ib,
                     float *A1, int lda1,
                     float *A2, int lda2,
                     float *T,  int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_sormlq(coreblas_enum_t side, coreblas_enum_t trans,
                     int m, int n, int k, int ib,
                     const float *A, int lda,
                     const float *T, int ldt,
                           float *C, int ldc,
                     coreblas_workspace_t work);

void coreblas_kernel_sormqr(coreblas_enum_t side, coreblas_enum_t trans,
                     int m, int n, int k, int ib,
                     const float *A, int lda,
                     const float *T, int ldt,
                           float *C, int ldc,
                     coreblas_workspace_t work);

#undef REAL

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // COREBLAS_CORE_BLAS_S_H
