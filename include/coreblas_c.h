/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from include/coreblas_z.h, normal z -> c, Mon Mar 18 06:35:29 2024
 *
 **/
#ifndef COREBLAS_CORE_BLAS_C_H
#define COREBLAS_CORE_BLAS_C_H

#include "coreblas_types.h"
#include "coreblas_workspace.h"

#ifdef __cplusplus
extern "C" {
#endif

#define COMPLEX

/******************************************************************************/
#ifdef COMPLEX
float coreblas_scabs1(coreblas_complex32_t alpha);
#endif

void coreblas_cgbtype1cb(coreblas_enum_t uplo, int n, int nb,
                      coreblas_complex32_t *A, int lda,
                      coreblas_complex32_t *VQ, coreblas_complex32_t *TAUQ,
                      coreblas_complex32_t *VP, coreblas_complex32_t *TAUP,
                      int st, int ed, int sweep, int Vblksiz, int WANTZ,
                      coreblas_complex32_t *work);
    
void coreblas_cgbtype2cb(coreblas_enum_t uplo, int n, int nb,
                      coreblas_complex32_t *A, int lda,
                      coreblas_complex32_t *VQ, coreblas_complex32_t *TAUQ,
                      coreblas_complex32_t *VP, coreblas_complex32_t *TAUP,
                      int st, int ed, int sweep, int Vblksiz, int WANTZ,
                      coreblas_complex32_t *work);
    
void coreblas_cgbtype3cb(coreblas_enum_t uplo, int n, int nb,
                      coreblas_complex32_t *A, int lda,
                      coreblas_complex32_t *VQ, coreblas_complex32_t *TAUQ,
                      coreblas_complex32_t *VP, coreblas_complex32_t *TAUP,
                      int st, int ed, int sweep, int Vblksiz, int WANTZ,
                      coreblas_complex32_t *work);
    
int coreblas_cgeadd(coreblas_enum_t transa,
                int m, int n,
                coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                coreblas_complex32_t beta,        coreblas_complex32_t *B, int ldb);

int coreblas_cgelqt(int m, int n, int ib,
                coreblas_complex32_t *A, int lda,
                coreblas_complex32_t *T, int ldt,
                coreblas_complex32_t *tau,
                coreblas_complex32_t *work);

void coreblas_cgemm(coreblas_enum_t transa, coreblas_enum_t transb,
                int m, int n, int k,
                coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                                          const coreblas_complex32_t *B, int ldb,
                coreblas_complex32_t beta,        coreblas_complex32_t *C, int ldc);

int coreblas_cgeqrt(int m, int n, int ib,
                coreblas_complex32_t *A, int lda,
                coreblas_complex32_t *T, int ldt,
                coreblas_complex32_t *tau,
                coreblas_complex32_t *work);

void coreblas_cgessq(int m, int n,
                 const coreblas_complex32_t *A, int lda,
                 float *scale, float *sumsq);

//void coreblas_cgetrf(coreblas_desc_t A, int *ipiv, int ib, int rank, int size,
//                 volatile int *max_idx, volatile coreblas_complex32_t *max_val,
//                 volatile int *info);

int coreblas_chegst(int itype, coreblas_enum_t uplo,
                int n,
                coreblas_complex32_t *A, int lda,
                coreblas_complex32_t *B, int ldb);

void coreblas_chemm(coreblas_enum_t side, coreblas_enum_t uplo,
                int m, int n,
                coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                                          const coreblas_complex32_t *B, int ldb,
                coreblas_complex32_t beta,        coreblas_complex32_t *C, int ldc);

void coreblas_cher2k(coreblas_enum_t uplo, coreblas_enum_t trans,
                 int n, int k,
                 coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                                           const coreblas_complex32_t *B, int ldb,
                 float beta,                    coreblas_complex32_t *C, int ldc);

void coreblas_cherk(coreblas_enum_t uplo, coreblas_enum_t trans,
                int n, int k,
                float alpha, const coreblas_complex32_t *A, int lda,
                float beta,        coreblas_complex32_t *C, int ldc);

void coreblas_chessq(coreblas_enum_t uplo,
                 int n,
                 const coreblas_complex32_t *A, int lda,
                 float *scale, float *sumsq);

void coreblas_csyssq(coreblas_enum_t uplo,
                 int n,
                 const coreblas_complex32_t *A, int lda,
                 float *scale, float *sumsq);

void coreblas_clacpy(coreblas_enum_t uplo, coreblas_enum_t transa,
                 int m, int n,
                 const coreblas_complex32_t *A, int lda,
                       coreblas_complex32_t *B, int ldb);

void coreblas_clacpy_lapack2tile_band(coreblas_enum_t uplo,
                                  int it, int jt,
                                  int m, int n, int nb, int kl, int ku,
                                  const coreblas_complex32_t *A, int lda,
                                        coreblas_complex32_t *B, int ldb);

void coreblas_clacpy_tile2lapack_band(coreblas_enum_t uplo,
                                  int it, int jt,
                                  int m, int n, int nb, int kl, int ku,
                                  const coreblas_complex32_t *B, int ldb,
                                        coreblas_complex32_t *A, int lda);

void coreblas_clange(coreblas_enum_t norm,
                 int m, int n,
                 const coreblas_complex32_t *A, int lda,
                 float *work, float *result);

void coreblas_clanhe(coreblas_enum_t norm, coreblas_enum_t uplo,
                 int n,
                 const coreblas_complex32_t *A, int lda,
                 float *work, float *value);

void coreblas_clansy(coreblas_enum_t norm, coreblas_enum_t uplo,
                 int n,
                 const coreblas_complex32_t *A, int lda,
                 float *work, float *value);

void coreblas_clantr(coreblas_enum_t norm, coreblas_enum_t uplo, coreblas_enum_t diag,
                 int m, int n,
                 const coreblas_complex32_t *A, int lda,
                 float *work, float *value);

int coreblas_clarfb_gemm(coreblas_enum_t side, coreblas_enum_t trans, int direct, int storev,
                     int M, int N, int K,
                     const coreblas_complex32_t *V, int LDV,
                     const coreblas_complex32_t *T, int LDT,
                     coreblas_complex32_t *C, int LDC,
                     coreblas_complex32_t *WORK, int LDWORK);

void coreblas_clascl(coreblas_enum_t uplo,
                 float cfrom, float cto,
                 int m, int n,
                 coreblas_complex32_t *A, int lda);

void coreblas_claset(coreblas_enum_t uplo,
                 int m, int n,
                 coreblas_complex32_t alpha, coreblas_complex32_t beta,
                 coreblas_complex32_t *A, int lda);
/*
void coreblas_cgeswp(coreblas_enum_t colrow,
                 coreblas_desc_t A, int k1, int k2, const int *ipiv, int incx);

void coreblas_cheswp(int rank, int num_threads,
                 int uplo, coreblas_desc_t A, int k1, int k2, const int *ipiv,
                 int incx);
*/
int coreblas_clauum(coreblas_enum_t uplo,
                int n,
                coreblas_complex32_t *A, int lda);

int coreblas_cpamm(coreblas_enum_t op, coreblas_enum_t side, coreblas_enum_t storev,
               int m, int n, int k, int l,
               const coreblas_complex32_t *A1, int lda1,
                     coreblas_complex32_t *A2, int lda2,
               const coreblas_complex32_t *V,  int ldv,
                     coreblas_complex32_t *W,  int ldw);

int coreblas_cparfb(coreblas_enum_t side, coreblas_enum_t trans, coreblas_enum_t direct,
                coreblas_enum_t storev,
                int m1, int n1, int m2, int n2, int k, int l,
                      coreblas_complex32_t *A1,   int lda1,
                      coreblas_complex32_t *A2,   int lda2,
                const coreblas_complex32_t *V,    int ldv,
                const coreblas_complex32_t *T,    int ldt,
                      coreblas_complex32_t *work, int ldwork);

int coreblas_cpemv(coreblas_enum_t trans, int storev,
               int m, int n, int l,
               coreblas_complex32_t alpha,
               const coreblas_complex32_t *A, int lda,
               const coreblas_complex32_t *X, int incx,
               coreblas_complex32_t beta,
               coreblas_complex32_t *Y, int incy,
               coreblas_complex32_t *work);

int coreblas_cpotrf(coreblas_enum_t uplo,
                int n,
                coreblas_complex32_t *A, int lda);

void coreblas_csymm(coreblas_enum_t side, coreblas_enum_t uplo,
                int m, int n,
                coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                                          const coreblas_complex32_t *B, int ldb,
                coreblas_complex32_t beta,        coreblas_complex32_t *C, int ldc);

void coreblas_csyr2k(
    coreblas_enum_t uplo, coreblas_enum_t trans,
    int n, int k,
    coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                              const coreblas_complex32_t *B, int ldb,
    coreblas_complex32_t beta,        coreblas_complex32_t *C, int ldc);

void coreblas_csyrk(coreblas_enum_t uplo, coreblas_enum_t trans,
                int n, int k,
                coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                coreblas_complex32_t beta,        coreblas_complex32_t *C, int ldc);

int coreblas_ctradd(coreblas_enum_t uplo, coreblas_enum_t transa,
                int m, int n,
                coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                coreblas_complex32_t beta,        coreblas_complex32_t *B, int ldb);

void coreblas_ctrmm(coreblas_enum_t side, coreblas_enum_t uplo,
                coreblas_enum_t transa, coreblas_enum_t diag,
                int m, int n,
                coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                                                coreblas_complex32_t *B, int ldb);

void coreblas_ctrsm(coreblas_enum_t side, coreblas_enum_t uplo,
                coreblas_enum_t transa, coreblas_enum_t diag,
                int m, int n,
                coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                                                coreblas_complex32_t *B, int ldb);

void coreblas_ctrssq(coreblas_enum_t uplo, coreblas_enum_t diag,
                 int m, int n,
                 const coreblas_complex32_t *A, int lda,
                 float *scale, float *sumsq);

int coreblas_ctrtri(coreblas_enum_t uplo, coreblas_enum_t diag,
                int n,
                coreblas_complex32_t *A, int lda);

int coreblas_ctslqt(int m, int n, int ib,
                coreblas_complex32_t *A1, int lda1,
                coreblas_complex32_t *A2, int lda2,
                coreblas_complex32_t *T,  int ldt,
                coreblas_complex32_t *tau,
                coreblas_complex32_t *work);

int coreblas_ctsmlq(coreblas_enum_t side, coreblas_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      coreblas_complex32_t *A1,   int lda1,
                      coreblas_complex32_t *A2,   int lda2,
                const coreblas_complex32_t *V,    int ldv,
                const coreblas_complex32_t *T,    int ldt,
                      coreblas_complex32_t *work, int ldwork);

int coreblas_ctsmqr(coreblas_enum_t side, coreblas_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      coreblas_complex32_t *A1,   int lda1,
                      coreblas_complex32_t *A2,   int lda2,
                const coreblas_complex32_t *V,    int ldv,
                const coreblas_complex32_t *T,    int ldt,
                      coreblas_complex32_t *work, int ldwork);

int coreblas_ctsqrt(int m, int n, int ib,
                coreblas_complex32_t *A1, int lda1,
                coreblas_complex32_t *A2, int lda2,
                coreblas_complex32_t *T,  int ldt,
                coreblas_complex32_t *tau,
                coreblas_complex32_t *work);

int coreblas_cttlqt(int m, int n, int ib,
                coreblas_complex32_t *A1, int lda1,
                coreblas_complex32_t *A2, int lda2,
                coreblas_complex32_t *T,  int ldt,
                coreblas_complex32_t *tau,
                coreblas_complex32_t *work);

int coreblas_cttmlq(coreblas_enum_t side, coreblas_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      coreblas_complex32_t *A1,   int lda1,
                      coreblas_complex32_t *A2,   int lda2,
                const coreblas_complex32_t *V,    int ldv,
                const coreblas_complex32_t *T,    int ldt,
                      coreblas_complex32_t *work, int ldwork);

int coreblas_cttmqr(coreblas_enum_t side, coreblas_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      coreblas_complex32_t *A1,   int lda1,
                      coreblas_complex32_t *A2,   int lda2,
                const coreblas_complex32_t *V,    int ldv,
                const coreblas_complex32_t *T,    int ldt,
                      coreblas_complex32_t *work, int ldwork);

int coreblas_cttqrt(int m, int n, int ib,
                coreblas_complex32_t *A1, int lda1,
                coreblas_complex32_t *A2, int lda2,
                coreblas_complex32_t *T,  int ldt,
                coreblas_complex32_t *tau,
                coreblas_complex32_t *work);

int coreblas_cunmlq(coreblas_enum_t side, coreblas_enum_t trans,
                int m, int n, int k, int ib,
                const coreblas_complex32_t *A,    int lda,
                const coreblas_complex32_t *T,    int ldt,
                      coreblas_complex32_t *C,    int ldc,
                      coreblas_complex32_t *work, int ldwork);

int coreblas_cunmqr(coreblas_enum_t side, coreblas_enum_t trans,
                int m, int n, int k, int ib,
                const coreblas_complex32_t *A,    int lda,
                const coreblas_complex32_t *T,    int ldt,
                      coreblas_complex32_t *C,    int ldc,
                      coreblas_complex32_t *work, int ldwork);

/******************************************************************************/
void coreblas_kernel_scamax(int colrow, int m, int n,
                     const coreblas_complex32_t *A, int lda,
                     float *values);

void coreblas_kernel_cgeadd(
    coreblas_enum_t transa, int m, int n,
    coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
    coreblas_complex32_t beta,        coreblas_complex32_t *B, int ldb);

void coreblas_kernel_cgelqt(int m, int n, int ib,
                     coreblas_complex32_t *A, int lda,
                     coreblas_complex32_t *T, int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_cgemm(
    coreblas_enum_t transa, coreblas_enum_t transb,
    int m, int n, int k,
    coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                              const coreblas_complex32_t *B, int ldb,
    coreblas_complex32_t beta,        coreblas_complex32_t *C, int ldc);

void coreblas_kernel_cgeqrt(int m, int n, int ib,
                     coreblas_complex32_t *A, int lda,
                     coreblas_complex32_t *T, int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_cgessq(int m, int n,
                     const coreblas_complex32_t *A, int lda,
                     float *scale, float *sumsq);

void coreblas_kernel_cgessq_aux(int n,
                         const float *scale, const float *sumsq,
                         float *value);

void coreblas_kernel_chegst(int itype, coreblas_enum_t uplo,
                     int n,
                     coreblas_complex32_t *A, int lda,
                     coreblas_complex32_t *B, int ldb);

void coreblas_kernel_chemm(
    coreblas_enum_t side, coreblas_enum_t uplo,
    int m, int n,
    coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                              const coreblas_complex32_t *B, int ldb,
    coreblas_complex32_t beta,        coreblas_complex32_t *C, int ldc);

void coreblas_kernel_cher2k(
    coreblas_enum_t uplo, coreblas_enum_t trans,
    int n, int k,
    coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                              const coreblas_complex32_t *B, int ldb,
    float beta,                    coreblas_complex32_t *C, int ldc);

void coreblas_kernel_cherk(coreblas_enum_t uplo, coreblas_enum_t trans,
                    int n, int k,
                    float alpha, const coreblas_complex32_t *A, int lda,
                    float beta,        coreblas_complex32_t *C, int ldc);

void coreblas_kernel_chessq(coreblas_enum_t uplo,
                     int n,
                     const coreblas_complex32_t *A, int lda,
                     float *scale, float *sumsq);

void coreblas_kernel_csyssq(coreblas_enum_t uplo,
                     int n,
                     const coreblas_complex32_t *A, int lda,
                     float *scale, float *sumsq);

void coreblas_kernel_csyssq_aux(int m, int n,
                         const float *scale, const float *sumsq,
                         float *value);

void coreblas_kernel_clacpy(coreblas_enum_t uplo, coreblas_enum_t transa,
                     int m, int n,
                     const coreblas_complex32_t *A, int lda,
                           coreblas_complex32_t *B, int ldb);

void coreblas_kernel_clacpy_lapack2tile_band(coreblas_enum_t uplo,
                                      int it, int jt,
                                      int m, int n, int nb, int kl, int ku,
                                      const coreblas_complex32_t *A, int lda,
                                            coreblas_complex32_t *B, int ldb);

void coreblas_kernel_clacpy_tile2lapack_band(coreblas_enum_t uplo,
                                      int it, int jt,
                                      int m, int n, int nb, int kl, int ku,
                                      const coreblas_complex32_t *B, int ldb,
                                            coreblas_complex32_t *A, int lda);

void coreblas_kernel_clange(coreblas_enum_t norm,
                     int m, int n,
                     const coreblas_complex32_t *A, int lda,
                     float *work, float *result);

void coreblas_kernel_clange_aux(coreblas_enum_t norm,
                         int m, int n,
                         const coreblas_complex32_t *A, int lda,
                         float *value);

void coreblas_kernel_clanhe(coreblas_enum_t norm, coreblas_enum_t uplo,
                     int n,
                     const coreblas_complex32_t *A, int lda,
                     float *work, float *value);

void coreblas_kernel_clanhe_aux(coreblas_enum_t norm, coreblas_enum_t uplo,
                         int n,
                         const coreblas_complex32_t *A, int lda,
                         float *value);

void coreblas_kernel_clansy(coreblas_enum_t norm, coreblas_enum_t uplo,
                     int n,
                     const coreblas_complex32_t *A, int lda,
                     float *work, float *value);

void coreblas_kernel_clansy_aux(coreblas_enum_t norm, coreblas_enum_t uplo,
                         int n,
                         const coreblas_complex32_t *A, int lda,
                         float *value);

void coreblas_kernel_clantr(coreblas_enum_t norm, coreblas_enum_t uplo, coreblas_enum_t diag,
                     int m, int n,
                     const coreblas_complex32_t *A, int lda,
                     float *work, float *value);

void coreblas_kernel_clantr_aux(coreblas_enum_t norm, coreblas_enum_t uplo,
                         coreblas_enum_t diag,
                         int m, int n,
                         const coreblas_complex32_t *A, int lda,
                         float *value);

void coreblas_kernel_clascl(coreblas_enum_t uplo,
                     float cfrom, float cto,
                     int m, int n,
                     coreblas_complex32_t *A, int lda);

void coreblas_kernel_claset(coreblas_enum_t uplo,
                     int mb, int nb,
                     int i, int j,
                     int m, int n,
                     coreblas_complex32_t alpha, coreblas_complex32_t beta,
                     coreblas_complex32_t *A);

void coreblas_kernel_clauum(coreblas_enum_t uplo,
                     int n,
                     coreblas_complex32_t *A, int lda);

void coreblas_kernel_cpotrf(coreblas_enum_t uplo,
                     int n,
                     coreblas_complex32_t *A, int lda,
                     int iinfo);

void coreblas_kernel_csymm(
    coreblas_enum_t side, coreblas_enum_t uplo,
    int m, int n,
    coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                              const coreblas_complex32_t *B, int ldb,
    coreblas_complex32_t beta,        coreblas_complex32_t *C, int ldc);

void coreblas_kernel_csyr2k(
    coreblas_enum_t uplo, coreblas_enum_t trans,
    int n, int k,
    coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                              const coreblas_complex32_t *B, int ldb,
    coreblas_complex32_t beta,        coreblas_complex32_t *C, int ldc);

void coreblas_kernel_csyrk(
    coreblas_enum_t uplo, coreblas_enum_t trans,
    int n, int k,
    coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
    coreblas_complex32_t beta,        coreblas_complex32_t *C, int ldc);

void coreblas_kernel_ctradd(
    coreblas_enum_t uplo, coreblas_enum_t transa,
    int m, int n,
    coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
    coreblas_complex32_t beta,        coreblas_complex32_t *B, int ldb);

void coreblas_kernel_ctrmm(
    coreblas_enum_t side, coreblas_enum_t uplo,
    coreblas_enum_t transa, coreblas_enum_t diag,
    int m, int n,
    coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                                    coreblas_complex32_t *B, int ldb);

void coreblas_kernel_ctrsm(
    coreblas_enum_t side, coreblas_enum_t uplo,
    coreblas_enum_t transa, coreblas_enum_t diag,
    int m, int n,
    coreblas_complex32_t alpha, const coreblas_complex32_t *A, int lda,
                                    coreblas_complex32_t *B, int ldb);

void coreblas_kernel_ctrssq(coreblas_enum_t uplo, coreblas_enum_t diag,
                     int m, int n,
                     const coreblas_complex32_t *A, int lda,
                     float *scale, float *sumsq);

void coreblas_kernel_ctrtri(coreblas_enum_t uplo, coreblas_enum_t diag,
                     int n,
                     coreblas_complex32_t *A, int lda,
                     int iinfo);

void coreblas_kernel_ctslqt(int m, int n, int ib,
                     coreblas_complex32_t *A1, int lda1,
                     coreblas_complex32_t *A2, int lda2,
                     coreblas_complex32_t *T,  int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_ctsmlq(coreblas_enum_t side, coreblas_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           coreblas_complex32_t *A1, int lda1,
                           coreblas_complex32_t *A2, int lda2,
                     const coreblas_complex32_t *V,  int ldv,
                     const coreblas_complex32_t *T,  int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_ctsmqr(coreblas_enum_t side, coreblas_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           coreblas_complex32_t *A1, int lda1,
                           coreblas_complex32_t *A2, int lda2,
                     const coreblas_complex32_t *V, int ldv,
                     const coreblas_complex32_t *T, int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_ctsqrt(int m, int n, int ib,
                     coreblas_complex32_t *A1, int lda1,
                     coreblas_complex32_t *A2, int lda2,
                     coreblas_complex32_t *T,  int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_cttlqt(int m, int n, int ib,
                     coreblas_complex32_t *A1, int lda1,
                     coreblas_complex32_t *A2, int lda2,
                     coreblas_complex32_t *T,  int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_cttmlq(coreblas_enum_t side, coreblas_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           coreblas_complex32_t *A1, int lda1,
                           coreblas_complex32_t *A2, int lda2,
                     const coreblas_complex32_t *V,  int ldv,
                     const coreblas_complex32_t *T,  int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_cttmqr(coreblas_enum_t side, coreblas_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           coreblas_complex32_t *A1, int lda1,
                           coreblas_complex32_t *A2, int lda2,
                     const coreblas_complex32_t *V, int ldv,
                     const coreblas_complex32_t *T, int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_cttqrt(int m, int n, int ib,
                     coreblas_complex32_t *A1, int lda1,
                     coreblas_complex32_t *A2, int lda2,
                     coreblas_complex32_t *T,  int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_cunmlq(coreblas_enum_t side, coreblas_enum_t trans,
                     int m, int n, int k, int ib,
                     const coreblas_complex32_t *A, int lda,
                     const coreblas_complex32_t *T, int ldt,
                           coreblas_complex32_t *C, int ldc,
                     coreblas_workspace_t work);

void coreblas_kernel_cunmqr(coreblas_enum_t side, coreblas_enum_t trans,
                     int m, int n, int k, int ib,
                     const coreblas_complex32_t *A, int lda,
                     const coreblas_complex32_t *T, int ldt,
                           coreblas_complex32_t *C, int ldc,
                     coreblas_workspace_t work);

#undef COMPLEX

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // COREBLAS_CORE_BLAS_C_H
