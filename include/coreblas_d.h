/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from include/coreblas_z.h, normal z -> d, Mon Mar 18 06:35:29 2024
 *
 **/
#ifndef COREBLAS_CORE_BLAS_D_H
#define COREBLAS_CORE_BLAS_D_H

#include "coreblas_types.h"
#include "coreblas_workspace.h"

#ifdef __cplusplus
extern "C" {
#endif

#define REAL

/******************************************************************************/
#ifdef COMPLEX
double fabs(double alpha);
#endif

void coreblas_dgbtype1cb(coreblas_enum_t uplo, int n, int nb,
                      double *A, int lda,
                      double *VQ, double *TAUQ,
                      double *VP, double *TAUP,
                      int st, int ed, int sweep, int Vblksiz, int WANTZ,
                      double *work);
    
void coreblas_dgbtype2cb(coreblas_enum_t uplo, int n, int nb,
                      double *A, int lda,
                      double *VQ, double *TAUQ,
                      double *VP, double *TAUP,
                      int st, int ed, int sweep, int Vblksiz, int WANTZ,
                      double *work);
    
void coreblas_dgbtype3cb(coreblas_enum_t uplo, int n, int nb,
                      double *A, int lda,
                      double *VQ, double *TAUQ,
                      double *VP, double *TAUP,
                      int st, int ed, int sweep, int Vblksiz, int WANTZ,
                      double *work);
    
int coreblas_dgeadd(coreblas_enum_t transa,
                int m, int n,
                double alpha, const double *A, int lda,
                double beta,        double *B, int ldb);

int coreblas_dgelqt(int m, int n, int ib,
                double *A, int lda,
                double *T, int ldt,
                double *tau,
                double *work);

void coreblas_dgemm(coreblas_enum_t transa, coreblas_enum_t transb,
                int m, int n, int k,
                double alpha, const double *A, int lda,
                                          const double *B, int ldb,
                double beta,        double *C, int ldc);

int coreblas_dgeqrt(int m, int n, int ib,
                double *A, int lda,
                double *T, int ldt,
                double *tau,
                double *work);

void coreblas_dgessq(int m, int n,
                 const double *A, int lda,
                 double *scale, double *sumsq);

//void coreblas_dgetrf(coreblas_desc_t A, int *ipiv, int ib, int rank, int size,
//                 volatile int *max_idx, volatile double *max_val,
//                 volatile int *info);

int coreblas_dsygst(int itype, coreblas_enum_t uplo,
                int n,
                double *A, int lda,
                double *B, int ldb);

void coreblas_dsymm(coreblas_enum_t side, coreblas_enum_t uplo,
                int m, int n,
                double alpha, const double *A, int lda,
                                          const double *B, int ldb,
                double beta,        double *C, int ldc);

void coreblas_dsyr2k(coreblas_enum_t uplo, coreblas_enum_t trans,
                 int n, int k,
                 double alpha, const double *A, int lda,
                                           const double *B, int ldb,
                 double beta,                    double *C, int ldc);

void coreblas_dsyrk(coreblas_enum_t uplo, coreblas_enum_t trans,
                int n, int k,
                double alpha, const double *A, int lda,
                double beta,        double *C, int ldc);

void coreblas_dsyssq(coreblas_enum_t uplo,
                 int n,
                 const double *A, int lda,
                 double *scale, double *sumsq);

void coreblas_dsyssq(coreblas_enum_t uplo,
                 int n,
                 const double *A, int lda,
                 double *scale, double *sumsq);

void coreblas_dlacpy(coreblas_enum_t uplo, coreblas_enum_t transa,
                 int m, int n,
                 const double *A, int lda,
                       double *B, int ldb);

void coreblas_dlacpy_lapack2tile_band(coreblas_enum_t uplo,
                                  int it, int jt,
                                  int m, int n, int nb, int kl, int ku,
                                  const double *A, int lda,
                                        double *B, int ldb);

void coreblas_dlacpy_tile2lapack_band(coreblas_enum_t uplo,
                                  int it, int jt,
                                  int m, int n, int nb, int kl, int ku,
                                  const double *B, int ldb,
                                        double *A, int lda);

void coreblas_dlange(coreblas_enum_t norm,
                 int m, int n,
                 const double *A, int lda,
                 double *work, double *result);

void coreblas_dlansy(coreblas_enum_t norm, coreblas_enum_t uplo,
                 int n,
                 const double *A, int lda,
                 double *work, double *value);

void coreblas_dlansy(coreblas_enum_t norm, coreblas_enum_t uplo,
                 int n,
                 const double *A, int lda,
                 double *work, double *value);

void coreblas_dlantr(coreblas_enum_t norm, coreblas_enum_t uplo, coreblas_enum_t diag,
                 int m, int n,
                 const double *A, int lda,
                 double *work, double *value);

int coreblas_dlarfb_gemm(coreblas_enum_t side, coreblas_enum_t trans, int direct, int storev,
                     int M, int N, int K,
                     const double *V, int LDV,
                     const double *T, int LDT,
                     double *C, int LDC,
                     double *WORK, int LDWORK);

void coreblas_dlascl(coreblas_enum_t uplo,
                 double cfrom, double cto,
                 int m, int n,
                 double *A, int lda);

void coreblas_dlaset(coreblas_enum_t uplo,
                 int m, int n,
                 double alpha, double beta,
                 double *A, int lda);
/*
void coreblas_dgeswp(coreblas_enum_t colrow,
                 coreblas_desc_t A, int k1, int k2, const int *ipiv, int incx);

void coreblas_dsyswp(int rank, int num_threads,
                 int uplo, coreblas_desc_t A, int k1, int k2, const int *ipiv,
                 int incx);
*/
int coreblas_dlauum(coreblas_enum_t uplo,
                int n,
                double *A, int lda);

int coreblas_dpamm(coreblas_enum_t op, coreblas_enum_t side, coreblas_enum_t storev,
               int m, int n, int k, int l,
               const double *A1, int lda1,
                     double *A2, int lda2,
               const double *V,  int ldv,
                     double *W,  int ldw);

int coreblas_dparfb(coreblas_enum_t side, coreblas_enum_t trans, coreblas_enum_t direct,
                coreblas_enum_t storev,
                int m1, int n1, int m2, int n2, int k, int l,
                      double *A1,   int lda1,
                      double *A2,   int lda2,
                const double *V,    int ldv,
                const double *T,    int ldt,
                      double *work, int ldwork);

int coreblas_dpemv(coreblas_enum_t trans, int storev,
               int m, int n, int l,
               double alpha,
               const double *A, int lda,
               const double *X, int incx,
               double beta,
               double *Y, int incy,
               double *work);

int coreblas_dpotrf(coreblas_enum_t uplo,
                int n,
                double *A, int lda);

void coreblas_dsymm(coreblas_enum_t side, coreblas_enum_t uplo,
                int m, int n,
                double alpha, const double *A, int lda,
                                          const double *B, int ldb,
                double beta,        double *C, int ldc);

void coreblas_dsyr2k(
    coreblas_enum_t uplo, coreblas_enum_t trans,
    int n, int k,
    double alpha, const double *A, int lda,
                              const double *B, int ldb,
    double beta,        double *C, int ldc);

void coreblas_dsyrk(coreblas_enum_t uplo, coreblas_enum_t trans,
                int n, int k,
                double alpha, const double *A, int lda,
                double beta,        double *C, int ldc);

int coreblas_dtradd(coreblas_enum_t uplo, coreblas_enum_t transa,
                int m, int n,
                double alpha, const double *A, int lda,
                double beta,        double *B, int ldb);

void coreblas_dtrmm(coreblas_enum_t side, coreblas_enum_t uplo,
                coreblas_enum_t transa, coreblas_enum_t diag,
                int m, int n,
                double alpha, const double *A, int lda,
                                                double *B, int ldb);

void coreblas_dtrsm(coreblas_enum_t side, coreblas_enum_t uplo,
                coreblas_enum_t transa, coreblas_enum_t diag,
                int m, int n,
                double alpha, const double *A, int lda,
                                                double *B, int ldb);

void coreblas_dtrssq(coreblas_enum_t uplo, coreblas_enum_t diag,
                 int m, int n,
                 const double *A, int lda,
                 double *scale, double *sumsq);

int coreblas_dtrtri(coreblas_enum_t uplo, coreblas_enum_t diag,
                int n,
                double *A, int lda);

int coreblas_dtslqt(int m, int n, int ib,
                double *A1, int lda1,
                double *A2, int lda2,
                double *T,  int ldt,
                double *tau,
                double *work);

int coreblas_dtsmlq(coreblas_enum_t side, coreblas_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      double *A1,   int lda1,
                      double *A2,   int lda2,
                const double *V,    int ldv,
                const double *T,    int ldt,
                      double *work, int ldwork);

int coreblas_dtsmqr(coreblas_enum_t side, coreblas_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      double *A1,   int lda1,
                      double *A2,   int lda2,
                const double *V,    int ldv,
                const double *T,    int ldt,
                      double *work, int ldwork);

int coreblas_dtsqrt(int m, int n, int ib,
                double *A1, int lda1,
                double *A2, int lda2,
                double *T,  int ldt,
                double *tau,
                double *work);

int coreblas_dttlqt(int m, int n, int ib,
                double *A1, int lda1,
                double *A2, int lda2,
                double *T,  int ldt,
                double *tau,
                double *work);

int coreblas_dttmlq(coreblas_enum_t side, coreblas_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      double *A1,   int lda1,
                      double *A2,   int lda2,
                const double *V,    int ldv,
                const double *T,    int ldt,
                      double *work, int ldwork);

int coreblas_dttmqr(coreblas_enum_t side, coreblas_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      double *A1,   int lda1,
                      double *A2,   int lda2,
                const double *V,    int ldv,
                const double *T,    int ldt,
                      double *work, int ldwork);

int coreblas_dttqrt(int m, int n, int ib,
                double *A1, int lda1,
                double *A2, int lda2,
                double *T,  int ldt,
                double *tau,
                double *work);

int coreblas_dormlq(coreblas_enum_t side, coreblas_enum_t trans,
                int m, int n, int k, int ib,
                const double *A,    int lda,
                const double *T,    int ldt,
                      double *C,    int ldc,
                      double *work, int ldwork);

int coreblas_dormqr(coreblas_enum_t side, coreblas_enum_t trans,
                int m, int n, int k, int ib,
                const double *A,    int lda,
                const double *T,    int ldt,
                      double *C,    int ldc,
                      double *work, int ldwork);

/******************************************************************************/
void coreblas_kernel_damax(int colrow, int m, int n,
                     const double *A, int lda,
                     double *values);

void coreblas_kernel_dgeadd(
    coreblas_enum_t transa, int m, int n,
    double alpha, const double *A, int lda,
    double beta,        double *B, int ldb);

void coreblas_kernel_dgelqt(int m, int n, int ib,
                     double *A, int lda,
                     double *T, int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_dgemm(
    coreblas_enum_t transa, coreblas_enum_t transb,
    int m, int n, int k,
    double alpha, const double *A, int lda,
                              const double *B, int ldb,
    double beta,        double *C, int ldc);

void coreblas_kernel_dgeqrt(int m, int n, int ib,
                     double *A, int lda,
                     double *T, int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_dgessq(int m, int n,
                     const double *A, int lda,
                     double *scale, double *sumsq);

void coreblas_kernel_dgessq_aux(int n,
                         const double *scale, const double *sumsq,
                         double *value);

void coreblas_kernel_dsygst(int itype, coreblas_enum_t uplo,
                     int n,
                     double *A, int lda,
                     double *B, int ldb);

void coreblas_kernel_dsymm(
    coreblas_enum_t side, coreblas_enum_t uplo,
    int m, int n,
    double alpha, const double *A, int lda,
                              const double *B, int ldb,
    double beta,        double *C, int ldc);

void coreblas_kernel_dsyr2k(
    coreblas_enum_t uplo, coreblas_enum_t trans,
    int n, int k,
    double alpha, const double *A, int lda,
                              const double *B, int ldb,
    double beta,                    double *C, int ldc);

void coreblas_kernel_dsyrk(coreblas_enum_t uplo, coreblas_enum_t trans,
                    int n, int k,
                    double alpha, const double *A, int lda,
                    double beta,        double *C, int ldc);

void coreblas_kernel_dsyssq(coreblas_enum_t uplo,
                     int n,
                     const double *A, int lda,
                     double *scale, double *sumsq);

void coreblas_kernel_dsyssq(coreblas_enum_t uplo,
                     int n,
                     const double *A, int lda,
                     double *scale, double *sumsq);

void coreblas_kernel_dsyssq_aux(int m, int n,
                         const double *scale, const double *sumsq,
                         double *value);

void coreblas_kernel_dlacpy(coreblas_enum_t uplo, coreblas_enum_t transa,
                     int m, int n,
                     const double *A, int lda,
                           double *B, int ldb);

void coreblas_kernel_dlacpy_lapack2tile_band(coreblas_enum_t uplo,
                                      int it, int jt,
                                      int m, int n, int nb, int kl, int ku,
                                      const double *A, int lda,
                                            double *B, int ldb);

void coreblas_kernel_dlacpy_tile2lapack_band(coreblas_enum_t uplo,
                                      int it, int jt,
                                      int m, int n, int nb, int kl, int ku,
                                      const double *B, int ldb,
                                            double *A, int lda);

void coreblas_kernel_dlange(coreblas_enum_t norm,
                     int m, int n,
                     const double *A, int lda,
                     double *work, double *result);

void coreblas_kernel_dlange_aux(coreblas_enum_t norm,
                         int m, int n,
                         const double *A, int lda,
                         double *value);

void coreblas_kernel_dlansy(coreblas_enum_t norm, coreblas_enum_t uplo,
                     int n,
                     const double *A, int lda,
                     double *work, double *value);

void coreblas_kernel_dlansy_aux(coreblas_enum_t norm, coreblas_enum_t uplo,
                         int n,
                         const double *A, int lda,
                         double *value);

void coreblas_kernel_dlansy(coreblas_enum_t norm, coreblas_enum_t uplo,
                     int n,
                     const double *A, int lda,
                     double *work, double *value);

void coreblas_kernel_dlansy_aux(coreblas_enum_t norm, coreblas_enum_t uplo,
                         int n,
                         const double *A, int lda,
                         double *value);

void coreblas_kernel_dlantr(coreblas_enum_t norm, coreblas_enum_t uplo, coreblas_enum_t diag,
                     int m, int n,
                     const double *A, int lda,
                     double *work, double *value);

void coreblas_kernel_dlantr_aux(coreblas_enum_t norm, coreblas_enum_t uplo,
                         coreblas_enum_t diag,
                         int m, int n,
                         const double *A, int lda,
                         double *value);

void coreblas_kernel_dlascl(coreblas_enum_t uplo,
                     double cfrom, double cto,
                     int m, int n,
                     double *A, int lda);

void coreblas_kernel_dlaset(coreblas_enum_t uplo,
                     int mb, int nb,
                     int i, int j,
                     int m, int n,
                     double alpha, double beta,
                     double *A);

void coreblas_kernel_dlauum(coreblas_enum_t uplo,
                     int n,
                     double *A, int lda);

void coreblas_kernel_dpotrf(coreblas_enum_t uplo,
                     int n,
                     double *A, int lda,
                     int iinfo);

void coreblas_kernel_dsymm(
    coreblas_enum_t side, coreblas_enum_t uplo,
    int m, int n,
    double alpha, const double *A, int lda,
                              const double *B, int ldb,
    double beta,        double *C, int ldc);

void coreblas_kernel_dsyr2k(
    coreblas_enum_t uplo, coreblas_enum_t trans,
    int n, int k,
    double alpha, const double *A, int lda,
                              const double *B, int ldb,
    double beta,        double *C, int ldc);

void coreblas_kernel_dsyrk(
    coreblas_enum_t uplo, coreblas_enum_t trans,
    int n, int k,
    double alpha, const double *A, int lda,
    double beta,        double *C, int ldc);

void coreblas_kernel_dtradd(
    coreblas_enum_t uplo, coreblas_enum_t transa,
    int m, int n,
    double alpha, const double *A, int lda,
    double beta,        double *B, int ldb);

void coreblas_kernel_dtrmm(
    coreblas_enum_t side, coreblas_enum_t uplo,
    coreblas_enum_t transa, coreblas_enum_t diag,
    int m, int n,
    double alpha, const double *A, int lda,
                                    double *B, int ldb);

void coreblas_kernel_dtrsm(
    coreblas_enum_t side, coreblas_enum_t uplo,
    coreblas_enum_t transa, coreblas_enum_t diag,
    int m, int n,
    double alpha, const double *A, int lda,
                                    double *B, int ldb);

void coreblas_kernel_dtrssq(coreblas_enum_t uplo, coreblas_enum_t diag,
                     int m, int n,
                     const double *A, int lda,
                     double *scale, double *sumsq);

void coreblas_kernel_dtrtri(coreblas_enum_t uplo, coreblas_enum_t diag,
                     int n,
                     double *A, int lda,
                     int iinfo);

void coreblas_kernel_dtslqt(int m, int n, int ib,
                     double *A1, int lda1,
                     double *A2, int lda2,
                     double *T,  int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_dtsmlq(coreblas_enum_t side, coreblas_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           double *A1, int lda1,
                           double *A2, int lda2,
                     const double *V,  int ldv,
                     const double *T,  int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_dtsmqr(coreblas_enum_t side, coreblas_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           double *A1, int lda1,
                           double *A2, int lda2,
                     const double *V, int ldv,
                     const double *T, int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_dtsqrt(int m, int n, int ib,
                     double *A1, int lda1,
                     double *A2, int lda2,
                     double *T,  int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_dttlqt(int m, int n, int ib,
                     double *A1, int lda1,
                     double *A2, int lda2,
                     double *T,  int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_dttmlq(coreblas_enum_t side, coreblas_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           double *A1, int lda1,
                           double *A2, int lda2,
                     const double *V,  int ldv,
                     const double *T,  int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_dttmqr(coreblas_enum_t side, coreblas_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           double *A1, int lda1,
                           double *A2, int lda2,
                     const double *V, int ldv,
                     const double *T, int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_dttqrt(int m, int n, int ib,
                     double *A1, int lda1,
                     double *A2, int lda2,
                     double *T,  int ldt,
                     coreblas_workspace_t work);

void coreblas_kernel_dormlq(coreblas_enum_t side, coreblas_enum_t trans,
                     int m, int n, int k, int ib,
                     const double *A, int lda,
                     const double *T, int ldt,
                           double *C, int ldc,
                     coreblas_workspace_t work);

void coreblas_kernel_dormqr(coreblas_enum_t side, coreblas_enum_t trans,
                     int m, int n, int k, int ib,
                     const double *A, int lda,
                     const double *T, int ldt,
                           double *C, int ldc,
                     coreblas_workspace_t work);

#undef REAL

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // COREBLAS_CORE_BLAS_D_H
