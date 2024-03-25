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
#include "coreblas_internal.h"
#include "core_lapack.h"

/***************************************************************************//**
 *
 * @ingroup core_lacpy
 *
 *  Copies all or part of a two-dimensional matrix A to another matrix B.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - CoreBlasGeneral: entire A,
 *          - CoreBlasUpper:   upper triangle,
 *          - CoreBlasLower:   lower triangle.
 *
 * @param[in] transa
 *          - CoreBlasNoTrans:   A is not transposed,
 *          - CoreBlasTrans:     A is transposed,
 *          - CoreBlasConjTrans: A is conjugate transposed.
 *
 * @param[in] m
 *          The number of rows of the matrices A and B.
 *          m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrices A and B.
 *          n >= 0.
 *
 * @param[in] A
 *          The m-by-n matrix to copy.
 *
 * @param[in] lda
 *          The leading dimension of the array A.
 *          lda >= max(1,m).
 *
 * @param[out] B
 *          The m-by-n copy of the matrix A.
 *          On exit, B = A ONLY in the locations specified by uplo.
 *
 * @param[in] ldb
 *          The leading dimension of the array B.
 *          ldb >= max(1,m).
 *
 ******************************************************************************/
__attribute__((weak))
void coreblas_zlacpy(coreblas_enum_t uplo, coreblas_enum_t transa,
                 int m, int n,
                 const coreblas_complex64_t *A, int lda,
                       coreblas_complex64_t *B, int ldb)
{
    if (transa == CoreBlasNoTrans) {
        #ifdef COREBLAS_USE_64BIT_BLAS
            LAPACKE_zlacpy_work64_(LAPACK_COL_MAJOR,
                            lapack_const(uplo),
                            m, n,
                            A, lda,
                            B, ldb);
        #else
            LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,
                            lapack_const(uplo),
                            m, n,
                            A, lda,
                            B, ldb);
        #endif

    }
    else if (transa == CoreBlasTrans) {
        switch (uplo) {
        case CoreBlasUpper:
            for (int i = 0; i < imin(m, n); i++)
                for (int j = i; j < n; j++)
                    B[j + i*ldb] = A[i + j*lda];
            break;
        case CoreBlasLower:
            for (int i = 0; i < m; i++)
                for (int j = 0; j <= imin(i, n); j++)
                    B[j + i*ldb] = A[i + j*lda];
            break;
        case CoreBlasGeneral:
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    B[j + i*ldb] = A[i + j*lda];
            break;
        }
    }
    else {
        switch (uplo) {
        case CoreBlasUpper:
            for (int i = 0; i < imin(m, n); i++)
                for (int j = i; j < n; j++)
                    B[j + i*ldb] = conj(A[i + j*lda]);
            break;
        case CoreBlasLower:
            for (int i = 0; i < m; i++)
                for (int j = 0; j <= imin(i, n); j++)
                    B[j + i*ldb] = conj(A[i + j*lda]);
            break;
        case CoreBlasGeneral:
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    B[j + i*ldb] = conj(A[i + j*lda]);
            break;
        }
    }
}