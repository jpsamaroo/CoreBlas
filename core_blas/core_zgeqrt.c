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
 * @ingroup core_geqrt
 *
 *  Computes a QR factorization of an m-by-n tile A:
 *  The factorization has the form
 *    \f[
 *        A = Q \times R
 *    \f]
 *  The tile Q is represented as a product of elementary reflectors
 *    \f[
 *        Q = H(1) H(2) ... H(k),
 *    \f]
 *  where \f$ k = min(m,n) \f$.
 *
 *  Each \f$ H(i) \f$ has the form
 *    \f[
 *        H(i) = I - \tau \times v \times v^H
 *    \f]
 *  where \f$ tau \f$ is a scalar, and \f$ v \f$ is a vector with
 *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
 *  and \f$ tau \f$ in tau(i).
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows of the tile A.  m >= 0.
 *
 * @param[in] n
 *         The number of columns of the tile A.  n >= 0.
 *
 * @param[in] ib
 *         The inner-blocking size.  ib >= 0.
 *
 * @param[in,out] A
 *         On entry, the m-by-n tile A.
 *         On exit, the elements on and above the diagonal of the array
 *         contain the min(m,n)-by-n upper trapezoidal tile R (R is
 *         upper triangular if m >= n); the elements below the diagonal,
 *         with the array tau, represent the unitary tile Q as a
 *         product of elementary reflectors (see Further Details).
 *
 * @param[in] lda
 *         The leading dimension of the array A.  lda >= max(1,m).
 *
 * @param[out] T
 *         The ib-by-n triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] ldt
 *         The leading dimension of the array T. ldt >= ib.
 *
 * @param tau
 *         Auxiliary workspace array of length n.
 *
 * @param work
 *         Auxiliary workspace array of length ib*n.
 *
 * @param[in] lwork
 *         Size of the array work. Should be at least ib*n.
 *
 *******************************************************************************
 *
 * @retval CoreBlasSuccess successful exit
 * @retval < 0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
__attribute__((weak))
int coreblas_zgeqrt(int m, int n, int ib,
                coreblas_complex64_t *A, int lda,
                coreblas_complex64_t *T, int ldt,
                coreblas_complex64_t *tau,
                coreblas_complex64_t *work)
{
    // Check input arguments.
    if (m < 0) {
        coreblas_error("illegal value of m");
        return -1;
    }
    if (n < 0) {
        coreblas_error("illegal value of n");
        return -2;
    }
    if ((ib < 0) || ( (ib == 0) && (m > 0) && (n > 0) )) {
        coreblas_error("illegal value of ib");
        return -3;
    }
    if (A == NULL) {
        coreblas_error("NULL A");
        return -4;
    }
    if (lda < imax(1, m) && m > 0) {
        coreblas_error("illegal value of lda");
        return -5;
    }
    if (T == NULL) {
        coreblas_error("NULL T");
        return -6;
    }
    if (ldt < imax(1, ib) && ib > 0) {
        coreblas_error("illegal value of ldt");
        return -7;
    }
    if (tau == NULL) {
        coreblas_error("NULL tau");
        return -8;
    }
    if (work == NULL) {
        coreblas_error("NULL work");
        return -9;
    }

    // quick return
    if (m == 0 || n == 0 || ib == 0)
        return CoreBlasSuccess;

    int k = imin(m, n);
    for (int i = 0; i < k; i += ib) {
        int sb = imin(ib, k-i);
        #ifdef COREBLAS_USE_64BIT_BLAS
            LAPACKE_zgeqr264_(LAPACK_COL_MAJOR,
                        m-i, sb,
                        &A[lda*i+i], lda,
                        &tau[i]);
        #else
            LAPACKE_zgeqr2_work(LAPACK_COL_MAJOR,
                        m-i, sb,
                        &A[lda*i+i], lda,
                        &tau[i], work);
        #endif 

        #ifdef COREBLAS_USE_64BIT_BLAS
            printf("[i:%d] Calling zlarft64\n", i);
            LAPACKE_zlarft_work64_(LAPACK_COL_MAJOR,
                       lapack_const(CoreBlasForward),
                       lapack_const(CoreBlasColumnwise),
                       m-i, sb,
                       &A[lda*i+i], lda,
                       &tau[i],
                       &T[ldt*i], ldt);
        #else
            LAPACKE_zlarft_work(LAPACK_COL_MAJOR,
                        lapack_const(CoreBlasForward),
                        lapack_const(CoreBlasColumnwise),
                        m-i, sb,
                        &A[lda*i+i], lda,
                        &tau[i],
                        &T[ldt*i], ldt);
        #endif 


        if (n > i+sb) {
        #ifdef COREBLAS_USE_64BIT_BLAS
            LAPACKE_zlarfb64_(LAPACK_COL_MAJOR,
                                lapack_const(CoreBlasLeft),
                                lapack_const(CoreBlas_ConjTrans),
                                lapack_const(CoreBlasForward),
                                lapack_const(CoreBlasColumnwise),
                                m-i, n-i-sb, sb,
                                &A[lda*i+i],      lda,
                                &T[ldt*i],        ldt,
                                &A[lda*(i+sb)+i], lda);
        #else
            LAPACKE_zlarfb_work(LAPACK_COL_MAJOR,
                                lapack_const(CoreBlasLeft),
                                lapack_const(CoreBlas_ConjTrans),
                                lapack_const(CoreBlasForward),
                                lapack_const(CoreBlasColumnwise),
                                m-i, n-i-sb, sb,
                                &A[lda*i+i],      lda,
                                &T[ldt*i],        ldt,
                                &A[lda*(i+sb)+i], lda,
                                work, n-i-sb);
        #endif 

        }
    }

    return CoreBlasSuccess;
}
