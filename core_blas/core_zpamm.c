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

static inline int coreblas_zpamm_a2(coreblas_enum_t side, coreblas_enum_t trans,
                                coreblas_enum_t uplo,
                                int m, int n, int k, int l, int vi2, int vi3,
                                      coreblas_complex64_t *A2, int lda2,
                                const coreblas_complex64_t *V,  int ldv,
                                      coreblas_complex64_t *W,  int ldw);

static inline int coreblas_zpamm_w(coreblas_enum_t side, coreblas_enum_t trans,
                               coreblas_enum_t uplo,
                               int m, int n, int k, int l, int vi2, int vi3,
                               const coreblas_complex64_t *A1, int lda1,
                                     coreblas_complex64_t *A2, int lda2,
                               const coreblas_complex64_t *V,  int ldv,
                                     coreblas_complex64_t *W,  int ldw);

/***************************************************************************//**
 *
 * @ingroup core_pamm
 *
 *  Performs one of the matrix-matrix operations
 *
 *                    CoreBlasLeft                CoreBlasRight
 *     OP CoreBlasW  :  W  = A1 + op(V) * A2  or  W  = A1 + A2 * op(V)
 *     OP CoreBlasA2 :  A2 = A2 - op(V) * W   or  A2 = A2 - W * op(V)
 *
 *  where  op( V ) is one of
 *
 *     op( V ) = V   or   op( V ) = V^H
 *
 *  A1, A2 and W are general matrices, and V is:
 *
 *        l = k: rectangle + triangle
 *        l < k: rectangle + trapezoid
 *        l = 0: rectangle
 *
 *  Size of V, both rowwise and columnwise, is:
 *
 *         ----------------------
 *          side   trans    size
 *         ----------------------
 *          left     N     M x K
 *                   T     K x M
 *          right    N     K x N
 *                   T     N x K
 *         ----------------------
 *
 *  CoreBlasLeft (columnwise and rowwise):
 *
 *              |    K    |                 |         M         |
 *           _  __________   _              _______________        _
 *              |    |    |                 |             | \
 *     V:       |    |    |           V^H:  |_____________|___\    K
 *              |    |    | M-L             |                  |
 *           M  |    |    |                 |__________________|   _
 *              |____|    |  _
 *              \    |    |                 |    M - L    | L  |
 *                \  |    |  L
 *           _      \|____|  _
 *
 *
 *  CoreBlasRight (columnwise and rowwise):
 *
 *          |         K         |                   |    N    |
 *          _______________        _             _  __________   _
 *          |             | \                       |    |    |
 *    V^H:  |_____________|___\    N        V:      |    |    |
 *          |                  |                    |    |    | K-L
 *          |__________________|   _             K  |    |    |
 *                                                  |____|    |  _
 *          |    K - L    | L  |                    \    |    |
 *                                                    \  |    |  L
 *                                               _      \|____|  _
 *
 *  Arguments
 *  ==========
 *
 * @param[in] op
 *
 *         OP specifies which operation to perform:
 *
 *         - CoreBlasW  : W  = A1 + op(V) * A2  or  W  = A1 + A2 * op(V)
 *         - CoreBlasA2 : A2 = A2 - op(V) * W   or  A2 = A2 - W * op(V)
 *
 * @param[in] side
 *
 *         SIDE specifies whether  op( V ) multiplies A2
 *         or W from the left or right as follows:
 *
 *         - CoreBlasLeft  : multiply op( V ) from the left
 *                            OP CoreBlasW  :  W  = A1 + op(V) * A2
 *                            OP CoreBlasA2 :  A2 = A2 - op(V) * W
 *
 *         - CoreBlasRight : multiply op( V ) from the right
 *                            OP CoreBlasW  :  W  = A1 + A2 * op(V)
 *                            OP CoreBlasA2 :  A2 = A2 - W * op(V)
 *
 * @param[in] storev
 *
 *         Indicates how the vectors which define the elementary
 *         reflectors are stored in V:
 *
 *         - CoreBlasColumnwise
 *         - CoreBlasRowwise
 *
 * @param[in] m
 *         The number of rows of the A1, A2 and W
 *         If SIDE is CoreBlasLeft, the number of rows of op( V )
 *
 * @param[in] n
 *         The number of columns of the A1, A2 and W
 *         If SIDE is CoreBlasRight, the number of columns of op( V )
 *
 * @param[in] k
 *         If SIDE is CoreBlasLeft, the number of columns of op( V )
 *         If SIDE is CoreBlasRight, the number of rows of op( V )
 *
 * @param[in] l
 *         The size of the triangular part of V
 *
 * @param[in] A1
 *         On entry, the m-by-n tile A1.
 *
 * @param[in] lda1
 *         The leading dimension of the array A1. lda1 >= max(1,m).
 *
 * @param[in,out] A2
 *         On entry, the m-by-n tile A2.
 *         On exit, if OP is CoreBlasA2 A2 is overwritten
 *
 * @param[in] lda2
 *         The leading dimension of the tile A2. lda2 >= max(1,m).
 *
 * @param[in] V
 *         The matrix V as described above.
 *         If SIDE is CoreBlasLeft : op( V ) is m-by-k
 *         If SIDE is CoreBlasRight: op( V ) is k-by-n
 *
 * @param[in] ldv
 *         The leading dimension of the array V.
 *
 * @param[in,out] W
 *         On entry, the m-by-n matrix W.
 *         On exit, W is overwritten either if OP is CoreBlasA2 or CoreBlasW.
 *         If OP is CoreBlasA2, W is an input and is used as a workspace.
 *
 * @param[in] ldw
 *         The leading dimension of array W.
 *
 *******************************************************************************
 *
 * @retval CoreBlasSuccess successful exit
 * @retval < 0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
__attribute__((weak))
int coreblas_zpamm(coreblas_enum_t op, coreblas_enum_t side, coreblas_enum_t storev,
               int m, int n, int k, int l,
               const coreblas_complex64_t *A1, int lda1,
                     coreblas_complex64_t *A2, int lda2,
               const coreblas_complex64_t *V,  int ldv,
                     coreblas_complex64_t *W,  int ldw)
{
    // Check input arguments.
    if ((op != CoreBlasW) && (op != CoreBlasA2)) {
        coreblas_error("illegal value of op");
        return -1;
    }
    if ((side != CoreBlasLeft) && (side != CoreBlasRight)) {
        coreblas_error("illegal value of side");
        return -2;
    }
    if ((storev != CoreBlasColumnwise) && (storev != CoreBlasRowwise)) {
        coreblas_error("illegal value of storev");
        return -3;
    }
    if (m < 0) {
        coreblas_error("illegal value of m");
        return -4;
    }
    if (n < 0) {
        coreblas_error("illegal value of n");
        return -5;
    }
    if (k < 0) {
        coreblas_error("illegal value of k");
        return -6;
    }
    if (l < 0) {
        coreblas_error("illegal value of l");
        return -7;
    }
    if (A1 == NULL) {
        coreblas_error("NULL A1");
        return -8;
    }
    if (lda1 < 0) {
        coreblas_error("illegal value of lda1");
        return -9;
    }
    if (A2 == NULL) {
        coreblas_error("NULL A2");
        return -10;
    }
    if (lda2 < 0) {
        coreblas_error("illegal value of lda2");
        return -11;
    }
    if (V == NULL) {
        coreblas_error("NULL V");
        return -12;
    }
    if (ldv < 0) {
        coreblas_error("illegal value of ldv");
        return -13;
    }
    if (W == NULL) {
        coreblas_error("NULL W");
        return -14;
    }
    if (ldw < 0) {
        coreblas_error("illegal value of ldw");
        return -15;
    }

    // quick return
    if (m == 0 || n == 0 || k == 0)
        return CoreBlasSuccess;

    // TRANS is set as:
    //
    //        -------------------------------------
    //         side   direct     CoreBlasW  CoreBlasA2
    //        -------------------------------------
    //         left   colwise       T        N
    //                rowwise       N        T
    //         right  colwise       N        T
    //                rowwise       T        N
    //        -------------------------------------

    coreblas_enum_t uplo;
    coreblas_enum_t trans;
    int vi2, vi3;

    //===================
    // CoreBlasColumnwise
    //===================
    if (storev == CoreBlasColumnwise) {
        uplo = CoreBlasUpper;
        if (side == CoreBlasLeft) {
            trans = op == CoreBlasA2 ? CoreBlasNoTrans : CoreBlas_ConjTrans;
            vi2 = trans == CoreBlasNoTrans ? m - l : k - l;
        }
        else {
            trans = op == CoreBlasW ? CoreBlasNoTrans : CoreBlas_ConjTrans;
            vi2 = trans == CoreBlasNoTrans ? k - l : n - l;
        }
        vi3 = ldv * l;
    }
    //================
    // CoreBlasRowwise
    //================
    else {
        uplo = CoreBlasLower;
        if (side == CoreBlasLeft) {
            trans = op == CoreBlasW ? CoreBlasNoTrans : CoreBlas_ConjTrans;
            vi2 = trans == CoreBlasNoTrans ? k - l : m - l;
        }
        else {
            trans = op == CoreBlasA2 ? CoreBlasNoTrans : CoreBlas_ConjTrans;
            vi2 = trans == CoreBlasNoTrans ? n - l : k - l;
        }
        vi2 *= ldv;
        vi3  = l;
    }

    if (op == CoreBlasW) {
        coreblas_zpamm_w(side, trans, uplo,
                     m, n, k, l, vi2, vi3,
                     A1, lda1,
                     A2, lda2,
                     V, ldv,
                     W, ldw);
    }
    else if (op == CoreBlasA2) {
        coreblas_zpamm_a2(side, trans, uplo,
                      m, n, k, l, vi2, vi3,
                      A2, lda2,
                      V,  ldv,
                      W,  ldw);
    }

    return CoreBlasSuccess;
}

/******************************************************************************/
static inline int coreblas_zpamm_w(
        coreblas_enum_t side, coreblas_enum_t trans, coreblas_enum_t uplo,
        int m, int n, int k, int l, int vi2, int vi3,
        const coreblas_complex64_t *A1, int lda1,
              coreblas_complex64_t *A2, int lda2,
        const coreblas_complex64_t *V,  int ldv,
              coreblas_complex64_t *W,  int ldw)
{
    // W = A1 + op(V) * A2  or  W = A1 + A2 * op(V)

    coreblas_complex64_t zone  = 1.0;
    coreblas_complex64_t zzero = 0.0;

    //=============
    // CoreBlasLeft
    //=============
    if (side == CoreBlasLeft) {
        if ((trans == CoreBlas_ConjTrans && uplo == CoreBlasUpper) ||
            (trans == CoreBlasNoTrans && uplo == CoreBlasLower)) {
            // W = A1 + V^H * A2

            // W = A2_2
            #ifdef COREBLAS_USE_64BIT_BLAS
                 LAPACKE_zlacpy64_(LAPACK_COL_MAJOR,
                            lapack_const(CoreBlasGeneral),
                            l, n,
                            &A2[k-l], lda2,
                            W,       ldw);
            #else
                 LAPACKE_zlacpy(LAPACK_COL_MAJOR,
                            lapack_const(CoreBlasGeneral),
                            l, n,
                            &A2[k-l], lda2,
                            W,       ldw);
            #endif 


            // W = V_2^H * W + V_1^H * A2_1 (ge+tr, top L rows of V^H)
            if (l > 0) {
                // W = V_2^H * W
            #ifdef COREBLAS_USE_64BIT_BLAS
                cblas_ztrmm64_(CblasColMajor,
                            CblasLeft, (CBLAS_UPLO)uplo,
                            (CBLAS_TRANSPOSE)trans, CblasNonUnit,
                            l, n,
                            CBLAS_SADDR(zone), &V[vi2], ldv,
                                                W,      ldw);
            #else
                cblas_ztrmm(CblasColMajor,
                            CblasLeft, (CBLAS_UPLO)uplo,
                            (CBLAS_TRANSPOSE)trans, CblasNonUnit,
                            l, n,
                            CBLAS_SADDR(zone), &V[vi2], ldv,
                                                W,      ldw);
            #endif 


                // W = W + V_1^H * A2_1
                if (k > l) {
                #ifdef COREBLAS_USE_64BIT_BLAS
                    cblas_zgemm64_(CblasColMajor,
                                (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                                l, n, k-l,
                                CBLAS_SADDR(zone), V,  ldv,
                                                   A2, lda2,
                                CBLAS_SADDR(zone), W,  ldw);

                #else
                    cblas_zgemm(CblasColMajor,
                                (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                                l, n, k-l,
                                CBLAS_SADDR(zone), V,  ldv,
                                                   A2, lda2,
                                CBLAS_SADDR(zone), W,  ldw);
                #endif 

                }
            }

            // W_2 = V_3^H * A2: (ge, bottom M-L rows of V^H)
            if (m > l) {
            #ifdef COREBLAS_USE_64BIT_BLAS
                cblas_zgemm64_(CblasColMajor,
                            (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                            (m-l), n, k,
                            CBLAS_SADDR(zone),  &V[vi3], ldv,
                                                 A2,     lda2,
                            CBLAS_SADDR(zzero), &W[l],   ldw);
            #else
                cblas_zgemm(CblasColMajor,
                            (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                            (m-l), n, k,
                            CBLAS_SADDR(zone),  &V[vi3], ldv,
                                                 A2,     lda2,
                            CBLAS_SADDR(zzero), &W[l],   ldw);
            #endif 

            }

            // W = A1 + W
            for (int j = 0; j < n; j++) {
            #ifdef COREBLAS_USE_64BIT_BLAS
                cblas_zaxpy64_(m, CBLAS_SADDR(zone),
                            &A1[lda1*j], 1,
                            &W[ldw*j], 1);
            #else
                cblas_zaxpy(m, CBLAS_SADDR(zone),
                            &A1[lda1*j], 1,
                            &W[ldw*j], 1);
            #endif 

            }
        }
        else {
            coreblas_error(
                "Left Upper/NoTrans & Lower/[Conj]Trans not implemented");
            return CoreBlasErrorNotSupported;
        }
    }
    //==============
    // CoreBlasRight
    //==============
    else {
        if ((trans == CoreBlas_ConjTrans && uplo == CoreBlasUpper) ||
            (trans == CoreBlasNoTrans && uplo == CoreBlasLower)) {
            coreblas_error(
                "Right Upper/[Conj]Trans & Lower/NoTrans not implemented");
            return CoreBlasErrorNotSupported;
        }
        else {
            // W = A1 + A2 * V
            if (l > 0) {
                // W = A2_2
                #ifdef COREBLAS_USE_64BIT_BLAS
                    LAPACKE_zlacpy64_(LAPACK_COL_MAJOR,
                                lapack_const(CoreBlasGeneral),
                                m, l,
                                &A2[lda2*(k-l)], lda2,
                                W,              ldw);
                #else
                    LAPACKE_zlacpy(LAPACK_COL_MAJOR,
                                lapack_const(CoreBlasGeneral),
                                m, l,
                                &A2[lda2*(k-l)], lda2,
                                W,              ldw);
                #endif 

                // W = W * V_2 --> W = A2_2 * V_2
                #ifdef COREBLAS_USE_64BIT_BLAS
                    cblas_ztrmm64_(CblasColMajor,
                            CblasRight, (CBLAS_UPLO)uplo,
                            (CBLAS_TRANSPOSE)trans, CblasNonUnit,
                            m, l,
                            CBLAS_SADDR(zone), &V[vi2], ldv,
                                                W,      ldw);
                #else
                    cblas_ztrmm(CblasColMajor,
                            CblasRight, (CBLAS_UPLO)uplo,
                            (CBLAS_TRANSPOSE)trans, CblasNonUnit,
                            m, l,
                            CBLAS_SADDR(zone), &V[vi2], ldv,
                                                W,      ldw);
                #endif 


                // W = W + A2_1 * V_1
                if (k > l) {
                #ifdef COREBLAS_USE_64BIT_BLAS
                        cblas_zgemm64_(CblasColMajor,
                                CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                                m, l, k-l,
                                CBLAS_SADDR(zone), A2, lda2,
                                                   V,  ldv,
                                CBLAS_SADDR(zone), W,  ldw);
                #else
                        cblas_zgemm(CblasColMajor,
                                CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                                m, l, k-l,
                                CBLAS_SADDR(zone), A2, lda2,
                                                   V,  ldv,
                                CBLAS_SADDR(zone), W,  ldw);
                #endif 
                }
            }

            // W = W + A2 * V_3
            if (n > l) {
            #ifdef COREBLAS_USE_64BIT_BLAS
                    cblas_zgemm64_(CblasColMajor,
                            CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                            m, n-l, k,
                            CBLAS_SADDR(zone),   A2,       lda2,
                                                &V[vi3],   ldv,
                            CBLAS_SADDR(zzero), &W[ldw*l], ldw);
            #else
                    cblas_zgemm(CblasColMajor,
                            CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                            m, n-l, k,
                            CBLAS_SADDR(zone),   A2,       lda2,
                                                &V[vi3],   ldv,
                            CBLAS_SADDR(zzero), &W[ldw*l], ldw);
            #endif 

            }

            // W = A1 + W
            for (int j = 0; j < n; j++) {
            #ifdef COREBLAS_USE_64BIT_BLAS
                    cblas_zaxpy64_(m, CBLAS_SADDR(zone),
                            &A1[lda1*j], 1,
                            &W[ldw*j],   1);
            #else
                    cblas_zaxpy(m, CBLAS_SADDR(zone),
                            &A1[lda1*j], 1,
                            &W[ldw*j],   1);
            #endif 

            }
        }
    }

    return CoreBlasSuccess;
}

/******************************************************************************/
static inline int coreblas_zpamm_a2(
        coreblas_enum_t side, coreblas_enum_t trans, coreblas_enum_t uplo,
        int m, int n, int k, int l, int vi2, int vi3,
              coreblas_complex64_t *A2, int lda2,
        const coreblas_complex64_t *V,  int ldv,
              coreblas_complex64_t *W,  int ldw)
{
    // A2 = A2 + op(V) * W  or  A2 = A2 + W * op(V)

    coreblas_complex64_t zone  =  1.0;
    coreblas_complex64_t zmone = -1.0;

    //=============
    // CoreBlasLeft
    //=============
    if (side == CoreBlasLeft) {
        if ((trans == CoreBlas_ConjTrans && uplo == CoreBlasUpper) ||
            (trans == CoreBlasNoTrans && uplo == CoreBlasLower)) {
            coreblas_error(
                "Left Upper/[Conj]Trans & Lower/NoTrans not implemented");
            return CoreBlasErrorNotSupported;
        }
        else {
            // A2 = A2 - V * W

            // A2_1 = A2_1 - V_1  * W_1
            if (m > l) {
            #ifdef COREBLAS_USE_64BIT_BLAS
                cblas_zgemm64_(CblasColMajor,
                            (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                            m-l, n, l,
                            CBLAS_SADDR(zmone), V,  ldv,
                                                W,  ldw,
                            CBLAS_SADDR(zone),  A2, lda2);
            #else
                cblas_zgemm(CblasColMajor,
                            (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                            m-l, n, l,
                            CBLAS_SADDR(zmone), V,  ldv,
                                                W,  ldw,
                            CBLAS_SADDR(zone),  A2, lda2);
            #endif 

            }

            // W_1 = V_2 * W_1
            #ifdef COREBLAS_USE_64BIT_BLAS
                cblas_ztrmm64_(CblasColMajor,
                        CblasLeft, (CBLAS_UPLO)uplo,
                        (CBLAS_TRANSPOSE)trans, CblasNonUnit,
                        l, n,
                        CBLAS_SADDR(zone), &V[vi2], ldv,
                                            W,      ldw);
            #else
                cblas_ztrmm(CblasColMajor,
                        CblasLeft, (CBLAS_UPLO)uplo,
                        (CBLAS_TRANSPOSE)trans, CblasNonUnit,
                        l, n,
                        CBLAS_SADDR(zone), &V[vi2], ldv,
                                            W,      ldw);
            #endif 


            // A2_2 = A2_2 - W_1
            for (int j = 0; j < n; j++) {
            #ifdef COREBLAS_USE_64BIT_BLAS
                cblas_zaxpy64_(l, CBLAS_SADDR(zmone),
                            &W[ldw*j], 1,
                            &A2[lda2*j+(m-l)], 1);
            #else
                cblas_zaxpy(l, CBLAS_SADDR(zmone),
                            &W[ldw*j], 1,
                            &A2[lda2*j+(m-l)], 1);
            #endif 

            }

            // A2 = A2 - V_3  * W_2
            if (k > l) {
                    #ifdef COREBLAS_USE_64BIT_BLAS
                        cblas_zgemm64_(CblasColMajor,
                            (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                            m, n, (k-l),
                            CBLAS_SADDR(zmone), &V[vi3], ldv,
                                                &W[l],   ldw,
                            CBLAS_SADDR(zone),   A2,     lda2);
                    #else
                        cblas_zgemm(CblasColMajor,
                            (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                            m, n, (k-l),
                            CBLAS_SADDR(zmone), &V[vi3], ldv,
                                                &W[l],   ldw,
                            CBLAS_SADDR(zone),   A2,     lda2);
                    #endif 

            }
        }
    }
    //==============
    // CoreBlasRight
    //==============
    else {
        if ((trans == CoreBlas_ConjTrans && uplo == CoreBlasUpper) ||
            (trans == CoreBlasNoTrans && uplo == CoreBlasLower)) {
            // A2 = A2 - W * V^H

            // A2 = A2 - W_2 * V_3^H
            if (k > l) {
            #ifdef COREBLAS_USE_64BIT_BLAS
                cblas_zgemm64_(CblasColMajor,
                            CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                            m, n, k-l,
                            CBLAS_SADDR(zmone), &W[ldw*l], ldw,
                                                &V[vi3],   ldv,
                            CBLAS_SADDR(zone),   A2,       lda2);
            #else
                cblas_zgemm(CblasColMajor,
                            CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                            m, n, k-l,
                            CBLAS_SADDR(zmone), &W[ldw*l], ldw,
                                                &V[vi3],   ldv,
                            CBLAS_SADDR(zone),   A2,       lda2);
            #endif 

            }

            // A2_1 = A2_1 - W_1 * V_1^H
            if (n > l) {
            #ifdef COREBLAS_USE_64BIT_BLAS
                cblas_zgemm64_(CblasColMajor,
                            CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                            m, n-l, l,
                            CBLAS_SADDR(zmone), W,  ldw,
                                                V,  ldv,
                            CBLAS_SADDR(zone),  A2, lda2);
            #else
                cblas_zgemm(CblasColMajor,
                            CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                            m, n-l, l,
                            CBLAS_SADDR(zmone), W,  ldw,
                                                V,  ldv,
                            CBLAS_SADDR(zone),  A2, lda2);
            #endif 

            }

            // A2_2 =  A2_2 -  W_1 * V_2^H
            if (l > 0) {
            #ifdef COREBLAS_USE_64BIT_BLAS
                cblas_ztrmm64_(CblasColMajor,
                            CblasRight, (CBLAS_UPLO)uplo,
                            (CBLAS_TRANSPOSE)trans, CblasNonUnit,
                            m, l,
                            CBLAS_SADDR(zmone), &V[vi2], ldv,
                                                 W,      ldw);
            #else
                cblas_ztrmm(CblasColMajor,
                            CblasRight, (CBLAS_UPLO)uplo,
                            (CBLAS_TRANSPOSE)trans, CblasNonUnit,
                            m, l,
                            CBLAS_SADDR(zmone), &V[vi2], ldv,
                                                 W,      ldw);
            #endif 


                for (int j = 0; j < l; j++) {
                #ifdef COREBLAS_USE_64BIT_BLAS
                    cblas_zaxpy64_(m, CBLAS_SADDR(zone),
                                &W[ldw*j], 1,
                                &A2[lda2*(n-l+j)], 1);
                #else
                    cblas_zaxpy(m, CBLAS_SADDR(zone),
                                &W[ldw*j], 1,
                                &A2[lda2*(n-l+j)], 1);
                #endif 

                }
            }
        }
        else {
            coreblas_error(
                "Right Upper/NoTrans & Lower/[Conj]Trans not implemented");
            return CoreBlasErrorNotSupported;
        }
    }

    return CoreBlasSuccess;
}