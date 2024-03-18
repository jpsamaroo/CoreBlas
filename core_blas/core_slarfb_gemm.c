/**
 *
 * @file
 *
 *  coreblas is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zlarfb_gemm.c, normal z -> s, Mon Mar 18 06:35:41 2024
 *
 **/

#include "coreblas.h"
#include "coreblas_types.h"
#include "core_lapack.h"

/***************************************************************************//**
 *
 * @ingroup CORE_coreblas_Complex64_t
 *
 *  CORE_slarfb_gemm applies a complex block reflector H or its transpose H'
 *  to a complex M-by-N matrix C, from either the left or the right.
 *  this kernel is similar to the lapack slarfb but it do a full gemm on the
 *  triangular Vs assuming that the upper part of Vs is zero and ones are on
 *  the diagonal. It is also based on the fact that a gemm on a small block of
 *  k reflectors is faster than a trmm on the triangular (k,k) + gemm below.
 *
 *  NOTE THAT:
 *  Only Columnwise/Forward cases are treated here.
 *
 *******************************************************************************
 *
 * @param[in] side
 *         @arg CoreBlasLeft  : apply Q or Q**H from the Left;
 *         @arg CoreBlasRight : apply Q or Q**H from the Right.
 *
 * @param[in] trans
 *         @arg CoreBlasNoTrans   : No transpose, apply Q;
 *         @arg CoreBlasConjTrans : ConjTranspose, apply Q**H.
 *
 * @param[in] direct
 *         Indicates how H is formed from a product of elementary
 *         reflectors
 *         @arg CoreBlasForward  : H = H(1) H(2) . . . H(k) (Forward)
 *         @arg CoreBlasBackward : H = H(k) . . . H(2) H(1) (Backward)
 *
 * @param[in] storev
 *         Indicates how the vectors which define the elementary
 *         reflectors are stored:
 *         @arg CoreBlasColumnwise
 *         @arg CoreBlasRowwise
 *
 * @param[in] M
 *         The number of rows of the matrix C.
 *
 * @param[in] N
 *         The number of columns of the matrix C.
 *
 * @param[in] K
 *          The order of the matrix T (= the number of elementary
 *          reflectors whose product defines the block reflector).
 *
 * @param[in] V
 *          COMPLEX*16 array, dimension
 *              (LDV,K) if storev = 'C'
 *              (LDV,M) if storev = 'R' and side = 'L'
 *              (LDV,N) if storev = 'R' and side = 'R'
 *          The matrix V. See further details.
 *
 * @param[in] LDV
 *         The leading dimension of the array V.
 *         If storev = 'C' and side = 'L', LDV >= max(1,M);
 *         if storev = 'C' and side = 'R', LDV >= max(1,N);
 *         if storev = 'R', LDV >= K.
 *
 * @param[in] T
 *         The triangular K-by-K matrix T in the representation of the
 *         block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= K.
 *
 * @param[in,out] C
 *         COMPLEX*16 array, dimension (LDC,N)
 *         On entry, the M-by-N matrix C.
 *         On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
 *
 * @param[in] LDC
 *         The leading dimension of the array C. LDC >= max(1,M).
 *
 * @param[in,out] WORK
 *         (workspace) COMPLEX*16 array, dimension (LDWORK,K).
 *
 * @param[in] LDWORK
 *         The dimension of the array WORK.
 *         If side = CoreBlasLeft,  LDWORK >= max(1,N);
 *         if side = CoreBlasRight, LDWORK >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval coreblas_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(coreblas_HAVE_WEAK)
#pragma weak CORE_slarfb_gemm = PCORE_slarfb_gemm
#define CORE_slarfb_gemm PCORE_slarfb_gemm
#endif
int coreblas_slarfb_gemm(coreblas_enum_t side, coreblas_enum_t trans, int direct, int storev,
                     int M, int N, int K,
                     const float *V, int LDV,
                     const float *T, int LDT,
                     float *C, int LDC,
                     float *WORK, int LDWORK)
{
    static float zzero =  0.0;
    static float zone  =  1.0;
    static float mzone = -1.0;

    // Quick return //
    if ((M == 0) || (N == 0) || (K == 0) )
        return CoreBlasSuccess;

    /* For Left case, switch the trans. noswitch for right case */
    if( side == CoreBlasLeft){
        if ( trans == CoreBlasNoTrans) {
            trans = CoreBlasConjTrans;
        }
        else {
            trans = CoreBlasNoTrans;
        }
    }

    // main code //
    if (storev == CoreBlasColumnwise )
    {
        if ( direct == CoreBlasForward )
        {
            /*
             * Let  V =  ( V1 )    (first K rows are unit Lower triangular)
             *           ( V2 )
             */
            if ( side == CoreBlasLeft )
            {
                /*
                 * Columnwise / Forward / Left
                 */
                /*
                 * Form  H * C  or  H' * C  where  C = ( C1 )
                 *                                     ( C2 )
                 *
                 * W := C' * V    (stored in WORK)
                 */
                cblas_sgemm(
                    CblasColMajor, CblasConjTrans, CblasNoTrans,
                    N, K, M,
                    (zone), C, LDC,
                    V, LDV,
                    (zzero), WORK, LDWORK );
                /*
                 * W := W * T'  or  W * T
                 */
                cblas_strmm(
                    CblasColMajor, CblasRight, CblasUpper,
                    (CBLAS_TRANSPOSE)trans, CblasNonUnit, N, K,
                    (zone), T, LDT, WORK, LDWORK );
                /*
                 * C := C - V * W'
                 */
                cblas_sgemm(
                    CblasColMajor, CblasNoTrans, CblasConjTrans,
                    M, N, K,
                    (mzone), V, LDV,
                    WORK, LDWORK,
                    (zone), C, LDC);
            }
            else {
                /*
                 * Columnwise / Forward / Right
                 */
                /*
                 * Form  C * H  or  C * H'  where  C = ( C1  C2 )
                 * W := C * V
                 */
                cblas_sgemm(
                    CblasColMajor, CblasNoTrans, CblasNoTrans,
                    M, K, N,
                    (zone), C, LDC,
                    V, LDV,
                    (zzero), WORK, LDWORK);
                /*
                 * W := W * T  or  W * T'
                 */
                cblas_strmm(
                    CblasColMajor, CblasRight, CblasUpper,
                    (CBLAS_TRANSPOSE)trans, CblasNonUnit, M, K,
                    (zone), T, LDT, WORK, LDWORK );
                /*
                 * C := C - W * V'
                 */
                cblas_sgemm(
                    CblasColMajor, CblasNoTrans, CblasConjTrans,
                    M, N, K,
                    (mzone), WORK, LDWORK,
                    V, LDV,
                    (zone), C, LDC);
            }
        }
        else {
            /*
             * Columnwise / Backward / Left or Right
             */
            //coreblas_error(3, "Not implemented (ColMajor / Backward / Left or Right)");
            return CoreBlasErrorNotSupported;
        }
    }
    else {
        if (direct == CoreBlasForward) {
            /*
             * Rowwise / Forward / Left or Right
             */
            //coreblas_error(3, "Not implemented (RowMajor / Backward / Left or Right)");
            return CoreBlasErrorNotSupported;;
        }
        else {
            /*
             * Rowwise / Backward / Left or Right
             */
            //coreblas_error(3, "Not implemented (RowMajor / Backward / Left or Right)");
            return CoreBlasErrorNotSupported;
        }
    }
    return CoreBlasSuccess;
}
