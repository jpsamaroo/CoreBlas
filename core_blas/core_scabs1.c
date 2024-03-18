/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_dcabs1.c, normal z -> c, Mon Mar 18 06:35:41 2024
 *
 **/

#include <coreblas.h>

#include <math.h>

/***************************************************************************//**
 *
 * @ingroup core_cabs1
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 *******************************************************************************
 *
 * @retval Complex 1-norm absolute value: abs(real(alpha)) + abs(imag(alpha)).
 *
 *******************************************************************************
 *
 * @sa coreblas_scabs1
 *
 ******************************************************************************/
float coreblas_scabs1(coreblas_complex32_t alpha)
{
    return fabsf(creal(alpha)) + fabsf(cimag(alpha));
}
