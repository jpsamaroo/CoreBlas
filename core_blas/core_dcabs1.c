/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @precisions normal z -> c
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
double coreblas_dcabs1(coreblas_complex64_t alpha)
{
    return fabs(creal(alpha)) + fabs(cimag(alpha));
}
