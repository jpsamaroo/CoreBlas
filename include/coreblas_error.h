/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 **/
#ifndef COREBLAS_ERROR_H
#define COREBLAS_ERROR_H

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************/
#define coreblas_warning(msg) \
    coreblas_warning_func_line_file(__func__, __LINE__, __FILE__, msg)

#define coreblas_error(msg) \
    coreblas_error_func_line_file(__func__, __LINE__, __FILE__, msg)

#define coreblas_error_with_code(msg, code) \
    coreblas_error_func_line_file_code(__func__, __LINE__, __FILE__, msg, \
                                         code)

#define coreblas_fatal_error(msg) \
    coreblas_fatal_error_func_line_file(__func__, __LINE__, __FILE__, msg)

/******************************************************************************/
static inline void coreblas_warning_func_line_file(
    char const *func, int line, const char *file, const char *msg)
{
    fprintf(stderr,
            "COREBLAS WARNING at %d of %s() in %s: %s\n",
            line, func, file, msg);
}

/******************************************************************************/
static inline void coreblas_error_func_line_file(
    char const *func, int line, const char *file, const char *msg)
{
    fprintf(stderr,
            "COREBLAS ERROR at %d of %s() in %s: %s\n",
            line, func, file, msg);
}

/******************************************************************************/
static inline void coreblas_error_func_line_file_code(
    char const *func, int line, const char *file, const char *msg, int code)
{
    fprintf(stderr,
            "COREBLAS ERROR at %d of %s() in %s: %s %d\n",
            line, func, file, msg, code);
}

/******************************************************************************/
static inline void coreblas_fatal_error_func_line_file(
    char const *func, int line, const char *file, const char *msg)
{
    fprintf(stderr,
            "COREBLAS FATAL ERROR at %d of %s() in %s: %s\n",
            line, func, file, msg);
    exit(EXIT_FAILURE);
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // COREBLAS_ERROR_H
