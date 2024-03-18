/**
 *
 * @file
 *
 *  COREBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 **/
#ifndef COREBLAS_WORKSPACE_H
#define COREBLAS_WORKSPACE_H

#include "coreblas_types.h"

#include <stdlib.h>
#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    void **spaces;      ///< array of nthread pointers to workspaces
    size_t lworkspace;  ///< length in elements of workspace on each core
    int nthread;        ///< number of threads
    coreblas_enum_t dtyp; ///< precision of the workspace
} coreblas_workspace_t;

/******************************************************************************/
int coreblas_workspace_create(coreblas_workspace_t *workspace, size_t lworkspace,
                           coreblas_enum_t dtyp);

int coreblas_workspace_destroy(coreblas_workspace_t *workspace);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // COREBLAS_WORKSPACE_H
