#include "coo_sparse.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

#define T coo_matrix_T

/*! \brief Matrix in triplet form using COO scheme.
 *
 * ========================================================================= */
struct T {
    double *val;
    size_t *col_ind;
    size_t *row_ind;

    size_t nz;
    size_t orig_size;
};


T *coo_matrix_new(size_t orig_size, size_t nz)
{
    T *ret = malloc(sizeof(T));

    register size_t bytes = nz * sizeof(size_t);
    ret->val = malloc(nz * sizeof(double));
    ret->col_ind = malloc(bytes);
    ret->row_ind = malloc(bytes);

    ret->nz = nz;
    ret->orig_size = orig_size;

    if (!ret || !ret->val || !ret->col_ind || !ret->row_ind) {
        perror("coo_serial");
        exit(EXIT_FAILURE);
    }

    return ret;
}


void coo_matrix_init_values(T *A, double *val)
{
    assert(A);
    assert(val);

    memcpy(A->val, val, A->nz * sizeof(double));
}


void coo_matrix_init_columns(T *A, size_t *col)
{
    assert(A);
    assert(col);

    memcpy(A->col_ind, col, A->nz * sizeof(size_t));
}


void coo_matrix_init_rows(T *A, size_t *row)
{
    assert(A);
    assert(row);

    memcpy(A->row_ind, row, A->nz * sizeof(size_t));
}


void coo_matrix_delete(T *A)
{
    assert(A);

    free(A->val);     A->val     = NULL;
    free(A->col_ind); A->col_ind = NULL;
    free(A->row_ind); A->row_ind = NULL;
    free(A);          A          = NULL;
}


void coo_matrix_vector_mul(const T *A,
                           const double *x,
                           double *y)
{
    assert(A);

    for (size_t i = 0; i < coo_matrix_non_zero_elemets(A); i++) {
        register size_t row_ind = coo_matrix_get_row_index(A, i);
        register size_t col_ind = coo_matrix_get_column_index(A, i);
        register double value   = coo_matrix_get_value(A, i);
        y[row_ind] += value * x[col_ind];
    }
}

#undef T
