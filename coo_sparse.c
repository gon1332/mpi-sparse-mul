#include "coo_sparse.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

#define T coo_matrix_T

/*!
 * A structure to represent sparse matrices in COO format.
 *
 */
struct T {
    double *val;    /*!< array of non-zero values */
    size_t *col_ind;/*!< array of row indices */
    size_t *row_ind;/*!< array of column indices */

    /*!
     * \name Supporting fields
     */
    size_t nz;          /*!< number of non-zero elements*/
    size_t orig_size;   /*!< dimension of the original sparse matrix*/
};


/*!
 * \details
 * Creates a new matrix by allocating memory and
 * initializing values of the structure fields.
 *
 * COO format was chosen for the needs of a project
 * for the lesson *Parallel and Distributed Computation*.
 *
 * \todo Embed the custom malloc checking mechanism
 * `Safecall`.
 *
 */
T coo_matrix_new(size_t orig_size, size_t nz)
{
    T ret = malloc(sizeof(struct T));

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


/*!
 * \details
 * Initializes the values array of the coo matrix by copying the
 * array `val` to the allocated memory from `coo_matrix_new`.
 *
 * \pre First create a coo matrix with `coo_matrix_new`.
 */
void coo_matrix_init_values(T A, double *val)
{
    assert(A);
    assert(val);

    memcpy(A->val, val, A->nz * sizeof(double));
}


/*!
 * \details
 * Initializes the columns array of the coo matrix by copying the
 * array `val` to the allocated memory from `coo_matrix_new`.
 *
 * \pre First create a coo matrix with `coo_matrix_new`.
 */
void coo_matrix_init_columns(T A, size_t *col)
{
    assert(A);
    assert(col);

    memcpy(A->col_ind, col, A->nz * sizeof(size_t));
}


/*!
 * \details
 * Initializes the rows array of the coo matrix by copying the
 * array `val` to the allocated memory from `coo_matrix_new`.
 *
 * \pre First create a coo matrix with `coo_matrix_new`.
 */
void coo_matrix_init_rows(T A, size_t *row)
{
    assert(A);
    assert(row);

    memcpy(A->row_ind, row, A->nz * sizeof(size_t));
}


/*!
 * \details
 * Frees the memory allocated by `coo_matrix_new` for the `A` coo
 * matrix. Also sets the corresponding pointers to `NULL`.
 *
 * \pre First create a coo matrix with `coo_matrix_new`.
 */
void coo_matrix_delete(T A)
{
    assert(A);

    free(A->val);     A->val     = NULL;
    free(A->col_ind); A->col_ind = NULL;
    free(A->row_ind); A->row_ind = NULL;
    free(A);          A          = NULL;
}


/*!
 * \details
 * Performs a matrix-vector multiplication that fits the
 * COO format of sparse matrices. If the coo matrix has
 * `nnz` elements and the names below for the arrays that
 * represent it:
 *
 * * `values` - *array of non-zero values from original matrix*
 * * `row_ind` - *array of row indices corresponding to the rows
 * of the original matrix*
 * * `col_ind` - *array of row indices corresponding to the columns
 * of the original matrix*
 *
 *
 * the algorithm used is:
 *
 * ~~~~~~~~~~~~~~~~~
 *     for i in 1..nnz:
 *         y[row_ind[i]] = y[row_ind[i]] + values[i] * x[col_ind[i]]
 * ~~~~~~~~~~~~~~~~~
 *
 *
 * \pre First create and initialize a coo matrix with `coo_matrix_new`.
 */
void coo_matrix_vector_mul(const T A,
                           const double *x,
                           double *y)
{
    assert(A);

    for (size_t i = 0; i < coo_matrix_get_non_zero_elements(A); i++) {
        register size_t row_ind = coo_matrix_get_row_index(A, i);
        register size_t col_ind = coo_matrix_get_column_index(A, i);
        register double value   = coo_matrix_get_value(A, i);
        y[row_ind] += value * x[col_ind];
    }
}

#undef T
