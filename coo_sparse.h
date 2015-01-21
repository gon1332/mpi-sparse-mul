#ifndef COO_SPARSE_H_7EUNGWRE
#define COO_SPARSE_H_7EUNGWRE
#include <stdlib.h> // for size_t

#define T coo_matrix_T
typedef struct T T;
T *ptr;

#define coo_matrix_original_size(A)         (A->orig_size)
#define coo_matrix_non_zero_elemets(A)      (A->nz)
#define coo_matrix_get_value(A, i)          (A->val[i])
#define coo_matrix_get_column_index(A, i)   (A->col_ind[i])
#define coo_matrix_get_row_index(A, i)      (A->row_ind[i])


/*! \brief Creates a matrix in a triplet form according to COO scheme.
 *
 * @param nz Number of non-zero elements
 * @return Pointer to the triplet form
 * ========================================================================= */
extern T *coo_matrix_new(size_t orig_size, size_t nz);


/*! \brief Initialize the value array of the triplet.
 *
 * @param A Sparse matrix in triplet COO scheme
 * @param val Pointer to the value array
 * ========================================================================= */
extern void coo_matrix_init_values(T *A, double *val);


/*! \brief Initialize the column index array of the triplet.
 *
 * @param A Sparse matrix in triplet COO scheme
 * @param val Pointer to the column index array
 * ========================================================================= */
extern void coo_matrix_init_columns(T *A, size_t *col);


/*! \brief Initialize the row index array of the triplet.
 *
 * @param A Sparse matrix in triplet COO scheme
 * @param val Pointer to the row index array
 * ========================================================================= */
extern void coo_matrix_init_rows(T *A, size_t *row);


/*! \brief Delete the given sparse matrix which uses the COO scheme.
 *
 * @param A Sparse matrix in triplet COO scheme
 * ========================================================================= */
extern void coo_matrix_delete(T *A);


/*! \brief Multiplies matrix A by vector x and stores the result in vector y.
 *
 * @param A Input matrix
 * @param x Input vactor
 * @param y Result of A*x
 * ========================================================================= */
extern void coo_matrix_vector_mul(const T *A,
                                  const double *x,
                                  double *y);

#undef T
#endif /* COO_SPARSE_H_7EUNGWRE */

