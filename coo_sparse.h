#ifndef COO_SPARSE_H_7EUNGWRE
#define COO_SPARSE_H_7EUNGWRE
#include <stdlib.h> // for size_t

/*!
 * \file coo_sparse.h
 * \brief Interface of COO matrix data type.
 */

#define T coo_matrix_T
typedef struct T *T;

#define coo_matrix_get_original_size(A)     (A->orig_size)
#define coo_matrix_get_non_zero_elements(A) (A->nz)
#define coo_matrix_get_value(A, i)          (A->val[i])
#define coo_matrix_get_column_index(A, i)   (A->col_ind[i])
#define coo_matrix_get_row_index(A, i)      (A->row_ind[i])


/*! \brief Creates a matrix in a triplet form according to COO scheme.
 *
 * \param orig_size The dimension of the original square matrix
 * \param nz Number of non-zero elements
 * \return Pointer to the coo_matrix created
 *
 * Usage:
 * ~~~~~~~~~~~~~~~~{.c}
 *     coo_matrix_T cm;
 *     cm = coo_matrix_new(50, 3);
 * ~~~~~~~~~~~~~~~~
 */
extern T coo_matrix_new(size_t orig_size, size_t nz);


/*! \brief Initializes the value array of the triplet.
 *
 * \param A     Sparse matrix in triplet COO scheme
 * \param val   Pointer to the value array
 *
 * Usage:
 * ~~~~~~~~~~~~~~~~{.c}
 *     // For an allocated matrix cm with 10 non-zero values:
 *
 *     double values[10] = { 1., 34., 21., 65.43, 0., 74.2, 0.3, 8., 19.2 };
 *
 *     coo_matrix_init_values(cm, values);
 * ~~~~~~~~~~~~~~~~
 */
extern void coo_matrix_init_values(T A, double *val);


/*! \brief Initialize the column index array of the triplet.
 *
 * \param A     Sparse matrix in triplet COO scheme
 * \param val   Pointer to the column index array
 *
 * Usage:
 * ~~~~~~~~~~~~~~~~{.c}
 *     // For an allocated matrix cm with 10 non-zero values:
 *
 *     double row_indices[10] = { 0, 0, 1, 15, 28, 34, 35, 48, 49 };
 *
 *     coo_matrix_init_values(cm, row_indices);
 * ~~~~~~~~~~~~~~~~
 */
extern void coo_matrix_init_columns(T A, size_t *col);


/*! \brief Initialize the row index array of the triplet.
 *
 * \param A     Sparse matrix in triplet COO scheme
 * \param val   Pointer to the row index array
 *
 * Usage:
 * ~~~~~~~~~~~~~~~~{.c}
 *     // For an allocated matrix cm with 10 non-zero values:
 *
 *     double col_indices[10] = { 1, 14, 21, 15, 20, 34, 0, 8, 19 };
 *
 *     coo_matrix_init_values(cm, col_indices);
 * ~~~~~~~~~~~~~~~~
 */
extern void coo_matrix_init_rows(T A, size_t *row);


/*! \brief Delete the given sparse matrix which uses the COO scheme.
 *
 * \param A     Sparse matrix in triplet COO scheme
 */
extern void coo_matrix_delete(T A);


/*! \brief Multiplies matrix A by vector x and stores the result in vector y.
 *
 * \param A[in]     Input matrix
 * \param x[in]     Input vactor
 * \param y[out]    Result of \f$ A*x \f$
 */
extern void coo_matrix_vector_mul(const T A,
                                  const double *x,
                                  double *y);

#undef T
#endif /* COO_SPARSE_H_7EUNGWRE */

