#include "test_suite.h"
#include "coo_sparse.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void test_suite(void)
{
    /* TEST 1 */
    {
        size_t size = 3;
        size_t nz = 3;
        coo_matrix_T *test = coo_matrix_new(size, nz);

        double val[3]     = {1., 2., 3.};
        size_t col_ind[3] = {0,  2,  2};
        size_t row_ind[3] = {0,  0,  2};
        coo_matrix_init_values(test, val);
        coo_matrix_init_columns(test, col_ind);
        coo_matrix_init_rows(test, row_ind);

        double x[3]     = {10., 20., 30.};
        double y_sol[3] = {70., 0.,  90.};

        double y_res[3] = {0.};

        coo_matrix_vector_mul(test, x, y_res);

        if (!memcmp(y_res, y_sol, size * sizeof(double))) {
            puts("[TEST 1] *** SUCCESS ***");
        } else {
            puts("[TEST 1] !!! FAILURE !!!");
        }

        coo_matrix_delete(test);
    }

    /* TEST 2 */
    {
        size_t size = 8;
        size_t nz = 14;
        coo_matrix_T *test = coo_matrix_new(size, nz);

        double val[14]     = {6., 9., 4., 4., 5., 3., 5., 8., 6., 5., 4., 3., 2., 2.};
        size_t col_ind[14] = {0,  2,  5,  5,  1,  2,  3,  4,  4,  5,  5,  6,  6,  7};
        size_t row_ind[14] = {0,  0,  0,  1,  2,  3,  3,  3,  4,  5,  6,  6,  7,  7};
        coo_matrix_init_values(test, val);
        coo_matrix_init_columns(test, col_ind);
        coo_matrix_init_rows(test, row_ind);

        double x[8]     = { 1., 3.,  6.,  2., 1., 0.,  5.,  3.};
        double y_sol[8] = {60., 0., 15., 36., 6., 0., 15., 16.};

        double y_res[8] = {0.};

        coo_matrix_vector_mul(test, x, y_res);

        if (!memcmp(y_res, y_sol, size * sizeof(double))) {
            puts("[TEST 2] *** SUCCESS ***");
        } else {
            puts("[TEST 2] !!! FAILURE !!!");
        }

        coo_matrix_delete(test);
    }
}

