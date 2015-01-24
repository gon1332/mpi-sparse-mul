#include "test_suite.h"
#include "coo_sparse.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"

void test_suite(void)
{
    /* TEST 1 */
    {
        size_t size = 3;
        size_t nz = 3;
        coo_matrix_T test = coo_matrix_new(size, nz);

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
        coo_matrix_T test = coo_matrix_new(size, nz);

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


void print_array(double *array, size_t size);


#define  MASTER		0
void parallel_tests(int numtasks, int argc, char *argv[])
{
    int taskid, rc, dest, offset, i, tag1, tag2, source, chunksize;

    MPI_Status status;

    /***** Initializations *****/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    printf("MPI task %d has started...\n", taskid);


    // ======================  TEST CASES  ====================================

    /***** Master task only ******/
    if (taskid == MASTER) {
        size_t size = 8;
        size_t nnz = 14;
        coo_matrix_T test = coo_matrix_new(size, nnz);

        double val[14]     = {6., 9., 4., 4., 5., 3., 5., 8., 6., 5., 4., 3., 2., 2.};
        size_t col_ind[14] = {0,  2,  5,  5,  1,  2,  3,  4,  4,  5,  5,  6,  6,  7};
        size_t row_ind[14] = {0,  0,  0,  1,  2,  3,  3,  3,  4,  5,  6,  6,  7,  7};
        coo_matrix_init_values(test, val);
        coo_matrix_init_columns(test, col_ind);
        coo_matrix_init_rows(test, row_ind);

        double x[8]     = { 1., 3.,  6.,  2., 1., 0.,  5.,  3.};
        double y_sol[8] = {60., 0., 15., 36., 6., 0., 15., 16.};

        double y_res[8] = {0.};


    // ======================  TEST CASES  ====================================

        if (size % numtasks != 0) {
            printf("Quitting. Number of MPI tasks must be divisible by 4.\n");
            MPI_Abort(MPI_COMM_WORLD, rc);
            exit(EXIT_FAILURE);
        }
        int num_rows_pp = size / numtasks;
        tag2 = 1;
        tag1 = 2;

        /* Send each task its portion of the array - master keeps 1st part */
        int dest = 0;
        int from = 0;
        for (size_t i = 0; i < nnz; i++) {
            if ((row_ind[i]+1) % (num_rows_pp+1) == 0) {
                if (dest != MASTER) {
                    int chunksize = i - from;
                    printf("MASTER has sent to process %d ", dest);
                    printf("array: ");
                    print_array(&val[from], chunksize);

                    MPI_Send(&chunksize, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
                    MPI_Send(&val[from], chunksize, MPI_DOUBLE, dest, tag2, MPI_COMM_WORLD);
                }
                dest++;
                from = i;
            }
        }


        printf("Sent %d elements to task %d offset= %d\n", chunksize, dest, offset);

        /* Master does its part of the work */
        offset = 0;
        //mysum = update(offset, chunksize, taskid);

        /* Wait to receive results from each task */
        /*
        for (i = 1; i < numtasks; i++) {
            source = i;
            MPI_Recv(&offset, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
            MPI_Recv(&ex[offset], chunksize, MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
            printf("offset recved from %d is %d\n"
                   "vector recved is [%d, %d]\n\n", i, offset, ex[offset], ex[offset+1]);
        }
        */

        /* Get final sum and print sample results */

    }  /* end of master section */



    /***** Non-master tasks only *****/

    if (taskid > MASTER) {

        int my_ex[2];
        /* Receive my portion of array from the master task */
        source = MASTER;
        MPI_Recv(&offset, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
        MPI_Recv(&my_ex, chunksize, MPI_FLOAT, source, tag2,
                    MPI_COMM_WORLD, &status);

        printf("process %d offset is %d: [%d, %d].\n", taskid, offset, my_ex[0], my_ex[1]);
        my_ex[0] = offset;
        my_ex[1] = offset;


        /* Send my results back to the master task */
        dest = MASTER;
        MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
        MPI_Send(&my_ex, chunksize, MPI_FLOAT, MASTER, tag2, MPI_COMM_WORLD);


    } /* end of non-master */

// FOR MASTER
    /*
    coo_matrix_vector_mul(test, x, y_res);

    if (!memcmp(y_res, y_sol, size * sizeof(double))) {
        puts("[TEST 2] *** SUCCESS ***");
    } else {
        puts("[TEST 2] !!! FAILURE !!!");
    }
    */

    // coo_matrix_delete(test);

    MPI_Finalize();
}

void print_array(double *array, size_t size)
{
    printf("[ ");
    for (size_t i = 0; i < size; i++) {
        printf("%lf ", array[i]);
    }
    printf("]\n");
}
