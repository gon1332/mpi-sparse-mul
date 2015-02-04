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



void print_array_int(int *array, size_t size);
void print_array_double(double *array, size_t size);


#define  MASTER		0
void parallel_tests(int numtasks, int argc, char *argv[])
{
    int taskid, rc, source, chunksize;

    MPI_Status status;

    /***** Initializations *****/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    printf("[info]: MPI task %d has started...\n", taskid);

	/* Tags are used to keep consistent and healthy flow of data during comm */
    int tag0 = 6;       // num_rows_pp
    int tag1 = 2;		// chuncksize
	int tag2 = 1;		// val_array
    int tag3 = 3;		// row_ind
    int tag4 = 4;		// col_ind
    int tag5 = 5;		// x multiplicator array

	/************************* Master task only *******************************/

    if (taskid == MASTER) {
    	// ======================  TEST CASES  ================================
        // ===================== 1 =========================================
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

        printf("[info]: Values are: ");
        print_array_double(val, nnz);
        printf("[info]: Vector is: ");
        print_array_double(x, size);

		// ======================== 2 =======================================
		/*size_t size = 8;
        size_t nnz = 10;
        coo_matrix_T test = coo_matrix_new(size, nnz);

        double val[14]     = {6., 9., 4., 4., 6., 5., 4., 3., 2., 2.};
        size_t col_ind[14] = {0,  2,  5,  5,  4,  5,  5,  6,  6,  7};
        size_t row_ind[14] = {0,  0,  0,  1,  4,  5,  6,  6,  7,  7};
        coo_matrix_init_values(test, val);
        coo_matrix_init_columns(test, col_ind);
        coo_matrix_init_rows(test, row_ind);

        double x[8]     = { 1., 3.,  6.,  2., 1., 0.,  5.,  3.};
        double y_sol[8] = {60., 0., 15., 36., 6., 0., 15., 16.};

        double y_res[8] = {0.};
		*/
		// ======================== 3 =======================================
		/*size_t size = 8;
        size_t nnz = 6;
        coo_matrix_T test = coo_matrix_new(size, nnz);

        double val[14]     = { 6., 5., 4., 3., 2., 2.};
        size_t col_ind[14] = {  4,  5,  5,  6,  6,  7};
        size_t row_ind[14] = {  4,  5,  6,  6,  7,  7};
        coo_matrix_init_values(test, val);
        coo_matrix_init_columns(test, col_ind);
        coo_matrix_init_rows(test, row_ind);

        double x[8]     = { 1., 3.,  6.,  2., 1., 0.,  5.,  3.};
        double y_sol[8] = {60., 0., 15., 36., 6., 0., 15., 16.};

        double y_res[8] = {0.};
		*/
    	// ======================  TEST CASES  ================================

        if (size % numtasks != 0) {
            printf("[info]: Quitting. Number of MPI tasks must be divisible by size.\n");
            MPI_Abort(MPI_COMM_WORLD, rc);
            exit(EXIT_FAILURE);
        }

        /* number of rows allocated to each process */
        int num_rows_pp = size / numtasks;
        printf("[info]: Number of rows per process : %d\n", num_rows_pp);
        printf("[info]: p0 is MASTER process.\n");


    	int dest = 0;
        int from = 0;
        int divisor = num_rows_pp;
        int product = 0;

        /* -----------   START OF PROCESS DATA INITIALIZATION   ----------- */

        /* Send each task its portion of the array - master keeps 1st part */
        for (size_t i = 0; i < nnz; i++) {

        	product = row_ind[i] / divisor;

        	/* 	if row_ind[i] mod divisor == 0 then the array recognized so far 
        		should be passed to the designated destination process and
        		re-initialize the counters for the next array chunk  */
            if ((row_ind[i]) % (divisor) == 0 || product > 1) {

            	// the rows with value 0 --> we want to keep them in MASTER array
            	if (row_ind[i] == 0)
            		continue;				// ayto mporei kai na fygei

	        	if (dest != MASTER) {
		            size_t chunksize = i - from;

                    /* send META info like chunksize and number of rows per process */
		            MPI_Send(&num_rows_pp, 1, MPI_INT, dest, tag0, MPI_COMM_WORLD);
                    printf("[comm@%s](p0 --> p%d) num_rows_pp : %d\n", __TIME__, dest, num_rows_pp);

		            MPI_Send(&chunksize, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
                    printf("[comm@%s](p0 --> p%d) chunksize : %zd\n", __TIME__, dest, chunksize);

		            if (chunksize > 0) {
		            	MPI_Send(&val[from], chunksize, MPI_DOUBLE, dest, tag2, MPI_COMM_WORLD);
                        printf("[comm@%s](p0 --> p%d) val : ", __TIME__, dest);
                        print_array_double(&val[from], chunksize);

		            	MPI_Send(&row_ind[from], chunksize, MPI_INT, dest, tag3, MPI_COMM_WORLD);
                        printf("[comm@%s](p0 --> p%d) row_ind : ", __TIME__, dest);
                        print_array_int(&row_ind[from], chunksize);

		            	MPI_Send(&col_ind[from], chunksize, MPI_INT, dest, tag4, MPI_COMM_WORLD);
                        printf("[comm@%s](p0 --> p%d) col_ind : ", __TIME__, dest);
                        print_array_int(&row_ind[from], chunksize);
		            }

		            /* distribute vector x */
                    MPI_Send(&x[divisor - num_rows_pp], num_rows_pp, MPI_INT,
                             dest, tag5, MPI_COMM_WORLD);

				}

                // ====== next iteration ======
                from = i;

                /* 	Fix divisor, destination if a long jump is performed in the array
        		It is used for arrays that do not have at least one nz elem
        		in their rows, therefore they might spoil the execution.*/
                if (product >= 1) {
                	// if there is a logic jump in row_ind array re-do the last iteration
                	i--;
        		}
                dest++;
                divisor += num_rows_pp;
            }
        }

        /* send to last process */
        if (dest < numtasks) {
            chunksize = nnz - from;

            MPI_Send(&num_rows_pp, 1, MPI_INT, dest, tag0, MPI_COMM_WORLD);
            printf("[comm@%s](p0 --> p%d) num_rows_pp : %d\n", __TIME__, dest, num_rows_pp);
            MPI_Send(&chunksize, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
            printf("[comm@%s](p0 --> p%d) chunksize : %d\n", __TIME__, dest, chunksize);

            if (chunksize > 0) {
                MPI_Send(&val[from], chunksize, MPI_DOUBLE, dest, tag2, MPI_COMM_WORLD);
                printf("[comm@%s](p0 --> p%d) val : ", __TIME__, dest);
                print_array_double(&val[from], chunksize);

                MPI_Send(&row_ind[from], chunksize, MPI_INT, dest, tag3, MPI_COMM_WORLD);
                printf("[comm@%s](p0 --> p%d) row_ind : ", __TIME__, dest);
                print_array_int(&row_ind[from], chunksize);

                MPI_Send(&col_ind[from], chunksize, MPI_INT, dest, tag4, MPI_COMM_WORLD);
                printf("[comm@%s](p0 --> p%d) col_ind : ", __TIME__, dest);
                print_array_int(&col_ind[from], chunksize);
            }

            /* distribute vector x */
            MPI_Send(&x[divisor - num_rows_pp], num_rows_pp, MPI_DOUBLE, dest, tag5, MPI_COMM_WORLD);
        }
        /* ------------   END OF PROCESS DATA INITIALIZATION   ------------ */


        /* -------   START OF REORDERING THE DATA TO LOCAL/GLOBAL   ------- */

        // Generally the process p will have locally the
        // [p*num_rows_pp, p*num_rows_pp+num_rows_pp) columns
        printf("[p%d]: My private area is the block [(%d,%d), (%d,%d)]\n",
               taskid, taskid*num_rows_pp, taskid*num_rows_pp,
               taskid*num_rows_pp+num_rows_pp, taskid*num_rows_pp+num_rows_pp);


        /* --------   END OF REORDERING THE DATA TO LOCAL/GLOBAL   -------- */



        /* Get final sum and print sample results */

    }  /* end of master section */

    /********************** Non-master tasks only *****************************/

    if (taskid > MASTER) {

        double *my_val;
        double *my_x;
        int *my_row_ind;
        int *my_col_ind;
        int num_rows_pp;

        /* Receive my portion of array from the master task */
        source = MASTER;
        MPI_Recv(&num_rows_pp, 1, MPI_INT, source, tag0, MPI_COMM_WORLD, &status);
        printf("[comm@%s](p%d <-- p%d) num_rows_pp: %d\n", __TIME__, taskid, source, num_rows_pp);

        MPI_Recv(&chunksize, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
        printf("[comm@%s](p%d <-- p%d) chunksize: %d\n", __TIME__, taskid, source, chunksize);

        if (chunksize > 0){
            // Receive values array
        	MPI_Recv(&my_val, chunksize, MPI_DOUBLE, source, tag2, MPI_COMM_WORLD, &status);
			printf("[comm@%s](p%d <-- p%d) val : ", __TIME__, taskid, source);
			print_array_double(&my_val, chunksize);

            // Receive row indeces array
            MPI_Recv(&my_row_ind, chunksize, MPI_INT, source, tag3, MPI_COMM_WORLD, &status);
			printf("[comm@%s](p%d <-- p%d) row_ind : ", __TIME__, taskid, source);
			print_array_int(&my_row_ind, chunksize);

            // Receive column indeces array
			MPI_Recv(&my_col_ind, chunksize, MPI_INT, source, tag4, MPI_COMM_WORLD, &status);
			printf("[comm@%s](p%d <-- p%d) col_ind : ", __TIME__, taskid, source);
			print_array_int(&my_col_ind, chunksize);

		}

        /* receive portion of array x */
		MPI_Recv(&my_x, num_rows_pp, MPI_DOUBLE, source, tag5, MPI_COMM_WORLD, &status);
		printf("[comm@%s](p%d <-- p%d) x : ", __TIME__, taskid, source);
		print_array_double(&my_x, num_rows_pp);


        /* -------   START OF REORDERING THE DATA TO LOCAL/GLOBAL   ------- */

        // Generally the process p will have locally the
        // [p*num_rows_pp, p*num_rows_pp+num_rows_pp) columns
        printf("[p%d]: My private area is the block [(%d,%d), (%d,%d)]\n",
               taskid, taskid*num_rows_pp, taskid*num_rows_pp,
               taskid*num_rows_pp+num_rows_pp, taskid*num_rows_pp+num_rows_pp);

        int private_from = taskid * num_rows_pp;
        int private_to   = taskid * num_rows_pp + num_rows_pp;

        /* --------   END OF REORDERING THE DATA TO LOCAL/GLOBAL   -------- */

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

void print_array_double(double *array, size_t size)
{
    printf("[ ");
    for (size_t i = 0; i < size; i++) {
        printf("%g, ", array[i]);
    }
    printf("]\n");
}


void print_array_int(int *array, size_t size)
{
    printf("[ ");
    for (size_t i = 0; i < size; i++) {
        printf("%i, ", array[i]);
    }
    printf("]\n");
}
