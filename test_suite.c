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
    int taskid, rc, dest, offset, i, tag1, tag2, tag3, tag4, tag5 , source, chunksize;

    MPI_Status status;

    /***** Initializations *****/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    printf("MPI task %d has started...\n", taskid);

	/* Tags are used to keep consistent and healthy flow of data during comm */
	tag2 = 1;		// val_array
    tag1 = 2;		// chuncksize
    tag3 = 3;		// row_ind
    tag4 = 4;		// col_ind
    tag5 = 5;		// x multiplicator array
    
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
            printf("Quitting. Number of MPI tasks must be divisible by 4.\n");
            MPI_Abort(MPI_COMM_WORLD, rc);
            exit(EXIT_FAILURE);
        }
        /* number of rows allocated to each process */
        int num_rows_pp = size / numtasks;		
        printf("num_rows_pp : %d\n", num_rows_pp);
        
        /* Send each task its portion of the array - master keeps 1st part */
    	size_t dest = 0;
        size_t from = 0;
        int divisor = num_rows_pp;
        int product = 0;
        for (size_t i = 0; i < nnz; i++) {
        
        	product = row_ind[i] / divisor;
        	
        	/* 	if row_ind[i] mod divisor == 0 then the array recognized so far 
        		should be passed to the designated destination process and 
        		re-initialize the counters for the next array chunk  */
            if ((row_ind[i]) % (divisor) == 0 || product > 1 ) {
            	
            	// the rows with value 0 --> we want to keep them in MASTER array
            	if (row_ind[i] == 0)
            		continue;				// ayto mporei kai na fygei
            	
	        	if (dest != MASTER){
		            size_t chunksize = i - from;
		            printf("MASTER has sent to process %zd the array: ", dest);
		            print_array(&val[from], chunksize);
		            printf("\n");

					/* number of elements to be sent */
		            MPI_Send(&chunksize, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
		            
		            if (chunksize > 0) {
		            	MPI_Send(&val[from], chunksize, MPI_DOUBLE, dest, tag2, MPI_COMM_WORLD);
		            }
		            
		            /* distribute array x */
		            
				}
				            
                // ====== next iteration ====== 
                from = i;
                
                /* 	Fix divisor, destination if a long jump is performed in the array 
        		It is used for arrays that do not have at least one nz elem
        		in their rows, therefore they might spoil the execution.*/
                if (product >= 1){
                	// if there is a logic jump in row_ind array re-do the last iteration 
                	i --; 
        		}
                dest++;
                divisor += num_rows_pp;
            }
        }
        /*** send to last process ***/
        chunksize = nnz - from;
        printf("MASTER has sent to process %zd the array: ", dest);
        print_array(&val[from], chunksize);
        MPI_Send(&chunksize, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
        
        if (chunksize > 0) {
        	MPI_Send(&val[from], chunksize, MPI_DOUBLE, dest, tag2, MPI_COMM_WORLD);
        }
        
        //MPI_Send(&val[from], chunksize, MPI_DOUBLE, dest, tag2, MPI_COMM_WORLD);
		
       
        printf("Sent %zd elements to task %d offset= %zd\n", chunksize, dest, offset);

        /* Master does its part of the work */
        //offset = 0;
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

    /********************** Non-master tasks only *****************************/

    if (taskid > MASTER) {

        double *my_ex;
        /* Receive my portion of array from the master task */
        
        source = MASTER;
        MPI_Recv(&chunksize, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
        printf("I process %d have received from process %d the chunksize %d\n",
        		taskid, source, chunksize);
        
        if  (chunksize > 0){
        	MPI_Recv(&my_ex, chunksize, MPI_DOUBLE, source, tag2, MPI_COMM_WORLD, &status);

			printf(" I process %d have received :", taskid);
			print_array(&my_ex, chunksize);
		}
        /* Send my results back to the master task */
        //dest = MASTER;
        //MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
        //MPI_Send(&my_ex, chunksize, MPI_FLOAT, MASTER, tag2, MPI_COMM_WORLD);


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
