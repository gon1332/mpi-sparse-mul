#include "test_suite.h"
#include "coo_sparse.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include "mpi.h"
#include "mmio.h"


#define ACCURACY 14
#define RRANGE 10

/* Returns true if they are equal.
 *
 */
static bool vec_compare(const double *const v1,
                        const double *const v2,
                        size_t size)
{
    size_t count = 0;
    for (size_t i = 0; i < size; i++) {
        count += fabs(v1[i] - v2[i]) <= ACCURACY;
    }
    printf("\ncount: %lu\n", count);
    return count == size;
}

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
    int taskid, rc, source;

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
    int tag6 = 7;		// y result from each process


	double *x;			// array to be multiplied
	//double x[8]     = { 1., 3.,  6.,  2., 1., 0.,  5.,  3.};
	/************************* Master task only *******************************/
    if (taskid == MASTER) {
    	// ======================  TEST CASES  ================================
        
        // ===================== .mtx format read =============================
		// input .mtx file name
        char filename[] = "./arrays/matrix_test.mtx";//"fidapm11.mtx" ;
        
		/* Information to be returned
			I: row_ind 
			J: col_ind
		 	val: array values
		 	nz: number of non-zero elements.
		*/
		int nnz;
		int *I, *J;
		double *vall;
		int size;

        nnz = readMtx(filename, &I, &J, &vall, &size);
        printf("Number of Non-zero elements %d\n",nnz);
        
        for (int i=0; i < nnz; i++){
    		fprintf(stdout, "P0 input : %d %d %20.19g\n", I[i] + 1, J[i] + 1, vall[i]);
    	}
    	fflush(stdout);

		/*	COO Format characteristics	*/ 
		printf ("The size of the matrix is: %d\n",size);
		
		int *row_ind, *col_ind; 
		double *val;
		row_ind = I;
		col_ind = J;
		val = vall;
		
		for (int i=0; i < nnz; i++){
			row_ind[i] = row_ind[i] + 1;
			col_ind[i] = col_ind[i] + 1;
    	}
		
		/*x = malloc(nnz*sizeof(double));
		
		// initialize array x
		srand(time(NULL));
		for (int i=0; i < nnz; i++){
			x[i] = rand() % RRANGE + 1.0;
		} */
		/*	COO Format characteristics	*/ 
		
		double x_[8]     = { 1., 3.,  6.,  2., 1., 0.,  5.,  3.};
		
		x = x_;
		
    	// ======================  TEST CASES  ================================

		/* N (size of array) must be divisible by number of tasks */
        if (size % numtasks != 0) {
            printf("[info]: Quitting. Number of MPI tasks must be divisible by size.\n");
            MPI_Abort(MPI_COMM_WORLD, rc);
            exit(EXIT_FAILURE);
        }

        /* number of rows allocated to each process */
        int num_rows_pp = size / numtasks;
        printf("[info]: Number of rows per process : %d\n", num_rows_pp);
        printf("[info]: p0 is MASTER process.\n");

		/* -----------   START OF PROCESS DATA INITIALIZATION   ----------- */
    	
    	int dest = 0;
        int from = 0;
        int divisor = num_rows_pp;
        int product = 0;
        size_t temp_chunksize = 0;

        /* Send each task its portion of the array - master keeps 1st part */
        for (size_t i = 0; i < nnz; i++) {

        	product = row_ind[i] / divisor;

        	/* if row_ind[i] mod divisor == 0 then the array recognized so far
        	   should be passed to the designated destination process and
        	   re-initialize the counters for the next array chunk
            */
            if ((row_ind[i]) % (divisor) == 0 || product > 1) {

                if (dest == MASTER) {
                    temp_chunksize = i;
                }

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
                        print_array_int(row_ind+from, chunksize);

		            	MPI_Send(&col_ind[from], chunksize, MPI_INT, dest, tag4, MPI_COMM_WORLD);
                        printf("[comm@%s](p0 --> p%d) col_ind : ", __TIME__, dest);
                        print_array_int(col_ind+from, chunksize);
		            }

		            /* distribute vector x ---->> EDIT distribute whole vector*/
                    MPI_Send(&x[0], size, MPI_DOUBLE,
                             dest, tag5, MPI_COMM_WORLD);
                    printf("[comm@%s](p0 --> p%d) x : ", __TIME__, dest);
                    print_array_double(&x[0], size);

				}
                
                // ====== next iteration ======
                from = i;

                /* Fix divisor, destination if a long jump is performed in the
                   array. It is used for arrays that do not have at least one nz
                   elem in their rows, therefore they might spoil the execution.
                */
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
            size_t chunksize = nnz - from;

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
                print_array_int(row_ind+from, chunksize);

                MPI_Send(&col_ind[from], chunksize, MPI_INT, dest, tag4, MPI_COMM_WORLD);
                printf("[comm@%s](p0 --> p%d) col_ind : ", __TIME__, dest);
                print_array_int(col_ind+from, chunksize);
            }

            /* distribute vector x ---->> EDITED */
            MPI_Send(&x[0], size, MPI_DOUBLE, dest, tag5, MPI_COMM_WORLD);
            printf("[comm@%s](p0 --> p%d) x : ", __TIME__, dest);
            print_array_double(&x[0], size);
        }

        size_t chunksize = temp_chunksize;
        /* ------------   END OF PROCESS DATA INITIALIZATION   ------------ */


        /* -------   START OF REORDERING THE DATA TO LOCAL/GLOBAL   ------- */

        // Generally the process p will have locally the
        // [p*num_rows_pp, p*num_rows_pp+num_rows_pp) columns
        printf("[p%d]: My private area is the block [(%d,%d), (%d,%d)]\n",
               taskid, taskid*num_rows_pp, taskid*num_rows_pp,
               taskid*num_rows_pp+num_rows_pp, taskid*num_rows_pp+num_rows_pp);

        int private_from = taskid * num_rows_pp;
        int private_to   = taskid * num_rows_pp + num_rows_pp;

        printf("[p%d]: private{ ", taskid);
        for (int i = 0; i < num_rows_pp; i++) {
            printf("%g(%d, %d), ", val[i], row_ind[i], col_ind[i]);
        }
        puts("}");

        /* --------   END OF REORDERING THE DATA TO LOCAL/GLOBAL   -------- */

        double y[size];
        memset(y, 0, size * sizeof(double));

        for (size_t i = 0; i < chunksize; i++) {
            y[row_ind[i]] += val[i] * x[col_ind[i]];
            printf(">> y[%d] = %g * %g\n", row_ind[i], val[i], x[col_ind[i]]);
        }
        printf("[p0]: y: ");
        print_array_double(y, num_rows_pp);


        /* Get final sum and print sample results */
        for (int i = 1; i < numtasks; i++) {
            MPI_Recv(&y[i*num_rows_pp], num_rows_pp, MPI_DOUBLE, i, tag6, MPI_COMM_WORLD, &status);
        }

        printf("[p0]: I'm master and I have the final result which lies in y: ");
        print_array_double(y, size);

		/*
        if (vec_compare(y, y_sol, size)) {
            puts("[TEST 1] *** SUCCESS ***");
        } else {
            puts("[TEST 1] !!! FAILURE !!!");
        }
        */
    }  /* end of master section */

    /********************** Non-master tasks only *****************************/

    if (taskid > MASTER) {

        int num_rows_pp;
        int chunksize;

        /* Receive my portion of array from the master task */
        source = MASTER;
        MPI_Recv(&num_rows_pp, 1, MPI_INT, source, tag0, MPI_COMM_WORLD, &status);
        //printf("[comm@%s](p%d <-- p%d) num_rows_pp: %d\n", __TIME__, taskid, source, num_rows_pp);

        MPI_Recv(&chunksize, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
        //printf("[comm@%s](p%d <-- p%d) chunksize: %d\n", __TIME__, taskid, source, chunksize);

        // Allocate memory for the arrays
        double my_val[chunksize];
        int my_row_ind[chunksize];
        int my_col_ind[chunksize];
        //double my_x[num_rows_pp];

        if (chunksize > 0){
            // Receive values array
        	MPI_Recv(&my_val, chunksize, MPI_DOUBLE, source, tag2, MPI_COMM_WORLD, &status);
			printf("[comm@%s](p%d <-- p%d) val : ", __TIME__, taskid, source);
			print_array_double(my_val, chunksize);

            // Receive row indeces array
            MPI_Recv(&my_row_ind, chunksize, MPI_INT, source, tag3, MPI_COMM_WORLD, &status);
			printf("[comm@%s](p%d <-- p%d) row_ind : ", __TIME__, taskid, source);
			print_array_int(my_row_ind, chunksize);

            // Receive column indeces array
			MPI_Recv(&my_col_ind, chunksize, MPI_INT, source, tag4, MPI_COMM_WORLD, &status);
			printf("[comm@%s](p%d <-- p%d) col_ind : ", __TIME__, taskid, source);
			print_array_int(my_col_ind, chunksize);

		}

        /* receive portion of array x */
        x = malloc(num_rows_pp*numtasks*sizeof(double));
		MPI_Recv(x, num_rows_pp*numtasks, MPI_DOUBLE, source, tag5, MPI_COMM_WORLD, &status);
		printf("[comm@%s](p%d <-- p%d) x : ", __TIME__, taskid, source);
		print_array_double(x, num_rows_pp*numtasks);
        //
        /* ------------   END OF PROCESS DATA INITIALIZATION   ------------ */


        /* -------   START OF REORDERING THE DATA TO LOCAL/GLOBAL   ------- */

        // Generally the process p will have locally the
        // [p*num_rows_pp, p*num_rows_pp+num_rows_pp) columns
        printf("[p%d]: My private area is the block [(%d,%d), (%d,%d)]\n",
               taskid, taskid*num_rows_pp, taskid*num_rows_pp,
               taskid*num_rows_pp+num_rows_pp, taskid*num_rows_pp+num_rows_pp);

        int private_from = taskid * num_rows_pp;
        int private_to   = taskid * num_rows_pp + num_rows_pp;

        printf("[p%d]: private{ ", taskid);
        for (int i = private_from-private_from; i < private_to-private_from; i++) {
            printf("%g(%d, %d), ", my_val[i], my_row_ind[i], my_col_ind[i]);
        }
        puts("}");



        /* --------   END OF REORDERING THE DATA TO LOCAL/GLOBAL   -------- */

        double my_y[num_rows_pp];
        memset(my_y, 0, chunksize * sizeof(double));

        for (size_t i = 0; i < (size_t)chunksize; i++) {
            my_y[my_row_ind[i] - num_rows_pp*taskid] += my_val[i] * x[my_col_ind[i]];
        }
        printf("[p%d]: y: ", taskid);
        print_array_double(my_y, num_rows_pp);

        MPI_Send(my_y, num_rows_pp, MPI_DOUBLE, MASTER, tag6, MPI_COMM_WORLD);

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

/* 	void print_array_*(..)
		Prints the arrays
 */
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
        printf("%d, ", array[i]);
    }
    printf("]\n");
}

/*	int readMtx( FILE *f )
		Read .mtx files 
		Returns the number of nonzero elems (directly)
		Returns the row_ind, col_ind, val 	(indirecly)
 */
int readMtx( char filename[], int **Ii, int **Ji, double **vali, int *size )
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    
    int M, N;
    int nz = 0;  
    int i;
    int *I, *J;
    double *val;
    
	// Read File
         
	if ((f = fopen(filename, "r")) == NULL) 
		exit(1);

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }
    
    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }
    
    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);

    /* reseve memory for matrices */

    I =  malloc(nz * sizeof(int));
    J =  malloc(nz * sizeof(int));
    val =  malloc(nz * sizeof(double));

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (i=0; i < nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }

    if (f !=stdin) 
    	fclose(f);

    /************************/
    /* now write out matrix */
    /************************/

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    //for (i=0; i< nz; i++){
    //    fprintf(stdout, "Lala : %d %d %20.19g\n", I[i]+1, J[i]+1, val[i]);
    //    fflush(stdout);
    //}
    
    *Ii = I;
    *Ji = J;
    *vali = val; 
	*size = N;
	 
	return nz;
}



