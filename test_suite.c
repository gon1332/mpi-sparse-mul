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
    int tagx = 8;		// size of matrix N


	double *x;			// array to be multiplied
	//double x[8]     = { 1., 3.,  6.,  2., 1., 0.,  5.,  3.};
	/************************* Master task only *******************************/
    if (taskid == MASTER) {
    	// ======================  TEST CASES  ================================
        
        // ===================== .mtx format read =============================
		// input .mtx file name
        //char filename[] = "./arrays/matrix_test.mtx";//"fidapm11.mtx" ;
        char filename[] = "./arrays/fidapm11_edited.mtx";
        
        /*
        // "conf6_0_00l8x8_8000.mtx"
        // "fidapm11_edited.mtx"
        */
        
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
        printf("--Number of Non-zero elements %d\n",nnz);
        
        /* Print output matrix (optional) */
        /*for (int i=0; i < nnz; i++){
    		fprintf(stdout, "P0 input : %d %d %20.19g\n", I[i] + 1, J[i] + 1, vall[i]);
    	}
    	fflush(stdout);*/

		/*	COO Format characteristics	*/ 
		printf ("--The size of the matrix is: %d\n",size);
		
		int *row_ind, *col_ind; 
		double *val;
		col_ind = I;
		row_ind = J;
		val = vall;
		
		for (int i=0; i < nnz; i++){
			row_ind[i] = row_ind[i] + 1;
			col_ind[i] = col_ind[i] + 1;
    	}
		
		// initialize array x
		x = malloc(size*sizeof(double));
		
		srand(time(NULL));
		for (int i=0; i < size; i++){
			x[i] = 1.0; //rand() % RRANGE + 1.0;
		}
		/*	COO Format characteristics	*/ 
		
		//double x_[8]     = { 1., 3.,  6.,  2., 1., 0.,  5.,  3.};
		//x = x_;
		
    	// ======================  TEST CASES  ================================

		/* N (size of array) must be divisible by number of tasks */
        if ( size % numtasks != 0 ) {
            printf("[info]: Quitting. Number of MPI tasks must be divisible by size.\n");
            MPI_Abort(MPI_COMM_WORLD, rc);
            exit(EXIT_FAILURE);
        }

		/* -----------   START OF PROCESS DATA INITIALIZATION   ------------- */
    	
		/******************* ex3  ARRAY DISTRIBUTION START ********************/
        
        int limit = nnz/numtasks;	/* NZ/p elements per process   */
        int cnz;					/* current num of non zero read*/
		int current_row = 0;
		int i = 0;
		int *proc = malloc(numtasks*sizeof(int));
	
		int *total_elements = malloc(numtasks*sizeof(int));
		for ( int j=0; j < numtasks; j++ ){
			total_elements[j] = 0;
			proc[j] = -1;
		}
		
		int k = 0;
		
		while (current_row < size-1) {
		
			current_row = row_ind[i];
			cnz = 0;
			while ((i < nnz) && (row_ind[i] == current_row)){
				cnz += 1;
				i += 1;
			}
			//printf("---Line %d has %d non zero elements ,\n",current_row, cnz);		//!!!!!!!!!!!!!!
		
			if (total_elements[k] + cnz <= limit ){
				total_elements[k] += cnz;
				proc[k] = current_row;
			}
			else{
			k += 1;
			if (k > numtasks - 1)
				k -= 1;
			
			total_elements[k] += cnz;
			proc[k] = current_row;
			}
		}
		
		// validate ---> print boundary_list (proc) and total_elements list
		printf("Boundary List (proc). List will be split by lines: ");
		print_array_int(proc, numtasks);
		printf("\n----------------------------------------------\n");
		printf("Total Elements List. Each process will have: ");
		print_array_int(total_elements, numtasks);
		printf(" elements \n -----------------------------------\n"); 

		/************** distribute elements to other processes ****************/
		int prv_line = proc[0];
		int elements_distributed = total_elements[0];
		int elem = 0;
		for (size_t i = 1; i < numtasks; i++){
		
			/* Send size N of matrix to other proc */
			MPI_Send(&size, 1, MPI_INT, i, tagx, MPI_COMM_WORLD);
			//printf("[comm@%s](p0 --> p%d) size N: %d\n", __TIME__, i, size);
			
			/* how many lines am i going to take? */
			elem = proc[i] - prv_line;
			MPI_Send(&elem, 1, MPI_INT, i, tag0, MPI_COMM_WORLD);
            //printf("[comm@%s](p0 --> p%d) num_rows_p : %d\n", __TIME__, i, (proc[i] - prv_line));
			
			/* how many elements am i sending? */
            MPI_Send(&total_elements[i], 1, MPI_INT, i, tag1, MPI_COMM_WORLD);
            //printf("[comm@%s](p0 --> p%d) chunksize : %zd\n", __TIME__, i, total_elements[i]);
			
			/* send the chunks of COO arrays (based on chunksize >> total_elements[i]) */
			if (total_elements[i] > 0) {
            	MPI_Send(&val[elements_distributed], total_elements[i], MPI_DOUBLE, i, tag2, MPI_COMM_WORLD);
                //printf("[comm@%s](p0 --> p%d) val : ", __TIME__, i);
                //print_array_double(&val[elements_distributed], total_elements[i]);

            	MPI_Send(&row_ind[elements_distributed], total_elements[i], MPI_INT, i, tag3, MPI_COMM_WORLD);
                //printf("[comm@%s](p0 --> p%d) row_ind : ", __TIME__, i);
                //print_array_int(&row_ind[elements_distributed], total_elements[i]);

            	MPI_Send(&col_ind[elements_distributed], total_elements[i], MPI_INT, i, tag4, MPI_COMM_WORLD);
                //printf("[comm@%s](p0 --> p%d) col_ind : ", __TIME__, i);
                //print_array_int(&col_ind[elements_distributed], total_elements[i]);
            }
            
            /* distribute vector x ---->> EDIT distribute whole vector*/
            MPI_Send(&x[0], size, MPI_DOUBLE, i, tag5, MPI_COMM_WORLD);
            //printf("[comm@%s](p0 --> p%d) x : ", __TIME__, i);
            //print_array_double(&x[0], size);
            
			/* for next iteration: */
			elements_distributed += total_elements[i];
			prv_line = proc[i];
		}
		/******************* ex3  ARRAY DISTRIBUTION END *********************/

		/* Compute multiplication result y */
		// locally
        double y[size];
        memset(y, 0, size * sizeof(double));

        for (size_t i = 0; i < total_elements[0]; i++) {
            y[row_ind[i]] += val[i] * x[col_ind[i]];
        
        	//printf(">> y[%d] = %g * %g\n", row_ind[i], val[i], x[col_ind[i]]);
            
        }
        //printf("[p0]: y: ");
        //print_array_double(y, proc[0] + 1);

        /* Get final sum and print sample results */
        prv_line = proc[0];
        /* Receiving strategy 
         	Since every process has a different chunksize and handles different
         	types of rows for every process i :
         	we start receiving on the line prv_line = proc[i-1] and we receive 
         	that many elements as the difference of the 1st line that p[i] has
         	substracted by the last line that p[i] has, namely proc[i+1] - proc[i]
         */
        for (int i = 1; i < numtasks; i++) {
        	
            MPI_Recv(&y[prv_line+1], (proc[i] - prv_line), MPI_DOUBLE, i, tag6, MPI_COMM_WORLD, &status);
            prv_line = proc[i];
        }

        printf("[p0]: I'm master and I have the final result which lies in y: ");
        print_array_double_result(y, size);

		/* Validate results */
		/*
        if (vec_compare(y, y_sol, size)) {
            puts("[TEST 1] *** SUCCESS ***");
        } else {
            puts("[TEST 1] !!! FAILURE !!!");
        }
        */
    }	/* end of master section */

    /********************** Non-master tasks only *****************************/

    if (taskid > MASTER) {

        int num_rows_pp;
        int chunksize;
		int size;
	
        /* Receive my portion of array from the master task */
        source = MASTER;
        
        MPI_Recv(&size, 1, MPI_INT, source, tagx, MPI_COMM_WORLD, &status);
        printf("[comm@%s](p%d <-- p%d) size N: %d\n", __TIME__, taskid, source, size);
        
        MPI_Recv(&num_rows_pp, 1, MPI_INT, source, tag0, MPI_COMM_WORLD, &status);
        printf("[comm@%s](p%d <-- p%d) num_rows_p: %d\n", __TIME__, taskid, source, num_rows_pp);

        MPI_Recv(&chunksize, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
        printf("[comm@%s](p%d <-- p%d) chunksize: %d\n", __TIME__, taskid, source, chunksize);

        // Allocate memory for the arrays
        double my_val[chunksize];
        int my_row_ind[chunksize];
        int my_col_ind[chunksize];

        if (chunksize > 0){
            // Receive values array
        	MPI_Recv(&my_val, chunksize, MPI_DOUBLE, source, tag2, MPI_COMM_WORLD, &status);
			//printf("[comm@%s](p%d <-- p%d) val : ", __TIME__, taskid, source);
			//print_array_double(my_val, chunksize);

            // Receive row indeces array
            MPI_Recv(&my_row_ind, chunksize, MPI_INT, source, tag3, MPI_COMM_WORLD, &status);
			//printf("[comm@%s](p%d <-- p%d) row_ind : ", __TIME__, taskid, source);
			//print_array_int(my_row_ind, chunksize);

            // Receive column indeces array
			MPI_Recv(&my_col_ind, chunksize, MPI_INT, source, tag4, MPI_COMM_WORLD, &status);
			//printf("[comm@%s](p%d <-- p%d) col_ind : ", __TIME__, taskid, source);
			//print_array_int(my_col_ind, chunksize);

		}

        /* receive portion of array x */
        x = malloc(size*sizeof(double));
		MPI_Recv(x, size, MPI_DOUBLE, source, tag5, MPI_COMM_WORLD, &status);
		//printf("[comm@%s](p%d <-- p%d) x : ", __TIME__, taskid, source);
		//print_array_double(x, size);
        //
        /* ------------   END OF PROCESS DATA INITIALIZATION   ------------ */
		
        double my_y[num_rows_pp];
        memset(my_y, 0, num_rows_pp * sizeof(double));
        
        for (int i = 0; i < chunksize; i++) {
 			my_y[my_row_ind[i] - my_row_ind[0]] += my_val[i] * x[my_col_ind[i]];    
 			
            //printf(">> y[%d] = %g * %g\n", my_row_ind[i] - my_row_ind[0], my_val[i], x[my_col_ind[i]]);
        }
        
        //printf("[p%d]: y: ", taskid);
        //print_array_double(my_y, num_rows_pp);

        MPI_Send(my_y, num_rows_pp, MPI_DOUBLE, MASTER, tag6, MPI_COMM_WORLD);

    } /* end of non-master */

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

void print_array_double_result(double *array, size_t size)
{
    printf("\n[ ");
    for (int i = 0; i < size; i++) {
        printf(" y[%d] %g,  \n", i ,array[i]);
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



