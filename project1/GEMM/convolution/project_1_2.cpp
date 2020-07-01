#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include"head.h"

typedef double dtype;
#define MPI_DTYPE MPI_DOUBLE



int main(int argc, char* argv[]) {
	int rank;
	int proc_num;

	int row_num = 0;
	int col_num = 0;
	int kernel_num = 0;
	int kernel_size = 0;
	int org_row_num = 0;
	int org_col_num = 0;

	dtype* mat_1d;
	dtype** mat_2d;
	dtype* kernel_1d;
	dtype** kernel_2d;
	dtype* vec;
	dtype* local_res;
	dtype* coll_res;  // concatenation of every local_res of each processor
	dtype* final_res = NULL;
	int bias;

	int mode;

	/*int* send_cnt;
	int* send_disp;*/
	/*int* recv_cnt;
	int* recv_disp;*/

	int i, j, k;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);

	if (argc != 4) {  // proc_num, matrix_input_file_name, kernel_input_file_name
		if (!rank) printf("argc should be 3 instead of %d\n", argc);
		MPI_Finalize();
		exit(1);
	}


	MPI_Comm grid_comm;
	int grid_rank;
	int size[2];
	int periodic[2];
	int coord[2];
	int local_col_num, local_row_num;

	size[0] = size[1] = 0;
	periodic[0] = periodic[1] = 0;
	coord[0] = coord[1] = 0;

	double time;
	MPI_Barrier(MPI_COMM_WORLD);
	time = -MPI_Wtime();

	MPI_Dims_create(proc_num, 2, size);
	MPI_Cart_create(MPI_COMM_WORLD, 2, size, periodic, 0, &grid_comm);
	MPI_Comm_rank(grid_comm, &grid_rank);

	if (rank == 0) {
		printf("Process %d: size[0]: %d, size[1]: %d\n", rank, size[0], size[1]);
	}

	read_block_kernel(argv[2], (void***)&kernel_2d, (void**)&kernel_1d, MPI_DTYPE, &kernel_size, &kernel_num, grid_comm);

	int kernel_local_row_num = BLOCK_SIZE(coord[1], size[1], kernel_size * kernel_size);
	int kernel_local_col_num = kernel_num;
	/*printf("Process %d: kernel_size: %d, kernel_num: %d\n",
		rank, kernel_size, kernel_num);
	printf("Process %d: kernel_local_row_num: %d, kernel_local_col_num: %d\n",
		rank, kernel_local_row_num, kernel_local_col_num);*/

	/*for (int i = 0; i < kernel_local_row_num; i++) {
		for (int j = 0; j < kernel_local_col_num; j++) {
			printf("%6.3f\n", kernel_2d[i][j]);
			fflush(stdout);
		}
	}*/

	for (i = 0; i < kernel_local_row_num; i++) {
		printf("Process %d (kernel): ", rank);
		print_vector(kernel_2d[i], MPI_DTYPE, kernel_local_col_num);
		fflush(stdout);
	}


	MPI_Cart_coords(grid_comm, grid_rank, 2, coord);

	printf("Process %d: coord[0]: %d, coord[1]: %d\n", rank, coord[0], coord[1]);


	mode = atoi(argv[3]);
	if (mode == 0) {
		read_kernel_block_matrix(argv[1], (void***)&mat_2d, (void**)&mat_1d, MPI_DTYPE,
			&org_row_num, &org_col_num, &row_num, &col_num, kernel_size, grid_comm);
	}
	else {
		generate_kernel_block_matrix((void***)&mat_2d, (void**)&mat_1d, MPI_DTYPE,
			&org_row_num, &org_col_num, &row_num, &col_num, kernel_size, grid_comm);
	}
	

	printf("Process %d: row_num: %d, col_num: %d\n", rank, row_num, col_num);

	local_row_num = BLOCK_SIZE(coord[0], size[0], row_num);
	local_col_num = BLOCK_SIZE(coord[1], size[1], col_num);
	printf("Process %d: local_row_num: %d, local_col_num: %d\n", rank, local_row_num, local_col_num);

	/*for (int i = 0; i < local_row_num; i++) {
		printf("Process %d: ", rank);
		print_vector(mat_2d[i], MPI_DTYPE, local_col_num);
		fflush(stdout);
	}

	print_block_matrix((void**)mat_2d, MPI_DTYPE, row_num, col_num, grid_comm);*/

	
	// GEMM
	int data_size = get_mpi_dtype_size(MPI_DTYPE);
	local_res = (dtype*)malloc_s(grid_rank, local_row_num * kernel_local_col_num * data_size);
	printf("Process %d : local_row_num: %d, kernel_local_col_num: %d, local_col_num: %d, kernel_local_row_num: %d",
		rank, local_row_num, kernel_local_col_num, local_col_num, kernel_local_row_num);
	fflush(stdout);
	for (i = 0; i < local_row_num; i++) {
		for (j = 0; j < kernel_local_col_num; j++) {
			local_res[i * kernel_local_col_num + j] = 0;
			for (k = 0; k < local_col_num; k++) {
				local_res[i * kernel_local_col_num + j] += mat_2d[i][k] * kernel_2d[k][j];
			}
		}
	}

	/*for (i = 0; i < local_row_num; i++) {
		printf("Process %d (res): ", rank);
		print_vector((void*)&(local_res[i]), MPI_DTYPE, kernel_local_col_num);
		fflush(stdout);
	}*/
	/*printf("Process %d (res): ", rank);
	print_vector((void*)local_res, MPI_DTYPE, kernel_local_col_num * local_row_num);
	fflush(stdout);*/


	// Reduce and collect
	MPI_Comm row_comm;
	MPI_Comm_split(grid_comm, coord[0], coord[1], &row_comm);
	MPI_Comm col_comm;
	MPI_Comm_split(grid_comm, coord[1], coord[0], &col_comm);

	dtype* reduce_buffer = NULL;
	if (coord[1] == 0) {
		reduce_buffer = (dtype*)malloc_s(grid_rank, local_row_num * kernel_local_col_num * data_size);
	}
	MPI_Reduce(local_res, reduce_buffer, local_row_num* kernel_local_col_num, MPI_DTYPE, MPI_SUM, 0, row_comm);

	/*if (coord[1] == 0) {
		printf("Process %d (reduce_buffer): ", rank);
		print_vector((void*)reduce_buffer, MPI_DTYPE, local_row_num * kernel_local_col_num);
		fflush(stdout);
	}*/

	if (coord[1] == 0) {
		if (grid_rank == 0) {
			final_res = (dtype*)malloc_s(grid_rank, row_num * kernel_local_col_num * data_size);
		}
		int* cnt_array;
		int* disp_array;
		get_arrays_by_diff_block_size(coord[1], size[1], row_num, &cnt_array, &disp_array);
		MPI_Gatherv(reduce_buffer, local_row_num * kernel_local_col_num, MPI_DTYPE, final_res,
			cnt_array, disp_array, MPI_DTYPE, 0, col_comm);
	}

	/*if (grid_rank == 0) {
		printf("Process %d (final_res): ", rank);
		print_vector((void*)final_res, MPI_DTYPE, row_num * kernel_local_col_num);
		fflush(stdout);
	}*/

	MPI_Barrier(MPI_COMM_WORLD);
	time += MPI_Wtime();
	
	// only display first 100 rows (if final_res has more than 100 rows)
	if (grid_rank == 0) {
		int cal_row_num = org_row_num - kernel_size + 1;
		int cal_col_num = org_col_num - kernel_size + 1;
		printf("cal_row_num: %d, cal_col_num: %d\n", cal_row_num, cal_col_num);
		printf("Final result:\n");
		int disp_row_num = MIN(100, cal_row_num);
		for (k = 0; k < kernel_local_col_num; k++) {
			printf("#%d channel\n", k);
			fflush(stdout);
			for (i = 0; i < disp_row_num; i++) {  // i < cal_row_num
				printf("row %d: ", i);
				for (j = 0; j < cal_col_num; j++) {
					printf("%6.3f ", final_res[(i * cal_col_num + j) * kernel_local_col_num + k]);
					fflush(stdout);
				}
				printf("\n");
				fflush(stdout);
			}
		}
	}



	printf("Process %d is done.\n", rank);
	fflush(stdout);

	if (!rank) {
		printf("Time: %f\n", time);
		fflush(stdout);
	}

	MPI_Finalize();
	return 0;

}
