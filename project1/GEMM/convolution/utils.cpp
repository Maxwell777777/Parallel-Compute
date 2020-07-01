#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<mpi.h>
#include"head.h"

int get_mpi_dtype_size(MPI_Datatype mpi_dtype) {
	switch (mpi_dtype)
	{
	case MPI_CHAR: return sizeof(char);
	case MPI_INT: return sizeof(int);
	case MPI_FLOAT: return sizeof(float);
	case MPI_DOUBLE: return sizeof(double);
	}
	printf("Invalid MPI data type.\n");
	fflush(stdout);
	return -1;
}

int terminate(int rank, const char* msg) {
	if (!rank) {
		printf("ERROR: %s\n", msg);
		fflush(stdout);
	}
	MPI_Finalize();
	exit(-1);
}

void* malloc_s(int rank, int num_bytes) {
	void* _buffer = malloc(num_bytes);
	if (_buffer == NULL) {
		printf("ERROR: function malloc failed in process %d\n", rank);
		fflush(stdout);
		MPI_Abort(MPI_COMM_WORLD, MALLOC_ERROR);
	}
	return _buffer;
}

void get_arrays_by_diff_block_size(int rank, int proc_num, int elem_num, int** cnt_array, int** disp_array) {
	*cnt_array = (int*)malloc_s(rank, proc_num * sizeof(int) * 2);
	*disp_array = (int*)malloc_s(rank, proc_num * sizeof(int) * 2);
	(*cnt_array)[0] = BLOCK_SIZE(0, proc_num, elem_num);
	(*disp_array)[0] = 0;
	for (int i = 1; i < proc_num; i++) {
		(*cnt_array)[i] = BLOCK_SIZE(i, proc_num, elem_num);
		(*disp_array)[i] = (*disp_array)[i - 1] + (*cnt_array)[i - 1];
	}
}


void get_arrays_by_curr_block_size(int rank, int proc_num, int elem_num, int** cnt_array, int** disp_array) {
	*cnt_array = (int*)malloc_s(rank, proc_num * sizeof(int) * 2);
	*disp_array = (int*)malloc_s(rank, proc_num * sizeof(int) * 2);
	(*cnt_array)[0] = BLOCK_SIZE(rank, proc_num, elem_num);
	(*disp_array)[0] = 0;
	for (int i = 1; i < proc_num; i++) {
		(*cnt_array)[i] = (*cnt_array)[0];
		(*disp_array)[i] = (*disp_array)[i - 1] + (*cnt_array)[i - 1];
	}
}


void read_vector(char* file_name, void** vector, MPI_Datatype mpi_dtype, int* p_len, int rank) {

	int data_size = get_mpi_dtype_size(mpi_dtype);
	FILE* p_input_file = NULL;
  p_input_file = fopen(file_name, "r");
	//fopen_s(&p_input_file, file_name, "r");
	if (p_input_file == NULL) terminate(rank, "Fail to open input file in function read_replicated_vector.");
	else fread(p_len, sizeof(int), 1, p_input_file);

	*vector = malloc_s(rank, *p_len * data_size * 2);
	fread(*vector, data_size, *p_len, p_input_file);

	fclose(p_input_file);
}


void read_replicated_vector(char* file_name, void** vector, MPI_Datatype mpi_dtype, int* p_len, MPI_Comm comm) {
	int rank;
	int proc_num;
	int data_size = get_mpi_dtype_size(mpi_dtype);
	FILE* p_input_file = NULL;

	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &proc_num);

	if (rank == proc_num - 1) {
     p_input_file = fopen(file_name, "r");
		//fopen_s(&p_input_file, file_name, "r");
		if (p_input_file == NULL) *p_len = 0;
		else fread(p_len, sizeof(int), 1, p_input_file);
	}
	MPI_Bcast(p_len, 1, MPI_INT, proc_num - 1, comm);
	if (!*p_len) terminate(rank, "Fail to open input file in function read_replicated_vector.");

	*vector = malloc_s(rank, *p_len * data_size * 2);

	if (rank == proc_num - 1) {
		fread(*vector, data_size, *p_len, p_input_file);
		fclose(p_input_file);
	}

	MPI_Bcast(*vector, *p_len, mpi_dtype, proc_num - 1, comm);
}

void read_matrix(const char* file_name, void*** p_matrix_2d, void** p_matrix_1d, MPI_Datatype mpi_dtype,
	int* p_row_num, int* p_col_num, int rank) {

	int data_size = get_mpi_dtype_size(mpi_dtype);
	FILE* p_input_file = NULL;
	char* tmp_for_matrix_1d;  // char for byte level operations
	char** tmp_for_matrix_2d;

  p_input_file = fopen(file_name, "r");
	//fopen_s(&p_input_file, file_name, "r");
	if (p_input_file == NULL) terminate(rank, "Fail to open input file in function read_col_striped_matrix.");
	fread(p_row_num, sizeof(int), 1, p_input_file);
	fread(p_col_num, sizeof(int), 1, p_input_file);

	*p_matrix_1d = malloc_s(rank, *p_row_num * *p_col_num * data_size * 2);
	*p_matrix_2d = (void**)malloc_s(rank, *p_row_num * sizeof(void*) * 2);

	tmp_for_matrix_1d = (char*)*p_matrix_1d;
	tmp_for_matrix_2d = (char**)*p_matrix_2d;
	for (int i = 0; i < *p_row_num; i++) {
		*(tmp_for_matrix_2d++) = tmp_for_matrix_1d + i * (*p_col_num * data_size);
	}

	// read by row
	tmp_for_matrix_1d = (char*)*p_matrix_1d;
	for (int i = 0; i < *p_row_num; i++) {
		fread(tmp_for_matrix_1d + i * *p_col_num * data_size,
			data_size, *p_col_num, p_input_file);
	}

	fclose(p_input_file);
}

void read_col_striped_matrix(char* file_name, void*** p_matrix_2d, void** p_matrix_1d,
	MPI_Datatype mpi_dtype,int* p_row_num, int* p_col_num, MPI_Comm comm) {

	int rank;
	int proc_num;
	int data_size = get_mpi_dtype_size(mpi_dtype);
	FILE* p_input_file = NULL;
	int local_col_num;
	char* tmp_for_matrix_1d;  // char for byte level operations
	char** tmp_for_matrix_2d;
	int* send_cnt;
	int* send_disp;
	char* read_buffer;

	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &proc_num);

	if (rank == proc_num - 1) {
     p_input_file = fopen(file_name, "r");
		//fopen_s(&p_input_file, file_name, "r");
		if (p_input_file == NULL) terminate(rank, "Fail to open input file in function read_col_striped_matrix.");
		fread(p_row_num, sizeof(int), 1, p_input_file);
		fread(p_col_num, sizeof(int), 1, p_input_file);
	}
	MPI_Bcast(p_row_num, 1, MPI_INT, proc_num - 1, comm);
	MPI_Bcast(p_col_num, 1, MPI_INT, proc_num - 1, comm);

	local_col_num = BLOCK_SIZE(rank, proc_num, *p_col_num);

	*p_matrix_1d = malloc_s(rank, *p_row_num * local_col_num * data_size * 2);
	*p_matrix_2d = (void**)malloc_s(rank, *p_row_num * sizeof(void*) * 2);

	tmp_for_matrix_1d = (char*)*p_matrix_1d;
	tmp_for_matrix_2d = (char**)*p_matrix_2d;
	for (int i = 0; i < *p_row_num; i++) {
		*(tmp_for_matrix_2d++) = tmp_for_matrix_1d + i * (local_col_num * data_size);
	}

	// read by row
	get_arrays_by_diff_block_size(rank, proc_num, *p_col_num, &send_cnt, &send_disp);
	tmp_for_matrix_1d = (char*)*p_matrix_1d;
	read_buffer = (char*)malloc_s(rank, *p_col_num * data_size * 2);
	for (int i = 0; i < *p_row_num; i++) {
		if (rank == proc_num - 1) {
			fread(read_buffer,
				data_size, *p_col_num, p_input_file);
		}
		MPI_Scatterv(read_buffer, send_cnt, send_disp, mpi_dtype, tmp_for_matrix_1d,
			local_col_num, mpi_dtype, proc_num - 1, comm);
		tmp_for_matrix_1d += local_col_num * data_size;
	}

	free(send_cnt);
	free(send_disp);
	if (rank == proc_num - 1) {
		free(read_buffer);
		fclose(p_input_file);
	}
}


void read_block_matrix(char* file_name, void*** p_matrix_2d, void** p_matrix_1d,
	MPI_Datatype mpi_dtype, int* p_row_num, int* p_col_num, MPI_Comm grid_comm) {

	int grid_rank;
	int proc_num;
	int data_size = get_mpi_dtype_size(mpi_dtype);
	FILE* p_input_file = NULL;
	int local_col_num;
	int local_row_num;
	char* tmp_for_matrix_1d;  // char for byte level operations
	char** tmp_for_matrix_2d;
	/*int* send_cnt;
	int* send_disp;*/
	char* read_buffer;
	char* send_buffer;

	int grid_coord[2];
	int grid_period[2];
	int grid_size[2];  // number of processors at each row/col of grid

	int dest_coord[2];
	int dest_rank;

	int i, j, k;
	MPI_Status* status = (MPI_Status*)malloc(sizeof(MPI_Status));

	MPI_Comm_rank(grid_comm, &grid_rank);
	MPI_Comm_size(grid_comm, &proc_num);

	if (grid_rank == 0) {
     p_input_file = fopen(file_name, "r");
		//fopen_s(&p_input_file, file_name, "r");
		if (p_input_file == NULL) terminate(grid_rank, "Fail to open input file in function read_col_striped_matrix.");
		fread(p_row_num, sizeof(int), 1, p_input_file);
		fread(p_col_num, sizeof(int), 1, p_input_file);
		//printf("utils: Process %d: *p_row_num: %d, *p_col_num: %d\n", grid_rank, *p_row_num, *p_col_num);

	}
	MPI_Bcast(p_row_num, 1, MPI_INT, 0, grid_comm);
	MPI_Bcast(p_col_num, 1, MPI_INT, 0, grid_comm);

	MPI_Cart_get(grid_comm, 2, grid_size, grid_period, grid_coord);
	//printf("utils: Process %d: grid_coord[0]: %d, grid_coord[1]: %d\n", grid_rank, grid_coord[0], grid_coord[1]);


	local_row_num = BLOCK_SIZE(grid_coord[0], grid_size[0], *p_row_num);
	local_col_num = BLOCK_SIZE(grid_coord[1], grid_size[1], *p_col_num);

	*p_matrix_1d = malloc_s(grid_rank, local_row_num * local_col_num * data_size * 2);
	*p_matrix_2d = (void**)malloc_s(grid_rank, local_row_num * sizeof(void*) * 2);

	tmp_for_matrix_1d = (char*)*p_matrix_1d;
	tmp_for_matrix_2d = (char**)*p_matrix_2d;
	for (int i = 0; i < local_row_num; i++) {
		*(tmp_for_matrix_2d++) = tmp_for_matrix_1d + i * (local_col_num * data_size);
	}

	// read by row
	//get_arrays_by_diff_block_size(rank, proc_num, *p_col_num, &send_cnt, &send_disp);
	tmp_for_matrix_1d = (char*)*p_matrix_1d;
	read_buffer = (char*)malloc_s(grid_rank, *p_col_num * data_size * 2);
	for (i = 0; i < grid_size[0]; i++) { // each row of grid
		dest_coord[0] = i;

		for (j = 0; j < BLOCK_SIZE(i, grid_size[0], *p_row_num); j++) { // each row of matrix
			if (grid_rank == 0) {
				fread(read_buffer,
					data_size, *p_col_num, p_input_file);
				print_vector(read_buffer, mpi_dtype, *p_col_num);

			}

			for (k = 0; k < grid_size[1]; k++) {  // each col of grid
				dest_coord[1] = k;
				send_buffer = read_buffer + BLOCK_LOW(k, grid_size[1], *p_col_num) * data_size;

				MPI_Cart_rank(grid_comm, dest_coord, &dest_rank);
				if (grid_rank == 0) {  // root sends data
					if (dest_rank == 0) {
						memcpy((*p_matrix_2d)[j], send_buffer, local_col_num * data_size);
					}
					else {
						MPI_Send(send_buffer, BLOCK_SIZE(k, grid_size[1], *p_col_num),
							mpi_dtype, dest_rank, 0, grid_comm);
					}
				}
				else if (grid_rank == dest_rank) {
					MPI_Recv((*p_matrix_2d)[j], local_col_num, mpi_dtype, 0, 0, grid_comm, status);
				}
			}
		}
	}

	if (grid_rank == 0) {
		free(read_buffer);
		fclose(p_input_file);
	}
	//printf("utils: Process %d: *p_row_num: %d, *p_col_num: %d\n", grid_rank, *p_row_num, *p_col_num);
}


void read_kernel_block_matrix(char* file_name, void*** p_matrix_2d, void** p_matrix_1d,
	MPI_Datatype mpi_dtype, int* p_org_row_num, int* p_org_col_num, int* p_row_num, int* p_col_num,
	int kernel_size, MPI_Comm grid_comm) {

	int grid_rank;
	int proc_num;
	int data_size = get_mpi_dtype_size(mpi_dtype);
	FILE* p_input_file = NULL;
	int local_col_num;
	int local_row_num;
	char* tmp_for_matrix_1d;  // char for byte level operations
	char** tmp_for_matrix_2d;

	int org_row_num, org_col_num;
	int row_cal_num;  // number of calculation, i.e. conv or pooling, along row (dim=0)
	int col_cal_num;  // number of calculation, i.e. conv or pooling, along col (dim=1)
	char* read_buffer = NULL;
	char* send_buffer = NULL;  // pointer to some place in read_buffer
	char** big_matrix_2d = NULL;  // only used in root processor
	char* big_matrix_1d;  // only used in root processor

	int grid_coord[2];
	int grid_period[2];
	int grid_size[2];  // number of processors at each row/col of grid

	int dest_coord[2];
	int dest_rank;

	int i, j, k;
	MPI_Status* status = (MPI_Status*)malloc(sizeof(MPI_Status));

	MPI_Comm_rank(grid_comm, &grid_rank);
	MPI_Comm_size(grid_comm, &proc_num);

	if (grid_rank == 0) {
     p_input_file = fopen(file_name, "r");
		//fopen_s(&p_input_file, file_name, "r");
		if (p_input_file == NULL) terminate(grid_rank, "Fail to open input file in function read_col_striped_matrix.");
		fread(p_org_row_num, sizeof(int), 1, p_input_file);
		fread(p_org_col_num, sizeof(int), 1, p_input_file);

		org_row_num = *p_org_row_num;
		org_col_num = *p_org_col_num;

		printf("utils: Process %d: org_row_num: %d, org_col_num: %d\n", grid_rank, org_row_num, org_col_num);
		fflush(stdout);

		// size of big matrix
		row_cal_num = (org_row_num - kernel_size + 1);
		col_cal_num = (org_col_num - kernel_size + 1);
		*p_row_num = row_cal_num * col_cal_num;
		*p_col_num = kernel_size * kernel_size;

		//printf("utils: Process %d: *p_row_num: %d, *p_col_num: %d\n", grid_rank, *p_row_num, *p_col_num);
		//fflush(stdout);

		big_matrix_1d = (char*)malloc_s(grid_rank, *p_row_num * *p_col_num * data_size * 2);
		big_matrix_2d = (char**)malloc_s(grid_rank, *p_row_num * sizeof(void*) * 2);

		tmp_for_matrix_1d = big_matrix_1d;
		tmp_for_matrix_2d = big_matrix_2d;
		for (i = 0; i < *p_row_num; i++) {
			*(tmp_for_matrix_2d++) = tmp_for_matrix_1d + i * (*p_col_num * data_size);
		}

		// read by row
		//get_arrays_by_diff_block_size(rank, proc_num, *p_col_num, &send_cnt, &send_disp);
		read_buffer = (char*)malloc_s(grid_rank, org_col_num * data_size * 2);
		//read_buffer = (char*)malloc_s(grid_rank, 16 * 16 * data_size);
		int start_id[2];  // start index of a "block" in big matrix; these two number are not used together

		/*fread(read_buffer,
			data_size, 16*16, p_input_file);
		printf("read ");
		print_vector(read_buffer, mpi_dtype, 16*16);
		fflush(stdout);
		return;*/
		for (i = 0; i < org_row_num; i++) { // read each row of original matrix and convert it into big matrix
			fread(read_buffer,
				data_size, org_col_num, p_input_file);
			printf("read ");
			print_vector(read_buffer, mpi_dtype, org_col_num);
			fflush(stdout);


			for (j = i; j >= 0; j--) { // j is for block id, i.e. (j, i-j)
				start_id[0] = j * col_cal_num;
				start_id[1] = (i - j) * kernel_size;

				//printf("utils: Process %d: i: %d, start_id[0]: %d, start_id[1]: %d\n", grid_rank, i, start_id[0], start_id[1]);
				//fflush(stdout);

				if (start_id[0] >= * p_row_num || start_id[1] >= * p_col_num) continue;

				for (k = 0; k < col_cal_num; k++) {
					//printf("utils: Process %d: copy ", grid_rank);
					//print_vector(read_buffer + k * data_size, mpi_dtype, kernel_size);
					//fflush(stdout);
					memcpy(big_matrix_2d[start_id[0] + k] + start_id[1] * data_size,
						read_buffer + k * data_size, kernel_size * data_size);
				}

			}

		}
		/*for (i = 0; i < *p_row_num; i++) {
			print_vector(big_matrix_2d[i], mpi_dtype, *p_col_num);
		}*/
		free(read_buffer);
		fclose(p_input_file);
	}
	
	MPI_Bcast(p_row_num, 1, MPI_INT, 0, grid_comm);
	MPI_Bcast(p_col_num, 1, MPI_INT, 0, grid_comm);

	MPI_Cart_get(grid_comm, 2, grid_size, grid_period, grid_coord);
	//printf("utils: Process %d: grid_coord[0]: %d, grid_coord[1]: %d\n", grid_rank, grid_coord[0], grid_coord[1]);

	local_row_num = BLOCK_SIZE(grid_coord[0], grid_size[0], *p_row_num);
	local_col_num = BLOCK_SIZE(grid_coord[1], grid_size[1], *p_col_num);

	*p_matrix_1d = malloc_s(grid_rank, local_row_num * local_col_num * data_size * 2);
	*p_matrix_2d = (void**)malloc_s(grid_rank, local_row_num * sizeof(void*) * 2);

	tmp_for_matrix_1d = (char*)*p_matrix_1d;
	tmp_for_matrix_2d = (char**)*p_matrix_2d;
	for (int i = 0; i < local_row_num; i++) {
		*(tmp_for_matrix_2d++) = tmp_for_matrix_1d + i * (local_col_num * data_size);
	}

	int cur_row_id = -1;
	for (i = 0; i < grid_size[0]; i++) { // each row of grid
		dest_coord[0] = i;

		for (j = 0; j < BLOCK_SIZE(i, grid_size[0], *p_row_num); j++) { // each row of matrix
			if (grid_rank == 0) {
				cur_row_id += 1;
				read_buffer = big_matrix_2d[cur_row_id];
			}

			for (k = 0; k < grid_size[1]; k++) {  // each col of grid
				dest_coord[1] = k;
				send_buffer = read_buffer + BLOCK_LOW(k, grid_size[1], *p_col_num) * data_size;

				MPI_Cart_rank(grid_comm, dest_coord, &dest_rank);
				if (grid_rank == 0) {  // root sends data
					if (dest_rank == 0) {
						memcpy((*p_matrix_2d)[j], send_buffer, local_col_num * data_size);
					}
					else {
						MPI_Send(send_buffer, BLOCK_SIZE(k, grid_size[1], *p_col_num),
							mpi_dtype, dest_rank, 0, grid_comm);
					}
				}
				else if (grid_rank == dest_rank) {
					MPI_Recv((*p_matrix_2d)[j], local_col_num, mpi_dtype, 0, 0, grid_comm, status);
				}
			}
		}
	}

	//printf("utils: Process %d: *p_row_num: %d, *p_col_num: %d\n", grid_rank, *p_row_num, *p_col_num);
}

void generate_kernel_block_matrix(void*** p_matrix_2d, void** p_matrix_1d,
	MPI_Datatype mpi_dtype, int* p_org_row_num, int* p_org_col_num, int* p_row_num, int* p_col_num,
	int kernel_size, MPI_Comm grid_comm) {

	int grid_rank;
	int proc_num;
	int data_size = get_mpi_dtype_size(mpi_dtype);
	//FILE* p_input_file = NULL;
	int local_col_num;
	int local_row_num;
	char* tmp_for_matrix_1d;  // char for byte level operations
	char** tmp_for_matrix_2d;

	int org_row_num, org_col_num;
	int row_cal_num;  // number of calculation, i.e. conv or pooling, along row (dim=0)
	int col_cal_num;  // number of calculation, i.e. conv or pooling, along col (dim=1)
	char* read_buffer = NULL;
	char* send_buffer = NULL;  // pointer to some place in read_buffer
	double* gen_buffer = NULL;

	char** big_matrix_2d = NULL;  // only used in root processor
	char* big_matrix_1d;  // only used in root processor

	int grid_coord[2];
	int grid_period[2];
	int grid_size[2];  // number of processors at each row/col of grid

	int dest_coord[2];
	int dest_rank;

	int seed = 1234;

	int i, j, k;
	MPI_Status* status = (MPI_Status*)malloc(sizeof(MPI_Status));

	MPI_Comm_rank(grid_comm, &grid_rank);
	MPI_Comm_size(grid_comm, &proc_num);

	if (grid_rank == 0) {
		/*fopen_s(&p_input_file, file_name, "r");
		if (p_input_file == NULL) terminate(grid_rank, "Fail to open input file in function read_col_striped_matrix.");
		fread(p_org_row_num, sizeof(int), 1, p_input_file);
		fread(p_org_col_num, sizeof(int), 1, p_input_file);*/
		*p_org_row_num = 10;//1024
		*p_org_col_num = 10;

		org_row_num = *p_org_row_num;
		org_col_num = *p_org_col_num;

		printf("utils: Process %d: org_row_num: %d, org_col_num: %d\n", grid_rank, org_row_num, org_col_num);
		fflush(stdout);

		// size of big matrix
		row_cal_num = (org_row_num - kernel_size + 1);
		col_cal_num = (org_col_num - kernel_size + 1);
		*p_row_num = row_cal_num * col_cal_num;
		*p_col_num = kernel_size * kernel_size;

		//printf("utils: Process %d: *p_row_num: %d, *p_col_num: %d\n", grid_rank, *p_row_num, *p_col_num);
		//fflush(stdout);

		big_matrix_1d = (char*)malloc_s(grid_rank, *p_row_num * *p_col_num * data_size * 2);
		big_matrix_2d = (char**)malloc_s(grid_rank, *p_row_num * sizeof(void*) * 2);

		tmp_for_matrix_1d = big_matrix_1d;
		tmp_for_matrix_2d = big_matrix_2d;
		for (i = 0; i < *p_row_num; i++) {
			*(tmp_for_matrix_2d++) = tmp_for_matrix_1d + i * (*p_col_num * data_size);
		}

		// read by row
		//get_arrays_by_diff_block_size(rank, proc_num, *p_col_num, &send_cnt, &send_disp);
		read_buffer = (char*)malloc_s(grid_rank, org_col_num * data_size * 2);
		gen_buffer = (double*)malloc_s(grid_rank, org_col_num * data_size * 2);

		int start_id[2];  // start index of a "block" in big matrix; these two number are not used together
		/*for (int l = 0; l < org_col_num; l++) {
			printf("LOOK!!!");
			gen_buffer[l] = 1.0;
		}
		printf("OVER!!!");
		print_vector(read_buffer, mpi_dtype, org_col_num);
		fflush(stdout);*/

		for (i = 0; i < org_row_num; i++) { // read each row of original matrix and convert it into big matrix
			for (int l = 0; l < org_col_num; l++) {
				gen_buffer[l] = (double(rand() % 100)) / 10;
				memcpy(read_buffer + l * data_size, gen_buffer + l, data_size);
			}

			/*printf("read ");
			print_vector(read_buffer, mpi_dtype, org_col_num);
			fflush(stdout);*/


			for (j = i; j >= 0; j--) { // j is for block id, i.e. (j, i-j)
				start_id[0] = j * col_cal_num;
				start_id[1] = (i - j) * kernel_size;

				//printf("utils: Process %d: i: %d, start_id[0]: %d, start_id[1]: %d\n", grid_rank, i, start_id[0], start_id[1]);
				//fflush(stdout);

				if (start_id[0] >= *p_row_num || start_id[1] >= *p_col_num) continue;

				for (k = 0; k < col_cal_num; k++) {
					//printf("utils: Process %d: copy ", grid_rank);
					//print_vector(read_buffer + k * data_size, mpi_dtype, kernel_size);
					//fflush(stdout);
					memcpy(big_matrix_2d[start_id[0] + k] + start_id[1] * data_size,
						read_buffer + k * data_size, kernel_size * data_size);
				}

			}

		}
		/*for (i = 0; i < *p_row_num; i++) {
			print_vector(big_matrix_2d[i], mpi_dtype, *p_col_num);
		}*/
		free(read_buffer);
		free(gen_buffer);
		//fclose(p_input_file);
	}
	MPI_Bcast(p_row_num, 1, MPI_INT, 0, grid_comm);
	MPI_Bcast(p_col_num, 1, MPI_INT, 0, grid_comm);

	MPI_Cart_get(grid_comm, 2, grid_size, grid_period, grid_coord);
	//printf("utils: Process %d: grid_coord[0]: %d, grid_coord[1]: %d\n", grid_rank, grid_coord[0], grid_coord[1]);
	

	local_row_num = BLOCK_SIZE(grid_coord[0], grid_size[0], *p_row_num);
	local_col_num = BLOCK_SIZE(grid_coord[1], grid_size[1], *p_col_num);

	*p_matrix_1d = malloc_s(grid_rank, local_row_num * local_col_num * data_size * 2);
	*p_matrix_2d = (void**)malloc_s(grid_rank, local_row_num * sizeof(void*) * 2);
	

	tmp_for_matrix_1d = (char*)*p_matrix_1d;
	tmp_for_matrix_2d = (char**)*p_matrix_2d;
	for (int i = 0; i < local_row_num; i++) {
		*(tmp_for_matrix_2d++) = tmp_for_matrix_1d + i * (local_col_num * data_size);
	}


	int cur_row_id = -1;
	for (i = 0; i < grid_size[0]; i++) { // each row of grid
		dest_coord[0] = i;

		for (j = 0; j < BLOCK_SIZE(i, grid_size[0], *p_row_num); j++) { // each row of matrix
			if (grid_rank == 0) {
				cur_row_id += 1;
				read_buffer = big_matrix_2d[cur_row_id];
			}

			for (k = 0; k < grid_size[1]; k++) {  // each col of grid
				dest_coord[1] = k;
				send_buffer = read_buffer + BLOCK_LOW(k, grid_size[1], *p_col_num) * data_size;

				MPI_Cart_rank(grid_comm, dest_coord, &dest_rank);
				if (grid_rank == 0) {  // root sends data
					if (dest_rank == 0) {
						memcpy((*p_matrix_2d)[j], send_buffer, local_col_num * data_size);
					}
					else {
						MPI_Send(send_buffer, BLOCK_SIZE(k, grid_size[1], *p_col_num),
							mpi_dtype, dest_rank, 0, grid_comm);
					}
				}
				else if (grid_rank == dest_rank) {
					MPI_Recv((*p_matrix_2d)[j], local_col_num, mpi_dtype, 0, 0, grid_comm, status);
				}
			}
		}
	}


	//printf("utils: Process %d: *p_row_num: %d, *p_col_num: %d\n", grid_rank, *p_row_num, *p_col_num);
}


void read_block_kernel(char* file_name, void*** p_matrix_2d, void** p_matrix_1d,
	MPI_Datatype mpi_dtype, int* p_kernel_size, int* p_kernel_num, MPI_Comm grid_comm) {
	/* kernels will form a new matrix [kernel_num, kernel_size * kernel_size],
	 * then transpose to [kernel_size * kernel_size, kernel_num]
	 */


	int grid_rank;
	int proc_num;
	int data_size = get_mpi_dtype_size(mpi_dtype);
	FILE* p_input_file = NULL;
	int local_col_num;
	char* tmp_for_matrix_1d;  // char for byte level operations
	char** tmp_for_matrix_2d;

	char* read_buffer = NULL;
	char* send_buffer = NULL;  // pointer to some place in read_buffer

	int grid_coord[2];
	int grid_period[2];
	int grid_size[2];  // number of processors at each row/col of grid
	int kernel_matrix_row_num, kernel_matrix_col_num;

	int dest_coord[2];
	int dest_rank;

	int i, j, k;
	MPI_Status* status = (MPI_Status*)malloc(sizeof(MPI_Status));

	MPI_Comm_rank(grid_comm, &grid_rank);
	MPI_Comm_size(grid_comm, &proc_num);

	if (grid_rank == 0) {
    p_input_file = fopen(file_name, "r");
    printf("hello");
		//fopen_s(&p_input_file, file_name, "r");
		if (p_input_file == NULL) terminate(grid_rank, "Fail to open input file in function read_col_striped_matrix.");
		fread(p_kernel_num, sizeof(int), 1, p_input_file);
		fread(p_kernel_size, sizeof(int), 1, p_input_file);
		printf("utils: Process %d: *p_kernel_num: %d, *p_kernel_size: %d\n", grid_rank, *p_kernel_num, *p_kernel_size);
		//fflush(stdout);
	}

	MPI_Bcast(p_kernel_num, 1, MPI_INT, 0, grid_comm);
	MPI_Bcast(p_kernel_size, 1, MPI_INT, 0, grid_comm);

	MPI_Cart_get(grid_comm, 2, grid_size, grid_period, grid_coord);
	//printf("utils: Process %d: grid_coord[0]: %d, grid_coord[1]: %d\n", grid_rank, grid_coord[0], grid_coord[1]);

	kernel_matrix_row_num = *p_kernel_num;
	kernel_matrix_col_num = *p_kernel_size * *p_kernel_size;
	local_col_num = BLOCK_SIZE(grid_coord[1], grid_size[1], kernel_matrix_col_num);

	*p_matrix_1d = malloc_s(grid_rank, kernel_matrix_row_num * local_col_num * data_size * 2);
	*p_matrix_2d = (void**)malloc_s(grid_rank, kernel_matrix_row_num * sizeof(void*) * 2);

	tmp_for_matrix_1d = (char*)*p_matrix_1d;
	tmp_for_matrix_2d = (char**)*p_matrix_2d;
	for (int i = 0; i < kernel_matrix_row_num; i++) {
		*(tmp_for_matrix_2d++) = tmp_for_matrix_1d + i * (local_col_num * data_size);
	}

	
	if (grid_rank == 0) {
		printf("Process %d: reading kernels...\n", grid_rank);
		fflush(stdout);
		read_buffer = (char*)malloc_s(grid_rank, kernel_matrix_col_num * data_size * 2);
	}
	//printf("Process %d: kernel_matrix_row_num is %d\n", grid_rank, kernel_matrix_row_num);
	//printf("Process %d: kernel_matrix_col_num is %d\n", grid_rank, kernel_matrix_col_num);
	//fflush(stdout);

	for (j = 0; j < kernel_matrix_row_num; j++) {
		// each row of matrix of a processor (every processor has same row number here)
		/*printf("Process %d: here in j loop\n", grid_rank);
		fflush(stdout);*/
		if (grid_rank == 0) {
			/*printf("Process %d: before read \n", grid_rank);
			fflush(stdout);*/
			fread(read_buffer,
				data_size, kernel_matrix_col_num, p_input_file);
			/*printf("Process %d: after read \n", grid_rank);
			fflush(stdout);*/
			//print_vector(read_buffer, mpi_dtype, kernel_matrix_col_num);
			//fflush(stdout);
		}

		for (k = 0; k < grid_size[1]; k++) {  // each col of grid
			dest_coord[1] = k;
			send_buffer = read_buffer + BLOCK_LOW(k, grid_size[1], kernel_matrix_col_num) * data_size;

			for (i = 0; i < grid_size[0]; i++) { // each row of grid
				dest_coord[0] = i;

				MPI_Cart_rank(grid_comm, dest_coord, &dest_rank);
				if (grid_rank == 0) {  // root sends data
					if (dest_rank == 0) {
						memcpy((*p_matrix_2d)[j], send_buffer, local_col_num * data_size);
					}
					else {
						MPI_Send(send_buffer, BLOCK_SIZE(k, grid_size[1], kernel_matrix_col_num),
							mpi_dtype, dest_rank, 0, grid_comm);
					}
				}
				else if (grid_rank == dest_rank) {
					MPI_Recv((*p_matrix_2d)[j], local_col_num, mpi_dtype, 0, 0, grid_comm, status);
				}
			}
		}
	}

	// transpose: [kernel_num, kernel_size * kernel_size] -> [kernel_size * kernel_size, kernel_num]
	printf("Process %d: transposing...", grid_rank);
	fflush(stdout);

	char** matrix_t_2d = (char**)malloc_s(grid_rank, local_col_num * sizeof(void**) * 2);
	char* matrix_t_1d = (char*)malloc_s(grid_rank, local_col_num * kernel_matrix_row_num * data_size * 2);

	tmp_for_matrix_1d = matrix_t_1d;
	tmp_for_matrix_2d = matrix_t_2d;
	for (i = 0; i < local_col_num; i++) {
		*(tmp_for_matrix_2d++) = tmp_for_matrix_1d + i * (kernel_matrix_row_num * data_size);
	}

	for (i = 0; i < local_col_num; i++) {
		for (j = 0; j < kernel_matrix_row_num; j++) {
			memcpy(matrix_t_1d + (i * kernel_matrix_row_num + j) * data_size,
				(char*)*p_matrix_1d + (j * kernel_matrix_row_num + i) * data_size, data_size);
		}
	}
	memcpy(*p_matrix_2d, matrix_t_2d, local_col_num * sizeof(void**));
	memcpy(*p_matrix_1d, matrix_t_1d, local_col_num * kernel_matrix_row_num * data_size);

	
	/*for (i = 0; i < local_col_num; i++) {
		printf("Process %d: ", grid_rank);
		print_vector((*p_matrix_2d)[i], mpi_dtype, kernel_matrix_row_num);
		fflush(stdout);
	}*/

	if (grid_rank == 0) {
		free(read_buffer);
		fclose(p_input_file);
	}
	//printf("utils: Process %d: *p_row_num: %d, *p_col_num: %d\n", grid_rank, *p_row_num, *p_col_num);
}


void print_vector(void* vector, MPI_Datatype mpi_dtype, int len) {
	for (int i = 0; i < len; i++) {
		switch (mpi_dtype)
		{
		case MPI_INT:
			printf("%6d ", ((int*)vector)[i]);
			break;
		case MPI_FLOAT:
			printf("%6.3f ", ((float*)vector)[i]);
			break;
		case MPI_DOUBLE:
			printf("%6.3f ", ((double*)vector)[i]);
			break;
		}
	}
	printf("\n");
}

void print_replicated_vector(void* vector, MPI_Datatype mpi_dtype, int len, MPI_Comm comm) {
	int rank;
	MPI_Comm_rank(comm, &rank);
	if (!rank) {
		print_vector(vector, mpi_dtype, len);
	}
}

void print_col_striped_matrix(void** matrix_2d, MPI_Datatype mpi_dtype, int row_num, int col_num, MPI_Comm comm) {
	int data_size = get_mpi_dtype_size(mpi_dtype);
	int rank;
	int proc_num;
	int* recv_cnt;
	int* recv_disp;
	void* coll_buffer = NULL;
	int local_col_num;

	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &proc_num);

	local_col_num = BLOCK_SIZE(rank, proc_num, col_num);
	get_arrays_by_diff_block_size(rank, proc_num, col_num, &recv_cnt, &recv_disp);
	/*print_vector(recv_cnt, MPI_INT, proc_num);
	print_vector(recv_disp, MPI_INT, proc_num);*/

	if (!rank) coll_buffer = malloc_s(rank, col_num * data_size * 2);

	for (int i = 0; i < row_num; i++) {
		MPI_Gatherv(matrix_2d[i], local_col_num, mpi_dtype, coll_buffer, recv_cnt, recv_disp, mpi_dtype, 0, comm);
		if (!rank) print_vector(coll_buffer, mpi_dtype, col_num);
	}

	free(recv_cnt);
	free(recv_disp);
	if (rank == proc_num - 1) {
		free(coll_buffer);
	}
}


void print_block_matrix(void** matrix_2d, MPI_Datatype mpi_dtype, int row_num, int col_num, MPI_Comm grid_comm) {
	int data_size = get_mpi_dtype_size(mpi_dtype);
	int grid_rank;
	int proc_num;
	char* coll_buffer = NULL;  // buffer for collection
	char* recv_buffer;  // pointer to some location of coll_buffer
	int local_col_num;

	int grid_coord[2];
	int grid_period[2];
	int grid_size[2];  // number of processors at each row/col of grid

	int source_coord[2];
	int source_rank;
	MPI_Status* status = (MPI_Status*)malloc(sizeof(MPI_Status));

	int i, j, k;

	MPI_Comm_rank(grid_comm, &grid_rank);
	MPI_Comm_size(grid_comm, &proc_num);

	MPI_Cart_get(grid_comm, 2, grid_size, grid_period, grid_coord);

	local_col_num = BLOCK_SIZE(grid_coord[1], grid_size[1], col_num);

	/*print_vector(recv_cnt, MPI_INT, proc_num);
	print_vector(recv_disp, MPI_INT, proc_num);*/

	if (!grid_rank) coll_buffer = (char*)malloc_s(grid_rank, col_num * data_size * 2);

	for (i = 0; i < grid_size[0]; i++) { // each row of grid
		source_coord[0] = i;

		for (j = 0; j < BLOCK_SIZE(i, grid_size[0], row_num); j++) { // each row of matrix

			if (grid_rank == 0) {  // root processor will receive data
				for (k = 0; k < grid_size[1]; k++) {  // each col of grid
					source_coord[1] = k;
					recv_buffer = coll_buffer + BLOCK_LOW(k, grid_size[1], col_num) * data_size;

					MPI_Cart_rank(grid_comm, source_coord, &source_rank);
					if (grid_rank == 0) {  // root receives data
						if (source_rank == 0) {
							memcpy(recv_buffer, matrix_2d[j], local_col_num * data_size);
						}
						else {
							MPI_Recv(recv_buffer, BLOCK_SIZE(k, grid_size[1], col_num),
								mpi_dtype, source_rank, 0, grid_comm, status);
						}
					}
				}
				print_vector(coll_buffer, mpi_dtype, col_num);
			}
			else if (grid_coord[0] == i) {  // others will send data
				MPI_Send(matrix_2d[j], local_col_num, mpi_dtype, 0, 0, grid_comm);
			}
			
		}
	}

	if (grid_rank == 0) {
		free(coll_buffer);
	}
}