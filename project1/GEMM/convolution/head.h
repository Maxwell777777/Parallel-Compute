#pragma once


#define DATA_MSG 0
#define PROMPT_MSG 1
#define RESPONSE_MSG 2

#define OPEN_FILE_ERROR -1
#define MALLOC_ERROR -2
#define TYPE_ERROR -3

#define PTR_SIZE sizeof(void*)

#define MAX(a,b) ((a)<(b)?(b):(a))
#define MIN(a,b) ((a)>(b)?(b):(a))

#define BLOCK_LOW(i,p,n) ((i)*(n)/(p))  // i for index of process
#define BLOCK_HIGH(i,p,n) (BLOCK_LOW(i+1,p,n)-1)
#define BLOCK_SIZE(i,p,n) (BLOCK_LOW(i+1,p,n)-BLOCK_LOW(i,p,n))
#define BLOCK_OWNER(i,p,n) (((p)*(i+1)-1)/n)  // i for index of block

int get_mpi_dtype_size(MPI_Datatype mpi_dtype);
int terminate(int rank, const char* msg);
void* malloc_s(int rank, int num_bytes);

void get_arrays_by_diff_block_size(int rank, int proc_num, int elem_num, int** cnt_array, int** disp_array);
void get_arrays_by_curr_block_size(int rank, int proc_num, int elem_num, int** cnt_array, int** disp_array);

void read_vector(char* file_name, void** vector, MPI_Datatype mpi_dtype, int* p_len, int rank);

void read_replicated_vector(char* file_name, void** vector, MPI_Datatype mpi_dtype, int* p_len, MPI_Comm comm);

void read_matrix(const char* file_name, void*** p_matrix_2d, void** p_matrix_1d, MPI_Datatype mpi_dtype,
	int* p_row_num, int* p_col_num, int rank);

void read_col_striped_matrix(char* file_name, void*** p_matrix_2d, void** p_matrix_1d, MPI_Datatype mpi_dtype,
	int* p_row_num, int* p_col_num, MPI_Comm comm);

void read_block_matrix(char* file_name, void*** p_matrix_2d, void** p_matrix_1d,
	MPI_Datatype mpi_dtype, int* p_row_num, int* p_col_num, MPI_Comm grid_comm);

void read_kernel_block_matrix(char* file_name, void*** p_matrix_2d, void** p_matrix_1d,
	MPI_Datatype mpi_dtype, int* p_org_row_num, int* p_org_col_num, int* p_row_num, int* p_col_num,
	int kernel_size, MPI_Comm grid_comm);

void generate_kernel_block_matrix(void*** p_matrix_2d, void** p_matrix_1d,
	MPI_Datatype mpi_dtype, int* p_org_row_num, int* p_org_col_num, int* p_row_num, int* p_col_num,
	int kernel_size, MPI_Comm grid_comm);

void read_block_kernel(char* file_name, void*** p_matrix_2d, void** p_matrix_1d, MPI_Datatype mpi_dtype, int* p_kernel_size, int* p_kernel_num, MPI_Comm grid_comm);

void print_vector(void* vector, MPI_Datatype mpi_dtype, int len);

void print_replicated_vector(void* vector, MPI_Datatype mpi_dtype, int len, MPI_Comm comm);

void print_col_striped_matrix(void** matrix_2d, MPI_Datatype mpi_dtype, int row_num, int col_num, MPI_Comm comm);

void print_block_matrix(void** matrix_2d, MPI_Datatype mpi_dtype, int row_num, int col_num, MPI_Comm grid_comm);
