#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <string.h>


void my_allgather(float* send_data, int send_count, MPI_Datatype send_datatype,
                    float* recv_data, int recv_count, MPI_Datatype recv_datatype,
                    MPI_Comm communicator) {
    int world_rank;
    MPI_Comm_rank(communicator, &world_rank);
    int world_size;
    MPI_Comm_size(communicator, &world_size);

   
    for (int i = 0; i < world_size; i++) {
        if (i != world_rank) {
            MPI_Send(send_data, send_count, send_datatype, i, 1, MPI_COMM_WORLD);
            MPI_Recv(recv_data + i * recv_count , recv_count, recv_datatype, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    // self data
    memcpy(recv_data + world_rank * recv_count, send_data, sizeof(send_datatype)*send_count);

    }
}



// Creates an array of random numbers. Each number has a value from 0 - 1
float *create_rand_nums(int num_elements) {
    float *rand_nums = (float *)malloc(sizeof(float) * num_elements);
    assert(rand_nums != NULL);
    int i;
    for (i = 0; i < num_elements; i++) {
        rand_nums[i] = (rand() / (float)RAND_MAX);
    }
    return rand_nums;
}

// Computes the average of an array of numbers
float compute_avg(float *array, int num_elements) {
    float sum = 0.f;
    int i;
    for (i = 0; i < num_elements; i++) {
        sum += array[i];
    }
    return sum / num_elements;
}


int main(int argc, char** argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: avg num_elements_per_proc\n");
        exit(1);
    }

    int num_elements_per_proc = atoi(argv[1]);

    // Seed the random number generator to get different results each time
    srand(time(NULL));

    MPI_Init(NULL, NULL);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Create a random array of elements on the root process. Its total
    // size will be the number of elements per process times the number
    // of processes
    float *rand_nums = NULL;
    if (world_rank == 0) {
        rand_nums = create_rand_nums(num_elements_per_proc * world_size);
    }

    // For each process, create a buffer that will hold a subset of the entire
    // array
    float *sub_rand_nums = (float *)malloc(sizeof(float) * num_elements_per_proc);
    assert(sub_rand_nums != NULL);

    // Scatter the random numbers from the root process to all processes in
    // the MPI world
    MPI_Scatter(rand_nums, num_elements_per_proc, MPI_FLOAT, sub_rand_nums,
                num_elements_per_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Compute the average of your subset
    float sub_avg = compute_avg(sub_rand_nums, num_elements_per_proc);

    // Gather all partial averages down to all the processes
    float *sub_avgs1 = (float *)malloc(sizeof(float) * world_size);
    assert(sub_avgs1 != NULL);
    float *sub_avgs2 = (float *)malloc(sizeof(float) * world_size);
    assert(sub_avgs2 != NULL);

    double total_my_allgather_time = 0.0;
    double total_mpi_allgather_time = 0.0;

    MPI_Barrier(MPI_COMM_WORLD);
    total_mpi_allgather_time -= MPI_Wtime();
    MPI_Allgather(&sub_avg, 1, MPI_FLOAT, sub_avgs1, 1, MPI_FLOAT,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    total_mpi_allgather_time += MPI_Wtime();

    total_my_allgather_time -= MPI_Wtime();
    my_allgather(&sub_avg, 1, MPI_FLOAT, sub_avgs2, 1, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    total_my_allgather_time += MPI_Wtime();


    // print results
    float avg1 = compute_avg(sub_avgs1, world_size);
    printf("[MPI_AllGather]Avg of all elements from processor %d is %f\n", world_rank, avg1);
    float avg2 = compute_avg(sub_avgs2, world_size);
    printf("[my_allgather]Avg of all elements from processor %d is %f\n", world_rank, avg2);


    if (world_rank == 0) {
        // print time
        printf("Total mpi_allgather time = %lf\n", total_mpi_allgather_time);
        printf("Total my_allgather time = %lf\n", total_my_allgather_time);
        // Clean up
        free(rand_nums);

    }
    free(sub_avgs1);
    free(sub_avgs2);
    free(sub_rand_nums);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
