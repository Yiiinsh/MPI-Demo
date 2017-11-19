#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#include "pgmio.h"
#include "arraytool.h"

/* Check master rank */
#define is_master(rank) (0 == rank)

/* Processing related */
#define MAX_LOOP 1500
#define THRESHOLD 0.05
/* MPI related */
#define RANK_MASTER 0
#define DEFAULT_TAG 1
#define LEFT_SEND_TAG 2
#define RIGHT_SEND_TAG 3
#define LEFT_FROM_TAG 3
#define RIGHT_FROM_TAG 2
#define NDIMS 2
#define DIM_X 0
#define DIM_Y 1

/* Check arguments for user input */
static void check_arg(int argc, char **argv, char **input_file_name, char **output_file_name);
/* Boundary value for top & buttom */
static double boundaryval(int i, int m);


char *input_file_name, *output_file_name; // file source & output target

int main(int argc, char **argv)
{
    check_arg(argc, argv, &input_file_name, &output_file_name);

    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // Get file size
    int length, width;
    pgmsize(input_file_name, &length, &width);
    if (is_master(rank))
    {
        fprintf(stdout, "Original file size: %d x %d.\n", length, width);
    }

    // Virtual topology
    MPI_Comm cart_comm;
    int *dims = (int *)malloc(NDIMS * sizeof(int));
    int *periods = (int *)malloc(NDIMS * sizeof(int));
    periods[DIM_X] = 1;
    periods[DIM_Y] = 0;
    MPI_Dims_create(size, NDIMS, dims);
    MPI_Cart_create(comm, NDIMS, dims, periods, 0, &cart_comm);
    int coords[NDIMS];
    MPI_Cart_coords(cart_comm, rank, NDIMS, coords);
    if (is_master(rank))
    {
        fprintf(stdout, "Mapping processes into %d x %d grids.\n", dims[DIM_X], dims[DIM_Y]);
    }
    fprintf(stdout, "Rank %d in coords (%d , %d)\n", rank, coords[DIM_X], coords[DIM_Y]);

    // Get neighbor
    int rank_up, rank_down, rank_left, rank_right;
    MPI_Cart_shift(cart_comm, 0, 1, &rank_left, &rank_right);
    MPI_Cart_shift(cart_comm, 1, 1, &rank_down, &rank_up);
    fprintf(stdout, "Neighbors of rank %d : Left %d Right %d Up %d Down %d\n", rank, rank_left, rank_right, rank_up, rank_down);

    // Start
    if (is_master(rank))
    {
        fprintf(stdout, "Start...\n");
    }
    MPI_Barrier(comm);
    double start = MPI_Wtime();

    // Memory allocation
    int x_step = length / dims[DIM_X];
    int y_step = width / dims[DIM_Y];
    int plength = x_step;
    if (coords[0] == dims[0] - 1)
    {
        plength = length / dims[DIM_X] + length % dims[DIM_X];
    }
    int pwidth = y_step;
    if (coords[1] == dims[1] - 1)
    {
        pwidth = width / dims[DIM_Y] + width % dims[DIM_Y];
    }
    fprintf(stdout, "rank %d plength %d, pwidth %d\n", rank, plength, pwidth);

    double **masterbuf, **img_new, **img_old, **img_edge;
    double *master_content, *new_content, *old_content, *edge_content;
    double_2d_array_allocation(&masterbuf, &master_content, length, width);
    double_2d_array_allocation(&img_new, &new_content, plength + 2, pwidth + 2);
    double_2d_array_allocation(&img_old, &old_content, plength + 2, pwidth + 2);
    double_2d_array_allocation(&img_edge, &edge_content, plength + 2, pwidth + 2);

    fprintf(stdout, "Rank %d reading from file %s...\n", rank, input_file_name);
    part_pgmread(input_file_name, &(img_edge[1][1]), plength, pwidth, coords[0] * x_step, coords[1] * y_step, 2);

    // Processing
    for (int i = 0; i != plength + 2; ++i)
    {
        for (int j = 0; j != pwidth + 2; ++j)
        {
            img_old[i][j] = 255.0;
        }
    }
    if (MPI_PROC_NULL == rank_up)
    {
        for (int i = 1; i <= plength; ++i)
        {
            double val = boundaryval(coords[0] * x_step + i, length);
            // img_old[i][0] = (int)(255.0 * val);
            img_old[i][pwidth + 1] = (int)(255.0 * (1.0 - val));
        }
    }
    if (MPI_PROC_NULL == rank_down)
    {
        for (int i = 1; i <= plength; ++i)
        {
            double val = boundaryval(coords[0] * x_step + i, length);
            // img_old[i][pwidth + 1] = (int)(255.0 * (1.0 - val));
            img_old[i][0] = (int)(255.0 * val);
        }
    }

    fprintf(stdout, "Start iteration on rank %d\n", rank);
    MPI_Datatype column_type, row_type;
    MPI_Type_vector(pwidth + 2, 1, 1, MPI_DOUBLE, &column_type);
    MPI_Type_commit(&column_type);
    MPI_Type_vector(plength + 2, 1, pwidth + 2, MPI_DOUBLE, &row_type);
    MPI_Type_commit(&row_type);

    MPI_Status status, left_recv_status, right_recv_status, down_recv_status, up_recv_status;
    MPI_Request left_send_request, right_send_request, down_send_request, up_send_request;
    for (int cnt = 0; cnt != MAX_LOOP; ++cnt)
    {
        MPI_Issend(&(img_old[1][0]), 1, column_type, rank_left, LEFT_SEND_TAG, comm, &left_send_request);
        MPI_Issend(&(img_old[plength][0]), 1, column_type, rank_right, RIGHT_SEND_TAG, comm, &right_send_request);
        MPI_Issend(&(img_old[0][1]), 1, row_type, rank_down, DEFAULT_TAG, comm, &down_send_request);
        MPI_Issend(&(img_old[0][pwidth]), 1, row_type, rank_up, DEFAULT_TAG, comm, &up_send_request);
        for (int i = 2; i <= plength - 1; ++i)
        {
            for (int j = 2; j <= pwidth - 1; ++j)
            {
                img_new[i][j] = 0.25 * (img_old[i - 1][j] + img_old[i + 1][j] + img_old[i][j - 1] + img_old[i][j + 1] - img_edge[i][j]);
            }
        }
        MPI_Recv(&(img_old[0][0]), 1, column_type, rank_left, LEFT_FROM_TAG, comm, &left_recv_status);
        MPI_Recv(&(img_old[plength + 1][0]), 1, column_type, rank_right, RIGHT_FROM_TAG, comm, &right_recv_status);
        MPI_Recv(&(img_old[0][0]), 1, row_type, rank_down, DEFAULT_TAG, comm, &down_recv_status);
        MPI_Recv(&(img_old[0][pwidth + 1]), 1, row_type, rank_up, DEFAULT_TAG, comm, &up_recv_status);
        for (int j = 1; j <= pwidth; ++j)
        {
            img_new[1][j] = 0.25 * (img_old[1 - 1][j] + img_old[1 + 1][j] + img_old[1][j - 1] + img_old[1][j + 1] - img_edge[1][j]);
            img_new[plength][j] = 0.25 * (img_old[plength - 1][j] + img_old[plength + 1][j] + img_old[plength][j - 1] + img_old[plength][j + 1] - img_edge[plength][j]);
        }
        for (int i = 1; i <= plength; ++i)
        {
            img_new[i][1] = 0.25 * (img_old[i - 1][1] + img_old[i + 1][1] + img_old[i][1 - 1] + img_old[i][1 + 1] - img_edge[i][1]);
            img_new[i][pwidth] = 0.25 * (img_old[i - 1][pwidth] + img_old[i + 1][pwidth] + img_old[i][pwidth - 1] + img_old[i][pwidth + 1] - img_edge[i][pwidth]);
        }
        MPI_Wait(&left_send_request, &status);
        MPI_Wait(&right_send_request, &status);
        MPI_Wait(&up_send_request, &status);
        MPI_Wait(&down_send_request, &status);

        double delta = 0, overall_delta = 0;
        for (int i = 1; i <= plength; ++i)
        {
            for (int j = 1; j <= pwidth; ++j)
            {
                delta = fabs(img_new[i][j] - img_old[i][j]) > delta ? fabs(img_new[i][j] - img_old[i][j]) : delta;
                img_old[i][j] = img_new[i][j];
            }
        }
        MPI_Allreduce(&delta, &overall_delta, 1, MPI_DOUBLE, MPI_MAX, comm);
        if (overall_delta <= THRESHOLD)
        {
            if (is_master(rank))
            {
                printf("Finish at iteration %d, with delta : %.5f\n", cnt, overall_delta);
            }
            break;
        }
    }

    MPI_Datatype send_type;
    MPI_Type_vector(plength, pwidth, pwidth + 2, MPI_DOUBLE, &send_type);
    MPI_Type_commit(&send_type);
    if (!is_master(rank))
    {
        // Other send to master
        MPI_Ssend(&(img_old[1][1]), 1, send_type, RANK_MASTER, DEFAULT_TAG, comm);
    }
    MPI_Type_free(&send_type);
    
    if (is_master(rank))
    {
        // Init grid of master
        for(int i = 0; i != plength; ++i) {
            for(int j = 0; j != pwidth; ++j) {
                masterbuf[i][j] = img_old[i + 1][j + 1];
            }
        }

        // Recv other grids
        for(int r = 1; r != size; ++r) {
            int recv_coords[NDIMS];
            MPI_Cart_coords(cart_comm, r, NDIMS, recv_coords);
            int recv_plength = (recv_coords[DIM_X] == dims[DIM_X] - 1) ? x_step + length % dims[DIM_X] : x_step;
            int recv_pwidth = (recv_coords[DIM_Y] == dims[DIM_Y] - 1) ? y_step + width % dims[DIM_Y] : y_step;

            MPI_Datatype recv_type;
            MPI_Type_vector(recv_plength, recv_pwidth, width, MPI_DOUBLE, &recv_type);
            MPI_Type_commit(&recv_type);
            MPI_Recv(&(masterbuf[recv_coords[DIM_X] * x_step][recv_coords[DIM_Y] * y_step]), 1, recv_type, r, DEFAULT_TAG, comm, MPI_STATUS_IGNORE);
            MPI_Type_free(&recv_type);
        }

        // Output 
        pgmwrite(output_file_name, master_content, length, width);
    }

    printf("End processing on rank %d.\n", rank);

    // End
    double end = MPI_Wtime();
    MPI_Barrier(comm);
    if (is_master(rank))
    {
        fprintf(stdout, "Execution time %.5f\n", end - start);
        fprintf(stdout, "End...\n");
    }

    // Free resources
    MPI_Type_free(&column_type);
    MPI_Type_free(&row_type);
    free(dims);
    free(periods);
    double_2d_array_deallocation(&masterbuf, &master_content);
    double_2d_array_deallocation(&img_new, &new_content);
    double_2d_array_deallocation(&img_old, &old_content);
    double_2d_array_deallocation(&img_edge, &edge_content);
    MPI_Finalize();
    return 0;
}

void check_arg(int argc, char **argv, char **input_file_name, char **output_file_name)
{
    fprintf(stdout, "Checking arguments.\n");
    if (argc < 3)
    {
        fprintf(stderr, "Input & Output file name required.\n");
        exit(1);
    }
    *input_file_name = argv[1];
    *output_file_name = argv[2];
    fprintf(stdout, "Argument checking finshed. Get input file name %s and output file name %s.\n", *input_file_name, *output_file_name);
}

double boundaryval(int i, int m)
{
    double val;

    val = 2.0 * ((double)(i - 1)) / ((double)(m - 1));
    if (i >= m / 2 + 1)
        val = 2.0 - val;

    return val;
}
