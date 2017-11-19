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
static void check_arg(int argc, char **argv);
/* Boundary value for top & buttom */
static double boundaryval(int i, int m);

char *input_file_name, *output_file_name; // file source & output target

int main(int argc, char **argv)
{
    check_arg(argc, argv);

    int size, rank;      // size of communicator & current rank
    int length, width;   // original file size
    int plength, pwidth; // partial file size of current process

    /* Init */
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    /* Topology */
    MPI_Comm cart_comm;
    int coords[NDIMS];
    int *dims = (int *)malloc(NDIMS * sizeof(int));
    int periods[NDIMS] = {1, 0};
    MPI_Dims_create(size, NDIMS, dims);
    MPI_Cart_create(comm, NDIMS, dims, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, NDIMS, coords);
    if (is_master(rank))
    {
        fprintf(stdout, "Mapping processes into %d x %d grids.\n", dims[DIM_X], dims[DIM_Y]);
    }
    fprintf(stdout, "Rank %d in coords (%d , %d)\n", rank, coords[DIM_X], coords[DIM_Y]);

    /* Get neighbor of current rank */
    int rank_up, rank_down, rank_left, rank_right;
    MPI_Cart_shift(cart_comm, 0, 1, &rank_left, &rank_right);
    MPI_Cart_shift(cart_comm, 1, 1, &rank_down, &rank_up);
    fprintf(stdout, "Neighbors of rank %d : Left %d Right %d Up %d Down %d\n", rank, rank_left, rank_right, rank_up, rank_down);

    /* Get image size */
    pgmsize(input_file_name, &length, &width);
    if (is_master(rank))
    {
        fprintf(stdout, "Original image size: %d x %d.\n", length, width);
    }

    /* Start */
    if (is_master(rank))
    {
        fprintf(stdout, "Start...\n");
    }
    MPI_Barrier(comm);
    double start = MPI_Wtime();

    /* Memory allocation */
    int x_step = length / dims[DIM_X];
    int y_step = width / dims[DIM_Y];
    plength = (coords[DIM_X] == dims[DIM_X] - 1) ? x_step + length % dims[DIM_X] : x_step;
    pwidth = (coords[DIM_Y] == dims[DIM_Y] - 1) ? y_step + width % dims[DIM_Y] : y_step;
    fprintf(stdout, "rank %d plength %d, pwidth %d\n", rank, plength, pwidth);

    double **new, **old, **edge;
    double *new_content, *old_content, *edge_content;
    double_2d_array_allocation(&new, &new_content, plength + 2, pwidth + 2);
    double_2d_array_allocation(&old, &old_content, plength + 2, pwidth + 2);
    double_2d_array_allocation(&edge, &edge_content, plength + 2, pwidth + 2);

    /* Read image */
    fprintf(stdout, "Rank %d reading from file %s...\n", rank, input_file_name);
    part_pgmread(input_file_name, &(edge[1][1]), plength, pwidth, coords[DIM_X] * x_step, coords[DIM_Y] * y_step, 2);

    /* Image reconstruction */
    for (int i = 0; i <= plength + 1; ++i)
    {
        for (int j = 0; j <= pwidth + 1; ++j)
        {
            old[i][j] = 255.0;
        }
    }

    /* Set fixed boundary conditions on the bottom and top sides */
    double val;
    if (MPI_PROC_NULL == rank_up)
    {
        for (int i = 1; i <= plength; ++i)
        {
            val = boundaryval(coords[DIM_X] * x_step + i, length);
            old[i][pwidth + 1] = (int)(255.0 * (1.0 - val));
        }
    }
    if (MPI_PROC_NULL == rank_down)
    {
        for (int i = 1; i <= plength; ++i)
        {
            val = boundaryval(coords[DIM_X] * x_step + i, length);
            old[i][0] = (int)(255.0 * val);
        }
    }

    /* Row & Column type declaration */
    MPI_Datatype column_type, row_type;
    MPI_Type_vector(pwidth + 2, 1, 1, MPI_DOUBLE, &column_type);
    MPI_Type_commit(&column_type);
    MPI_Type_vector(plength + 2, 1, pwidth + 2, MPI_DOUBLE, &row_type);
    MPI_Type_commit(&row_type);

    /* Start iteration */
    fprintf(stdout, "Start iteration on rank %d\n", rank);
    MPI_Request send_request[4];
    MPI_Status recv_status[4];
    for (int cnt = 0; cnt != MAX_LOOP; ++cnt)
    {
        MPI_Issend(&(old[1][0]), 1, column_type, rank_left, LEFT_SEND_TAG, comm, &send_request[0]);
        MPI_Issend(&(old[plength][0]), 1, column_type, rank_right, RIGHT_SEND_TAG, comm, &send_request[1]);
        MPI_Issend(&(old[0][1]), 1, row_type, rank_down, DEFAULT_TAG, comm, &send_request[2]);
        MPI_Issend(&(old[0][pwidth]), 1, row_type, rank_up, DEFAULT_TAG, comm, &send_request[3]);
        for (int i = 2; i <= plength - 1; ++i)
        {
            for (int j = 2; j <= pwidth - 1; ++j)
            {
                new[i][j] = 0.25 * (old[i - 1][j] + old[i + 1][j] + old[i][j - 1] + old[i][j + 1] - edge[i][j]);
            }
        }
        MPI_Recv(&(old[0][0]), 1, column_type, rank_left, LEFT_FROM_TAG, comm, &recv_status[0]);
        MPI_Recv(&(old[plength + 1][0]), 1, column_type, rank_right, RIGHT_FROM_TAG, comm, &recv_status[1]);
        MPI_Recv(&(old[0][0]), 1, row_type, rank_down, DEFAULT_TAG, comm, &recv_status[2]);
        MPI_Recv(&(old[0][pwidth + 1]), 1, row_type, rank_up, DEFAULT_TAG, comm, &recv_status[3]);
        for (int j = 1; j <= pwidth; ++j)
        {
            new[1][j] = 0.25 * (old[0][j] + old[2][j] + old[1][j - 1] + old[1][j + 1] - edge[1][j]);
            new[plength][j] = 0.25 * (old[plength - 1][j] + old[plength + 1][j] + old[plength][j - 1] + old[plength][j + 1] - edge[plength][j]);
        }
        for (int i = 1; i <= plength; ++i)
        {
            new[i][1] = 0.25 * (old[i - 1][1] + old[i + 1][1] + old[i][0] + old[i][2] - edge[i][1]);
            new[i][pwidth] = 0.25 * (old[i - 1][pwidth] + old[i + 1][pwidth] + old[i][pwidth - 1] + old[i][pwidth + 1] - edge[i][pwidth]);
        }
        MPI_Wait(&send_request[0], MPI_STATUS_IGNORE);
        MPI_Wait(&send_request[1], MPI_STATUS_IGNORE);
        MPI_Wait(&send_request[2], MPI_STATUS_IGNORE);
        MPI_Wait(&send_request[3], MPI_STATUS_IGNORE);

        double delta = 0, overall_delta = 0;
        for (int i = 1; i <= plength; ++i)
        {
            for (int j = 1; j <= pwidth; ++j)
            {
                delta = fabs(new[i][j] - old[i][j]) > delta ? fabs(new[i][j] - old[i][j]) : delta;
                old[i][j] = new[i][j];
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

    /* Gather partial image from ranks */
    MPI_Datatype send_type;
    MPI_Type_vector(plength, pwidth, pwidth + 2, MPI_DOUBLE, &send_type);
    MPI_Type_commit(&send_type);

    // Other ranks send to master
    if (!is_master(rank))
    {
        MPI_Ssend(&(old[1][1]), 1, send_type, RANK_MASTER, DEFAULT_TAG, comm);
    }
    MPI_Type_free(&send_type);
    // Master recv partial image & output
    if (is_master(rank))
    {
        double **masterbuf;
        double *master_content;
        double_2d_array_allocation(&masterbuf, &master_content, length, width);
        // Fill in image from rank master
        for (int i = 0; i <= plength - 1; ++i)
        {
            for (int j = 0; j <= pwidth - 1; ++j)
            {
                masterbuf[i][j] = old[i + 1][j + 1];
            }
        }

        // Recv image from other ranks
        for (int r = 1; r <= size - 1; ++r)
        {
            int recv_coords[NDIMS]; // coords of source rank
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

        // Memory free
        double_2d_array_deallocation(&masterbuf, &master_content);
    }
    printf("End processing on rank %d.\n", rank);

    /* End */
    double end = MPI_Wtime();
    MPI_Barrier(comm);
    if (is_master(rank))
    {
        fprintf(stdout, "Execution time %.5f\n", end - start);
        fprintf(stdout, "End...\n");
    }

    /* Free resources */
    MPI_Type_free(&column_type);
    MPI_Type_free(&row_type);
    free(dims);
    double_2d_array_deallocation(&new, &new_content);
    double_2d_array_deallocation(&old, &old_content);
    double_2d_array_deallocation(&edge, &edge_content);
    MPI_Finalize();
    return 0;
}

void check_arg(int argc, char **argv)
{
    fprintf(stdout, "Checking arguments.\n");
    if (argc < 3)
    {
        fprintf(stderr, "Input & Output file name required.\n");
        exit(1);
    }
    input_file_name = argv[1];
    output_file_name = argv[2];
    fprintf(stdout, "Argument checking finshed. Get input file name %s and output file name %s.\n", input_file_name, output_file_name);
}

double boundaryval(int i, int m)
{
    double val;

    val = 2.0 * ((double)(i - 1)) / ((double)(m - 1));
    if (i >= m / 2 + 1)
        val = 2.0 - val;

    return val;
}
