/* 
 * Coursework image processing parallel implementation with MPI.
 * 
 * Note that all the functions should be called in the order (init -> preprocessing -> processing -> postprocessing -> finalize)
 * to achieve the goal.
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#include "imgprocessing.h"
#include "defs.h"
#include "pgmio.h"
#include "arraytool.h"
#include "boundary.h"

/* Macros for MPI */
#define is_master(rank) (0 == rank) // function to check master

#define RANK_MASTER 0
/* MPI Comm tag */
#define DEFAULT_TAG 1
#define LEFT_SEND_TAG 2
#define RIGHT_SEND_TAG 3
#define LEFT_FROM_TAG 3
#define RIGHT_FROM_TAG 2
/* MPI Topology dimension */
#define NDIMS 2
#define DIM_X 0
#define DIM_Y 1
/* MPI Topology neighbor direction */
#define LEFT 0
#define RIGHT 1
#define UP 2
#define DOWN 3

MPI_Comm comm, cart_comm;    // MPI communicator
int size, rank;              // size of communicator & current rank
int neighbor[4];             // neighbors of current rank
int coords[NDIMS];           // coords of current rank in topology
int dims[NDIMS];             // dims of topolog
int length_buf, width_buf;   // buf for length & width of original image
int plength_buf, pwidth_buf; //buf for plength & pwidth
int x_step, y_step;          // step for partial file coords calculation
double start, end;           // execution time measurement

/* 
 * Initialization Task including MPI initialization
 * 
 * in : input_file_name , source image file name/path
 * out : length , length of the source image
 * out : width , width of the source image
 * out : plength , processing length
 * out : pwidth , processing width
 * 
 */
void init(char *input_file_name, int *length, int *width, int *plength, int *pwidth)
{
    /* MPI Init */
    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    /* Topology */
    int periods[NDIMS] = {1, 0};
    MPI_Dims_create(size, NDIMS, dims);
    MPI_Cart_create(comm, NDIMS, dims, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, NDIMS, coords);
    if (is_master(rank))
    {
        fprintf(stdout, "Mapping processes into %d x %d grids.\n", dims[DIM_X], dims[DIM_Y]);
    }
    fprintf(stdout, "Rank %d in coords (%d , %d)\n", rank, coords[DIM_X], coords[DIM_Y]);

    /* Neighbor */
    MPI_Cart_shift(cart_comm, 0, 1, &neighbor[LEFT], &neighbor[RIGHT]);
    MPI_Cart_shift(cart_comm, 1, 1, &neighbor[DOWN], &neighbor[UP]);

    /* Source image size */
    pgmsize(input_file_name, length, width);
    length_buf = *length;
    width_buf = *width;
    if (is_master(rank))
    {
        fprintf(stdout, "Original image size: %d x %d.\n", *length, *width);
        fprintf(stdout, "Max iteration:%d\n", MAX_LOOP);
    }

    /* Processing size */
    x_step = *length / dims[DIM_X];
    y_step = *width / dims[DIM_Y];
    *plength = (coords[DIM_X] == dims[DIM_X] - 1) ? x_step + *length % dims[DIM_X] : x_step;
    *pwidth = (coords[DIM_Y] == dims[DIM_Y] - 1) ? y_step + *width % dims[DIM_Y] : y_step;
    plength_buf = *plength;
    pwidth_buf = *pwidth;
    fprintf(stdout, "Rank %d plength %d, pwidth %d\n", rank, *plength, *pwidth);
}

/* 
 * Preprocessing Task, partial file preprocessing in MPI
 * 
 * in : input_file_name, source image file name/path
 * inout : edge, edge buf of source image edge
 * inout : old, old buf for processing
 * inout : new, new buf for processing
 * in : plength, processing length
 * in : pwidth, processing width
 */
void preprocessing(char *input_file_name, double **edge, double **old, double **new, int plength, int pwidth)
{
    /* Read Image */
    part_pgmread(input_file_name, &(edge[1][1]), plength, pwidth, coords[DIM_X] * x_step, coords[DIM_Y] * y_step, 2);

    /* Old buf init */
    for (int i = 0; i <= plength + 1; ++i)
    {
        for (int j = 0; j <= pwidth + 1; ++j)
        {
            old[i][j] = 255.0;
        }
    }

    /* Set fixed boundary conditions on the bottom and top sides */
    double val;
    if (MPI_PROC_NULL == neighbor[UP])
    {
        for (int i = 1; i <= plength; ++i)
        {
            val = boundaryval(coords[DIM_X] * x_step + i, length_buf);
            old[i][pwidth + 1] = (int)(255.0 * (1.0 - val));
        }
    }
    if (MPI_PROC_NULL == neighbor[DOWN])
    {
        for (int i = 1; i <= plength; ++i)
        {
            val = boundaryval(coords[DIM_X] * x_step + i, length_buf);
            old[i][0] = (int)(255.0 * val);
        }
    }
}

/* 
 * Processing Task , partial file processing & halo swap in MPI
 * 
 * in : edge, edge buf of source image edge
 * inout : old, old buf for processing
 * inout : new, new buf for processing
 * int : plength, processing length
 * int : pwidth, processing width
 */
void processing(double **edge, double **old, double **new, int plength, int pwidth)
{
    /* MPI Datatype declaration */
    MPI_Datatype column_type, row_type;
    MPI_Type_contiguous(pwidth + 2, MPI_DOUBLE, &column_type);
    MPI_Type_commit(&column_type);
    MPI_Type_vector(plength + 2, 1, pwidth + 2, MPI_DOUBLE, &row_type);
    MPI_Type_commit(&row_type);

    /* Iteration */
    MPI_Barrier(comm);
    start = MPI_Wtime();
    fprintf(stdout, "Rank %d start iteration...\n", rank);
    MPI_Request send_request[4];
    MPI_Status recv_status[4];
    for (int iter = 1; iter <= MAX_LOOP; ++iter)
    {
        /* Issend halo */
        MPI_Issend(&(old[1][0]), 1, column_type, neighbor[LEFT], LEFT_SEND_TAG, comm, &send_request[LEFT]);
        MPI_Issend(&(old[plength][0]), 1, column_type, neighbor[RIGHT], RIGHT_SEND_TAG, comm, &send_request[RIGHT]);
        MPI_Issend(&(old[0][1]), 1, row_type, neighbor[DOWN], DEFAULT_TAG, comm, &send_request[DOWN]);
        MPI_Issend(&(old[0][pwidth]), 1, row_type, neighbor[UP], DEFAULT_TAG, comm, &send_request[UP]);

        /* Inner part calculation */
        for (int i = 2; i <= plength - 1; ++i)
        {
            for (int j = 2; j <= pwidth - 1; ++j)
            {
                new[i][j] = 0.25 * (old[i - 1][j] + old[i + 1][j] + old[i][j - 1] + old[i][j + 1] - edge[i][j]);
            }
        }

        /* Recv halo */
        MPI_Recv(&(old[0][0]), 1, column_type, neighbor[LEFT], LEFT_FROM_TAG, comm, &recv_status[LEFT]);
        MPI_Recv(&(old[plength + 1][0]), 1, column_type, neighbor[RIGHT], RIGHT_FROM_TAG, comm, &recv_status[RIGHT]);
        MPI_Recv(&(old[0][0]), 1, row_type, neighbor[DOWN], DEFAULT_TAG, comm, &recv_status[DOWN]);
        MPI_Recv(&(old[0][pwidth + 1]), 1, row_type, neighbor[UP], DEFAULT_TAG, comm, &recv_status[UP]);

        /* Halo part calculation */
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

        /* Wait for Issend */
        MPI_Wait(&send_request[LEFT], MPI_STATUS_IGNORE);
        MPI_Wait(&send_request[RIGHT], MPI_STATUS_IGNORE);
        MPI_Wait(&send_request[DOWN], MPI_STATUS_IGNORE);
        MPI_Wait(&send_request[UP], MPI_STATUS_IGNORE);

        /* Copy new to old */
        double delta = 0, overall_delta = 0;
        double sum = 0, overall_sum = 0;
        for (int i = 1; i <= plength; ++i)
        {
            for (int j = 1; j <= pwidth; ++j)
            {
                delta = fabs(new[i][j] - old[i][j]) > delta ? fabs(new[i][j] - old[i][j]) : delta;
                sum += new[i][j];
                old[i][j] = new[i][j];
            }
        }
        MPI_Allreduce(&delta, &overall_delta, 1, MPI_DOUBLE, MPI_MAX, comm);
        if (overall_delta <= THRESHOLD)
        {
            if (is_master(rank))
            {
                printf("Finish at iteration %d, with delta : %.5f\n", iter, overall_delta);
            }
            break;
        }
        if (0 == iter % INTERVAL)
        {
            MPI_Allreduce(&sum, &overall_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
            if (is_master(rank))
            {
                printf("Iter %d : average pixels %.3f\n", iter, overall_sum / (length_buf * width_buf));
            }
        }
    }
    end = MPI_Wtime();
    MPI_Barrier(comm);

    /* Resource free */
    MPI_Type_free(&column_type);
    MPI_Type_free(&row_type);

    /* End processing */
    if (is_master(rank))
    {
        fprintf(stdout, "Iteration time %.5f\n", end - start);
    }
}

/* 
 * Postprocessing Task & Image output
 * 
* in : output_file_name , output file name/path
 * in : old, old buf for processing
 * in : plength, processing length
 * in : pwidth, processing width
 */
void postprocessing(char *output_file_name, double **old, int plength, int pwidth)
{
    /* Gather partial image */
    if (!is_master(rank))
    {
        /* Others send to master */
        MPI_Datatype send_type;
        MPI_Type_vector(plength, pwidth, pwidth + 2, MPI_DOUBLE, &send_type);
        MPI_Type_commit(&send_type);

        MPI_Ssend(&(old[1][1]), 1, send_type, RANK_MASTER, DEFAULT_TAG, comm);

        MPI_Type_free(&send_type);
    }
    if (is_master(rank))
    {
        /* Master recv partial image & output */
        /* Buf allocation */
        double **masterbuf;
        double *master_content;
        double_2d_array_allocation(&masterbuf, &master_content, length_buf, width_buf);

        /* Init partial image of master */
        for (int i = 0; i <= plength - 1; ++i)
        {
            for (int j = 0; j <= pwidth - 1; ++j)
            {
                masterbuf[i][j] = old[i + 1][j + 1];
            }
        }

        /* Recv */
        for (int r = 1; r <= size - 1; ++r)
        {
            int recv_coords[NDIMS]; // coords of source rank
            MPI_Cart_coords(cart_comm, r, NDIMS, recv_coords);
            int recv_plength = (recv_coords[DIM_X] == dims[DIM_X] - 1) ? x_step + length_buf % dims[DIM_X] : x_step;
            int recv_pwidth = (recv_coords[DIM_Y] == dims[DIM_Y] - 1) ? y_step + width_buf % dims[DIM_Y] : y_step;

            MPI_Datatype recv_type;
            MPI_Type_vector(recv_plength, recv_pwidth, width_buf, MPI_DOUBLE, &recv_type);
            MPI_Type_commit(&recv_type);

            MPI_Recv(&(masterbuf[recv_coords[DIM_X] * x_step][recv_coords[DIM_Y] * y_step]), 1, recv_type, r, DEFAULT_TAG, comm, MPI_STATUS_IGNORE);

            MPI_Type_free(&recv_type);
        }

        /* Output */
        pgmwrite(output_file_name, master_content, length_buf, width_buf);

        /* Memory deallocation */
        double_2d_array_deallocation(&masterbuf, &master_content);
    }
}

/* 
 * Finalization Task & resources deallocation
 * 
 * in : edge , edge buf
 * in : edge_content , continuous memory for edge
 * in : old , old buf
 * in : old_content , continuous memory for old
 * in : new , new buf
 * in : new_content , continuous memory for new
 */
void finalize(double **edge, double *edge_content, double **old, double *old_content, double **new, double *new_content)
{
    /* Memory deallocation */
    double_2d_array_deallocation(&new, &new_content);
    double_2d_array_deallocation(&old, &old_content);
    double_2d_array_deallocation(&edge, &edge_content);

    /* MPI Finalize */
    MPI_Finalize();
}