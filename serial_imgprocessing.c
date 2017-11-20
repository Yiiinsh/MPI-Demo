#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "imgprocessing.h"
#include "defs.h"
#include "boundary.h"
#include "arraytool.h"
#include "pgmio.h"

clock_t start, end; // time measurement

/* 
 * Initialization Task
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
    /* Read file */
    pgmsize(input_file_name, length, width);
    *plength = *length;
    *pwidth = *width;
    fprintf(stdout, "Original file size %d x %d\n", *length, *width);
    fprintf(stdout, "Max iteration:%d\n", MAX_LOOP);
}

/* 
 * Preprocessing Task
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
    /* Read image */
    double **buf;
    double *buf_content;
    double_2d_array_allocation(&buf, &buf_content, plength, pwidth);
    pgmread(input_file_name, buf_content, plength, pwidth);

    /* Copy image to edge */
    for (int i = 1; i < plength + 1; ++i)
    {
        for (int j = 1; j < pwidth + 1; ++j)
        {
            edge[i][j] = buf[i - 1][j - 1];
        }
    }

    /* Old buf init */
    for (int i = 0; i < plength + 2; ++i)
    {
        for (int j = 0; j < pwidth + 2; ++j)
        {
            old[i][j] = 255.0;
        }
    }

    /* Set fixed boundary conditions on the bottom and top sides */
    for (int i = 1; i < plength + 1; ++i)
    {
        double val = boundaryval(i, plength);
        old[i][0] = (int)(255.0 * val);
        old[i][pwidth + 1] = (int)(255.0 * (1.0 - val));
    }

    double_2d_array_deallocation(&buf, &buf_content);
}

/* 
 * Processing Task
 * 
 * in : edge, edge buf of source image edge
 * inout : old, old buf for processing
 * inout : new, new buf for processing
 * int : plength, processing length
 * int : pwidth, processing width
 */
void processing(double **edge, double **old, double **new, int plength, int pwidth)
{
    start = clock();
    for (int iter = 1; iter <= MAX_LOOP; ++iter)
    {
        /* Implement periodic boundary conditions on left and right sides */
        for (int j = 1; j < pwidth + 1; j++)
        {
            old[0][j] = old[plength][j];
            old[plength + 1][j] = old[1][j];
        }

        for (int i = 1; i < plength + 1; i++)
        {
            for (int j = 1; j < pwidth + 1; j++)
            {
                new[i][j] = 0.25 * (old[i - 1][j] + old[i + 1][j] + old[i][j - 1] + old[i][j + 1] - edge[i][j]);
            }
        }

        double delta = 0;
        double sum = 0;
        for (int i = 1; i < plength + 1; i++)
        {
            for (int j = 1; j < pwidth + 1; j++)
            {
                delta = fabs(new[i][j] - old[i][j]) > delta ? fabs(new[i][j] - old[i][j]) : delta;
                sum += new[i][j];
                old[i][j] = new[i][j];
            }
        }
        if (delta <= THRESHOLD)
        {
            printf("Finish at iteration %d, with delta : %.5f\n", iter, delta);
            break;
        }
        if (0 == iter % INTERVAL)
        {
            printf("Iter %d : average pixels %.3f\n", iter, sum / (plength * pwidth));
        }
    }
    end = clock();
    fprintf(stdout, "Iteration time %.5f\n", (double)(end - start) / CLOCKS_PER_SEC);
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
    double **buf;
    double *buf_content;
    double_2d_array_allocation(&buf, &buf_content, plength, pwidth);

    for (int i = 1; i < plength + 1; ++i)
    {
        for (int j = 1; j < pwidth + 1; ++j)
        {
            buf[i - 1][j - 1] = old[i][j];
        }
    }

    pgmwrite(output_file_name, buf_content, plength, pwidth);

    double_2d_array_deallocation(&buf, &buf_content);
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
}