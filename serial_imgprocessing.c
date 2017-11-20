#include <stdio.h>
#include <stdlib.h>

#include "imgprocessing.h"
#include "defs.h"
#include "boundary.h"
#include "arraytool.h"
#include "pgmio.h"

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
    /* Init */
    fprintf(stdout, "Init...\n");

    /* Read file */
    pgmsize(input_file_name, length, width);
    *plength = *length;
    *pwidth = *width;
    fprintf(stdout, "Original file size %d x %d\n", *length, *width);
    fprintf(stdout, "Max iteration:%d\n", MAX_LOOP);

    /* End init */
    fprintf(stdout, "Init finished...\n");
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
    /* Preprocessing */
    fprintf(stdout, "Preprocessing...\n");

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
    /* End preprocessing */
    fprintf(stdout, "End preprocessing...\n");
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
    /* Processing */
    fprintf(stdout, "Processing...\n");

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

        for (int i = 1; i < plength + 1; i++)
        {
            for (int j = 1; j < pwidth + 1; j++)
            {
                old[i][j] = new[i][j];
            }
        }
    }

    /* End processing */
    fprintf(stdout, "Processing finished...\n");
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
    /* Postprocessing */
    fprintf(stdout, "Postprocessing...\n");

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
    /* End postprocessing */
    fprintf(stdout, "End postprocessing...\n");
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
    /* Finalize */
    fprintf(stdout, "Finalize...\n");

    /* Memory deallocation */
    double_2d_array_deallocation(&new, &new_content);
    double_2d_array_deallocation(&old, &old_content);
    double_2d_array_deallocation(&edge, &edge_content);

    /* End Finalization */
    fprintf(stdout, "End finalization...\n");
}