/* Tools for c continuous array allocation & deallocation */
#include <stdio.h>
#include <stdlib.h>

#include "arraytool.h"

/* 
 *  Allocate memory for a continuous 2d double array 
 *  output: arr is the pointer to the 2d array 
 *  output: content is the pointer to the continuous memory
 *  input: length is the length of array
 *  input: width is the width of array
*/
void double_2d_array_allocation(double ***arr, double **content, int length, int width)
{
    if (length < 0 || width < 0)
    {
        fprintf(stderr, "Invalid array size.");
        exit(-1);
    }

    if (NULL == (*arr = (double **)malloc(length * sizeof(double *))))
    {
        fprintf(stderr, "Allocation Failed.");
        exit(-1);
    }
    if (NULL == (*content = (double *)malloc(length * width * sizeof(double))))
    {
        fprintf(stderr, "Allocation Failed.");
        exit(-1);
    }

    for (int i = 0; i != length; ++i)
    {
        (*arr)[i] = *content + width * i;
    }
}

/* 
 * Deallocate memory for a continuous 2d double array created by double_2d_array_allocation
 * input: arr is the pointer to the 2d array
 * input: content is the pointer to the continuous memory
 */
void double_2d_array_deallocation(double ***arr, double **content)
{
    if (NULL != content)
    {
        free(*content);
    }
    if (NULL != arr)
    {
        free(*arr);
    }
}