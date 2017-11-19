/* Tools for c continuous array allocation & deallocation */ 
#ifndef __ARRAY_TOOL_H
#define __ARRAY_TOOL_H

#include <stdlib.h>

void double_2d_array_allocation(double ***arr, double **content, int length, int width);
void double_2d_array_deallocation(double ***arr, double **content);

#endif