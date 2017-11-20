/* 
 * Coursework image processing header file.
 * 
 * Note that all the functions should be called in the order (init -> preprocessing -> processing -> postprocessing -> finalize)
 * to achieve the goal.
 */
#ifndef __IMGPROCESSING_H
#define __IMGPROCESSING_H

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
void init(char *input_file_name, int *length, int *width, int *plength, int *pwidth);

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
void preprocessing(char *input_file_name, double **edge, double **old, double **new, int plength, int pwidth);

/* 
 * Processing Task
 * 
 * in : edge, edge buf of source image edge
 * inout : old, old buf for processing
 * inout : new, new buf for processing
 * int : plength, processing length
 * int : pwidth, processing width
 */
void processing(double **edge, double **old, double **new, int plength, int pwidth);

/* 
 * Postprocessing Task & Image output
 * 
 * in : output_file_name , output file name/path
 * in : old, old buf for processing
 * in : plength, processing length
 * in : pwidth, processing width
 */
void postprocessing(char *output_file_name, double **old, int plength, int pwidth);

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
void finalize(double **edge, double *edge_content, double **old, double *old_content, double **new, double *new_content);

#endif