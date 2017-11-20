#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#include "defs.h"
#include "pgmio.h"
#include "arraytool.h"
#include "boundary.h"
#include "imgprocessing.h"

/* Check arguments for user input */
static void check_arg(int argc, char **argv);

char *input_file_name, *output_file_name; // file source & output target

int main(int argc, char **argv)
{
    int length, width;                                // original image size
    int plength, pwidth;                              // processing image size
    double **edge, **old, **new;                      // image buf
    double *edge_content, *old_content, *new_content; // continuous memory for image buf

    /* Argument checking */
    check_arg(argc, argv);

    /* Init */
    init(input_file_name, &length, &width, &plength, &pwidth);

    /* Memory allocation */
    double_2d_array_allocation(&new, &new_content, plength + 2, pwidth + 2);
    double_2d_array_allocation(&old, &old_content, plength + 2, pwidth + 2);
    double_2d_array_allocation(&edge, &edge_content, plength + 2, pwidth + 2);

    /* Preprocessing */
    preprocessing(input_file_name, edge, old, new, plength, pwidth);

    /* Processing */
    processing(edge, old, new, plength, pwidth);

    /* Postprocessing */
    postprocessing(output_file_name, old, plength, pwidth);

    /* Finalize */
    finalize(edge, edge_content, old, old_content, new, new_content);

    return 0;
}

void check_arg(int argc, char **argv)
{
    if (argc < 3)
    {
        fprintf(stderr, "Input & Output file name required.\n");
        exit(1);
    }
    input_file_name = argv[1];
    output_file_name = argv[2];
}