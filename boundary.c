#include "boundary.h"

double boundaryval(int i, int m)
{
    double val;

    val = 2.0 * ((double)(i - 1)) / ((double)(m - 1));
    if (i >= m / 2 + 1)
        val = 2.0 - val;

    return val;
}