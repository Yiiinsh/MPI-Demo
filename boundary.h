#ifndef __BOUNDARY_H
#define __BOUNDARY_H

/* 
 * Boundary condition for top & buttom
 * 
 * in : i, current i index
 * in : m, length of i
 * out : sawtooth value of input boundary condition
 */
double boundaryval(int i, int m);

#endif