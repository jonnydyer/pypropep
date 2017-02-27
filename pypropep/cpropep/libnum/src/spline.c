/* spline.c  -  Cubic spline interpolation
 * $Id: spline.c,v 1.1 2001/06/10 21:06:00 antoine Exp $
 * Copyright (C) 2000
 *    Antoine Lefebvre <antoine.lefebvre@polymtl.ca>
 *
 *
 * Licensed under the GPLv2
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "num.h"

#define OUT_OF_RANGE -1

/* Computation of natural spline.
   
 */

int create_spline(double *data, int n_point, double *spline)
{
  int i;
  int col;
  double hi0, hi1, hi2;
  double fi0, fi1, fi2;
  double *matrix;
  
  col = 2;
  
  matrix = (double *) calloc (n_point*(n_point+1), sizeof(double));
 
  /* create the linear system */
  matrix[0 + 0] = 1;
  
  for (i = 1; i < n_point - 1; i++)
  {
    hi0 = data[0 + (i+0)*col] - data[0 + (i-1)*col];
    hi1 = data[0 + (i+1)*col] - data[0 + (i-0)*col];
    hi2 = data[0 + (i+1)*col] - data[0 + (i-1)*col];
    
    fi0 = data[1 + (i-1)*col]; /* f(x_(i-1)) */
    fi1 = data[1 + (i-0)*col]; /* f(x_i)     */
    fi2 = data[1 + (i+1)*col]; /* f(x_(i+1)) */
   
    matrix[i + (i-1)*n_point] = hi0/(hi0+hi1);
    matrix[i + (i-0)*n_point] = 2;
    matrix[i + (i+1)*n_point] = hi1/(hi0+hi1);

    matrix[i + n_point*n_point] = 6*( (fi2-fi1)/hi1 - (fi1-fi0)/hi0)/hi2;
  }

  matrix[n_point - 1 + (n_point - 1)*n_point] = 1;

  NUM_print_matrix(matrix, n_point);  
  NUM_lu(matrix, spline, n_point);
  NUM_print_vec(spline, n_point);
  
  free(matrix);
  return 0;
}



int eval_spline(double *data, double *spline, int n_point, double x, double *y)
{

  int i;

  int col = 2;

  double hi0;
  double xi0;
  double xi1;
  double fi0;
  double fi1;
  double d2f0;
  double d2f1;
  
  /* we must find the interval in which x is located */

  i = 0;

  if ((x < data[0]) || (x > data[0 + (n_point - 1)*col]))
  {
    printf("Out of range: %f\n", x);
    return OUT_OF_RANGE;
  }
  
  for (i = 0; x > data[0 + i*col]; i++);

  xi0  = data[0 + (i-1)*col];
  xi1  = data[0 + i*col];
  hi0  = xi1 - xi0;
  fi0  = data[1 + (i-1)*col];
  fi1  = data[1 + i*col];
  d2f0 = spline[i-1];
  d2f1 = spline[i];

  /* evaluation of the spline at the specified point */
  *y = -d2f0*pow(x-xi1, 3)/(6*hi0) + d2f1*pow(x-xi0, 3)/(6*hi0)
    -(fi0/hi0 - hi0*d2f0/6)*(x-xi1) + (fi1/hi0 - hi0*d2f1/6)*(x-xi0);
  
  return 0;
}


