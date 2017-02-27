/* sysnewton.c - Solve system of non-linear equations with newton's method
 * $Id: sysnewton.c,v 1.1 2000/10/20 20:17:20 antoine Exp $
 * Copyright (C) 2000
 *    Antoine Lefebvre <antoine.lefebvre@polymtl.ca>
 *
 * Licensed under the GPL
 */
   

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "num.h"


double norme(double *x, int n);

/* Jac : Pointer to a matrix of pointer to function (Jacobian matrix)
 * R   : Pointer to a vector of pointer to function (Residue vector)
 * x   : Initial estimate of the solution (will be overwrite by the answer)
 * nvar: Number of variable in the system
 * nmax: Maximal number of iterations
 * eps : Precision on the answer
 */
int NUM_sysnewton(func_t *Jac, func_t *R, double *x, int nvar,
                  int nmax, double eps)
{
  int i, j, l;
  
  double *r;
  double *dx; /* solution of the system */
  double *matrix; /* the matrix to be solve */

  /* Allow space for the left hand side */
  matrix = (double *) malloc (sizeof(double) * nvar * (nvar + 1));

  r = (double *) malloc (sizeof(double) * nvar);
  dx = (double *) malloc(sizeof(double) * nvar);

  l = 0;
  do
  {
    /* set up the matrix */
    for (i = 0; i < nvar; i++) /* line */
    {
      for (j = 0; j < nvar; j++) /* column */
      {
        matrix[i + nvar*j] = Jac[i + nvar*j](x);
      }
      r[i] = matrix[i + nvar*nvar] = -R[i](x);
    }

    /* if the matrix is singular */
    if (NUM_lu(matrix, dx, nvar))
      return NO_CONVERGENCE;
    
    for (i = 0; i < nvar; i++)
    {
      x[i] = x[i] + dx[i];
    }
    
    if ((norme(dx, nvar)/norme(x, nvar) < eps) &&
        (norme(r, nvar) <= eps))
    {
      /* the solution converged */
      return 0;
    }

    l++;
  } while (l < nmax); 

  return NO_CONVERGENCE;
  
}


double norme(double *x, int n)
{
  int i;
  double a = 0.0;
  for (i = 0; i < n; i++)
  {
    a += pow(x[i], 2);
  }
  return sqrt(a);
}
