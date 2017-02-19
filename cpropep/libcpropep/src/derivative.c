/* derivative.c  -  Fill the mattrix to compute thermochemical derivative
                    relative to logarithm of pressure and temperature */
/* $Id: derivative.c,v 1.2 2001/02/22 19:49:28 antoine Exp $ */
/* Copyright (C) 2000                                                  */
/*    Antoine Lefebvre <antoine.lefebvre@polymtl.ca>                   */
/*    Mark Pinese <pinese@cyberwizards.com.au>                         */
/*                                                                     */
/* Licensed under the GPLv2                                            */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "num.h"

#include "derivative.h"
#include "equilibrium.h"
#include "performance.h"
#include "thermo.h"

#include "print.h"
#include "compat.h"
#include "return.h"

int fill_temperature_derivative_matrix(double *matrix, equilibrium_t *e);
int fill_pressure_derivative_matrix(double *matrix, equilibrium_t *e);

/* Compute the specific_heat of the mixture using thermodynamics
   derivative with respect to logarithm of temperature */
double mixture_specific_heat(equilibrium_t *e, double *sol)
{
  short i, j;
  double cp, tmp;

  product_t       *p  = &(e->product);
  equilib_prop_t  *pr = &(e->properties);
  
  cp = 0.0;
  /* Compute Cp/R */
  for (i = 0; i < p->n_element; i++)
  {
    tmp = 0.0;
    for (j = 0; j < p->n[GAS]; j++)
      tmp += p->A[i][j] * p->coef[GAS][j] *
        enthalpy_0(p->species[GAS][j], pr->T);
    
    cp += tmp * sol[i];
    
  }
  
  for (i = 0; i < p->n[CONDENSED]; i++)
  {
    cp += enthalpy_0(p->species[CONDENSED][i], pr->T)
      * sol[i + p->n_element];
  }
  
  tmp = 0.0;
  for (i = 0; i < p->n[GAS]; i++)
  {
    tmp += p->coef[GAS][i] * enthalpy_0(p->species[GAS][i], pr->T);
  }
  cp += tmp * sol[p->n_element + p->n[CONDENSED]];
  
  cp += mixture_specific_heat_0(e, pr->T);
  
  for (i = 0; i < p->n[GAS]; i++)
  {
    cp += p->coef[GAS][i] * pow(enthalpy_0(p->species[GAS][i], pr->T), 2);
  }

  return cp;
}

int derivative(equilibrium_t *e)
{
  short size;
  double *matrix;
  double *sol;

  product_t      *p    = &(e->product);
  equilib_prop_t *prop = &(e->properties);
  
  
  /* the size of the coefficient matrix */
  size = p->n_element + p->n[CONDENSED] + 1;

  matrix = (double *) malloc (size * (size + 1) * sizeof(double));
  
  /* allocate the memory for the solution vector */
  sol = (double *) calloc (size, sizeof(double));

  fill_temperature_derivative_matrix(matrix, e);
  
  if (NUM_lu(matrix, sol, size) == -1)
  {
    fprintf(outputfile, "The matrix is singular.\n");
  }
  else
  {
    if (global_verbose > 2)
    {
      fprintf(outputfile, "Temperature derivative results.\n");
      NUM_print_vec(sol, size);
    }
    
    prop->Cp   = mixture_specific_heat(e, sol)*R;
    prop->dV_T = 1 + sol[e->product.n_element + e->product.n[CONDENSED]];  
  }

  fill_pressure_derivative_matrix(matrix, e);

  if (NUM_lu(matrix, sol, size) == -1)
  {
    fprintf(outputfile, "The matrix is singular.\n");
  }
  else
  {
    if (global_verbose > 2)
    {
      fprintf(outputfile, "Pressure derivative results.\n");
      NUM_print_vec(sol, size);
    }
    prop->dV_P = sol[e->product.n_element + e->product.n[CONDENSED]] - 1;
    
  }

  prop->Cv    = prop->Cp + e->itn.n * R * pow(prop->dV_T, 2)/prop->dV_P;
  prop->Isex  = -(prop->Cp / prop->Cv) / prop->dV_P;
  prop->Vson  = sqrt(1000 * e->itn.n * R * e->properties.T * prop->Isex);
  
  free(matrix);
  free(sol);
  return 0;
}


/* Fill the matrix with the coefficient for evaluating derivatives with
   respect to logarithm of temperature at constant pressure */
int fill_temperature_derivative_matrix(double *matrix, equilibrium_t *e)
{
  
  short j, k, size;
  double tmp;

  short idx_cond, idx_n, idx_T;

  product_t       *p  = &(e->product);
  equilib_prop_t  *pr = &(e->properties);

  idx_cond  = p->n_element;
  idx_n     = p->n_element + p->n[CONDENSED];
  idx_T     = p->n_element + p->n[CONDENSED] + 1;
  
  /* the size of the coefficient matrix */
  size      = p->n_element + p->n[CONDENSED] + 1;
  
  /* fill the common part */
  fill_matrix(matrix, e, TP);
  
  /* del ln(n)/ del ln(T) */
  matrix[idx_n + size * idx_n] = 0.0;
  
  /* right side */
  for (j = 0; j < p->n_element; j++)
  {
    tmp = 0.0;
    for (k = 0; k < p->n[GAS]; k++)
      tmp -= p->A[j][k] * p->coef[GAS][k] * enthalpy_0(p->species[GAS][k],
                                                       pr->T);
    matrix[j + size * idx_T] = tmp;
  }

  for (j = 0; j < p->n[CONDENSED]; j++) /* row */
    matrix[j + idx_cond + size * idx_T] =
      -enthalpy_0( p->species[CONDENSED][j], pr->T);
  
  tmp = 0.0;
  for (k = 0; k < p->n[GAS]; k++)
    tmp -= p->coef[GAS][k]*enthalpy_0(p->species[GAS][k], pr->T); 

  matrix[idx_n + size * idx_T] = tmp;
  
  return 0;
}

/* Fill the matrix with the coefficient for evaluating derivatives with
   respect to logarithm of pressure at constant temperature */
int fill_pressure_derivative_matrix(double *matrix, equilibrium_t *e)
{
  
  short j, k, size;
  double tmp;

  short idx_cond, idx_n, idx_T;

  product_t       *p  = &(e->product);

  idx_cond  = p->n_element;
  idx_n     = p->n_element + p->n[CONDENSED];
  idx_T     = p->n_element + p->n[CONDENSED] + 1;
  
  /* the size of the coefficient matrix */
  size = p->n_element + p->n[CONDENSED] + 1;
  
  /* fill the common part */
  fill_matrix(matrix, e, TP);
  
  /* del ln(n)/ del ln(T) */
  matrix[idx_n + size * idx_n] = 0.0;
  
  /* right side */
  for (j = 0; j < p->n_element; j++)
  {
    tmp = 0.0;
    for (k = 0; k < p->n[GAS]; k++)
      tmp += p->A[j][k] * p->coef[GAS][k];

    matrix[j + size * idx_T] = tmp;
  }

  for (j = 0; j < p->n[CONDENSED]; j++) /* row */
    matrix[j + idx_cond + size * idx_T] = 0.0;
  
  tmp = 0.0;
  for (k = 0; k < p->n[GAS]; k++)
    tmp += p->coef[GAS][k]; 

  matrix[idx_n + size * idx_T] = tmp;
  
  return 0;
}









