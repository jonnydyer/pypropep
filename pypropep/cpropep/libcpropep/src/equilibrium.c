/* equilibrium.c  -  Responsible of the chemical equilibrium          */
/* $Id: equilibrium.c,v 1.2 2001/02/22 19:49:28 antoine Exp $ */
/* Copyright (C) 2000                                                  */
/*    Antoine Lefebvre <antoine.lefebvre@polymtl.ca>                   */
/*    Mark Pinese <pinese@cyberwizards.com.au>                         */
/*                                                                     */
/* Licensed under the GPLv2                                            */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "num.h" /* matrix solution */

#include "print.h"
#include "equilibrium.h"

#include "conversion.h"
#include "compat.h"
#include "return.h"

#include "thermo.h" /* thermodynamics function */

/* Initial temperature estimate for problem with not-fixed temperature */
#define ESTIMATED_T 3800

#define CONC_TOL       1.0e-8
#define LOG_CONC_TOL -18.420681
#define CONV_TOL       0.5e-5

#define ITERATION_MAX 100


/* 1 for verbose, 0 for non-verbose */
int global_verbose = 0;

double product_molar_mass(equilibrium_t *e)
{
  return (1/e->itn.n);
}

int list_element(equilibrium_t *e)
{
  int n = 0;
  int t = 0;
  int i, j, k;

  composition_t *prop = &(e->propellant);
  product_t     *prod = &(e->product);
  
  /* reset the lement vector to -1 */
  reset_element_list(e);

  for (i = 0; i < e->propellant.ncomp; i++)
  {
    /* maximum of 6 different atoms in the composition */
    for (j = 0; j < 6; j++)
    {	       
      if (!( (propellant_list + prop->molecule[i])->coef[j] == 0))
      {
        /* get the element */
        t = (propellant_list + prop->molecule[i])->elem[j];
        
        for (k = 0; k <= n; k++)
        {
          /* verify if the element was not already in the list */
          if (prod->element[k] == t)
            break;
          /* if we have check each element, add it to the list */
          if (k == n)
          {
            if (n == MAX_ELEMENT)
            {
              fprintf(errorfile, "Maximum of %d elements. Abort.\n",
                      MAX_ELEMENT);
            }
            prod->element[n] = t;
            n++;
            break;
          }
        }
      }
    }
  }
  prod->n_element      = n;
  prod->element_listed = 1;
  
  return n;
}

/************************************************************
FUNCTION: This function search in thermo_list for all molecule
          that could be form with one or more of the element
          in element_list. The function fill product_list with
          the corresponding number of these molecule.

PARAMETER: e is a pointer to an equilibrium_t structure

COMMENTS: The return value is the number of elements found

DATE: February 6, 2000

AUTHOR: Antoine Lefebvre
**************************************************************/
int list_product(equilibrium_t *e)
{
  int i, j, k;

  int n = 0;   /* global counter (number of species found) */
  int st;      /* temporary variable to hold the state of one specie */
  int ok = 1;

  product_t    *prod = &(e->product);
  
  /* reset the product to zero */
  prod->n[GAS]       = 0;
  prod->n[CONDENSED] = 0;
  
  for (j = 0; j < num_thermo; j++)
  {
    /* for each of the five possible element of a species */
    for (k = 0; k < 5; k++)
    {
      if (!((thermo_list + j)->coef[k] == 0))
      {
        for (i = 0; i < prod->n_element; i++)
        {
          if (prod->element[i] == (thermo_list + j)->elem[k])
            break;
          else if (i == (prod->n_element - 1) )
            ok = 0;
        }    
        if (!ok)
          break;
      }
    }
    if (ok) /* add to the list */
    {
      st = (thermo_list + j)->state;

      prod->species[st][ prod->n[st] ] = j;
      prod->n[st]++;
      n++;
      
      if ((prod->n[GAS] > MAX_PRODUCT) || (prod->n[CONDENSED] > MAX_PRODUCT))
      {
        fprintf(errorfile,
                "Error: Maximum of %d differents product reach.\n",
                MAX_PRODUCT);
        fprintf(errorfile, "       Change MAX_PRODUCT and recompile!\n");
        return ERR_TOO_MUCH_PRODUCT;
      }
       
    }
    ok = 1;
  }

  prod->n_condensed = prod->n[CONDENSED];

  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    move it to the equilibrium function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  
  /* initialize tho mol number to 0.1mol/(nb of gazeous species) */
  e->itn.n    = e->itn.sumn = 0.1;
  e->itn.ln_n = log(e->itn.n);
  
  for (i = 0; i < e->product.n[GAS]; i++)
  {
    e->product.coef[GAS][i] = 0.1 / e->product.n[GAS];
    e->itn.ln_nj[i]       = log( e->product.coef[GAS][i] );
  }
  
  /* initialize condensed to zero */
  for (i = 0; i < e->product.n[CONDENSED]; i++)
    e->product.coef[CONDENSED][i] = 0;

  e->product.product_listed = 1;
  
  return n;


}

/* Initialisation of the product_t structure */
int initialize_product(product_t *p)
{
  int i, j;
  
  for (i = 0; i < STATE_LAST; i++)
    p->n[i] = 0;

  /* initialize the list to -1 */
  for (j = 0; j < STATE_LAST; j++)
  {
    for (i = 0; i < MAX_PRODUCT; i++)
      p->species[j][i] = -1;
  }

  p->n_condensed = 0;
  p->product_listed = 0;
  return 0;
}


int initialize_equilibrium(equilibrium_t *e)
{ 

  /* the composition have not been set */
  e->propellant.ncomp = 0;
  
  e->product.isequil        = false;
  e->product.element_listed = 0; /* the element haven't been listed */
  
  /* initialize the product */
  return initialize_product(&(e->product));

}


int copy_equilibrium(equilibrium_t *dest, equilibrium_t *src)
{
  memcpy(dest, src, sizeof(equilibrium_t));
  return 0;
}

int reset_element_list(equilibrium_t *e)
{
  int i;
  for (i = 0; i < MAX_ELEMENT; i++)
    e->product.element[i] = -1;
  return 0;
}


/* Probably not need */
int reset_equilibrium(equilibrium_t *e)
{
  int i;

  e->equilibrium_ok = false;
  e->itn.n = 0.1;
  /* initialize tho mol number to 0.1mol/(nb of gazeous species) */
  for (i = 0; i < e->product.n[GAS]; i++)
    e->product.coef[GAS][i] = 0.1/e->product.n[GAS];
  
  /* initialize condensed to zero */
  for (i = 0; i < e->product.n[CONDENSED]; i++)
    e->product.coef[CONDENSED][i] = 0;

  return 0;
}

int set_state(equilibrium_t *e, double T, double P)
{
  e->properties.T = T;
  e->properties.P = P;
  return 0;
}

int add_in_propellant(equilibrium_t *e, int sp, double mol)
{
  composition_t *c = &(e->propellant);
  c->molecule[ c->ncomp ] = sp;
  c->coef[ c->ncomp ]     = mol;
  c->ncomp++;
  return 0;
}


int product_element_coef(int element, int molecule)
{
  int i;
  for (i = 0; i < 5; i++)
  {
    if ((thermo_list + molecule)->elem[i] == element)
      return (thermo_list + molecule)->coef[i];
  }
  return 0;
}

int propellant_element_coef(int element, int molecule)
{
  int i;
  for (i = 0; i < 6; i++)
  {
    if ((propellant_list + molecule)->elem[i] == element)
      return (propellant_list + molecule)->coef[i];
  }
  return 0;
}

int compute_thermo_properties(equilibrium_t *e)
{
  equilib_prop_t  *pr = &(e->properties);
  /* Compute equilibrium properties */
  pr->H = product_enthalpy(e) * R * pr->T;
  pr->U = (product_enthalpy(e) - e->itn.n) * R * pr->T;
  pr->G = (product_enthalpy(e) - product_entropy(e)) * R * pr->T;
  pr->S = product_entropy(e) * R;  
  pr->M = product_molar_mass(e);
  e->properties_ok = true;
  return 0;
}

/* Compute an initial estimate of the product composition using
   a method develop by G. Eriksson */

/* not use for the moment */
/*
int initial_estimate(equilibrium_t *e)
{
  int i, j, mol;
  int components = e->product.n[GAS] + e->product.n[CONDENSED];
  double energy[components];

  for (i = 0; i < e->product.n[CONDENSED]; i++)
  {
    mol = 0;
    for (j = 0; j < 5; j++)
      mol += (thermo_list + e->product.species[CONDENSED][i])->coef[j];

    energy[i] = (enthalpy_0(e->product.species[CONDENSED][i], e->properties.T) -
                 entropy_0(e->product.species[CONDENSED][i], e->properties.T))*R*e->properties.T/mol;
    printf("%s \t %f \t %i\n",
           (thermo_list + e->product.species[CONDENSED][i])->name,
           energy[i], mol);
  }
  for (i = 0; i < e->product.n[GAS]; i++)
  {
    mol = 0;
    for (j = 0; j < 5; j++)
      mol += (thermo_list + e->product.species[GAS][i])->coef[j];
    
    energy[i + e->product.n[CONDENSED]] = (enthalpy_0(e->product.species[GAS][i],
                                               e->properties.T) -
                                      entropy_0(e->product.species[GAS][i],
                                              e->properties.T))*R*e->properties.T/mol;
    
    printf("%s \t %f \t %i\n",
           (thermo_list + e->product.species[GAS][i])->name,
           energy[i + e->product.n[CONDENSED]], mol);
  }
  
  return 0;
}
*/

int fill_equilibrium_matrix(double *matrix, equilibrium_t *e, problem_t P)
{

  short i, j, k;
  double tmp, mol;

  /* position of the right side dependeing on the type of problem */
  short roff = 2, size;
  
  double Mu[STATE_LAST][MAX_PRODUCT]; /* gibbs free energy for gases */
  double Ho[STATE_LAST][MAX_PRODUCT]; /* enthalpy in the standard state */

  /* The matrix is separated in five parts
     1- lagrangian multiplier (start at zero)
     2- delta(nj) for condensed (start at n_element)
     3- delta(ln n) (start at n_element + n[CONDENSED])
     4- delta(ln T) (start at n_element + n[CONDENSED] + 1)
     5- right side (start at n_element + n[CONDENSED] + roff)

     we defined one index for 2, 3 and 4
     the first start to zero and the last is the matrix size
  */
  short idx_cond, idx_n, idx_T;

  product_t       *p  = &(e->product);
  equilib_prop_t  *pr = &(e->properties);
  iteration_var_t *it = &(e->itn);
  
  if (P == TP)
    roff = 1;

  idx_cond  = p->n_element;
  idx_n     = p->n_element + p->n[CONDENSED];
  idx_T     = p->n_element + p->n[CONDENSED] + 1;

  size      = p->n_element + p->n[CONDENSED] + roff;
    
  mol = it->sumn;

  for (k = 0; k < p->n[GAS]; k++)
  {
    Mu[GAS][k] = gibbs(p->species[GAS][k], GAS, it->ln_nj[k] - it->ln_n,
                       pr->T, pr->P);
    Ho[GAS][k] = enthalpy_0( p->species[GAS][k], pr->T);
  }

  for (k = 0; k < p->n[CONDENSED]; k++)
  {
    Mu[CONDENSED][k] = gibbs(p->species[CONDENSED][k], CONDENSED, 0,
                             pr->T, pr->P);
    Ho[CONDENSED][k] = enthalpy_0(p->species[CONDENSED][k], pr->T);
  }
  
  /* fill the common part of the matrix */
  fill_matrix(matrix, e, P);
  
  /* delta ln(T) (for SP and HP only) */
  if (P != TP)
  {
    for (j = 0; j < p->n_element; j++)
    {
      tmp = 0.0;
      for (k = 0; k < p->n[GAS]; k++) {
        tmp += p->A[j][k] * p->coef[GAS][k] * Ho[GAS][k];
      }

      matrix[j + size * idx_T] = tmp;
    }
  }

  /* right side */
  for (j = 0; j < p->n_element; j++)
  {
    tmp = 0.0;
    
    for (k = 0; k < p->n[GAS]; k++)
      tmp += p->A[j][k] * p->coef[GAS][k] * Mu[GAS][k];
    
    /* b[i] */
    for (k = 0; k < STATE_LAST; k++)
      for (i = 0; i < p->n[k]; i++)
        tmp -= product_element_coef(p->element[j], p->species[k][i]) *
          p->coef[k][i];
    
    /* b[i]o */
    /* 04/06/2000 - division by propellant_mass(e) */
    for (i = 0; i < e->propellant.ncomp; i++)
      tmp += propellant_element_coef(p->element[j],e->propellant.molecule[i]) *
        e->propellant.coef[i] / propellant_mass(e);

    matrix[j + size * size] = tmp;
  }

  /* delta ln(T) */
  if (P != TP)
  {
    for (j = 0; j < p->n[CONDENSED]; j++) /* row */
      matrix[j + idx_cond + size * idx_T] = Ho[CONDENSED][j];
  }
  
  /* right side */
  for (j = 0; j < p->n[CONDENSED]; j++) /* row */
  {
    matrix[j + idx_cond + size * size] = Mu[CONDENSED][j]; 
  }

  /* delta ln(n) */
  matrix[idx_n + size * idx_n] = mol - it->n;
  
  /* delta ln(T) */
  if (P != TP)
  {
    tmp = 0.0;
    for (k = 0; k < p->n[GAS]; k++)
      tmp += p->coef[GAS][k] * Ho[GAS][k];

    matrix[idx_n + size * idx_T] = tmp;
    
  }
  
  /* right side */
  tmp = 0.0;
  for (k = 0; k < p->n[GAS]; k++)
  {
    tmp += p->coef[GAS][k] * Mu[GAS][k];
  }

  matrix[idx_n + size * size] = it->n - mol + tmp;
    
  /* for enthalpy/pressure problem */
  if (P == HP)
  {
    /* part with lagrangian multipliers */
    for (i = 0; i < p->n_element; i++) /* each column */
    {   
      tmp = 0.0;
      for (k = 0; k < p->n[GAS]; k++)
        tmp += p->A[i][k] * p->coef[GAS][k] * Ho[GAS][k];

      matrix[idx_T + size * i] = tmp;
    }

    /* Delta n */
    for (i = 0; i < p->n[CONDENSED]; i++)
      matrix[idx_T + size * ( i + idx_cond)] = Ho[CONDENSED][i];


    /* Delta ln(n) */
    tmp = 0.0;
    for (k = 0; k < p->n[GAS]; k++)
      tmp += p->coef[GAS][k] * Ho[GAS][k];

    matrix[idx_T + size * idx_n] = tmp;

    /* Delta ln(T) */
    tmp = 0.0;
    for (k = 0; k < p->n[GAS]; k++)
      tmp += p->coef[GAS][k]*specific_heat_0( p->species[GAS][k], pr->T );

    for (k = 0; k < p->n[CONDENSED]; k++)
      tmp += p->coef[CONDENSED][k] * specific_heat_0(p->species[CONDENSED][k],
                                                     pr->T);

    for (k = 0; k < p->n[GAS]; k++)
      tmp += p->coef[GAS][k] * Ho[GAS][k] * Ho[GAS][k];

    matrix[idx_T + size * idx_T] = tmp;

    
    /* right side */
    tmp = 0.0;
    tmp = propellant_enthalpy(e)/(R*pr->T) - product_enthalpy(e);
    
    for (k = 0; k < p->n[GAS]; k++)
      tmp += p->coef[GAS][k] * Ho[GAS][k] * Mu[GAS][k];

    matrix[idx_T + size * size] = tmp;
    
  } /* for entropy/pressure problem */
  else if (P == SP)
  {
    /* part with lagrangian multipliers */
    for (i = 0; i < p->n_element; i++) /* each column */
    {   
      tmp = 0.0;
      for (k = 0; k < p->n[GAS]; k++)
        tmp += p->A[i][k] * p->coef[GAS][k] *
          entropy(p->species[GAS][k], GAS, it->ln_nj[k] - it->ln_n, pr->T,
                  pr->P);
      
      matrix[idx_T + size * i] = tmp;
    }
    
    /* Delta n */
    for (i = 0; i < p->n[CONDENSED]; i++)
      matrix[idx_T + size * (i + p->n_element)] =
        entropy_0( p->species[CONDENSED][i], pr->T);
    
    /* Delta ln(n) */
    tmp = 0.0;
    for (k = 0; k < p->n[GAS]; k++)
      tmp += p->coef[GAS][k] * entropy(p->species[GAS][k], GAS, it->ln_nj[k]
                                       - it->ln_n, pr->T, pr->P);

    matrix[idx_T + size * idx_n] = tmp;
    
    tmp = 0.0;
    for (k = 0; k < p->n[GAS]; k++)
      tmp += p->coef[GAS][k]*specific_heat_0( p->species[GAS][k], pr->T );

    for (k = 0; k < p->n[CONDENSED]; k++)
      tmp += p->coef[CONDENSED][k]*
        specific_heat_0( p->species[CONDENSED][k], pr->T);

    for (k = 0; k < p->n[GAS]; k++)
      tmp += p->coef[GAS][k]*Ho[GAS][k]*
        entropy(p->species[GAS][k], GAS, it->ln_nj[k]-it->ln_n, pr->T, pr->P);
    
    matrix[idx_T + size * idx_T] = tmp;    
    
    /* entropy of reactant */
    tmp = e->entropy; /* assign entropy */
    tmp -= product_entropy(e);
    tmp += it->n;

    for (k = 0; k < p->n[GAS]; k++)
      tmp -= p->coef[GAS][k];

    for (k = 0; k < p->n[GAS]; k++)
      tmp += p->coef[GAS][k]
        * Mu[GAS][k]
        * entropy(p->species[GAS][k], GAS, it->ln_nj[k] - it->ln_n,
                  pr->T, pr->P);

    matrix[idx_T + size * size] = tmp;    
  }

  return 0;
}

/* This part of the matrix is the same for equilibrium and derivative */
int fill_matrix(double *matrix, equilibrium_t *e, problem_t P)
{

  short i, j, k, size;
  double tmp;

  product_t *p  = &(e->product);

  short idx_cond, idx_n, idx_T;
  
  int roff = 2;

  if (P == TP)
    roff = 1;

  idx_cond  = p->n_element;
  idx_n     = p->n_element + p->n[CONDENSED];
  idx_T     = p->n_element + p->n[CONDENSED] + 1;
  size      = p->n_element + p->n[CONDENSED] + roff;
  
  /* fill the matrix (part with the Lagrange multipliers) */
  for (i = 0; i < p->n_element; i++)  /* each column */
  {
    for (j = 0; j < p->n_element; j++) /* each row */
    {
      tmp = 0.0;
      for (k = 0; k < p->n[GAS]; k++)
      {
        tmp += p->A[j][k] * p->A[i][k] * p->coef[GAS][k]; 
      }

      matrix[j + size * i] = tmp;
    }
  }
  
  /* Delta n */
  for (i = 0; i < p->n[CONDENSED]; i++) /* column */
  {
    for (j = 0; j < p->n_element; j++) /* row */
    {
      matrix[j + size * (i + idx_cond)] = 
        product_element_coef(p->element[j], p->species[CONDENSED][i]); 
    }
  } 

  /* delta ln(n) */
  for (j = 0; j < p->n_element; j++)
  {
    tmp = 0.0;
    for (k = 0; k < p->n[GAS]; k++)
    {
      tmp += p->A[j][k] * p->coef[GAS][k];
    }
    matrix[j + size * idx_n] = tmp;
  }
   
  /* second row */
  for (i = 0; i < p->n_element; i++) /* column */
  {
    for (j = 0; j < p->n[CONDENSED]; j++) /* row */
    {
      /* copy the symetric part of the matrix */
      matrix[j + idx_cond + size * i] = matrix[i + size * (j + idx_cond)];
    }
  }
  
  /* set to zero */
  for (i = 0; i < p->n[CONDENSED] + 1; i++) /* column */
  {
    for (j = 0; j < p->n[CONDENSED]; j++) /* row */
    {
      matrix[j + idx_cond + size * (i + idx_cond)] = 0.0;
    }
  }
  
  /* third row */
  for (i = 0; i < p->n_element; i++) /* each column */
  {   
    /* copy the symetric part of the matrix */
    matrix[idx_n + size * i] = matrix[i + size * idx_n];
  }

  /* set to zero */
  for (i = 0; i < p->n[CONDENSED]; i++) /* column */
  {
    matrix[idx_n + size * (i + idx_cond)] = 0.0;
  }
  
  return 0;
}

int remove_condensed(short *size, short *n, equilibrium_t *e)
{

  int i, j, k, pos;
  int r = 0; /* something have been replace, 0=false, 1=true */

  int ok = 1;

  product_t       *p  = &(e->product);
  equilib_prop_t  *pr = &(e->properties);
  
  for (i = 0; i < p->n[CONDENSED]; i++)
  {

    /* if a condensed have negative coefficient, we should remove it */
    if (p->coef[CONDENSED][i] <= 0.0)
    {
      if (global_verbose > 1)
      {
        fprintf(outputfile,
                "%s should be remove, negative concentration.\n\n", 
                (thermo_list + p->species[CONDENSED][i])->name );
      }
      
      /* remove from the list ( put it at the end for later use )*/
      pos = p->species[CONDENSED][i];
      
      for (j = i; j < p->n[CONDENSED] - 1; j++)
      {
        p->species[CONDENSED][j] = p->species[CONDENSED][j + 1];
      }
      p->species[CONDENSED][ p->n[CONDENSED] - 1 ] = pos;
        
      (p->n[CONDENSED])--;
      
      //(*size)--; /* reduce the size of the matrix */
      r = 1;
    }
    else if ( !(temperature_check(p->species[CONDENSED][i], pr->T)) )
    {
      /* if the condensed species is present outside of the temperature
         range at which it could exist, we should either replace it by
         an other phase or add the other phase. If the difference between
         the melting point and the temperature is over 50 k, we replace,
         else we add the other phase. */

      /* Find the new molecule */
      for (j = p->n[CONDENSED]; j < (*n); j++)
      {
        /* if this is the same molecule and temperature_check is true,
           than it is the good molecule */

        for (k = 0; k < 5; k++)
        {
          if (!( ((thermo_list + p->species[CONDENSED][i])->coef[k] ==
                  (thermo_list + p->species[CONDENSED][j])->coef[k] ) &&
                 ((thermo_list + p->species[CONDENSED][i])->elem[k] ==
                  (thermo_list + p->species[CONDENSED][j])->elem[k] ) &&
                 (p->species[CONDENSED][i] != p->species[CONDENSED][j])));
          //temperature_check(p->species[CONDENSED][j], pr->T) ))
          {
            ok = 0;
          }
        }

        /* replace or add the molecule */
        if (ok)
        {

          if (fabs(pr->T - transition_temperature(p->species[CONDENSED][j],
                                                  pr->T)) > 50.0)
          {
            /* replace the molecule */
            if (global_verbose > 1)
            {
              fprintf(outputfile, "%s should be replace by %s\n\n",
                      (thermo_list + p->species[CONDENSED][i])->name,
                      (thermo_list + p->species[CONDENSED][j])->name);
            }
            
            pos = p->species[CONDENSED][i];
            p->species[CONDENSED][i] = p->species[CONDENSED][j];
            p->species[CONDENSED][j] = pos;
            
          }
          else
          {
            /* add the molecule */
            if (global_verbose > 1)
            {
              fprintf(outputfile, "%s should be add with %s\n\n",
                      (thermo_list + p->species[CONDENSED][i])->name,
                      (thermo_list + p->species[CONDENSED][j])->name);
            }

            /* to include the species, exchange the value */
            pos = p->species[CONDENSED][ p->n[CONDENSED] ];
            p->species[CONDENSED][ p->n[CONDENSED] ] = 
              p->species[CONDENSED][i];
            p->species[CONDENSED][i] = pos;
    
            p->n[CONDENSED]++;

          }
          

          r = 1; /* A species have been replace */
          
          /* we do not need to continue searching so we break */
          break;
        }

        ok = 1;
      }
    }
  } /* for each condensed */

  
  /* 0 if none remove */
  return r;
}

int include_condensed(short *size, short *n, equilibrium_t *e, 
                      double *sol)
{
  double tmp;
  double temp;
  int    i, j, k;
  int    pos;

  product_t       *p  = &(e->product);
  equilib_prop_t  *pr = &(e->properties);
  
  tmp = 0.0;
  j   = -1;

  /* We include a condensed if it minimize the gibbs free energy and
     if it could exist at the chamber temperature */
  for (i = p->n[CONDENSED] ; i < (*n); i++)
  {
    if (temperature_check(p->species[CONDENSED][i], pr->T))
    {
      temp = 0.0;
      for (k = 0; k < p->n_element; k++)
        temp += sol[k]*product_element_coef(p->element[k], 
                                            p->species[CONDENSED][i]);
      
      if ( gibbs_0(p->species[CONDENSED][i], pr->T) - temp < tmp )
      {
        tmp = gibbs_0(p->species[CONDENSED][i], pr->T) - temp;
        j = i; 
      }
    }
  }

  /* In the case we found a species that minimize the gibbs energy,
     we should include it */
  if (!(j == -1))
  {
    
    if (global_verbose > 1)
    { 
      fprintf(outputfile, "%s should be include\n\n", 
              (thermo_list + e->product.species[CONDENSED][j])->name );
    } 
    
    
    /* to include the species, exchange the value */
    pos = p->species[CONDENSED][ p->n[CONDENSED] ];
    p->species[CONDENSED][ p->n[CONDENSED] ] = 
      p->species[CONDENSED][j];
    p->species[CONDENSED][j] = pos;
    
    p->n[CONDENSED]++;
  
    return 1;
  }
  return 0;
}


int new_approximation(equilibrium_t *e, double *sol, problem_t P)
{
  int i, j;

  /* control factor */
  double lambda1, lambda2, lambda;
  
  double temp;

  product_t       *p  = &(e->product);
  equilib_prop_t  *pr = &(e->properties);
  iteration_var_t *it = &(e->itn);
  
  /* compute the values of delta ln(nj) */
  it->delta_ln_n = sol[ p->n_element + p->n[CONDENSED] ];

  if  (P != TP)
    it->delta_ln_T = sol[p->n_element + p->n[CONDENSED] + 1];
  else
    it->delta_ln_T = 0.0;

  
  for (i = 0; i < p->n[GAS]; i++)
  {
    temp = 0.0;
    for (j = 0; j < p->n_element; j++)
    {
      temp += p->A[j][i] * sol[j];
    }
    
    it->delta_ln_nj[i] =
      - gibbs(p->species[GAS][i], GAS, it->ln_nj[i] - it->ln_n, pr->T, pr->P)
      + temp + it->delta_ln_n
      + enthalpy_0(p->species[GAS][i], pr->T)*it->delta_ln_T;     
  }
  

  lambda2 = 1.0;
  lambda1 = __max(fabs(it->delta_ln_T), fabs(it->delta_ln_n));
  lambda1 = 5 * lambda1;
  
  for (i = 0; i < p->n[GAS]; i++)
  {
    if (it->delta_ln_nj[i] > 0.0)
    {
      if (it->ln_nj[i] - it->ln_n <= LOG_CONC_TOL)
      {
        lambda2 = __min(lambda2,
                        fabs( ((- it->ln_nj[i] + it->ln_n - 9.2103404)
                               /(it->delta_ln_nj[i] - it->delta_ln_n))) );
      }
      else if (it->delta_ln_nj[i] > lambda1) 
      {
        lambda1 = it->delta_ln_nj[i];
      }
    }
  }
  
  lambda1 = 2.0 / lambda1;
  
  lambda = _min(1.0, lambda1, lambda2);
  
  if (global_verbose > 3)
  {
    fprintf(outputfile,
            "lambda  = %.10f\nlambda1 = %.10f\nlambda2 = %.10f\n\n",
            lambda, lambda1, lambda2);
    fprintf(outputfile, "%-19s  nj \t\t  ln_nj_n \t  Delta ln(nj)\n", "");
    
    for (i = 0; i < p->n[GAS]; i++)
    {
      fprintf(outputfile, "%-19s % .4e \t % .4e \t % .4e\n", 
              (thermo_list + p->species[GAS][i])->name, 
              p->coef[GAS][i], it->ln_nj[i], it->delta_ln_nj[i]);
    }
  }

  it->sumn = 0.0;
  
  /* compute the new value for nj (gazeous) and ln_nj */
  for (i = 0; i < p->n[GAS]; i++)
  {
    it->ln_nj[i] = it->ln_nj[i] + lambda * it->delta_ln_nj[i];

    if (it->ln_nj[i] - it->ln_n <= LOG_CONC_TOL)
    {
      p->coef[GAS][i] = 0.0;
    }
    else
    {
      p->coef[GAS][i] = exp(it->ln_nj[i]);
      it->sumn += p->coef[GAS][i];
    }
    
  }
  
  /* compute the new value for nj (condensed) */
  for (i = 0; i < p->n[CONDENSED]; i++)
  {
    p->coef[CONDENSED][i] = p->coef[CONDENSED][i] +
      lambda*sol[p->n_element + i];     
  }

  if (global_verbose > 3)
  {
    for (i = 0; i < p->n[CONDENSED]; i++)
    {
      fprintf(outputfile, "%-19s % .4e\n", 
              (thermo_list + p->species[CONDENSED][i])->name, 
              p->coef[CONDENSED][i]);
    }
  }
    
  /* new value of T */
  if (P != TP)
    pr->T = exp( log(pr->T) + lambda * it->delta_ln_T);
      
  if (global_verbose > 2)
    fprintf(outputfile, "Temperature: %f\n", pr->T);
      
  /* new value of n */
  it->ln_n = it->ln_n + lambda * it->delta_ln_n;
  it->n = exp(it->ln_n);
  
  return SUCCESS;
}
    
bool convergence(equilibrium_t *e, double *sol)
{
  int i;
  double mol;

  /* for convergence test, mol is the sum of all mol
     even condensed */
  
  mol = e->itn.sumn;
      
  /* check for convergence */ 
  for (i = 0; i < e->product.n[GAS]; i++)
  {
    if (!(e->product.coef[GAS][i]*fabs(e->itn.delta_ln_nj[i])/mol <= CONV_TOL))
      return false; /* haven't converge yet */
  }
      
  for ( i = 0; i < e->product.n[CONDENSED]; i++ )
  {
    /* test for the condensed phase */
    if (!(sol[e->product.n_element+1]/mol <= CONV_TOL))
      return false; /* haven't converge yet */
  }
      
  if (!(e->itn.n*fabs(e->itn.delta_ln_n)/mol <= CONV_TOL))
    return false; /* haven't converge yet */

  if (!(fabs(e->itn.delta_ln_T) <= 1.0e-4))
    return false;
  
  return true;
}

int equilibrium(equilibrium_t *equil, problem_t P)
{
  int err_code;
  
  short   i, j, k;
  short   size;     /* size of the matrix */
  double *matrix;
  double *sol;
  
  bool convergence_ok;
  bool stop           = false;
  bool gas_reinserted = false;
  bool solution_ok    = false;

  product_t *p  = &(equil->product);
  
  /* position of the right side of the matrix dependeing on the
     type of problem */
  int roff = 2;

  if (P == TP)
    roff = 1;

  /* initial temperature for assign enthalpy, entropy/pressure */
  if ( P != TP)
    equil->properties.T = ESTIMATED_T;

  
  if (!(equil->product.element_listed))
    /* if the element and the product haven't  been listed */
  {
    list_element(equil);
  }

  if (!(equil->product.product_listed))
  {
    if ((err_code = list_product(equil)) < 0)
    {
      return err_code;
    }
    /* equil->product.n_condensed = equil->product.n[CONDENSED]; */
  }


  /* build up the coefficient matrix */
  for (i = 0; i < p->n_element; i++)
    for (j = 0; j < p->n[GAS]; j++)
      p->A[i][j] = product_element_coef(p->element[i], p->species[GAS][j]);
  
  
  /* First determine an initial estimate of the composition
     to accelerate the convergence */
  /* initial_estimate(equil); */
  
  /* For the first equilibrium, we do not consider the condensed
     species. */
  if (!(equil->product.isequil))
  {
//    equil->product.n_condensed = equil->product.n[CONDENSED];
    equil->product.n[CONDENSED] = 0;
    equil->itn.n = 0.1; /* initial estimate of the mol number */
  }
  
  /* the size of the coefficient matrix */
  size = equil->product.n_element + equil->product.n[CONDENSED] + roff;
  
  /* allocate the memory for the matrix */
  matrix = (double *) malloc (size*(size+1)*sizeof(double));
  
  /* allocate the memory for the solution vector */
  sol = (double *) calloc (size, sizeof(double));

  /* main loop */
  for (k = 0; k < ITERATION_MAX; k++)
  {
    /* Initially we haven't a good solution */
    solution_ok = false;

    while (!solution_ok)
    {      
      fill_equilibrium_matrix(matrix, equil, P);
      
      if (global_verbose > 2)
      {
        fprintf(outputfile, "Iteration %d\n", k+1);
        NUM_print_matrix(matrix, size);
      }
      if (NUM_lu(matrix, sol, size) == -1) /* solve the matrix */
      {
        /* the matrix have no unique solution */
        fprintf(outputfile,
                "The matrix is singular, removing excess condensed.\n");
          
        /* Try removing excess condensed */
        if (!remove_condensed(&size, &(equil->product.n_condensed), equil))
        {
          if (gas_reinserted)
          {
            fprintf(errorfile, "ERROR: No convergence, don't trust results\n");
            /* finish the main loop */
            stop = true;
            break;
          }
          fprintf(errorfile, "None remove. Try reinserting remove gaz\n");
          for (i = 0; i < equil->product.n[GAS]; i++)
          {
            /* It happen that some species were eliminated in the
               process even if they should be present in the equilibrium.
               In such case, we have to reinsert them */
            if (equil->product.coef[GAS][i] == 0.0)
              equil->product.coef[GAS][i] = 1e-6;
          }
          gas_reinserted = true;
        }
        else
        {
          gas_reinserted = false;
        }
        
        /* Restart the loop counter to zero for a new loop */
        k = -1;
      }
      else /* There is a solution */
      {
        solution_ok = true;
      }
    }
      
    if (global_verbose > 2)
    {
      NUM_print_vec(sol, size);    /* print the solution vector */
      fprintf(outputfile, "\n");
    }
    
    /* compute the new approximation */
    new_approximation(equil, sol, P);

    convergence_ok = false;

    /* verify the convergence */
    if (convergence(equil, sol))
    {
      convergence_ok = true;

      if (global_verbose > 0)
      {
        fprintf(outputfile,
                "The solution converge in %-2d iterations (%.2f degK)\n",
                k+1, equil->properties.T);
        //fprintf(outputfile, "T = %f\n", equil->T);
      }
      gas_reinserted = false;


      
      /* print the list of condensed */
      /*
      for (i = 0; i < p->n[CONDENSED]; i++)
      {
        printf("%d ", p->species[CONDENSED][i]);
      }
      printf("\n");
      */
      
      
      /* find if a new condensed species should be include or remove */
      if (remove_condensed(&size, &(equil->product.n_condensed), equil) ||
          include_condensed(&size, &(equil->product.n_condensed), equil, sol))
      {

        free(matrix);

        
        free(sol);
        /* new size */
        size = equil->product.n_element + equil->product.n[CONDENSED] + roff;

        /* allocate the memory for the matrix */
        matrix = (double *) malloc(size*(size+1)*sizeof(double));
                
        /* allocate the memory for the solution vector */
        sol = (double *) malloc (sizeof(double) * size);
          
        /* haven't converge yet */
        convergence_ok = false;    
      }
        
      /* reset the loop counter to compute a new equilibrium */
      k = -1;
    }
    else if (global_verbose > 2)
    {
      fprintf(outputfile, "The solution doesn't converge\n\n");
      /* ?? */
      /*remove_condensed(&size, &n_condensed, equil); */
    }

    if (convergence_ok || stop)
    {
      /* when the solution have converge, we could get out of the
         main loop */
      /* if there was problem, the stop flag is set and we also get out */
      break;
    }
   
    /* suppose that it will converge the next time */
    convergence_ok = true;
    
  } /* end of main loop */

  free (sol);
  free (matrix);
  
  if (k == ITERATION_MAX)
  {
    //fprintf(outputfile, "\n");
    //fprintf(outputfile, "Maximum number of %d iterations attain\n",
    //        ITERATION_MAX);
    //fprintf(outputfile, "Don't thrust results.\n"); 
    return ERR_TOO_MANY_ITER;
  }
  else if (stop)
  {
    //fprintf(outputfile, "\n");
    //fprintf(outputfile, "Problem computing equilibrium...aborted.\n");
    //fprintf(outputfile, "Don't thrust results.\n");
    return ERR_EQUILIBRIUM;
  }

  equil->product.isequil = true;
  equil->equilibrium_ok = true;
  compute_thermo_properties(equil); 
  derivative(equil);
  
  return SUCCESS;
}


