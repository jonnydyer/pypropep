/* thermo.c  -  Compute thermodynamic properties of individual
                species and composition of species           */
/* $Id: thermo.c,v 1.2 2001/02/22 19:48:44 antoine Exp $ */
/* Copyright (C) 2000                                                  */
/*    Antoine Lefebvre <antoine.lefebvre@polymtl.ca>                   */
/*    Mark Pinese <pinese@cyberwizards.com.au>                         */
/*                                                                     */
/* Licensed under the GPLv2                                            */


#include <math.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>

#include "thermo.h"
#include "compat.h"
#include "conversion.h"

/**************************************************************
These variables hold the number of records for propellant and thermo data
***************************************************************/
unsigned long num_propellant, num_thermo;

/* global variable containing the information about chemical species */
propellant_t	*propellant_list;
thermo_t	    *thermo_list;


/****************************************************************
VARIABLE: Contain the molar mass of element by atomic number
          molar_mass[0] contain hydrogen and so on.
          Data come from Sargent-Welch 1996
*****************************************************************/
const float molar_mass[N_SYMB] = { 
  1.00794,   4.002602, 6.941,      9.012182, 10.811,    12.0107,
  14.00674,  15.9994,  18.9984032, 20.11797, 22.989770, 24.305, 
  26.981538, 28.0855,  30.973761,  32.066,   35.4527,   39.948, 
  39.0983,   40.078,   44.95591,   47.88,    50.9415,   51.996,
  54.938,    55.847,   58.9332,    58.6934,  63.546,    65.39,
  69.723,    72.61,    74.9216,    78.96,    79.904,    83.80,
  85.4678,   87.62,    88.9059,    91.224,   92.9064,   95.94,
  98.0,      101.07,   102.9055,   106.42,   107.868,   112.41,
  114.82,    118.71,   121.757,    127.60,   126.9045,  131.29,
  132.9054,  137.33,   138.9055,   140.12,   140.9077,  144.24,
  145.,      150.36,   151.965,    157.25,   158.9253,  162.50,
  164.9303,  167.26,   168.9342,   173.04,   174.967,   178.49,
  180.9479,  183.85,   186.207,    190.2,    192.22,    195.08,
  196.9665,  200.59,   204.383,    207.2,    208.9804,  209.,
  210.,      222.,     223.,       226.0254, 227.,      232.0381,
  231.0359,  238.029,  237.0482,   244.,     12.011,    9.01218,
  10.811,    24.305,   26.98154,   257.0,    0,         2};


/****************************************************************
VARIABLE: Contain the symbol of the element in the same way as
          for the molar mass.

COMMENTS: It is use in the loading of the data file to recognize
          the chemical formula.
*****************************************************************/
const char symb[N_SYMB][3] = {
  "H ","HE","LI","BE","B ","C ","N ","O ",
  "F ","NE","NA","MG","AL","SI","P ","S ","CL","AR","K ","CA",
  "SC","TI","V ","CR","MN","FE","CO","NI","CU","ZN","GA","GE",
  "AS","SE","BR","KR","RB","SR","Y ","ZR","NB","MO","TC","RU",
  "RH","PD","AG","CD","IN","SN","SB","TE","I ","XE","CS","BA",
  "LA","CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER",
  "TM","YB","LU","HF","TA","W ","RE","OS","IR","PT","AU","HG","TL",
  "PB","BI","PO","AT","RN","FR","RA","AC","TH","PA","U ","NP",
  "U6","U5","U1","U2","U3","U4","FM",
  "E ", "D " }; /* the E stand for electron and D for deuterium*/


/* Enthalpy in the standard state (Dimensionless) */
double enthalpy_0(int sp, float T)
{
  thermo_t *s = (thermo_list + sp);

  double val;
  int    pos = 0, i;
  
  if (T < s->range[0][0]) /* Temperature below the lower range */
  {
    pos = 0;
  }       /*Temperature above the higher range */
  else if (T >= s->range[s->nint-1][1]) 
  {
    pos = s->nint - 1;
  }
  else
  {
    for (i = 0; i < s->nint; i++) /* Find the range */
    {
      if ((T >= s->range[i][0]) && (T < s->range[i][1]))
        pos = i;
    }
  }
  
  /* parametric equation for dimentionless enthalpy */
  val = -s->param[pos][0]*pow(T, -2) + s->param[pos][1]*pow(T, -1)*log(T)
    + s->param[pos][2] + s->param[pos][3]*T/2 + s->param[pos][4]*pow(T, 2)/3
    + s->param[pos][5]*pow(T, 3)/4 + s->param[pos][6]*pow(T, 4)/5
    + s->param[pos][7]/T;

  return val; /* dimensionless enthalpy */
}

/* Entropy in the standard state (Dimensionless)*/
double entropy_0(int sp, float T)
{
  thermo_t *s = (thermo_list + sp);
  double val;
  int    pos = 0, i;

  if (T < s->range[0][0])
  {
    pos = 0;
  }
  else if (T >= s->range[s->nint-1][1])
  {
    pos = s->nint - 1;
  }
  else
  {
    for (i = 0; i < s->nint; i++)
    {
      if ((T >= s->range[i][0]) && (T < s->range[i][1]))
        pos = i;
    }
  }
  
  /* parametric equation for dimentionless entropy */
  val = -s->param[pos][0]*pow(T, -2)/2 - s->param[pos][1]*pow(T, -1)
    + s->param[pos][2]*log(T) + s->param[pos][3]*T
    + s->param[pos][4]*pow(T, 2)/2
    + s->param[pos][5]*pow(T, 3)/3 + s->param[pos][6]*pow(T, 4)/4 
    + s->param[pos][8];
  
  return val;
}

/* Specific heat in the standard state (Dimensionless) */
double specific_heat_0(int sp, float T)
{
  thermo_t *s = (thermo_list + sp);
  double val;
  int    pos = 0, i;

  if (T < s->range[0][0])
  {
    pos = 0;
  }
  else if (T >= s->range[s->nint-1][1])
  {
    pos = s->nint - 1;
  }
  else
  {
    for (i = 0; i < s->nint; i++)
    {
      if ((T >= s->range[i][0]) && (T < s->range[i][1]))
        pos = i;
    }
  }
  
  /* parametric equation for dimentionless specific_heat */
  val = s->param[pos][0]*pow(T, -2) + s->param[pos][1]*pow(T, -1)
    + s->param[pos][2] + s->param[pos][3]*T + s->param[pos][4]*pow(T, 2)
    + s->param[pos][5]*pow(T, 3) + s->param[pos][6]*pow(T, 4);

  return val;
}

/* Dimensionless Gibbs free energy in the standard state */
double gibbs_0(int sp, float T)
{
  return enthalpy_0(sp, T) - entropy_0(sp, T); /* dimensionless */
}

/* Check if the species is in its range of definition
   0 if out of range, 1 if ok */
int temperature_check(int sp, float T)
{
  thermo_t *s = (thermo_list + sp);

  if ((T > s->range[s->nint-1][1]) || (T < s->range[0][0]))
    return 0;

  return 1;
}

/* This function return the transition temperature of the species
   considered which is nearest of the temperature T */
double transition_temperature(int sp, float T)
{
  thermo_t *s = (thermo_list + sp);

  /* first assume that the lowest temperature is the good one */
  double transition_T = s->range[0][0];

  /* verify if we did the good bet */
  if (fabs(transition_T - T) > fabs(s->range[s->nint-1][1] - T))
  {
    transition_T = s->range[s->nint-1][1];
  }

  return transition_T;
}

double entropy(int sp, state_t st, double ln_nj_n, float T, float P)
{
  double s;
  
  switch (st)
  {
    case GAS:
        /* The thermodynamic data are based on a standard state pressure
           of 1 bar (10^5 Pa) */
        s = entropy_0(sp, T) - ln_nj_n - log(P * ATM_TO_BAR);
        break;
    case CONDENSED:
        s = entropy_0(sp, T);
        break;
    default:
        s = 0;
  }  
  return s;
}


/* J/mol T is in K, P is in atm */
double gibbs(int sp, state_t st, double ln_nj_n, float T, float P)
{
  double g;
  
  switch (st)
  {
    case GAS:    
        g = gibbs_0(sp, T) + ln_nj_n + log(P * ATM_TO_BAR);
        break;
    case CONDENSED:
        g = gibbs_0(sp, T);
        break;
    default:
        g = 0;
  }  
  return g;
}

double propellant_molar_mass(int molecule)
{     
  int i = 0, coef;
  double ans = 0;

  while ((coef = (propellant_list + molecule)->coef[i]))
	{
		ans += coef * molar_mass[(propellant_list + molecule)->elem[i]];
		i++;
	}
	return ans;
}

/* J/mol */
double heat_of_formation(int molecule)
{
  double hf = (propellant_list + molecule)->heat * 
    propellant_molar_mass(molecule);
  return hf;
}


/* should not be in thermo.c */
double propellant_enthalpy(equilibrium_t *e)
{
  int i;
  double h = 0.0;
  for (i = 0; i < e->propellant.ncomp; i++)
  {
    h += e->propellant.coef[i] * heat_of_formation (e->propellant.molecule[i])
      / propellant_mass (e);
  }
  return h;
}

/* should not be in thermo.c */
double product_enthalpy(equilibrium_t *e)
{
  int i;
  double h = 0.0;
  
  for (i = 0; i < e->product.n[GAS]; i++)
  {
    h += e->product.coef[GAS][i] * enthalpy_0(e->product.species[GAS][i], e->properties.T);
  }
  
  for (i = 0; i < e->product.n[CONDENSED]; i++)
  {
    h += e->product.coef[CONDENSED][i] * enthalpy_0(e->product.species[CONDENSED][i], e->properties.T);
  }
  return h;
}

/* should not be in thermo.c */
double product_entropy(equilibrium_t *e)
{
  int i;
  double ent = 0.0;
  for (i = 0; i < e->product.n[GAS]; i++)
  {
    ent += e->product.coef[GAS][i]*entropy(e->product.species[GAS][i], GAS,
                                     e->itn.ln_nj[i] - e->itn.ln_n,
                                     e->properties.T, e->properties.P);
  }
  for (i = 0; i < e->product.n[CONDENSED]; i++)
  {
    ent += e->product.coef[CONDENSED][i]*entropy(e->product.species[CONDENSED][i],
                                           CONDENSED, 0, e->properties.T, e->properties.P);
  }
  return ent;
}

/* should not be in thermo.c */
/* The specific heat of the mixture for frozen performance */
double mixture_specific_heat_0(equilibrium_t *e, double temp)
{
  int i;
  double cp = 0.0;
  /* for gases */
  for (i = 0; i < e->product.n[GAS]; i++)
  {
    cp += e->product.coef[GAS][i]*specific_heat_0(e->product.species[GAS][i], temp);
  }
  /* for condensed */
  for (i = 0; i < e->product.n[CONDENSED]; i++)
  {
    cp += e->product.coef[CONDENSED][i]*
      specific_heat_0(e->product.species[CONDENSED][i], temp);
  }
  return cp;
}


int thermo_search(char *str)
{
  int i;
  int last = -1;
  
  for (i = 0; i < num_thermo; i++)
  {
    if (!(STRNCASECMP(str, (thermo_list + i)->name, strlen(str))))
    {
      last = i;
      printf("%-5d %s\n", i, (thermo_list + i)->name);
    }
  }
  return last;
}

int propellant_search(char *str)
{
  int i;
  int last = -1;
  
  for (i = 0; i < num_propellant; i++)
  {
    if (!(STRNCASECMP(str, (propellant_list + i)->name, strlen(str))))
    {
      last = i;
      printf("%-5d %s\n", i, (propellant_list + i)->name);
    }
  }
  return last; 
}


int atomic_number(char *symbole)
{
  int i;
  int element = -1;
  
  /* find the atomic number of the element */
  for (i = 0; i < N_SYMB; i++)
  {
    if (!STRCASECMP(symbole, symb[i]))
    {
      element = i;
      break;
    }
  }
  return element;
}

int compute_density(composition_t *c)
{
  short i;
  double mass = 0;

  c->density = 0.0;
  
  for (i = 0; i < c->ncomp; i++)
  {
    mass += c->coef[i] * propellant_molar_mass(c->molecule[i]);
  }
  
  for (i = 0; i < c->ncomp; i++)
  {
    if ((propellant_list + c->molecule[i])->density != 0.0)
    {
      c->density += c->coef[i] * propellant_molar_mass(c->molecule[i])
        / (mass * (propellant_list + c->molecule[i])->density);
    }
  }
  
  if (c->density != 0.0)
  {
    c->density = 1/c->density;
  }

  return 0;
}

/* This fonction return the offset of the molecule in the propellant_list
   the argument is the chemical formula of the molecule */
int propellant_search_by_formula(char *str)
{
  int i = 0, j ;
  
  char  tmp[5];
  char *ptr;
  
  int   elem[6] = {0, 0, 0, 0, 0, 1};   
  int   coef[6] = {0, 0, 0, 0, 0, 0}; 

  int molecule = -1;
  
  ptr = str; /* beginning of the string */

  while ( (i < 6) && ((ptr - str) < strlen(str)) )
  {    
    if (isupper(*ptr) && islower(*(ptr+1)) && (isupper(*(ptr+2)) ||
                                               iscntrl(*(ptr+2))) )
    {
      tmp[0] = *ptr;
      tmp[1] = toupper(*(ptr+1));
      tmp[2] = '\0';
      /* find the atomic number of the element */
      elem[i] = atomic_number(tmp);
      coef[i] = 1;
      i++;   
      ptr += 2;
    }
    else if (isupper(*ptr) && (isupper(*(ptr+1)) ||
                               iscntrl(*(ptr+1))) )
    {
      tmp[0] = *ptr;
      tmp[1] = ' ';
      tmp[2] = '\0';
      elem[i] = atomic_number(tmp);
      coef[i] = 1;
      i++;
      ptr++;
    }
    else if (isupper(*ptr) && isdigit(*(ptr+1)))
    {
      tmp[0] = *ptr;
      tmp[1] = ' ';
      tmp[2] = '\0';
      elem[i] = atomic_number(tmp);
      
      j = 0;
      do
      {
        tmp[j] = *(ptr + 1 + j);
        j++;
      } while (isdigit(*(ptr + 1 + j)));

      tmp[j] = '\0';
      
      coef[i] = atoi(tmp);
      i++;
      
      ptr = ptr + j + 1;
    }
    else if (isupper(*ptr) && islower(*(ptr+1)) && isdigit(*(ptr+2)))
    {
      tmp[0] = *ptr;
      tmp[1] = toupper(*(ptr+1));
      tmp[2] = '\0';
      elem[i] = atomic_number(tmp);
      
      j = 0;
      while (isdigit(*(ptr + 2 + j)))
      {
        tmp[j] = *(ptr + 1 + j);
        j++;
      }
      tmp[j] = '\0';
      
      coef[i] = atoi(tmp);
      i++;
      
      ptr = ptr + j + 2;
    }
  }

  /*
  for (i = 0; i < 6; i++)
  {
    if (elem[i] != -1)
      printf("%s %d\n", symb[elem[i]], coef[i]);
  }
  */

  for (i = 0; i < num_propellant; i++)
  {
    for (j = 0; j < 6; j++)
    {
      /* set to the same value as the previous one if the same */
      if (!( ((propellant_list+i)->coef[j] == coef[j]) &&
             ((propellant_list+i)->elem[j] == elem[j]) ))
        break;
    }
    
  
  /* Now search in propellant list for this molecule */
/*
  for (j = 0; j < num_propellant; j++)
  {
    for (i = 0; i < 6; i++)
    {
      if ( (coef[i] != propellant_element_coef(elem[i], j)) &&
           (propellant_list + i)
        break;
    }
*/  

    if (j == 5) /* we found the molecule ! */
    {

      /* check if the inverse is true */
      molecule = i;
      break;
    }
  }
  
  return molecule;
}


/* Mass of propellant in gram */
double propellant_mass(equilibrium_t *e)
{
  int i;
  double mass = 0.0;
  for (i = 0; i < e->propellant.ncomp; i++)
  {
    mass += e->propellant.coef[i] *
      propellant_molar_mass(e->propellant.molecule[i]);
  }
  return mass;
}
