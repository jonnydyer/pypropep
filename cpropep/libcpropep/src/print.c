/* print.c  -  Output functions           */
/* $Id: print.c,v 1.2 2001/02/22 19:49:28 antoine Exp $ */
/* Copyright (C) 2000                                                  */
/*    Antoine Lefebvre <antoine.lefebvre@polymtl.ca>                   */
/*    Mark Pinese <pinese@cyberwizards.com.au>                         */
/*                                                                     */
/* Licensed under the GPLv2                                            */

#include <stdio.h>

#include "print.h"
#include "performance.h"
#include "equilibrium.h"
#include "conversion.h"
#include "thermo.h"
#include "const.h"

char header[][32] = {
  "CHAMBER",
  "THROAT",
  "EXIT" };

char err_message[][64] = {
  "Error allocating memory",
  "Error opening file",
  "Error EOF",
  "Error memory not allocated",
  "Error too much product",
  "Error in equilibrium",
  "Error bad aera ratio",
  "Error bad aera ratio type"};

FILE * errorfile;
FILE * outputfile;

int print_error_message(int error_code)
{
  fprintf(errorfile, "%s\n", err_message[-error_code - 1]);
  return 0;
}

int print_propellant_info(int sp)
{
  int j;

  if (sp > num_propellant || sp < 0)
    return -1;
  
  fprintf(outputfile, "Code %-35s Enthalpy  Density  Composition\n", "Name");
  fprintf(outputfile, "%d  %-35s % .4f % .2f", sp,
          (propellant_list + sp)->name,
          (propellant_list + sp)->heat,
          (propellant_list + sp)->density);
  
  fprintf(outputfile, "  ");

  /* print the composition */
  for (j = 0; j < 6; j++)
  {
    if (!((propellant_list + sp)->coef[j] == 0))
      fprintf(outputfile, "%d%s ", (propellant_list + sp)->coef[j],
             symb[ (propellant_list + sp)->elem[j] ]);
  }
  fprintf(outputfile, "\n");
  return 0;
}

int print_thermo_info(int sp)
{
  int   i, j;
  thermo_t *s;

  if (sp > num_thermo || sp < 0)
    return -1;

  s = (thermo_list + sp);
  
  fprintf(outputfile, "---------------------------------------------\n");
  fprintf(outputfile, "Name: \t\t\t%s\n", s->name);
  fprintf(outputfile, "Comments: \t\t%s\n", s->comments);
  fprintf(outputfile, "Id: \t\t\t%s\n", s->id);
  fprintf(outputfile, "Chemical formula:\t");
  
  for (i = 0; i < 5; i++)
  {
    if (!(s->coef[i] == 0))
      fprintf(outputfile, "%d%s", s->coef[i], symb[ s->elem[i]]);
  }
  fprintf(outputfile, "\n");
  fprintf(outputfile, "State:\t\t\t");
  switch (s->state)
  {
    case GAS:
        fprintf(outputfile, "GAZ\n");
        break;
    case CONDENSED:
        fprintf(outputfile, "CONDENSED\n");
        break;
    default:
        printf("UNRECOGNIZE\n");
  }
  
  fprintf(outputfile, "\n");
  fprintf(outputfile, "Molecular weight: \t\t% f g/mol\n", s->weight);
  fprintf(outputfile, "Heat of formation at 298.15 K : % f J/mol\n", s->heat);
  fprintf(outputfile, "Assign enthalpy               : % f J/mol\n", s->enth);
  fprintf(outputfile, "HO(298.15) - HO(0): \t\t% f J/mol\n", s->dho);
  fprintf(outputfile, "Number of temperature range: % d\n\n", s->nint);
  
  for (i = 0; i < s->nint; i++)
  {
    fprintf(outputfile, "Interval: %f - %f \n", s->range[i][0],
            s->range[i][1]);
    for (j = 0; j < 9; j++)
      fprintf(outputfile, "% .9e ", s->param[i][j]);
    fprintf(outputfile, "\n\n");
  }
  fprintf(outputfile, "---------------------------------------------\n");
  return 0;
}


int print_thermo_list(void)
{
  int i;
  for (i = 0; i < num_thermo; i++)
    fprintf(outputfile, "%-4d %-15s % .2f\n", i, (thermo_list + i)->name,
            (thermo_list + i)->heat);
  
  return 0;
}

int print_propellant_list(void)
{
  int i;
  for (i = 0; i < num_propellant; i++)
    fprintf(outputfile, "%-4d %-30s %5f\n", i, (propellant_list + i)->name,
            (propellant_list +i)->heat);
 
  return 0;
}


int print_condensed(product_t p)
{
  int i;
  for (i = 0; i < p.n[CONDENSED]; i ++)
    fprintf(outputfile, "%s ",
            (thermo_list + p.species[CONDENSED][i])->name );
  fprintf(outputfile, "\n");
  return 0;
}

int print_gazeous(product_t p)
{
  int i;
  for (i = 0; i < p.n[GAS]; i++)
    fprintf(outputfile, "%s ", (thermo_list + p.species[GAS][i])->name);
  fprintf(outputfile, "\n");
  return 0;
}

int print_product_composition(equilibrium_t *e, short npt)
{
  int i, j, k;

  double mol_g = e->itn.n;

  /* we have to build a list of all condensed species present
     in the three equilibrium */
  int n = 0;
  int condensed_list[MAX_PRODUCT];

  /* ok become false if the species already exist in the list */
  int ok = 1;

  double qt;
  
  for (i = 0; i < e->product.n[CONDENSED]; i++)
    mol_g += e->product.coef[CONDENSED][i];
  
  fprintf(outputfile, "\nMolar fractions\n\n");
  for (i = 0; i < e->product.n[GAS]; i++)
  {
    if (e->product.coef[GAS][i]/e->itn.n > 0.0)
    {
      fprintf(outputfile, "%-20s",
              (thermo_list + e->product.species[GAS][i])->name);

      for (j = 0; j < npt; j++)
        fprintf(outputfile, " %11.4e", (e+j)->product.coef[GAS][i]/mol_g);
      fprintf(outputfile,"\n");
      
    }
  }

  /* build the list of condensed */
  for (i = 0; i < npt; i++)
  {
    for (j = 0; j < (e+i)->product.n[CONDENSED]; j++)
    {
      for (k = 0; k < n; k++)
      {
        /* check if the condensed are to be include in the list */
        if (condensed_list[k] == (e+i)->product.species[CONDENSED][j])
        {
          /* we do not have to include the species */
          ok = 0;
          break;
        }
      } /* k */

      if (ok)
      {
        condensed_list[n] = (e+i)->product.species[CONDENSED][j];
        n++;
      }
      
      /* reset the flag */
      ok = 1;
    } /* j */
  } /* i */

  
  if (n > 0)
  {
    fprintf(outputfile, "Condensed species\n");
    for (i = 0; i < n; i++)  
    {
      fprintf(outputfile,   "%-20s",
              (thermo_list + condensed_list[i])->name);

      for (j = 0; j < npt; j++)
      {
        /* search in the product of each equilibrium if the
           condensed is present */
        
        qt = 0.0;
        
        for (k = 0; k < (e+j)->product.n[CONDENSED]; k++)
        {
          if (condensed_list[i] == (e+j)->product.species[CONDENSED][k])
          {
            qt = (e+j)->product.coef[CONDENSED][k];
            break;
          }
        }
          
        fprintf(outputfile, " %11.4e", qt/mol_g);
      }
      fprintf(outputfile,"\n");
      
    }
  }
  fprintf(outputfile, "\n");
  return 0;
}


int print_propellant_composition(equilibrium_t *e)
{
  int i, j;
  
  fprintf(outputfile, "Propellant composition\n");
  fprintf(outputfile, "Code  %-35s mol    Mass (g)  Composition\n", "Name");
  for (i = 0; i < e->propellant.ncomp; i++)
  {
    fprintf(outputfile, "%-4d  %-35s %.4f %.4f ", e->propellant.molecule[i],
            PROPELLANT_NAME(e->propellant.molecule[i]), e->propellant.coef[i], 
            e->propellant.coef[i] *
            propellant_molar_mass(e->propellant.molecule[i]));
    
    fprintf(outputfile, "  ");
    /* print the composition */
    for (j = 0; j < 6; j++)
    {
      if (!((propellant_list + e->propellant.molecule[i])->coef[j] == 0))
        fprintf(outputfile, "%d%s ",
                (propellant_list + e->propellant.molecule[i])->coef[j],
                symb[(propellant_list + e->propellant.molecule[i])->elem[j]]);
    }
    fprintf(outputfile, "\n");
  }
  fprintf(outputfile, "Density : % .3f g/cm^3\n", e->propellant.density); 

  if (e->product.element_listed)
  {
    fprintf(outputfile, "%d different elements\n", e->product.n_element);
    /* Print those elements */
    for (i = 0; i < e->product.n_element; i++)
      fprintf(outputfile, "%s ", symb[e->product.element[i]] );
    fprintf(outputfile, "\n");
  }
  
  fprintf(outputfile, "Total mass: % f g\n", propellant_mass(e));
  
  fprintf(outputfile, "Enthalpy  : % .2f kJ/kg\n",
          propellant_enthalpy(e));
  
  fprintf(outputfile, "\n");

  if (e->product.product_listed)
  {
    fprintf(outputfile, "%d possible gazeous species\n", e->product.n[GAS]);
    if (global_verbose > 1)
      print_gazeous(e->product);
    fprintf(outputfile, "%d possible condensed species\n\n",
            e->product.n_condensed);
    if (global_verbose > 1)
      print_condensed(e->product);
  }
  
  return 0;
}

int print_performance_information (equilibrium_t *e, short npt)
{
  short i;
  
  fprintf(outputfile, "Ae/At            :            ");
  for (i = 1; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->performance.ae_at);
  fprintf(outputfile, "\n");
  
  fprintf(outputfile, "A/dotm (m/s/atm) :            ");
  for (i = 1; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->performance.a_dotm);
  fprintf(outputfile, "\n");

  fprintf(outputfile, "C* (m/s)         :            ");
  for (i = 1; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->performance.cstar);
  fprintf(outputfile, "\n");

  fprintf(outputfile, "Cf               :            ");
  for (i = 1; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->performance.cf);
  fprintf(outputfile, "\n");

  fprintf(outputfile, "Ivac (m/s)       :            ");
  for (i = 1; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->performance.Ivac);
  fprintf(outputfile, "\n");

  fprintf(outputfile, "Isp (m/s)        :            ");
  for (i = 1; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->performance.Isp);
  fprintf(outputfile, "\n");

  fprintf(outputfile, "Isp/g (s)        :            ");
  for (i = 1; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->performance.Isp/Ge);
  fprintf(outputfile, "\n");

  return 0;
}


int print_product_properties(equilibrium_t *e, short npt)
{
  short i;

  fprintf(outputfile, "                  ");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " %11s", header[i]);
  fprintf(outputfile, "\n");
  
  fprintf(outputfile, "Pressure (atm)   :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.3f", (e+i)->properties.P);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "Temperature (K)  :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.3f", (e+i)->properties.T);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "H (kJ/kg)        :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.3f", (e+i)->properties.H);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "U (kJ/kg)        :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.3f", (e+i)->properties.U);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "G (kJ/kg)        :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.3f", (e+i)->properties.G);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "S (kJ/(kg)(K)    :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.3f", (e+i)->properties.S);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "M (g/mol)        :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.3f", (e+i)->properties.M);
  fprintf(outputfile, "\n");
  
  fprintf(outputfile, "(dLnV/dLnP)t     :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->properties.dV_P);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "(dLnV/dLnT)p     :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->properties.dV_T);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "Cp (kJ/(kg)(K))  :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->properties.Cp);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "Cv (kJ/(kg)(K))  :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->properties.Cv);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "Cp/Cv            :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->properties.Cp/(e+i)->properties.Cv);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "Gamma            :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->properties.Isex);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "Vson (m/s)       :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->properties.Vson);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "\n");
  return 0;
}
