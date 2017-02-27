/* performance.c  -  Compute performance caracteristic of a motor
                     considering equilibrium                      */
/* $Id: performance.c,v 1.2 2001/02/22 19:49:28 antoine Exp $ */
/* Copyright (C) 2000                                                  */
/*    Antoine Lefebvre <antoine.lefebvre@polymtl.ca>                   */
/*    Mark Pinese  <pinese@cyberwizards.com.au>                        */
/*                                                                     */
/* Licensed under the GPLv2                                            */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "performance.h"
#include "derivative.h"
#include "print.h"
#include "equilibrium.h"

#include "const.h"
#include "compat.h"
#include "return.h"
#include "thermo.h"

#define TEMP_ITERATION_MAX  8
#define PC_PT_ITERATION_MAX 5
#define PC_PE_ITERATION_MAX 6


double compute_temperature(equilibrium_t *e, double pressure,
                           double p_entropy);


/* Entropy of the product at the exit pressure and temperature */
double product_entropy_exit(equilibrium_t *e, double pressure, double temp)
{
  double ent;
  double t = e->properties.T;
  double p = e->properties.P;
  e->properties.T = temp;
  e->properties.P = pressure;
  ent = product_entropy(e); /* entropy at the new pressure */
  e->properties.T = t;
  e->properties.P = p;
  return ent;
}

/* Enthalpy of the product at the exit temperature */
double product_enthalpy_exit(equilibrium_t *e, double temp)
{
  double enth;
  double t = e->properties.T;
  e->properties.T = temp;
  enth = product_enthalpy(e);
  e->properties.T = t;
  return enth;
}
    
/* The temperature could be found by entropy conservation with a
   specified pressure */
double compute_temperature(equilibrium_t *e, double pressure, double p_entropy)
{
  int i = 0;

  double delta_lnt;
  double temperature;

  /* The first approximation is the chamber temperature */
  temperature = e->properties.T;
  
  do
  {
    delta_lnt = (p_entropy  - product_entropy_exit(e, pressure, temperature))
      / mixture_specific_heat_0(e, temperature);

    temperature = exp (log(temperature) + delta_lnt);
        
    i++;

  } while (fabs(delta_lnt) >= 0.5e-4 && i < TEMP_ITERATION_MAX);

  if (i == TEMP_ITERATION_MAX)
  {
    fprintf(errorfile,
       "Temperature do not converge in %d iterations. Don't thrust results.\n",
            TEMP_ITERATION_MAX);
  }
  
  return temperature;
}

int frozen_performance(equilibrium_t *e, exit_condition_t exit_type,
                       double value)
{
  int err_code;
  
  short i;

  double sound_velocity;
  double flow_velocity;
  double pc_pt;            /* Chamber pressure / Throat pressure */
  double pc_pe;            /* Chamber pressure / Exit pressure   */
  double log_pc_pe;        /* log(pc_pe)                         */
  double ae_at;            /* Exit aera / Throat aera            */
  double cp_cv;
  double chamber_entropy;
  double exit_pressure = 0;
  
  equilibrium_t *t  = e + 1; /* throat equilibrium */
  equilibrium_t *ex = e + 2; /* exit equilibrium   */
  
  /* find the equilibrium composition in the chamber */
  if (!(e->product.isequil))
  {
    if ((err_code = equilibrium(e, HP)) != SUCCESS)
    {
      fprintf(outputfile,
              "No equilibrium, performance evaluation aborted.\n");
      return err_code;
    }
  }

  /* Simplification due to frozen equilibrium */
  e->properties.dV_T =  1.0;
  e->properties.dV_P = -1.0;
  e->properties.Cp   = mixture_specific_heat_0(e, e->properties.T) * R;
  e->properties.Cv   = e->properties.Cp - e->itn.n * R;
  e->properties.Isex = e->properties.Cp/e->properties.Cv;

  compute_thermo_properties(e);
  
  chamber_entropy  = product_entropy(e);
  
  /* begin computation of throat caracteristic */
  copy_equilibrium(t, e);
  
  cp_cv = e->properties.Cp/e->properties.Cv;

  /* first estimate of Pc/Pt */
  pc_pt = pow (cp_cv/2 + 0.5, cp_cv/(cp_cv - 1));

  i = 0;
  do
  {
   
    t->properties.T = compute_temperature(t, e->properties.P/pc_pt,
                                          chamber_entropy);

    /* Cp of the combustion point assuming frozen */
    t->properties.Cp   = mixture_specific_heat_0(e, t->properties.T) * R;
    /* Cv = Cp - nR  (for frozen) */
    t->properties.Cv   = t->properties.Cp - t->itn.n * R;
    t->properties.Isex = t->properties.Cp/t->properties.Cv;

    compute_thermo_properties(t);
    
    sound_velocity = sqrt(1000 * e->itn.n * R * t->properties.T *
                          t->properties.Isex);
    
    flow_velocity = sqrt(2000*(product_enthalpy(e)*R*e->properties.T -
                               product_enthalpy_exit(e, t->properties.T)*
                               R*t->properties.T));
    
    pc_pt = pc_pt / ( 1 + ((pow(flow_velocity, 2) - pow(sound_velocity, 2))
                           /(1000*(t->properties.Isex + 1)*
                             t->itn.n * R *t->properties.T)));
    i++;
  } while ((fabs((pow(flow_velocity, 2) - pow(sound_velocity, 2))
                 /pow(flow_velocity, 2)) > 0.4e-4) &&
           (i < PC_PT_ITERATION_MAX));

  if (i == PC_PT_ITERATION_MAX)
  {
    fprintf(errorfile,
    "Throat pressure do not converge in %d iterations. Don't thrust results\n",
            PC_PT_ITERATION_MAX);
  }
  
  //printf("%d iterations to evaluate throat pressure.\n", i);
  
  t->properties.P    = e->properties.P/pc_pt;
  t->performance.Isp = t->properties.Vson = sound_velocity;

  /* Now compute exit properties */ 
  copy_equilibrium(ex, e);

  if (exit_type == PRESSURE)
  {
    exit_pressure = value;
  }
  else 
  {
    ae_at = value;
    
    /* Initial estimate of pressure ratio */
    if (exit_type == SUPERSONIC_AREA_RATIO)
    {   
      if ((ae_at > 1.0) && (ae_at < 2.0))
      {
        log_pc_pe = log(pc_pt) + sqrt (3.294*pow(ae_at,2) + 1.535*log(ae_at));
      }
      else if (ae_at >= 2.0)
      {
        log_pc_pe = t->properties.Isex + 1.4 * log(ae_at);
      }
      else
      { 
        printf("Aera ratio out of range ( < 1.0 )\n");
        return ERR_AERA_RATIO;
      }
    }
    else if (exit_type == SUBSONIC_AREA_RATIO)
    {
      if ((ae_at > 1.0) && (ae_at < 1.09))
      {
        log_pc_pe = 0.9 * log(pc_pt) /
          (ae_at + 10.587 * pow(log(ae_at), 3) + 9.454 * log(ae_at));
      }
      else if (ae_at >= 1.09)
      {
        log_pc_pe = log(pc_pt) /
          (ae_at + 10.587 * pow(log(ae_at), 3) + 9.454 * log(ae_at));
      }
      else
      { 
        printf("Aera ratio out of range ( < 1.0 )\n");
        return ERR_AERA_RATIO;
      }
    }
    else
    {
      return ERR_RATIO_TYPE;
    }
    

    /* Improved the estimate */
    i = 0;
    do
    {      
      pc_pe            = exp(log_pc_pe);
      ex->properties.P = exit_pressure   = e->properties.P/pc_pe;
      ex->properties.T = compute_temperature(e, exit_pressure,
                                             chamber_entropy);
      /* Cp of the combustion point assuming frozen */
      ex->properties.Cp   = mixture_specific_heat_0(e, ex->properties.T) * R;
      /* Cv = Cp - nR  (for frozen) */
      ex->properties.Cv   = ex->properties.Cp - ex->itn.n * R;
      ex->properties.Isex = ex->properties.Cp/ex->properties.Cv;
      
      compute_thermo_properties(ex);
    
      sound_velocity = sqrt(1000 * ex->itn.n * R * ex->properties.T *
                            ex->properties.Isex);
    
      ex->performance.Isp =
        flow_velocity = sqrt(2000*(product_enthalpy(e)*R*e->properties.T -
                                   product_enthalpy_exit(e, ex->properties.T)*
                                   R*ex->properties.T));
      
      ex->performance.ae_at =
        (ex->properties.T * t->properties.P * t->performance.Isp) /
        (t->properties.T * ex->properties.P * ex->performance.Isp);

      log_pc_pe = log_pc_pe +
        (ex->properties.Isex*pow(flow_velocity, 2)/
         (pow(flow_velocity, 2) - pow(sound_velocity,2))) *
        (log(ae_at) - log(ex->performance.ae_at));

      i++;
      
    } while ( (fabs((log_pc_pe - log(pc_pe))) > 0.00004) &&
              (i < PC_PE_ITERATION_MAX) );

    if (i == PC_PE_ITERATION_MAX)
    {
      fprintf(errorfile,
    "Exit pressure do not converge in %d iterations. Don't thrust results\n",
              PC_PE_ITERATION_MAX);
    }
    
    //printf("%d iterations to evaluate exit pressure.\n", i);
    
    pc_pe            = exp(log_pc_pe);
    exit_pressure    = e->properties.P/pc_pe;
    
  }
      
  ex->properties.T = compute_temperature(e, exit_pressure,
                                         chamber_entropy);
  /* We must check if the exit temperature is more than 50 K lower
     than any transition temperature of condensed species.
     In this case the results are not good and must be reject. */

  ex->properties.P    = exit_pressure;
  ex->performance.Isp = sqrt(2000*(product_enthalpy(e)*R*e->properties.T -
                                    product_enthalpy_exit(e, ex->properties.T)*
                                    R*ex->properties.T));


  /* units are (m/s/atm) */
  ex->performance.a_dotm = 1000 * R * ex->properties.T * ex->itn.n /
    (ex->properties.P * ex->performance.Isp);

  /* Cp of the combustion point assuming frozen */
  ex->properties.Cp   = mixture_specific_heat_0(e, ex->properties.T) * R;
  /* Cv = Cp - nR  (for frozen) */
  ex->properties.Cv   = ex->properties.Cp - ex->itn.n * R;
  ex->properties.Isex = ex->properties.Cp/ex->properties.Cv;

  compute_thermo_properties(ex);
  
  ex->properties.Vson = sqrt(1000 * e->itn.n * R * ex->properties.T *
                             e->properties.Isex);

  
  t->performance.a_dotm = 1000 * R * t->properties.T
    * t->itn.n / (t->properties.P * t->performance.Isp);
  t->performance.ae_at = 1.0;
  t->performance.cstar = e->properties.P * t->performance.a_dotm;
  t->performance.cf    = t->performance.Isp /
    (e->properties.P * t->performance.a_dotm);
  t->performance.Ivac  = t->performance.Isp + t->properties.P
    * t->performance.a_dotm;

  ex->performance.ae_at =
    (ex->properties.T * t->properties.P * t->performance.Isp) /
    (t->properties.T * ex->properties.P * ex->performance.Isp);
  ex->performance.cstar = e->properties.P * t->performance.a_dotm;
  ex->performance.cf    = ex->performance.Isp /
    (e->properties.P * t->performance.a_dotm);
  ex->performance.Ivac  = ex->performance.Isp + ex->properties.P
    * ex->performance.a_dotm;

  return SUCCESS;
}


int shifting_performance(equilibrium_t *e, exit_condition_t exit_type,
                         double value)
{
  int err_code;
  short i;
  double sound_velocity = 0.0;
  double flow_velocity;
  double pc_pt;
  double pc_pe;
  double log_pc_pe;
  double ae_at;
  double chamber_entropy;
  double exit_pressure = 0;
  
  equilibrium_t *t  = e + 1; /* throat equilibrium */
  equilibrium_t *ex = e + 2; /* throat equilibrium */
  
  /* find the equilibrium composition in the chamber */
  if (!(e->product.isequil))
    /* if the equilibrium have not already been compute */
  {
    if ((err_code = equilibrium(e, HP)) < 0)
    {
      fprintf(outputfile, "No equilibrium, performance evaluation aborted.\n");
      return err_code;
    }
  }

  /* Begin by first aproximate the new equilibrium to be
     the same as the chamber equilibrium */

  copy_equilibrium(t, e);
  
  chamber_entropy = product_entropy(e);
  
  /* Computing throat condition */
  /* Approximation of the throat pressure */
  pc_pt = pow(t->properties.Isex/2 + 0.5,
              t->properties.Isex/(t->properties.Isex - 1) );

  t->entropy = chamber_entropy;
    
  i = 0;
  do
  { 
    t->properties.P = e->properties.P/pc_pt;

    /* We must compute the new equilibrium each time */
    if ((err_code = equilibrium(t, SP)) < 0)
    {
      fprintf(outputfile, "No equilibrium, performance evaluation aborted.\n");
      return err_code;
    }

    sound_velocity = sqrt (1000*t->itn.n*R*t->properties.T*
                           t->properties.Isex);
    
    flow_velocity = sqrt (2000*(product_enthalpy(e)*R*e->properties.T -
                                product_enthalpy(t)*R*t->properties.T));

    pc_pt = pc_pt / ( 1 + ((pow(flow_velocity, 2) - pow(sound_velocity, 2))
                           /(1000*(t->properties.Isex + 1)*t->itn.n*R*
                             t->properties.T)));
    i++;
  } while ((fabs((pow(flow_velocity, 2) - pow(sound_velocity, 2))
                 /pow(flow_velocity, 2)) > 0.4e-4) &&
           (i < PC_PT_ITERATION_MAX));

  if (i == PC_PT_ITERATION_MAX)
  {
    fprintf(errorfile, "Throat pressure do not converge in %d iterations."
            " Don't thrust results.\n", PC_PT_ITERATION_MAX);
  }
  
  //printf("%d iterations to evaluate throat pressure.\n", i);

  t->properties.P    = e->properties.P/pc_pt;
  t->properties.Vson = sound_velocity;
  t->performance.Isp = sound_velocity;


  t->performance.a_dotm = 1000 * R *
    t->properties.T * t->itn.n /
    (t->properties.P * t->performance.Isp);
  
  copy_equilibrium(ex, e);

  if (exit_type == PRESSURE)
  {
    exit_pressure = value;
  }
  else
  {
    ae_at = value;

    /* Initial estimate of pressure ratio */
    if (exit_type == SUPERSONIC_AREA_RATIO)
    {
      if ((ae_at > 1.0) && (ae_at < 2.0))
      {
        log_pc_pe = log(pc_pt) + sqrt (3.294*pow(ae_at,2) + 1.535*log(ae_at));
      }
      else if (ae_at >= 2.0)
      {
        log_pc_pe = t->properties.Isex + 1.4 * log(ae_at);
      }
      else
      { 
        printf("Aera ratio out of range ( < 1.0 )\n");
        return ERR_AERA_RATIO;
      }
    }
    else if (exit_type == SUBSONIC_AREA_RATIO)
    {
      if ((ae_at > 1.0) && (ae_at < 1.09))
      {
        log_pc_pe = 0.9 * log(pc_pt) /
          (ae_at + 10.587 * pow(log(ae_at), 3) + 9.454 * log(ae_at));
      }
      else if (ae_at >= 1.09)
      {
        log_pc_pe = log(pc_pt) /
          (ae_at + 10.587 * pow(log(ae_at), 3) + 9.454 * log(ae_at));
      }
      else
      { 
        printf("Aera ratio out of range ( < 1.0 )\n");
        return ERR_AERA_RATIO;
      }
    }
    else
    {
      return ERR_RATIO_TYPE;
    }

    
    /* Improved the estimate */
    ex->entropy      = chamber_entropy;
    i = 0;
    do
    {
      pc_pe            = exp(log_pc_pe);
      ex->properties.P = exit_pressure    = e->properties.P/pc_pe;

      
      /* Find the exit equilibrium */
      if ((err_code = equilibrium(ex, SP)) < 0)
      {
        fprintf(outputfile,
                "No equilibrium, performance evaluation aborted.\n");
        return err_code;
      }
      
      sound_velocity = ex->properties.Vson;
     
      ex->performance.Isp =
        flow_velocity = sqrt(2000*(product_enthalpy(e)*R*e->properties.T -
                                   product_enthalpy(ex)*R*ex->properties.T));
      
      ex->performance.ae_at =
        (ex->properties.T * t->properties.P * t->performance.Isp) /
        (t->properties.T * ex->properties.P * ex->performance.Isp);

      log_pc_pe = log_pc_pe +
        (ex->properties.Isex*pow(flow_velocity, 2)/
         (pow(flow_velocity, 2) - pow(sound_velocity,2))) *
        (log(ae_at) - log(ex->performance.ae_at));
      i++;
    } while ((fabs((log_pc_pe - log(pc_pe))) > 0.00004) &&
             (i < PC_PE_ITERATION_MAX));

    if (i == PC_PE_ITERATION_MAX)
    {
      fprintf(errorfile, "Exit pressure do not converge in %d iteration."
              " Don't thrust results.\n", PC_PE_ITERATION_MAX);
    }
    
    //printf("%d iterations to evaluate exit pressure.\n", i);
    
    pc_pe            = exp(log_pc_pe);
    exit_pressure    = e->properties.P/pc_pe;
  }
  
  ex->entropy      = chamber_entropy;
  ex->properties.P = exit_pressure;

  /* Find the exit equilibrium */
  if ((err_code = equilibrium(ex, SP)) < 0)
  {
    fprintf(outputfile, "No equilibrium, performance evaluation aborted.\n");
    return err_code;
  }
  
  flow_velocity = sqrt(2000*(product_enthalpy(e)*R*e->properties.T -
                             product_enthalpy(ex)*R*ex->properties.T));

  
  ex->performance.Isp = flow_velocity;

  ex->performance.a_dotm = 1000 * R *
    ex->properties.T * ex->itn.n /
    (ex->properties.P * ex->performance.Isp);
 
  t->performance.ae_at = 1.0;
  t->performance.cstar = e->properties.P * t->performance.a_dotm;
  t->performance.cf    = t->performance.Isp /
    (e->properties.P * t->performance.a_dotm);
  t->performance.Ivac  = t->performance.Isp + t->properties.P
    * t->performance.a_dotm;
  
  ex->performance.ae_at =
    (ex->properties.T * t->properties.P * t->performance.Isp) /
    (t->properties.T * ex->properties.P * ex->performance.Isp);
  ex->performance.cstar = e->properties.P * t->performance.a_dotm;
  ex->performance.cf    = ex->performance.Isp /
    (e->properties.P * t->performance.a_dotm);
  ex->performance.Ivac  = ex->performance.Isp + ex->properties.P
    * ex->performance.a_dotm;
  
  return SUCCESS;
}





