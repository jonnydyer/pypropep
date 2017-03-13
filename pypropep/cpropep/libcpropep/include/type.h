#ifndef type_h
#define type_h

#define MAX_PRODUCT 400 /* Maximum species in product */
#define MAX_ELEMENT  15 /* Maximum different element  */
#define MAX_COMP     20 /* Maximum different ingredient in
                           composition */

#include "compat.h"

/****************************************************************
TYPE:  Enumeration of the possible state of a substance
*****************************************************************/
typedef enum 
{
  GAS,
  CONDENSED,
  STATE_LAST
} state_t;

typedef enum
{
  TP,          /* assign temperature and pressure */
  HP,          /* assign enthalpy and pressure */
  SP           /* assign entropy and pressure */
} problem_t;

typedef enum
{
  SUBSONIC_AREA_RATIO,
  SUPERSONIC_AREA_RATIO,
  PRESSURE
} exit_condition_t;


/********************************************
Note: Specific impulse have unit of m/s
      Ns/kg = (kg m / s^2) * (s / kg)
            = (m / s)

      It is habitual to found in literature
      specific impulse in units of second.
      It is in reality Isp/g where g is
      the earth acceleration.
**********************************************/
typedef struct _performance_prop
{
  double ae_at;   /* Exit aera / Throat aera              */   
  double a_dotm;  /* Exit aera / mass flow rate (m/s/atm) */
  double cstar;   /* Characteristic velocity              */
  double cf;      /* Coefficient of thrust                */
  double Ivac;    /* Specific impulse (vacuum)            */
  double Isp;     /* Specific impulse                     */
  
} performance_prop_t;


/***************************************************************
TYPE: Hold the composition of a specific propellant
      ncomp is the number of component
      molecule[ ] hold the number in propellant_list corresponding
                  to the molecule
      coef[ ] hold the stochiometric coefficient

NOTE: It should be great to allocate the memory of the array in 
      function of the number of element

DATE: February 6, 2000
****************************************************************/
typedef struct _composition
{
  short  ncomp;              /* Number of different component */
  short  molecule[MAX_COMP]; /* Molecule code                 */
  double coef[MAX_COMP];     /* Moles of molecule             */ 
  double density;            /* Density of propellant         */
} composition_t;


/*****************************************************************
TYPE: Hold the composition of the combustion product. The molecule
      are separate between their different possible state.

NOTE: This structure should be initialize with the function 
      initialize_product.

DATE: February 13, 2000
******************************************************************/
typedef struct _product
{
  int   element_listed;                 /* true if element have been listed */
  int   product_listed;                 /* true if product have been listed */
  int   isequil;                        /* true if equilibrium is ok        */

  /* coefficient matrix for the gases */ 
  unsigned short A[MAX_ELEMENT][MAX_PRODUCT];
  
  short  n_element;                        /* n. of different element        */
  short  element[MAX_ELEMENT];             /* element list                   */
  short  n[STATE_LAST];                    /* n. of species for each state   */
  short  n_condensed;                      /* n. of total possible condensed */
  short  species[STATE_LAST][MAX_PRODUCT]; /* possible species in each state */
  double coef[STATE_LAST][MAX_PRODUCT];    /* coef. of each molecule         */
  
} product_t;

/* Structure to hold information during the iteration procedure */
typedef struct _iteration_var
{
  double n;                        /* mol/g of the mixture                  */
  double ln_n;                     /* ln(n)                                 */
  double sumn;                     /* sum of all the nj                     */
  double delta_ln_n;               /* delta ln(n) in the iteration process  */
  double delta_ln_T;               /* delta ln(T) in the iteration process  */
  double delta_ln_nj[MAX_PRODUCT]; /* delta ln(nj) in the iteration process */
  double ln_nj[MAX_PRODUCT];       /* ln(nj) nj are the individual mol/g    */

} iteration_var_t;

/**********************************************
Hold information on equilibrium properties once
it have been compute

June 14, 2000
***********************************************/
typedef struct _equilib_prop
{
  double P;    /* Pressure (atm)              */
  double T;    /* Temperature (K)             */
  double H;    /* Enthalpy (kJ/kg)            */
  double U;    /* Internal energy (kJ/kg)     */
  double G;    /* Gibbs free energy (kJ/kg)   */
  double S;    /* Entropy (kJ/(kg)(K))        */
  double M;    /* Molar mass (g/mol)          */
  double dV_P; /* (d ln(V) / d ln(P))t        */
  double dV_T; /* (d ln(V) / d ln(T))p        */
  double Cp;   /* Specific heat (kJ/(kg)(K))  */
  double Cv;   /* Specific heat (kJ/(kg)(K))  */
  double Isex; /* Isentropic exponent (gamma) */
  double Vson; /* Sound speed (m/s)           */
} equilib_prop_t;


typedef struct _new_equilibrium
{  
  int equilibrium_ok;  /* true if the equilibrium have been compute */
  int properties_ok;   /* true if the properties have been compute  */
  int performance_ok;  /* true if the performance have been compute */

  //temporarily
  double entropy;
  
  iteration_var_t    itn;
  composition_t      propellant;
  product_t          product;
  equilib_prop_t     properties;
  performance_prop_t performance;
  
} equilibrium_t;


#endif
