#ifndef thermo_h
#define thermo_h

#include "equilibrium.h"
#include "const.h"

/* MACRO: Number of symbol in the symbol table */
#define N_SYMB      102

/***************************************************************
TYPE: Structure to hold information of species contain in the
      thermo data file
****************************************************************/
typedef struct _thermo
{
  char    name[19];
  char    comments[57];
  int     nint;         /* number of different temperature interval */
  char    id[7];        /* identification code */
  int     elem[5]; 
  int     coef[5];
  state_t state;
  double  weight;       /* molecular weight */
  float   heat;         /* heat of formation at 298.15 K  (J/mol)  */
  double  dho;          /* HO(298.15) - HO(0) */
  float   range[4][2];  /* temperature range */
  int     ncoef[4];     /* number of coefficient for Cp0/R   */
  int     ex[4][8];     /* exponent in empirical equation */
  
  double param[4][9];
  
  /* for species with data at only one temperature */
  /* especially condensed                          */
  float temp;
  float enth;
  
} thermo_t;

/***************************************************************
TYPE: Structure to hold information of species contain in the
      propellant data file
****************************************************************/
typedef struct _propellant
{
  char  name[120]; /* name of the propellant */
  int   elem[6];   /* element in the molecule (atomic number) max 6 */
  int   coef[6];   /* stochiometric coefficient of this element 
		                  (0 for none) */
  float heat;      /* heat of formation in Joule/gram */
  float density;   /* density in g/cubic cm */
  
} propellant_t;


extern propellant_t	*propellant_list;
extern thermo_t	    *thermo_list;

extern const float molar_mass[];
extern const char symb[][3];

extern unsigned long num_thermo;
extern unsigned long num_propellant;

/*************************************************************
FUNCTION: Search in the field name of thermo_list and return
          the value of the found item.

PARAMETER: A string corresponding to what we search, 
           example: "CO2"

COMMENTS: If nothing is found, it return -1

AUTHOR: Antoine Lefebvre
        modification bye Mark Pinese
**************************************************************/
int thermo_search(char *str);

int propellant_search(char *str);

int atomic_number(char *symbole);

int propellant_search_by_formula(char *str);

/*************************************************************
FUNCTION: Return the enthalpy of the molecule in thermo_list[sp]
          at the temperature T in K. (Ho/RT)

PARAMETER: sp is the position in the array of the molecule
           T is the temperature in K

COMMENTS: It use the parametric form explain in the documentation
          to compute the value from the data read in the file
	  thermo.dat

AUTHOR: Antoine Lefebvre
**************************************************************/
double enthalpy_0(int sp, float T);

/*************************************************************
FUNCTION: Return the entropy of the molecule in thermo_list[sp]
          at the temperature T in K. (So/RT)

PARAMETER: sp is the position in the array of the molecule
           T is the temperature in K

COMMENTS: It use the parametric form explain in the documentation
          to compute the value from the data read in the file
	  thermo.dat

AUTHOR: Antoine Lefebvre
**************************************************************/
double entropy_0(int sp, float T);

double entropy(int sp, state_t st, double ln_nj_n, float T, float P);

/*************************************************************
FUNCTION: Return the specific heat (Cp) of the molecule in 
          thermo_list[sp] at the temperature T in K. (Cp/RT)

PARAMETER: sp is the position in the array of the molecule
           T is the temperature in K

COMMENTS: It use the parametric form explain in the documentation
          to compute the value from the data read in the file
	  thermo.dat

AUTHOR: Antoine Lefebvre
**************************************************************/
double specific_heat_0(int sp, float T);

double mixture_specific_heat_0(equilibrium_t *e, double temp);

/*************************************************************
FUNCTION: Return true if the thermochemical data are define for
          this temperature.

PARAMETER: The same as for entropy

COMMENTS:  It is useful to determine if a specie is present at
           a given temperature.

AUTHOR: Antoine Lefebvre
**************************************************************/
int temperature_check(int sp, float T);

double transition_temperature(int sp, float T);

double propellant_enthalpy(equilibrium_t *e);
double product_enthalpy(equilibrium_t *e);
double product_entropy(equilibrium_t *e);
double propellant_mass(equilibrium_t *e);

int compute_density(composition_t *c);


/*************************************************************
FUNCTION: Return the gibbs free energy of the molecule in 
          thermo_list[sp] at temperature T. (uo/RT)

PARAMETER: sp is the position in the array of the molecule
           T is the temperature in K

COMMENTS: g = H - ST where H is the enthalpy, T the temperature
          and S the entropy.
**************************************************************/
double gibbs_0(int sp, float T);


/*************************************************************
FUNCTION: Return the gibbs free energy of the molecule in 
          thermo_list[sp] at temperature T, pressure P. (u/RT)

PARAMETER: sp is the position in the array of the molecule
           T is the temperature in K

COMMENTS: g = uo + ln(nj/n) + ln(P) for gazes
          g = uo for condensed

AUTHOR: Antoine Lefebvre
**************************************************************/
//double gibbs(int sp, state_t st, double nj, double n, float T, float P);
double gibbs(int sp, state_t st, double nj_n_n, float T, float P);


/***************************************************************
FUNCTION: Return the heat of formation of a propellant in kJ/mol
****************************************************************/
double heat_of_formation(int molecule);

/*************************************************************
FUNCTION: Return the molar mass of a propellant (g/mol)

PARAMETER: molecule is the number in propellant_list
**************************************************************/
double propellant_molar_mass(int molecule);



double propellant_mass(equilibrium_t *e);


#endif
