#ifndef equilibrium_h
#define equilibrium_h
/* equilibrium.h  -  Calculation of Complex Chemical Equilibrium       */
/* $Id: equilibrium.h,v 1.1 2000/10/13 19:24:31 antoine Exp $ */
/* Copyright (C) 2000                                                  */
/*    Antoine Lefebvre <antoine.lefebvre@polymtl.ca>                   */
/*    Mark Pinese <pinese@cyberwizards.com.au>                         */
/*                                                                     */
/* Licensed under the GPLv2                                            */

#include "compat.h"
#include "type.h"

#define GRAM_TO_MOL(g, sp)   g/propellant_molar_mass(sp)

#define _min(a, b, c) __min( __min(a, b), c)
#define _max(a, b, c) __max( __max(a, b), c)

extern int global_verbose;


/***************************************************************
FUNCTION PROTOTYPE SECTION
****************************************************************/

int set_verbose(equilibrium_t *e, int v);

/************************************************************
FUNCTION: This function search for all elements present in
          the composition and fill the list with the 
	  corresponding number.

PARAMETER: e is of type equilibrium_t and hold the information
           about the propellant composition.

COMMENTS: It fill the member element in equilibrium_t

DATE: February 6, 2000

AUTHOR: Antoine Lefebvre
**************************************************************/
int list_element(equilibrium_t *e);
int reset_element_list(equilibrium_t *e);

int list_product(equilibrium_t *e);

/***************************************************************
FUNCTION: This function initialize the equilibrium structure.
          The function allocate memory for all the structure
	  it need. It is important to call dealloc_equilibrium
	  after.

AUTHOR:   Antoine Lefebvre

DATE: February 27, 2000
****************************************************************/
int initialize_equilibrium(equilibrium_t *e);


/***************************************************************
FUNCTION: Dealloc what have been allocated by 
          initialize_equilibrium
***************************************************************/
int dealloc_equilibrium(equilibrium_t *e);

int reset_equilibrium(equilibrium_t *e);

int copy_equilibrium(equilibrium_t *dest, equilibrium_t *src);

int compute_thermo_properties(equilibrium_t *e);

/***************************************************************
FUNCTION: Set the state at which we want to compute the 
          equilibrium.

PARAMETER: e is a pointer to an equilibrium_t structure
           T is the temperature in deg K
	   P is the pressure in atm

AUTHOR:    Antoine Lefebvre
****************************************************************/
int set_state(equilibrium_t *e, double T, double P);


/***************************************************************
FUNCTION: Add a new molecule in the propellant

PARAMETER: e is a pointer to the equilibrium_t structure
           sp is the number of the molecule in the list
	   mol is the quantity in mol

AUTHOR:    Antoine Lefebvre
****************************************************************/
int add_in_propellant(equilibrium_t *e, int sp, double mol);

/***************************************************************
FUNCTION: Return the stochiometric coefficient of an element
          in a molecule. If the element isn't present, it return 0.

COMMENTS: There is a different function for the product and for the
          propellant.

AUTHOR:   Antoine Lefebvre
****************************************************************/
int product_element_coef(int element, int molecule);
//int propellant_element_coef(int element, int molecule);



/***************************************************************
FUNCTION: This function fill the matrix in function of the data
          store in the structure equilibrium_t. The solution
	  of this matrix give corresction to initial estimate.

COMMENTS: It use the theory explain in 
          "Computer Program for Calculation of Complex Chemical
	  Equilibrium Compositions, Rocket Performance, Incident
	  and Reflected Shocks, and Chapman-Jouguet Detonations"
	  by Gordon and McBride

AUTHOR:   Antoine Lefebvre
****************************************************************/


//#ifdef TRUE_ARRAY
//int fill_equilibrium_matrix(double *matrix, equilibrium_t *e, problem_t P);
int fill_matrix(double *matrix, equilibrium_t *e, problem_t P);
//#else
//int fill_equilibrium_matrix(double **matrix, equilibrium_t *e, problem_t P);
//int fill_matrix(double **matrix, equilibrium_t *e, problem_t P);
//#endif


/****************************************************************
FUNCTION: This function compute the equilibrium composition at
          at specific pressure/temperature point. It use fill_matrix
	  to obtain correction to initial estimate. It correct the 
	  value until equilibrium is obtain.

AUTHOR:   Antoine Lefebvre
******************************************************************/
int equilibrium(equilibrium_t *equil, problem_t P);


double product_molar_mass(equilibrium_t *e);

#endif









