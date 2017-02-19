#ifndef print_h
#define print_h

#include "type.h"

#include "equilibrium.h"
#include "derivative.h"
#include "performance.h"

#define PROPELLANT_NAME(sp) (propellant_list + sp)->name

extern FILE * errorfile;
extern FILE * outputfile;

int print_error_message(int error_code);

/***************************************************************
FUNCTION: Print the information of a specie in the thermo_list

PARAMETER: an integer corresponding to the molecule

AUTHOR: Antoine Lefebvre
***************************************************************/
int print_propellant_info(int sp);
int print_thermo_info(int sp);


/*************************************************************
FUNCTION: Print the content of the respective list with the
          number which refer to the molecule

AUTHOR: Antoine Lefebvre
        modification by Mark Pinese
**************************************************************/
int print_thermo_list(void);
int print_propellant_list(void);

/*************************************************************
FUNCTION: Print the list of condensed species in the product

PARAMETER: a structure of tupe product_t

AUTHOR: Antoine Lefebvre
**************************************************************/
int print_condensed(product_t p);



/*************************************************************
FUNCTION: Print the list of gazeous species in the product

PARAMETER: a structure of tupe product_t

AUTHOR: Antoine Lefebvre
**************************************************************/
int print_gazeous(product_t p);

int print_product_composition(equilibrium_t *e, short npt);

int print_product_properties(equilibrium_t *e, short npt);

int print_propellant_composition(equilibrium_t *e);

//int print_derivative_results(deriv_t *d);
//int print_performance_information(performance_t *p);

int print_performance_information(equilibrium_t *e, short npt);

#endif
