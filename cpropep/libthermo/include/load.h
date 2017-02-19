#ifndef load_h
#define load_h

/***************************************************************
FUNCTION: Load the propellant data contain in filename

PARAMETER: filename should be the path of the file containing
           the information on propellant. It should be exactly
           in the format of the file propellant.dat include
	   in the distribution.

COMMENTS: It load the information in the global variable
          propellant_list[MAX_PROPELLANT] that is of type 
	  propellant_t

AUTHOR: Antoine Lefebvre
        modification bye Mark Pinese
****************************************************************/
int load_propellant(char *filename);

/***************************************************************
FUNCTION: Load the thermo data contain in filename

PARAMETER: filename should be the path of the file containing
           the thermo information . It should be exactly
           in the format of the file thermo.dat include
	   in the distribution.

COMMENTS: It load the information in the global variable
          thermo_list[MAX_THERMO] that is of type 
	  thermo_t

AUTHOR: Antoine Lefebvre
        modification bye Mark Pinese
****************************************************************/
int load_thermo(char *filename);

/***************************************************************
Removes trailing ' ' in str.  If str is all ' ', removes all
but the first.
str - pointer to a character array (not necessarily a string)
len - length of str.

AUTHOR: Mark Pinese
****************************************************************/
void trim_spaces(char *str, unsigned int len);


#endif







