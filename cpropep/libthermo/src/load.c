#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "equilibrium.h"
#include "load.h"
#include "thermo.h"

#include "conversion.h"
#include "return.h"


/***************************************************************************
Initial format of thermo.dat:
interval   variable   type	size	description
-----------------------------------------------------------------------------
(0, 18)    name	      string	18	compound name
(18, 73)   comments   string	55	comment
(73, 75)   nint	      int	2	the number of temperature intervals
(75, 81)   id	      string	6	the material id
81	   state      int	1	0 - GAS, else CONDENSED
(82, 95)   weight     float	13	molecular weight
(95, 108)  enth/heat  float	13	enthaply if nint == 0 
                                        else heat of formation
			...
			rest of file
			...
***************************************************************************/

int load_thermo(char *filename)
{
  FILE *fd;
  
  int i = 0;
  int j, k, l;

  bool ok;
  
  char buf[88], *buf_ptr, tmp[32], *tmp_ptr;
  buf_ptr = &buf[0];
  tmp_ptr = &tmp[0];
  
  /* open the file for reading */
  if ((fd = fopen(filename, "r")) == NULL )
    return ERR_FOPEN;

  if (global_verbose)
  {
    printf("Scanning thermo data file...");
		fflush(stdout);
  }

  num_thermo = 0;


  /* Scan thermo.dat to find the number of positions in thermo_list
     to allocate */
	while (fgets(buf_ptr, 88, fd))
	{
    /*
      All that is required is to count the number of lines not
      starting with ' ', '!' or '-'
    */
		if (*buf_ptr != ' ' && *buf_ptr != '!' && *buf_ptr != '-')
			num_thermo++;
	}
  
	/* Reset the file pointer */
	fseek(fd, 0, SEEK_SET);

	if (global_verbose)
	{
		printf("\nScan complete.  %ld records found.  Allocating memory...",
           num_thermo);
	}

	if ((thermo_list = (thermo_t *)malloc (sizeof(thermo_t) * num_thermo)) ==
      NULL)
	{
		printf("\n\nMemory allocation error with thermo_t thermo_list[%ld], %ld bytes required", num_thermo, sizeof(thermo_t) * num_thermo);
		return ERR_MALLOC;
	}

	if (global_verbose)
	{
		printf("\nSuccessful.  Loading thermo data file...");
		fflush(stdout);
	}

	for (i = 0; i < num_thermo; i++)
	{
		/* Read in the next line and check for EOF */
		if (!fgets(buf_ptr, 88, fd))
		{
			fclose(fd);
			free(thermo_list);
			return ERR_EOF;
		}

		/* Skip commented lines */
		while (*buf_ptr == '!')
		{
			if (!fgets(buf_ptr, 88, fd))
			{
				fclose(fd);
				free(thermo_list);
				return ERR_EOF;
			}
		}

		/* Read in the name and the comments */
		strncpy((thermo_list + i)->name, buf_ptr, 18);
		trim_spaces((thermo_list + i)->name, 18);
        
		strncpy((thermo_list + i)->comments, buf_ptr + 18, 55);
		trim_spaces((thermo_list + i)->comments, 55);
      
		// Read in the next line and check for EOF
		if (!fgets(buf_ptr, 88, fd))
		{
			fclose(fd);
			free(thermo_list);
			return ERR_EOF;
		}
      
		strncpy(tmp_ptr, buf_ptr, 3);
		(thermo_list + i)->nint = atoi(tmp_ptr);
      
		strncpy((thermo_list + i)->id, buf_ptr + 3, 6);
		trim_spaces((thermo_list + i)->id, 6);
      
		/* get the chemical formula and coefficient */
		/* grep the elements (5 max) */
		for (k = 0; k < 5; k++)
		{
			tmp[0] = buf[k * 8 + 10];
			tmp[1] = buf[k * 8 + 11];
			tmp[2] = '\0';

			/* Check for an empty place (no more atoms) */
			if (strcmp(tmp, "  "))
			{
				/* Atoms still to be processed */
		    
				/* find the atomic number of the element */
        (thermo_list + i)->elem[k] = atomic_number(tmp);
		    
				/* And the number of atoms */
				strncpy(tmp_ptr, buf_ptr + k * 8 + 13, 6);
				tmp[6] = '\0';

				/* Should this be an int?  If so, why is it stored in x.2 format? */
				(thermo_list + i)->coef[k] = (int) atof(tmp_ptr);
			}
			else
			{
				/* No atom here */
				(thermo_list + i)->coef[k] = 0;
			}
		}
	       
		/* grep the state */
		if (buf[51] == '0')
			(thermo_list + i)->state = GAS;
		else
			(thermo_list + i)->state = CONDENSED;
      
		/* grep the molecular weight */
		strncpy(tmp_ptr, buf_ptr + 52, 13);
		tmp[13] = '\0';
		(thermo_list + i)->weight = atof(tmp_ptr);
      
		/* grep the heat of formation (J/mol) or enthalpy if condensed */
		/* The values are assigned in the if block following */
		strncpy(tmp_ptr, buf_ptr + 65, 15);
		tmp[15] = '\0';
      
		/* now get the data */
		/* there is '(thermo_list + i)->nint' set of data */
		if ((thermo_list + i)->nint == 0)
		{
			/* Set the enthalpy */
			(thermo_list + i)->enth = atof(tmp_ptr);
          
			/* condensed phase, different info */
			/* Read in the next line and check for EOF */
			if (!fgets(buf_ptr, 88, fd))
			{
				fclose(fd);
				free(thermo_list);
				return ERR_EOF;
			}
			  
			/* treat the line */
			/* get the temperature of the assigned enthalpy */
			strncpy(tmp_ptr, buf_ptr + 1, 10);
			tmp[10] = '\0';

			(thermo_list + i)->temp = atof(tmp_ptr);
		}
		else 
		{ 
			/* Set the heat of formation */
			(thermo_list + i)->heat = atof(tmp_ptr);


			/* I'm not quite sure this is necessary */
			/* if the value is 0 and this is the same substance as
			the previous one but in a different state ... */
			if ((thermo_list + i)->heat == 0 && i != 0)
			{
        ok = true;
				for (j = 0; j < 5; j++)
				{
					/* set to the same value as the previous one if the same */
					if (!((thermo_list+i)->coef[j] == (thermo_list+i-1)->coef[j] &&
                (thermo_list+i)->elem[j] == (thermo_list+i-1)->elem[j]))
            ok = false;
						 
				}
        if (ok)
          (thermo_list+i)->heat = (thermo_list+i-1)->heat;
			}
            
			for (j = 0; j < (thermo_list + i)->nint; j++)
			{
				/* Get the first line of three */
				/* Read in the line and check for EOF */
				if (!fgets(buf_ptr, 88, fd))
				{
					fclose(fd);
					free(thermo_list);
					return ERR_EOF;
				}
              
				/* low */
				strncpy(tmp_ptr, buf_ptr + 1, 10);
				tmp[10] = '\0';
				(thermo_list + i)->range[j][0] = atof(tmp_ptr);
	  
				/* high */
				strncpy(tmp_ptr, buf_ptr + 11, 10);
				tmp[10] = '\0';
				(thermo_list + i)->range[j][1] = atof(tmp_ptr);
	  
				tmp[0] = buf[22];
				tmp[1] = '\0';
				(thermo_list + i)->ncoef[j] = atoi(tmp_ptr);
	  
				/* grep the exponent */
				for (l = 0; l < 8; l++)
				{
					strncpy(tmp_ptr, buf_ptr + l * 5 + 23, 5);
					tmp[5] = '\0';					     
					(thermo_list + i)->ex[j][l] = atoi(tmp_ptr);
				}
	  
				/* HO(298.15) -HO(0) */
				strncpy(tmp_ptr, buf_ptr + 65, 15);
				tmp[15] = '\0';
				(thermo_list + i)->dho = atof(tmp);
	  
				/* Get the second line of three */
				/* Read in the line and check for EOF */
				if (!fgets(buf_ptr, 88, fd))
				{
					fclose(fd);
					free(thermo_list);
					return ERR_EOF;
				}
			       
				/* grep the first data line */
				/* there are 5 coefficients */
				for (l = 0; l < 5; l++)
				{
					strncpy(tmp_ptr, buf_ptr + l * 16, 16);
					tmp[16] = '\0';
	    
					(thermo_list + i)->param[j][l] = atof(tmp_ptr);
          //(thermo_list + i)->param[j][l] = strtod(tmp_ptr, NULL);
        }
	  
				/* Get the third line of three */
				/* Read in the line and check for EOF */
				if (!fgets(buf_ptr, 88, fd))
				{
					fclose(fd);
					free(thermo_list);
					return ERR_EOF;
				}
	  
				/* grep the second data line */
				for (l = 0; l < 2; l++)
				{
					strncpy(tmp_ptr, buf_ptr + l * 16, 16);
					tmp[16] = '\0';
	    
					(thermo_list + i)->param[j][l + 5] = atof(tmp_ptr);
				}
	  
				for (l = 0; l < 2; l++)
				{
					strncpy(tmp_ptr, buf_ptr + l * 16 + 48, 16);
					tmp[16] = '\0';
	    
					(thermo_list + i)->param[j][l + 7] = atof(tmp_ptr);	    
				}
			}
		}
	}
	
	fclose(fd);
  
	if (global_verbose)
		printf("%d species loaded.\n", i);
  
	return i;
}


int load_propellant(char *filename) 
{
  
  FILE *fd;
  
  int i = 0, j, len, name_start, name_end, name_len;
  
  /* temporary string to store string in order to treat the informations */
  char buf[88], *buf_ptr, tmp[70], *tmp_ptr;
  buf_ptr = &buf[0];
  tmp_ptr = &tmp[0];
  
  /* open the file for reading */
  if ((fd = fopen(filename, "r")) == NULL )
		return ERR_FOPEN;

  if (global_verbose)
  {
    printf("Scanning propellant data file...");
    fflush(stdout);
  }

  num_propellant = 0;
  
  /* Scan propellant.dat to find the number of positions in propellant_list
     to allocate */
	while (fgets(buf_ptr, 88, fd))
	{
		/* All that is required is to count the number of lines not starting
       with '*' or '+' */
		if (*buf_ptr != '*' && *buf_ptr != '+')
			num_propellant++;
	}

	/* Reset the file pointer */
	fseek(fd, 0, SEEK_SET);

	if (global_verbose)
	{
		printf("\nScan complete.  %ld records found.  Allocating memory...",
           num_propellant);
		fflush(stdout);
	}

	if ((propellant_list = (propellant_t *) malloc(sizeof(propellant_t) *
                                                 num_propellant)) == NULL)
	{
		printf ("\n\nMemory allocation error with propellant_t propellant_list[%ld], %ld bytes required", num_propellant, sizeof(propellant_t) * num_propellant);
		return ERR_MALLOC;
	}

	if (global_verbose)
	{
		printf("\nSuccessful.  Loading propellant data file...");
		fflush(stdout);
	}

	if (!fgets(buf_ptr, 88, fd))
	{
		fclose(fd);
		free(propellant_list);
		return ERR_EOF;
	}

  
	for (i = 0; i < num_propellant; i++)
	{
		/* Skip commented code */
		do
		{
			if (!fgets(buf_ptr, 88, fd))
			{
				fclose(fd);
				free(propellant_list);
				return ERR_EOF;
			}
		}
		while (*buf_ptr == '*');

		/* Check for a continued name */
		while (*buf_ptr == '+')
		{
			/* A continued name found */
			strncpy(tmp_ptr, buf_ptr + 9, 70);

 
			/* Find the end of the whitespaces.  name_start + 1 is used to leave
         one space. */
			for (name_start = 0; name_start < 70; name_start++)
			{
				if (*(tmp_ptr + name_start + 1) != ' ')
					break;
			}

			/* Find the end of the name.  > 0 is used to be consistent with the
         one space left */
			/* when finding name_start */
			for (name_end = 69; name_end > 0; name_end--)
			{
				if (*(tmp_ptr + name_end) != ' ')
					break;
			}

			name_len = name_end - name_start + 1;
			len = strlen((propellant_list + i - 1)->name);
      
			/* Check for room in the destination string.  Take into account
         the possibility of a
         multiple line continuation
         TODO - > 120 or >= 120 */
			if (len + name_len >= 120)
			{
				/* Not enough room - copy as much as possible and leave the
           name alone */
				strncpy((propellant_list + i - 1)->name + len,
                tmp_ptr + name_start, 119 - len);
				*((propellant_list + i - 1)->name + 119) = '\x0';
			}
			else
			{
				/* Concatenate the entire string */
				strncpy((propellant_list + i - 1)->name + len,
                tmp_ptr + name_start, name_len);
				*((propellant_list + i - 1)->name + len + name_len) = '\x0';
			}

      
			/* Processing of this line is done, so get the next one */
			if (!fgets(buf_ptr, 88, fd))
			{
				fclose(fd);
				free(propellant_list);
				return ERR_EOF;
			}
		}
		
		/* grep the name */
		strncpy((propellant_list + i)->name, buf_ptr + 9, 30);
		trim_spaces((propellant_list + i)->name, 30);
      
		for (j = 0; j < 6; j++)
		{
			tmp[0] = buf[j * 5 + 39];
			tmp[1] = buf[j * 5 + 40];
			tmp[2] = buf[j * 5 + 41];
			tmp[3] = '\0';
		
			(propellant_list + i)->coef[j] = atoi(tmp);
        
			tmp[0] = buf[j * 5 + 42];
			tmp[1] = buf[j * 5 + 43];
			tmp[2] = '\0';
        
			/* find the atomic number of the element */
/*
      for (k = 0; k < N_SYMB; k++)
			{
				if (!(strcmp(tmp, symb[k]))) 
				{
					(propellant_list + i)->elem[j] = k;
					break;
				}
			}
*/
      (propellant_list + i)->elem[j] = atomic_number(tmp);
		}
      
		strncpy(tmp_ptr, buf_ptr + 69, 5);
		tmp[5] = '\0';		    
		propellant_list[i].heat = atof(tmp) * CAL_TO_JOULE;
      
		strncpy(tmp_ptr, buf_ptr + 75, 5);
		tmp[5] = '\0';
		propellant_list[i].density = atof(tmp) *  LBS_IN3_TO_G_CM3;
      
	} 
  
	fclose(fd);

	if (global_verbose)
		printf("%d species loaded.\n", i);

	return i;
}


void trim_spaces(char *str, unsigned int len)
{
  unsigned int i;
  
  for (i = len - 1; i > 0; i--)
  {
    if (*(str + i) != ' ' && *(str + i) != '\t')
    {
      *(str + i + 1) = '\0';
      return;
    }
  }
  *(str + 1) = '\0';
}
