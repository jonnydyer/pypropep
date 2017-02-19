/* num.h - Library of numerical method
 * $Id: num.h,v 1.3 2001/02/22 19:47:37 antoine Exp $
 * Copyright (C) 2000
 *    Antoine Lefebvre <antoine.lefebvre@polymtl.ca>
 *
 * Licensed under the GPL
 */

#ifndef num_h
#define num_h

/* NOTE on matrix representation
 * -----------------------------
 * All matrix should be allocate the following way
 *
 * matrix = (double *) malloc (sizeof(double) * line * column)
 *
 * to access the value at line 2, column 4, you do it like that
 *
 * matrix[2 + 4*line]
 */


/* This type is used as function pointer for the
 * sysnewton algorithm.
 */
typedef double (*func_t)(double *x);

typedef struct status
{

  int itn;
  
} status_t;

#define NO_CONVERGENCE 1
#define NO_SOLUTION 2


/* Find the solution of linear system of equation using the
 * LU factorisation method as explain in
 * 'Advanced engineering mathematics' bye Erwin Kreyszig.
 *
 * This algorithm will also do column permutation to find
 * the larger pivot.
 *
 * ARGUMENTS
 * ---------
 * matrix: the augmented matrix of coefficient in the system
 *         with right hand side value.
 *
 * solution: the solution vector
 *
 * neq: number of equation in the system
 *
 * Antoine Lefebvre
 *    february 6, 2000 Initial version
 *    october 20, 2000 revision of the permutation method
 */
int NUM_lu(double *matrix, double *solution, int neq);
//int old_lu(double *matrix, double *solution, int neq);

/* This function print the coefficient of the matrix to
 * the screen. 
 *
 * matrix: should be an augmented matrix
 * neq   : number of equation
 */
int NUM_print_matrix(double *matrix, int neq);

/* Print the coefficient of the square matrix to the screen
 */
int NUM_print_square_matrix(double *matrix, int neq);


/* This function print the contents of the vector to
 * the screen. 
 *
 * vec: vector containing neq element
 */
int NUM_print_vec(double *vec, int neq);


/**************************************************************
FUNCTION: This function solve systems of ODE of the first order
          with the Runge-Kutta method of the fourth order.

PARAMETER: The first parameter is a pointer to the function we
           we want to solve. This function take five parameter
	   neq is the number of equations in the system
	   time is the time at which we want to evaluate
	   y is an array containing an initial value
	   dy will store the result of the function
	   ierr any error field

	   step is the time variation
	   duration is the total time of the simulation
	   ic are the initial conditions
	   **y is an array containing all the data
	   
	   void * is nay user data

COMMENTS: **y must be properly allocated, [number of points]X[neq]

          It could be interesting to add a tolerance and use 
	  variable step to reach our tolerance instead of using a 
	  fixed step.

AUTHOR: Antoine Lefebvre

DATE: February 11
*****************************************************************/
int NUM_rk4(int (*f)(int neq, double time, double *y, double *dy, void *data), 
            int neq, double step, double duration, double *ic, 
            double **y, void *data );

int NUM_rkf(int (*f)(int neq, double time, double *y, double *dy, void *data), 
            int neq, double step, double duration, double *ic, 
            double **y, double epsil, void *data);

/* this function return the nearest integer to a */
/* it is a replacement of rint which is not ANSI complient */
int Round(double a);

double epsilon(void);

int NUM_sec(double (*f)(double x), double x0, double x1, int nmax,
            double epsilon, double *ans);

int NUM_newton(double (*f)(double x), double (*df)(double x), double x0,
               int nmax, double epsilon, double *ans);

int NUM_ptfix(double (*f)(double x), double x0,
              double nmax, double epsilon, double *ans);

int NUM_sysnewton(func_t *Jac, func_t *R, double *x, int nvar,
                  int nmax, double eps);


int trapeze(double *data, int n_point, int col, int off, double *integral);

int simpson(double *data, int n_point, int col, int off, double *integral);

int create_spline(double *data, int n_point, double *spline);

#endif



