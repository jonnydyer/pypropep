/* test.c - Testing the functionnality of the various method
 * $Id: test.c,v 1.5 2001/06/10 21:06:00 antoine Exp $
 * Copyright (C) 2000
 *    Antoine Lefebvre <antoine.lefebvre@polymtl.ca>
 *
 * Licensed under the GPL
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include <time.h>

#include "num.h"

FILE * errorfile;
FILE * outputfile;

int test_rk4(void);
int test_lu(void);
int test_sysnewton(void);
int test_sec(void);
int test_newton(void);
int test_ptfix(void);
int test_spline(void);

/* g1(x) = x + 1 - ln(x) */
double g1(double x) {
  return x + 1 - log(x); 
}
/* f1(x) = ln(x) - 1 */
double f1(double x) {
  return log(x) - 1; 
}
/* f1'(x) = 1/x */
double df1(double x) {
  return 1/x; 
}
double dg1(double x) {
  return 1 - 1/x;
}
/* f(x) = (x-1)^3 */
double f(double x) {
  return pow(x-1, 3); 
}

int function(int neq, double time, double *y, double *dy, 
             void *data)
{
  
  dy[0] = y[1];
  dy[1] = -9.8;
  dy[2] = y[3];
  dy[3] = 0;
  return 0;
}

/* functions to test the sysnewton algorythm
 *
 * The system to be solve is the following
 * r1(x1, x2) = e^x1 - x2 = 0
 * r2(x1, x2) = x1^2 + x2^2 - 16 = 0
 *
 */
double r1(double *x)
{
  return exp(x[0]) - x[1];
}
double r2(double *x)
{
  return pow(x[0], 2) + pow(x[1], 2) - 16;
}
/* We need partial derivative of these function with each variable
 * dr1_dx1(x1, x2) = e^x1
 * dr1_dx2(x1, x2) = -1
 * dr2_dx1(x1, x2) = 2*x1
 * dr2_dx2(x1, x2) = 2*x2
 */
double dr1_dx1(double *x)
{
  return exp(x[0]);
}
double dr1_dx2(double *x)
{
  return -1;
}
double dr2_dx1(double *x)
{
  return 2*x[0];
}
double dr2_dx2(double *x)
{
  return 2*x[1];
}


int main(void)
{

  
  test_lu();
  test_spline();
 
  test_rk4();
  test_sysnewton();

  test_sec();
  test_newton();
  test_ptfix();
 
  return 0;
}

int test_sec(void)
{
  double ans;
  printf("Testing Secante method\n");
  if (NUM_sec(f1, 1, 2, 100, 0.0001, &ans))
    printf("No solution: error in the method.\n\n");
  else
    printf("Solution: %f\n\n",  ans);
  return 0;
}

int test_newton(void)
{
  double ans;
  printf("Testing Newton method\n");
  if (NUM_newton(f1, df1, 1, 100, 0.0001, &ans))
    printf("No solution: error in the method.\n\n");
  else
    printf("Solution: %f\n\n", ans);
  return 0;
}

int test_ptfix(void)
{
  double ans;
  printf("Testing fixed point method\n");
  if (NUM_ptfix(g1, 1, 100, 0.0001, &ans))
     printf("No solution: error in the method.\n");
  else
    printf("Solution: %f\n\n", ans);
  return 0;
}

int test_sysnewton(void)
{
  func_t *jac;
  func_t *r;
  int nvar = 2;

  double x[2]; /* solution vector, contain initially the initial conditions */
  
  printf("Testing newton method for solving"
         " non-linear system of equations.\n");

  jac = (func_t *) malloc (sizeof(func_t) * nvar * nvar);
  r = (func_t *) malloc (sizeof(func_t) * nvar);

  /* Initialize the functions pointers array */
  r[0] = r1;
  r[1] = r2;

  jac[0] = dr1_dx1;
  jac[1] = dr2_dx1;
  jac[2] = dr1_dx2;
  jac[3] = dr2_dx2;

  /* set the initial estimate */
  x[0] = 2.8;
  x[1] = 2.8;

  /* call the sysnewton function */

  NUM_sysnewton(jac, r, x, 2, 100, 1e-8);

  /* print the solution */

  printf("Solution: x1 = %f, x2 = %f\n", x[0], x[1]);

  printf("\n");
  return 0;
  
}

int test_lu(void)
{
  int i;
  double *matrix;
  double *solution;
  int size = 8;
  
  printf("Testing the LU factorisation algotythm.\n");
  matrix = (double *) malloc (sizeof(double)*size*(size+1));
  solution = (double *) malloc (sizeof(double)*size);
  
  matrix[0] = 4.77088e-02; matrix[8]  = 1.17204e-01; matrix[16] = 1.88670e-02;
  matrix[1] = 1.17204e-01; matrix[9]  = 4.07815e-01; matrix[17] = 1.25752e-02;
  matrix[2] = 1.88670e-02; matrix[10] = 1.25752e-02; matrix[18] = 4.40215e-02;
  matrix[3] = 0.00000e+00; matrix[11] = 0.00000e+00; matrix[19] = 0.00000e+00;
  matrix[4] = 0.00000e+00; matrix[12] = 0.00000e+00; matrix[20] = 0.00000e+00;
  matrix[5] = 1.00000e+00; matrix[13] = 0.00000e+00; matrix[21] = 3.00000e+00;
  matrix[6] = 2.97603e-02; matrix[14] = 8.67031e-02; matrix[22] = 2.51546e-02;
  matrix[7] = 5.15962e+01; matrix[15] = 1.67940e+02; matrix[23] = -7.46975e+01;

  matrix[24] = 0.00000e+00;
  matrix[25] = 0.00000e+00;
  matrix[26] = 0.00000e+00;
  matrix[27] = 1.28681e-02;
  matrix[28] = 0.00000e+00;
  matrix[29] = 0.00000e+00;
  matrix[30] = 6.43406e-03;
  matrix[31] = -7.23674e-01;
  
  matrix[32] = 0.0000e+00; matrix[40] = 1.00000e+00; matrix[48] = 2.97603e-02;
  matrix[33] = 0.0000e+00; matrix[41] = 0.00000e+00; matrix[49] = 8.67031e-02;
  matrix[34] = 0.0000e+00; matrix[42] = 3.00000e+00; matrix[50] = 2.51546e-02;
  matrix[35] = 0.0000e+00; matrix[43] = 0.00000e+00; matrix[51] = 6.43406e-03;
  matrix[36] = 0.0000e+00; matrix[44] = 2.00000e+00; matrix[52] = 0.00000e+00;
  matrix[37] = 2.0000e+00; matrix[45] = 0.00000e+00; matrix[53] = 0.00000e+00;
  matrix[38] = 0.0000e+00; matrix[46] = 0.00000e+00; matrix[54] = 1.45616e-02;
  matrix[39] = 0.0000e+00; matrix[47] =-9.68727e+03; matrix[55] =-3.06747e+01; 
  
  matrix[56] = 5.15962e+01; matrix[64] =-1.12677e+01;
  matrix[57] = 1.67940e+02; matrix[65] = 8.29437e+00;
  matrix[58] =-7.46975e+01; matrix[66] =-7.39145e+01;
  matrix[59] =-7.23674e-01; matrix[67] =-5.99249e-01;
  matrix[60] = 0.00000e+00; matrix[68] = 1.58804e-02;
  matrix[61] =-9.68727e+03; matrix[69] =-9.62706e+03;
  matrix[62] =-3.06747e+01; matrix[70] =-4.63904e+01;
  matrix[63] = 4.68590e+05; matrix[71] = 2.37298e+05;
  
  //NUM_matscale(matrix, size);
  NUM_print_matrix(matrix, size);

  if (NUM_lu(matrix, solution, size))
    printf("No solution: Error in the numerical method,\n");
  else
    NUM_print_vec(solution, size);
/*
  for (i = 0; i < size; i++)
  {
    if (solution[i] != 1.0)
    {
      printf("Error found in the solution.\n");
    }
  }
*/
  
  printf("\n");
  return 0;
}


int test_rk4(void)
{
  int i, n;
  double *ans;
  double *ic;
  
  printf("\nTesting the RK4 and RKF algorythm.\n");
  
  ic = (double *) malloc(sizeof(double) * 4);
  
  ic[0] = 0;
  ic[1] = 100;
  ic[2] = 0;
  ic[3] = 10;

  /* it return the length of the answer vector */
  //n = NUM_rk4 (function, 4, 0.1, 10, ic, &ans, NULL);

  //for (i = 0; i < n; i++)
  //{
  //  printf("%f %f %f %f \n", ans[4*i], ans[1 + 4*i], ans[2+4*i], ans[3+4*i]);
  //}

  n = NUM_rkf (function, 4, 0.1, 20, ic, &ans, 1e-4, NULL);

  printf("n = %i\n", n);

  for( i = 0; i < n; i++)
  {
    printf("%f  %f  %f  %f  %f\n", ans[4 + 5*i], ans[5*i],
           ans[1+5*i], ans[2 + 5*i], ans[3 + 5*i]);  
  }

  free(ic);
  free(ans);

  return 0;
}

int test_spline(void)
{
  int i;
  
  double *data;

  double *spline;

  int size = 10;

  double ans;
  
  printf("Testing the spline function.\n\n");
  data = malloc (2 * size * sizeof(double));
  spline = malloc (size * sizeof(double));
  
  data[0] = 0; /* x0 */
  data[1] = 55; /* y0 */
  data[2] = 5; /* x1 */
  data[3] = 60; /* y1 */
  data[4] = 10; /* ... */
  data[5] = 58;
  data[6] = 15;
  data[7] = 54;
  data[8] = 20;
  data[9] = 55;
  data[10] = 25;
  data[11] = 60;
  data[12] = 30;
  data[13] = 54;
  data[14] = 35;
  data[15] = 57;
  data[16] = 40;
  data[17] = 52;
  data[18] = 45;
  data[19] = 49;
  
  create_spline(data, size, spline);

  for (i = 0; i < 45; i++)
  {
    eval_spline(data, spline, size, (double)i, &ans);
    printf("%d %f\n", i, ans);
  }
    
  printf("Spline test finish\n");
  return 0;
}

