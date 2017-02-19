#include <math.h>
#include <stdio.h>

#include "num.h"

int NUM_newton(double (*f)(double x), double (*df)(double x), double x0,
               int nmax, double eps, double *ans)
{

  int i = 0;
  double x1;
  double delta;
  
  do
  {
    if (i >= nmax)
    {
      return NO_CONVERGENCE;
    }
    
    x1 = x0 - f(x0)/df(x0);
    delta = fabs(x1 - x0)/fabs(x1);
    
    x0 = x1;
    i++;
 
  }  while (delta > eps);

  *ans = x1;
  return 0;
}
