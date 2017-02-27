#include <math.h>
#include <stdio.h>

#include "num.h"

int NUM_ptfix(double (*f)(double x), double x0, double nmax,
              double epsilon, double *ans)
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
    
    x1 = f(x0);

    delta = fabs(x1 - x0)/fabs(x1);
    x0 = x1;
    i++;

  } while (delta > epsilon);

  *ans = x1;
  return 0;
}
