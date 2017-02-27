#include <math.h>
#include <stdio.h>

#include "num.h"

/* Resolution of non-linear equation of the for
   f(x) = 0 with the secante method which is a modified
   newton method using secante instead of derivative */
int NUM_sec(double (*f)(double x), double x0, double x1, int nmax,
            double epsilon, double *ans)
{
  int i = 0;
  double x2;
  double delta;
  
  do
  {
    if (i >= nmax)
    {
      return NO_CONVERGENCE;
    }
    
    x2 = x1 - f(x1)*(x1 - x0)/(f(x1) - f(x0));

    delta = fabs(x2 - x1)/fabs(x2);   
    x0 = x1;
    x1 = x2;
    i++;
  } while (delta > epsilon);

  *ans = x2;
  return 0;
}
