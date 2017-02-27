
#include "num.h"

int Round(double a)
{
  int t = a;
  
  if (a - (double)t < 0.5)
    return t;
  else if (a - (double)t > 0.5)
    return t + 1;
  else
    return (t + (t % 2));
}


double epsilon(void)
{	
  double epsilon;
  double temp;

  epsilon = 1.0;
  temp = epsilon + 1.0;
  
  while (temp > 1.0)
  {
    epsilon = epsilon / 2.0;
    temp = epsilon + 1.0;
  } 
  return (epsilon * 2.0);

}
