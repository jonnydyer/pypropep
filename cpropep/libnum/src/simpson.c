
#include <math.h>

int simpson(double *data, int n_point, int col, int off, double *integral)
{

  int i;
  
  double val = 0;

  double x0, x1, x2, fx0, fx1, fx2;

  double dd_x0x1;
  double dd_x1x2;
  double dd_x0x1x2;

  int beg = 0;
  if ((n_point%2) != 1)
  {
    beg = 1;
    x0 = data[0 + 0*col];
    x1 = data[0 + 1*col];
    fx0 = data[off + 0*col];
    fx1 = data[off + 1*col];
    val += (fx0+fx1)*(x1-x0)/2;
  }
  for (i = beg; i < n_point - 2; i += 2)
  {
    x0 = data[0 + (i+0)*col];
    x1 = data[0 + (i+1)*col];
    x2 = data[0 + (i+2)*col];

    fx0 = data[off + (i+0)*col];
    fx1 = data[off + (i+1)*col];
    fx2 = data[off + (i+2)*col];

    dd_x0x1 = (fx1 - fx0)/(x1 - x0);
    dd_x1x2 = (fx2 - fx1)/(x2 - x1);

    dd_x0x1x2 = (dd_x1x2 - dd_x0x1)/(x2 - x0);
    
    val += (fx0 - dd_x0x1*x0 + dd_x0x1x2*x0*x1)*(x2 - x0) +
      (dd_x0x1 - dd_x0x1x2*(x0+x1))*(pow(x2, 2)-pow(x0, 2))/2 +
      dd_x0x1x2*(pow(x2, 3) - pow(x0,3))/3;
  }

  *integral = val;
  return 0;
}
