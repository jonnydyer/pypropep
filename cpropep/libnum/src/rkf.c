#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "num.h"

/* neq:      number of equations
 * step:     initial time step
 * duration: total duration of the simulation
 * ic:       initial conditions vector
 * y:        solution vector
 * epsil:    precision on the data
 * data:     pointer to data structure
 */

/* Runge-Kutta-Fehlberg 4-5 order */
 
int NUM_rkf(int (*f)(int neq, double time, double *y, double *dy, void *data), 
            int neq, double step, double duration, double *ic, 
            double **y, double epsil, void *data)
{
  int i;
  int n;

  int col;
  
  double h;
  double t = 0.0;

  double *a;
  
  double *tmp;
  double *dy;
  double *K1, *K2, *K3, *K4, *K5, *K6;   

  double *E; /* Error vector of the error */
  
  double err; /* maximum error */
  double beta;
  
  tmp = (double *) malloc(sizeof(double) * neq);
  dy  = (double *) malloc(sizeof(double) * neq);

  E   = (double *) malloc(sizeof(double) * neq);
  
  K1  = (double *) malloc(sizeof(double) * neq);
  K2  = (double *) malloc(sizeof(double) * neq);
  K3  = (double *) malloc(sizeof(double) * neq);
  K4  = (double *) malloc(sizeof(double) * neq);
  K5  = (double *) malloc(sizeof(double) * neq);
  K6  = (double *) malloc(sizeof(double) * neq);

  
  h = step;
  n = 0;

  col = neq + 1;
  
  a = *y = (double *) malloc(sizeof(double) * col * (n + 1));
 
  
  for (i = 0; i < neq; i++)
  {
    tmp[i] = a[i + col*0] = ic[i]; /* initial conditions */
  }
  
  while (t < duration)
  {  
    a = *y = (double *) realloc(*y, sizeof(double) * col * (n + 2));
    
    for (i = 0; i < neq; i++)
    {
      f(neq, t, tmp, dy, data);
      K1[i] = h * dy[i];
      
      tmp[i] = a[i + col*n] + K1[i]/4.0;  /* for the next step */          
    }
      
    for (i = 0; i < neq; i++)
    {
      f(neq, t + h/4.0, tmp, dy, data);
      K2[i] = h * dy[i];

      tmp[i] = a[i + col*n] + 3.0*K1[i]/32.0 + 9.0*K2[i]/32.0;
        
    }

    for (i = 0; i < neq; i++)
    {
      f(neq, t + 3.0*h/8.0, tmp, dy, data);
      K3[i] = h * dy[i];

      tmp[i] = a[i + col*n] + 1932.0*K1[i]/2197.0 - 7200.0*K2[i]/2197.0 +
        7296.0*K3[i]/2197.0;
    }

    
    for (i = 0; i < neq; i++)
    {
      f(neq, t + 12.0*h/13.0, tmp, dy, data);
      K4[i] = h * dy[i];
      
      tmp[i] = a[i + col*n] + 439.0*K1[i]/216.0 - 8.0*K2[i] +
        3680.0*K3[i]/513.0 - 845.0*K4[i]/4104.0;
    }

    for (i = 0; i < neq; i++)
    {
      f(neq, t + h, tmp, dy, data);
      K5[i] = h * dy[i];
      
      tmp[i] = a[i + col*n] - 8.0*K1[i]/27.0 + 2.0*K2[i] -
        3544.0*K3[i]/2565.0 + 1859.0*K4[i]/4104.0 - 11.0*K5[i]/40.0;
    }
    
    for (i = 0; i < neq; i++)
    {
      f(neq, t + h/2.0, tmp, dy, data);
      K6[i] = h * dy[i];
    }

    err = 0.0;
    for (i = 0; i <  neq; i++)
    {
      E[i] = fabs(K1[i]/360.0 - 128.0*K3[i]/4275.0 - 2197.0*K4[i]/75240.0 +
                   K5[i]/50.0 + 2.0*K6[i]/55.0); /* /h ?? */
      //printf("E[%d] = %f\n", i, E[i]);
      err = ((E[i] > err) ? E[i] : err);      
    }

    //printf("err = %e\n", err);

    if ((err < epsil) || (h <= step/1000) )
    {
      t += h;
      
      for (i = 0; i < neq; i++)
      {
        a[i + col*(n+1)] = a[i + col*n] + 25.0*K1[i]/216.0 +
          1408.0*K3[i]/2565.0 + 2197.0*K4[i]/4104.0 - 0.2*K5[i];
      }
      /* store the time */
      a[neq + col*(n+1)] = t;
      n++;
    }
    
    beta = pow(epsil/(2*err), 0.25);
    
    if (beta < 0.1)
    {
      h = 0.1* h;
    }
    else if (beta > 4)
    {
      h = 4.0*h;
    }
    else
    {
      h = beta * h;
    }

    /* we have to prevent too little h*/
    if (h < step/1000)
      h = step/1000;

    /* prevent too big */
    if (h > duration/16)
      h = duration/16;
    
    if (t + h > duration)
    {
      h =  duration - t;
    }
    
  }

  return n+1;
}

