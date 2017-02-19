#include <stdlib.h>
#include <math.h>
#include "num.h"


int NUM_rk4(int (*f)(int neq, double time, double *y, double *dy, void *data), 
            int neq, double step, double duration, double *ic, 
            double **y, void *data)
{
  int i;
  int n;

  int status;
  
  int length;
  
  double t = 0.0;

  double *tmp;
  double *dy;
  double *K1, *K2, *K3, *K4;   

  double *a;
  
  tmp = (double *) malloc(sizeof(double) * neq);
  dy  = (double *) malloc(sizeof(double) * neq);

  K1  = (double *) malloc(sizeof(double) * neq);
  K2  = (double *) malloc(sizeof(double) * neq);
  K3  = (double *) malloc(sizeof(double) * neq);
  K4  = (double *) malloc(sizeof(double) * neq);

  /* allocation of the answer vector */
  length = (int)ceil(duration/step) + 1;
  a = *y = (double *) malloc(sizeof(double) * neq * length);
  
  for (i = 0; i < neq; i++)
  {
    tmp[i] = a[i + neq*0] = ic[i]; /* initials conditions */
  }
 
  for (n = 0; n < Round(duration/step); n++)
  {

    for (i = 0; i < neq; i++)
    {
      status = f(neq, t, tmp, dy, data);

      if (status)
        return status;
        
      K1[i] = step * dy[i];
      
      tmp[i] = a[i + neq*n] + K1[i]/2;  /* for the next step */           
    }

    /* FIXME: verify the coefficient t + step/2 */
    for (i = 0; i < neq; i++)
    {
      status = f(neq, t + step/2, tmp, dy, data);

      if (status)
        return -1;
        
      K2[i] = step * dy[i];
      
      tmp[i] = a[i + neq*n] + K2[i]/2;
    }
    
    for (i = 0; i < neq; i++)
    {
      status = f(neq, t + step/2, tmp, dy, data);

      if (status)
        return status;
        
      K3[i] = step * dy[i];
      
      tmp[i] = a[i + neq*n] + K3[i];
    }
    
    for (i = 0; i < neq; i++)
    {
      status = f(neq, t + step, tmp, dy, data);
        if (status)
          return -1;
          
      K4[i] = step * dy[i];
    }
    
    for (i = 0; i < neq; i++)
      a[i + neq*(n+1)] = a[i + neq*n] +
        (1.0/6.0)*(K1[i] + 2.0*K2[i] + 2.0*K3[i] + K4[i]);

    t = t + step;
  }
  
  free(tmp);
  free(dy);
  free(K1);
  free(K2);
  free(K3);
  free(K4);
  
  return length;
}

