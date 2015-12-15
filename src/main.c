/*******************************************************************
 *** Runge-Kutta 4th order method for the solution of 
 *** 2nd order ODEs: main source file.
 *** 
 *** Author: Nikos Tryfonidis, November 2015
 *** The MIT License (MIT)
 *** Copyright (c) 2015 Nikos Tryfonidis
 *** See LICENSE.txt
 *******************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Headers ("headers" directory): */
#include "../headers/struct.h"
#include "../headers/RK4.h"
#include "../headers/io.h"
#include "../headers/memory.h"

int main(int argc, char * argv[])
{
  int i, j, totalSteps, nOutput, interval;
  double h, length;
  double *t, *x, *v;
  struct xv old;
  
  /* Get command line arguments */
  getArgs(argc, argv, &h, &length, &interval);
  
  /* Calculate total timesteps */
  totalSteps = length/h;
  printf("Total steps: %d\n", totalSteps);
  nOutput = totalSteps/interval+1;
  printf("Output size: %d points\n", nOutput);
  
  /* Allocate dependent and independent variable arrays */
  t = allocateDoubleArray1D(nOutput);
  x = allocateDoubleArray1D(nOutput);
  v = allocateDoubleArray1D(nOutput);
  
  /*********************************/
  /*** ODE Solution              ***/
  /*********************************/
  
  /* Set Initial conditions here!  */
  t[0] = 0; x[0] = 0; v[0] = 1.0;  
  
  double t_old;
  t_old = t[0];
  old.x = x[0];
  old.v = v[0];
  
  /* Write solution every [interval] steps */
  for (i=1; i<nOutput; i++) {
    for(j=0; j<interval; j++) {
      t_old += h;
      old = RK4(t_old, old, h, rhs);
    }
    t[i] = t_old;
    x[i] = old.x;
    v[i] = old.v; 
  }

  /* Write File Output */
  writeFileOutput (t, x, v, nOutput, "output/out.txt");
  
  /* Free Memory */
  free(t); free(x); free(v);
  
  return 0;
}

