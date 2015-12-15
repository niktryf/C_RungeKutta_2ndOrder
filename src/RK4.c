/*******************************************************************
 *** Runge-Kutta 4th order method for the solution of 
 *** 2nd order ODEs: 
 *** Runge-Kutta method function and "source" function.
 ***
 ***  Author: Nikos Tryfonidis, November 2015
 *** The MIT License (MIT)
 *** Copyright (c) 2015 Nikos Tryfonidis
 *** See LICENSE.txt
 *******************************************************************/
#include <stdio.h>
#include <math.h>

#include "../headers/struct.h"

/* 
  Force function (enter desired function here).
  Contains k, l, m coefficients, for defining
  general functions (i.e. F(t, x, v) = kt + lx + mv),
  that may depend on time, position and velocity.
*/
double F(double t, double x, double v) 
{
  double k, l, m;
  k = 0.0;
  l = 1.0;
  m = 0.1;

  return -l*x -m*v;
}

/* Right-hand side function of the 2nd order ODE. 
   This version sticks to Newton's 2nd law formalism
   (x''(t) = Force/mass)
*/
double rhs(double t, double x, double v) 
{
  double mass;
  mass = 1.0;

  return F(t, x, v)/mass;
}

/*
  Runge - Kutta 4th order for 2nd Order ODEs.
  (Follows Newton's 2nd Law formalism)
  t : independent variable
  x : "position" variable
  v : "velocity" variable
  
  Returns new, calculated solutions for x and v,
  nicely packed into a structure.  
*/

struct xv RK4(double t, struct xv old, double h, 
              double (*rhs)(double, double, double))
{
  int i;
  double  k[4], l[4];
  struct xv new;
  
  k[0] = h*old.v;
  l[0] = h*rhs(t, old.x, old.v);

  k[1] = h*(old.v + 0.5*l[0]);
  l[1] = h*rhs(t + 0.5*h, old.x + 0.5*h, old.v + 0.5*k[0]);

  k[2] = h*(old.v + 0.5*l[1]);
  l[2] = h*rhs(t + 0.5*h, old.x + 0.5*h, old.v + 0.5*k[1]);

  k[3] = h*((old.v + l[2]));
  l[3] = h*rhs(t + h, old.x + h, old.v + k[2]);
  
  new.x = old.x + (1.0/6.0)*(k[0] + 2*k[1] + 2*k[2] + k[3]);
  new.v = old.v + (1.0/6.0)*(l[0] + 2*l[1] + 2*l[2] + l[3]);  

  return new;
}
