/*******************************************************************
 *** RK4 Header file: 
 *** Functions: Runge-Kutta 4th order method for the solution of 
 *** 2nd order ODEs.
 *** 
 *******************************************************************/

struct xv RK4(double t, struct xv old, double h, 
              double (*rhs)(double, double, double));
double rhs(double t, double x, double v);
