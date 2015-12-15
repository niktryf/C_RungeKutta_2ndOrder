/*******************************************************************
 *** Runge-Kutta 4th order method for the solution of 
 *** 2nd order ODEs: 
 *** Structure definition
 ***
 ***  Author: Nikos Tryfonidis, November 2015
 *** The MIT License (MIT)
 *** Copyright (c) 2015 Nikos Tryfonidis
 *** See LICENSE.txt
 *******************************************************************/

/* Structure for passing and returning x and v 
   (independent and dependent variables), packed nice and pretty.
*/
struct xv{
  double x;
  double v;
};
