// C++ file CK_RungeKutta.h
//
// Cash-Karp Runge-Kutta solution of ordinary diffrential equations
// (call routine odeint)
//
// The original functions are taken from the book "NUMERICAL RECIPES in C".
//  
//  Version 1.0, T.Pospischil, 12/97: 
//   Modification of the original routines:  
//    - memory handling (nrutil routines replaced by "new" statements)
//    - step storage is written into a file instead of an array.
//    - value of TINY raised from 1e-30 to 0.001.
//

#ifndef __CK_RungeKutta_h__
#define __CK_RungeKutta_h__

/* ------- rkck: Take a Cash-Karp Runge-Kutta step:

   Given values for n variables y[1..n] and their derivatives dydx[1..n] 
   known at x, use the fifth-order Cash-Karp Runge-Kutta method to advance 
   the solution over an interval h and return the incremented variables as 
   yout[1..n]. Also return an estimate of the local truncation error in 
   yout using the embedded fourth-order method. 
   The user supplies the routine derivs(x,y,dydx), which returns dervatives 
   dydx at x.
*/

void rkck(double y[], 
	  double dydx[],
	  int n, 
	  double x,
	  double h,
	  double yout[],
	  double yerr[],
	  void (*derivs)(double, double [], double []) );


/* ----------------------- rkqs: Fifth-order Runge-Kutta stepper

  Fifth-order Runge-Kutta step with monitoring of local truncation error 
  to ensure accuracy and adjust stepsize. Input are the dependent variable 
  vector y[1..n] and jts derivative dydx[1..n] at the starting value of the 
  independent variable x. Also input are the stepsize to be attempted htry, 
  the required accuracy eps, and the vector yscal [1..n] against which the 
  error is scaled. 
  On output, y and x are replaced by their new values, hdid is the stepsize 
  that was actually accomplished, and hnext is the estimated next stepsize. 
  derivs is the user-supplied routine that computes the right-hand side 
  derivatives.
*/

void rkqs(double y[], 
	  double dydx[],
	  int n, 
	  double *x, 
	  double htry, 
	  double eps,
	  double yscal[], 
	  double *hdid,
	  double *hnext,
	  void (*derivs)(double, double [], double []));


/* ---------------------------- odeint -------------------------------------
   Runge-Kutta driver with adaptive stepsize control. Integrate starting 
   values ystart [1..nvar] from x1 to x2 with accuracy eps, storing inter-
   mediate results in global variables. h1 should be set as a guessed first 
   stepsize, hmin as the minimum allowed stepsize (can be zero). 
   On output nok and nbad are the number of good and bad (but retried and 
   fixed) steps taken, and ystart is replaced by values at the end of the 
   integration interval. derivs is the user-supplied routine for calculating
   the right-hand side derivative, while rkqs is the name of the stepper
   routine to be used. */

int odeint(double ystart[],  // Start values
	    int    nvar,      // number of ODEs (= dimension of the arrays)
	    double x1,        // start value for the independent variable
	    double x2,        // stop at this value of the independent variable
	    double y1min,     // stop if y[1] smaller than y1min
	    double eps,       // required accuracy
	    double h1,        // guessed first stepsize
	    double hmin,      // min. allowed stepsize (can be zero)
	    int    *nok,      // number of good steps taken
	    int    *nbad,     // number of bad (retried and fixed) steps taken
	    double dxsav,     // interval size for step storing
	    int    kmax,      // max number of steps to store
	    int    &kount,    // number of stored steps
	    void (*stepstore)(double, double [], int, int),// func stores steps
	    void (*derivs)(double, double [], double []),  // derivatives
	    double *xact );   // actual value for the independent variable

int odeint2(double ystart[],  // Start values
	    int    nvar,      // number of ODEs (= dimension of the arrays)
	    double x1,        // start value for the independent variable
	    double x2,        // stop at this value of the independent variable
	    double y1min,     // stop if y[1] smaller than y1min
	    double eps,       // required accuracy
	    double h1,        // guessed first stepsize
	    double hmin,      // min. allowed stepsize (can be zero)
	    int    *nok,      // number of good steps taken
	    int    *nbad,     // number of bad (retried and fixed) steps taken
	    double dxsav,     // interval size for step storing
	    int    kmax,      // max number of steps to store
	    int    &kount,    // number of stored steps
	    void (*stepstore)(double, double [], int, int),// func stores steps
	    void (*derivs)(double, double [], double []),  // derivatives
	    double *xact,	// actual value for the independent variable
	    int (*checkPosition)(double []) ) ; 

int odeint3(double ystart[],  // Start values
	    int    nvar,      // number of ODEs (= dimension of the arrays)
	    double x1,        // start value for the independent variable
	    double x2,        // stop at this value of the independent variable
	    double y1min,     // stop if y[1] smaller than y1min
	    double eps,       // required accuracy
	    double h1,        // guessed first stepsize
	    double hmin,      // min. allowed stepsize (can be zero)
	    int    *nok,      // number of good steps taken
	    int    *nbad,     // number of bad (retried and fixed) steps taken
	    double dxsav,     // interval size for step storing
	    int    kmax,      // max number of steps to store
	    int    &kount,    // number of stored steps
	    void (*stepstore)(double, double [], int, int),// func stores steps
	    void (*derivs)(double, double [], double []),  // derivatives
	    double *xact,	// actual value for the independent variable
	    int (*checkPosition)(double []),
	    int (*checkInside)(double x, double y, double z)) ; 


#endif

