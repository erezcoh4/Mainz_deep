#ifndef CK_RUNGEKUTTA_CXX
#define CK_RUNGEKUTTA_CXX

#include "CK_RungeKutta.h"
// C++ file CK_RungeKutta.cc
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

using namespace std;

#define SIGN(a,b) ((b) >= 0 ? fabs(a) : -fabs(a))

static double dmaxarg1, dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))

static double dminarg1, dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))


/* ------- rkck: Take a Cash-Karp Runge-Kutta step:
 
 Given values for n variables y[1..n] and their derivatives dydx[1..n]
 known at x, use the fifth-order Cash-Karp Runge-Kutta method to advance
 the solution over an interval h and return the incremented variables as
 yout[1..n]. Also return an estimate of the local truncation error in
 yout using the embedded fourth-order method.
 The user supplies the routine derivs(x,y,dydx), which returns dervatives
 dydx at x.
 */

void CK_RungeKutta::rkck(double y[],
          double dydx[],
          int n,
          double x,
          double h,
          double yout[],
          double yerr[],
          void (*derivs)(double, double [], double []) )
{
    int i;
    
    // Carsh-Karp Parameters for embedded Runge-Kutta Method:
    static double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875, b21=0.2,
    b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42=-0.9, b43=1.2,
    b51=-11.0/54.0, b52=2.5, b53=-70.0/27.0, b54=35.0/27.0,
    b61=1631.0/55296.0, b62=175.0/512.0, b63=575.0/13824.0,
    b64=44275.0/110592.0, b65=253.0/4096.0, c1=37.0/378.0,
    c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
    double dc1=c1-2825.0/27648.0, dc3=c3-18578.0/48384.0,
    dc4=c4-13525.0/55296.0, dc6=c6-0.25;
    
    // alloc temp vectors:
    double *ak2 = new double[n+1];
    double *ak3 = new double[n+1];
    double *ak4 = new double[n+1];
    double *ak5 = new double[n+1];
    double *ak6 = new double[n+1];
    double *ytemp = new double[n+1];
    
    // ------------------- Fifth order Runge-Kutta (with embedded fourth order):
    // First + second step:
    for (i=1; i<=n; i++)
        ytemp[i] = y[i] + b21*h*dydx[i];
    (*derivs)(x+a2*h,ytemp,ak2);
    
    // Third step:
    for (i=1; i<=n; i++)
        ytemp[i] = y[i] + h*(b31*dydx[i]+b32*ak2[i]);
    (*derivs)(x+a3*h,ytemp,ak3);
    
    // Fourth step:
    for (i=1; i<=n; i++)
        ytemp[i] = y[i] + h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
    (*derivs)(x+a4*h,ytemp,ak4);
    
    // Fifth step:
    for (i=1; i<=n; i++)
        ytemp[i] = y[i] + h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
    (*derivs)(x+a5*h,ytemp,ak5);
    
    // Sixth step:
    for (i=1; i<=n; i++)
        ytemp[i] = y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
    (*derivs)(x+a6*h,ytemp,ak6);
    
    // Accumulate increments with proper weights:
    for (i=1; i<=n; i++)
        yout[i] = y[i] + h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
    // Estimate error as difference between fourth and fifth order methods:
    for (i=1; i<=n; i++)
        yerr[i] = h * (dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
    
    delete[] ytemp;
    delete[] ak6;
    delete[] ak5;
    delete[] ak4;
    delete[] ak3;
    delete[] ak2;
    
}


#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
// ERRCON equals (5/SAFETY) raised to the power (1/PGROW), see use below.

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

void CK_RungeKutta::rkqs(double y[],
          double dydx[],
          int n,
          double *x,
          double htry,
          double eps,
          double yscal[],
          double *hdid,
          double *hnext,
          void (*derivs)(double, double [], double []))
{
    
    // void rkck(double y[], double dydx[], int n, double x, double h,
    //	    double yout[], double yerr[],
    //	    void (*derivs)(double, double [], double []));
    
    int i;
    double errmax,h,htemp,xnew;
    double *yerr  = new double[n+1];
    double *ytemp = new double[n+1];
    
    // Set stepsize to the initial trial value:
    h=htry;
    
    for (;;){
        
        rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);  // Take a step.
        
        // Evaluate accuracy:
        errmax = 0.0;
        for (i=1; i<=n; i++) errmax = DMAX(errmax, fabs(yerr[i]/yscal[i]));
        errmax /= eps;	        // Scale relative to required tolerance.
        
        // Step succeeded -> Compute size of next step:
        if (errmax <= 1.0) break;
        
        // Truncation error too large -> reduce stepsize:
        htemp = SAFETY * h * pow(errmax,PSHRNK);
        h = (h >= 0.0 ? DMAX(htemp, 0.1*h) : DMIN(htemp, 0.1*h));
        // Not more than a factor of 10.
        
        xnew = (*x) + h;
        if ( xnew == *x )
        {
            //Israel
            //cerr << "ERROR: stepsize underflow (" << h;
            //cerr << ") in rkqs" << endl;
            //exit(1);
            delete[] ytemp;
            delete[] yerr;
            return;
        }
    }
    
    // Compute size of next step:
    if (errmax > ERRCON) *hnext = SAFETY * h * pow(errmax,PGROW);
    else *hnext = 5.0*h;	              // No more than a factor of 5 increase.
    
    *x += (*hdid=h);
    for (i=1; i<=n; i++) y[i]=ytemp[i];
    
    delete[] ytemp;
    delete[] yerr;
}

#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON


#define MAXSTP 500000
#define TINY 0.001

int CK_RungeKutta::odeint(double ystart[],  // Start values
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
           double *xact )    // actual value for the independent variable
{
    /* Runge-Kutta driver with adaptive stepsize control. Integrate starting
     values ystart [1..nvar] from x1 to x2 with accuracy eps, storing inter-
     mediate results in global variables. h1 should be set as a guessed first
     stepsize, hmin as the minimum allowed stepsize (can be zero).
     On output nok and nbad are the number of good and bad (but retried and
     fixed) steps taken, and ystart is replaced by values at the end of the
     integration interval. derivs is the user-supplied routine for calculating
     the right-hand side derivative, while rkqs is the name of the stepper
     routine to be used. */
    
    int nstp,i;
    double xsav=0, x, hnext, hdid, h;
    
    double *yscal = new double[nvar+1];
    double *y     = new double[nvar+1];
    double *dydx  = new double[nvar+1];
    
    x = x1;
    h = SIGN(h1,x2-x1);
    *nok = (*nbad) = kount = 0;
    for (i=1; i<=nvar; i++) y[i] = ystart[i];
    if (kmax > 0) xsav=x-dxsav*2.0;	// Assures storage of first step.
    
    for (nstp=1; nstp<=MAXSTP; nstp++) {  // Take at most MAXSTP steps.
        //  for (nstp=1; nstp<=500000; nstp++) {  // Take at most MAXSTP steps.
        (*derivs)(x,y,dydx);
        for(i=1; i<=nvar; i++)              // Scaling used to monitor accuracy.
            yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + TINY;
        // This general-purpose choice can be modified if need be.
        
        // Store intermediate results:
        if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
            (*stepstore)(x, y, nvar, ++kount);
            xsav=x;
        }
        // lf stepsize can overshoot, decrease:
        if ((x+h-x2)*(x+h-x1) > 0.0) h = x2-x;
        
        rkqs(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext, derivs);
        if (hdid == h) ++(*nok); else ++(*nbad);
        
        // Are we done?
        if ( ((x-x2)*(x2-x1) >= 0.0) ) {
            for (i=1; i<=nvar; i++) ystart[i] = y[i];
            // Save final step:
            if (kmax) (*stepstore)(x, y, nvar, ++kount);
            
            delete[] dydx;
            delete[] y;
            delete[] yscal;
            *xact = x;
            return 1;	// Normal exit.
        }
        // Are we done?
        if ( (y[1] < y1min) ) {
            for (i=1; i<=nvar; i++) ystart[i] = y[i];
            // Save final step:
            if (kmax) (*stepstore)(x, y, nvar, ++kount);
            
            delete[] dydx;
            delete[] y;
            delete[] yscal;
            *xact = x;
            return 2;	// Normal exit.
        }
        
        if (fabs(hnext) <= hmin) {
//            cerr << "Step size (" << hnext << ") too small in odeint" << endl;
            return 0;
        }
        h = hnext;
    }
    
    delete[] dydx;
    delete[] y;
    delete[] yscal;
    
    //cerr << "Too many steps in routine odeint" << endl;
    return -1;
}

#undef MAXSTP
#undef TINY

#define MAXSTP 500000
#define TINY 0.001

int CK_RungeKutta::odeint2(double ystart[],  // Start values
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
            int (*checkPosition)(double []) )
{
    /* Runge-Kutta driver with adaptive stepsize control. Integrate starting
     values ystart [1..nvar] from x1 to x2 with accuracy eps, storing inter-
     mediate results in global variables. h1 should be set as a guessed first
     stepsize, hmin as the minimum allowed stepsize (can be zero).
     On output nok and nbad are the number of good and bad (but retried and
     fixed) steps taken, and ystart is replaced by values at the end of the
     integration interval. derivs is the user-supplied routine for calculating
     the right-hand side derivative, while rkqs is the name of the stepper
     routine to be used. */
    
    int nstp,i;
    double xsav=0, x, hnext, hdid, h;
    
    double *yscal = new double[nvar+1];
    double *y     = new double[nvar+1];
    double *dydx  = new double[nvar+1];
    
    x = x1;
    h = SIGN(h1,x2-x1);
    *nok = (*nbad) = kount = 0;
    for (i=1; i<=nvar; i++) y[i] = ystart[i];
    if (kmax > 0) xsav=x-dxsav*2.0;	// Assures storage of first step.
    
    for (nstp=1; nstp<=MAXSTP; nstp++) {  // Take at most MAXSTP steps.
        (*derivs)(x,y,dydx);
        for(i=1; i<=nvar; i++)              // Scaling used to monitor accuracy.
            yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + TINY;
        // This general-purpose choice can be modified if need be.
        
        // Store intermediate results:
        if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
            (*stepstore)(x, y, nvar, ++kount);
            xsav=x;
        }
        // lf stepsize can overshoot, decrease:
        if ((x+h-x2)*(x+h-x1) > 0.0) h = x2-x;
        
        rkqs(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext, derivs);
        if (hdid == h) ++(*nok); else ++(*nbad);
        
        // Are we done?
        if ( ((x-x2)*(x2-x1) >= 0.0) ) {
            for (i=1; i<=nvar; i++) ystart[i] = y[i];
            // Save final step:
            if (kmax) (*stepstore)(x, y, nvar, ++kount);
            
            delete[] dydx;
            delete[] y;
            delete[] yscal;
            *xact = x;
            return 1;	// Normal exit.
        }
        // Are we done?
        if ( (*checkPosition)(y) ) {
            for (i=1; i<=nvar; i++) ystart[i] = y[i];
            // Save final step:
            if (kmax) (*stepstore)(x, y, nvar, ++kount);
            
            delete[] dydx;
            delete[] y;
            delete[] yscal;
            *xact = x;
            return 2;	// Normal exit.
        }
        
        if (fabs(hnext) <= hmin) {
            cerr << "Step size (" << hnext << ") too small in odeint" << endl;
            return 0;
        }
        h = hnext;
    }
    //cerr << "Too many steps in routine odeint" << endl;
    delete[] dydx;
    delete[] y;
    delete[] yscal;
    return -1;
}

#undef MAXSTP
#undef TINY

#define MAXSTP 500000
#define TINY 0.001

int CK_RungeKutta::odeint3(double ystart[],  // Start values
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
            int (*checkInside)(double x, double y, double z)  )
{
    /* Runge-Kutta driver with adaptive stepsize control. Integrate starting
     values ystart [1..nvar] from x1 to x2 with accuracy eps, storing inter-
     mediate results in global variables. h1 should be set as a guessed first
     stepsize, hmin as the minimum allowed stepsize (can be zero).
     On output nok and nbad are the number of good and bad (but retried and
     fixed) steps taken, and ystart is replaced by values at the end of the
     integration interval. derivs is the user-supplied routine for calculating
     the right-hand side derivative, while rkqs is the name of the stepper
     routine to be used. */
    
    int nstp,i;
    double xsav=0, x, hnext, hdid, h;
    
    double *yscal = new double[nvar+1];
    double *y     = new double[nvar+1];
    double *dydx  = new double[nvar+1];
    
    x = x1;
    h = SIGN(h1,x2-x1);
    *nok = (*nbad) = kount = 0;
    for (i=1; i<=nvar; i++) y[i] = ystart[i];
    if (kmax > 0) xsav=x-dxsav*2.0;	// Assures storage of first step.
    
    for (nstp=1; nstp<=MAXSTP; nstp++) {  // Take at most MAXSTP steps.
        (*derivs)(x,y,dydx);
        for(i=1; i<=nvar; i++)              // Scaling used to monitor accuracy.
            yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + TINY;
        // This general-purpose choice can be modified if need be.
        
        // Store intermediate results:
        if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
            (*stepstore)(x, y, nvar, ++kount);
            xsav=x;
        }
        // lf stepsize can overshoot, decrease:
        if ((x+h-x2)*(x+h-x1) > 0.0) h = x2-x;
        
        rkqs(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext, derivs);
        if (hdid == h) ++(*nok); else ++(*nbad);
        
        // Are we done?
        if ( ((x-x2)*(x2-x1) >= 0.0) ) {
            for (i=1; i<=nvar; i++) ystart[i] = y[i];
            // Save final step:
            if (kmax) (*stepstore)(x, y, nvar, ++kount);
            
            delete[] dydx; 
            delete[] y;
            delete[] yscal;
            *xact = x;
            return 1;	// Normal exit.
        }
        // Are we done?
        if ( (*checkPosition)(y) ) {
            for (i=1; i<=nvar; i++) ystart[i] = y[i];
            // Save final step:
            if (kmax) (*stepstore)(x, y, nvar, ++kount);
            
            delete[] dydx; 
            delete[] y;
            delete[] yscal;
            *xact = x;
            return 2;	// Normal exit.
        }
        // preverimo, ali smo med integracijo padli ven iz detektroja!!!
        if ( !(*checkInside)(y[1], y[2], y[3]) ) {
            for (i=1; i<=nvar; i++) ystart[i] = y[i];
            // Save final step:
            if (kmax) (*stepstore)(x, y, nvar, ++kount);
            delete[] dydx; 
            delete[] y;
            delete[] yscal;
            *xact = x;
            return 3;	// Normal exit.
        }
        
        if (fabs(hnext) <= hmin) { 
            cerr << "Step size (" << hnext << ") too small in odeint" << endl;
            return 0;
        }
        h = hnext;
    }
    //cerr << "Too many steps in routine odeint" << endl;
    delete[] dydx; 
    delete[] y;
    delete[] yscal;
    return -1;
}

#undef MAXSTP
#undef TINY


#endif
