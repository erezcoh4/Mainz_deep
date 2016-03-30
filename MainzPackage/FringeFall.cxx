//								      -*-c++-*-
// FringeFall.cc - class FringeFallFunc, 
//                 calculation of multipol fringe field fall off,
//                 its Gradients and their derivatives
//
// derivatives determined by H. Merkel with mathematica in Jan. 98
//   
// V1.0, T.Pospischil, 22.01.98
//
// Copyright (c) 1998-2001
//
// Institut für Kernphysik, Universität Mainz	tel. +49 6131 39-25802
// 55099 Mainz, Germany				fax  +49 6131 39-22964
//
// $Id: FringeFall.cc 2216 2008-06-13 21:13:47Z distler $
//

#include "FringeFall.h"
#include <math.h>

// -------------------------- init:
void
FringeFallFunc::init(double CX[6], // Fall off Coeff., as given in RAYTRACE
		     double radius, // Multipole specific fall off radius
		     double delta_z)// Multipole specific inward displacem.
{ 
  R=radius; dz=delta_z; 
  // The coeffincients are recalculated that s=z-dz (RAYTRACE: s=(y-dz)/(2R) )
  for(int i=0; i<=5; i++) c[i] = CX[i]/pow((2.0*R),i);
}
  
double FringeFallFunc::d0(double z) 
{
  double s = z + dz;
  return 1 / (1 + exp(S(s)));
}
double FringeFallFunc::d1(double z) 
{
  double s = z + dz;
  return -((exp(S(s))*dev1(s))/pow(1 + exp(S(s)),2));
}

double FringeFallFunc::d2(double z)
{ 
  double s = z + dz;
  return (exp(S(s))*((-1 + exp(S(s)))*pow(dev1(s),2) -
		     (1 + exp(S(s)))*dev2(s)))/pow(1 + exp(S(s)),3);
}

double FringeFallFunc::d3(double z) 
{
  double s = z + dz;
  return -((exp(S(s))*((1 - 4*exp(S(s)) + exp(2*S(s)))*
		       pow(dev1(s),3) - 
		       3*(-1 + exp(2*S(s)))*dev1(s)* dev2(s) + 
		       pow(1 + exp(S(s)),2)*dev3(s)))/pow(1 + exp(S(s)),4));
}

double FringeFallFunc::d4(double z) 
{ 
  double s = z + dz;
  return (exp(S(s))*((-1 + 11*exp(S(s)) - 11*exp(2*S(s)) + 
		      exp(3*S(s)))*pow(dev1(s),4) - 
		     6*(1 - 3*exp(S(s)) - 3*exp(2*S(s)) + exp(3*S(s)))*
		     pow(dev1(s),2)*dev2(s) + 
		     4*(-1 + exp(S(s)))*pow(1 + exp(S(s)),2)*
		     dev1(s)*dev3(s) + 
		     pow(1 + exp(S(s)),2)*
		     (3*(-1 + exp(S(s)))*pow(dev2(s),2)
		      - (1 + exp(S(s)))*dev4(s))))/pow(1 + exp(S(s)),5);
}

double FringeFallFunc::d5(double z) 
{
  double s = z + dz;
  return -((exp(S(s))*
	    ((1 - 26*exp(S(s)) + 66*exp(2*S(s)) - 
	      26*exp(3*S(s)) + exp(4*S(s))) * pow(dev1(s),5) - 
	     10*(-1 + 10*exp(S(s))-10*exp(3*S(s))+exp(4*S(s)))*pow(dev1(s),3)*
	     dev2(s) + 10*pow(1 + exp(S(s)),2)*(1 - 4*exp(S(s)) + exp(2*S(s)))*
	     pow(dev1(s),2)*dev3(s) + 5*pow(1 + exp(S(s)),2)* dev1(s)*
	     (3*(1 - 4*exp(S(s)) + exp(2*S(s)))*
	      pow(dev2(s),2) - (-1 + exp(2*S(s)))*dev4(s)) - 
	     pow(1 + exp(S(s)),3)*
	     (10*(-1 + exp(S(s)))*dev2(s) * dev3(s) -
	      (1 + exp(S(s)))*dev5(s))))/pow(1 + exp(S(s)),6));
}
    
double FringeFallFunc::d6(double z) 
{
  double s = z + dz;
  return (exp(S(s)) * 
	  ((-1 + 57*exp(S(s)) - 302*exp(2*S(s)) + 
	    302*exp(3*S(s)) - 57*exp(4*S(s)) + exp(5*S(s))) * pow(dev1(s),6) - 
	   15*(1 - 25*exp(S(s)) + 40*exp(2*S(s)) + 
	       40*exp(3*S(s)) - 25*exp(4*S(s)) + exp(5*S(s)))*pow(dev1(s),4)*
	   dev2(s) + 20*pow(1 + exp(S(s)), 2)
	   *(-1 + 11*exp(S(s)) - 11*exp(2*S(s)) + 
	     exp(3*S(s)))*pow(dev1(s),3)*
	   dev3(s) + 15 * pow(1 + exp(S(s)),2) * pow(dev1(s),2)*
	   (3*(-1 + 11*exp(S(s)) - 11*exp(2*S(s)) + 
	       exp(3*S(s)))*pow(dev2(s),2) -
	    (1 - 3*exp(S(s)) - 3*exp(2*S(s)) + 
	     exp(3*S(s)))*dev4(s)) - 6*pow(1 + exp(S(s)),3) * dev1(s) 
	   *(10*(1 - 4*exp(S(s)) + exp(2*S(s)))*
	     dev2(s)*dev3(s) - (-1 + exp(2*S(s)))*dev5(s)) - 
	   pow(1 + exp(S(s)),3)* 
	   (15*(1 - 4*exp(S(s)) + exp(2*S(s)))*
	    pow(dev2(s),3) -	15*(-1 + exp(2*S(s)))*dev2(s)*
	    dev4(s) - (1 + exp(S(s)))* 
	    (-10*pow(dev3(s),2) + 
	     10*exp(S(s))*pow(dev3(s),2)))))/pow(1 + exp(S(s)),7);
}










