//								      -*-c++-*-
// FringeFall.h - class FringeFallFunc, 
//                calculation of multipol fringe field fall off,
//                its Gradients and their derivatives
//
// V1.0, T.Pospischil, 22.01.98
//
// Copyright (c) 1998-2001
//
// Institut für Kernphysik, Universität Mainz	tel. +49 6131 39-25802
// 55099 Mainz, Germany				fax  +49 6131 39-22964
//
// $Id: FringeFall.h 2216 2008-06-13 21:13:47Z distler $
//

#ifndef __FRINGEFALL_H__
#define __FRINGEFALL_H__

class FringeFallFunc {
  
private:
  double S(double s) { 
    return c[0]+s*(c[1]+s*(c[2]+s*(c[3]+s*(c[4]+s*c[5]))));  }
  double dev1(double s) { 
    return (c[1]+s*(2*c[2]+s*(3*c[3]+s*(4*c[4]+s*5*c[5]))));  };
  double dev2(double s) { 
    return 2*c[2]+s*(2*3*c[3]+s*(3*4*c[4]+s*4*5*c[5]));  };
  double dev3(double s) { 
    return (2*3*c[3]+s*(2*3*4*c[4]+s*3*4*5*c[5])); };
  double dev4(double s) { 
    return 2*3*4*c[4]+s*2*3*4*5*c[5]; };
  double dev5(double s) { 
    return (2*3*4*5*c[5] + 0.0*s); };
  
public:
  double R;     // Multipole specific fringe field fall off radius
  double dz;    // Multipole specific inward displacement of EFB 
  double c[6];  // Multipole specific fall off coeffincients (=C0n/(2*R))

  // Constructor:
  FringeFallFunc() { R=1; dz=0; for(int i=0; i<=5; i++) c[i]=0; }
  
  // Functions:
  void init(double CX[6],    // Fall off Coefficients, as given in RAYTRACE
	    double radius,   // Multipole specific fall off radius
	    double delta_z); // Multipole specific inward displacement

  double d0(double z);
  double d1(double z);
  double d2(double z);
  double d3(double z);
  double d4(double z);
  double d5(double z);
  double d6(double z);
};

#endif
