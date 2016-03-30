//
// $Id: DipolB.h 2216 2008-06-13 21:13:47Z distler $
//

#ifndef __DIPOLB_H__
#define __DIPOLB_H__

#include "Matrix3D.h"
#include "FringeFall.h"
#include "Magnets.h"
// +++++++++++++++++++++++++++++++++ Dipol ++++++++++++++++++++++++++++++++++++

class DipolB : public Magnet {
public:

  double B_Nom; // nominal homogenous B field in Tesla
  double B_Act; // actual homogenous B field in Tesla, B_Act = Fract*B_Nom;
  double R;     // dipole (reference ray) radius (mm)
  double D;     // Gap width (mm)
  double S0[9]; // Coefficients in polynom description of entrance boundary
  double S1[9]; // Coefficients in polynom description of entrance boundary

  // --------------------------------------Magnet virtual functions to implement:

  // check if a ACS position X_A is between the entrance and the exit boundary:
  int inside(Vector3D X_A);
  // Calculate Magnetic Field at a given position:
  Vector3D GetField(int Inside,  // result of inside(r), if unknown set it < 0
		    Vector3D r); // position in the ACS
  // Set Magnetic Field to fraction (factor Fract) of nominal Field:
  int ScaleField(double Fract);   // i.e. Fract=0.5 sets half of nom. F.

  // write data file with boundary curvature:
  int showBoundary(char *FileName, // Name of output data file 
  		   int StepNum);   // Number of points on boundary
  // write data file with entrance boundary curvature:
  int showEntranceBoundary(char *FileName, // Name of output data file 
			   int StepNum);   // Number of points on boundary
  // write data file with exit boundary curvature:
  int showExitBoundary(char *FileName, // Name of output data file 
		       int StepNum);   // Number of points on boundary

  // ----------------- Constructors:
  DipolB();      // Normal construktor, everything is set to some initial value

  		     
private:
  int BField(const double x0[], double B[]);
  // calculate B_y in the Midplane of the Spectrometer at a position r
  double B_yMidPlane(int Inside,  // result of inside(r), if <0 value is det.
		     Vector3D r); // position in the ACS
  double dz_En(double x); // calculate Entrance boundary dz for a given EnCS x
  double dz_Ex(double x); // calculate Exit boundary dz for a given ExCS x

};

#endif

