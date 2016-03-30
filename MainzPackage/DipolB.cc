//
// $Id: DipolB.cc 2348 2009-10-09 10:07:29Z merkel $
//

// implementation of the magnet class DipolB for Spectrometer B,
// based on code by T.Pospischil, 22.01.98 for SpectrometerA
// but with much simpler parameterization of the magnetic field



//********************************************************************
//     This is fixed dipole for Spec B
//********************************************************************


#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "DipolB.h"

// Normal constructor, everything is set to some initial value
DipolB::DipolB(){

  int i;

  B_Nom = 1.5;
  B_Act = B_Nom;
  R = 1000;
  D = 100;
  En = Vector3D(0.0, 0.0, 0.0);
  En_phi = 0.0;
  En_xmin = -1000;
  En_xmax = 1000;
  En_zmin = -1000;
  En_zmax = 1000;
  Ex = Vector3D(0.0, 0.0, 0.0);
  Ex_phi = 0.0;
  Ex_xmin = -1000;
  Ex_xmax = 1000;
  Ex_zmin = -1000;
  Ex_zmax = 1000;
  for(i=0; i<=8; i++) { S0[i] = 0.0; S1[i] = 0.0; }
  for(i=0; i<=5; i++) { C0[i] = 0.0; C1[i] = 0.0; }

  E_y = Vector3D(0.0, 1.0, 0.0);       // Einheitsvektor in y-Richtung
}



int DipolB::BField(const double x0[], double B[]) {
  const double tilt = 3.495000 * M_PI/180;            // pole piece tilt angle
  const double gap = 20;                              // gap width at 1.5T line
  const double entrance = tan((90-55+51.8)*M_PI/180); // entrance angle
  const double xe = 450;                              // entrance point
  const double T15Line = 55 * M_PI/180;               // angle of 1.5T line
  const double s55=sin(T15Line), c55=cos(T15Line);
  const double D = gap/tilt;
  //  const double T15length = 350, x1r=320;
  double x[3];
	
  x[0] = -x0[0]; // Rotating coordinate system
  x[1] = -x0[1]; // Rotating coordinate system
  x[2] = x0[2];  

  B[0] = B[1] = B[2] = 0;
  if((x[0]-234.5)*(x[0]-234.5)+(x[2]-452.0)*(x[2]-452.0) >297.0*297.0) return 0;
  if((x[0]-213.5)*(x[0]-213.5)+(x[2]-429.5)*(x[2]-429.5) <158.5*158.5) return 0;
  if((x[0]-592.0)*(x[0]-592.0)+(x[2]-559.5)*(x[2]-559.5) <320.0*320.0) return 0;
  if (x[2] < xe - x[0]/entrance) return 0 ;

  double xr[] = {c55*x[0]-s55*(x[2]-xe)-D, x[1], s55*x[0]+c55*(x[2]-xe)};
  double dr   = sqrt(xr[0]*xr[0] + xr[1]*xr[1]);
  double phi  = atan2(xr[1], -xr[0]);
  if (fabs(phi)>tilt/2.0) return 0;
 
  double B2 = B_Act*D/dr;
  double B1 = sin(phi)*B2;
  
  B[0] =  -(c55 * B1)*1.0;
  B[1] =  -(cos(phi)*B2);
  B[2] =  -s55 * B1; 
  
  return 1;
}

 // calculate Entrance boundary dz for a given EnCS x
double 
DipolB::dz_En(double x)
{
  return 0.0;
}

// calculate Exit boundary dz for a given ExCS x
double 
DipolB::dz_Ex(double x)
{
 
  return 0.0;
}

// ----------------------------- inside --------------------------------------

int DipolB::inside(Vector3D X_A)
{
  //  int inval = 0;
  double B[3];
  double X[3];

  // Method BField requires inpud data in cm not im mm. We must also consider transformation X[0]-> -X[0];	
  X[0] = X_A[1]/10.0;
  X[1] = X_A[2]/10.0;
  X[2] = X_A[3]/10.0; 
  return DipolB::BField(X,B);
}


// ---------------------------- B_yMidPlane ----------------------------------
// calculate B_y in the Midplane of the Spectrometer at a postion r
double 
DipolB::B_yMidPlane(int Inside,  // result of inside(r), if <0 value is det.
		   Vector3D r)  // position in the ACS
{

  return 0.0;
}

// --------------------------- GetField ----------------------------------------
// Calculate Magnetic Field at a given position:
Vector3D 
DipolB::GetField(int Inside, Vector3D r)
{
  Vector3D B;
  double X[3];	
  double Bi[3];

  X[0] = r[1]/10.0;
  X[1] = r[2]/10.0;
  X[2] = r[3]/10.0;
 
  DipolB::BField(X,Bi); 

  for(int i = 0; i<3; i++) B[i+1] = Bi[i];
  
  return B;
}

// Set Magnetic Field to fraction (factor Fract) of nominal Field:
int 
DipolB::ScaleField(double Fract)   // i.e. Fract=0.5 sets half of nom. F.
{
  B_Act = Fract*B_Nom;
  return 0;
}


// write data file with boundary curvature:
int 
DipolB::showBoundary(char *FileName, // Name of output data file 
		     int StepNum)    // Number of points on boundary
{
  // open file:
  FILE *outfile;
  if ( NULL == (outfile = fopen(FileName,"w")) ) {
    fprintf(stderr," Cannot open file %s for output!\n",FileName);
    return 2;
  }

  double x, z;
  
  // levi rob
  double entrance = tan((90-55+51.8)*M_PI/180);
  
  for(x = -550.0; x<=630.0; x = x+10.0){

    z = 4500.0 + x/entrance;
    fprintf(outfile,"%8.3f %8.3f\n",z, -1.0*x); 
  }
  
  // Spodnji rob
  for(double phi = 0.01; phi<=M_PI/2.0 + 0.30; phi = phi+0.05){
    x = 2970.0*cos(phi) - 2340.0;
    z = 2970.0*sin(phi) + 4520.0;
    fprintf(outfile,"%8.3f %8.3f\n",z, -1.0*x);
  }
  
  // desni zgorni rob
  for(double phi = 0.58; phi>= 0.055; phi = phi-0.1){

    x = 3200.0*cos(phi) - 5920.0;
    z = 3200.0*sin(phi) + 5595.0;
    fprintf(outfile,"%8.3f %8.3f\n",z, -1.0*x);
  }
  
  // zgornji rob
  for(double phi = M_PI/2.0 + 0.4; phi>=0.071; phi = phi - 0.1){

    x = 1585.0*cos(phi) - 2135.0;
    z = 1585.0*sin(phi) + 4295.0;
    fprintf(outfile,"%8.3f %8.3f\n",z, -1.0*x); 
  }
  
  // sticisce
  x = -550.0;
  z = 4500.0 + -550.0/entrance;
  fprintf(outfile,"%8.3f %8.3f\n",z, -1.0*x);

  fclose(outfile);
  return 0;
}

// ------------------------- showEntranceBoundary -----------------------------
// write data file with entrance boundary curvature:
int 
DipolB::showEntranceBoundary(char *FileName, // Name of output data file 
			     int StepNum)    // Number of points on boundary
{
  return 0;
}

// ------------------------- showExitBoundary -----------------------------
// write data file with exit boundary curvature:
int 
DipolB::showExitBoundary(char *FileName, // Name of output data file 
			 int StepNum)    // Number of points on boundary
{
  return 0;
}

