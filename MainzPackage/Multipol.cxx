//   Dipol.cc - implementation of the magnet class Multipol, 
//              field calculation decribed in the RAYTRACE manual
//              WARNING: may be not accurate for very large quads, ...
//
//   V1.0, T.Pospischil, 22.01.98
//
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "Magnets.h"

// Normal constructor, everything is set to some initial value
Multipol::Multipol(){

  int i;

  R = 1000;
  L = 100;
  for(i=0;i<=5;i++) { B_Nom[i] = 0.0; G[i]=0.0; }

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
  for(i=0; i<=5; i++) { C0[i] = 0.0; C1[i] = 0.0; }

  E_y = Vector3D(0.0, 1.0, 0.0);       // Einheitsvektor in y-Richtung
}

// Constructor which initialises everything:
//  EnCS and ExCS are assumed to be antiparallel with the same width
//  (xmin to xmax) and a distance L between the origins
Multipol::Multipol(double iB_Nom[6],// Multipole field strength (Tesla)
		   double iR,       // Radius of aperture (mm)
		   double iL,       // Effective length of the multipole (mm)
		   Vector3D iMCS,   // Origin of the MCS in the ACS (mm)
		   double iphi,     // Angle of MCS in ACS zx-plane (deg)
		   double ixmin,    // Effective x-range min of the multipole
		   double ixmax,    // Effective x-range max of the multipole
		   double iC0[6],   // Coeff. of entrance fringe field fall off
		   double iEn_zmin, // Min. of z-Range of Entrance Boundary
		   double iEn_zmax, // Max. of z-Range of Entrance Boundary
		   double iC1[6],   // Coeff. of exit fringe field fall off
		   double iEx_zmin, // Min. of z-Range of Exit Boundary
		   double iEx_zmax, // Max. of z-Range of Exit Boundary
		   double FRH,      // ff fall off hexapole radius rel. to R
		   double FRO,      // ff fall off octapole radius rel. to R
		   double FRD,      // ff fall off decapole radius rel. to R
		   double FRDD,     // ff fall off dodecapole rad. rel. to R
		   double DSH,      // hexa EFB inward displacem (DSH = dz/2R) 
		   double DSO,      // rel. octapole EFB inward displacement 
		   double DSD,      // rel. decapole EFB inward displacement 
		   double DSDD)     // rel. dodecapole EFB inward displacement
{
  int i;

  R = iR;
  L = iL;

  B_Nom[0] = 0.0;

  if ( (B_Nom[1]=iB_Nom[1]) != 0.0 ) Quad  = 1; else Quad  = 0;
  EmQ.init(iC0,-R,0.0);  
  ExQ.init(iC1,R,0.0);
  if ( (B_Nom[2]=iB_Nom[2]) != 0.0 ) Hexa  = 1; else Hexa  = 0;
  EmH.init(iC0,(-R*FRH),(DSH*2.0*(-R*FRH)));
  ExH.init(iC1,(R*FRH),(DSH*2.0*(R*FRH)));
  if ( (B_Nom[3]=iB_Nom[3]) != 0.0 ) Octa  = 1; else Octa  = 0;
  EmO.init(iC0,(-R*FRO),(DSO*2.0*(-R*FRO)));
  ExO.init(iC1,(R*FRO),(DSO*2.0*(R*FRO)));
  if ( (B_Nom[4]=iB_Nom[4]) != 0.0 ) Deca  = 1; else Deca  = 0;
  EmD.init(iC0,(-R*FRD),(DSD*2.0*(-R*FRD)));
  ExD.init(iC1,(R*FRD),(DSD*2.0*(R*FRD)));
  if ( (B_Nom[5]=iB_Nom[5]) != 0.0 ) Dodeca = 1; else Dodeca = 0;
  EmDD.init(iC0,(-R*FRDD),(DSDD*2.0*(-R*FRDD)));
  ExDD.init(iC1,(R*FRDD),(DSDD*2.0*(R*FRDD)));

  // Set Gradients:
  for(i=1; i<=5; i++) G[i] = B_Nom[i]/pow(R,i);

  En = iMCS;
  En_phi = iphi + 180.0;
  En_xmin = -1.0*ixmax;
  En_xmax = -1.0*ixmin;
  En_zmin = iEn_zmin;
  En_zmax = iEn_zmax;
  E_y = Vector3D(0.0, 1.0, 0.0);       // Einheitsvektor in y-Richtung
  rot_A_to_En = Matrix3D(E_y, En_phi*deg);
  rot_En_to_A = Matrix3D(E_y, -1.0*En_phi*deg);

  // origin of ExCS is at (0, 0, -L) in the En_CS: almost -mod-
  Ex[1] = En[1]-sin(iphi*deg)*L; Ex[2] = En[2]; Ex[3] = En[3]+cos(iphi*deg)*L;

  Ex_phi = iphi;
  Ex_xmin = ixmin;
  Ex_xmax = ixmax;
  Ex_zmin = iEx_zmin;
  Ex_zmax = iEx_zmax;
  rot_A_to_Ex = Matrix3D(E_y, Ex_phi*deg);
  rot_Ex_to_A = Matrix3D(E_y, -1.0*Ex_phi*deg);
  // cout << rot_Ex_to_A << endl; exit(1);

  for(i=0; i<=5; i++) { C0[i] = iC0[i]; C1[i] = iC1[i]; }
 
  bEn_xmin = bfwef*En_xmin;
  bEn_xmax = bfwef*En_xmax;
  bEx_xmin = bfwef*Ex_xmin;
  bEx_xmax = bfwef*Ex_xmax;
   
}

// --------------------------- showBoundary ----------------------------------
// write data file with boundary curvature:
int 
Multipol::showBoundary(char *FileName, // Name of output data file 
                       int StepNum)    // Number of points on boundary
{
  // open file:
  FILE *outfile;
  if ( NULL == (outfile = fopen(FileName,"w")) ) {
    fprintf(stderr," Cannot open file %s for output!\n",FileName);
    return 2;
  }

  Vector3D X_E, X_A;
  X_E[1] = StepNum; // dummy, but supresses warning

  X_E[1]=En_xmin;
  X_E[3]=0.0;
  X_A = rot_En_to_A*X_E + En;  // trafo to ACS 
  fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
  X_E[1]=En_xmax;
  X_E[3] = 0.0;
  X_A = rot_En_to_A*X_E + En;  // trafo to ACS 
  fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
  X_E[1]=Ex_xmin;
  X_E[3]=0.0;
  X_A = rot_Ex_to_A*X_E + Ex;  // trafo to ACS 
  fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
  X_E[1]=Ex_xmax;
  X_E[3] = 0.0;
  X_A = rot_Ex_to_A*X_E + Ex;  // trafo to ACS 
  fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
  X_E[1]=En_xmin;
  X_E[3]=0.0;
  X_A = rot_En_to_A*X_E + En;  // trafo to ACS 
  fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);


  fclose(outfile);
  return 0;
}

// ------------------------- showEntranceBoundary -----------------------------
// write data file with entrance boundary curvature:
int 
Multipol::showEntranceBoundary(char *FileName, // Name of output data file 
			       int StepNum)    // Number of points on boundary
{
  // open file:
  FILE *outfile;
  if ( NULL == (outfile = fopen(FileName,"w")) ) {
    fprintf(stderr," Cannot open file %s for output!\n",FileName);
    return 2;
  }

  Vector3D X_E, X_A;
  X_E[1] = StepNum; // dummy, but supresses warning

  X_E[1]=bEn_xmin;
  X_E[3]=En_zmin;
  X_A = rot_En_to_A*X_E + En;  // trafo to ACS 
  fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
  X_E[1]=bEn_xmax;
  X_E[3] =En_zmin;
  X_A = rot_En_to_A*X_E + En;  // trafo to ACS 
  fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
  X_E[1]=bEn_xmax;
  X_E[3]=En_zmax;
  X_A = rot_En_to_A*X_E + En;  // trafo to ACS 
  fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
  X_E[1]=bEn_xmin;
  X_E[3] =En_zmax;
  X_A = rot_En_to_A*X_E + En;  // trafo to ACS 
  fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
  X_E[1]=bEn_xmin;
  X_E[3]=En_zmin;
  X_A = rot_En_to_A*X_E + En;  // trafo to ACS 
  fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);

  fclose(outfile);
  return 0;
}

// ------------------------- showExitBoundary -----------------------------
// write data file with exit boundary curvature:
int 
Multipol::showExitBoundary(char *FileName, // Name of output data file 
			   int StepNum)    // Number of points on boundary
{
  // open file:
  FILE *outfile;
  if ( NULL == (outfile = fopen(FileName,"w")) ) {
    fprintf(stderr," Cannot open file %s for output!\n",FileName);
    return 2;
  }

  Vector3D X_E, X_A;
  X_E[1] = StepNum; // dummy, but supresses warning

  X_E[1]=bEx_xmin;
  X_E[3]=Ex_zmin;
  X_A = rot_Ex_to_A*X_E + Ex;  // trafo to ACS 
  fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
  X_E[1]=bEx_xmax;
  X_E[3] =Ex_zmin;
  X_A = rot_Ex_to_A*X_E + Ex;  // trafo to ACS 
  fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
  X_E[1]=bEx_xmax;
  X_E[3]=Ex_zmax;
  X_A = rot_Ex_to_A*X_E + Ex;  // trafo to ACS 
  fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
  X_E[1]=bEx_xmin;
  X_E[3] =Ex_zmax;
  X_A = rot_Ex_to_A*X_E + Ex;  // trafo to ACS 
  fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
  X_E[1]=bEx_xmin;
  X_E[3]=Ex_zmin;
  X_A = rot_Ex_to_A*X_E + Ex;  // trafo to ACS 
  fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);

  fclose(outfile);
  return 0;
}


// ----------------------------- inside --------------------------------------
// check if a ACS position is between the entrance and the exit boundary
//  possible return values of inside() consists of 3 bits:
//   const int InHomogenField  = 1;
//   const int InEntranceField = 2;
//   const int InExitField     = 4;
int
Multipol::inside(Vector3D X_A)
{
  int inval = 0;
  Vector3D X_En, X_Ex;

  X_Ex = rot_A_to_Ex * (X_A - Ex);   // transform X_A to ExCS

  // transform X_A to EnCS  
  X_En[1] = -1.0*X_Ex[1]; X_En[2] = X_Ex[2]; X_En[3] = -1.0*(X_Ex[3]+L);   

  if ((X_En[1] < bEn_xmax) && (X_En[1] > bEn_xmin)){    // En side enough
    // Are we inside the homogenous field?
    if ( (X_En[3] < 0) && (X_Ex[3] < 0) ) inval += Magnet::HomogenField;
    // Are we inside the entrance field?
    if ( (X_En[3] > En_zmin) && (X_En[3]<En_zmax) ) inval += Magnet::EntranceField;
    // Are we inside the exit field?
    if ( (X_Ex[3] > Ex_zmin) && (X_Ex[3]<Ex_zmax) ) inval += Magnet::ExitField;
  }

  return inval;
}

// ------------------------- GetQuadFringeField ------------------------------
Vector3D
Multipol::GetQuadFringeField(Vector3D r,         // Position in fringe field
			     FringeFallFunc FF)  // fringe field description
{
  Vector3D B;

  double x  = r[1];     double y  = r[2];     double G1_0 = G[1]*FF.d0(r[3]); 
  double x2 = x*x;      double y2 = y*y;      double G1_1 = G[1]*FF.d1(r[3]); 
  double x3 = x*x2;     double y3 = y*y2;     double G1_2 = G[1]*FF.d2(r[3]); 
  double x4 = x*x3;     double y4 = y*y3;     double G1_3 = G[1]*FF.d3(r[3]); 
  double x5 = x*x4;     double y5 = y*y4;     double G1_4 = G[1]*FF.d4(r[3]); 
  double x6 = x*x5;     double y6 = y*y5;     double G1_5 = G[1]*FF.d5(r[3]); 
  double x7 = x*x6;     double y7 = y*y6;     double G1_6 = G[1]*FF.d6(r[3]); 
  
  B[1] = G1_0*y  -  1.0/12.0 * G1_2 * (3*x2*y + y3) 
    +  1.0/384.0   * G1_4 * (5*x4*y + 6*x2*y3 + y5)
    -  1.0/23040.0 * G1_6 * (7*x6*y + 15*x4*y3 + 9*x2*y5 + y7);

  B[2] = G1_0*x  -  1.0/12.0 * G1_2 * (x3 + 3*x*y2) 
    +  1.0/384.0   * G1_4 * (x5 + 6*x3*y2 + 5*x*y4)
    -  1.0/23040.0 * G1_6 * (x7 + 9*x5*y2 + 15*x3*y4 + 7*x*y6);

  B[3] = G1_1*x*y  -  1.0/12.0 * G1_3 * (x3*y + x*y3)
    +  1.0/384.0   * G1_5 * (x5*y + 2*x3*y3 + x*y5);

  return B;
}

// ------------------------- GetHexaFringeField ------------------------------
Vector3D
Multipol::GetHexaFringeField(Vector3D r,         // Position in fringe field
			     FringeFallFunc FF)  // fringe field description
{
  Vector3D B;

  double x  = r[1];     double y  = r[2];     double G2_0 = G[2]*FF.d0(r[3]); 
  double x2 = x*x;      double y2 = y*y;      double G2_1 = G[2]*FF.d1(r[3]); 
  double x3 = x*x2;     double y3 = y*y2;     double G2_2 = G[2]*FF.d2(r[3]); 
  double x4 = x*x3;     double y4 = y*y3;     double G2_3 = G[2]*FF.d3(r[3]); 
                        double y5 = y*y4;                                     
  
  B[1] = G2_0 * (2*x*y)  -  1.0/48.0 * G2_2 * (12*x3*y + 4*x*y3); 
  B[2] = G2_0 * (x2-y2)  -  1.0/48.0 * G2_2 * (3*x4 + 6*x2*y2 - 5*y4); 
  B[3] = G2_1 * (x2*y - y3/3.0)  -  1.0/48.0 * G2_3 * (3*x4*y + 2*x2*y3 - y5);

  return B;
}

// ------------------------- GetOctaFringeField ------------------------------
Vector3D
Multipol::GetOctaFringeField(Vector3D r,         // Position in fringe field
			     FringeFallFunc FF)  // fringe field description
{
  Vector3D B;

  double x  = r[1];     double y  = r[2];     double G3_0 = G[3]*FF.d0(r[3]); 
  double x2 = x*x;      double y2 = y*y;      double G3_1 = G[3]*FF.d1(r[3]); 
  double x3 = x*x2;     double y3 = y*y2;     double G3_2 = G[3]*FF.d2(r[3]); 
  double x4 = x*x3;     double y4 = y*y3;
  double x5 = x*x4;     double y5 = y*y4;                                     
  
  B[1] = G3_0 * (3*x2*y - y3)  -  1.0/80.0 * G3_2 * (20*x4*y - 4*y5); 
  B[2] = G3_0 * (x3 - 3*x*y2)  -  1.0/80.0 * G3_2 * (4*x5 - 20*x*y4); 
  B[3] = G3_1 * (x3*y - x*y3);

  return B;
}

// ------------------------- GetDecaFringeField ------------------------------
Vector3D
Multipol::GetDecaFringeField(Vector3D r,         // Position in fringe field
			     FringeFallFunc FF)  // fringe field description
{
  Vector3D B;

  double x  = r[1];     double y  = r[2];     double G4_0 = G[4]*FF.d0(r[3]); 
  double x2 = x*x;      double y2 = y*y;      double G4_1 = G[4]*FF.d1(r[3]); 
  double x3 = x*x2;     double y3 = y*y2;
  double x4 = x*x3;     double y4 = y*y3;
                        double y5 = y*y4;                                     
  
  B[1] = G4_0 * 4.0*(x3*y - x*y3);
  B[2] = G4_0 * (x4 - 6*x2*y2 + y4);
  B[3] = G4_1 * (x4*y - 2*x2*y3 + y5/5.0);

  return B;
}

// ------------------------- GetDodecaFringeField ------------------------------
Vector3D
Multipol::GetDodecaFringeField(Vector3D r,       // Position in fringe field
			       FringeFallFunc FF)// fringe field description
{
  Vector3D B;

  double x  = r[1];     double y  = r[2];     double G5_0 = G[5]*FF.d0(r[3]); 
  double x2 = x*x;      double y2 = y*y;
  double x3 = x*x2;     double y3 = y*y2;
  double x4 = x*x3;     double y4 = y*y3;
  double x5 = x*x4;     double y5 = y*y4;                                     
  
  B[1] = G5_0 * (5*x4*y - 10*x2*y3 + y5);
  B[2] = G5_0 * (x5 - 10*x3*y2 + 5*x*y4);
  B[3] = 0.0;

  return B;
}

// --------------------------- GetField --------------------------------------
// Calculate Magnetic Field at a given position:
Vector3D 
Multipol::GetField(int Inside, // result of inside(r), if <0 it is determin.
		   Vector3D r) // position in the ACS
{
  Vector3D B(0.0,0.0,0.0);     // Magnetic field

  // Inside value not known -> determine
  if (Inside<0) Inside = inside(r);
  if (!Inside) return B;

  Vector3D X_Em, X_Ex;

  // Rotation to EnCS would be wrong here (<= no field description in EnCS)
  X_Em = rot_A_to_Ex * (r - En);               // transform r to MCS
  X_Ex = rot_A_to_Ex * (r - Ex);               // transform r to ExCS

  // --- Add up the Multipole contributions to the B-Field (use MCS/ExCS):

  // ------------------------------ Quadrupole: ------------------------------
  static Vector3D B_quad;
  B_quad=Vector3D(0.0, 0.0, 0.0);
  if (Quad) {
    // Calculate Uniform field:
    Vector3D B_quad_Uni( (G[1]*X_Ex[2]), (G[1]*X_Ex[1]), 0.0); 

    if (Inside & Magnet::EntranceField) {
      B_quad = GetQuadFringeField(X_Em, EmQ);            // Calculate Field
    }
    if (Inside & Magnet::ExitField) {
      B_quad = B_quad + GetQuadFringeField(X_Ex, ExQ);   // Calculate Field
    }
    if ((Inside & Magnet::ExitField) && (Inside & Magnet::EntranceField)) {
      B_quad = B_quad - B_quad_Uni;
    }
    if (Inside == Magnet::HomogenField) {
      B_quad = B_quad_Uni;
    }
    B = B + B_quad;
  }

  // ------------------------------ Sextupole: --------------------------------
  static Vector3D B_hexa;
  B_hexa=Vector3D(0.0, 0.0, 0.0);
  if (Hexa) {
    // Calculate Uniform field:
    Vector3D B_hexa_Uni( (2.0*G[2]*X_Em[1]*X_Em[2]), 
			 (G[2]*(X_Em[1]*X_Em[1] - X_Em[2]*X_Em[2])), 0.0);
    // Where are we?
    int hexaInside = 0;    // Own Inside cause EFB can be shifted against Quad
    if ( ((-X_Em[3]-EmH.dz) > En_zmin) && ((-X_Em[3]-EmH.dz) < En_zmax) )
      hexaInside += Magnet::EntranceField;
    if ( ((X_Ex[3]+ExH.dz) > Ex_zmin) && ((X_Ex[3]+ExH.dz) < Ex_zmax) )
      hexaInside += Magnet::ExitField;
    if ( ((X_Em[3]+EmH.dz) > 0) && ((X_Ex[3]+ExH.dz) < 0) )
      hexaInside += Magnet::HomogenField;

    if (hexaInside & Magnet::EntranceField) {
      B_hexa = GetHexaFringeField(X_Em, EmH);           // Calculate Field
    }
    if (hexaInside & Magnet::ExitField) {
      B_hexa = B_hexa + GetHexaFringeField(X_Ex, ExH);  // Calculate Field
    }
    if ((hexaInside & Magnet::ExitField) &&
	(hexaInside & Magnet::EntranceField)) {
      B_hexa = B_hexa - B_hexa_Uni;
    }
    if (Inside == Magnet::HomogenField) {
      B_hexa = B_hexa_Uni;
    }
    B = B + B_hexa;
  }

  // ------------------------------ Octapole: --------------------------------
  static Vector3D B_octa;
  B_octa=Vector3D(0.0, 0.0, 0.0);
  if (Octa) {
    // Calculate Uniform field:
    Vector3D B_octa_Uni( (G[3]*(3*X_Em[1]*X_Em[1]*X_Em[2] - 
			       X_Em[2]*X_Em[2]*X_Em[2])), 
			 (G[3]*(X_Em[1]*X_Em[1]*X_Em[1] -
				3*X_Em[1]*X_Em[2]*X_Em[2])), 0.0);
    // Where are we?
    int octaInside = 0;    // Own Inside cause EFB can be shifted against Quad
    if ( ((-X_Em[3]-EmO.dz) > En_zmin) && ((-X_Em[3]-EmO.dz) < En_zmax) )
      octaInside += Magnet::EntranceField;
    if ( ((X_Ex[3]+ExO.dz) > Ex_zmin) && ((X_Ex[3]+ExO.dz) < Ex_zmax) )
      octaInside += Magnet::ExitField;
    if ( ((X_Em[3]+EmO.dz) > 0) && ((X_Ex[3]+ExO.dz) < 0) )
      octaInside += Magnet::HomogenField;

    if (octaInside & Magnet::EntranceField) {
      B_octa = GetOctaFringeField(X_Em, EmO);           // Calculate Field
    }
    if (octaInside & Magnet::ExitField) {
      B_octa = B_octa + GetOctaFringeField(X_Ex, ExO);  // Calculate Field
    }
    if ((octaInside & Magnet::ExitField) &&
	(octaInside & Magnet::EntranceField)) {
      B_octa = B_octa - B_octa_Uni;
    }
    if (Inside == Magnet::HomogenField) {
      B_octa = B_octa_Uni;
    }
    B = B + B_octa;
  }

  // ------------------------------ Decapole: --------------------------------
  static Vector3D B_deca;
  B_deca=Vector3D(0.0, 0.0, 0.0);
  if (Deca) {
    // Calculate Uniform field:
    Vector3D B_deca_Uni( (G[4]*4*(X_Em[1]*X_Em[1]*X_Em[1]*X_Em[2] - 
				  X_Em[1]*X_Em[2]*X_Em[2]*X_Em[2])), 
			 (G[4]*(X_Em[1]*X_Em[1]*X_Em[1]*X_Em[1] -
				6*X_Em[1]*X_Em[1]*X_Em[2]*X_Em[2] +
				X_Em[2]*X_Em[2]*X_Em[2]*X_Em[2])), 0.0);
    // Where are we?
    int decaInside = 0;    // Own Inside cause EFB can be shifted against Quad
    if ( ((-X_Em[3]-EmD.dz) > En_zmin) && ((-X_Em[3]-EmD.dz) < En_zmax) )
      decaInside += Magnet::EntranceField;
    if ( ((X_Ex[3]+ExD.dz) > Ex_zmin) && ((X_Ex[3]+ExD.dz) < Ex_zmax) )
      decaInside += Magnet::ExitField;
    if ( ((X_Em[3]+EmD.dz) > 0) && ((X_Ex[3]+ExD.dz) < 0) )
      decaInside += Magnet::HomogenField;

    if (decaInside & Magnet::EntranceField) {
      B_deca = GetDecaFringeField(X_Em, EmD);           // Calculate Field
    }
    if (decaInside & Magnet::ExitField) {
      B_deca = B_deca + GetDecaFringeField(X_Ex, ExD);  // Calculate Field
    }
    if ((decaInside & Magnet::ExitField) &&
	(decaInside & Magnet::EntranceField)) {
      B_deca = B_deca - B_deca_Uni;
    }
    if (Inside == Magnet::HomogenField) {
      B_deca = B_deca_Uni;
    }
    B = B + B_deca;
  }

  // ------------------------------ Dodecapole: --------------------------------
  static Vector3D B_dodeca;
  B_dodeca=Vector3D(0.0,0.0,0.0);
  if (Dodeca) {
    // Calculate Uniform field:
    Vector3D B_dodeca_Uni( (G[5]*(5*X_Em[1]*X_Em[1]*X_Em[1]*X_Em[1]*X_Em[2] -
				  10*X_Em[1]*X_Em[1]*X_Em[2]*X_Em[2]*X_Em[2] +
				  X_Em[2]*X_Em[2]*X_Em[2]*X_Em[2]*X_Em[2])), 
			   (G[5]*(X_Em[1]*X_Em[1]*X_Em[1]*X_Em[1]*X_Em[1] -
				  10*X_Em[1]*X_Em[1]*X_Em[1]*X_Em[2]*X_Em[2] +
				  5*X_Em[1]*X_Em[2]*X_Em[2]*X_Em[2]*X_Em[2])), 
			   0.0);
    // Where are we?
    int DodecaInside = 0;    // Own Inside cause EFB can be shifted against Quad
    if ( ((-X_Em[3]-EmDD.dz) > En_zmin) && ((-X_Em[3]-EmDD.dz) < En_zmax) )
      DodecaInside += Magnet::EntranceField;
    if ( ((X_Ex[3]+ExDD.dz) > Ex_zmin) && ((X_Ex[3]+ExDD.dz) < Ex_zmax) )
      DodecaInside += Magnet::ExitField;
    if ( ((X_Em[3]+EmDD.dz) > 0) && ((X_Ex[3]+ExDD.dz) < 0) )
      DodecaInside += Magnet::HomogenField;

    if (DodecaInside & Magnet::EntranceField) {
      B_dodeca = GetDodecaFringeField(X_Em, EmDD);           // Calculate Field
    }
    if (DodecaInside & Magnet::ExitField) {
      B_dodeca = B_dodeca + GetDodecaFringeField(X_Ex, ExDD);// Calculate Field
    }
    if ((DodecaInside & Magnet::ExitField) &&
	(DodecaInside & Magnet::EntranceField)) {
      B_dodeca = B_dodeca - B_dodeca_Uni;
    }
    if (Inside == Magnet::HomogenField) {
      B_dodeca = B_dodeca_Uni;
    }
    B = B + B_dodeca;
  }
  // rotate B back to ACS:
  return (rot_Ex_to_A * B);
}

// Set Magnetic Field to fraction (factor Fract) of nominal Field:
int 
Multipol::ScaleField(double Fract)   // i.e. Fract=0.5 sets half of nom. F.
{
  int i;
  for(i=1; i<=5; i++) G[i] = Fract*B_Nom[i]/pow(R,i);
  return 0;
}
