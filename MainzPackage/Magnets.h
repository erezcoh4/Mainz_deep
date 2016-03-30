//								      -*-c++-*-
// Magnet.h - abstract class magnet and its derivates
//
// V1.0, T.Pospischil, 22.01.98
//
// Copyright (c) 1998-2001
//
// Institut für Kernphysik, Universität Mainz	tel. +49 6131 39-25802
// 55099 Mainz, Germany				fax  +49 6131 39-22964
//
// $Id: Magnets.h 2216 2008-06-13 21:13:47Z distler $
//

#ifndef __MAGNETS_H__
#define __MAGNETS_H__

#include "Matrix3D.h"
#include "FringeFall.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

// ACS: Absolute Coordinate System
// (z: start direction of the central ray, x downwards there)
// identical to Spek A target coordinate system

// Boundary Field width Extension factor:
const double bfwef = 1.4;

class Magnet{
public:
    virtual ~Magnet() { ; }
    enum Inside { HomogenField=1, EntranceField=2, ExitField=4 };
    
    // ---- Entrance boundary, EnCS: Entrance Coordinate System:
    Vector3D En;  // Origin of the EnCS in the ACS (mm)
    double En_phi;// Angle of EnCS z-axis against ACS z-axis in ACS zx-plane (deg)
    double C0[6]; // Coefficients of entrance fringing field fall off
    double En_xmin, En_xmax; // x width of Entrance Boundary in EnCS (mm)
    double En_zmin, En_zmax; // z length of Entrance Boundary in EnCS (mm)
    
    // ---- Exit boundary, ExCS: Exit Coordinate System:
    Vector3D Ex;  // Origin of the ExCS in the ACS (mm)
    double Ex_phi;// Angle of ExCS z-axis against ACS z-axis in ACS zx-plane (deg)
    double C1[6]; // Coefficients of exit fringing field fall off
    double Ex_xmin, Ex_xmax; // x width of Exit Boundary in ExCS (mm)
    double Ex_zmin, Ex_zmax; // z length of Entrance Boundary in ExCS (mm)
    
    
    // +++++++++++++++ Magnet is an abstract class, functions are 'pure virtual':
    
    // check if a ACS position X_A is between the entrance and the exit boundary:
    virtual int inside(Vector3D X_A)=0;
    // possible return values of inside() consists of 3 bits:
    
    // Calculate Magnetic Field at a given position:
    virtual Vector3D GetField(int Inside,  // result inside(r), if unkn. set it <0
                              Vector3D r)=0;// position in the ACS
    
    // Set Magnetic Field to fraction (factor Fract) of nominal Field:
    virtual int ScaleField(double Fract)=0; // i.e. Fract=0.5 sets half of nom. F.
    
    // write data file with boundary curvature:
    virtual int showBoundary(char *FileName, // Name of output data file
                             int StepNum)=0; // Number of points on boundary
    // write data file with entrance boundary curvature:
    virtual int showEntranceBoundary(char *FileName, // Name of output data file
                                     int StepNum)=0; // Number of points on bound.
    // write data file with exit boundary curvature:
    virtual int showExitBoundary(char *FileName,  // Name of output data file
                                 int StepNum)=0;  // Number of points on boundary
    
protected:
    Vector3D E_y;                      // Einheitsvektor in y-Richtung
    Matrix3D rot_A_to_En;              // rotation matrix from ACS to EnCS
    Matrix3D rot_En_to_A;              // rotation matrix from EnCS to ACS
    Matrix3D rot_A_to_Ex;              // rotation matrix from ACS to ExCS
    Matrix3D rot_Ex_to_A;              // rotation matrix from ExCS to ACS
    
    // ... extended boundary field widths:
    double bEn_xmin, bEn_xmax, bEx_xmin, bEx_xmax;
};











// +++++++++++++++++++++++++++++++++ Dipol ++++++++++++++++++++++++++++++++++++

class Dipol : public Magnet {
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
     int ScaleField(double Fract);   // i.e. Fract=0.5 sets half of nom. F. // virtual
    
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
    Dipol();      // Normal construktor, everything is set to some initial value
    // Constructor which initialises everything:
    Dipol(double iB_Nom,   // Homogenous Field B_y (Tesla)
          double iR,       // dipole (reference ray) radius (mm)
          double iD,       // Gap width (mm)
          Vector3D iEn,    // Origin of the EnCS in the ACS (mm)
          double iEn_phi,  // Angle of EnCS against ACS in zx-plane (deg)
          double iS0[9],   // polynom description entrance boundary, only 2-8 used
          double iC0[6],   // Coefficients of entrance fringing field fall off
          double iEn_xmin, // Min. of x-Range of Entrance Boundary
          double iEn_xmax, // Max. of x-Range of Entrance Boundary
          double iEn_zmin, // Min. of z-Range of Entrance Boundary
          double iEn_zmax, // Max. of z-Range of Entrance Boundary
          Vector3D iEx,    // Origin of the ExCS in the ACS (mm)
          double iEx_phi,  // Angle of ExCS against ACS in zx-plane (deg)
          double iS1[9],   // polynom description of exit boundary, only 2-8 used
          double iC1[6],   // Coefficients of exit fringing field fall off
          double iEx_xmin, // Min. of x-Range of Exit Boundary
          double iEx_xmax, // Max. of x-Range of Exit Boundary
          double iEx_zmin, // Min. of z-Range of Exit Boundary
          double iEx_zmax);// Max. of z-Range of Exit Boundary
    ~Dipol() { ; } // virtual
    
private:
    
    // calculate B_y in the Midplane of the Spectrometer at a position r
    double B_yMidPlane(int Inside,  // result of inside(r), if <0 value is det.
                       Vector3D r); // position in the ACS
    double dz_En(double x); // calculate Entrance boundary dz for a given EnCS x
    double dz_Ex(double x); // calculate Exit boundary dz for a given ExCS x
    
};












// +++++++++++++++++++++++++++++++ Multipol ++++++++++++++++++++++++++++++++++++

class Multipol : public Magnet {
public:
    
    // MCS: Multipole Coordinate System, directions like ExCS, but EnCS origin
    
    // B_Nom[1]: quad, [2]:sext, [3]:octa, [4]:deca, [5]:dodeca
    double B_Nom[6]; // nominal B field at r=R in Tesla,
    double G[6];     // Grandients in field calc: G[n] = Fract*B_Nom[n]/(R^n)
    
    double R;        // Radius of aperure (mm)
    double L;        // Effective length of the multipole (mm)
    
    // --------------------------------------Magnet virtual functions to implement:
    // check if a ACS position X_A is between the entrance and the exit boundary:
    int inside(Vector3D X_A);
    // Calculate Magnetic Field at a given position:
    Vector3D GetField(int Inside,  // result of inside(r), if unknown set it < 0
                      Vector3D r); // position in the ACS
    // Set Magnetic Field to fraction (factor Fract) of nominal Field:
    virtual int ScaleField(double Fract);   // i.e. Fract=0.5 sets half of nom. F.
    
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
    Multipol();   // Normal construktor, everything is set to some initial value
    // Constructor which initialises everything:
    //  EnCS and ExCS are assumed to be antiparallel with the same width
    //  (xmin to xmax) and a distance L between the origins
    Multipol(double iB_Nom[6],// Multipole field strength (Tesla)
             double iR,       // Radius of (quad) aperture (mm)
             double iL,       // Effective length of the multipole (mm)
             Vector3D iMCS,   // Origin of the MCS in the ACS (mm)
             double iphi,     // Angle of MCS against ACS in zx-plane (deg)
             double ixmin,    // Effective x-range min of the multipole
             double ixmax,    // Effective x-range max of the multipole
             double iC0[6],   // Coefficients of entrance fringing field fall off
             double iEn_zmin, // Min. of z-Range of Entrance Boundary
             double iEn_zmax, // Max. of z-Range of Entrance Boundary
             double iC1[6],   // Coefficients of exit fringing field fall off
             double iEx_zmin, // Min. of z-Range of Exit Boundary
             double iEx_zmax, // Max. of z-Range of Exit Boundary
             double FRH,      // fringe field fall off hexapole radius rel. to R
             double FRO,      // fringe field fall off octapole radius rel. to R
             double FRD,      // fringe field fall off decapole radius rel. to R
             double FRDD,     // fringe field fall off dodecapole rad. rel. to R
             double DSH,      // hexapole EFB inward displacement (DSH = dz/2R)
             double DSO,      // rel. octapole EFB inward displacement
             double DSD,      // rel. decapole EFB inward displacement
             double DSDD);    // rel. dodecapole EFB inward displacement
    ~Multipol() { ; } // virtual
    
private:
    
    int Quad, Hexa, Octa, Deca, Dodeca; // Booleans, !=0 if B_Nom[i] !=0
    
    // --- Fall off functions of Entrance fringing field (work in MCS):
    FringeFallFunc EmQ, EmH, EmO, EmD, EmDD;
    // --------- Fall off functions of Exit fringing field (work in ExCS):
    FringeFallFunc ExQ, ExH, ExO, ExD, ExDD;
    
    Vector3D GetQuadFringeField(Vector3D r,          // Position in fringe field
                                FringeFallFunc FF);  // fringe field func
    Vector3D GetHexaFringeField(Vector3D r,          // Position in fringe field
                                FringeFallFunc FF);  // fringe field func
    Vector3D GetOctaFringeField(Vector3D r,          // Position in fringe field
                                FringeFallFunc FF);  // fringe field func
    Vector3D GetDecaFringeField(Vector3D r,          // Position in fringe field
                                FringeFallFunc FF);  // fringe field func
    Vector3D GetDodecaFringeField(Vector3D r,        // Position in fringe field
                                  FringeFallFunc FF);// fringe field func
};

#endif

