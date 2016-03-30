#ifndef DIPOL_CXX
#define DIPOL_CXX

#include "Magnets.h"


//								      -*-c++-*-
// Dipol.cc - implementation of the magnet class Dipol,
//            field calculation decribed in the RAYTRACE manual
//
// V1.0, T.Pospischil, 22.01.98
//
// Copyright (c) 1998-2001
//
// Institut f¸r Kernphysik, Universit‰t Mainz	tel. +49 6131 39-25802
// 55099 Mainz, Germany				fax  +49 6131 39-22964
//
// $Id: Dipol.cc 2216 2008-06-13 21:13:47Z distler $
//



// Normal construktor, everything is set to some initial value
Dipol::Dipol(){
    
    int i;
    
    B_Nom = 0.0;
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

// Constructor which initialises everything:
Dipol::Dipol(double iB_Nom,// Homogenous Field B_y (Tesla)
             double iR,      // dipole (reference ray) radius (mm)
             double iD,      // Gap width (mm)
             Vector3D iEn,   // Origin of the EnCS in the ACS (mm)
             double iEn_phi, // Angle of EnCS against ACS in zx-plane (deg)
             double iS0[9],  // polynom entrance boundary, only 2-8 used
             double iC0[6],   // Coeff. of entrance fringing field fall off
             double iEn_xmin,// Min. of x-Range of Entrance Boundary
             double iEn_xmax,// Max. of x-Range of Entrance Boundary
             double iEn_zmin, // Min. of z-Range of Entrance Boundary
             double iEn_zmax, // Max. of z-Range of Entrance Boundary
             Vector3D iEx,   // Origin of the ExCS in the ACS (mm)
             double iEx_phi, // Angle of ExCS against ACS in zx-plane (deg)
             double iS1[9],  // polynom exit boundary, only 2-8 used
             double iC1[6],   // Coefficients of exit fringing field fall off
             double iEx_xmin,// Min. of x-Range of Exit Boundary
             double iEx_xmax,// Max. of x-Range of Exit Boundary
             double iEx_zmin,// Min. of z-Range of Exit Boundary
             double iEx_zmax)// Max. of z-Range of Exit Boundary
{
    int i;
    
    B_Nom = iB_Nom;
    B_Act = B_Nom;
    R = iR;
    D = iD;
    
    En = iEn;
    En_phi = iEn_phi;
    En_xmin = iEn_xmin;
    En_xmax = iEn_xmax;
    En_zmin = iEn_zmin;
    En_zmax = iEn_zmax;
    
    Ex = iEx;
    Ex_phi = iEx_phi;
    Ex_xmin = iEx_xmin;
    Ex_xmax = iEx_xmax;
    Ex_zmin = iEx_zmin;
    Ex_zmax = iEx_zmax;
    
    for(i=0; i<=8; i++) { S0[i] = iS0[i]; S1[i] = iS1[i]; }
    for(i=0; i<=5; i++) { C0[i] = iC0[i]; C1[i] = iC1[i]; }
    
    E_y = Vector3D(0.0, 1.0, 0.0);       // Einheitsvektor in y-Richtung
    rot_A_to_En = Matrix3D(E_y, En_phi*deg);
    rot_En_to_A = Matrix3D(E_y, -1.0*En_phi*deg);
    rot_A_to_Ex = Matrix3D(E_y, Ex_phi*deg);
    rot_Ex_to_A = Matrix3D(E_y, -1.0*Ex_phi*deg);
    
    bEn_xmin = bfwef*iEn_xmin;
    bEn_xmax = bfwef*iEn_xmax;
    bEx_xmin = bfwef*iEx_xmin;
    bEx_xmax = bfwef*iEx_xmax;
    
}

// calculate Entrance boundary dz for a given EnCS x
double
Dipol::dz_En(double x)
{
    double dz= 0.0;
    double w = x/R;
    dz = w*w*(S0[2]+w*(S0[3]+w*(S0[4]+w*(S0[5]+w*(S0[6]+w*(S0[7]+w*S0[8]))))));
    return (-dz*R);
}

// calculate Exit boundary dz for a given ExCS x
double
Dipol::dz_Ex(double x)
{
    double dz= 0.0;
    double w = x/R;
    dz = w*w*(S1[2]+w*(S1[3]+w*(S1[4]+w*(S1[5]+w*(S1[6]+w*(S1[7]+w*S1[8]))))));
    return (-dz*R);
}

// ----------------------------- inside --------------------------------------
// check if a ACS position is between the entrance and the exit boundary
//  possible return values of inside() consists of 3 bits:
//   const int InHomogenField  = 1;
//   const int InEntranceField = 2;
//   const int InExitField     = 4;
int
Dipol::inside(Vector3D X_A)
{
    int inval = 0;
    Vector3D X_En, X_Ex;
    
    X_En = rot_A_to_En * (X_A - En);   // transform X_A to EnCS
    X_Ex = rot_A_to_Ex * (X_A - Ex);   // transform X_A to ExCS
    
    double dzEn = dz_En(X_En[1]);
    double dzEx = dz_Ex(X_Ex[1]);
    
    // Are we inside the homogenous field?
    if ( ((X_En[3]-dzEn) < 0) && ( (X_Ex[3]-dzEx) < 0) &&
        ( ((X_En[1] < bEn_xmax) && (X_En[1] > bEn_xmin)) ||         // En side
         ((X_Ex[1] < bEx_xmax) && (X_Ex[1] > bEx_xmin)) )   )      // Ex side
        inval += Magnet::HomogenField;
    // Are we inside the entrance field?
    if ( ((X_En[3]-dzEn)>En_zmin) && ((X_En[3]-dzEn)<En_zmax) &&
        ((X_En[1] < bEn_xmax)    && (X_En[1] > bEn_xmin))       )    // En side
        inval += Magnet::EntranceField;
    // Are we inside the exit field?
    if ( ((X_Ex[3]-dzEx)>Ex_zmin) && ((X_Ex[3]-dzEx)<Ex_zmax) &&
        ((X_Ex[1] < bEx_xmax) && (X_Ex[1] > bEx_xmin))          )    // Ex side
        inval += Magnet::ExitField;
    
    return inval;
}


// ---------------------------- B_yMidPlane ----------------------------------
// calculate B_y in the Midplane of the Spectrometer at a postion r
double
Dipol::B_yMidPlane(int Inside,  // result of inside(r), if <0 value is det.
                   Vector3D r)  // position in the ACS
{
    // Inside value not known -> determine
    if (Inside<0) Inside = inside(r);
    if (!Inside) return 0.0;
    
    // Fringe Field regions:
    double B_y_En=0.0, B_y_Ex=0.0;
    if (Inside & Magnet::EntranceField) {
        
        // trafo from ACS to the EnCS
        Vector3D X_En = rot_A_to_En * (r - En);
        // calculate field:
        const double s = (X_En[3] - dz_En(X_En[1]))/D;
        const double S = C0[0]+s*(C0[1]+s*(C0[2]+s*(C0[3]+s*(C0[4]+s*(C0[5])))));
        B_y_En = (B_Act/(1+exp(S)));
        
        if ( !(Inside & Magnet::ExitField) ) return B_y_En;
    }
    if (Inside & Magnet::ExitField) {
        
        // trafo from ACS to the EnCS
        Vector3D X_Ex = rot_A_to_Ex * (r - Ex);
        // calculate field:
        const double s = (X_Ex[3] - dz_Ex(X_Ex[1]))/D;
        const double S = C1[0]+s*(C1[1]+s*(C1[2]+s*(C1[3]+s*(C1[4]+s*(C1[5])))));
        B_y_Ex = (B_Act/(1+exp(S)));
        
        if ( !(Inside & Magnet::EntranceField) ) return B_y_Ex;
    }
    // We are in Entrance- AND Exit Fringe Field
    if ( (Inside & Magnet::EntranceField) && (Inside & Magnet::ExitField) )
        return (B_y_En + B_y_Ex - B_Act);
    
    // We are only in the homogenous region of the Dipol:
    if (Inside == Magnet::HomogenField) return B_Act;
    
    fprintf(stderr,"WARNING: Undefined Dipol.inside value (%d)!\n",Inside);
    return 0.0;
}

// --------------------------- GetField ----------------------------------------
// Calculate Magnetic Field at a given position:
Vector3D
Dipol::GetField(int Inside,  // result of inside(r), if <0 value it is determin.
                Vector3D r)  // position in the ACS
{
    Vector3D B;                // Magnetic field
    double Delta=0.01*D;       // grid distance for numerical derivatives
    Vector3D rg;               // grid positions
    
    double B_m2_0, B_m1_m1, B_m1_0, B_m1_p1, B_0_p2, B_0_p1, B_0_0;
    double B_0_m1, B_0_m2, B_p1_m1, B_p1_0, B_p1_p1, B_p2_0;
    
    // Calculate B_y on a 13 point grid in the median plane:
    rg = r;
    rg[3] -= 2*Delta;  B_m2_0  = B_yMidPlane(Inside,rg);
    rg[3] += Delta;
    rg[1] -= Delta;    B_m1_m1 = B_yMidPlane(Inside,rg);
    rg[1] += Delta;    B_m1_0  = B_yMidPlane(Inside,rg);
    rg[1] += Delta;    B_m1_p1 = B_yMidPlane(Inside,rg);
    rg[3] += Delta;
    rg[1] += Delta;    B_0_p2  = B_yMidPlane(Inside,rg);
    rg[1] -= Delta;    B_0_p1  = B_yMidPlane(Inside,rg);
    rg[1] -= Delta;    B_0_0   = B_yMidPlane(Inside,rg);
    rg[1] -= Delta;    B_0_m1  = B_yMidPlane(Inside,rg);
    rg[1] -= Delta;    B_0_m2  = B_yMidPlane(Inside,rg);
    rg[3] += Delta;
    rg[1] += Delta;    B_p1_m1 = B_yMidPlane(Inside,rg);
    rg[1] += Delta;    B_p1_0  = B_yMidPlane(Inside,rg);
    rg[1] += Delta;    B_p1_p1 = B_yMidPlane(Inside,rg);
    rg[3] += Delta;
    rg[1] -= Delta;    B_p2_0  = B_yMidPlane(Inside,rg);
    
    // Calculate the 3 components of the magnetic field B by Taylor expansion
    // in y out of the mid plane (<= Maxwell's equations):
    
    double yD = r[2]/Delta;
    double yD2 = yD*yD;
    double yD3 = yD*yD2;
    double yD4 = yD*yD3;
    
    double B00_4     = 4.0 * B_0_0;
    double Term1     = B_p1_0 + B_m1_0 + B_0_p1 + B_0_m1 - B00_4;
    double Term2     = B_p2_0 + B_m2_0 + B_0_p2 + B_0_m2 - B00_4;
    
    // y component:
    B[2] = B_0_0
    - yD2 * (  2.0/3.0 * Term1 - 1.0/24.0 * Term2 )
    + yD4 * ( -1.0/6.0 * Term1 + 1.0/24.0 * Term2
             +1.0/12.0* (B_p1_p1+B_m1_p1+B_p1_m1+B_m1_m1 - 2*Term1 - B00_4) );
    // z component:
    B[3] = yD * ( 2.0/3.0*(B_p1_0-B_m1_0) - 1.0/12.0*(B_p2_0-B_m2_0) )
    + yD3 * ( 1.0/6.0*(B_p1_0-B_m1_0)
             -1.0/12.0*(B_p2_0-B_m2_0 + B_p1_p1+B_p1_m1
                        - B_m1_p1-B_m1_m1 + 2*(B_m1_0-B_p1_0)) );
    // x component:
    B[1] = yD * ( 2.0/3.0*(B_0_p1-B_0_m1) - 1.0/12.0*(B_0_p2-B_0_m2) )
    + yD3 * ( 1.0/6.0*(B_0_p1-B_0_m1)
             -1.0/12.0*(B_0_p2-B_0_m2 + B_p1_p1+B_m1_p1
                        - B_p1_m1-B_m1_m1 + 2*(B_0_m1-B_0_p1)) );
    
    return B;
}

// Set Magnetic Field to fraction (factor Fract) of nominal Field:
int
Dipol::ScaleField(double Fract)   // i.e. Fract=0.5 sets half of nom. F.
{
    B_Act = Fract*B_Nom;
    return 0;
}

// --------------------------- showBoundary ------------------------------------
// write data file with boundary curvature:
int
Dipol::showBoundary(char *FileName, // Name of output data file
                    int StepNum)    // Number of points on boundary
{
    // open file:
    FILE *outfile;
    if ( NULL == (outfile = fopen(FileName,"w")) ) {
        fprintf(stderr," Cannot open file %s for output!\n",FileName);
        return 2;
    }
    // run through entrance boundary:
    double xstep = (En_xmax - En_xmin) / StepNum;
    Vector3D X_E, X_A;
    for(X_E[1]=En_xmin; X_E[1]<En_xmax; X_E[1] += xstep){
        X_E[3] = dz_En(X_E[1]);
        X_A = rot_En_to_A*X_E + En;  // trafo to ACS
        fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
    }
    // run through exit boundary:
    xstep = (Ex_xmax - Ex_xmin) / StepNum;
    for(X_E[1]=Ex_xmin; X_E[1]<Ex_xmax; X_E[1] += xstep){
        X_E[3] = dz_Ex(X_E[1]);
        X_A = rot_Ex_to_A*X_E + Ex;  // trafo to ACS
        fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
    }
    // close path:
    X_E[1]=En_xmin; X_E[3] = dz_En(X_E[1]);
    X_A = rot_En_to_A*X_E + En;  // trafo to ACS
    fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
    
    fclose(outfile);
    return 0;
}

// ------------------------- showEntranceBoundary -----------------------------
// write data file with entrance boundary curvature:
int
Dipol::showEntranceBoundary(char *FileName, // Name of output data file
                            int StepNum)    // Number of points on boundary
{
    // open file:
    FILE *outfile;
    if ( NULL == (outfile = fopen(FileName,"w")) ) {
        fprintf(stderr," Cannot open file %s for output!\n",FileName);
        return 2;
    }
    // run forwards through zmin boundary:
    double xstep = bfwef * (En_xmax - En_xmin) / StepNum;
    Vector3D X_E, X_A;
    for(X_E[1]=bfwef*En_xmin; X_E[1]<bfwef*En_xmax; X_E[1] += xstep){
        X_E[3] = dz_En(X_E[1]) + En_zmin;
        X_A = rot_En_to_A*X_E + En;  // trafo to ACS
        fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
    }
    // run backwards through zmax boundary:
    for(X_E[1]=bfwef*En_xmax; X_E[1]>bfwef*En_xmin; X_E[1] -= xstep){
        X_E[3] = dz_En(X_E[1]) + En_zmax;
        X_A = rot_En_to_A*X_E + En;  // trafo to ACS
        fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
    }
    // close path:
    X_E[1]=bfwef*En_xmin;
    X_E[3] = dz_En(X_E[1]) + En_zmin;
    X_A = rot_En_to_A*X_E + En;  // trafo to ACS
    fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
    
    fclose(outfile);
    return 0;
}

// ------------------------- showExitBoundary -----------------------------
// write data file with exit boundary curvature:
int
Dipol::showExitBoundary(char *FileName, // Name of output data file
                        int StepNum)    // Number of points on boundary
{
    // open file:
    FILE *outfile;
    if ( NULL == (outfile = fopen(FileName,"w")) ) {
        fprintf(stderr," Cannot open file %s for output!\n",FileName);
        return 2;
    }
    // run forwards through zmin boundary:
    double xstep = bfwef * (Ex_xmax - Ex_xmin) / StepNum;
    Vector3D X_E, X_A;
    for(X_E[1]=bfwef*Ex_xmin; X_E[1]<bfwef*Ex_xmax; X_E[1] += xstep){
        X_E[3] = dz_Ex(X_E[1]) + Ex_zmin; 
        X_A = rot_Ex_to_A*X_E + Ex;  // trafo to ACS 
        fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
    }
    // run backwards through zmax boundary:
    for(X_E[1]=bfwef*Ex_xmax; X_E[1]>bfwef*Ex_xmin; X_E[1] -= xstep){
        X_E[3] = dz_Ex(X_E[1]) + Ex_zmax; 
        X_A = rot_Ex_to_A*X_E + Ex;  // trafo to ACS 
        fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
    }
    // close path:
    X_E[1]=bfwef*Ex_xmin;
    X_E[3] = dz_Ex(X_E[1]) + Ex_zmin;
    X_A = rot_Ex_to_A*X_E + Ex;  // trafo to ACS 
    fprintf(outfile,"%8.3f %8.3f\n",X_A[3], -1.0*X_A[1]);
    
    fclose(outfile);
    return 0;
}



#endif
