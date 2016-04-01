/**
 * \file Tqspin.h
 *
 * \ingroup MainzPackage
 *
 * \brief Class def header for a class Tqspin
 *
 * @author erezcohen
 */

/** \addtogroup MainzPackage
 
 @{*/
#ifndef TQSPIN_H
#define TQSPIN_H

#include <iostream>
using namespace std;
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "TqspinSpectrometer.h"
#include "CK_RungeKutta.h"
#include "TVector3.h"

// some constants:
#define SpeedOfLight 299.792458 // c in mm/nsec
#define M_p 938.27231 // MeV
#define g_p 5.58569478


/**
 \class Tqspin
 User defined class Tqspin ... these comments are used to generate
 doxygen documentation!
 */
class Tqspin{
    
public:
    
    
    // from qspin.cc
    
    
    TqspinSpectrometer *spec = new qspinSpecA();
    
    
    // some nasty global (!!) variables:
    double gammaParticle;
    double K_p;
    int multitrack = 0;   // Flag, assume we want to calc only one track
    static const int odenum = 9;    // Number of ODEs (dimension of the arrays)
    CK_RungeKutta rk;
    
    
    // ---------------------------- StoreStep ------------------------------------
    FILE   *storefile=NULL;  // filehandle of open file for step storage
    FILE   *storefile2=NULL; // filehandle for storage of magnetic field strength
    Vector3D r_old;
    double act_length;
    
    
    Vector3D O_D_vdc_x1; // unsicher !! (un certain)
    double phi_bend_ref , tma_x , tma_theta , tma_y , tma_phi , alpha;
    Vector3D Ez_vdc , Ex_vdc , Ey_vdc;
    Matrix3D M_x;  // Phi-Offset der VDC
    Vector3D O_vdc;
    double p, dp, th, ph, p_ref, y0;       // Spectrometer-Target-coordinates
    double cmd_Sx, cmd_Sy, cmd_Sz;         // Spin comp. read from command line
    int dp_c; // dp-counter
    int nok, nbad, nstored;
    double acttime , starttime , endtime;
    double betaParticle;
    
    
    
    TVector3    Position , Velocity;
    Float_t     pGamma , pBeta;
    /// Default constructor
    Tqspin();
    
    /// Default destructor
    ~Tqspin(){}
    
    
    
    
    // methods from qspin.cc
    TVector3 pSpinPrecessionSpecA (float, float, float, float, float, TVector3 );
    
};

#endif
/** @} */ // end of doxygen group

