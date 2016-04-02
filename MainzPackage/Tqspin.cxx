#ifndef TQSPIN_CXX
#define TQSPIN_CXX

#include "Tqspin.h"


void StoreStep(double t, double y[], int ysize, int stepN) {
    // function to store contents of y[] and magnetic Field in files:
    return;  // do not store steps for all 10290 tracks
}


void magmove(double t,  double y[], double dydx[]){
    // ------------------------------ magmove -------------------------------------
    // routine which calculates the derivatives dydx[] of y[] with the ODEs
    // for the numerical Cash Karp Runge-Kutta solution:
    
    static Vector3D r , v , s , dr, dv, ds;     // position, velocity and spin (and time derivatives)
    static Vector3D B, E_v , B_L, B_T;          // Fields longitudinal and transversal to v
    
    // Fill position, velocity and spin-Vector:
    r[1] = y[1]; r[2] = y[2]; r[3] = y[3];
    v[1] = y[4]; v[2] = y[5]; v[3] = y[6];
    s[1] = y[7]; s[2] = y[8]; s[3] = y[9];
    
    B[1] = t;         // Dummy, makes warning 't unused' disappear
    // Calculate magnetic field at position r:
    TqspinSpectrometer *spec = new qspinSpecA();
    spec->getMagField(B,r);
    
    // Split calulated B in a longitudinal and a transversal part to v:
    E_v = v/v.abs();            // Vector in v-direction with length 1
    B_L = (E_v * B) * E_v;
    B_T = B - B_L;
    
    // ++++++++++ ODEs: Ordinary Differential Equations to integrate
    TVector3 Velocity( y[4] , y[5] , y[6] );
    double betaParticle  = Velocity.Mag() / SpeedOfLight ;
    double gammaParticle = 1./sqrt( 1. - betaParticle*betaParticle );
    double K_p = 0.09578831/gammaParticle;           // [K_p] = 1/(Tesla*nsec)
    dr = v;
    dv = K_p * ( v && B );
    ds = K_p * ( s && (g_p/2 * B_L + (1+gammaParticle*(g_p-2)/2) * B_T) );
    
    // Fill derivatives of position, velocity and spin-Vector to dydx:
    dydx[1]=dr[1]; dydx[2]=dr[2]; dydx[3]=dr[3];
    dydx[4]=dv[1]; dydx[5]=dv[2]; dydx[6]=dv[3];
    dydx[7]=ds[1]; dydx[8]=ds[2]; dydx[9]=ds[3];
    
    return;
}




Tqspin::Tqspin(){
    act_length = 0.0;
    r_old = Vector3D(0.0,0.0,0.0);
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++ define VDC-system:
    O_D_vdc_x1= Vector3D(-6081.59, 0.0, 4949.3); // unsicher !! (un certain)
    phi_bend_ref = 100.0238 *deg;
    // Tma-Offsets   (qsdda665.mxl):         (qsdda495.mxl2):
    tma_x     =   46.7; // mm          =  37.4; // mm
    tma_theta =  787.9; // mrad        = 785.0; // mrad
    tma_y     =   -5.1; // mm          =   1.2; // mm
    tma_phi   =  -1.15; // mrad        =   4.2; // mrad
    // Orientierung des VDC-Systems:
    alpha= (phi_bend_ref-90.0*deg+tma_theta/1000.0); // main Ang of VDC
    Ez_vdc = Vector3D( -cos(alpha), 0.0, -sin(alpha) );
    Ex_vdc = Vector3D( -sin(alpha), 0.0,  cos(alpha) );
    Ey_vdc = Ez_vdc && Ex_vdc;
    M_x = Matrix3D(Ex_vdc, (tma_phi/1000.0) );  // Phi-Offset der VDC
    Ey_vdc = M_x * Ey_vdc;  Ez_vdc = M_x * Ez_vdc;
    // Urspungslage des VDC-Systems:
    O_vdc = O_D_vdc_x1 -tma_x*Ex_vdc -tma_y*Ey_vdc;
    dp_val[0] = -5.7794;
    dp_val[1] = 0.0327;
    dp_val[2] = 5.6206;
    dp_val[3] = 11.0182;
    dp_val[4] = 16.2521;
    winkel_hdc=10.0*deg;                        // main Angle of HDC
    Ez_hdc = Vector3D( -cos(winkel_hdc), 0.0, -sin(winkel_hdc) );
    Ex_hdc = Vector3D( -sin(winkel_hdc), 0.0,  cos(winkel_hdc) );
    Ey_hdc = Ez_hdc && Ex_hdc ;

}



TVector3 Tqspin::pSpinPrecessionSpecA(float fdp , float fth , float fph , float fy0 , float fp_ref, TVector3 fS){
    // <dp/%c> <th/mrad> <y0/mm> <ph/mrad> <p_ref/MeV/c> <S_tg/%>
    
    double *ystart  = new double[odenum+1];
    starttime       = nstored = 0;
    endtime         = 1000.0;
    dp  = fdp;      th = fth;   ph = fph; y0 = fy0; p_ref = fp_ref;
    TVector3 Spin_tg(fS.X(),fS.Y(),fS.Z());
    p               = (1.0 + dp/100)*p_ref; // MeV/c
    spec            ->  scaleField(p_ref);
    
    // start position and velocity:
    Position.SetXYZ(0, y0, 0); // im mm !
    pGamma      = sqrt(p*p+M_p*M_p)/M_p;
    K_p         = 0.09578831/pGamma;           // [K_p] = 1/(Tesla*nsec)
    pBeta       = sqrt(1.0 - 1.0/(pGamma*pGamma));
    TVector3 Velocity( TVector3(tan(th/1000.0), tan(ph/1000.0), 1.0).Unit() *  SpeedOfLight * pBeta ); // Richtung

    // Vector containing start values for y:
    ystart[1] = Position.X(); // x-Coordinate in mm
    ystart[2] = Position.Y(); // y-Coordinate in mm
    ystart[3] = Position.Z(); // z-Coordinate in mm
    ystart[4] = Velocity.X(); // x-Velocity in mm/nsec
    ystart[5] = Velocity.Y(); // y-Velocity in mm/nsec
    ystart[6] = Velocity.Z(); // z-Velocity in mm/nsec
    ystart[7] = Spin_tg.X();  // x component of the spin at the target (%)
    ystart[8] = Spin_tg.Y();  // y component of the spin at the target (%)
    ystart[9] = Spin_tg.Z();  // z component of the spin at the target (%)
    
    
    
    // +++++++++++++++++++++++++++++++++++  integrate ODEs:
    rk.odeint(ystart,    // Start values
              odenum,    // number of ODEs (= dimension of the arrays)
              starttime, // start value for the independent variable
              endtime,   // stop at this value of the independent variable
              -8800,     // stop if y[1] smaller than y1min (mm)
              5e-6,      // required accuracy
              1,         // guessed first stepsize
              1e-20,     // min. allowed stepsize (can be zero)
              &nok,      // number of good steps taken
              &nbad,     // number of bad (retried and fixed) steps taken
              0.2,       // interval size for step storing (nsec)
              100000,    // max number of steps to store
              nstored,   // number of stored steps
              StoreStep, // function which writes the step to file
              magmove,   // derivatives function
              &acttime); // actual value for the independent variable
    
    // --- Get END position, velocity and spin: out of ystart
    Position.SetXYZ(ystart[1] , ystart[2] , ystart[3]);
    Velocity.SetXYZ(ystart[4] , ystart[5] , ystart[6]);
    Spin_end.SetXYZ(ystart[7],ystart[8],ystart[9]);
    // Calculate End Spin in particle system (rot to HDC-System)
    Vector3D DAxis_to_HDC;
    SEAR_H = rotVec(Vector3D(Spin_end.X(),Spin_end.Y(),Spin_end.Z()), Vector3D(Velocity.X(),Velocity.Y(),Velocity.Z()), Ez_hdc, DAxis_to_HDC, Dangle_to_HDC);
    TVector3 Spin_end_particle(SEAR_H*Ex_hdc, SEAR_H*Ey_hdc, SEAR_H*Ez_hdc);
    return Spin_end_particle;
}


#endif
