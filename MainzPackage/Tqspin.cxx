#ifndef TQSPIN_CXX
#define TQSPIN_CXX

#include "Tqspin.h"


// function to store contents of y[] and magnetic Field in files:
void StoreStep(double t,   // independent variable
               double y[], // fields to store
               int ysize,  // size of y (y[1..ysize] must exist!)
               int stepN)  // number of the stored step
{
    return;  // do not store steps for all 10290 tracks
}


void magmove(double t,
             double y[],      // INPUT:  variables
             double dydx[])   // OUTPUT: their derivatives
{
    // ------------------------------ magmove -------------------------------------
    // routine which calculates the derivatives dydx[] of y[] with the ODEs
    // for the numerical Cash Karp Runge-Kutta solution:
    
    static Vector3D r;          // position
    static Vector3D v;          // velocity
    static Vector3D s;          // spin
    static Vector3D dr, dv, ds; // time derivates of position, velocity and spin
    static Vector3D B;          // Magnetic Field
    static Vector3D E_v;        // Vector in v-direction with length 1
    static Vector3D B_L, B_T;   // Fields longitudinal and transversal to v
    
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
    
    // Set long. B-Field to zero for test purpose:
    // B_L = Vector3D(0.0,0.0,0.0);
    // B_T = Vector3D(0.0,0.0,0.0);
    
    // ++++++++++ ODEs: Ordinary Differential Equations to integrate
    // (K_p = e/(M_p*gamma)):
    Vector3D Velocity( y[4] , y[5] , y[6] );
    double betaParticle  = Velocity.abs() / SpeedOfLight ;
    double gammaParticle = 1./sqrt( 1. - betaParticle*betaParticle );
    double K_p = 0.09578831/gammaParticle;           // [K_p] = 1/(Tesla*nsec)
    //    printf("in MagMove beta = %f, gamma = %f , K_p = %f\n",betaParticle,gammaParticle,K_p);
    dr = v;
    dv = K_p * ( v && B );
    ds = K_p * ( s && (g_p/2 * B_L + (1+gammaParticle*(g_p-2)/2) * B_T) );
    
    // Fill derivatives of position, velocity and spin-Vector to dydx:
    dydx[1]=dr[1]; dydx[2]=dr[2]; dydx[3]=dr[3];
    dydx[4]=dv[1]; dydx[5]=dv[2]; dydx[6]=dv[3];
    dydx[7]=ds[1]; dydx[8]=ds[2]; dydx[9]=ds[3];
    
//    delete r;
//    delete v;
//    delete s;
//    delete dr;
//    delete dv;
//    delete ds;
//    delete B;          // Magnetic Field
//    delete E_v;        // Vector in v-direction with length 1
//    delete B_L;
//    delete B_T;   // Fields longitudinal and transversal to v
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
}





TVector3 Tqspin::pSpinPrecessionSpecA(float fdp , float fth , float fph , float fy0 , float fp_ref, float fSx , float fSy , float fSz){
    // <dp/%c> <th/mrad> <y0/mm> <ph/mrad> <p_ref/MeV/c> <Sx_tg> <Sy_tg> <Sz_tg>
    
    // Vector containing start values for y:
    double *ystart = new double[odenum+1];
    
    double dp_val[5] = {-5.7794, 0.0327, 5.6206, 11.0182, 16.2521};
    
    // Loop for multitrack (is left after first track for single track calc):
    for (dp_c = 0; dp_c < 5; dp_c++){
        dp = dp_val[dp_c];
        for (y0 = -30.0; y0 <= 30.0; y0+=10.0)
            for (ph = -105.0; ph <= 105.0; ph+=35.0)
                for (p_ref = 480.0; p_ref <= 630.0; p_ref+=30.0)
                    for (th = -75.0; th <= 75.0; th+=25.0){
                        
                        // particle coordinates in Spectrometer-Target-System from command line:
                        dp = fdp;    // %
                        th = fth;    // th_tg in mrad
                        ph = fph;    // ph_tg in mrad
                        y0 = fy0;    // y0_tg in mm
                        p_ref = fp_ref; // MeV/c
                        cmd_Sx = fSx;
                        cmd_Sy = fSy;
                        cmd_Sz = fSz;
                        
                        Vector3D Spin_tg(cmd_Sx,cmd_Sy,cmd_Sz);
                        p = (1.0 + dp/100)*p_ref; // MeV/c
                        
                        spec->scaleField(p_ref);
                        
                        // calc start position and velocity:
                        Vector3D Position(0, y0, 0); // im mm !
                        Vector3D Velocity(tan(th/1000.0), tan(ph/1000.0), 1.0); // Richtung
                        Velocity = Velocity/Velocity.abs(); // Normierung des Betrags auf 1.0
                        gammaParticle = sqrt(p*p+M_p*M_p)/M_p;
                        betaParticle = sqrt(1.0 - 1.0/(gammaParticle*gammaParticle));
                        Velocity = Velocity * SpeedOfLight * betaParticle;
                        
                        // store Start_Velocity:
                        Vector3D Start_Velocity(tan(th/1000.0), tan(ph/1000.0), 1.0); // Richtung
                        // Normierung des Betrags auf 1.0
                        Start_Velocity = Start_Velocity/Start_Velocity.abs();
                        
                        // Vector containing start values for y:
                        ystart[1] = Position[1]; // x-Coordinate in mm
                        ystart[2] = Position[2]; // y-Coordinate in mm
                        ystart[3] = Position[3]; // z-Coordinate in mm
                        ystart[4] = Velocity[1]; // x-Velocity in mm/nsec
                        ystart[5] = Velocity[2]; // y-Velocity in mm/nsec
                        ystart[6] = Velocity[3]; // z-Velocity in mm/nsec
                        ystart[7] = Spin_tg[1];  // x component of the spin at the target (%)
                        ystart[8] = Spin_tg[2];  // y component of the spin at the target (%)
                        ystart[9] = Spin_tg[3];  // z component of the spin at the target (%)
                        
                        // K_p = e/(M_p*gamma)
                        K_p = 0.09578831/gammaParticle;           // [K_p] = 1/(Tesla*nsec)
                        
                        starttime = 0;
                        endtime = 1000.0;
                        nstored=0;
                        
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
                        Position[1] = ystart[1];
                        Position[2] = ystart[2];
                        Position[3] = ystart[3];
                        Velocity[1] = ystart[4];
                        Velocity[2] = ystart[5];
                        Velocity[3] = ystart[6];
                        betaParticle = Velocity.abs()/SpeedOfLight;
                        Vector3D Spin_end(ystart[7],ystart[8],ystart[9]);
                        
                        // Calculate End Spin in particle system (rot to HDC-System)
                        double winkel_hdc=10.0*deg;                        // main Angle of HDC
                        Vector3D Ez_hdc( -cos(winkel_hdc), 0.0, -sin(winkel_hdc) );
                        Vector3D Ex_hdc( -sin(winkel_hdc), 0.0,  cos(winkel_hdc) );
                        Vector3D Ey_hdc = Ez_hdc && Ex_hdc;
                        Vector3D DAxis_to_HDC;
                        double Dangle_to_HDC;
                        Vector3D SEAR_H = rotVec(Spin_end, Velocity, Ez_hdc,
                                                 DAxis_to_HDC, Dangle_to_HDC);
                        Vector3D Spin_end_particle(SEAR_H*Ex_hdc, SEAR_H*Ey_hdc, SEAR_H*Ez_hdc);

                        
                        return TVector3(Spin_end_particle[1],Spin_end_particle[2],Spin_end_particle[3]);
                      //                        printf("START SPIN:\n \\(S{_x}, S{_y}, S{_z}\\){_ACS}\n");
                        //                        printf("( %6.2f, %6.2f, %6.2f )\n", Spin_tg[1],Spin_tg[2],Spin_tg[3]);
                        //                        printf("\nEND SPIN:\n \\(S{_x}, S{_y}, S{_z}\\){_particle}\n");
                        //                        printf("( %6.2f, %6.2f, %6.2f )\n",Spin_end_particle[1],Spin_end_particle[2],Spin_end_particle[3]);
                        // close some files:
                        
                    }
    }
    
    
    // close some files:
    return TVector3();
}



#endif
