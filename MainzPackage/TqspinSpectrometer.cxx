#ifndef TQSPINSPECTROMETER_CXX
#define TQSPINSPECTROMETER_CXX

#include "TqspinSpectrometer.h"



qspinSpecA::qspinSpecA(const char *basename)
:  TqspinSpectrometer(basename)
{
    // Quadrupole
    double Q_FringeFall[6] = {0.10120, 6.5287, -1.5947, 2.2773, 0.496, 0.0};
    double Q_Fields[6] = {0.0, -0.51615, 0.0, 0.04405, 0.0, 0.00339};
    Quad = new Multipol(Q_Fields,// Multipole field strength (Tesla)
                        200.0,           // Radius of (quad) aperture (mm)
                        1016.0,          // Effective length of the multipole (mm)
                        Vector3D(0,0,792), // Origin of the MCS in the ACS (mm)
                        0.0,             // Angle of MCS against ACS in zx-plane (deg)
                        -400.0,          // Effective x-range min of the multipole
                        400.0,           // Effective x-range max of the multipole
                        Q_FringeFall,    // Coeff of entrance fringing field fall off
                        -500.0,          // Min. of z-Range of Entrance Boundary
                        500.0,           // Max. of z-Range of Entrance Boundary
                        Q_FringeFall,    // Coefficients of exit fringing field fall off
                        -500.0,          // Min. of z-Range of Exit Boundary
                        500.0,           // Max. of z-Range of Exit Boundary
                        1.0,             // ff fall off hexapole radius rel. to R
                        1.0,             // ff fall off octapole radius rel. to R
                        1.0,             // ff fall off decapole radius rel. to R
                        1.0,             // ff fall off dodecapole rad. rel. to R
                        0.0,             // hexapole EFB inward displacem. (DSH = dz/2R)
                        0.254,           // rel. octapole EFB inward displacement
                        0.0,             // rel. decapole EFB inward displacement
                        0.5);            // rel. dodecapole EFB inward displacement
    
    // Sextupole
    double S_Fields[6] = {0.0, 0.0, 0.33179, 0.0, 0.02651, 0.0};
    double S_FringeFall[6] = {-0.03358,9.2438,-1.9596,-5.2344,-12.2224,60.3104};
    Sext = new Multipol(S_Fields,// Multipole field strength (Tesla)
                        400.0,           // Radius of (quad) aperture (mm)
                        470.0,           // Effective length of the multipole (mm)
                        Vector3D(0,0,2115), // Origin of the MCS in the ACS (mm)
                        0.0,             // Angle of MCS against ACS in zx-plane (deg)
                        -600,            // Effective x-range min of the multipole
                        600,             // Effective x-range max of the multipole
                        S_FringeFall,    // Coeff of entrance fringing field fall off
                        -600.0,          // Min. of z-Range of Entrance Boundary
                        500.0,           // Max. of z-Range of Entrance Boundary
                        S_FringeFall,    // Coefficients of exit fringing field fall off
                        -600.0,          // Min. of z-Range of Exit Boundary
                        500.0,           // Max. of z-Range of Exit Boundary
                        0.8,             // ff fall off hexapole radius rel. to R
                        0.8,             // ff fall off octapole radius rel. to R
                        0.8,             // ff fall off decapole radius rel. to R
                        0.8,             // ff fall off dodecapole rad. rel. to R
                        0.0,             // hexapole EFB inward displacem. (DSH = dz/2R)
                        0.0,             // rel. octapole EFB inward displacement
                        0.0,             // rel. decapole EFB inward displacement
                        0.0);            // rel. dodecapole EFB inward displacement
    
    // Dipole 1
    double D_FringeFall[6] = {0.148, 1.8674, -0.3337, 0.6745, -0.1445, 0.0426};
    double D1EnBound[9] = {0, 0, 0.416, -0.502, 0.03, 0.225, 0, 0, 0};
    double D1ExBound[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    D1 = new Dipol(1.50101,  // Homogenious Field By (Tesla)
                   1400,     // Radius of curvature (mm)
                   200,      // Gap width (mm)
                   Vector3D(0,0,3050),        // Origin of the EnCS in the ACS (mm)
                   169,      // Angle of EnCS against ACS in zx-plane ( 180 + 11 deg)
                   D1EnBound,// entrance boundary curvature polynom coefficients
                   D_FringeFall, // Coefficients of entrance fringing field fall off
                   -750,     // Min. of x-Range of Entrance Boundary
                   1000,     // Max. of x-Range of Entrance Boundary
                   -550,     // Min. of z-Range of Entrance Boundary
                   500,     // Max. of z-Range of Entrance Boundary
                   Vector3D(-597,0,4196.8),   // Origin of the ExCS in the ACS (mm)
                   33,       // Angle of ExCS against ACS in zx-plane (      -33 deg)
                   D1ExBound,// exit boundary curvature polynom coefficients
                   D_FringeFall, // Coefficients of exit fringing field fall off
                   -800,     // Min. of x-Range of Exit Boundary
                   900,      // Max. of x-Range of Exit Boundary
                   -550,     // Min. of z-Range of Exit Boundary
                   750);     // Max. of z-Range of Exit Boundary
    
    // Dipole 2
    double D2EnBound[9] = {0, 0, -0.062, 0.063, 0.095, 0.01, 0, 0, 0};
    double D2ExBound[9] = {0, 0, -0.086, 0.152, -0.18, 0.15, -0.085, 0.015, 0};
    D2 = new Dipol(1.50101,  // Homogenious Field By (Tesla)
                   1400,     // Radius of curvature (mm)
                   200,      // Gap width (mm)
                   Vector3D(-2071.4,0,5229.2), // Origin of the EnCS in the ACS (mm)
                   -83.4,     // Angle of EnCS against ACS in zx-plane ( 90 - 6.6 deg)
                   D2EnBound, // entrance boundary curvature polynom coefficients
                   D_FringeFall, // Coefficients of entrance fringing field fall off
                   -1600,     // Min. of x-Range of Entrance Boundary
                   1000,      // Max. of x-Range of Entrance Boundary
                   -550,      // Min. of z-Range of Entrance Boundary
                   750,       // Max. of z-Range of Entrance Boundary
                   Vector3D(-3120,0,5460.7),   // Origin of the ExCS in the ACS (mm)
                   114,       // Angle of ExCS against ACS in zx-plane (-90 - 24 deg)
                   D2ExBound, // exit boundary curvature polynom coefficients
                   D_FringeFall, // Coefficients of exit fringing field fall off
                   -800,      // Min. of x-Range of Exit Boundary
                   1800,      // Max. of x-Range of Exit Boundary
                   -550,      // Min. of z-Range of Exit Boundary
                   750);      // Max. of z-Range of Exit Boundary
}

int qspinSpecA::getMagField(Vector3D &B, Vector3D &r)
{
    int bval=0; // last determined inside-value
    
    B[1]=0.0; B[2]=0.0; B[3]=0.0; // Set B initially to zero
    //B[1]=-0.47e-4; B[2]=-0.17e-4; B[3]=0.0; // Test: Earth-Field ~0.5e-4 Tesla
    
    // +++++++  Here, one should put in all defined magnetic elements:
    
    // normally 'if ( (bval=Element.inside(r)) ) B=B+Element.GetField(bval,r);'
    // is enough, but adding if-statements about the rough position makes
    // the program running faster.
    
    // Quadrupol:
    if ((r[3] < 2310)){
        if ( (bval=Quad->inside(r)) ) B=B+Quad->GetField(bval,r);
    }
    // Sextupol:
    if ((r[3] > 1600) && (r[3] < 3100)){
        if ( (bval=Sext->inside(r)) ) B=B+Sext->GetField(bval,r);
    }
    // Dipol 1:
    if ((r[3] > 2400) && (r[1]>-2000))
        if ( (bval=D1->inside(r)) ) B=B+D1->GetField(bval,r);
    // Dipol 2:
    if ((r[3] > 3500) && (r[1]<-1000))
        if ( (bval=D2->inside(r)) ) B=B+D2->GetField(bval,r);
    
    return 0;
}

int qspinSpecA::showBoundary()
{
    // --- write magnet boundary data files:
    Quad->showBoundary((char*)"dat/quad.dat", 50);
    Quad->showEntranceBoundary((char*)"dat/quad_En.dat", 50);
    Quad->showExitBoundary((char*)"dat/quad_Ex.dat", 50);
    Sext->showBoundary((char*)"dat/sext.dat", 50);
    Sext->showEntranceBoundary((char*)"dat/sext_En.dat", 50);
    Sext->showExitBoundary((char*)"dat/sext_Ex.dat", 50);
    D1->showBoundary((char*)"dat/d1.dat", 50);
    D1->showEntranceBoundary((char*)"dat/d1_En.dat", 50);
    D1->showExitBoundary((char*)"dat/d1_Ex.dat", 50);
    D2->showBoundary((char*)"dat/d2.dat", 50);
    D2->showEntranceBoundary((char*)"dat/d2_En.dat", 50);
    D2->showExitBoundary((char*)"dat/d2_Ex.dat", 50);
    
    return 0;
}

int qspinSpecA::checkBoundary()
{
    FILE *inquad_En, *inquad, *inquad_Ex;
    FILE *insext_En, *insext, *insext_Ex;
    FILE *ind1_En, *ind1, *ind1_Ex;
    FILE *ind2_En, *ind2, *ind2_Ex;
    FILE *outside;
    inquad_En=fopen("ps/in_quad_En.dat","w");
    inquad   =fopen("ps/in_quad.dat","w");
    inquad_Ex=fopen("ps/in_quad_Ex.dat","w");
    insext_En=fopen("ps/in_sext_En.dat","w");
    insext   =fopen("ps/in_sext.dat","w");
    insext_Ex=fopen("ps/in_sext_Ex.dat","w");
    ind1_En=fopen("ps/in_d1_En.dat","w");
    ind1   =fopen("ps/in_d1.dat","w");
    ind1_Ex=fopen("ps/in_d1_Ex.dat","w");
    ind2_En=fopen("ps/in_d2_En.dat","w");
    ind2   =fopen("ps/in_d2.dat","w");
    ind2_Ex=fopen("ps/in_d2_Ex.dat","w");
    outside=fopen("ps/outside.dat","w");
    if ( ind1_En && ind1 && ind1_Ex && ind2_En && ind2 && ind2_Ex &&
        inquad_En && inquad && inquad_Ex && insext_En && insext && insext_Ex
        && outside ) {
        Vector3D X;
        int inval, j;
        std::cout << " please wait ..... " << std::endl;
        for(j=0; j<=20000; j++){
            X[3] = 8000.0*(rand()/(RAND_MAX+1.0));
            X[1] = -5000.0 + 7000.0*(rand()/(RAND_MAX+1.0));
            if ( (inval=Quad->inside(X)) ) {
                if ((inval & Magnet::HomogenField) && ((rand()/(RAND_MAX+1.0))>0.5))
                    fprintf(inquad,"%8.3f %8.3f \n",X[3], -1.0*X[1]);
                if ((inval & Magnet::EntranceField) && ((rand()/(RAND_MAX+1.0))>0.5))
                    fprintf(inquad_En,"%8.3f %8.3f \n",X[3], -1.0*X[1]);
                if ((inval & Magnet::ExitField) && ((rand()/(RAND_MAX+1.0))>0.5))
                    fprintf(inquad_Ex,"%8.3f %8.3f \n",X[3], -1.0*X[1]);
            }
            if ( (inval=Sext->inside(X)) ) {
                if ((inval & Magnet::HomogenField) && ((rand()/(RAND_MAX+1.0))>0.5))
                    fprintf(insext,"%8.3f %8.3f \n",X[3], -1.0*X[1]);
                if ((inval & Magnet::EntranceField) && ((rand()/(RAND_MAX+1.0))>0.5))
                    fprintf(insext_En,"%8.3f %8.3f \n",X[3], -1.0*X[1]);
                if ((inval & Magnet::ExitField) && ((rand()/(RAND_MAX+1.0))>0.5))
                    fprintf(insext_Ex,"%8.3f %8.3f \n",X[3], -1.0*X[1]);
            }
            if ( (inval=D1->inside(X)) ) {
                if ((inval & Magnet::HomogenField) && ((rand()/(RAND_MAX+1.0))>0.5))
                    fprintf(ind1,"%8.3f %8.3f \n",X[3], -1.0*X[1]);
                if ((inval & Magnet::EntranceField) && ((rand()/(RAND_MAX+1.0))>0.5))
                    fprintf(ind1_En,"%8.3f %8.3f \n",X[3], -1.0*X[1]);
                if ((inval & Magnet::ExitField) && ((rand()/(RAND_MAX+1.0))>0.5))
                    fprintf(ind1_Ex,"%8.3f %8.3f \n",X[3], -1.0*X[1]);
            }
            if ( (inval=D2->inside(X)) ) {
                if ((inval & Magnet::HomogenField) && ((rand()/(RAND_MAX+1.0))>0.5))
                    fprintf(ind2,"%8.3f %8.3f \n",X[3], -1.0*X[1]);
                if ((inval & Magnet::EntranceField) && ((rand()/(RAND_MAX+1.0))>0.5))
                    fprintf(ind2_En,"%8.3f %8.3f \n",X[3], -1.0*X[1]);
                if ((inval & Magnet::ExitField) && ((rand()/(RAND_MAX+1.0))>0.5))
                    fprintf(ind2_Ex,"%8.3f %8.3f \n",X[3], -1.0*X[1]);
            }
            else if (!D1->inside(X))
                fprintf(outside,"%8.3f %8.3f \n",X[3], -1.0*X[1]);
            
        }
        fclose(ind1); fclose(ind1_En); fclose(ind1_Ex);
        fclose(ind2); fclose(ind2_En); fclose(ind2_Ex);
        fclose(insext); fclose(insext_En); fclose(insext_Ex);
        fclose(inquad); fclose(inquad_En); fclose(inquad_Ex);
        fclose(outside);
    } else {
        fprintf(stderr," Couldn't open some boundary check data file!\n");
        return 1;
    }
    return 0;
}




#endif
