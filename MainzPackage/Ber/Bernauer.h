/**
 * \file Bernauer.h
 *
 * \ingroup MainzPackage
 * 
 * \brief Class def header for a class Bernauer
 *
 * @author erezcohen
 */

/** \addtogroup MainzPackage

    @{*/
#ifndef BERNAUER_H
#define BERNAUER_H

#include <iostream>
#define m_proton   0.938272046   ;
#define m_neutron  0.93956563    ;
#define m_deuteron 1.87561339    ;
#define m_e        0.510998928e-3;
#define m_C12      11.174866     ;
#define PI         3.14159265    ;
#define mu_p       2.79284739    ;

#define k_v 3.706;
#define k_s -0.12;
#define m_ro 0.776;
#define m_w 0.784;
#define m_phi 1.019;
#define m_ro_prime 1.45;

////PhysRevC.69.068201 (2004)
#define mu_n -1.913;
#define gama 0.515;
#define m_ro9 0.776;
#define m_phi 1.019;

/**
   \class Bernauer
   User defined class Bernauer ... these comments are used to generate
   doxygen documentation!
 */
class Bernauer{

public:

    
//    const double m_proton,  m_neutron , m_deuteron  , m_e ,  m_C12  , PI  , mu_p ;

    //VDM FF
    //PhysRevC.64.035204
    //Table 1 Title
//    double k_v , k_s, m_ro , m_w , m_phi , m_ro_prime;
    
    //Table 1 Content
    int currentModel;
    double g_roPrimeTof_roPrime[7] , k_ro_prime[7]  , g_wTof_w[7] , k_w[7] , g_phiTof_phi[7];
    double k_phi[7]  , mu_phy[7] ,  lamda1[7]  ,  lamdaD[7]  , lamda2[7]  , lamdaQCD[7] , N[7];
    
    double g_roPrimeTof_roPrime(){return g_roPrimeTof_roPrime[currentModel];};
    double k_ro_prime(){return k_ro_prime[currentModel];};
    double g_wTof_w(){return g_wTof_w[currentModel];};
    double k_w(){return k_w[currentModel];};
    double g_phiTof_phi(){return g_phiTof_phi[currentModel];};
    double k_phi(){return k_phi[currentModel];};
    double mu_phy(){return mu_phy[currentModel];};
    double lamda1(){return lamda1[currentModel];};
    double lamdaD(){return lamdaD[currentModel];};
    double lamda2(){return lamda2[currentModel];};
    double lamdaQCD(){return lamdaQCD[currentModel];};
    inline double sqr ( double x ){    return x * x; };
    
    
  
    
  /// Default constructor
    Bernauer();

  /// Default destructor
  ~Bernauer(){}

    
    //DipoleFit
    double dipole_fit(double q2) const;
    double G_ED(double q2) const;
    double G_MD(double q2) const;
    inline double Bernauer::Q2Tilda( const double Q2 );
    inline double g ( const double Q2 );
    inline double alph ( const double Q2,const double mass);
    inline double beta ( const double Q2,const double mass);
    inline double massRatio ( const double Q2,const double mass);
    inline double F1D ( const double Q2 );
    inline double F2D ( const double Q2 );
    inline double F1_10 ( const double Q2 );
    inline double F2_10 ( const double Q2 );
    inline double F1RO ( const double Q2 );
    inline double F2RO ( const double Q2 );
    inline double F1W ( const double Q2 );
    inline double F2W ( const double Q2 );
    inline double F1PHI ( const double Q2 );
    inline double F2PHI ( const double Q2 );
    inline double F1V ( const double Q2 );
    inline double F2V ( const double Q2 );
    inline double F1S ( const double Q2 );
    inline double F2S ( const double Q2 );
    inline double F1S9 ( const double Q2 );
    inline double F2S9 ( const double Q2 );
    inline double F1V9 ( const double Q2 );
    inline double F2V9 ( const double Q2 );
    inline double G_E6 ( const double q2 );
    inline double G_M6 ( const double q2 );
    inline double TanTh2O2 ( double q2 = -.4,double E = 0.855 );
    inline double G_E10 ( const double q2 );
    inline double G_M10 ( const double q2 );
    inline double G_E8 ( const double q2 );
    inline double G_E9 ( const double q2 );
    inline double G_M9 ( const double q2 );
    inline double Px_Pz_6 ( double q2 = -0.4 ,double E = 0.855 );
    inline double Px_Pz8 ( double q2 = -0.4 ,double E = 0.855 );
    inline double Px_Pz_10 ( double q2 = -0.4 ,double E = 0.855 );
    inline double PxPzBernauer_To_PxPzModifiedFF ( double q2 = -0.4 ,double E = 0.855 );
};

#endif
/** @} */ // end of doxygen group 

