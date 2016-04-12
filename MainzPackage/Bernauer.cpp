#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreePlayer.h"
#include "TMath.h"
#include "TCut.h"
#include "THStack.h"
#include "TProfile.h"
#include "TStyle.h"
#include <iostream>
#include <fstream>
#include "TCanvas.h"
#include "TChain.h"
#include "TH2F.h"
#include "TRandom2.h"
#include "TDirectory.h"


using namespace std;
using namespace TMath;

inline double sqr ( double x )
{    return x * x; };

const double m_proton   = 0.938272046   ;
const double m_neutron  = 0.93956563    ;
const double m_deuteron = 1.87561339    ;
const double m_e        = 0.510998928e-3;
const double m_C12      = 11.174866     ;
const double PI         = 3.14159265    ;
const double mu_p       = 2.79284739    ;

//DipoleFit
double dipole_fit(double q2) const{
    return 1/sqr(1 - q2/0.71);
}
double G_ED(double q2) const{
    return dipole_fit(q2);
}
double G_MD(double q2) const{
    return  mu_p * dipole_fit(q2);
}


//VDM FF
//PhysRevC.64.035204
//Table 1 Title
double k_v = 3.706;
double k_s = -0.12;
double m_ro = 0.776;
double m_w = 0.784;
double m_phi = 1.019;
double m_ro_prime = 1.45;


//Table 1 Content
int currentModel = 0;
double g_roPrimeTof_roPrime[7] = {0.4466,0.0514,0.3223,0.1013,0.0808,0.0625,0.0636};
double k_ro_prime[7] = {4.3472,23.533,4.982,-15.87,-17.993,0.9397,-0.4175};
double g_wTof_w[7] = {0.4713,0.0588,0.344,0.6604,0.8038,0.8029,0.7918};
double k_w[7] = {21.762,18.934,40.661,8.847,4.0526,5.5225,5.1109};
double g_phiTof_phi[7] = {-0.8461,-0.5283,-0.9315,-0.4054,-0.2336,-0.3070,-0.3011};
double k_phi[7] = {11.849,1.2236,14.6805,13.6415,13.5963,14.4123,13.4385};
double mu_phy[7] = {1.1498,1.1670,1.1411,1.127,1.1218,1.2379,1.1915};
double lamda1[7] = {0.9006,0.5902,0.8956,0.89361,0.9295,0.9916,0.9660};
double lamdaD[7] = {1.7038,0.7273,1.7038,1.0454,1.2207,1.2589,1.3406};
double lamda2[7] = {1.1336,1.9368,0.9551,2.1614,3.9736,2.1327,2.1382};
double lamdaQCD[7] = {0.0312,0.1377,0.0604,0.2452,0.4394,0.1377,0.1163};
double N[7] = {0.0,0.0,0.0,0.7838,1.0,1.0,1.0};

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


inline double Q2Tilda( const double Q2 )
{
    double res1 = Q2*Log( (sqr(lamda2())+ Q2)/sqr(lamdaQCD()) ) / Log( sql(lamda2())/sqr(lamdaQCD()) );
    double res2 = Q2*Log( (sqr(lamdaD())+ Q2)/sqr(lamdaQCD()) ) / Log( sql(lamdaD())/sqr(lamdaQCD()) );
    double res = res1;
    if (currentModel == 2)
    res = res2;
    return res;
}




////PhysRevC.69.068201 (2004)
double mu_n = -1.913;
double gama = 0.515;
double m_ro9 = 0.776;
double m_w = 0.783;
double m_phi = 1.019;
//double beta_ro = 0.512;
//double beta_w = 1.129;
//double beta_phi = -0.263;
//double alfa_ro = 2.675;
//double alfa_phi = -0.200;



inline double g ( const double Q2 ){
    double res = pow(1+gama*Q2,-2);
    return res;
}
inline double alph ( const double Q2,const double mass){
    double res = 2.0/PI*pow( (4.0*sqr(mass)+sqrt(Q2))/sqrt(Q2), 0.5);
    res = res*Log( (sqrt(4.0*sqr(mass)+Q2) + sqrt(Q2))/(2.0*mass) );
    return res;
}
inline double beta ( const double Q2,const double mass){
    double tau1 = Q2/4.0/sqr(mass);
    double res = sqrt(1.0+1.0/tau1);
    return res;
}
inline double massRatio ( const double Q2,const double mass){
    double res2 = sqr(mass)/(sqr(mass) + Q2);
    if (mass == m_ro)
    {
        double res1 = () / ( sqr(m_ro) +Q2 + (4.0*sqr(m_poin)+Q2)*0.112*alph(Q2,m_poin)/m_poin );
        return res1;
    }
    return res2;
}




// F1 / F2
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// Dipole
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double F1D ( const double Q2 ){
    double res = ( sqr(lamdaD())/(sqr(lamdaD()) + Q2Tilda) ) * ( sqr(lamda2())/(sqr(lamda2()) + Q2Tilda) );
    return res;
}
inline double F2D ( const double Q2 ){
    double res1 = sqr( sqr(lamdaD())/(sqr(lamdaD()) + Q2Tilda) ) * ( sqr(lamda2())/(sqr(lamda2()) + Q2Tilda) );
    double res2 = ( sqr(lamdaD())/(sqr(lamdaD()) + Q2Tilda) ) * sqr( sqr(lamda2())/(sqr(lamda2()) + Q2Tilda) );
    double res = res1;
    if (currentModel < 2)
    res = res2;
    
    return res;
    
}

// Modified Form Factor
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double F1_10 ( const double Q2 ){
    double p0= 8.73123;//   +/-   4.01161
    double p1=-130.677;//   +/-   103.246
    double p2= 1150.05;//   +/-   1119.3
    double p3=-6348.41;//   +/-   6723.79
    double p4= 22536.4;//   +/-   24736
    double p5=-51965.3;//   +/-   57981.4
    double p6= 77093.3;//   +/-   86885.2
    double p7=-70797.1;//   +/-   80529.7
    double p8= 36553;//   +/-   42026
    double p9=-8102.82;//   +/-   9437.31
    double res = p0+p1*Q2+p2*pow(Q2,2)+p3*pow(Q2,3)+p4*pow(Q2,4)+p5*pow(Q2,5)+p6*pow(Q2,6)+p7*pow(Q2,7)+p8*pow(Q2,8)+p9*pow(Q2,9);
    return res;
}
inline double F2_10 ( const double Q2 ){
    double p0                        =       33.059;//   +/-   92.9574
    double p1                        =     -518.566;//   +/-   2049.21
    double p2                        =      3944.47;//   +/-   19578.3
    double p3                        =     -17964.7;//   +/-   106424
    double p4                        =      53062.8;//   +/-   362885
    double p5                        =      -104946;//   +/-   805450
    double p6                        =       138845;//   +/-   1.16465e+06
    double p7                        =      -118444;//   +/-   1.05887e+06
    double p8                        =      59063.1;//   +/-   549813
    double p9                        =       -13095;//   +/-   124354
    double res = p0+p1*Q2+p2*pow(Q2,2)+p3*pow(Q2,3)+p4*pow(Q2,4)+p5*pow(Q2,5)+p6*pow(Q2,6)+p7*pow(Q2,7)+p8*pow(Q2,8)+p9*pow(Q2,9);
    return res;
}



// R0
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double F1RO ( const double Q2 ){
    double res = ( sqr(lamda1())/(sqr(lamda1()) + Q2Tilda) ) * ( sqr(lamda2())/(sqr(lamda2()) + Q2Tilda) );
    return res;
    
}
inline double F2RO ( const double Q2 ){
    double res1 = sqr( sqr(lamda1())/(sqr(lamda1()) + Q2Tilda) ) * ( sqr(lamda2())/(sqr(lamda2()) + Q2Tilda) );
    double res2 = ( sqr(lamda1())/(sqr(lamda1()) + Q2Tilda) ) * sqr( sqr(lamda2())/(sqr(lamda2()) + Q2Tilda) );
    double res = res1;
    if (currentModel < 2)
    res = res2;
    return res;
    
}


// W....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double F1W ( const double Q2 ){
    double res = ( sqr(lamda1())/(sqr(lamda1()) + Q2Tilda) ) * ( sqr(lamda2())/(sqr(lamda2()) + Q2Tilda) );
    return res;
    
}
inline double F2W ( const double Q2 ){
    double res1 = sqr( sqr(lamda1())/(sqr(lamda1()) + Q2Tilda) ) * ( sqr(lamda2())/(sqr(lamda2()) + Q2Tilda) );
    double res2 = ( sqr(lamda1())/(sqr(lamda1()) + Q2Tilda) ) * sqr( sqr(lamda2())/(sqr(lamda2()) + Q2Tilda) );
    double res = res1;
    if (currentModel < 2)
    res = res2;
    return res;
    
}



// PHI...
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double F1PHI ( const double Q2 ){
    double res = F1RO(Q2)*pow((Q2/(sqr(lamda1())+Q2)),1.5);
    return res;
}
inline double F2PHI ( const double Q2 ){
    double ration1 = sqr(lamda1())/sqr(mu_phy())*(Q2+sqr(mu_phy())/(sqr(lamda1())+Q2));
    double res = F2RO(Q2)*pow(ration1,1.5);
    return res;
}

// V...
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double F1V ( const double Q2 ){
    double m_roRatio1 = sqr(m_ro)/(sqr(m_ro)+Q2);
    double m_roRatio2 = sqr(m_ro_prime)/(sqr(m_ro_prime)+Q2);
    double m_roRatio = m_roRatio1;
    if (currentModel > 2)
    {
        m_roRatio = m_roRatio2;
        double res1 = N[currentModel]*(1.0317 + 0.0875*pow(1+Q2/0.3176,-2))/(1+Q2/0.5496)*F1RO(Q2);
        res1 += g_roPrimeTof_roPrime()*m_roRatio*F1RO(Q2);
        res1 += (1.0-1.1192*N[currentModel]-g_roPrimeTof_roPrime())*F1D(Q2);
        return res1;
    }
    double res = g_roPrimeTof_roPrime()*m_roRatio*F1RO(Q2) + (1.0-g_roPrimeTof_roPrime())*F1D(Q2);
    return res;
}
inline double F2V ( const double Q2 ){
    double m_roRatio1 = sqr(m_ro)/(sqr(m_ro)+Q2);
    double m_roRatio2 = sqr(m_ro_prime)/(sqr(m_ro_prime)+Q2);
    double m_roRatio = m_roRatio1;
    if (currentModel > 2)
    {
        m_roRatio = m_roRatio2;
        double res1 = N[currentModel]*(5.7824 + 0.3907*pow(1+Q2/0.1422,-1))/(1+Q2/0.5362)*F2RO(Q2);
        res1 += k_ro_prime()*g_roPrimeTof_roPrime()*m_roRatio*F2RO(Q2);
        res1 += (k_v - 6.1731*N[currentModel]-k_ro_prime()*g_roPrimeTof_roPrime())*F2D(Q2);
        return res1;
    }
    
    double res = k_ro_prime()*g_roPrimeTof_roPrime()*m_roRatio*F2RO(Q2) + (k_v-k_ro_prime()*g_roPrimeTof_roPrime())*F2D(Q2);
    return res;
}


// S...
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double F1S ( const double Q2 ){
    double m_wRatio = sqr(m_w)/(sqr(m_w)+Q2);
    double m_phiRatio = sqr(m_phi)/(sqr(m_phi)+Q2);
    double res = g_wTof_w()*m_wRatio*F1W(Q2) + g_phiTof_phi()*m_phiRatio*F1PHI(Q2) + (1.0-g_wTof_w())*F1D(Q2);
    return res;
    
}
inline double F2S ( const double Q2 ){
    double m_wRatio = sqr(m_w)/(sqr(m_w)+Q2);
    double m_phiRatio = sqr(m_phi)/(sqr(m_phi)+Q2);
    double res = k_w()*g_wTof_w()*m_wRatio*F2W(Q2) + k_phi()*g_phiTof_phi()*m_phiRatio*F2PHI(Q2) + (k_s-k_w()*g_wTof_w()-k_phi()*g_phiTof_phi())*F2D(Q2);
    return res;
    
}



// S9...
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double F1S9 ( const double Q2 ){
    double res = 0.5*g(Q2)*(1-beta(Q2,m_w)-beta(Q2,m_phi)+beta(Q2,m_w)*massRatio(Q2,m_w)+beta(Q2,m_phi)*massRatio(Q2,m_phi));
    return res;
}
inline double F2S9 ( const double Q2 ){
    double res = 0.5*g(Q2)*( (mu_p + mu_n - 1.0 - alph(Q2,m_phi))*massRatio(Q2,m_w) +alph(Q2,m_phi)*massRatio(Q2,m_phi));
    return res;
}


// V9...
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double F1V9 ( const double Q2 ){
    double res = 0.5*g(Q2)*( 1-beta(Q2,m_ro9)+beta(Q2,m_ro9)*massRatio(Q2,m_ro9) );
    return res;
}
inline double F2V9 ( const double Q2 ){
    double res = 0.5*g(Q2)*( (mu_p - mu_n - 1.0 - alph(Q2,m_ro9))/(1+gama*Q2) +alph(Q2,m_ro9)*massRatio(Q2,m_ro9));
    return res;
}














// G(E) / G(M)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// Bernauer Form Factor
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double G_E6 ( const double q2 ){
    double Q2 = TMath::Abs(q2);
    return 1+(-3.3686*Q2)+(14.5606*pow(Q2,2))+(-88.1912*pow(Q2,3))+(453.6244*pow(Q2,4))+(-1638.7911*pow(Q2,5))+(3980.7174*pow(Q2,6))
    +(-6312.6333*pow(Q2,7))+(6222.3646*pow(Q2,8))+(-3443.2251*pow(Q2,9))+(814.4112*pow(Q2,10));
}
inline double G_M6 ( const double q2 ){
    double Q2 = TMath::Abs(q2);
    return ( 1+(-2.5952*Q2)+(1.0222*pow(Q2,2))+(23.4945*pow(Q2,3))+(-93.0372*pow(Q2,4))+(140.7984*pow(Q2,5))
            + (-0.3656*pow(Q2,6))+(-305.6759*pow(Q2,7))+(444.6251*pow(Q2,8))+(-273.6688*pow(Q2,9))+(64.5811*pow(Q2,10)) ) * 2.7928474;
}
inline double TanTh2O2 ( double q2 = -.4,double E){
    double Q2 = TMath::Abs(q2);
    double eprime = E - 0.5 * Q2 / m_proton;
    double sinTh2O2 = Q2 / ( 4 * E * eprime ); // Q2 = 4 E E' Sin(th/2)^2 ==> Sin(th/2)^2 = 4 E E' / Q2
    return sinTh2O2 / ( 1.0 - sinTh2O2 ); // ( tan(arcsin(x)) )^2 = x^2 / (1-x^2)
}


// Modified Form Factor
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double G_E10 ( const double q2 ){
    double Q2 = TMath::Abs(q2);
    double tau = Q2 / 4.0 / sqr(m_proton);
    double res = F1_10(Q2) - tau*F2_10(Q2);
    return res;
}
inline double G_M10 ( const double q2 ){
    double Q2 = TMath::Abs(q2);
    double res = F1_10(Q2) + F2_10(Q2);
    return res;
}



// 8(?)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double G_E8 ( const double q2 ){
    double Q2 = TMath::Abs(q2);
    double tau = Q2 / 4.0 / sqr(m_proton);
    double res = F1S(Q2) + F1V(Q2) - tau*(F2S(Q2) + F2V(Q2));
    return res;
}
inline double G_M8 ( const double q2 ){
    double Q2 = TMath::Abs(q2);
    double res = (F1S(Q2) + F1V(Q2)) + (F2S(Q2) + F2V(Q2));
    return res;
}


// 9(?)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double G_E9 ( const double q2 ){
    double Q2 = TMath::Abs(q2);
    double tau = Q2 / 4.0 / sqr(m_proton);
    double res = F1S9(Q2) + F1V9(Q2) - tau*(F2S9(Q2) + F2V9(Q2));
    return res;
}
inline double G_M9 ( const double q2 ){
    double Q2 = TMath::Abs(q2);
    double res = (F1S9(Q2) + F1V9(Q2)) + (F2S9(Q2) + F2V9(Q2));
    return res;
}








// Px / Pz ratios
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Bernauer Form Factor
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double Px_Pz_6 ( double q2 = -0.4 ,double E) {
    double Q2 = TMath::Abs(q2);
    double eprime = E - 0.5 * Q2 / m_proton;
    double res = -2 * m_proton / ( eprime + E ) / Sqrt(TanTh2O2(q2,E)) * G_E6(q2) / G_M6(q2); // 6 is bernauer
    if (-5 < res && res < -0.005)
    return res;
    return -1.0;
}


// 8(?)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double Px_Pz8 ( double q2 = -0.4 ,double E){
    double Q2 = TMath::Abs(q2);
    double eprime = E - 0.5 * Q2 / m_proton;
    double res = -2 * m_proton / ( eprime + E ) / Sqrt(TanTh2O2(q2,E)) * G_E8(q2) / G_M8(q2);
    if (-5 < res && res < -0.005)
    return res;
    return -1.0;
}


// Modified Form Factor
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double Px_Pz_10 ( double q2 = -0.4 ,double E) { // Modified Form Factor
    double Q2 = TMath::Abs(q2);
    double eprime = E - 0.5 * Q2 / m_proton;
    double res = -2 * m_proton / ( eprime + E ) / Sqrt(TanTh2O2(q2,E)) * G_E10(q2) / G_M10(q2); // 6 is bernauer
    if (-5 < res && res < -0.005)
    return res;
    return -1.0;
}


// super ratio
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double PxPzBernauer_To_PxPzModifiedFF ( double q2 = -0.4 ,double E){
    double Px_Pz_Bernauer   = Px_Pz_6(q2,E);
    double Px_Pz_ModifiedFF = Px_Pz_10(q2,E);
    
    double res = Px_Pz_Bernauer / Px_Pz_ModifiedFF;
    if (-50 < res && res < 50){
        Printf("q2=%f , E = %f, Px_Pz_Bernauer = %f , Px_Pz_ModifiedFF = %f, ratio = %f",q2,E,Px_Pz_Bernauer,Px_Pz_ModifiedFF,res);
        return res;
    }
    return -1.0;
}

/*
 double I0(double q2,double tanTh2O2)
 {
	double ge = G_E6(q2);
	double gm = G_M6(q2);
	
	double tau    = -q2 / 4.0 / sqr(m_proton);
	double epsilon = 1.0/(1.0+2.0*(1+tau)*tanTh2O2);
	double temp = sqr(ge)+ tau*sqr(gm)*(1.0+2.0*(1.0+tau)*tanTh2O2);
	return temp;
 }
 double Px(double q2,double E)
 {
	double Q2 = TMath::Abs(q2);
	double eprime = E - 0.5 * Q2 / m_proton;
	double tanTh2O2 = TanTh2O2(q2,E);
	
	double ge = G_E6(q2);
	double gm = G_M6(q2);
	
	double tau    = -q2 / 4.0 / sqr(m_proton);
	double epsilon = 1.0/(1.0+2.0*(1.0+tau)*tanTh2O2);
	double res = -2.0*sqrt(tau*(tau+1))*gm*ge*sqrt(tanTh2O2);
	res = res/I0(q2,tanTh2O2);
	if (-5 < res && res < -0.005)
 return res;
	return -1.0;
 }
 
 double Pz(double q2,double E)
 {
	double Q2 = TMath::Abs(q2);
	double eprime = E - 0.5 * Q2 / m_proton;
	double tanTh2O2 = TanTh2O2(q2,E);
	
	double ge = G_E6(q2);
	double gm = G_M6(q2);
 
	double tau    = -q2 / 4.0 / sqr(m_proton);
	double epsilon = 1.0/(1.0+2.0*(1.0+tau)*tanTh2O2);
	double res = ((E+eprime)/m_proton)*sqrt(tau*(tau+1.0))*gm*gm*tanTh2O2;
	res = res/I0(q2,tanTh2O2);
	if (-5 < -res && -res < -0.005)
 return res;
	return 1.0;
 }
 */
void Bernauer()
{
    /*
     double q2=-0.4;
     double E=0.855;
     cout << "E=0.855 q2= -0.40 Px/Pz=" << Px_Pz(-0.40,E) << " Px=" << Px(-0.40,E) << " Pz=" << Pz(-0.40,E) << endl;
     cout << "E=0.855 q2= -0.45 Px/Pz=" << Px_Pz(-0.45,E) << " Px=" << Px(-0.45,E) << " Pz=" << Pz(-0.45,E) << endl;
     cout << "E=0.855 q2= -0.50 Px/Pz=" << Px_Pz(-0.50,E) << " Px=" << Px(-0.50,E) << " Pz=" << Pz(-0.50,E) << endl;
     
     cout << "E=0.600 q2= -0.40 Px/Pz=" << Px_Pz(-0.40,0.600) << " Px=" << Px(-0.40,0.600) << " Pz=" << Pz(-0.40,0.600) << endl;
     
     cout << "E=0.650 q2= -0.50 Px/Pz=" << Px_Pz(-0.50,0.650) << " Px=" << Px(-0.50,0.650) << " Pz=" << Pz(-0.50,0.650) << endl;
     */
    
    //cout << "Px/Pz= " << Px_Pz(-1.0*0.4,0.600)<< endl;
    //cout << "Px/Pz= " << Px_Pz8(-1.0*0.4,0.600)<< endl;
    //currentModel = 6;
    //TF1 *fa1 = new TF1("fa1","Px_Pz(-1.0*x,0.600)/Px_Pz8(-1.0*x,0.600)",0.1,0.5);
    //TF1 *fa1 = new TF1("fa1","G_M8(-1.0*x)/G_M6(-1.0*x)",0.1,0.5); 
    //fa1->Draw(); 	
    //TF1 *fa2 = new TF1("fa2","Px_Pz(-1.0*x,0.600)/Px_Pz8(-1.0*x,0.600)",0.1,0.5); 
    //fa2->Draw(); 	
    //TF1 *fa2 = new TF1("fa2","2.79284739*G_E8(-1.0*x)/G_M8(-1.0*x)",0.1,10.0); 
    //fa2->Draw(); 	

    
    TCanvas * c = new TCanvas("PxPzBernauer_To_PxPzModifiedFF","PxPzBernauer_To_PxPzModifiedFF",1000,1000);
    c -> Divide(3,3);
    
    float MinQ2 = 0.1, MaxQ2 = 0.9 , MinE = 0.3 , MaxE = 1.0;
    
    c -> cd(2);
    TF1 *fG_E6 = new TF1("fG_E6","G_E6(x)",MinQ2,MaxQ2);
    fG_E6->SetTitle("Bernauer G(E) (black) / G(M) (red)");
    fG_E6->GetXaxis()->SetTitle("Q2");
    fG_E6->SetLineColor(1);
    fG_E6->Draw();
    TF1 *fG_M6 = new TF1("fG_M6","G_M6(x)",MinQ2,MaxQ2);
    fG_M6->SetLineColor(2);
    fG_M6->Draw("same");
    
    c -> cd(4);
    TF1 *f1_10 = new TF1("fa1","F1_10(x)",MinQ2,MaxQ2);
    f1_10->SetTitle("Modified FF F1 (black) / F2 (red)");
    f1_10->GetXaxis()->SetTitle("Q2");
    f1_10->SetLineColor(1);
    f1_10->Draw();
    TF1 *f2_10 = new TF1("fa2","F2_10(x)",MinQ2,MaxQ2);
    f2_10->SetLineColor(2);
    f2_10->Draw("same");
    
    c -> cd(5);
    TF1 *fG_E10 = new TF1("fG_E10","G_E10(x)",MinQ2,MaxQ2);
    fG_E10->SetTitle("Modified FF G(E) (black) / G(M) (red)");
    fG_E10->GetXaxis()->SetTitle("Q2");
    fG_E10->SetLineColor(1);
    fG_E10->Draw();
    TF1 *fG_M10 = new TF1("fG_M10","G_M10(x)",MinQ2,MaxQ2);
    fG_M10->SetLineColor(2);
    fG_M10->Draw("same");
    
    c -> cd(3);
    TF2 *PxPzBernauer = new TF2("PxPzBernauer","Px_Pz_6(x,y)",MinQ2,MaxQ2,MinE,MaxE);
    PxPzBernauer->SetNpx(200);
    PxPzBernauer->SetNpy(200);
    PxPzBernauer->SetTitle("Bernauer Px/Pz ratio");
    PxPzBernauer->GetXaxis()->SetTitle("Q ^{2} (GeV/c) ^{2}");
    PxPzBernauer->GetYaxis()->SetTitle("E (GeV)");
    PxPzBernauer->Draw("colz");
    
    c -> cd(6);
    TF2 *PxPzModified = new TF2("PxPzBernauer","Px_Pz_10(x,y)",MinQ2,MaxQ2,MinE,MaxE);
    PxPzModified->SetNpx(200);
    PxPzModified->SetNpy(200);
    PxPzModified->SetTitle("Modified FF Px/Pz ratio");
    PxPzModified->GetXaxis()->SetTitle("Q ^{2} (GeV/c) ^{2}");
    PxPzModified->GetYaxis()->SetTitle("E (GeV)");
    PxPzModified->Draw("colz");
    
    
    c -> cd(9);
    TF2 *PxPzBernauer_To_PxPzModifiedFF = new TF2("PxPzBernauer_To_PxPzModifiedFF","PxPzBernauer_To_PxPzModifiedFF(x,y)",MinQ2,MaxQ2,MinE,MaxE);
    PxPzBernauer_To_PxPzModifiedFF->SetNpx(200);
    PxPzBernauer_To_PxPzModifiedFF->SetNpy(200);
    PxPzBernauer_To_PxPzModifiedFF->SetTitle("Bernauer Px/Pz / Modified FF Px/Pz");
    PxPzBernauer_To_PxPzModifiedFF->GetXaxis()->SetTitle("Q ^{2} (GeV/c) ^{2}");
    PxPzBernauer_To_PxPzModifiedFF->GetYaxis()->SetTitle("E (GeV)");
    PxPzBernauer_To_PxPzModifiedFF->Draw("colz");
}

