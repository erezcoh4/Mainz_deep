#ifndef BERNAUER_CXX
#define BERNAUER_CXX

#include "Bernauer.h"

Bernauer::Bernauer(){
    //    //Table 1 Content
    currentModel = 0;
    double Tmp_g_roPrimeTof_roPrime[7] = {0.4466,0.0514,0.3223,0.1013,0.0808,0.0625,0.0636};
    double Tmp_k_ro_prime[7] = {4.3472,23.533,4.982,-15.87,-17.993,0.9397,-0.4175};
    double Tmp_g_wTof_w[7] = {0.4713,0.0588,0.344,0.6604,0.8038,0.8029,0.7918};
    double Tmp_k_w[7] = {21.762,18.934,40.661,8.847,4.0526,5.5225,5.1109};
    double Tmp_g_phiTof_phi[7] = {-0.8461,-0.5283,-0.9315,-0.4054,-0.2336,-0.3070,-0.3011};
    double Tmp_k_phi[7] = {11.849,1.2236,14.6805,13.6415,13.5963,14.4123,13.4385};
    double Tmp_mu_phy[7] = {1.1498,1.1670,1.1411,1.127,1.1218,1.2379,1.1915};
    double Tmp_lamda1[7] = {0.9006,0.5902,0.8956,0.89361,0.9295,0.9916,0.9660};
    double Tmp_lamdaD[7] = {1.7038,0.7273,1.7038,1.0454,1.2207,1.2589,1.3406};
    double Tmp_lamda2[7] = {1.1336,1.9368,0.9551,2.1614,3.9736,2.1327,2.1382};
    double Tmp_lamdaQCD[7] = {0.0312,0.1377,0.0604,0.2452,0.4394,0.1377,0.1163};
    double Tmp_N[7] = {0.0,0.0,0.0,0.7838,1.0,1.0,1.0};
    for (int i = 0; i < 7; i++) {
        g_roPrimeTof_roPrime[i] = Tmp_g_roPrimeTof_roPrime[i];
        k_ro_prime[i] = Tmp_k_ro_prime[i];
        g_wTof_w[i] = Tmp_g_wTof_w[i];
        k_w[i] = Tmp_k_w[i];
        g_phiTof_phi[i] = Tmp_g_phiTof_phi[i];
        k_phi[i] = Tmp_k_phi [i];
        mu_phy[i] = Tmp_mu_phy[i];
        lamda1[i] = Tmp_lamda1[i];
        lamdaD[i] = Tmp_lamdaD[i];
        lamda2[i] = Tmp_lamda2[i];
        lamdaQCD[i] = Tmp_lamdaQCD[i];
        N[i] = Tmp_N[i];

    }
}

//DipoleFit
double Bernauer::dipole_fit(double q2) const{
    return 1/((1 - q2/0.71)*(1 - q2/0.71));
}
double Bernauer::G_ED(double q2) const{
    return dipole_fit(q2);
}
double Bernauer::G_MD(double q2) const{
    return  mu_p * dipole_fit(q2);
}


inline double Bernauer::Q2Tilda( const double Q2 ){
    double res1[i] = Q2*Log( (sqr(lamda2())+ Q2)/sqr(lamdaQCD()) ) / Log( sql(lamda2())/sqr(lamdaQCD()) );
    double res2 = Q2*Log( (sqr(lamdaD())+ Q2)/sqr(lamdaQCD()) ) / Log( sql(lamdaD())/sqr(lamdaQCD()) );
    double res = res1;
    if (currentModel == 2)
    res = res2;
    return res;
}






inline double Bernauer::g ( const double Q2 ){
    double res = pow(1+gama*Q2,-2);
    return res;
}
inline double Bernauer::alph ( const double Q2,const double mass){
    double res = 2.0/PI*pow( (4.0*sqr(mass)+sqrt(Q2))/sqrt(Q2), 0.5);
    res = res*Log( (sqrt(4.0*sqr(mass)+Q2) + sqrt(Q2))/(2.0*mass) );
    return res;
}
inline double Bernauer::beta ( const double Q2,const double mass){
    double tau1 = Q2/4.0/sqr(mass);
    double res = sqrt(1.0+1.0/tau1);
    return res;
}
inline double Bernauer::massRatio ( const double Q2,const double mass){
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
inline double Bernauer::F1D ( const double Q2 ){
    double res = ( sqr(lamdaD())/(sqr(lamdaD()) + Q2Tilda) ) * ( sqr(lamda2())/(sqr(lamda2()) + Q2Tilda) );
    return res;
}
inline double Bernauer::F2D ( const double Q2 ){
    double res1 = sqr( sqr(lamdaD())/(sqr(lamdaD()) + Q2Tilda) ) * ( sqr(lamda2())/(sqr(lamda2()) + Q2Tilda) );
    double res2 = ( sqr(lamdaD())/(sqr(lamdaD()) + Q2Tilda) ) * sqr( sqr(lamda2())/(sqr(lamda2()) + Q2Tilda) );
    double res = res1;
    if (currentModel < 2)
    res = res2;
    
    return res;
    
}

// Modified Form Factor
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double Bernauer::F1_10 ( const double Q2 ){
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
inline double Bernauer::F2_10 ( const double Q2 ){
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
inline double Bernauer::F1RO ( const double Q2 ){
    double res = ( sqr(lamda1())/(sqr(lamda1()) + Q2Tilda) ) * ( sqr(lamda2())/(sqr(lamda2()) + Q2Tilda) );
    return res;
    
}
inline double Bernauer::F2RO ( const double Q2 ){
    double res1 = sqr( sqr(lamda1())/(sqr(lamda1()) + Q2Tilda) ) * ( sqr(lamda2())/(sqr(lamda2()) + Q2Tilda) );
    double res2 = ( sqr(lamda1())/(sqr(lamda1()) + Q2Tilda) ) * sqr( sqr(lamda2())/(sqr(lamda2()) + Q2Tilda) );
    double res = res1;
    if (currentModel < 2)
    res = res2;
    return res;
    
}


// W....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double Bernauer::F1W ( const double Q2 ){
    double res = ( sqr(lamda1())/(sqr(lamda1()) + Q2Tilda) ) * ( sqr(lamda2())/(sqr(lamda2()) + Q2Tilda) );
    return res;
    
}
inline double Bernauer::F2W ( const double Q2 ){
    double res1 = sqr( sqr(lamda1())/(sqr(lamda1()) + Q2Tilda) ) * ( sqr(lamda2())/(sqr(lamda2()) + Q2Tilda) );
    double res2 = ( sqr(lamda1())/(sqr(lamda1()) + Q2Tilda) ) * sqr( sqr(lamda2())/(sqr(lamda2()) + Q2Tilda) );
    double res = res1;
    if (currentModel < 2)
    res = res2;
    return res;
    
}



// PHI...
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double Bernauer::F1PHI ( const double Q2 ){
    double res = F1RO(Q2)*pow((Q2/(sqr(lamda1())+Q2)),1.5);
    return res;
}
inline double Bernauer::F2PHI ( const double Q2 ){
    double ration1 = sqr(lamda1())/sqr(mu_phy())*(Q2+sqr(mu_phy())/(sqr(lamda1())+Q2));
    double res = F2RO(Q2)*pow(ration1,1.5);
    return res;
}

// V...
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double Bernauer::F1V ( const double Q2 ){
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
inline double Bernauer::F2V ( const double Q2 ){
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
inline double Bernauer::F1S ( const double Q2 ){
    double m_wRatio = sqr(m_w)/(sqr(m_w)+Q2);
    double m_phiRatio = sqr(m_phi)/(sqr(m_phi)+Q2);
    double res = g_wTof_w()*m_wRatio*F1W(Q2) + g_phiTof_phi()*m_phiRatio*F1PHI(Q2) + (1.0-g_wTof_w())*F1D(Q2);
    return res;
    
}
inline double Bernauer::F2S ( const double Q2 ){
    double m_wRatio = sqr(m_w)/(sqr(m_w)+Q2);
    double m_phiRatio = sqr(m_phi)/(sqr(m_phi)+Q2);
    double res = k_w()*g_wTof_w()*m_wRatio*F2W(Q2) + k_phi()*g_phiTof_phi()*m_phiRatio*F2PHI(Q2) + (k_s-k_w()*g_wTof_w()-k_phi()*g_phiTof_phi())*F2D(Q2);
    return res;
    
}



// S9...
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double Bernauer::F1S9 ( const double Q2 ){
    double res = 0.5*g(Q2)*(1-beta(Q2,m_w)-beta(Q2,m_phi)+beta(Q2,m_w)*massRatio(Q2,m_w)+beta(Q2,m_phi)*massRatio(Q2,m_phi));
    return res;
}
inline double Bernauer::F2S9 ( const double Q2 ){
    double res = 0.5*g(Q2)*( (mu_p + mu_n - 1.0 - alph(Q2,m_phi))*massRatio(Q2,m_w) +alph(Q2,m_phi)*massRatio(Q2,m_phi));
    return res;
}


// V9...
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double Bernauer::F1V9 ( const double Q2 ){
    double res = 0.5*g(Q2)*( 1-beta(Q2,m_ro9)+beta(Q2,m_ro9)*massRatio(Q2,m_ro9) );
    return res;
}
inline double Bernauer::F2V9 ( const double Q2 ){
    double res = 0.5*g(Q2)*( (mu_p - mu_n - 1.0 - alph(Q2,m_ro9))/(1+gama*Q2) +alph(Q2,m_ro9)*massRatio(Q2,m_ro9));
    return res;
}














// G(E) / G(M)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// Bernauer Form Factor
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double Bernauer::G_E6 ( const double q2 ){
    double QQ = -1.0*q2;// *pow(0.197327,2);
    return 1+(-3.3686*QQ)+(14.5606*pow(QQ,2))+(-88.1912*pow(QQ,3))+(453.6244*pow(QQ,4))+(-1638.7911*pow(QQ,5))+(3980.7174*pow(QQ,6))
    +(-6312.6333*pow(QQ,7))+(6222.3646*pow(QQ,8))+(-3443.2251*pow(QQ,9))+(814.4112*pow(QQ,10));
}
inline double Bernauer::G_M6 ( const double q2 ){
    double QQ = -1.0*q2;// *pow(0.197327,2);
    return ( 1+(-2.5952*QQ)+(1.0222*pow(QQ,2))+(23.4945*pow(QQ,3))+(-93.0372*pow(QQ,4))+(140.7984*pow(QQ,5))
            + (-0.3656*pow(QQ,6))+(-305.6759*pow(QQ,7))+(444.6251*pow(QQ,8))+(-273.6688*pow(QQ,9))+(64.5811*pow(QQ,10)) ) * 2.7928474;
}
inline double Bernauer::TanTh2O2 ( double q2 ,double E){
    double Q2 = TMath::Abs(q2);
    double eprime = E - 0.5 * Q2 / m_proton;
    double sinTh2O2 = Q2 / ( 4 * E * eprime ); // Q2 = 4 E E' Sin(th/2)^2 ==> Sin(th/2)^2 = 4 E E' / Q2
    return sinTh2O2 / ( 1.0 - sinTh2O2 ); // ( tan(arcsin(x)) )^2 = x^2 / (1-x^2)
}


// Modified Form Factor
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double Bernauer::G_E10 ( const double q2 ){
    double Q2 = TMath::Abs(q2);
    double tau = Q2 / 4.0 / sqr(m_proton);
    double res = F1_10(Q2) - tau*F2_10(Q2);
    return res;
}
inline double Bernauer::G_M10 ( const double q2 ){
    double Q2 = TMath::Abs(q2);
    double res = F1_10(Q2) + F2_10(Q2);
    return res;
}



// 8(?)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double Bernauer::G_E8 ( const double q2 ){
    double Q2 = TMath::Abs(q2);
    double tau = Q2 / 4.0 / sqr(m_proton);
    double res = F1S(Q2) + F1V(Q2) - tau*(F2S(Q2) + F2V(Q2));
    return res;
}
inline doublBernauer::G_M8 ( const double q2 ){
    double Q2 = TMath::Abs(q2);
    double res = (F1S(Q2) + F1V(Q2)) + (F2S(Q2) + F2V(Q2));
    return res;
}


// 9(?)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double Bernauer::G_E9 ( const double q2 ){
    double Q2 = TMath::Abs(q2);
    double tau = Q2 / 4.0 / sqr(m_proton);
    double res = F1S9(Q2) + F1V9(Q2) - tau*(F2S9(Q2) + F2V9(Q2));
    return res;
}
inline double Bernauer::G_M9 ( const double q2 ){
    double Q2 = TMath::Abs(q2);
    double res = (F1S9(Q2) + F1V9(Q2)) + (F2S9(Q2) + F2V9(Q2));
    return res;
}








// Px / Pz ratios
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Bernauer Form Factor
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double Bernauer::Px_Pz_6 ( double q2 ,double E) {
    double Q2 = TMath::Abs(q2);
    double eprime = E - 0.5 * Q2 / m_proton;
    double res = -2 * m_proton / ( eprime + E ) / Sqrt(TanTh2O2(q2,E)) * G_E6(q2) / G_M6(q2); // 6 is bernauer
    if (-5 < res && res < -0.005)
    return res;
    return -1.0;
}


// 8(?)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double Bernauer::Px_Pz8 ( double q2 ,double E){
    double Q2 = TMath::Abs(q2);
    double eprime = E - 0.5 * Q2 / m_proton;
    double res = -2 * m_proton / ( eprime + E ) / Sqrt(TanTh2O2(q2,E)) * G_E8(q2) / G_M8(q2);
    if (-5 < res && res < -0.005)
    return res;
    return -1.0;
}


// Modified Form Factor
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double Bernauer::Px_Pz_10 ( double q2 ,double E) { // Modified Form Factor
    double Q2 = TMath::Abs(q2);
    double eprime = E - 0.5 * Q2 / m_proton;
    double res = -2 * m_proton / ( eprime + E ) / Sqrt(TanTh2O2(q2,E)) * G_E10(q2) / G_M10(q2); // 6 is bernauer
    if (-5 < res && res < -0.005)
    return res;
    return -1.0;
}


// super ratio
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline double Bernauer::PxPzBernauer_To_PxPzModifiedFF ( double q2,double E){
    double Q2 = TMath::Abs(q2);
    double Px_Pz_Bernauer =
    double res = -2 * m_proton / ( eprime + E ) / Sqrt(TanTh2O2(q2,E)) * G_E10(q2) / G_M10(q2); // 6 is bernauer
    if (-5 < res && res < -0.005)
    return res;
    return -1.0;
}

//
//void Bernauer()
//{
//    TF1 *fa2 = new TF1("fa2","F1_10(x)",0.13,0.9); ///G_M10(-1.0*x)
//    fa2->Draw();
//}


#endif
