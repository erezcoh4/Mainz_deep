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
#include "TLine.h"
#include "TText.h"

#include "TGaxis.h"
//#define NOAH
#define CREATE_DAT_FILES
#define PKINEM
//#define K_PRIME_PLANE

using namespace std;
using namespace TMath;

inline double sqr ( double x ) {
  return x * x;
};

const double m_proton   = 0.938272046   ;
const double m_neutron  = 0.93956563    ;
const double m_deuteron = 1.87561339    ;
const double m_e        = 0.510998928e-3;
const double m_C12	= 11.174866     ;
const double PI         = 3.14159265    ;

//Bernauer FF
inline double G_E6 ( const double q2 ) {
  double QQ = -1.0*q2;// *pow(0.197327,2);
  return 1+(-3.3686*QQ)+(14.5606*pow(QQ,2))+(-88.1912*pow(QQ,3))+(453.6244*pow(QQ,4))+(-1638.7911*pow(QQ,5))+(3980.7174*pow(QQ,6))
    +(-6312.6333*pow(QQ,7))+(6222.3646*pow(QQ,8))+(-3443.2251*pow(QQ,9))+(814.4112*pow(QQ,10));
}
inline double G_M6 ( const double q2 ) {
  double QQ = -1.0*q2;// *pow(0.197327,2);
  return ( 1+(-2.5952*QQ)+(1.0222*pow(QQ,2))+(23.4945*pow(QQ,3))+(-93.0372*pow(QQ,4))+(140.7984*pow(QQ,5))
    + (-0.3656*pow(QQ,6))+(-305.6759*pow(QQ,7))+(444.6251*pow(QQ,8))+(-273.6688*pow(QQ,9))+(64.5811*pow(QQ,10)) ) * 2.7928474;
}

inline double TanTh2O2 ( double q2 = -.4 ) {
  double Q2 = TMath::Abs(q2);
  double E = ( Q2 < 0.29 ) ? 0.630 : 0.600;
  double eprime = E - 0.5 * Q2 / m_proton;
  double sinTh2O2 = Q2 / ( 4 * E * eprime ); // Q2 = 4 E E' Sin(th/2)^2 ==> Sin(th/2)^2 = 4 E E' / Q2
  return sinTh2O2 / ( 1.0 - sinTh2O2 ); // ( tan(arcsin(x)) )^2 = x^2 / (1-x^2)
}

inline double Px_Pz ( double q2 = -.4 ) {
  double Q2 = TMath::Abs(q2);
  double E = ( Q2 < 0.29 ) ? 0.630 : 0.600;
  double eprime = E - 0.5 * Q2 / m_proton;
//   double sinTh2O2 = sqrt ( Q2 / ( 4 * E * eprime ) );
//   double tanTh2O2 = sinTh2O2 / sqrt ( 1 - sinTh2O2 * sinTh2O2 );
  double res = -2 * m_proton / ( eprime + E ) / Sqrt(TanTh2O2(q2)) * G_E6(q2) / G_M6(q2);
   if (-5 < res && res < -0.005)
    return res;
   return -1.0;
}

double I0(double q2,double tanTh2O2)
{
	
	double ge = G_E6(q2);
	double gm = G_M6(q2);
	
	double tau    = -q2 / 4.0 / sqr(m_proton);
	double epsilon = 1.0/(1.0+2.0*(1+tau)*tanTh2O2);
	double temp = sqr(ge)+ tau*sqr(gm)*(1.0+2.0*(1.0+tau)*tanTh2O2);
	return temp;
}
double Px(double q2)
{	
	double Q2 = TMath::Abs(q2);
	double E = ( Q2 < 0.29 ) ? 0.630 : 0.600;
	double eprime = E - 0.5 * Q2 / m_proton;
	double tanTh2O2 = TanTh2O2(q2);
	
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

double Pz(double q2)
{		
	double Q2 = TMath::Abs(q2);
	double E = ( Q2 < 0.29 ) ? 0.630 : 0.600;
	double eprime = E - 0.5 * Q2 / m_proton;
	double tanTh2O2 = TanTh2O2(q2);
	
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

//virtuality of a proton in a deteron
inline double virtF ( double p ) {
  const double p2 = 1e-6 * sqr(p);
  const double emiss = m_deuteron - sqrt ( p2 + sqr(m_neutron) );
  return sqr(emiss) - p2 - sqr(m_proton);
}

TRandom2 rnd;
double virtCF ( double p, int shell ) {
  const double bmmissm[3]={10.28008,10.25143,10.25735};
  const double bmmisss[3]={0.01171,0.00081,0.00314};
  const double p2 = 1e-6 * sqr(p);  
  const double emiss = m_C12 - sqrt ( p2 + sqr(rnd.Gaus(bmmissm[shell],bmmisss[shell])) );
  return sqr(emiss) - p2 - sqr(m_proton);
}

inline double virtPrimeF ( double p ) {
  const double p2 = 1e-6 * sqr(p);
  const double emiss = m_deuteron/2;
  return sqr(emiss) - p2 - sqr(m_proton);
}

double virtPrimeCF ( double p, int shell ) {
  //const double bmmissm[3]={10.28008,10.25143,10.25735};
  //const double bmmisss[3]={0.01171,0.00081,0.00314};
  const double p2 = 1e-6 * sqr(p);  
  const double emiss = m_C12/12;
  return sqr(emiss) - p2 - sqr(m_proton);
}




const int nSets = 5;
TString setNames[nSets] = {"All", "K30", "K18", "G90", "G91A"};
#ifdef NOAH
TCut cuts = "q2<0";
#else
TCut cutx = "-1 < Px0 && Px0 < 0";
TCut cuty = "-1 < Py0 && Py0 < 1";
TCut cutz = " 0 < Pz0 && Pz0 < 1";
TCut cuts = cutx && cuty && cutz;
#endif

void findBadEvents ( const int set = 3 ) {
  TString fileName = TString("~/deep/")+setNames[set-1];
#ifdef NOAH
  fileName += "_NOAH.root";
#else
  fileName += "_D/his_2015_11_09_K_plane/Pol_A_12_6.root";//his_2015_11_09_K_plane
#endif
  TFile *file1 = new TFile(fileName);
  TTree *tree = (TTree*)file1->Get("ntuple");
  ((TTreePlayer*)(tree->GetPlayer()))->SetScanRedirect(true);
  ((TTreePlayer*)(tree->GetPlayer()))->SetScanFileName(setNames[set]+".bad");
  tree->Scan("*",!cuts);
}

  

TString like ( const int set = 0, const unsigned nBins = 0, 
	       TString var = "TMath::Sign(virtF(Pm),1+q2/(2*0.938272046*(E-Eprime/1000.0)))", 
	       bool divide = false, bool smallUnc = false ) {
  TString buff = "";
  
  TString fileName = TString("~/deep/")+setNames[set];
#ifdef NOAH
  fileName += "_NOAH.root";
#else
  if (set)
    fileName += "_D/his_2015_11_09_K_plane/Pol_A_12_6.root";//his_2015_11_09_K_plane
  else
    fileName += "DeepSets.root";
#endif
  TFile *file1 = new TFile(fileName);
  TTree *tree = (TTree*)file1->Get("ntuple");
  
  TCut binCut[nBins+1];
  binCut[0]=cuts;
  if (nBins) {
    const int nEntries = tree->Draw(var.Data(), cuts, "goff");
    Double_t *varVec = tree->GetV1();
    Int_t *index = new Int_t[nEntries];
    TMath::Sort(nEntries, varVec, index, false);
    double binX[nBins+1];
    binX[0] = varVec[index[0]];
    for (unsigned bin = 1; bin <= nBins; bin++) {
      binX[bin] = varVec[index[ nEntries * bin / nBins - 1 ]];
      binCut[bin]  = cuts;
      binCut[bin] += Form ( "%f<%s", binX[bin-1], var.Data() );
      binCut[bin] += Form ( "%f>%s", binX[bin],   var.Data() );
    }
  }

	
  TString cx = "cos(euler1)", sx = "sin(euler1)";
  TString cy = "cos(euler2)", sy = "sin(euler2)";
  TString cz = "cos(euler3)", sz = "sin(euler3)";
  TString s11 = TString("(") +cy+"*"+cz+")"; // cy*cz
  TString s12 = TString("(-")+cy+"*"+sz+")"; //-cy*sz
  TString s13 = TString("(-")+sy+")";        //-sy
  TString s21 = TString("(-")+sx+"*"+sy+"*"+cz+"+"+cx+"*"+sz+")"; //-sx*sy*cz+cx*sz
  TString s22 = TString("(") +sx+"*"+sy+"*"+sz+"+"+cx+"*"+cz+")"; // sx*sy*sz+cx*cz
  TString s23 = TString("(-")+sx+"*"+cy+")";                      //-sx*cy
  TString gXX = TString("A * ( -") + s11 + "*sin(phi) + " + s21 + "*cos(phi) )";
  TString gYX = TString("A * ( -") + s12 + "*sin(phi) + " + s23 + "*cos(phi) )";
  TString gX = gXX;// + "- sqrt(1-cosphi**2)/cosphi *" + gYX;
  TString gY = gYX;//TString("TMath::Abs(A)/0.85*( -") + s12 + "*sin(phi) + " + s23 + "*cos(phi) )";
  TString gZ = TString("A * ( -") + s13 + "*sin(phi) + " + s23 + "*cos(phi) )";
  
  const Long64_t MAX_DRAW = 10e5;

  for (unsigned bin=0; bin <= nBins; bin++) {
    buff += Form("%4i%4i", set, bin);
    double b1=0, b2=0, b3=0;
    double i11=0, i12=0, i13=0, i23=0, i22=0, i33=0;
    double X_sum=0, X_sum2=0, X, dX, M7=0, H=0, M7Py=0;
    double HPx=0,HPz=0,DPx=0,DPz=0;
    Long64_t n = tree->GetEntries(binCut[bin]), nPerReg;
    for (Long64_t first_event = 0; first_event < n; first_event+=MAX_DRAW) {
      if (divide)//Py0*
#ifdef K_PRIME_PLANE
	nPerReg = tree->Draw((gX+"*(Px):"+gY+":"+gZ+"*Pz:"+var+":1:1:1").Data(), binCut[bin], "para goff", MAX_DRAW, first_event);
#else
	nPerReg = tree->Draw((gX+"*(PxRot):"+gY+":"+gZ+"*PzRot:"+var+":1:1:1").Data(), binCut[bin], "para goff", MAX_DRAW, first_event);
#endif //K_PRIME_PLANE
      else
#ifdef K_PRIME_PLANE
	nPerReg = tree->Draw((gX+":"+gY+":"+gZ+":"+var+":Px/Pz:Px_Pz(q2):PyRot").Data(), binCut[bin], "para goff", MAX_DRAW, first_event);
#else
	nPerReg = tree->Draw((gX+":"+gY+":"+gZ+":"+var+":PxRot/PzRot:Px_Pz(q2):PyRot:PxRot:PzRot:Px(q2):Pz(q2)").Data(), binCut[bin], "para goff", MAX_DRAW, first_event);
#endif //K_PRIME_PLANE
      double *v1 = tree->GetV1();
      double *v2 = tree->GetV2();
      double *v3 = tree->GetV3();
      double *v4 = tree->GetV4();
      double *v5 = tree->GetVal(4);
      double *v6 = tree->GetVal(5);
      double *v7 = tree->GetVal(6);
      double *v8 = tree->GetVal(7);
      double *v9 = tree->GetVal(8);
      double *v10 = tree->GetVal(9);
      double *v11 = tree->GetVal(10);
      double *last = v1 + nPerReg;
      while ( v1 != last ) {
	b1     += (*v1);
	b2     += (*v2);
	b3     += (*v3);
	X_sum  += (*v4);
	M7     += (*v5);
	H      += (*v6);
	M7Py   += (*v7);

	DPx   += (*v8);
	DPz   += (*v9);
	HPx   += (*v10);
	HPz   += (*v11);

	i11    += (*v1) * (*v1);
	i12    += (*v1) * (*v2);
	i13    += (*v1) * (*v3);
	i22    += (*v2) * (*v2);
	i23    += (*v2) * (*v3);
	i33    += (*v3) * (*v3);
	X_sum2 += (*v4) * (*v4);
	++v1;
	++v2;
	++v3;
	++v4;
	++v5;
	++v6;
	++v7;
	++v8;
	++v9;
	++v10;
	++v11;
      }
    }
    
    double n_double = (double) n;
    X = X_sum/n_double;
    dX = sqrt ( X_sum2/n_double - X*X );
    M7   /= n_double;
    H    /= n_double;
    M7Py /= n_double;
    HPx /= n_double;
    HPz /= n_double;
    DPx /= n_double;
    DPz /= n_double;

    double det2D = i11*i33 - i13*i13;
    double u12D = b3*i13 - b1*i33;
    double u32D = b1*i13 - b3*i11;
    double  R2D = u12D/u32D;
    double dR2D = Abs(det2D) / sqr(u32D) * Sqrt( b1*b1*i33 - 2*b1*b3*i13 + b3*b3*i11 );
    double dRSmall2D = Abs(det2D/u32D) / Sqrt(i11);
    
    double det3D = Power(b3*(Power(i12,2) - i11*i22) + b2*(-(i12*i13) + i11*i23) + b1*(i13*i22 - i12*i23),2)
                   / Abs(Power(i13,2)*i22 - 2*i12*i13*i23 + Power(i12,2)*i33 + i11*(Power(i23,2) - i22*i33));
    double up3D  = Sqrt ( Power(b3,2)*i22*(-Power(i12,2) + i11*i22) + 2*b3*i22*(i13*(b2*i12 - b1*i22) + (-(b2*i11) + b1*i12)*i23) 
			+ i23*(2*b2*i13*(-(b2*i12) + b1*i22) + (Power(b2,2)*i11 - Power(b1,2)*i22)*i23) + Power(b2*i12 - b1*i22,2)*i33 );
    double  R3D1 = (b3*(i13*i22 - i12*i23) + b2*(-(i13*i23) + i12*i33) + b1*(Power(i23,2) - i22*i33));
    double  R3D3 = (b3*(Power(i12,2) - i11*i22) + b2*(-(i12*i13) + i11*i23) + b1*(i13*i22 - i12*i23));
    double  R3D = R3D1/R3D3;
    double dR3D = up3D / Abs(det3D);
    //cout << "Px=" << u12D/n_double << " Pz=" << u32D/n_double << " Px/Pz= "<< R3D << endl;
    
    double  PxRes = (b3*(-(i12*i23) + i22*i13) + b1*(Power(i23,2) - i22*i33) + b2*(-(i23*i13) + i12*i33))/(Power(i23,2)*i11 - 1*i12*i23*i13 + Power(i12,2)*i33 + i22*(Power(i13,2) - i11*i33));
    double dPxRes = Sqrt((Power(i23,2) - i22*i33)/(Power(i23,2)*i11 - 1*i12*i23*i13 + Power(i12,2)*i33 + i22*(Power(i13,2) - i11*i33)));
    double  PyRes = (b3*(-(i12*i13) + i11*i23) + b2*(Power(i13,2) - i11*i33) + b1*(-(i13*i23) + i12*i33))/(Power(i13,2)*i22 - 2*i12*i13*i23 + Power(i12,2)*i33 + i11*(Power(i23,2) - i22*i33));
    double dPyRes = Sqrt((Power(i13,2) - i11*i33)/(Power(i13,2)*i22 - 2*i12*i13*i23 + Power(i12,2)*i33 + i11*(Power(i23,2) - i22*i33)));
    double  PzRes = (b2*(-(i13*i12) + i11*i23) + b3*(Power(i12,2) - i11*i22) + b1*(-(i12*i23) + i13*i22))/(Power(i12,2)*i33 - 3*i13*i12*i23 + Power(i13,2)*i22 + i11*(Power(i23,2) - i33*i22));
    double dPzRes = Sqrt((Power(i12,2) - i11*i22)/(Power(i12,2)*i33 - 3*i13*i12*i23 + Power(i13,2)*i22 + i11*(Power(i23,2) - i33*i22)));

#ifdef PRINT_FOR_SPLINE
    const int nFields = 7;
#else
    const int nFields = 10;// + 4;
#endif
    int maxField = nFields;
    double fields3DWithPy[9] = {X, dX, R3D, smallUnc ? dRSmall2D : dR3D, M7, det3D/n/n/sqrt(n), H, PyRes, dPyRes};
    if (*fields3DWithPy>0.0)cout<<"";// NOTE: This line does not do anyting. It just prevents compilation warnings.
    double fields2DWithPy[10] = {X, dX, R2D, smallUnc ? dRSmall2D : dR2D, PyRes, dPyRes, M7, det2D/n_double/n_double, H, M7Py};
    if (*fields2DWithPy>0.0)cout<<"";// NOTE: This line does not do anyting. It just prevents compilation warnings.
#ifdef PRINT_FOR_SPLINE
    double fields[nFields] = {X, dX, -u12D/det2D, sqrt(i33/det2D), -u32D/det2D, sqrt(i11/det2D), -i12/sqrt(i11*i22)};
#else
    double fields[nFields] = {X, dX, R3D, dR3D, PyRes, dPyRes, M7, 1, H, M7Py};
    double fields22[nFields+4] = {X, dX, R3D, dR3D, PxRes, dPxRes, M7, PzRes, H, dPzRes,HPx,HPz,DPx,DPz};
    if (divide)
      maxField = 6;
#endif
    for(int fld=0; fld < maxField; fld++)
      buff += Form("%17.5f", fields[fld]);
    buff += "\n";
    TString buff1 = "";
    for(int fld1=0; fld1 < maxField+4; fld1++)
      buff1 += Form("%17.5f", fields22[fld1]);
    buff1 += "\n";
    //cout << buff1;

  }
  return buff;
}

TString kinem ( const int set = 0, const unsigned nBins = 0 ) {
  TString buff = "";
	
  const unsigned nPars = 7;
  TString pars[nPars] = { "-TMath::Sign(Pm,1+q2/(2*0.938272046*(E-Eprime/1000.0)))",
			  "-q2/(2*0.938272046*(E-Eprime/1000.0))", "Eprime", "theta", "theta1", 
			  "phi_pq-TMath::Sign(180,phi_pq)+180", "-q2" };
			  
  TString var= pars[0];
  
  TString fileName = TString("~/deep/")+setNames[set];
#ifdef NOAH
  fileName += "_NOAH.root";
#else
  if (set)
    fileName += "_D/his_2015_11_09_K_plane/Pol_A_12_6.root";//his_2015_11_09_K_plane
  else
    fileName += "DeepSets.root";
#endif
  TFile *file1 = new TFile(fileName);
  TTree *tree = (TTree*)file1->Get("ntuple");
  
  TCut binCut[nBins+1];
  binCut[0]=cuts;
  if (nBins) {
    const int nEntries = tree->Draw(var.Data(), cuts, "goff");
    Double_t *varVec = tree->GetV1();
    Int_t *index = new Int_t[nEntries];
    TMath::Sort(nEntries, varVec, index, false);
    double binX[nBins+1];
    binX[0] = varVec[index[0]];
    for (unsigned bin = 1; bin <= nBins; bin++) {
      binX[bin] = varVec[index[ nEntries * bin / nBins - 1 ]];
      binCut[bin]  = cuts;
      binCut[bin] += Form ( "%f<%s", binX[bin-1], var.Data() );
      binCut[bin] += Form ( "%f>%s", binX[bin],   var.Data() );
    }
  }
  
//   double parMeans[nBins][nPars];
//   double parRMSs [nBins][nPars];
  
  for (unsigned bin=0; bin <= nBins; bin++) {
    Long64_t n = tree->GetEntries(binCut[bin]);
    tree->SetEstimate(n);
    TH1F *h1;
    buff += Form("%4i%4i", set, bin);
    for (unsigned par = 0; par < nPars; par++) {
      tree->Draw((pars[par]+">>h1(1)").Data(), binCut[bin], "goff");
      h1 = (TH1F*)gDirectory->Get("h1");
//       parMeans[par] = h1 -> GetMean();
//       parRMSs [par] = h1 -> GetRMS ();
      buff += Form("%17.5f%17.5f", h1 -> GetMean(), h1 -> GetRMS());
    }      
    buff += "\n";
  }
  return buff;
}



void PlotPm() {
	TCanvas canvas = new TCanvas();
	gPad->SetLogz(1);
	gStyle->SetOptStat(0);
	gStyle->SetPalette(53);

	TChain* cn = new TChain("ntuple");
	cn->Add("K30_D/his_2015_11_09_K_plane/Pol_A_12_6.root");
	cn->Add("K18_D/his_2015_11_09_K_plane/Pol_A_12_6.root");
	cn->Add("G90_D/his_2015_11_09_K_plane/Pol_A_12_6.root");
	cn->Add("G91A_D/his_2015_11_09_K_plane/Pol_A_12_6.root");
	cn->GetEntries();
	//virtuality dependance
	//TH2F *h2 = new TH2F("h2", "Calculations as function of #nu;Sign(x_{B}-1) #times #nu;P_{x}/P_{z} M7/H",1600,-.1,.15,1600,.25,1.15);
	//cn->Draw("PxRot/PzRot/Px_Pz(q2):TMath::Sign(virtF(Pm),1+q2/(2*0.938272046*(E-Eprime/1000.0)))>>h2","Pm<300","colz");
	//h2r2 = h2->Rebin2D (8,8,"h2r2");
	//hpx = h2r2->ProfileX("_pfx",1,-1,"d hist l same");
	// gPad->SaveAs("nu.pdf")
	
	//Pm dependance
	//TH2F *h3 = new TH2F("h3", "Calculations as function of P_{m};P_{m};P_{x}/P_{z} M7/H",1600,-250,.15,1600,250,1.15);
	//cn->Draw("PxRot/PzRot/Px_Pz(q2):TMath::Sign(Pm,1+q2/(2*0.938272046*(E-Eprime/1000.0)))>>h3","abs(Pm)<400","colz");
	//gPad->SaveAs("Pm.pdf");
	
	//W
	//cn->Draw("sqrt(q2+2*1.875613*(E-Eprime/1000.0)+1.875613**2)-1.875613:TMath::Sign(virtF(Pm),1+q2/(2*0.938272046*(E-Eprime/1000.0)))>>h2(1000,-.15,.2,1000)","Pm<300","colz")
}

void Concat() {
  TFile f("AllDeepSets.root", "recreate");

  TChain* cn = new TChain("ntuple");
  cn->Add("K30_D/his_2015_11_09_K_plane/Pol_A_12_6.root");
  cn->Add("K18_D/his_2015_11_09_K_plane/Pol_A_12_6.root");
  cn->Add("G90_D/his_2015_11_09_K_plane/Pol_A_12_6.root");
  cn->Add("G91A_D/his_2015_11_09_K_plane/Pol_A_12_6.root");
  cn->GetEntries();
  cn->Write();
  f.Close();
}


//void compModels(TString xPar = "Pm"){//"virtF(Pm)"
//    const int nMod = 7;
//    TString tits[nMod] = {"PWBA (NR)", "DWIA (NR)", "PWBA (RC)", "DWIA+MEC (NR)", "DWIA+MEC+IC (NR)", "DWIA (RC)", "DWIA+MEC+IC (RC)"};
//    TCanvas *canvas = new TCanvas();
//    TChain *cn, *chain[nMod];
//    TProfile *prof[nMod];
//    TH1 *pr;
//    THStack *hs;
//    if (xPar == "virtF(Pm)")
//    hs = new THStack ("hsNiceName", ";Sign(x_{B}-1)#times#nu [GeV^{2}/c^{2}];(P_{t}/P_{l})_{A}/(P_{t}/P_{l})_{H}");
//    if (xPar == "Pm")
//    hs = new THStack ("hsNiceName", ";p_{m} [MeV/c];(P_{t}/P_{l})_{A}/(P_{t}/P_{l})_{H}");
//    for (int mod = 0; mod < nMod; mod++) {
//        chain[mod] = new TChain("ntuple");
//        cn = chain[mod];
//        for (int set = 2; set < nSets; set++) {
//            //      cn -> Add ( Form("/home/israely/deep/%s_D/his/Pol_A_12_%i.root", setNames[set-1].Data(), mod) );
//            cn -> Add ( Form("/home/israely/deep/%s_D/his_2015_11_09_K_plane/Pol_A_12_%i.root", setNames[set-1].Data(), mod) );
//        }
//        TH2F *h2;
//        if (xPar == "virtF(Pm)")
//        h2 = new TH2F(Form("h2%i",mod), "Calculations as function of #nu;Sign(x_{B}-1) #times #nu;P_{x}/P_{z} HA/H",200,-.1,.15,200,.25,1.15);
//        if (xPar == "Pm")
//        h2 = new TH2F(Form("h2%i",mod), "Calculations as function of #nu;Sign(x_{B}-1) #times #nu;P_{x}/P_{z} HA/H",100,-150,250,200,.25,1.15);
//        cn->Draw(Form("PxRot/PzRot/Px_Pz(q2):TMath::Sign(%s,1+q2/(2*0.938272046*(E-Eprime/1000.0)))>>h2%i",xPar.Data(),mod),"Pm<300","goff");
//        prof[mod] = h2 -> ProfileX("_pfx",1,-1,Form("d hist l %s",mod?"":"same"));
//        pr = /*static_cast<TH1*>*/ (prof[mod]);
//        pr -> SetLineColor(mod+1);
//        pr -> SetTitle(tits[mod].Data());
//        hs -> Add (pr);
//    }
//    
//    //  new TCanvas();
//    hs -> Draw("hist l nostack");
//    //   hs -> SaveAs ("curvesHStack.root");
//    gPad -> BuildLegend();
//    gPad -> Update();
//    TCanvas *canvasCol = new TCanvas("canvasCol","canvasCol");
//    cn->Draw("PxRot/PzRot/Px_Pz(q2):TMath::Sign(Pm,1+q2/(2*0.938272046*(E-Eprime/1000.0))):1e3*sqrt(q2+2*1.875613*(E-Eprime/1000.0)+1.875613**2)-1875.613>>h2temp","Pm<300","colz");
//}
//



// re-edited by Erez Cohen, April 12, 2016, only for virtuatily! (virtF(Pm))
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void compModels(TString xPar = "virtF(Pm)"){

    //{"PWBA (NR)", "DWIA (NR)", "PWBA (RC)", "DWIA+MEC (NR)", "DWIA+MEC+IC (NR)", "DWIA (RC)", "DWIA+MEC+IC (RC)"};
    std::vector<TString> models;
    models.push_back("PWBA (RC)"); // 3
    models.push_back("DWIA+MEC+IC (NR)"); // 5
    models.push_back("DWIA (RC)"); // 6
    models.push_back("DWIA+MEC+IC (RC)"); // 7

    TFile * OutFile = new TFile("/home/erez/VirtualityModels.root","recreate");
    TCanvas * canvas = new TCanvas("models","models");
    canvas -> Divide(2,1);
    
    TChain * chain = new TChain("ntuple");
    TProfile * prof[models.size()];
    TH1 *pr;
    THStack *hs = new THStack ("hsNiceName", "Relativistic corrections ; #nu [GeV ^{2}] ; (P_{x}/P_{z})_{A}/(P_{x}/P_{z})_{H}");

    TString xVar = "-TMath::Sign(virtF(Pm),1+q2/(2*0.938272046*(E-Eprime/1000.0)))";
    TString yVar = "PxRot/PzRot/Px_Pz(q2)";

    
    for (size_t i = 0; i < models.size(); i++) {
        Printf("processing %s...",models[i].Data());
        for (int set = 2; set < nSets; set++) {
            chain -> Add ( Form("/home/israely/deep/%s_D/his_2015_11_09_K_plane/Pol_A_12_%i.root", setNames[set-1].Data(), i) );
        }
        TH2F * h2 = new TH2F(Form("h2%i",i), "Calculations as function of #nu;Sign(x_{B}-1) #times #nu;P_{x}/P_{z} HA/H",200,-.13,.062,200,.25,1.15);
        chain -> Draw(Form("%s:%s>>h2%i",yVar.Data(),xVar.Data(),i),Form("Pm<300 && (%s<-0.011 || 0.014<%s)",xVar.Data(),xVar.Data()),"goff");
        prof[i] = h2 -> ProfileX("_pfx",1,-1,Form("d hist l %s",i?"":"same"));
        pr = prof[i];
        pr -> SetLineColor(i+1);
        pr -> SetTitle(models[i].Data());
        hs -> Add (pr);
        
        h2 -> Write();

    }
    // horizontal line at y = 1
    TLine * line = new TLine(-.13,1,0.062,1);
    line -> SetLineStyle( 2 );
    line -> SetLineColor( 1 );
    line -> SetLineWidth( 2 );
    line -> Draw();
    
    // add centeral axis
    TGaxis *CentralAxis = new TGaxis(0,0.57,0,1.3,"");
    CentralAxis -> SetName("CentralAxis");
    CentralAxis -> SetLabelSize(0.03);
    CentralAxis -> SetTextFont(72);
    CentralAxis -> SetLabelOffset(0.025);
    CentralAxis -> Draw();
    
    
    // add Pmiss texts
    TText * TLeft   = new TText(-0.06,0.6,"P_{miss} < 0");
    TLeft->Draw();
    TText * TRight  = new TText( 0.02,0.6,"P_{miss} > 0");
    TRight->Draw();

    hs -> Draw("hist l nostack");
    gPad -> BuildLegend();
    gPad -> SetTicks(1,1);
    gPad -> Update();
    
    
    hs -> Write();
    OutFile -> Write();
    OutFile -> Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......







void printKinem(){
  #ifdef CREATE_DAT_FILES
  ofstream *outDat;
  outDat = new ofstream("kinem.dat");
  #endif
  unsigned bins[nSets] = {0,6,0,7,0};
  for(int set = 1; set < nSets; set++) {
    TString resString = kinem(set, bins[set]);
    cout << resString;
    #ifdef CREATE_DAT_FILES
    (*outDat) << resString;
    #endif
  }
  cout << endl;
  #ifdef CREATE_DAT_FILES
  outDat->close();
  #endif
}



void LoopLike(){
#ifdef CREATE_DAT_FILES
  ofstream *outDat;
#endif
  const int nFigs = 3;
#ifdef PRINT_FOR_SPLINE
  unsigned bins[nSets]={0, 100,75,200,0};
//   unsigned bins[nSets]={0,0,0,4};
  int fg = 2;
#else
  unsigned bins[nSets] = {0,6,0,7,0};  
  int fg = 0;
#endif
  TString paras[nFigs] = {"-TMath::Sign(Pm,1+q2/(2*0.938272046*(E-Eprime/1000.0)))", "TMath::Sign(virtPrimeF(Pm),1+q2/(2*0.938272046*(E-Eprime/1000.0)))", "virtF(Pm)"};
  paras[2]=paras[1];
  bool dividers[nFigs] = {false, false, true};
  for (; fg < nFigs; fg++) {
#ifdef CREATE_DAT_FILES
    outDat = new ofstream(Form("fig%i.dat", fg + 1));
#endif
    for(int set = 1; set < nSets; set++) {
      TString resString = like(set, bins[set], paras[fg], dividers[fg], false);
      cout << resString;
#ifdef CREATE_DAT_FILES
      (*outDat) << resString;
#endif
    }
    cout << endl;
#ifdef CREATE_DAT_FILES
    outDat->close();
#endif
  }
}

void GuyMethodNOHA(){
#ifdef PKINEM
  printKinem();
#else
  LoopLike();
#endif
}
