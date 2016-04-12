#include "TString.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TColor.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TPaveStats.h"
#include "TThread.h"
#include "TMatrixT.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TF1.h"
#include "TText.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPavesText.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TFrame.h"
#include "TAttMarker.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

//#define CURVES
// #define MOD4
//#define ROOT6
// #define ALLBINS
// #define VAN_ORDEN_COLS
#define C12 0

const double m_proton   = 0.938272046;
const double m_neutron  = 0.93956563;
const double m_deuteron = 1.87561339;
double virtF(double p){
    const double p2=1e-6*p*p;
    const double emiss=m_deuteron-sqrt(p2+m_neutron*m_neutron);
    return emiss*emiss-p2-m_proton*m_proton;
}


TMultiGraph *mg, *mgu, *mgd, *mg3, *mg4, *mg5, *mg6;
const int nMod=2*2,nSet=4,nGED=4,nJLab=3,nFigs=6,nC12=1, nnC12=C12*nC12,nSubPads=2,nAxs=2;
const int nTot=nMod+nSet+nJLab+nnC12,nMG=nFigs+nSubPads-1,nTrees=nFigs-1;
TTree *t1, *t2, *t3, *t4, *t4_C12;
TTree *t[nTrees]={t1,t2,t3,t4,t4_C12};
TGraph*gr[nTot][nFigs];
int ln[nSet][nFigs];
int ln_C12[nC12];
double*dat[nSet][nGED][nFigs];
double *dat_C12;
TCanvas *c1, *c2, *c3, *c4, *c5, *c6;
TCanvas *cc[nFigs]={c1,c2,c3,c4,c5,c6};
TLegend *leg1, *leg2, *leg21, *leg22, *legF;
TBox *box2;
TLatex*tl;

void drawFigSplit(){
    EColor cols[nTot]={kRed,kRed,kBlue,kBlue,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack};//,kGray}; .
    int styles[nTot]={3,3,28,28,23,21,20,22,33,25,25};//21
    EColor cols_C12[nC12] = {kBlue};
    int style_C12 = 29;
    TString tits1[nTot]={"Calculation [5]","","#splitline{Non-relativistic}{calculation [5]}",
        "","^{2}H, setup A","^{2}H, setup B","^{2}H, setup D","^{2}H, setup C","d","#alpha","#alpha"};
#ifdef VAN_ORDEN_COLS
    TString modPar22 = "DWBA/proton:v";
    //   TString modPar[nFigs][nMod/2]={{"M7:Pm","DWBA:Pm"},{"M7/H:v",modPar22},{"1:1","1:1"}};
    TString tits2[nTot]={"DW+MEC+IC+RC HA","","DWBA WVO","",
        "","","^{2}H, This work","",
        "^{2}H, JLab Q^{2}=1.0 [GeV^{2}/c^{2}]","^{4}He, JLab Q^{2}=0.8 [GeV^{2}/c^{2}]","^{4}He, JLab Q^{2}=0.8 [GeV^{2}/c^{2}]"};
#else
    TString tits2[nTot]={"^{2}H Calculation","","Non-relativistic calculation","",
        "","","^{2}H, This work","",
        "^{2}H, JLab","^{4}He, JLab","^{4}He, JLab"
    };//"^{2}H, JLab Q^{2}=1.0 [GeV^{2}/c^{2}]","^{4}He, JLab Q^{2}=0.8 [GeV^{2}/c^{2}]"
    TString modPar22 = "M4/H:v";
#ifdef MOD4
    tits1[0]="Full "+tits1[0];
    tits2[0]="Full "+tits2[0];
#endif // MOD4
#endif // VAN_ORDEN_COLS
    if (C12)
    for (int i = 0; i < nnC12; i++) {
        cols[nTot-nnC12+i] = cols_C12[i];
        styles[nTot-nnC12+i] = style_C12;
        tits2[nTot-nnC12+i] = "";
    }
    TString modPar[nFigs][nMod/2]={{"M7:Pm","M4:Pm"},{"M7/H:-v",modPar22},{"1:1","1:1"},{"1:1","1:1"},{"M7Py:v","1:1"},{"1:1","1:1"}};
    TString expPar[nFigs]={"Pm:dPm:dat:ddat","-v:dv:dat/H:-1.0*ddat/H","-v:dv:dat:1.0*ddat","-v:dv:dat:1.0*ddat","-v:dv:Py:dPy","-v:dv:Py:dPy"};
    
    
    //Create a canvas for the temporary graphs
    c1=new TCanvas("c1","c1");
    //This multigraph will hold both the data and the calculations
    mg  = new TMultiGraph("mg", ";#nu [GeV^{2}/c^{2}];(P_{x}/P_{z})_{A}/(P_{x}/P_{z})_{H}");
    mg3 = new TMultiGraph("mg3",";#nu [GeV^{2}/c^{2}];(P_{x}/P_{z})^{exp}_{d}/(P_{x}/P_{z})^{calc}_{d}");
    mg4 = new TMultiGraph("mg4",";#nu [GeV^{2}/c^{2}];(P_{x}/P_{z})^{exp}_{d}/(P_{x}/P_{z})^{calc}_{d}");
    mg5 = new TMultiGraph("mg5",";#nu [GeV^{2}/c^{2}];P_{n}");
    mg6 = new TMultiGraph("mg6",";#nu [GeV^{2}/c^{2}];(P_{n})^{exp}_{d}/(P_{n})^{calc}_{d}");
    mgu = new TMultiGraph("mgu",";;P_{x}/P_{z}");
    mgd = new TMultiGraph("mgd",";p_{miss} [MeV/c]");
    TString fig1Tits[nSubPads]={"(a)   Q^{2} = 0.40 [GeV^{2}/c^{2}]","(b)   Q^{2} = 0.18 [GeV^{2}/c^{2}]"};
#ifdef VAN_ORDEN_COLS
    TString treeKey2 = "set/I:bin:v/F:dv:dat:ddat:M7:M4:H:PWIA:PWBA:DWBA:proton";
#else
    TString treeKey2 = "set/I:bin:v/F:dv:dat:ddat:Py:dPy:M7:M4:H:M7Py";
#endif
    TString treeKey[nTrees]={"set/I:bin:Pm/F:dPm:dat:ddat:Py:dPy:M7:M4:H:M7Py",treeKey2,"set/I:bin:v/F:dv:dat:ddat:Py:dPy","set/I:bin:v/F:dv:dat:ddat","set/I:bin:v/F:dv:dat:ddat"};
    //Read input files to tree
    for ( int tt=0; tt<nTrees; tt++ ) {
        if (tt == nTrees-1) {
            t[tt]=new TTree("t4_C12","t4_C12");
            t[tt]->ReadFile("fig4_C.dat",treeKey[tt]);
        }
        else {
            t[tt]=new TTree(Form("t%i",tt+1),Form("t%i",tt+1));
            t[tt]->ReadFile(Form("fig%i.dat",tt+1),treeKey[tt]);
        }
        //     t[tt]->Scan("*");
    }
    cerr << "curr t4_C12\t" << (t4_C12=t[nTrees-1]) << endl;
    
    
    //====================================================
    // This section loops over the different experimental settings and creates a seperate graph for each
    //====================================================
    for ( int set = 1; set <= nSet; set++ ){
        for ( int fg = 0; fg < nFigs; fg++  ){
            int treeNum=fg;
            if (fg>3)
            treeNum -= 3;
            TCut setCut=Form("set==%i",set);
            TCut binCut;
            //if the setting is odd  (K30 and G90) draw the binned data
            //if the setting is even (K18 and G91) draw the single bin data
#ifndef ALLBINS
            if(!(set%2))
            binCut="bin==0";
            else
#endif
            binCut="bin>0";
            ln[set-1][fg]=t[treeNum]->Draw(expPar[fg].Data(),setCut+binCut,"goff");
            for (int clmn = 0; clmn < nGED; clmn++)
            dat[set-1][clmn][fg]=t[treeNum]->GetVal(clmn);
            gr[nMod+set-1][fg]=new TGraphErrors(ln[set-1][fg],dat[set-1][0][fg],dat[set-1][2][fg],dat[set-1][1][fg],dat[set-1][3][fg]);
        }
    }
    
    //====================================================
    // This section loops over the different calculation models and creates a seperate graph for each
    //====================================================
#ifdef ALLBINS
    TCut binCut="bin>0",subFigCut;
#else
    TCut binCut="(set%2)*(bin>0)+(1-set%2)*(bin==0)",subFigCut;
#endif
    for(int mod=0;mod<nMod;mod++){
        if(!(mod%2))
        subFigCut="set<3";
        else
        subFigCut="set>2";
        for(int fg=0;fg<nFigs;fg++){
            int treeNum=fg;
            if (fg>3)
            treeNum -= 3;
            t[treeNum]->Draw(modPar[fg][mod/2].Data(),binCut+subFigCut);
            gr[mod][fg]=(TGraph*)(gPad->GetPrimitive("Graph")->Clone());
        }
    }
    c1->Close();
    
    //====================================================
    //This section creates graphs for the carbon data
    //====================================================
    
    if(C12) {
        t4_C12 = t[nTrees-1];
        ln_C12[0]=t4_C12->Draw(expPar[2].Data(),"bin>0 && set<3","goff");
        gr[nMod+nSet+nJLab+0][1] = new TGraphErrors(ln_C12[0], t4_C12->GetV1(),t4_C12->GetV3(),t4_C12->GetV2(),t4_C12->GetV4());
    }
    
    //====================================================
    //This section calculates JLab d values
    //====================================================
    
    const int nJLabD=6;
    double  pmJLabD[nJLabD]={  -135.,  -170.,   -26.,   -56.,  26.,  57.};
    double  PRJLabD[nJLabD]={-0.485,-0.420,-0.545,-0.553,-0.520,-0.539};
    double dPRJLabD[nJLabD]={ 0.052, 0.050, 0.019, 0.025, 0.018, 0.019};
    double   vJLabD[nJLabD]={0.};
    double  dvJLabD[nJLabD]={0.};
    double  DRJLabD[nJLabD]={0.};
    double dDRJLabD[nJLabD]={0.};
    
    const double R_h = -0.537;//Hydrogen Px/Pz at Q2=1
    for(int point=0;point<nJLabD;point++){
        vJLabD[point]   = TMath::Sign(virtF(pmJLabD[point]),pmJLabD[point]);//virtF(pmJLabD[point])
        DRJLabD[point]  =  PRJLabD[point]/R_h;
        dDRJLabD[point] =-dPRJLabD[point]/R_h;
    }
    gr[nMod+nSet+0][1] = new TGraphErrors(nJLabD,vJLabD,DRJLabD,dvJLabD,dDRJLabD);
    
    //====================================================
    // This section creates a graph for JLab He values
    //====================================================
    
    //   const int nJLabHe=5;
    //   double x1[nJLabHe], y1[nJLabHe], ex1[nJLabHe], ey1[nJLabHe];
    //
    //   x1[0]=-0.0541785; y1[0]=0.867779; ex1[0]=0; ey1[0]=0.0534269;
    //   x1[1]=-0.0474135; y1[1]=0.917724; ex1[1]=0; ey1[1]=0.0327462;
    //   x1[2]=-0.0426928; y1[2]=0.893016; ex1[2]=0; ey1[2]=0.0244888;
    //   x1[3]=-0.0393226; y1[3]=0.887166; ex1[3]=0; ey1[3]=0.0188161;
    //   x1[4]=-0.0370113; y1[4]=0.921571; ex1[4]=0; ey1[4]=0.0194034;
    
    //  -114.257	-46.154	53.255	3.5505	0.8035505	0.0497045
    // -69.002	56.36	99.85	78.105	0.878105	0.021745
    // -33.264	42.16	84.319	63.2395	0.8632395	0.0210795
    // 33.357	114.498	158.871	136.6845	0.9366845	0.0221865
    // 68.51	103.848	151.329	127.5885	0.9275885	0.0237405
    // 112.988	79.439	222.337	150.888	0.950888	0.071449
    
    //   const int nJLabHe=6;
    //   double x1[nJLabHe], y1[nJLabHe], ex1[nJLabHe], ey1[nJLabHe];
    //   //NOTE minus corrected
    //   x1[0]= 114.257; y1[0]=0.8035505; ex1[0]=0; ey1[0]=0.0497045;
    //   x1[1]=  69.002; y1[1]=0.878105 ; ex1[1]=0; ey1[1]=0.021745 ;
    //   x1[2]=  33.264; y1[2]=0.8632395; ex1[2]=0; ey1[2]=0.0210795;
    //   x1[3]= -33.357; y1[3]=0.9366845; ex1[3]=0; ey1[3]=0.0221865;
    //   x1[4]= -68.510; y1[4]=0.9275885; ex1[4]=0; ey1[4]=0.0237405;
    //   x1[5]=-112.988; y1[5]=0.950888 ; ex1[5]=0; ey1[5]=0.071449 ;
    const int nJLabHe=6;
    double x1[nJLabHe] = { -0.0569448991032313, -0.0459520254010585, -0.0411057848090918,0.0411140013247896,0.0458623238666614,0.0565627565350318},
    ex1[nJLabHe] = {0.00574508723804499,0.00575736851809157,0.00576278684923525,0.00576277766065254,0.0057574687859793,0.00574551395312394},
    y1[nJLabHe] = {0.8035505,0.878105,0.8632395,0.9366845,0.9275885,0.950888},
    ey1[nJLabHe] = {0.0497045,0.021745,0.0210795,0.0221865,0.0237405,0.071449};
    
    gr[nMod+nSet+1][1] = new TGraphErrors(nJLabHe, x1, y1, ex1, ey1);
    
    
    
    //   const int nJLabHePy=2;
    //   double x2[nJLabHePy] ={TMath::Mean(nJLabHe,x1),TMath::Mean(nJLabHe,x1)};
    //   double ex2[nJLabHePy]={TMath::RMS (nJLabHe,x1),TMath::RMS (nJLabHe,x1)};
    //   double y2[nJLabHePy] ={-.0415,-.0415};
    //   double ey2[nJLabHePy]={ .0050, .0072};
    
    const int nJLabHePy=6;
    double x2[nJLabHePy] ={-0.052550677,-0.042059884,-0.037270806,-0.037238672,-0.041828492,-0.052441323};
    double ex2[nJLabHePy]={0.0};
    double y2[nJLabHePy] ={-0.0623179,-0.0506718,-0.0351829,-0.0426834,-0.0340251,-0.0396349};
    double ey2[nJLabHePy]={0.022805,0.010671,0.0104271,0.0103658,0.0117074,0.027561};
    
    gr[nMod+nSet+2][4] = new TGraphErrors(nJLabHePy, x2, y2, ex2, ey2);
    
    
    
    //====================================================
    // Global Settings
    //====================================================
    TMultiGraph*mg1[nSubPads] = {mgu,mgd};
    TMultiGraph*mgA[nMG]      = {mgu,mgd,mg,mg3,mg4,mg5,mg6};
    
    int figureNum = 1;
    int nLocalMax = nTot;//;-nnC12
    for(int i=nLocalMax-nJLab-nnC12;i<nLocalMax;i++){
        if (i==nMod+nSet+2)
        figureNum = 4;
        gr[i][figureNum]->SetTitle(tits2[i].Data());
        gr[i][figureNum]->SetMarkerColor(cols[i]);
        gr[i][figureNum]->SetMarkerStyle(styles[i]);
        gr[i][figureNum]->SetMarkerSize (2);
        mgA[figureNum+1] -> Add(gr[i][figureNum]);
        if (i==nMod+nSet+2)
        figureNum = 1;
    }
    
    
    
    for(int fg=0;fg<nFigs;fg++){
        for(int i=0;i<nTot-nJLab-nnC12;i++){
#ifndef MOD4
            //Skip the graphs of model 4
            if(1<i&&i<nMod)
            continue;
#endif
            //Skip the calculation graphs
            if(fg>1&&fg!=4&&i<nMod)
            continue;
            
            gr[i][fg]->SetTitle(tits2[i].Data());
            if ( fg == 3 ) {
                gr[i][fg] -> SetMarkerColor(kBlue);
                gr[i][fg] -> SetLineColor  (kBlue);
            }
            else
            gr[i][fg]->SetMarkerColor(cols[i]);
#ifdef ROOT6
            //makes the markers semi- see through
            //works only with ROOT v. 6+
            gr[i][fg]->SetMarkerColorAlpha(0.8);
#endif
            gr[i][fg]->SetMarkerStyle(styles[i]);
            gr[i][fg]->SetMarkerSize (2);
            if ( fg )
            mgA[fg+1] -> Add(gr[i][fg]);
            else if((i<nMod&&!(i%2))||(nMod<=i&&i<nMod+2)) //Selects the graphs to go to the top and bottom frame
            mgu-> Add(gr[i][0]);
            else
            mgd-> Add(gr[i][0]);
        }
    }
    
    //====================================================
    // Draw
    //====================================================
    
    TAxis *xaxis, *yaxis;
    tl=new TLatex();
    tl->SetTextAlign(12);
    double labelSizeFig1=.09;
    double labelOffSFig1=.01;
    double titleSizeFig1=.11;
    double titleOffSFig1=.55;
    double botMarg=.25;
    double topMarg=botMarg/2.;
    double padSizeFactor, sizeFactor[nSubPads]={1.1,1.0};
    gStyle->SetTitleSize(0.73);
    cc[0] = c1 = new TCanvas ( "c1", "c1" );
    c1->Divide(1,nSubPads,0,0);
    for(int pad=0;pad<nSubPads;pad++){
        c1->cd(pad+1);
        padSizeFactor=sizeFactor[pad];
        gPad->SetRightMargin(0.02);
        gPad->SetLeftMargin(0.15);
        if(pad){
            gPad->SetBottomMargin(botMarg);
        }
        else{
            gPad->SetTopMargin(topMarg);
            
        }
        mg1[pad]->Draw("ap");
        gPad->SetFillColor(0);
        gPad->SetFillStyle(0);
        gPad->GetFrame()->SetFillColor(0);
        gPad->GetFrame()->SetFillStyle(0);
        tl->SetTextSize(0.095*padSizeFactor);
        //     tl->DrawTextNDC(0.2,0.9-(pad?0:topMarg),fig1Tits[pad].Data());
        tl->DrawLatexNDC(0.2,0.9-(pad?0:topMarg),fig1Tits[pad].Data());
        gPad->SetTicks(1,1);
        gPad->Update();
        xaxis=mg1[pad]->GetXaxis();
        yaxis=mg1[pad]->GetYaxis();
        xaxis->SetLimits(-240,151);
        xaxis->SetTitleSize(titleSizeFig1);
        xaxis->SetTitleOffset(titleOffSFig1*1.75);
        xaxis->SetLabelSize(labelSizeFig1);
        xaxis->SetLabelOffset(labelOffSFig1);
        xaxis->SetNdivisions(6,6,1);
        yaxis->SetTitleSize(titleSizeFig1*padSizeFactor);
        yaxis->SetTitleOffset(titleOffSFig1);
        yaxis->SetLabelSize(labelSizeFig1*padSizeFactor);
        yaxis->SetLabelOffset(labelOffSFig1*padSizeFactor);
        yaxis->SetNdivisions(6,6,1);
    }
    
    mgu->GetYaxis()->SetRangeUser(-1.09,-0.41);
    mgd->GetYaxis()->SetRangeUser(-1.44,-0.76);
    
    
#ifdef MOD4
    const int nLeg1 = nMod/2 + nSet,     nLeg2 = 3 + nJLab - 1, nLeg21 = nSet, nLeg22 = 2 + nJLab - 1;
    int posLeg1 [nLeg1]  = { 0+nMod, 1+nMod,      3+nMod,      2+nMod, 2, 0 };
    int posLeg2 [nLeg2]  = { 0+nMod, 0+nMod+nSet, 1+nMod+nSet,         2, 0 };
    int posLeg21[nLeg21] = { 0+nMod, 1+nMod,      3+nMod,      2+nMod       };
    int posLeg22[nLeg22] = {         0+nMod+nSet, 1+nMod+nSet,         2, 0 };
#else
    const int nLeg1 = nMod/2 - 1 + nSet, nLeg2 = 2 + nJLab - 1, nLeg21 = nSet, nLeg22 = 1 + nJLab - 1;
    int posLeg1 [nLeg1]  = { 0+nMod, 1+nMod,      3+nMod,      2+nMod, 0 };
    int posLeg2 [nLeg2]  = { 0+nMod, 0+nMod+nSet, 1+nMod+nSet,         0 };
    int posLeg21[nLeg21] = { 0+nMod, 1+nMod,      3+nMod,      2+nMod    };
    int posLeg22[nLeg22] = {         0+nMod+nSet, 1+nMod+nSet,         0 };
#endif
    
    //   leg1 =  new TLegend(0.675, 0.05, 0.975, 0.95-topMarg, 0);
    leg1 =  new TLegend(0.65, 0.275, 0.95, 0.975, 0);
#ifdef MOD4
    leg2   = new TLegend(0.110, 0.570, 0.475, 0.895, 0);
    leg21  = new TLegend(0.110, 0.850, 0.475, 0.895, 0);
    leg22  = new TLegend(0.110, 0.570, 0.475, 0.850, 0);
    box2   = new TBox   (0.110, 0.570, 0.475, 0.895   );
    leg1  -> SetTextSize(0.08);
    leg2  -> SetTextSize(0.04);
    leg21 -> SetTextSize(0.04);
    leg22 -> SetTextSize(0.04);
#else
    leg2   = new TLegend(0.107, 0.575, 0.539, 0.887, 0);
    if (C12) {
        legF   = new TLegend(0.547, 0.575, 0.895, 0.895, 0);
        leg21  = new TLegend(0.547, 0.815, 0.895, 0.895, 0);
        leg22  = new TLegend(0.547, 0.575, 0.895, 0.815, 0);
    }
    else {
        legF   = new TLegend(0.135, 0.575, 0.455, 0.895, 0);
        leg21  = new TLegend(0.135, 0.815, 0.455, 0.895, 0);
        leg22  = new TLegend(0.135, 0.575, 0.455, 0.815, 0);
    }
    box2   = new TFrame (0.110, 0.575, 0.539, 0.875   );
    
    leg1  -> SetTextSize(0.10);
    leg2  -> SetTextSize(0.05);
    leg21 -> SetTextSize(0.05);
    leg22 -> SetTextSize(0.05);
#endif
    box2  -> SetFillStyle(0);
    box2  -> SetLineWidth(3);
    box2  -> SetLineColor(kBlack);
    legF  -> SetLineColor(kBlack);
    legF  -> SetFillStyle(0);
    legF  -> SetLineWidth(1);
    leg2  -> SetLineColor(0);
    leg21 -> SetLineColor(0);
    leg22 -> SetLineColor(0);
    leg1  -> SetFillColor(0);
    leg2  -> SetFillColor(0);
    legF  -> SetFillColor(kGreen);
    leg21 -> SetFillColor(0);
    leg22 -> SetFillColor(0);
    //   leg1  -> SetFillStyle(0);
    leg21 -> SetLineWidth(0);
    leg22 -> SetLineWidth(0);
    leg2  -> SetMargin(0.1);
    leg21 -> SetMargin(0.4);
    leg22 -> SetMargin(0.1);
    leg21 -> SetNColumns(4);
    leg21 -> SetColumnSeparation(0);
    for(int k=0;k<nLeg1; k++){
        leg1 ->  AddEntry(gr[posLeg1 [k]][0],tits1[posLeg1 [k]].Data(),"P");
    }
    for(int k=0;k<nLeg2; k++){
        leg2 ->  AddEntry(gr[posLeg2 [k]][1],tits2[posLeg2 [k]].Data(),"P");
    }
    for(int k=0;k<nLeg21;k++){
        leg21 -> AddEntry(gr[posLeg21[k]][1],tits2[posLeg21[k]].Data(),"P");
    }
    for(int k=0;k<nLeg22;k++){
        leg22 -> AddEntry(gr[posLeg22[k]][1],tits2[posLeg22[k]].Data(),"P");
    }
    c1->cd(/*1*/2);
    leg1 -> Draw("same");
    
    TLine *line;
    TGaxis *ax[2], *ay[2];
    double yUmin[2]={.544691858010134977, .887377999350428603};
    
    for ( int fg = 1; fg < /*nFigs*/3; fg++ ){
        TString fgname = Form ( "c%i", fg+1 );
        cc [fg] = new TCanvas ( fgname.Data(), fgname.Data() );
        mgA[fg+1] -> Draw ( "ap" );
        if ( fg == 3 )
        mgA[3] -> Draw ( "p" );
        gPad->SetTicks(1,1);
        xaxis = mgA[fg+1] -> GetXaxis();
        yaxis = mgA[fg+1] -> GetYaxis();
        TAxis * axs[nAxs] = { xaxis , yaxis };
        
        xaxis -> SetLimits(-.12,.07);
        xaxis -> SetBit(TAxis::kLabelsHori);
        const int nXLab=10;
        const int xLabBin[nXLab] = {1,11,21,32,42,53,63,74,84,95};
        double xLabs[nXLab];
        
        for (int lab=1; lab<nXLab; lab++){
            double v=-.12 + .02 * lab;
            xaxis -> SetBinLabel(xLabBin[lab], Form("%.2f", v>0 ? -v : v));
        }
        ax[fg-1]= new TGaxis (-.12,yUmin[fg-1],.07,yUmin[fg-1],-.12,.07,/*510*/xaxis->GetNdivisions(),"B");
        //ax[fg-1] -> SetNDivisions(6,6,1);
        ax[fg-1] -> SetLabelSize(0.);
        ax[fg-1] -> Draw();
        
        for ( int j = 0; j < nAxs; j++ ){
            axs[j] -> SetTitleSize(.045);
            axs[j] -> SetLabelSize(j?.045:.0705);
            if ( fg == 2 )
            axs[j] -> SetNdivisions(6,6,1);
        }
        xaxis -> SetBit(TAxis::kLabelsHori);
        
        ay[fg-1]=new TGaxis(0, yaxis->GetXmin(), 0, yaxis->GetXmax(), yaxis->GetXmin(),yaxis->GetXmax(), /*10606*/yaxis->GetNdivisions(), "+-");
        ay[fg-1] -> SetLabelSize(0.);
        ay[fg-1]->Draw();
        
        tl->SetTextSize(0.08);
        tl->SetTextAlign(31);
        tl->DrawLatexNDC(0.565,0.175,"p_{miss}<0");
        tl->DrawLatexNDC(0.865,0.175,"p_{miss}>0");
        
        if ( fg == 2 ) {
            c3 = cc[fg];
        }
        
        if ( fg == 3 ) {
            c4 = cc[fg];
        }
        
        if ( fg == 4 ) {
            c5 = cc[fg];
        }
        
        if ( fg == 5 ) {
            c6 = cc[fg];
        }
        
        if ( fg == 1 ) {
            c2 = cc[fg];
            
            
            
            
            THStack*curves=(THStack*)((new TFile("~/deep/output/hs.root"))->Get("hs")->Clone());
            double minx = -.103, maxx = -.0045;
            double miny =  .615, maxy = 1.3500;
#ifdef CURVES
            curves->Draw("SAME HIST NOSTACK C");
            xaxis->SetRangeUser(minx,maxx);
            yaxis->SetRangeUser(miny,maxy);
            curves->GetXaxis()->SetRangeUser(minx,maxx);
            curves->GetYaxis()->SetRangeUser(miny,maxy);
            line = new TLine(minx,1,maxx,1);
#else
            TH1 * h6 = (TH1*) curves -> GetHists() -> FindObject("_pfx_1_1_5");
#ifndef MOD4
            h6->SetLineColor(kBlue);
#endif
            h6->GetXaxis()->SetRangeUser(minx,maxx);
            h6->GetYaxis()->SetRangeUser(miny,maxy);
            TString curveLeg = "PWIA";
#ifdef VAN_ORDEN_COLS
            curveLeg += " HA";
#endif
            //       leg2  -> AddEntry(h6,curveLeg.Data(),"L");
            //       leg22 -> AddEntry(h6,curveLeg.Data(),"L");
            //       h6->Draw("SAME HIST NOSTACK C");
            line = new TLine(xaxis->GetXmin(),1,xaxis->GetXmax(),1);
#endif
            
            leg21 -> Draw("same");
            leg22 -> Draw("same");
            legF  -> Draw("same");
            //       ax[0]->Draw();
            //       
            //       box2  -> Draw("L");
            //       box2  -> Paint();
            
        }
        
        //====================================================
        // Draw a dashed line
        //====================================================
        line->SetLineStyle(2);
        line->Draw("same");
    }
    
    //====================================================
    // Save the figures
    //====================================================
    for ( int fg = 0; fg < /*nFigs*/3; fg++ ) {
        cc[fg] -> Update();
        cc[fg] -> SaveAs(Form("output/fig%i.pdf",fg+1));
        cc[fg] -> SaveAs(Form("output/fig%i.png",fg+1));
    }
}




