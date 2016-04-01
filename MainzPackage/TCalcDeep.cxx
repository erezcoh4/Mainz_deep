#ifndef TCALCDEEP_CXX
#define TCALCDEEP_CXX

#include "TCalcDeep.h"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TCalcDeep::TCalcDeep( TTree * fInTree, TTree * fOutTree ){
    SetInTree       (fInTree);
    SetOutTree      (fOutTree);
    InitGlobals     ();
    InitInputTree   ();
    InitOutputTree  ();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TCalcDeep::InitGlobals(){
    p_ref = 630;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TCalcDeep::InitInputTree(){
    
//    InTree -> SetBranchAddress("Theta_pq"       , &Theta_pq);
    InTree -> SetBranchAddress("Theta_cm"       , &Theta_cm);
    InTree -> SetBranchAddress("Thet_prq"       , &Thet_prq);
    InTree -> SetBranchAddress("theta"          , &Th_e_i);
    InTree -> SetBranchAddress("phi_e"          , &Ph_e_i);
    InTree -> SetBranchAddress("p_e"            , &Pf_e_i);
    InTree -> SetBranchAddress("theta_p"        , &Th_p_i);
    InTree -> SetBranchAddress("phi_p"          , &Ph_p_i);
    InTree -> SetBranchAddress("Px0"            , &Px);
    InTree -> SetBranchAddress("Py0"            , &Py);
    InTree -> SetBranchAddress("Pz0"            , &Pz);
    InTree -> SetBranchAddress("Qmag"           , &Qmag);
    InTree -> SetBranchAddress("Theq"           , &Theq);
//    InTree -> SetBranchAddress("Qmu2"           , &Q2);
    InTree -> SetBranchAddress("X"              , &Xb);
    InTree -> SetBranchAddress("Emiss"          , &Emiss);
    InTree -> SetBranchAddress("Prmag"          , &Prmag);
    InTree -> SetBranchAddress("W"              , &W);
    InTree -> SetBranchAddress("Prec_x"         , &Prec_x);
    InTree -> SetBranchAddress("PxRot"          , &PxRot);
    InTree -> SetBranchAddress("PyRot"          , &PyRot);
    InTree -> SetBranchAddress("PzRot"          , &PzRot);

    
    
    Nentries    = InTree -> GetEntries();
    std::cout << "Initialized Input InTree TCalcDeep for " << InTree -> GetName() <<", Nentries = " <<  Nentries << std::endl;
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TCalcDeep::InitOutputTree(){
    
    
    // Float_t branches
    OutTree -> Branch("Xb"                  ,&Xb                    , "Xb/F");
    OutTree -> Branch("Q2"                  ,&Q2                    , "Q2/F");
    
    
    // TVector3 branches
    OutTree -> Branch("e"                   ,"TLorentzVector"       ,&e);
    OutTree -> Branch("p"                   ,"TLorentzVector"       ,&p);
    OutTree -> Branch("prec"                ,"TLorentzVector"       ,&prec);
    
    OutTree -> Branch("S_tg"                ,"TVector3"             ,&S_tg);    // spin in target position
    OutTree -> Branch("S_HDC"               ,"TVector3"             ,&S_HDC);    // spin in HDC

    
    std::cout << "Initialized Output Tree TCalcPhysVarsEG2 on " << OutTree -> GetTitle() << std::endl;
    
    
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TCalcDeep::ComputePhysVars(int entry){
    
    InTree -> GetEntry(entry);

    Q2 = (Qmag/1000)*(Qmag/1000);
    e.SetXYZM( Pf_e_i*sin(Th_e_i)*cos(Ph_e_i) , Pf_e_i*sin(Th_e_i)*sin(Ph_e_i)  , Pf_e_i*cos(Th_e_i) , 0.511 );
    p.SetXYZM( 1000*Px , 1000*Py , 1000*Pz  , 938 );
    prec.SetXYZM( Prec_x , Prec_y , Prec_z, 938 );
    
    S_tg = TVector3(PxRot , PyRot , PzRot);
    
    dp = 100.*(p.P() / p_ref - 1.);
    
    // <dp/%c> <th/mrad> <y0/mm> <ph/mrad> <p_ref/MeV/c> <Sx_tg> <Sy_tg> <Sz_tg>
    S_HDC = qspin.pSpinPrecessionSpecA( dp , Th_p_i  , 0 , Ph_p_i , p_ref , S_tg );
    
    OutTree -> Fill();

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TCalcDeep::PrintData(int entry){

    SHOW(entry);
    SHOW(Xb);
    SHOW(Q2);
    SHOWTLorentzVector(e);
    SHOWTLorentzVector(p);
    SHOWTLorentzVector(prec);
    SHOWTVector3(S_tg);
    SHOWTVector3(S_HDC);
    PrintLine();

}

#endif







