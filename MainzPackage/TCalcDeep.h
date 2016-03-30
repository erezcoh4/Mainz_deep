/**
 * \file TCalcDeep.h
 *
 * \ingroup MainzPackage
 * 
 * \brief Class def header for a class TCalcDeep
 *
 * @author erezcohen
 */

/** \addtogroup MainzPackage

    @{*/
#ifndef TCALCDEEP_H
#define TCALCDEEP_H

#include <iostream>
#include "MySoftwarePackage/TPlots.h"
#include "MySoftwarePackage/TAnalysis.h"
#include "Tqspin.h"
#define r2d TMath::RadToDeg()

/**
   \class TCalcDeep
   User defined class TCalcDeep ... these comments are used to generate
   doxygen documentation!
 */
class TCalcDeep{

public:
    
    TString DataType    , Path;
    TString InFileName  , InTreeName;
    TString OutFileName , OutTreeName;
    
    
    TFile   * InFile    , * OutFile;
    TTree   * InTree    , * OutTree;
    
    int     Nentries    , Entry;
    TString FrameName   ;          // prefered frame of axes to work in....
    
    
    Float_t Numer_wt    , Theta_pq  , Theta_cm  , Thet_prq  , Th_e_i    , Ph_e_i    , Pf_e_i;
    Float_t Th_p_i      , Ph_p_i    , Pf_p_i    , Qmag      , Theq      , Q2      , Xb     , Emiss , Prmag , W;
    Float_t Prec_x      , Prec_y    , Prec_z    ;
    Float_t PxRot       , PyRot     , PzRot    ;
    Float_t Px          , Py        , Pz;
    
    Float_t p_ref       , dp;

    
    
    TVector3 Prec       , S_tg      , S_HDC;
    
    TLorentzVector      e , p , prec;
    
    Tqspin qspin;
    
  /// Default constructor
    TCalcDeep (){InitGlobals();}
    TCalcDeep (TTree * fInTree, TTree * fOutTree);
  /// Default destructor
  ~TCalcDeep(){}

    
    
    
    // GETs
    
    
    // SETs
    void          SetPath (TString path){Path = path;};
    void        SetInTree (TTree * tree){InTree = tree;};
    void       SetOutTree (TTree * tree){OutTree = tree;};
    
    // initializations
    void      InitGlobals ();
    void    InitInputTree ();
    void   InitOutputTree ();
    
    
    void  ComputePhysVars (int entry);
    
    
    void       PrintData (int);

};

#endif
/** @} */ // end of doxygen group 

