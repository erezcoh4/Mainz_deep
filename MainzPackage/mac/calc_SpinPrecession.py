import ROOT
from ROOT import TCalcDeep , Tqspin

DoMCEEPtreeCalculation  = False
DoSingeEventCalculation = True


if DoMCEEPtreeCalculation:
    Path        = "/Users/erezcohen/Desktop/MAINZ"
    FileName    = "Pol181_A_12_0"
    InFile      = ROOT.TFile(Path + "/simulation_data/" + FileName+".root")
    InTree      = InFile.Get("ntuple")
    Nentries    = InTree.GetEntries()
    OutFile     = ROOT.TFile(Path + "/AnaFiles/"+"Ana_"+FileName+".root","recreate")
    OutTree     = ROOT.TTree("anaTree","deep after spin precession")
    calc        = TCalcDeep(InTree , OutTree)
    for entry in range(0, 1000): # (int)(0.01*Nentries)
        calc.ComputePhysVars( entry );
        if entry%10==0:
            calc.PrintData( entry );
    print "done filling %d events" % OutTree.GetEntries()
    OutTree.Write()
    OutFile.Close()





if DoSingeEventCalculation:
    #   E0[GeV]	Q2[GeV2]    E'[GeV]	ThetaE0[Deg]	Pp[GeV/c]	ThetaP[Deg]	Theta[Rad]	Sin^2(Theta)    Mott[uBarn]	Xsec[uBarn]	Ge          Gm      I0      Px      Pz      Px/Pz	Ge/Gm
    
    #   0.855	0.4         0.641	50.5            0.667       47.91       0.882       0.182           0.131       0.349       0.409       1.142   0.389	-0.403	0.424	-0.950	0.358

    #   0.855	0.45        0.615	55.1            0.712       45.1        0.961       0.214           0.087       0.324       0.374       1.046   0.366	-0.424	0.484	-0.875	0.358

    #   0.855	0.5         0.588	59.8            0.755       42.3        1.043       0.248           0.059       0.305       0.344       0.962   0.349	-0.439	0.542	-0.809	0.358

    #   0.855	0.55        0.562	64.6            0.797       39.5        1.129       0.286           0.040       0.292       0.317       0.886   0.337	-0.448	0.599	-0.748	0.358
    qspin = Tqspin()
    Q2 = 0.5
    dp = 0
    th = 1043
    ph = 0
    y0 = 0
    p_ref = 755 
    S_tg_dipole = ROOT.TVector3( -0.439 , 0.0 , 0.542 )
    S_tg = S_tg_dipole
    S_HDC = qspin.pSpinPrecessionSpecA( dp , th , ph , y0 , p_ref, S_tg )
    # input: <dp/%c> <th_tg/mrad> <y0_tg/mm> <ph_tg/mrad> <p_ref/MeV/c> <Sx_tg> <Sy_tg> <Sz_tg>
    print "calculated single event precession with p_ref=%.1f MeV/c, dp = %.1f%%, ph = %.2f mrad, th = %.2f mrad, y0 = %.1fmm"%(p_ref,dp,ph,th,y0)
    S_tg.Print()
    S_HDC.Print()






