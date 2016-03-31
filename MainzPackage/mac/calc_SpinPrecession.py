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
    qspin = Tqspin()
    Q2 = 0.445
    Xv = 1.006
    dp = 0
    th = 45.5 * 3.14 / 180;
    ph = 0
    y0 = 0
    p_ref = 523
    # E = 854.2 MeV electron
    # electron at -54.4 deg., momentum 619 MeV
    # proton at 51.4 deg., momentum 523 MeV
    # proton at 45.5 deg., momentum 705 MeV
    S_tg = ROOT.TVector3( 0.5 , 0.0 , 0.5 )
    S_HDC = qspin.pSpinPrecessionSpecA( dp , th , ph , y0 , p_ref, S_tg.X() , S_tg.Y() , S_tg.Z() )
    # input: <dp/%c> <th_tg/mrad> <y0_tg/mm> <ph_tg/mrad> <p_ref/MeV/c> <Sx_tg> <Sy_tg> <Sz_tg>
    print "calculated single event precession with p_ref=%.1f MeV/c, dp = %.1f%%, ph = %.2f mrad, th = %.2f mrad, y0 = %.1fmm"%(p_ref,dp,ph,th,y0)
    S_tg.Print()
    S_HDC.Print()






