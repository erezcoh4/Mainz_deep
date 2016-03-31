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
    dp = 10
    th = 0.1
    ph = 0.1
    y0 = 0
    p_ref = 523
    # E = 854.2 MeV electron
    # electron at -54.4 deg., momentum 619 MeV
    S_tg = ROOT.TVector3( 0.5 , 0.0 , 0.5 )
    S_HDC = qspin.pSpinPrecessionSpecA( dp , th , ph , y0 , p_ref, S_tg.X() , S_tg.Y() , S_tg.Z() )
    # input: <dp/%c> <th/mrad> <y0/mm> <ph/mrad> <p_ref/MeV/c> <Sx_tg> <Sy_tg> <Sz_tg>
    print "calculated single event precession with p_ref=%.1f MeV/c, dp = %.1f%%, ph = %.2f mrad, th = %.2f mrad, y0 = %.1fmm"%(p_ref,dp,ph,th,y0)
    S_tg.Print()
    S_HDC.Print()






