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
    th  = 0 # target coordinate
    ph  = 0 # target coordinate
    y0  = 0
    p   = 667
    p_ref = 630
    dp = 100*((p - p_ref)/p_ref)
    S_tg_dipole = ROOT.TVector3( -40.3 , 0 , 42.4 )
    S_tg = S_tg_dipole
    S_HDC = qspin.pSpinPrecessionSpecA( dp , th , ph , y0 , p_ref, S_tg )
    # input: <dp/%c> <th_tg/mrad> <y0_tg/mm> <ph_tg/mrad> <p_ref/MeV/c> <S_tg/%>
    print "p_ref=%.0f MeV/c, dp = %.1f%%, ph = %.1f mrad, th = %.1f mrad, y0 = %.1f mm"%(p_ref,dp,ph,th,y0)
    print "polarization at target (%.3f,%.3f,%.3f)"%(S_tg.X()/100,S_tg.Y()/100,S_tg.Z()/100)
    print "polarization at HDC (%.3f,%.3f,%.3f)"%(S_HDC.X()/100,S_HDC.Y()/100,S_HDC.Z()/100)






