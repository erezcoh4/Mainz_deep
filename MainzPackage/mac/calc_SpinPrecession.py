import ROOT
from ROOT import TCalcDeep

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






