import ROOT , sys , os , time
from ROOT import TPlots
from rootpy.interactive import wait
import time
ROOT.gStyle.SetOptStat(0000)
ROOT.gStyle.SetOptFit(1111)



AnalysisType = "F1/F2 from paper"


dirfmt = "~/Desktop/%4d-%02d-%02d"
dirname = dirfmt % time.localtime()[0:3]
try:
    os.makedirs(dirname)
except OSError, e:
    if e.errno != 17:
        raise # This was not a "directory exist" error..


if AnalysisType == "F1/F2 from paper":
    anaF1 = TPlots("/Users/erezcohen/Desktop/MAINZ/Software/MediumModifications/F1P_Solid.root","ntuple","F1")
    anaF2 = TPlots("/Users/erezcohen/Desktop/MAINZ/Software/MediumModifications/F2P_Solid.root","ntuple","F2")
    c = anaF1.CreateCanvas("F1 and F2","Divide",2,1)
    c.cd(1)
    hF1 = anaF1.H2("q2","f",ROOT.TCut(),"",100,0.1,0.9,100,0,2.0,"F1","Q^{2}  (GeV/c) ^{2}","F _{1}")
    hF1.Fit("pol9","","R",0.133,0.88)
    c.cd(2)
    hF2 = anaF2.H2("q2","f",ROOT.TCut(),"",100,0.1,0.9,100,0,2.0,"F2","Q^{2} (GeV/c) ^{2}","F _{2}")
    hF2.Fit("pol9","","R",0.2,0.88)
    c.Update()
    wait()
    c.SaveAs(dirname+"/F1F2.pdf")