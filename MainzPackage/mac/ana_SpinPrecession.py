import ROOT
from ROOT import TPlots
from rootpy.interactive import wait
import time
ROOT.gStyle.SetOptStat(0000)


DoSpinBeforeAfter   = True


Path        = "/Users/erezcohen/Desktop/MAINZ"
FileName    = "Pol181_A_12_0"
ana         = TPlots(Path + "/AnaFiles/"+"Ana_"+FileName+".root","anaTree","SpinPrecession")
Nbins       = 1000
Smin        = -1
Smax        = 1



if DoSpinBeforeAfter :
    c = ana.CreateCanvas("c.m. momentum","Divide",3,2)
    c.cd(1)
    ana.H1("S_tg.X()",ROOT.TCut(),"HIST",Nbins,Smin,Smax,"at target position","S(tg) ^{x}","",1,46)
    c.cd(2)
    ana.H1("S_tg.Y()",ROOT.TCut(),"HIST",Nbins,Smin,Smax,"at target position","S(tg) ^{y}","",1,46)
    c.cd(3)
    ana.H1("S_tg.Z()",ROOT.TCut(),"HIST",Nbins,Smin,Smax,"at target position","S(tg) ^{z}","",1,46)
    c.cd(4)
    ana.H1("S_HDC.X()",ROOT.TCut(),"HIST",Nbins,Smin,Smax,"at HDC","S(HDC) ^{x}","",1,46)
    c.cd(5)
    ana.H1("S_HDC.Y()",ROOT.TCut(),"HIST",Nbins,Smin,Smax,"at HDC","S(HDC) ^{y}","",1,46)
    c.cd(6)
    ana.H1("S_HDC.Z()",ROOT.TCut(),"HIST",Nbins,Smin,Smax,"at HDC","S(HDC) ^{z}","",1,46)
    c.Update()
    wait()

