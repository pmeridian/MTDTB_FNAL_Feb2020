import ROOT as R


def Map(tf):
    m = {}
    for k in tf.GetListOfKeys():
        n = k.GetName()
        m[n] = tf.Get(n)
    return m

def gausSum(x,par):
    return par[0]*R.TMath.Gaus(x[0],par[1],par[2])+par[3]*R.TMath.Gaus(x[0],par[1]+par[4],par[5]*par[2])+par[6]*R.TMath.Gaus(x[0],par[1]-par[4],par[5]*par[2])

def fitGausSum(h):
    f=R.TF1(h.GetName()+"_gausSumFit",gausSum,-1000,1000,7)
    hmax=h.GetXaxis().GetBinCenter(h.GetMaximumBin())
    hrms=h.GetRMS()
    f.SetParameter(0,100)
    f.SetParLimits(0,0,999999)
    f.SetParameter(1,hmax)
    f.SetParLimits(1,hmax-50,hmax+50)
    f.SetParameter(2,100)
    f.SetParLimits(2,50,180)
    f.SetParameter(3,50)
    f.SetParLimits(3,0,999999)
    f.SetParameter(4,350)
    f.SetParLimits(4,50.,500)
    f.SetParameter(5,2.)
    f.SetParLimits(5,1.,5.)
    f.SetParameter(6,50)
    f.SetParLimits(6,0,999999)

    f.Print()
    h.Fit(f,"RB","", -1000,1000)
    return f
    #    c1.SaveAs(h.GetName()+"_gausSumFit.png")


barID = { '3.10':'067c3.10', '8.0': '067c8.0', '8.2': '067c8.2', }

ov= [ "4.0", "6.0", "8.0" ]
th = [ 10, 20, 40 ]

R.gROOT.SetBatch(1)

c1=R.TCanvas("c1","c1",800,600)
R.gStyle.SetOptFit(111111)
R.gStyle.SetOptTitle(0)

l=R.TLatex()
l.SetTextSize(0.05)

histosOut={}
a=R.TH2F("a","a",10,0,50,10,0,100)
a.GetXaxis().SetTitle("Threshold [DAC counts]")
a.GetYaxis().SetTitle("#sigma_{bar} [ps]")
a.SetStats(0)

histosOut={}
histosOut['ctr_ByConf_FNAL']=R.TGraphErrors()
histosOut['ctr_ByConf_FNAL'].SetName('ctr_ByProd_FNAL')
    
for iprod,prod in enumerate(barID.keys()):
    f=R.TFile("timeAnalysis_prod1"+'_bar'+barID[prod]+".root")
    histos=Map(f)

    for v in ov:
        histosOut["prod1"+'_bar'+barID[prod]+'_tDiff_ov'+v]=R.TGraphErrors()
        histosOut["prod1"+'_bar'+barID[prod]+'_tDiff_ov'+v].SetName("prod"+'_bar'+barID[prod]+'_tDiff_ov'+v)
        for it,t in enumerate(th):
            key="ov{}_th{}".format(v,str(t))
            ff=fitGausSum(histos["tDiff_"+key])
            histosOut["prod1"+'_bar'+barID[prod]+'_tDiff_ov'+v].SetPoint(it,t,ff.GetParameter(2)/2.)
            histosOut["prod1"+'_bar'+barID[prod]+'_tDiff_ov'+v].SetPointError(it,0,ff.GetParError(2)/2.)
            c1.SaveAs("prod1"+'_bar'+barID[prod]+"_tDiff_"+key+".png")

    a.Draw()
    leg=R.TLegend(0.15,0.65,0.4,0.85)
    leg.SetBorderSize(0)
    leg.SetFillColorAlpha(0,0)
    leg.SetTextSize(0.05)
    for iv,v in enumerate(ov):
        histosOut["prod1"+'_bar'+barID[prod]+'_tDiff_ov'+v].SetMarkerColor(R.kBlack+iv)
        histosOut["prod1"+'_bar'+barID[prod]+'_tDiff_ov'+v].SetLineColor(R.kBlack+iv)
        histosOut["prod1"+'_bar'+barID[prod]+'_tDiff_ov'+v].SetMarkerStyle(20)
        histosOut["prod1"+'_bar'+barID[prod]+'_tDiff_ov'+v].SetMarkerSize(1.2)
        histosOut["prod1"+'_bar'+barID[prod]+'_tDiff_ov'+v].Draw("PLSAME")
        leg.AddEntry(histosOut["prod1"+'_bar'+barID[prod]+'_tDiff_ov'+v],"OV=%s V"%v,"LP")
    leg.Draw()
    l.DrawLatexNDC(0.12,0.93,"prod1"+'_bar'+barID[prod]+' - S12572-015C - protons 120 GeV')
    c1.SaveAs("prod1"+'_bar'+barID[prod]+"_tDiff.png")

    xF,yF=R.Double(0),R.Double(0)
    thVals=[1,2]
    sigmaT,sigmaT_err=0.,0.

    for t in thVals:
        histosOut["prod1"+'_bar'+barID[prod]+'_tDiff_ov6.0'].GetPoint(t,xF,yF)
        yE=histosOut["prod1"+'_bar'+barID[prod]+'_tDiff_ov6.0'].GetErrorY(t)
        sigmaT=sigmaT+yF
        sigmaT_err=sigmaT_err+yE
    sigmaT=sigmaT/len(thVals)
    sigmaT_err=sigmaT_err/len(thVals)

    histosOut['ctr_ByConf_FNAL'].SetPoint(iprod,iprod+1,sigmaT)
    histosOut['ctr_ByConf_FNAL'].SetPointError(iprod,0,sigmaT_err)


c1.SetGridy(1)

l=R.TLatex()
l.SetTextSize(0.05)

a1=R.TH2F("a","a",len(barID.keys()),0.5,0.5+len(barID.keys()),10,30,70)
a1.GetXaxis().SetTitle("Conf ID")
a1.GetYaxis().SetTitle("#sigma_{bar} [ps]")
a1.SetStats(0)
a1.Draw()
histosOut['ctr_ByConf_FNAL'].SetMarkerStyle(20)
histosOut['ctr_ByConf_FNAL'].SetMarkerSize(1.2)
histosOut['ctr_ByConf_FNAL'].Draw("PSAME")
for ic,conf in enumerate(barID. keys()):
    a1.GetXaxis().SetBinLabel(ic+1,conf)

a1.GetXaxis().SetLabelSize(0.05)
l.DrawLatexNDC(0.12,0.92,"bar# 67(prod1) - protons 120 GeV")

l.SetTextSize(0.04)
l.DrawLatexNDC(0.13,0.85,"LYSO:Ce 3x3x57 mm^{3}")
l.DrawLatexNDC(0.6,0.85,"S125712-015C OV=6V")

c1.SaveAs("ctr_ByConf_FNAL.png")

fOut=R.TFile("timingConfAnalysisFNAL.root","RECREATE")
for hN,h in histosOut.items():
    h.Write()
fOut.Write()
fOut.Close()
