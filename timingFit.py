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


barID = { 1:'027', 2: '045', 3: '059', 4: '095' , 5: '081', 6: '110', 7: '166', 8:'127', 9:'150'}

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
histosOut['ctr_ByProd_FNAL']=R.TGraphErrors()
histosOut['ctr_ByProd_FNAL'].SetName('ctr_ByProd_FNAL')
histosOut['ctrRel_ByProd_FNAL']=R.TGraphErrors()
histosOut['ctrRel_ByProd_FNAL'].SetName('ctrRel_ByProd_FNAL')
    
for prod in range(1,10):
    f=R.TFile("timeAnalysis_prod"+str(prod)+'_bar'+barID[prod]+".root")
    histos=Map(f)

    for v in ov:
        histosOut["prod"+str(prod)+'_bar'+barID[prod]+'_tDiff_ov'+v]=R.TGraphErrors()
        histosOut["prod"+str(prod)+'_bar'+barID[prod]+'_tDiff_ov'+v].SetName("prod"+str(prod)+'_bar'+barID[prod]+'_tDiff_ov'+v)
        for it,t in enumerate(th):
            key="ov{}_th{}".format(v,str(t))
            ff=fitGausSum(histos["tDiff_"+key])
            histosOut["prod"+str(prod)+'_bar'+barID[prod]+'_tDiff_ov'+v].SetPoint(it,t,ff.GetParameter(2)/2.)
            histosOut["prod"+str(prod)+'_bar'+barID[prod]+'_tDiff_ov'+v].SetPointError(it,0,ff.GetParError(2)/2.)
            c1.SaveAs("prod"+str(prod)+'_bar'+barID[prod]+"_tDiff_"+key+".png")

    a.Draw()
    leg=R.TLegend(0.15,0.65,0.4,0.85)
    leg.SetBorderSize(0)
    leg.SetFillColorAlpha(0,0)
    leg.SetTextSize(0.05)
    for iv,v in enumerate(ov):
        histosOut["prod"+str(prod)+'_bar'+barID[prod]+'_tDiff_ov'+v].SetMarkerColor(R.kBlack+iv)
        histosOut["prod"+str(prod)+'_bar'+barID[prod]+'_tDiff_ov'+v].SetLineColor(R.kBlack+iv)
        histosOut["prod"+str(prod)+'_bar'+barID[prod]+'_tDiff_ov'+v].SetMarkerStyle(20)
        histosOut["prod"+str(prod)+'_bar'+barID[prod]+'_tDiff_ov'+v].SetMarkerSize(1.2)
        histosOut["prod"+str(prod)+'_bar'+barID[prod]+'_tDiff_ov'+v].Draw("PLSAME")
        leg.AddEntry(histosOut["prod"+str(prod)+'_bar'+barID[prod]+'_tDiff_ov'+v],"OV=%s V"%v,"LP")
    leg.Draw()
    l.DrawLatexNDC(0.12,0.93,"prod"+str(prod)+'_bar'+barID[prod]+' - S12572-015C - protons 120 GeV')
    c1.SaveAs("prod"+str(prod)+'_bar'+barID[prod]+"_tDiff.png")

    xF,yF=R.Double(0),R.Double(0)
    thVals=[1,2]
    sigmaT,sigmaT_err=0.,0.

    for t in thVals:
        histosOut["prod"+str(prod)+'_bar'+barID[prod]+'_tDiff_ov6.0'].GetPoint(t,xF,yF)
        yE=histosOut["prod"+str(prod)+'_bar'+barID[prod]+'_tDiff_ov6.0'].GetErrorY(t)
        sigmaT=sigmaT+yF
        sigmaT_err=sigmaT_err+yE
    sigmaT=sigmaT/len(thVals)
    sigmaT_err=sigmaT_err/len(thVals)

    if prod==1:
        ref=sigmaT

    histosOut['ctr_ByProd_FNAL'].SetPoint(prod-1,prod,sigmaT)
    histosOut['ctr_ByProd_FNAL'].SetPointError(prod-1,0,sigmaT_err)
    histosOut['ctrRel_ByProd_FNAL'].SetPoint(prod-1,prod,sigmaT/ref)
    histosOut['ctrRel_ByProd_FNAL'].SetPointError(prod-1,0,sigmaT_err/ref)

l.SetTextSize(0.05)

c1.SetGridy(1)
a1=R.TH2F("a","a",10,0,10,10,40,70)
a1.GetXaxis().SetTitle("Vendor ID")
a1.GetYaxis().SetTitle("#sigma_{bar} [ps]")
a1.SetStats(0)
a1.Draw()
histosOut['ctr_ByProd_FNAL'].SetMarkerStyle(20)
histosOut['ctr_ByProd_FNAL'].SetMarkerSize(1.2)
histosOut['ctr_ByProd_FNAL'].Draw("PSAME")

l.SetTextSize(0.05)
l.DrawLatexNDC(0.12,0.92,"protons 120 GeV")

l.SetTextSize(0.04)
l.DrawLatexNDC(0.13,0.85,"LYSO:Ce 3x3x57 mm^{3}")
l.DrawLatexNDC(0.6,0.85,"S125712-015C OV=6V")

c1.SaveAs("ctr_ByProd_FNAL.png")

a2=R.TH2F("a2","a2",10,0,10,10,0.8,1.3)
a2.GetXaxis().SetTitle("Vendor ID")
a2.GetYaxis().SetTitle("#sigma_{bar} [ps]")
a2.SetStats(0)
a2.Draw()
histosOut['ctrRel_ByProd_FNAL'].SetMarkerStyle(20)
histosOut['ctrRel_ByProd_FNAL'].SetMarkerSize(1.2)
histosOut['ctrRel_ByProd_FNAL'].Draw("PSAME")

l.SetTextSize(0.05)
l.DrawLatexNDC(0.12,0.92,"protons 120 GeV")

l.SetTextSize(0.04)
l.DrawLatexNDC(0.13,0.85,"LYSO:Ce 3x3x57 mm^{3}")
l.DrawLatexNDC(0.6,0.85,"S125712-015C OV=6V")

c1.SaveAs("ctrRel_ByProd_FNAL.png")

c1.SetGridy(0)
f1=R.TFile("romeCernCorrelation.root")
histosRomeCern=Map(f1)

histosOut['ctr_RomeVsFNAL']=R.TGraphErrors()
histosOut['ctr_RomeVsFNAL'].SetName('ctr_RomeVsFNAL')
histosOut['ctr_CERNVsFNAL']=R.TGraphErrors()
histosOut['ctr_CERNVsFNAL'].SetName('ctr_CERNVsFNAL')

histosOut['ctrRel_RomeVsFNAL']=R.TGraphErrors()
histosOut['ctrRel_RomeVsFNAL'].SetName('ctrRel_RomeVsFNAL')
histosOut['ctrRel_CERNVsFNAL']=R.TGraphErrors()
histosOut['ctrRel_CERNVsFNAL'].SetName('ctrRel_CERNVsFNAL')

for prod in range(1,10):
    xF,yF,xR,yR,xC,yC=R.Double(0),R.Double(0),R.Double(0),R.Double(0),R.Double(0),R.Double(0)
    histosOut['ctr_ByProd_FNAL'].GetPoint(prod-1,xF,yF)
    histosRomeCern['ctr_ByProd_ROME'].GetPoint(prod-1,xR,yR)
    histosRomeCern['ctr_ByProd_CERN'].GetPoint(prod-1,xC,yC)

    if prod==1:
        relF=yF
        relC=yC
        relR=yR

    histosOut['ctr_RomeVsFNAL'].SetPoint(prod-1,yR,yF)
    histosOut['ctr_RomeVsFNAL'].SetPointError(prod-1,histosRomeCern['ctr_ByProd_ROME'].GetErrorY(prod-1),histosOut['ctr_ByProd_FNAL'].GetErrorY(prod-1))
    histosOut['ctr_CERNVsFNAL'].SetPoint(prod-1,yC,yF)
    histosOut['ctr_CERNVsFNAL'].SetPointError(prod-1,histosRomeCern['ctr_ByProd_CERN'].GetErrorY(prod-1),histosOut['ctr_ByProd_FNAL'].GetErrorY(prod-1))

    histosOut['ctrRel_RomeVsFNAL'].SetPoint(prod-1,yR/relR,yF/relF)
    histosOut['ctrRel_RomeVsFNAL'].SetPointError(prod-1,histosRomeCern['ctr_ByProd_ROME'].GetErrorY(prod-1)/relR,histosOut['ctr_ByProd_FNAL'].GetErrorY(prod-1)/relF)
    histosOut['ctrRel_CERNVsFNAL'].SetPoint(prod-1,yC/relC,yF/relF)
    histosOut['ctrRel_CERNVsFNAL'].SetPointError(prod-1,histosRomeCern['ctr_ByProd_CERN'].GetErrorY(prod-1)/relC,histosOut['ctr_ByProd_FNAL'].GetErrorY(prod-1)/relF)

histosOut['ctr_RomeVsFNAL'].SetMarkerStyle(20)
histosOut['ctr_RomeVsFNAL'].SetMarkerSize(1.2)
histosOut['ctr_RomeVsFNAL'].GetXaxis().SetLimits(80,200)
histosOut['ctr_RomeVsFNAL'].GetXaxis().SetRangeUser(80,200)
histosOut['ctr_RomeVsFNAL'].GetXaxis().SetTitle("#sigma_{t}^{ROME}");
histosOut['ctr_RomeVsFNAL'].GetYaxis().SetLimits(50,70)
histosOut['ctr_RomeVsFNAL'].GetYaxis().SetRangeUser(50,70)
histosOut['ctr_RomeVsFNAL'].GetYaxis().SetTitle("#sigma_{t}^{FNAL}");
histosOut['ctr_RomeVsFNAL'].Draw("AP")
c1.SaveAs("ctr_RomeVsFNAL.png")

histosOut['ctr_CERNVsFNAL'].SetMarkerStyle(20)
histosOut['ctr_CERNVsFNAL'].SetMarkerSize(1.2)
histosOut['ctr_CERNVsFNAL'].GetXaxis().SetLimits(75,110)
histosOut['ctr_CERNVsFNAL'].GetXaxis().SetRangeUser(75,110)
histosOut['ctr_CERNVsFNAL'].GetXaxis().SetTitle("#sigma_{t}^{CERN}");
histosOut['ctr_CERNVsFNAL'].GetYaxis().SetLimits(50,70)
histosOut['ctr_CERNVsFNAL'].GetYaxis().SetRangeUser(50,70)
histosOut['ctr_CERNVsFNAL'].GetYaxis().SetTitle("#sigma_{t}^{FNAL}");
histosOut['ctr_CERNVsFNAL'].Draw("AP")
c1.SaveAs("ctr_CERNVsFNAL.png")

histosOut['ctrRel_CERNVsFNAL'].SetMarkerStyle(20)
histosOut['ctrRel_CERNVsFNAL'].SetMarkerSize(1.2)
histosOut['ctrRel_CERNVsFNAL'].GetXaxis().SetLimits(0.85,1.25)
histosOut['ctrRel_CERNVsFNAL'].GetXaxis().SetRangeUser(0.85,1.25)
histosOut['ctrRel_CERNVsFNAL'].GetXaxis().SetTitle("#sigma_{t}^{CERN}");
histosOut['ctrRel_CERNVsFNAL'].GetYaxis().SetLimits(0.85,1.3)
histosOut['ctrRel_CERNVsFNAL'].GetYaxis().SetRangeUser(0.85,1.25)
histosOut['ctrRel_CERNVsFNAL'].GetYaxis().SetTitle("#sigma_{t}^{FNAL}");
histosOut['ctrRel_CERNVsFNAL'].Draw("AP")
l=R.TF1("l","x",0.7,1.3)
l.Draw("SAME")
c1.SaveAs("ctrRel_CERNVsFNAL.png")

histosOut['ctrRel_RomeVsFNAL'].SetMarkerStyle(20)
histosOut['ctrRel_RomeVsFNAL'].SetMarkerSize(1.2)
histosOut['ctrRel_RomeVsFNAL'].GetXaxis().SetLimits(0.7,1.8)
histosOut['ctrRel_RomeVsFNAL'].GetXaxis().SetRangeUser(0.7,1.8)
histosOut['ctrRel_RomeVsFNAL'].GetXaxis().SetTitle("#sigma_{t}^{ROME}");
histosOut['ctrRel_RomeVsFNAL'].GetYaxis().SetLimits(0.85,1.2)
histosOut['ctrRel_RomeVsFNAL'].GetYaxis().SetRangeUser(0.85,1.2)
histosOut['ctrRel_RomeVsFNAL'].GetYaxis().SetTitle("#sigma_{t}^{FNAL}");
histosOut['ctrRel_RomeVsFNAL'].Draw("AP")
c1.SaveAs("ctrRel_RomeVsFNAL.png")

fOut=R.TFile("timingAnalysisFNAL.root","RECREATE")
for hN,h in histosOut.items():
    h.Write()
fOut.Write()
fOut.Close()
