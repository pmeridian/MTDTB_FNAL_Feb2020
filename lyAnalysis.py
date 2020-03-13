import ROOT as R

R.gROOT.SetBatch(1)

def Map(tf):
    """                                                                                                                  
    Maps objets as dict[obj_name][0] using a TFile (tf) and TObject to browse.                                           
    """
    m = {}
    for k in tf.GetListOfKeys():
        n = k.GetName()
        m[n] = tf.Get(n)
    return m

f=R.TFile("histos.root")
histos=Map(f)

barID = { 1:'027', 2: '045', 3: '059', 4: '095' , 5: '081', 6: '110', 7: '166', 8:'127', 9:'150'}

refLY=-9999

histosOut={}
histosOut['lyS_ByProd_FNAL']=R.TGraphErrors()
histosOut['lyS_ByProd_FNAL'].SetName('lyByProd_FNAL')
for prod in range(1,10):
    if (prod==1):
        ref=histos['prod'+str(prod)+'_bar'+barID[prod]+'_energyTot_ov6.0_lanGausFit'].GetParameter(1)/1.05
    lyRel=histos['prod'+str(prod)+'_bar'+barID[prod]+'_energyTot_ov6.0_lanGausFit'].GetParameter(1)/ref
    histosOut['lyS_ByProd_FNAL'].SetPoint(prod-1,prod,lyRel)

f1=R.TFile("romeCernCorrelation.root")
histosRomeCern=Map(f1)

histosOut['lyS_RomeVsFNAL']=R.TGraphErrors()
histosOut['lyS_RomeVsFNAL'].SetName('lyS_RomeVsFNAL')
histosOut['lyS_RomeFNALRatio']=R.TH1F('lyS_RomeFNALRatio','lyS_RomeFNALRatio',15,0.8,1.2)

histosOut['lyS_CERNVsFNAL']=R.TGraphErrors()
histosOut['lyS_CERNVsFNAL'].SetName('lyS_CERNVsFNAL')
histosOut['lyS_CERNFNALRatio']=R.TH1F('lyS_CERNFNALRatio','lyS_CERNFNALRatio',15,0.8,1.2)

for prod in range(1,10):
    xF,yF,xR,yR,xC,yC=R.Double(0),R.Double(0),R.Double(0),R.Double(0),R.Double(0),R.Double(0)
    histosOut['lyS_ByProd_FNAL'].GetPoint(prod-1,xF,yF)
    histosRomeCern['lyS_ByProd_ROME'].GetPoint(prod-1,xR,yR)
    histosRomeCern['lyS_ByProd_CERN'].GetPoint(prod-1,xC,yC)

    histosOut['lyS_RomeVsFNAL'].SetPoint(prod-1,yR,yF)
    histosOut['lyS_RomeVsFNAL'].SetPointError(prod-1,histosRomeCern['lyS_ByProd_ROME'].GetErrorY(prod-1),0.035)
    histosOut['lyS_RomeFNALRatio'].Fill(yR/yF)

    histosOut['lyS_CERNVsFNAL'].SetPoint(prod-1,yC,yF)
    histosOut['lyS_CERNVsFNAL'].SetPointError(prod-1,histosRomeCern['lyS_ByProd_CERN'].GetErrorY(prod-1),0.035)
    histosOut['lyS_CERNFNALRatio'].Fill(yC/yF)

c1=R.TCanvas("c1","c1",800,600)

R.gStyle.SetOptTitle(0)

a=R.TH2F("a","a",10,0.7,1.3,10,0.7,1.3)
a.GetXaxis().SetTitle('LY_{ROME}')
a.GetYaxis().SetTitle('LY_{FNAL}')
a.SetStats(0)
a.Draw()
histosOut['lyS_RomeVsFNAL'].SetMarkerStyle(20)
histosOut['lyS_RomeVsFNAL'].SetMarkerSize(0.8)
histosOut['lyS_RomeVsFNAL'].SetMarkerColor(R.kBlack)
histosOut['lyS_RomeVsFNAL'].Draw("PSAME")
l=R.TF1("l","x",0.7,1.3)
l.Draw("SAME")
c1.SaveAs("lyS_RomeVsFNAL.png")

a.GetXaxis().SetTitle('LY_{CERN}')
a.GetYaxis().SetTitle('LY_{FNAL}')
a.SetStats(0)
a.Draw()
histosOut['lyS_CERNVsFNAL'].SetMarkerStyle(20)
histosOut['lyS_CERNVsFNAL'].SetMarkerSize(0.8)
histosOut['lyS_CERNVsFNAL'].SetMarkerColor(R.kBlack)
histosOut['lyS_CERNVsFNAL'].Draw("PSAME")
l=R.TF1("l","x",0.7,1.3)
l.Draw("SAME")
c1.SaveAs("lyS_CERNVsFNAL.png")

R.gStyle.SetOptStat(0)
R.gStyle.SetOptFit(111)

histosOut['lyS_RomeFNALRatio'].GetXaxis().SetTitle('LY_{ROME}/LY_{FNAL}')
histosOut['lyS_RomeFNALRatio'].Fit("gaus","L")
histosOut['lyS_RomeFNALRatio'].SetMarkerStyle(20)
histosOut['lyS_RomeFNALRatio'].SetMarkerSize(1.2)
histosOut['lyS_RomeFNALRatio'].SetMarkerColor(R.kBlack)
histosOut['lyS_RomeFNALRatio'].SetLineColor(R.kBlack)
histosOut['lyS_RomeFNALRatio'].Draw("PE")
c1.SaveAs("lyS_RomeFNALRatio.png")

histosOut['lyS_CERNFNALRatio'].GetXaxis().SetTitle('LY_{CERN}/LY_{FNAL}')
histosOut['lyS_CERNFNALRatio'].Fit("gaus","L")
histosOut['lyS_CERNFNALRatio'].SetMarkerStyle(20)
histosOut['lyS_CERNFNALRatio'].SetMarkerSize(1.2)
histosOut['lyS_CERNFNALRatio'].SetMarkerColor(R.kBlack)
histosOut['lyS_CERNFNALRatio'].SetLineColor(R.kBlack)
histosOut['lyS_CERNFNALRatio'].Draw("PE")
c1.SaveAs("lyS_CERNFNALRatio.png")

fOut=R.TFile("lyAnalysisFNAL.root","RECREATE")
for hN,h in histosOut.items():
    h.Write()
fOut.Write()
fOut.Close()
