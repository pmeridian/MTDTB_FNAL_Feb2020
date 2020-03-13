#!/usr/bin/env python
# coding: utf-8

# In[2]:


import ROOT as R
#R.ROOT.EnableImplicitMT()

run_ranges = {}
run_ranges['prod2_bar045']=[ 20912,21019 ]
run_ranges['prod3_bar059']=[ 21033,21081 ]
run_ranges['prod4_bar095']=[ 21086,21143 ]
run_ranges['prod5_bar081']=[ 21150,21378 ]
run_ranges['prod6_bar110']=[ 22014,22040 ]
run_ranges['prod7_bar166']=[ 22240,22284 ]
run_ranges['prod8_bar127']=[ 22299,22344 ]
run_ranges['prod9_bar150']=[ 22361,22406 ]
run_ranges['prod1_bar027']=[ 22445,22508 ]
run_ranges['prod1_bar067c3.10']=[ 22512,22538 ]
run_ranges['prod1_bar067c8.0']=[ 23657,23876 ]
run_ranges['prod1_bar067c8.2']=[ 23980,24200 ]

channels = {}
channels['prod2_bar045']=[ 213,219 ]
channels['prod3_bar059']=[ 213,219 ]
channels['prod4_bar095']=[ 213,219 ]
channels['prod5_bar081']=[ 213,219 ]
channels['prod6_bar110']=[ 213,219 ]
channels['prod7_bar166']=[ 213,219 ]
channels['prod8_bar127']=[ 213,219 ]
channels['prod9_bar150']=[ 213,219 ]
channels['prod1_bar027']=[ 213,219 ]
channels['prod1_bar067c3.10']=[ 213,219 ]
channels['prod1_bar067c8.0']=[ 195,217 ]
channels['prod1_bar067c8.2']=[ 195,217 ]


# In[3]:

chains={}
dir='/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_FNAL_Feb2020/TOFHIR/RecoData/v1/RecoWithTracks'
for crys,rg in run_ranges.iteritems():
    chains[crys]=R.TChain("data")
    for run in range(rg[0],rg[1]+1):
        chains[crys].Add(dir+'/run'+str(run)+'_events.root')


# In[4]:


histos={}

for crys,c in chains.iteritems():
    
    histos[crys+'_'+'STEP1']=R.TH1F(crys+'_'+'STEP1',crys+'_'+'STEP1',10,0.5,10.5)
    histos[crys+'_'+'TH1']=R.TH1F(crys+'_'+'TH1',crys+'_'+'TH1',100,0.5,100.5)
    for ch in channels[crys]:
        histos[crys+'_'+'qfineVsToT_ch'+str(ch)]=R.TProfile(crys+'_'+'qfineVsToT_ch'+str(ch),crys+'_'+'qfineVsToT_ch'+str(ch),1000,0.5,1000.5)

#first loop to identify various steps and cleaning cuts 
    for ev in range(0,c.GetEntries()/4):
        c.GetEntry(ev)
        th1=(c.step2/10000)%100-1
        histos[crys+'_'+'STEP1'].Fill(c.step1)
        histos[crys+'_'+'TH1'].Fill(th1)
        for ch in channels[crys]:
           histos[crys+'_'+'qfineVsToT_ch'+str(ch)].Fill(c.tot[ch]/1000,c.qfine[ch]) 


# In[5]:


c1=R.TCanvas("c1","c1",600,400)


# In[6]:


ranges={}

for crys in chains.keys():
    ranges[crys]={}
    for ch in channels[crys]:
        ranges[crys+'_ch'+str(ch)]={}
        histos[crys+'_qfineVsToT_ch'+str(ch)].Fit("pol0","R","",50.5,196.5)
        ranges[crys+'_ch'+str(ch)]['qfineMin']=int(histos[crys+'_qfineVsToT_ch'+str(ch)].GetFunction("pol0").GetParameter(0))+2
        histos[crys+'_qfineVsToT_ch'+str(ch)].Fit("pol0","R","",650,800)
        ranges[crys+'_ch'+str(ch)]['qfineMax']=int(histos[crys+'_qfineVsToT_ch'+str(ch)].GetFunction("pol0").GetParameter(0))-20
        histos[crys+'_qfineVsToT_ch'+str(ch)].Fit("pol1","R","",250,600)
        ranges[crys+'_ch'+str(ch)]['qfineSlope']=histos[crys+'_qfineVsToT_ch'+str(ch)].GetFunction("pol1").GetParameter(1)


# In[7]:


print ranges


# In[8]:


ov={}
th1={}
for crys in chains.keys():
    ov[crys]=[]
    th1[crys]=[]
    for i in range(1,histos[crys+'_'+'STEP1'].GetNbinsX()+1):
        if (histos[crys+'_'+'STEP1'].GetBinContent(i)>0):
            ov[crys].append(histos[crys+'_'+'STEP1'].GetXaxis().GetBinCenter(i))
    for i in range(1,histos[crys+'_'+'TH1'].GetNbinsX()+1):
        if (histos[crys+'_'+'TH1'].GetBinContent(i)>0):
            th1[crys].append(histos[crys+'_'+'TH1'].GetXaxis().GetBinCenter(i))
   
    for v in ov[crys]:
        histos[crys+'_energyTot'+'_ov%.1f'%v]=R.TH1F(crys+'_energyTot'+'_ov%.1f'%v,crys+'_energyTot'+'_ov%.1f'%v,400,0.5,100.5)    
        for ch in channels[crys]:
            histos[crys+'_energy_ch'+str(ch)+'_ov%.1f'%v]=R.TH1F(crys+'_energy_ch'+str(ch)+'_ov%.1f'%v,crys+'_energy_ch'+str(ch)+'_ov%.1f'%v,400,0.5,100.5)
            histos[crys+'_tot_ch'+str(ch)+'_ov%.1f'%v]=R.TH1F(crys+'_tot_ch'+str(ch)+'_ov%.1f'%v,crys+'_tot_ch'+str(ch)+'_ov%.1f'%v,200,0.5,1000.5)


# In[ ]:


for crys,c in chains.iteritems():
#second loop to fill the energy histograms
    for ev in range(0,c.GetEntries()/4):
#    for ev in range(0,10):
        c.GetEntry(ev)
        o='_ov%.1f'%(c.step1)
        
        etot=0
        ech=0
        for ch in channels[crys]:
            if c.qfine[ch]<ranges[crys+'_ch'+str(ch)]['qfineMin']:
                continue
            if c.qfine[ch]>ranges[crys+'_ch'+str(ch)]['qfineMax']:
                continue
            if c.energy[ch]<0:
                continue
            ech=ech+1  
            etot=etot+c.energy[ch]
            histos[crys+'_energy_ch'+str(ch)+o].Fill(c.energy[ch])
            histos[crys+'_energy_ch'+str(ch)+o].Fill(c.energy[ch])
        if(ech>1):
            histos[crys+'_energyTot'+o].Fill(etot)


# In[ ]:


def langaufun(x,par):
    invsq2pi = 0.3989422804014
    mpshift  = -0.22278298
    
    np = 100.0
    sc =   5.0    
    ssum = 0.0


    mpc = par[1] - mpshift * par[0]
 
    xlow = x[0] - sc * par[3]
    xupp = x[0] + sc * par[3]
 
    step = (xupp-xlow) / np

    for i in range(1,int(np/2)+1):
        xx = xlow + (float(i)-0.5) * step
        fland = R.TMath.Landau(xx,mpc,par[0]) / par[0]
        ssum = ssum + fland * R.TMath.Gaus(x[0],xx,par[3])

        xx = xupp - (float(i)-0.5) * step
        fland = R.TMath.Landau(xx,mpc,par[0]) / par[0]
        ssum = ssum + fland * R.TMath.Gaus(x[0],xx,par[3])
     
    return (par[2] * step * ssum * invsq2pi / par[3])

def fitLanGaus(h):
    histos[h.GetName()+"_lanGausFit"]=R.TF1(h.GetName()+"_lanGausFit",langaufun,0,100,4)
    h.GetXaxis().SetRangeUser(0,100)
    hmax=h.GetXaxis().GetBinCenter(h.GetMaximumBin())
    hrms=h.GetRMS()
    histos[h.GetName()+"_lanGausFit"].SetParameter(0,1)
    histos[h.GetName()+"_lanGausFit"].SetParLimits(0,0.,2.)
    histos[h.GetName()+"_lanGausFit"].SetParameter(1,hmax)
    histos[h.GetName()+"_lanGausFit"].SetParameter(2,histos[crys+'_energy_ch'+str(ch)+ovS].GetMaximum())
    histos[h.GetName()+"_lanGausFit"].SetParameter(3,0.001)
    histos[h.GetName()+"_lanGausFit"].SetParLimits(3,0.001,2.)
    print(hmax,hrms)
    h.Fit(histos[h.GetName()+"_lanGausFit"],"RB","", hmax-hmax*0.15,hmax+0.8*hrms)
    h.GetXaxis().SetRangeUser(hmax-hmax*0.15,hmax+0.8*hrms)
    c1.SaveAs(h.GetName()+"_lanGausFit.png")


# In[ ]:


R.gStyle.SetOptFit(111111)
for crys in chains.keys():
    for v in ov[crys]:
        ovS='_ov%.1f'%(v)
        for ch in channels[crys]:
            fitLanGaus(histos[crys+'_energy_ch'+str(ch)+ovS])
        fitLanGaus(histos[crys+'_energyTot'+ovS])


# In[ ]:


for crys in chains.keys():
    histos[crys+'_energyVsOV']=R.TGraphErrors()
    histos[crys+'_energyVsOV'].SetName(crys+'_energyVsOV')
    for ch in channels[crys]:
        histos[crys+'_energy_ch'+str(ch)+'VsOV']=R.TGraphErrors()
        histos[crys+'_energy_ch'+str(ch)+'VsOV'].SetName(crys+'_energy_ch'+str(ch)+'VsOV')

    for v in ov[crys]:
        ovS='_ov%.1f'%(v)
        i=histos[crys+'_energyVsOV'].GetN()
        histos[crys+'_energyVsOV'].SetPoint(i,v,histos[crys+'_energyTot'+ovS+'_lanGausFit'].GetParameter(1))
        for ch in channels[crys]:
            i=histos[crys+'_energy_ch'+str(ch)+'VsOV'].GetN()
            histos[crys+'_energy_ch'+str(ch)+'VsOV'].SetPoint(i,v,histos[crys+'_energy_ch'+str(ch)+ovS+'_lanGausFit'].GetParameter(1))


# In[ ]:


R.gStyle.SetOptStat(0)
R.gStyle.SetOptTitle(0)

a=R.TH2F("a","a",10,2,10,10,0,50)
a.GetXaxis().SetTitle("OverVoltage [V]")
a.GetYaxis().SetTitle("MIP Peak [ADC]")
a.SetStats(0)

colors={0:R.kRed, 1:R.kBlue}
labels={0:'Left',1:'Right'}

l=R.TLatex()
l.SetTextSize(0.05)

for crys in chains.keys():
    a.Draw()
    leg=R.TLegend(0.15,0.65,0.4,0.85)
    leg.SetBorderSize(0)
    leg.SetFillColorAlpha(0,0)
    leg.SetTextSize(0.05)
#    leg.Set(crys)
    
    histos[crys+'_energyVsOV'].SetLineColor(R.kBlack)
    histos[crys+'_energyVsOV'].SetMarkerColor(R.kBlack)
    histos[crys+'_energyVsOV'].SetMarkerSize(1.2)
    histos[crys+'_energyVsOV'].SetMarkerStyle(20)
    histos[crys+'_energyVsOV'].Draw("LPSAME")
    leg.AddEntry(histos[crys+'_energyVsOV'],"Bar Energy SUM","LP")
    for ich,ch in enumerate(channels[crys]):
        histos[crys+'_energy_ch'+str(ch)+'VsOV'].SetLineColor(colors[ich])
        histos[crys+'_energy_ch'+str(ch)+'VsOV'].SetMarkerColor(colors[ich])
        histos[crys+'_energy_ch'+str(ch)+'VsOV'].SetMarkerSize(1.2)
        histos[crys+'_energy_ch'+str(ch)+'VsOV'].SetMarkerStyle(20)
        histos[crys+'_energy_ch'+str(ch)+'VsOV'].Draw("LPSAME")
        histos[crys+'_energy_ch'+str(ch)+'VsOV'].Draw("LPSAME")
        leg.AddEntry(histos[crys+'_energy_ch'+str(ch)+'VsOV'],"Bar Energy "+labels[ich],"LP")
    leg.Draw()
    l.DrawLatexNDC(0.12,0.93,crys+' - S12572-015C - protons 120 GeV')
    c1.SaveAs(crys+"_energyVsOV.png")


# In[ ]:


for crys in chains.keys():
    for ch in channels[crys]:
        for v in ov[crys]:
            ovS='_ov%.1f'%(v)
            ranges[crys+'_ch'+str(ch)]['eMin'+ovS]=histos[crys+'_energy_ch'+str(ch)+ovS+'_lanGausFit'].GetParameter(1)*0.9
            ranges[crys+'_ch'+str(ch)]['eMax'+ovS]=histos[crys+'_energy_ch'+str(ch)+ovS+'_lanGausFit'].GetParameter(1)*2
print ranges


# In[ ]:

R.gROOT.ProcessLine(".L timeAnalysis.C+")

def Map(tf):
    m = {}
    for k in tf.GetListOfKeys():
        n = k.GetName()
        m[n] = tf.Get(n)
    return m

for crys,c in chains.items():
    t=R.timeAnalysis(c)
    t.outFileName="timeAnalysis_%s.root"%crys
    t.xMin=19.5-3;
    t.xMax=19.5+3;
    t.energyMin.resize(len(channels[crys]))
    t.energyMax.resize(len(channels[crys]))
    for ich,ch in enumerate(channels[crys]):
        t.channels.push_back(ch)
        t.qfineMin.push_back(ranges[crys+'_ch'+str(ch)]['qfineMin'])
        t.qfineMax.push_back(ranges[crys+'_ch'+str(ch)]['qfineMax'])
        for v in ov[crys]:
            ovS='%.1f'%(v)
            t.energyMin[ich][ovS]=ranges[crys+'_ch'+str(ch)]['eMin_ov'+ovS]
            t.energyMax[ich][ovS]=ranges[crys+'_ch'+str(ch)]['eMax_ov'+ovS]
    t.Loop()

    fTime=R.TFile("timeAnalysis_%s.root"%crys)
                       
#    for v in ov[crys]:
#        ovS='_ov%.1f'%(v)
#        for th in th1[crys]:
#            thS='_th%d'%th
#                       histos[crys+'_tDiff'+ovS+thS]=R.TH1F('tDiff'+ovS+thS,'tDiff'+ovS+thS).Clone(crys+'_tDiff'+ovS+thS)


# In[ ]:


fOut=R.TFile("histos.root","RECREATE")
for hN,h in histos.iteritems():
    h.Write()
fOut.Write()
fOut.Close()

