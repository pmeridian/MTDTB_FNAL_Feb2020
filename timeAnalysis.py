import ROOT as R

R.gROOT.ProcessLine(".L timeAnalysis.C++")

c=R.TChain('data')
c.Add('/afs/cern.ch/user/m/meridian/eosmtd/comm_mtd/TB/MTDTB_FNAL_Feb2020/TOFHIR/RecoData/v1/RecoWithTracks/run2092*_events.root')

t=R.timeAnalysis(c)
t.outFileName="timeAnalysis.root"

t.xMin=19.5-2
t.xMax=19.5+2

t.channels.push_back(213)
t.channels.push_back(219)

t.qfineMin.push_back(15)
t.qfineMin.push_back(15)

t.qfineMax.push_back(565)
t.qfineMax.push_back(565)

t.energyMin.resize(t.channels.size())
t.energyMin[0]['4.0']=6
t.energyMin[0]['6.0']=6
t.energyMin[0]['8.0']=6
t.energyMin[1]['4.0']=6
t.energyMin[1]['6.0']=6
t.energyMin[1]['8.0']=6

t.energyMax.resize(t.channels.size())
t.energyMax[0]['4.0']=15
t.energyMax[0]['6.0']=15
t.energyMax[0]['8.0']=15
t.energyMax[1]['4.0']=15
t.energyMax[1]['6.0']=15
t.energyMax[1]['8.0']=15

t.Loop()
