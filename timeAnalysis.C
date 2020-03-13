#define timeAnalysis_cxx
#include "timeAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>

void timeAnalysis::Loop()
{
//   In a ROOT session, you can do:
//      root> .L timeAnalysis.C
//      root> timeAnalysis t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   int thresholds[] = { 10, 20, 40 };
   float overvoltages[] = { 4.0, 6.0, 8.0 };

   for (auto& th: thresholds)
     for (auto& ov: overvoltages)
	 histos[Form("tDiff_ov%.1f_th%d",ov,th)]=new TH1F(Form("tDiff_ov%.1f_th%d",ov,th),Form("tDiff_ov%.1f_th%d",ov,th),200,-2000,2000);

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //   for (Long64_t jentry=0; jentry<nentries/4;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      if (jentry%1000==0)
	std::cout << "Event #"<< jentry << std::endl;

      nb = fChain->GetEntry(jentry);   nbytes += nb;
      int ich=-1;
      int nch=0;
      int th1=int(step2/10000)%100-1;

      if (step1!=4 && step1!=6 && step1!=8)
	continue;

      if (th1!=10 && th1!=20 && th1!=40)
	continue;

      if (ntracks!=1)
	continue;

      if (x_dut<xMin)
	continue;

      if (x_dut>xMax)
	continue;
      
      for (auto& ch: channels)
	{
	  ++ich;
	  if (qfine[ch]<qfineMin[ich])
	    continue;
	  if (qfine[ch]>qfineMax[ich])
	    continue;
	  if (energy[ch]<energyMin[ich][Form("%.1f",step1)])
	    continue;
	  if (energy[ch]>energyMax[ich][Form("%.1f",step1)])
	    continue;
	  ++nch;
	}
      
      if (nch<2) //requires at 2 channels passing cuts
	continue;

      histos[Form("tDiff_ov%.1f_th%d",step1,th1)]->Fill((time[channels[0]]-time[channels[1]]));
   }

   fOut=new TFile(outFileName.c_str(),"RECREATE");
   for (auto& h: histos)
     h.second->Write();
   fOut->Write();
   fOut->Close();
     
}
