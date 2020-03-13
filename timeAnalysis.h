//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 11 11:52:31 2020 by ROOT version 6.18/04
// from TTree data/data
// found on file: /afs/cern.ch/user/m/meridian/eosmtd/comm_mtd/TB/MTDTB_FNAL_Feb2020/TOFHIR/RecoData/v1/RecoWithTracks/run20912_events.root
//////////////////////////////////////////////////////////

#ifndef timeAnalysis_h
#define timeAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>

#include <vector>

// Header file for the classes stored in the TTree if any.

class timeAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   std::string    outFileName;
   TFile          *fOut;

   float xMin;
   float xMax;

   std::vector<int> channels;
   std::vector<float> qfineMin;
   std::vector<float> qfineMax;
   std::vector<std::map<std::string,float>> energyMin;
   std::vector<std::map<std::string,float>> energyMax;

   std::map<std::string,TH1F*> histos;


// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Float_t         step1;
   Float_t         step2;
   Long64_t        time[256];
   Float_t         tot[256];
   Float_t         energy[256];
   Float_t         qfine[256];
   Float_t         xIntercept;
   Float_t         yIntercept;
   Float_t         xSlope;
   Float_t         ySlope;
   Float_t         x_dut;
   Float_t         y_dut;
   Float_t         chi2;
   Int_t           ntracks;
   Int_t           nplanes;
   Int_t           matchEff;
   Int_t           SlowTriggerTag;
   Long64_t        tDiffTrigger[150];
   Float_t         TimeDiff;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_step1;   //!
   TBranch        *b_step2;   //!
   TBranch        *b_time;   //!
   TBranch        *b_tot;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_qfine;   //!
   TBranch        *b_xIntercept;   //!
   TBranch        *b_yIntercept;   //!
   TBranch        *b_xSlope;   //!
   TBranch        *b_ySlope;   //!
   TBranch        *b_x_dut;   //!
   TBranch        *b_y_dut;   //!
   TBranch        *b_chi2;   //!
   TBranch        *b_ntracks;   //!
   TBranch        *b_nplanes;   //!
   TBranch        *b_matchEff;   //!
   TBranch        *b_SlowTriggerTag;   //!
   TBranch        *b_tDiffTrigger;   //!
   TBranch        *b_TimeDiff;   //!

   timeAnalysis(TTree *tree=0);
   virtual ~timeAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef timeAnalysis_cxx
timeAnalysis::timeAnalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/afs/cern.ch/user/m/meridian/eosmtd/comm_mtd/TB/MTDTB_FNAL_Feb2020/TOFHIR/RecoData/v1/RecoWithTracks/run20912_events.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/afs/cern.ch/user/m/meridian/eosmtd/comm_mtd/TB/MTDTB_FNAL_Feb2020/TOFHIR/RecoData/v1/RecoWithTracks/run20912_events.root");
      }
      f->GetObject("data",tree);

   }
   Init(tree);
}

timeAnalysis::~timeAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t timeAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t timeAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void timeAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("step1", &step1, &b_step1);
   fChain->SetBranchAddress("step2", &step2, &b_step2);
   fChain->SetBranchAddress("time", time, &b_time);
   fChain->SetBranchAddress("tot", tot, &b_tot);
   fChain->SetBranchAddress("energy", energy, &b_energy);
   fChain->SetBranchAddress("qfine", qfine, &b_qfine);
   fChain->SetBranchAddress("xIntercept", &xIntercept, &b_xIntercept);
   fChain->SetBranchAddress("yIntercept", &yIntercept, &b_yIntercept);
   fChain->SetBranchAddress("xSlope", &xSlope, &b_xSlope);
   fChain->SetBranchAddress("ySlope", &ySlope, &b_ySlope);
   fChain->SetBranchAddress("x_dut", &x_dut, &b_x_dut);
   fChain->SetBranchAddress("y_dut", &y_dut, &b_y_dut);
   fChain->SetBranchAddress("chi2", &chi2, &b_chi2);
   fChain->SetBranchAddress("ntracks", &ntracks, &b_ntracks);
   fChain->SetBranchAddress("nplanes", &nplanes, &b_nplanes);
   fChain->SetBranchAddress("matchEff", &matchEff, &b_matchEff);
   fChain->SetBranchAddress("SlowTriggerTag", &SlowTriggerTag, &b_SlowTriggerTag);
   fChain->SetBranchAddress("tDiffTrigger", tDiffTrigger, &b_tDiffTrigger);
   fChain->SetBranchAddress("TimeDiff", &TimeDiff, &b_TimeDiff);
   Notify();
}

Bool_t timeAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void timeAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t timeAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef timeAnalysis_cxx
