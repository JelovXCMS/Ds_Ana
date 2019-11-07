//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 27 17:43:54 2018 by ROOT version 6.02/13
// from TTree t_DsMassData_pt20to40/t_DsMassData_pt20to40
// found on file: FitFile_pp.root
//////////////////////////////////////////////////////////

#ifndef t_DsMassData_pt20to40_h
#define t_DsMassData_pt20to40_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class t_DsMassData_pt20to40 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         Dmass;
   Float_t         Ddca;

   // List of branches
   TBranch        *b_Dmass;   //!
   TBranch        *b_Ddca;   //!

   t_DsMassData_pt20to40(TTree *tree=0);
   virtual ~t_DsMassData_pt20to40();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef t_DsMassData_pt20to40_cxx
t_DsMassData_pt20to40::t_DsMassData_pt20to40(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("FitFile_pp.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("FitFile_pp.root");
      }
      f->GetObject("t_DsMassData_pt20to40",tree);

   }
   Init(tree);
}

t_DsMassData_pt20to40::~t_DsMassData_pt20to40()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t t_DsMassData_pt20to40::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t t_DsMassData_pt20to40::LoadTree(Long64_t entry)
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

void t_DsMassData_pt20to40::Init(TTree *tree)
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

   fChain->SetBranchAddress("Dmass", &Dmass, &b_Dmass);
   fChain->SetBranchAddress("Ddca", &Ddca, &b_Ddca);
   Notify();
}

Bool_t t_DsMassData_pt20to40::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void t_DsMassData_pt20to40::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t t_DsMassData_pt20to40::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef t_DsMassData_pt20to40_cxx
