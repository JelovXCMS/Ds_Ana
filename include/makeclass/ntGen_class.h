//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb 26 17:53:53 2018 by ROOT version 6.02/13
// from TTree ntGen/
// found on file: Dntuple_finder_PbPb_1111.root_1104.root
//////////////////////////////////////////////////////////

#ifndef ntGen_class_h
#define ntGen_class_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ntGen_class {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         GPVx;
   Float_t         GPVy;
   Float_t         GPVz;
   Int_t           Gsize;
   Float_t         Gy[1];   //[Gsize]
   Float_t         Geta[1];   //[Gsize]
   Float_t         Gphi[1];   //[Gsize]
   Float_t         Gpt[1];   //[Gsize]
   Int_t           GpdgId[1];   //[Gsize]
   Int_t           GcollisionId[1];   //[Gsize]
   Int_t           GisSignal[1];   //[Gsize]
   Float_t         GBAncestorpt[1];   //[Gsize]
   Int_t           GBAncestorpdgId[1];   //[Gsize]
   Int_t           GfromgenPV[1];   //[Gsize]
   Float_t         GprodvtxX[1];   //[Gsize]
   Float_t         GprodvtxY[1];   //[Gsize]
   Float_t         GprodvtxZ[1];   //[Gsize]
   Float_t         GdecayvtxX[1];   //[Gsize]
   Float_t         GdecayvtxY[1];   //[Gsize]
   Float_t         GdecayvtxZ[1];   //[Gsize]
   Float_t         Gtk1pt[1];   //[Gsize]
   Float_t         Gtk1eta[1];   //[Gsize]
   Float_t         Gtk1y[1];   //[Gsize]
   Float_t         Gtk1phi[1];   //[Gsize]
   Int_t           Gtk1pdgId[1];   //[Gsize]
   Float_t         Gtk2pt[1];   //[Gsize]
   Float_t         Gtk2eta[1];   //[Gsize]
   Float_t         Gtk2y[1];   //[Gsize]
   Float_t         Gtk2phi[1];   //[Gsize]
   Int_t           Gtk2pdgId[1];   //[Gsize]
   Float_t         Gtk3pt[1];   //[Gsize]
   Float_t         Gtk3eta[1];   //[Gsize]
   Float_t         Gtk3y[1];   //[Gsize]
   Float_t         Gtk3phi[1];   //[Gsize]
   Float_t         Gtk4pt[1];   //[Gsize]
   Float_t         Gtk4eta[1];   //[Gsize]
   Float_t         Gtk4y[1];   //[Gsize]
   Float_t         Gtk4phi[1];   //[Gsize]
   Float_t         GRestk1pt[1];   //[Gsize]
   Float_t         GRestk1eta[1];   //[Gsize]
   Float_t         GRestk1y[1];   //[Gsize]
   Float_t         GRestk1phi[1];   //[Gsize]
   Int_t           GRestk1pdgId[1];   //[Gsize]
   Float_t         GRestk2pt[1];   //[Gsize]
   Float_t         GRestk2eta[1];   //[Gsize]
   Float_t         GRestk2y[1];   //[Gsize]
   Float_t         GRestk2phi[1];   //[Gsize]
   Int_t           GRestk2pdgId[1];   //[Gsize]
   Float_t         GRestk3pt[1];   //[Gsize]
   Float_t         GRestk3eta[1];   //[Gsize]
   Float_t         GRestk3y[1];   //[Gsize]
   Float_t         GRestk3phi[1];   //[Gsize]
   Float_t         GRestk4pt[1];   //[Gsize]
   Float_t         GRestk4eta[1];   //[Gsize]
   Float_t         GRestk4y[1];   //[Gsize]
   Float_t         GRestk4phi[1];   //[Gsize]

   // List of branches
   TBranch        *b_GPVx;   //!
   TBranch        *b_GPVy;   //!
   TBranch        *b_GPVz;   //!
   TBranch        *b_Gsize;   //!
   TBranch        *b_Gy;   //!
   TBranch        *b_Geta;   //!
   TBranch        *b_Gphi;   //!
   TBranch        *b_Gpt;   //!
   TBranch        *b_GpdgId;   //!
   TBranch        *b_GcollisionId;   //!
   TBranch        *b_GisSignal;   //!
   TBranch        *b_GBAncestorpt;   //!
   TBranch        *b_GBAncestorpdgId;   //!
   TBranch        *b_GfromgenPV;   //!
   TBranch        *b_GprodvtxX;   //!
   TBranch        *b_GprodvtxY;   //!
   TBranch        *b_GprodvtxZ;   //!
   TBranch        *b_GdecayvtxX;   //!
   TBranch        *b_GdecayvtxY;   //!
   TBranch        *b_GdecayvtxZ;   //!
   TBranch        *b_Gtk1pt;   //!
   TBranch        *b_Gtk1eta;   //!
   TBranch        *b_Gtk1y;   //!
   TBranch        *b_Gtk1phi;   //!
   TBranch        *b_Gtk1pdgId;   //!
   TBranch        *b_Gtk2pt;   //!
   TBranch        *b_Gtk2eta;   //!
   TBranch        *b_Gtk2y;   //!
   TBranch        *b_Gtk2phi;   //!
   TBranch        *b_Gtk2pdgId;   //!
   TBranch        *b_Gtk3pt;   //!
   TBranch        *b_Gtk3eta;   //!
   TBranch        *b_Gtk3y;   //!
   TBranch        *b_Gtk3phi;   //!
   TBranch        *b_Gtk4pt;   //!
   TBranch        *b_Gtk4eta;   //!
   TBranch        *b_Gtk4y;   //!
   TBranch        *b_Gtk4phi;   //!
   TBranch        *b_GRestk1pt;   //!
   TBranch        *b_GRestk1eta;   //!
   TBranch        *b_GRestk1y;   //!
   TBranch        *b_GRestk1phi;   //!
   TBranch        *b_GRestk1pdgId;   //!
   TBranch        *b_GRestk2pt;   //!
   TBranch        *b_GRestk2eta;   //!
   TBranch        *b_GRestk2y;   //!
   TBranch        *b_GRestk2phi;   //!
   TBranch        *b_GRestk2pdgId;   //!
   TBranch        *b_GRestk3pt;   //!
   TBranch        *b_GRestk3eta;   //!
   TBranch        *b_GRestk3y;   //!
   TBranch        *b_GRestk3phi;   //!
   TBranch        *b_GRestk4pt;   //!
   TBranch        *b_GRestk4eta;   //!
   TBranch        *b_GRestk4y;   //!
   TBranch        *b_GRestk4phi;   //!

   ntGen_class(TTree *tree=0);
   virtual ~ntGen_class();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ntGen_class_cxx
ntGen_class::ntGen_class(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Dntuple_finder_PbPb_1111.root_1104.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Dntuple_finder_PbPb_1111.root_1104.root");
      }
      f->GetObject("ntGen",tree);

   }
   Init(tree);
}

ntGen_class::~ntGen_class()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ntGen_class::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ntGen_class::LoadTree(Long64_t entry)
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

void ntGen_class::Init(TTree *tree)
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

   fChain->SetBranchAddress("GPVx", &GPVx, &b_GPVx);
   fChain->SetBranchAddress("GPVy", &GPVy, &b_GPVy);
   fChain->SetBranchAddress("GPVz", &GPVz, &b_GPVz);
   fChain->SetBranchAddress("Gsize", &Gsize, &b_Gsize);
   fChain->SetBranchAddress("Gy", &Gy, &b_Gy);
   fChain->SetBranchAddress("Geta", &Geta, &b_Geta);
   fChain->SetBranchAddress("Gphi", &Gphi, &b_Gphi);
   fChain->SetBranchAddress("Gpt", &Gpt, &b_Gpt);
   fChain->SetBranchAddress("GpdgId", &GpdgId, &b_GpdgId);
   fChain->SetBranchAddress("GcollisionId", &GcollisionId, &b_GcollisionId);
   fChain->SetBranchAddress("GisSignal", &GisSignal, &b_GisSignal);
   fChain->SetBranchAddress("GBAncestorpt", &GBAncestorpt, &b_GBAncestorpt);
   fChain->SetBranchAddress("GBAncestorpdgId", &GBAncestorpdgId, &b_GBAncestorpdgId);
   fChain->SetBranchAddress("GfromgenPV", &GfromgenPV, &b_GfromgenPV);
   fChain->SetBranchAddress("GprodvtxX", &GprodvtxX, &b_GprodvtxX);
   fChain->SetBranchAddress("GprodvtxY", &GprodvtxY, &b_GprodvtxY);
   fChain->SetBranchAddress("GprodvtxZ", &GprodvtxZ, &b_GprodvtxZ);
   fChain->SetBranchAddress("GdecayvtxX", &GdecayvtxX, &b_GdecayvtxX);
   fChain->SetBranchAddress("GdecayvtxY", &GdecayvtxY, &b_GdecayvtxY);
   fChain->SetBranchAddress("GdecayvtxZ", &GdecayvtxZ, &b_GdecayvtxZ);
   fChain->SetBranchAddress("Gtk1pt", &Gtk1pt, &b_Gtk1pt);
   fChain->SetBranchAddress("Gtk1eta", &Gtk1eta, &b_Gtk1eta);
   fChain->SetBranchAddress("Gtk1y", &Gtk1y, &b_Gtk1y);
   fChain->SetBranchAddress("Gtk1phi", &Gtk1phi, &b_Gtk1phi);
   fChain->SetBranchAddress("Gtk1pdgId", &Gtk1pdgId, &b_Gtk1pdgId);
   fChain->SetBranchAddress("Gtk2pt", &Gtk2pt, &b_Gtk2pt);
   fChain->SetBranchAddress("Gtk2eta", &Gtk2eta, &b_Gtk2eta);
   fChain->SetBranchAddress("Gtk2y", &Gtk2y, &b_Gtk2y);
   fChain->SetBranchAddress("Gtk2phi", &Gtk2phi, &b_Gtk2phi);
   fChain->SetBranchAddress("Gtk2pdgId", &Gtk2pdgId, &b_Gtk2pdgId);
   fChain->SetBranchAddress("Gtk3pt", &Gtk3pt, &b_Gtk3pt);
   fChain->SetBranchAddress("Gtk3eta", &Gtk3eta, &b_Gtk3eta);
   fChain->SetBranchAddress("Gtk3y", &Gtk3y, &b_Gtk3y);
   fChain->SetBranchAddress("Gtk3phi", &Gtk3phi, &b_Gtk3phi);
   fChain->SetBranchAddress("Gtk4pt", &Gtk4pt, &b_Gtk4pt);
   fChain->SetBranchAddress("Gtk4eta", &Gtk4eta, &b_Gtk4eta);
   fChain->SetBranchAddress("Gtk4y", &Gtk4y, &b_Gtk4y);
   fChain->SetBranchAddress("Gtk4phi", &Gtk4phi, &b_Gtk4phi);
   fChain->SetBranchAddress("GRestk1pt", &GRestk1pt, &b_GRestk1pt);
   fChain->SetBranchAddress("GRestk1eta", &GRestk1eta, &b_GRestk1eta);
   fChain->SetBranchAddress("GRestk1y", &GRestk1y, &b_GRestk1y);
   fChain->SetBranchAddress("GRestk1phi", &GRestk1phi, &b_GRestk1phi);
   fChain->SetBranchAddress("GRestk1pdgId", &GRestk1pdgId, &b_GRestk1pdgId);
   fChain->SetBranchAddress("GRestk2pt", &GRestk2pt, &b_GRestk2pt);
   fChain->SetBranchAddress("GRestk2eta", &GRestk2eta, &b_GRestk2eta);
   fChain->SetBranchAddress("GRestk2y", &GRestk2y, &b_GRestk2y);
   fChain->SetBranchAddress("GRestk2phi", &GRestk2phi, &b_GRestk2phi);
   fChain->SetBranchAddress("GRestk2pdgId", &GRestk2pdgId, &b_GRestk2pdgId);
   fChain->SetBranchAddress("GRestk3pt", &GRestk3pt, &b_GRestk3pt);
   fChain->SetBranchAddress("GRestk3eta", &GRestk3eta, &b_GRestk3eta);
   fChain->SetBranchAddress("GRestk3y", &GRestk3y, &b_GRestk3y);
   fChain->SetBranchAddress("GRestk3phi", &GRestk3phi, &b_GRestk3phi);
   fChain->SetBranchAddress("GRestk4pt", &GRestk4pt, &b_GRestk4pt);
   fChain->SetBranchAddress("GRestk4eta", &GRestk4eta, &b_GRestk4eta);
   fChain->SetBranchAddress("GRestk4y", &GRestk4y, &b_GRestk4y);
   fChain->SetBranchAddress("GRestk4phi", &GRestk4phi, &b_GRestk4phi);
   Notify();
}

Bool_t ntGen_class::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ntGen_class::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ntGen_class::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ntGen_class_cxx
