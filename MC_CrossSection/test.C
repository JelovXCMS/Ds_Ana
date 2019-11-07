using namespace std;
#include <memory>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <list>
#include <iomanip>
#include <cmath>

#include <TLorentzVector.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include <TString.h>
#include <TNtuple.h>
#include <TVector3.h>


int CSProject(TString infile="", TString outfile="", Bool_t REAL=false, Int_t isPbPb=0){

  TString ifname;
	ifname = infile;
  cout<<"ifname = "<<ifname<<endl;
  if (!TFile::Open(ifname))   { cout << " fail to open file" << endl; return 0;}
  TFile* fin = TFile::Open(ifname);


	TTree *t_evt=(TTree*)fin->Get("Events");

	int nbin=10;
	double bins[nbin+1]={0,2,3,4,5,6,8,10,20,40,200};

	
  TFile* fout = TFile::Open(Form("%s", outfile.Data()),"recreate");
	fout->cd();	
	TH1D *hGpt_D0=new TH1D("hGpt_D0","hGpt_D0",nbin,bins); hGpt_D0->Sumw2();
	TH1D *hGpt_Ds=new TH1D("hGpt_Ds","hGpt_Ds",nbin,bins); hGpt_Ds->Sumw2();
	

	t_evt->Project("hGpt_D0","recoGenParticles_genParticles__GEN.obj.pt()","abs(recoGenParticles_genParticles__GEN.obj.pdgId())==421&&abs(recoGenParticles_genParticles__GEN.obj.y())<1");
	t_evt->Project("hGpt_Ds","recoGenParticles_genParticles__GEN.obj.pt()","abs(recoGenParticles_genParticles__GEN.obj.pdgId())==431&&abs(recoGenParticles_genParticles__GEN.obj.y())<1");

// adding cut: recoGenParticles_genParticles__GEN.obj.mother().pdgid()<500  , and need to consider its mother's mother...

//root [4] Events->Draw("recoGenParticles_genParticles__GEN.obj.pdgId()","abs(recoGenParticles_genParticles__GEN.obj.pdgId())==421")
// (Long64_t) 2998
//root [5] Events->Draw("recoGenParticles_genParticles__GEN.obj.mother().pdgId()","abs(recoGenParticles_genParticles__GEN.obj.p^[[A^[[B^[[B^[[B(Long64_t) 2998
//root [6] Events->Draw("recoGenParticles_genParticles__GEN.obj.mother().pdgId()","abs(recoGenParticles_genParticles__GEN.obj.pdgId())==431")
//(Long64_t) 513
//root [7] Events->Draw("recoGenParticles_genParticles__GEN.obj.mother().mother().pdgId()","abs(recoGenParticles_genParticles__^[[A^[[B^[[B^[[A^[[B(Long64_t) 513
//root [8] Events->Draw("recoGenParticles_genParticles__GEN.obj.mother().mother().pdgId()","abs(recoGenParticles_genParticles__GEN.obj.pdgId())==431 && abs(recoGenParticles_genParticles__GEN.obj.mother().mother().pdgId())<600")
//(Long64_t) 508


	hGpt_D0->Write();
	hGpt_Ds->Write();

	fout->Close();

	return 0;
}

int main(int argc, char *argv[])
{
  if(argc==3)
  {
    CSProject(argv[1], argv[2]);
  }
  else if(argc==5)
  {
    cout<<"real = "<<atoi(argv[3])<<" ,isPbPb = "<<atoi(argv[4])<<endl;
    CSProject(argv[1], argv[2],atoi(argv[3]), atoi(argv[4]));
  }
  else
  {
    std::cout << "Usage: loop.exe <input_collection> <output_file> <isReal> <isPbPb>" << std::endl;
    return 1;
  }

  return 0;
}

