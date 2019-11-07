#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"

#include <TString.h>
#include <TCanvas.h>
#include "TBranch.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLatex.h"
#include <TH2F.h>
#include <TCut.h>
#include <TStyle.h>

#include <math.h>

#include <string>
#include <TChain.h>
#include <TFile.h>
#include "TH1.h"
#include "TROOT.h"
#include "TPad.h"
#include "TRandom.h"
#include "THStack.h"
#include "TH2.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TAxis.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <TCut.h>

using namespace std;


int bfraction_CMSGen7TeV(){

	TFile *fout=new TFile("bfraction7TeV.root","recreate");


	TChain *ch_B7TeV=new TChain("Events");
	
	ch_B7TeV->Add("/mnt/hadoop/store/group/hi/chengchi/Btest/Btest7TeV/Btest7TeV_step0_GEN/Btest7TeV_step0_GEN/180717_041845/0000/*.root");

	fout->cd();

	TH1D *h_B0_7TeV=new TH1D("h_B0_7TeV","h_B0_7TeV",18,-4.5,4.5); h_B0_7TeV->Sumw2();
	TH1D *h_Bp_7TeV=new TH1D("h_Bp_7TeV","h_B0_7TeV",18,-4.5,4.5); h_Bp_7TeV->Sumw2();
	TH1D *h_Bs_7TeV=new TH1D("h_Bs_7TeV","h_B0_7TeV",18,-4.5,4.5); h_Bs_7TeV->Sumw2();

	TH1D *h_BAll_7TeV=new TH1D("h_BAll_7TeV","h_BAll_7TeV",18,-4.5,4.5); h_BAll_7TeV->Sumw2();

	TH1D *h_B0_7TeV_frac=new TH1D("h_B0_7TeV_frac","h_B0_7TeV_frac",18,-4.5,4.5); h_B0_7TeV_frac->Sumw2();
	TH1D *h_Bp_7TeV_frac=new TH1D("h_Bp_7TeV_frac","h_B0_7TeV_frac",18,-4.5,4.5); h_Bp_7TeV_frac->Sumw2();
	TH1D *h_Bs_7TeV_frac=new TH1D("h_Bs_7TeV_frac","h_B0_7TeV_frac",18,-4.5,4.5); h_Bs_7TeV_frac->Sumw2();

	ch_B7TeV->Draw("recoGenParticles_genParticles__GEN.obj.y()>>h_B0_7TeV","abs(recoGenParticles_genParticles__GEN.obj.pdgId())==511");
	ch_B7TeV->Draw("recoGenParticles_genParticles__GEN.obj.y()>>h_Bp_7TeV","abs(recoGenParticles_genParticles__GEN.obj.pdgId())==521");
	ch_B7TeV->Draw("recoGenParticles_genParticles__GEN.obj.y()>>h_Bs_7TeV","abs(recoGenParticles_genParticles__GEN.obj.pdgId())==531");

	h_BAll_7TeV->Add(h_B0_7TeV);
	h_BAll_7TeV->Add(h_Bp_7TeV);
	h_BAll_7TeV->Add(h_Bs_7TeV);

	h_B0_7TeV_frac->Add(h_B0_7TeV);
	h_Bp_7TeV_frac->Add(h_Bp_7TeV);
	h_Bs_7TeV_frac->Add(h_Bs_7TeV);

	h_B0_7TeV_frac->Divide(h_B0_7TeV_frac,h_BAll_7TeV,1,1,"B");
	h_Bp_7TeV_frac->Divide(h_Bp_7TeV_frac,h_BAll_7TeV,1,1,"B");
	h_Bs_7TeV_frac->Divide(h_Bs_7TeV_frac,h_BAll_7TeV,1,1,"B");


	TCanvas *c_Bfr_7TeV=new TCanvas("c_Bfr_7TeV","c_Bfr_7TeV",800,600);
	c_Bfr_7TeV->cd();
	h_B0_7TeV_frac->SetLineColor(1);
	h_B0_7TeV_frac->Draw("same");
	h_Bp_7TeV_frac->SetLineColor(2);
	h_Bp_7TeV_frac->Draw("same");
	h_Bs_7TeV_frac->SetLineColor(4);
	h_Bs_7TeV_frac->Draw("same");



	h_B0_7TeV->Write();
	h_Bp_7TeV->Write();
	h_Bs_7TeV->Write();

	h_B0_7TeV_frac->Write();
	h_Bp_7TeV_frac->Write();
	h_Bs_7TeV_frac->Write();





	return 0;

}
