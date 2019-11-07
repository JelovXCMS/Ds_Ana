// adding singleRatio 

#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"


#include <TString.h>
#include <TCanvas.h>
#include "TBranch.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLatex.h"
#include <TH2F.h>
#include <TCut.h>
#include <TStyle.h>
#include <TLine.h>
#include <TGraph.h>
#include <TGraphErrors.h>

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

#include "TFitter.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"


// using namespace RooFit;
using namespace std;

// Int_t StartBin=0;

TCanvas *cav[500];
TCanvas *cav1[500];
Int_t c_count=0;

double shiftY=0;
double oneshift=0.075;


int CutScanSys_FixShape(Int_t isPbPb=0, Int_t startbin=0, Int_t start_var=0){
  if(isPbPb==3){startbin=2;}
  // if(isPbPb==3){startbin=4;}
  if(isPbPb==0){startbin=0;}
  // StartBin=startbin;

	TFitResultPtr fitR;
	double x0[1]={0};
	double x0err[1];


	bool doDauPtCut=0;

  InitStyle();
  initParameter();

	double RatioMax=2.5;
	double RatioMin=0;


	gStyle->SetOptStat(0);

  TString str_PbPb="pp";
	TString str_PbPbtext="pp";
  int nbin_pt=nbin_pt_pp;
  double *bins_pt=bins_pt_pp;

  double LumiNevt=LumiSum;
  TString str_eff_Prompt_phi=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_pp_Prompt_phikkpi.root",s_CutSet.Data());
  TString str_eff_NonPrompt_phi=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_pp_NonPrompt_phikkpi.root",s_CutSet.Data());
  TString str_eff_Prompt_f0=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_pp_Prompt_f0kkpi.root",s_CutSet.Data());
  TString str_eff_NonPrompt_f0=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_pp_NonPrompt_f0kkpi.root",s_CutSet.Data());

  if(isPbPb==3){
    LumiNevt=NevtPbPb3;
    nbin_pt=nbin_pt_PbPb3;
    bins_pt=bins_pt_PbPb3;
    str_eff_Prompt_phi=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_PbPb3_Prompt_phikkpi.root",s_CutSet.Data());
    str_eff_NonPrompt_phi=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_PbPb3_NonPrompt_phikkpi.root",s_CutSet.Data());
    str_eff_Prompt_f0=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_PbPb3_Prompt_f0kkpi.root",s_CutSet.Data());
    str_eff_NonPrompt_f0=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_PbPb3_NonPrompt_f0kkpi.root",s_CutSet.Data());
    str_PbPb="PbPb3";
  }

  // reading fit raw yield

	cout<<__LINE__<<endl;
  TString s_f_raw=Form("../SignalFit/output%s/RawFitYield_FixShapeTrkPtScan_%s.root",s_CutSet.Data(),str_PbPb.Data());

  TFile *f_raw=TFile::Open(s_f_raw.Data(),"read");
	TH1F *h_RawFitYield=(TH1F*)f_raw->Get("h_RawBinFitYield");


	cout<<__LINE__<<endl;
 // for Efficiency histogram
  TFile *f_Prompt_phikkpi= TFile::Open(str_eff_Prompt_phi.Data());
  TFile *f_NonPrompt_phikkpi= TFile::Open(str_eff_NonPrompt_phi.Data());
  TFile *f_Prompt_f0kkpi= TFile::Open(str_eff_Prompt_f0.Data());
  TFile *f_NonPrompt_f0kkpi= TFile::Open(str_eff_NonPrompt_f0.Data());

	TString str_effH="h_RecoNormEff";

	TH1D *h_Eff_Prompt_phikkpi=(TH1D*)f_Prompt_phikkpi->Get(str_effH.Data());
	TH1D *h_Eff_NonPrompt_phikkpi=(TH1D*)f_NonPrompt_phikkpi->Get(str_effH.Data());
	TH1D *h_Eff_Prompt_f0kkpi=(TH1D*)f_Prompt_f0kkpi->Get(str_effH.Data());
	TH1D *h_Eff_NonPrompt_f0kkpi=(TH1D*)f_NonPrompt_f0kkpi->Get(str_effH.Data());

	TH1D *h_Eff_Prompt_AllBR=(TH1D*)h_Eff_Prompt_phikkpi->Clone("h_Eff_Prompt_AllBR");
	h_Eff_Prompt_AllBR->Sumw2();
	h_Eff_Prompt_AllBR->Add(h_Eff_Prompt_phikkpi,h_Eff_Prompt_f0kkpi,BRphi,BRf0);

	TH1D *h_Eff_NonPrompt_AllBR=(TH1D*)h_Eff_NonPrompt_phikkpi->Clone("h_Eff_NonPrompt_AllBR");
	h_Eff_NonPrompt_AllBR->Sumw2();
	h_Eff_NonPrompt_AllBR->Add(h_Eff_NonPrompt_phikkpi,h_Eff_NonPrompt_f0kkpi,BRphi,BRf0);


	// read efficiency for cut scan

  TH1D *h_Eff_Prompt_phikkpi_Dchi2clMinScan[nbin_Dchi2clMinScan];
  TH1D *h_Eff_Prompt_f0kkpi_Dchi2clMinScan[nbin_Dchi2clMinScan];
  TH1D *h_Eff_NonPrompt_phikkpi_Dchi2clMinScan[nbin_Dchi2clMinScan];
  TH1D *h_Eff_NonPrompt_f0kkpi_Dchi2clMinScan[nbin_Dchi2clMinScan];
	TH1D *h_Eff_Prompt_AllBR_Dchi2clMinScan[nbin_Dchi2clMinScan];
	TH1D *h_Eff_NonPrompt_AllBR_Dchi2clMinScan[nbin_Dchi2clMinScan];
	str_effH="h_RecoNormEff_Dchi2cl";

	TH1F *h_RawFitYield_Dchi2clMinScan[nbin_Dchi2clMinScan];

	cout<<__LINE__<<endl;
	for(int i=0; i<nbin_Dchi2clMinScan; i++){
		h_Eff_Prompt_phikkpi_Dchi2clMinScan[i]=(TH1D*)f_Prompt_phikkpi->Get(Form("%s_%i",str_effH.Data(),i));	
		h_Eff_Prompt_f0kkpi_Dchi2clMinScan[i]=(TH1D*)f_Prompt_f0kkpi->Get(Form("%s_%i",str_effH.Data(),i));	
		h_Eff_NonPrompt_phikkpi_Dchi2clMinScan[i]=(TH1D*)f_NonPrompt_phikkpi->Get(Form("%s_%i",str_effH.Data(),i));	
		h_Eff_NonPrompt_f0kkpi_Dchi2clMinScan[i]=(TH1D*)f_NonPrompt_f0kkpi->Get(Form("%s_%i",str_effH.Data(),i));
		h_Eff_Prompt_AllBR_Dchi2clMinScan[i]=(TH1D*)h_Eff_Prompt_phikkpi_Dchi2clMinScan[i]->Clone(Form("h_Eff_Prompt_AllBR_Dchi2clMinScan_%i",i));
		h_Eff_Prompt_AllBR_Dchi2clMinScan[i]->Sumw2();
		h_Eff_Prompt_AllBR_Dchi2clMinScan[i]->Add( h_Eff_Prompt_phikkpi_Dchi2clMinScan[i],h_Eff_Prompt_f0kkpi_Dchi2clMinScan[i],BRphi,BRf0);

		h_Eff_NonPrompt_AllBR_Dchi2clMinScan[i]=(TH1D*)h_Eff_NonPrompt_phikkpi_Dchi2clMinScan[i]->Clone(Form("h_Eff_NonPrompt_AllBR_Dchi2clMinScan_%i",i));
		h_Eff_NonPrompt_AllBR_Dchi2clMinScan[i]->Sumw2();
		h_Eff_NonPrompt_AllBR_Dchi2clMinScan[i]->Add( h_Eff_NonPrompt_phikkpi_Dchi2clMinScan[i],h_Eff_NonPrompt_f0kkpi_Dchi2clMinScan[i],BRphi,BRf0);

		h_RawFitYield_Dchi2clMinScan[i]=(TH1F*)f_raw->Get(Form("h_RawRooFitYield_Dchi2clMinScan_%i",i));

	}

	cout<<__LINE__<<endl;

	// read efficiency for cut scan for Ddls
  TH1D *h_Eff_Prompt_phikkpi_DdlsMinScan[nbin_DdlsMinScan];
  TH1D *h_Eff_Prompt_f0kkpi_DdlsMinScan[nbin_DdlsMinScan];
  TH1D *h_Eff_NonPrompt_phikkpi_DdlsMinScan[nbin_DdlsMinScan];
  TH1D *h_Eff_NonPrompt_f0kkpi_DdlsMinScan[nbin_DdlsMinScan];
	TH1D *h_Eff_Prompt_AllBR_DdlsMinScan[nbin_DdlsMinScan];
	TH1D *h_Eff_NonPrompt_AllBR_DdlsMinScan[nbin_DdlsMinScan];
	str_effH="h_RecoNormEff_Ddls";
	TH1F *h_RawFitYield_DdlsMinScan[nbin_DdlsMinScan];

	for(int i=0; i<nbin_DdlsMinScan; i++){
		h_Eff_Prompt_phikkpi_DdlsMinScan[i]=(TH1D*)f_Prompt_phikkpi->Get(Form("%s_%i",str_effH.Data(),i));	
		h_Eff_Prompt_f0kkpi_DdlsMinScan[i]=(TH1D*)f_Prompt_f0kkpi->Get(Form("%s_%i",str_effH.Data(),i));	
		h_Eff_NonPrompt_phikkpi_DdlsMinScan[i]=(TH1D*)f_NonPrompt_phikkpi->Get(Form("%s_%i",str_effH.Data(),i));	
		h_Eff_NonPrompt_f0kkpi_DdlsMinScan[i]=(TH1D*)f_NonPrompt_f0kkpi->Get(Form("%s_%i",str_effH.Data(),i));
		h_Eff_Prompt_AllBR_DdlsMinScan[i]=(TH1D*)h_Eff_Prompt_phikkpi_DdlsMinScan[i]->Clone(Form("h_Eff_Prompt_AllBR_DdlsMinScan_%i",i));
		h_Eff_Prompt_AllBR_DdlsMinScan[i]->Sumw2();
		h_Eff_Prompt_AllBR_DdlsMinScan[i]->Add( h_Eff_Prompt_phikkpi_DdlsMinScan[i],h_Eff_Prompt_f0kkpi_DdlsMinScan[i],BRphi,BRf0);
		h_Eff_NonPrompt_AllBR_DdlsMinScan[i]=(TH1D*)h_Eff_NonPrompt_phikkpi_DdlsMinScan[i]->Clone(Form("h_Eff_NonPrompt_AllBR_DdlsMinScan_%i",i));
		h_Eff_NonPrompt_AllBR_DdlsMinScan[i]->Sumw2();
		h_Eff_NonPrompt_AllBR_DdlsMinScan[i]->Add( h_Eff_NonPrompt_phikkpi_DdlsMinScan[i],h_Eff_NonPrompt_f0kkpi_DdlsMinScan[i],BRphi,BRf0);
    h_RawFitYield_DdlsMinScan[i]=(TH1F*)f_raw->Get(Form("h_RawRooFitYield_DdlsMinScan_%i",i));

	}
	// end efficiency for cut scan for Ddls


	// read efficiency for cut scan for KaonPt
  TH1D *h_Eff_Prompt_phikkpi_KaonPtScan[nbin_DauPtScan];
  TH1D *h_Eff_Prompt_f0kkpi_KaonPtScan[nbin_DauPtScan];
  TH1D *h_Eff_NonPrompt_phikkpi_KaonPtScan[nbin_DauPtScan];
  TH1D *h_Eff_NonPrompt_f0kkpi_KaonPtScan[nbin_DauPtScan];
	TH1D *h_Eff_Prompt_AllBR_KaonPtScan[nbin_DauPtScan];
	TH1D *h_Eff_NonPrompt_AllBR_KaonPtScan[nbin_DauPtScan];
	str_effH="h_RecoNormEff_Ddls";
	TH1F *h_RawFitYield_KaonPtScan[nbin_DauPtScan];

  if(doDauPtCut){

	for(int i=0; i<nbin_DauPtScan; i++){
		h_Eff_Prompt_phikkpi_KaonPtScan[i]=(TH1D*)f_Prompt_phikkpi->Get(Form("%s_%i",str_effH.Data(),i));	
		h_Eff_Prompt_f0kkpi_KaonPtScan[i]=(TH1D*)f_Prompt_f0kkpi->Get(Form("%s_%i",str_effH.Data(),i));	
		h_Eff_NonPrompt_phikkpi_KaonPtScan[i]=(TH1D*)f_NonPrompt_phikkpi->Get(Form("%s_%i",str_effH.Data(),i));	
		h_Eff_NonPrompt_f0kkpi_KaonPtScan[i]=(TH1D*)f_NonPrompt_f0kkpi->Get(Form("%s_%i",str_effH.Data(),i));
		h_Eff_Prompt_AllBR_KaonPtScan[i]=(TH1D*)h_Eff_Prompt_phikkpi_KaonPtScan[i]->Clone(Form("h_Eff_Prompt_AllBR_KaonPtScan_%i",i));
		h_Eff_Prompt_AllBR_KaonPtScan[i]->Sumw2();
		h_Eff_Prompt_AllBR_KaonPtScan[i]->Add( h_Eff_Prompt_phikkpi_KaonPtScan[i],h_Eff_Prompt_f0kkpi_KaonPtScan[i],BRphi,BRf0);
		h_Eff_NonPrompt_AllBR_KaonPtScan[i]=(TH1D*)h_Eff_NonPrompt_phikkpi_KaonPtScan[i]->Clone(Form("h_Eff_NonPrompt_AllBR_KaonPtScan_%i",i));
		h_Eff_NonPrompt_AllBR_KaonPtScan[i]->Sumw2();
		h_Eff_NonPrompt_AllBR_KaonPtScan[i]->Add( h_Eff_NonPrompt_phikkpi_KaonPtScan[i],h_Eff_NonPrompt_f0kkpi_KaonPtScan[i],BRphi,BRf0);
    h_RawFitYield_KaonPtScan[i]=(TH1F*)f_raw->Get(Form("h_RawRooFitYield_KaonPtScan_%i",i));

	}
	// end efficiency for cut scan for KaonPt
	} // end doDauPtCut

	// read efficiency for cut scan for PionPt
  TH1D *h_Eff_Prompt_phikkpi_PionPtScan[nbin_DauPtScan];
  TH1D *h_Eff_Prompt_f0kkpi_PionPtScan[nbin_DauPtScan];
  TH1D *h_Eff_NonPrompt_phikkpi_PionPtScan[nbin_DauPtScan];
  TH1D *h_Eff_NonPrompt_f0kkpi_PionPtScan[nbin_DauPtScan];
	TH1D *h_Eff_Prompt_AllBR_PionPtScan[nbin_DauPtScan];
	TH1D *h_Eff_NonPrompt_AllBR_PionPtScan[nbin_DauPtScan];
	str_effH="h_RecoNormEff_Ddls";
	TH1F *h_RawFitYield_PionPtScan[nbin_DauPtScan];

  if(doDauPtCut){

	for(int i=0; i<nbin_DauPtScan; i++){
		h_Eff_Prompt_phikkpi_PionPtScan[i]=(TH1D*)f_Prompt_phikkpi->Get(Form("%s_%i",str_effH.Data(),i));	
		h_Eff_Prompt_f0kkpi_PionPtScan[i]=(TH1D*)f_Prompt_f0kkpi->Get(Form("%s_%i",str_effH.Data(),i));	
		h_Eff_NonPrompt_phikkpi_PionPtScan[i]=(TH1D*)f_NonPrompt_phikkpi->Get(Form("%s_%i",str_effH.Data(),i));	
		h_Eff_NonPrompt_f0kkpi_PionPtScan[i]=(TH1D*)f_NonPrompt_f0kkpi->Get(Form("%s_%i",str_effH.Data(),i));
		h_Eff_Prompt_AllBR_PionPtScan[i]=(TH1D*)h_Eff_Prompt_phikkpi_PionPtScan[i]->Clone(Form("h_Eff_Prompt_AllBR_PionPtScan_%i",i));
		h_Eff_Prompt_AllBR_PionPtScan[i]->Sumw2();
		h_Eff_Prompt_AllBR_PionPtScan[i]->Add( h_Eff_Prompt_phikkpi_PionPtScan[i],h_Eff_Prompt_f0kkpi_PionPtScan[i],BRphi,BRf0);
		h_Eff_NonPrompt_AllBR_PionPtScan[i]=(TH1D*)h_Eff_NonPrompt_phikkpi_PionPtScan[i]->Clone(Form("h_Eff_NonPrompt_AllBR_PionPtScan_%i",i));
		h_Eff_NonPrompt_AllBR_PionPtScan[i]->Sumw2();
		h_Eff_NonPrompt_AllBR_PionPtScan[i]->Add( h_Eff_NonPrompt_phikkpi_PionPtScan[i],h_Eff_NonPrompt_f0kkpi_PionPtScan[i],BRphi,BRf0);
    h_RawFitYield_PionPtScan[i]=(TH1F*)f_raw->Get(Form("h_RawRooFitYield_PionPtScan_%i",i));

	}
	// end efficiency for cut scan for PionPt
	}

	// read efficiency for cut scan for AllDauPt
  TH1D *h_Eff_Prompt_phikkpi_AllDauPtScan[nbin_DauPtScan];
  TH1D *h_Eff_Prompt_f0kkpi_AllDauPtScan[nbin_DauPtScan];
  TH1D *h_Eff_NonPrompt_phikkpi_AllDauPtScan[nbin_DauPtScan];
  TH1D *h_Eff_NonPrompt_f0kkpi_AllDauPtScan[nbin_DauPtScan];
	TH1D *h_Eff_Prompt_AllBR_AllDauPtScan[nbin_DauPtScan];
	TH1D *h_Eff_NonPrompt_AllBR_AllDauPtScan[nbin_DauPtScan];
	str_effH="h_RecoNormEff_Ddls";
	TH1F *h_RawFitYield_AllDauPtScan[nbin_DauPtScan];

	if(doDauPtCut){

	for(int i=0; i<nbin_DauPtScan; i++){
		h_Eff_Prompt_phikkpi_AllDauPtScan[i]=(TH1D*)f_Prompt_phikkpi->Get(Form("%s_%i",str_effH.Data(),i));	
		h_Eff_Prompt_f0kkpi_AllDauPtScan[i]=(TH1D*)f_Prompt_f0kkpi->Get(Form("%s_%i",str_effH.Data(),i));	
		h_Eff_NonPrompt_phikkpi_AllDauPtScan[i]=(TH1D*)f_NonPrompt_phikkpi->Get(Form("%s_%i",str_effH.Data(),i));	
		h_Eff_NonPrompt_f0kkpi_AllDauPtScan[i]=(TH1D*)f_NonPrompt_f0kkpi->Get(Form("%s_%i",str_effH.Data(),i));
		h_Eff_Prompt_AllBR_AllDauPtScan[i]=(TH1D*)h_Eff_Prompt_phikkpi_AllDauPtScan[i]->Clone(Form("h_Eff_Prompt_AllBR_AllDauPtScan_%i",i));
		h_Eff_Prompt_AllBR_AllDauPtScan[i]->Sumw2();
		h_Eff_Prompt_AllBR_AllDauPtScan[i]->Add( h_Eff_Prompt_phikkpi_AllDauPtScan[i],h_Eff_Prompt_f0kkpi_AllDauPtScan[i],BRphi,BRf0);
		h_Eff_NonPrompt_AllBR_AllDauPtScan[i]=(TH1D*)h_Eff_NonPrompt_phikkpi_AllDauPtScan[i]->Clone(Form("h_Eff_NonPrompt_AllBR_AllDauPtScan_%i",i));
		h_Eff_NonPrompt_AllBR_AllDauPtScan[i]->Sumw2();
		h_Eff_NonPrompt_AllBR_AllDauPtScan[i]->Add( h_Eff_NonPrompt_phikkpi_AllDauPtScan[i],h_Eff_NonPrompt_f0kkpi_AllDauPtScan[i],BRphi,BRf0);
    h_RawFitYield_AllDauPtScan[i]=(TH1F*)f_raw->Get(Form("h_RawRooFitYield_AllDauPtScan_%i",i));

	}
	// end efficiency for cut scan for AllDauPt
	}




	// read efficiency for cut scan Reschi2cl
  TH1D *h_Eff_Prompt_phikkpi_Reschi2clScan[nbin_Reschi2clScan];
  TH1D *h_Eff_Prompt_f0kkpi_Reschi2clScan[nbin_Reschi2clScan];
  TH1D *h_Eff_NonPrompt_phikkpi_Reschi2clScan[nbin_Reschi2clScan];
  TH1D *h_Eff_NonPrompt_f0kkpi_Reschi2clScan[nbin_Reschi2clScan];
	TH1D *h_Eff_Prompt_AllBR_Reschi2clScan[nbin_Reschi2clScan];
	TH1D *h_Eff_NonPrompt_AllBR_Reschi2clScan[nbin_Reschi2clScan];
	str_effH="h_RecoNormEff_Reschi2cl";

	TH1F *h_RawFitYield_Reschi2clScan[nbin_Reschi2clScan];

	cout<<__LINE__<<endl;
	for(int i=0; i<nbin_Reschi2clScan; i++){
		h_Eff_Prompt_phikkpi_Reschi2clScan[i]=(TH1D*)f_Prompt_phikkpi->Get(Form("%s_%i",str_effH.Data(),i));	
		h_Eff_Prompt_f0kkpi_Reschi2clScan[i]=(TH1D*)f_Prompt_f0kkpi->Get(Form("%s_%i",str_effH.Data(),i));	
		h_Eff_NonPrompt_phikkpi_Reschi2clScan[i]=(TH1D*)f_NonPrompt_phikkpi->Get(Form("%s_%i",str_effH.Data(),i));	
		h_Eff_NonPrompt_f0kkpi_Reschi2clScan[i]=(TH1D*)f_NonPrompt_f0kkpi->Get(Form("%s_%i",str_effH.Data(),i));
		h_Eff_Prompt_AllBR_Reschi2clScan[i]=(TH1D*)h_Eff_Prompt_phikkpi_Reschi2clScan[i]->Clone(Form("h_Eff_Prompt_AllBR_Reschi2clScan_%i",i));
		h_Eff_Prompt_AllBR_Reschi2clScan[i]->Sumw2();
		h_Eff_Prompt_AllBR_Reschi2clScan[i]->Add( h_Eff_Prompt_phikkpi_Reschi2clScan[i],h_Eff_Prompt_f0kkpi_Reschi2clScan[i],BRphi,BRf0);

		h_Eff_NonPrompt_AllBR_Reschi2clScan[i]=(TH1D*)h_Eff_NonPrompt_phikkpi_Reschi2clScan[i]->Clone(Form("h_Eff_NonPrompt_AllBR_Reschi2clScan_%i",i));
		h_Eff_NonPrompt_AllBR_Reschi2clScan[i]->Sumw2();
		h_Eff_NonPrompt_AllBR_Reschi2clScan[i]->Add( h_Eff_NonPrompt_phikkpi_Reschi2clScan[i],h_Eff_NonPrompt_f0kkpi_Reschi2clScan[i],BRphi,BRf0);

		h_RawFitYield_Reschi2clScan[i]=(TH1F*)f_raw->Get(Form("h_RawRooFitYield_Reschi2clScan_%i",i));

	}

	// end efficiency for cut scan Reschi2cl


  TString outfile="output/CutScanSys_FixShape_pp.root";
  TString infile="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/CrossSection_dNdpt/output/PromptDsCrossSectionPP_FixShape.root";

  double *DalphaMaxScan_bins=DalphaMaxScan_bins_pp;
  double *DdlsMinScan_bins=DdlsMinScan_bins_pp;
  double *Dchi2clMinScan_bins=Dchi2clMinScan_bins_pp;

  double *bins_DalphaMaxScan=bins_DalphaMaxScan_pp;
  double *bins_DdlsMinScan=bins_DdlsMinScan_pp;
  double *bins_Dchi2clMinScan=bins_Dchi2clMinScan_pp;

	double *DdlsMin_bins=DdlsMin_bins_pp;
	double *DalphaMax_bins=DalphaMax_bins_pp;
	double *Dchi2clMin_bins=Dchi2clMin_bins_pp;
	
	double *DtrkPtMin_bins=Dtrk1PtMin_bins_pp;


  if(isPbPb==3){
    str_PbPb="PbPb3";
		str_PbPbtext="PbPb";
    nbin_pt=nbin_pt_PbPb3;
    bins_pt=bins_pt_PbPb3;
    outfile="output/CutScanSys_FixShape_PbPb3.root";
    infile="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/CrossSection_dNdpt/output/PromptDsdNdptPbPb3_FixShape.root";

    DalphaMaxScan_bins=DalphaMaxScan_bins_PbPb3;
    DdlsMinScan_bins=DdlsMinScan_bins_PbPb3;
    Dchi2clMinScan_bins=Dchi2clMinScan_bins_PbPb3;

    bins_DalphaMaxScan=bins_DalphaMaxScan_PbPb3;
    bins_DdlsMinScan=bins_DdlsMinScan_PbPb3;
    bins_Dchi2clMinScan=bins_Dchi2clMinScan_PbPb3;

		DdlsMin_bins=DdlsMin_bins_PbPb3;
		DalphaMax_bins=DalphaMax_bins_PbPb3;
		Dchi2clMin_bins=Dchi2clMin_bins_PbPb3;

		DtrkPtMin_bins=Dtrk1PtMin_bins_PbPb3;

  }

	// for plot start from x=0, set new bin
	cout<<"nbin = "<<nbin_DalphaMaxScan<<endl;

	double bins_DdlsMinScan_new[nbin_DdlsMinScan+2];
	bins_DdlsMinScan_new[0]=0;
	for(int i=0; i<nbin_DdlsMinScan+1; i++){
		bins_DdlsMinScan_new[i+1]=bins_DdlsMinScan[i];
		cout<<"ibin = "<<i<<" , bins_DdlsMinScan_new = "<<bins_DdlsMinScan_new[i+1]<<endl;
	}




	cout<<"start marco"<<endl;
	TFile *fin=TFile::Open(infile.Data(),"READ");
	TFile *fout=TFile::Open(outfile.Data(),"RECREATE");
	fout->cd();

	// saving results
	TH1D *h_PromptDs_CutScanAll_SysRel=new TH1D("h_PromptDs_CutScanAll_SysRel","h_PromptDs_CutScanAll_SysRel",nbin_pt,bins_pt); h_PromptDs_CutScanAll_SysRel->Sumw2();
	TH1D *h_PromptDs_DalphaMaxScan_SysRel=new TH1D("h_PromptDs_DalphaMaxScan_SysRel","h_PromptDs_DalphaMaxScan_SysRel",nbin_pt,bins_pt); h_PromptDs_DalphaMaxScan_SysRel->Sumw2();
	TH1D *h_PromptDs_Dchi2clMinScan_SysRel=new TH1D("h_PromptDs_Dchi2clMinScan_SysRel","h_PromptDs_Dchi2clMinScan_SysRel",nbin_pt,bins_pt); h_PromptDs_Dchi2clMinScan_SysRel->Sumw2();
	TH1D *h_PromptDs_DdlsMinScan_SysRel=new TH1D("h_PromptDs_DdlsMinScan_SysRel","h_PromptDs_DdlsMinScan_SysRel",nbin_pt,bins_pt); h_PromptDs_DdlsMinScan_SysRel->Sumw2();

	TH1D *h_PromptDs_PhiMassScan_SysRel=new TH1D("h_PromptDs_PhiMassScan_SysRel","h_PromptDs_PhiMassScan_SysRel",nbin_pt,bins_pt); h_PromptDs_PhiMassScan_SysRel->Sumw2();
	TH1D *h_PromptDs_Reschi2clScan_SysRel=new TH1D("h_PromptDs_Reschi2clScan_SysRel","h_PromptDs_Reschi2clScan_SysRel",nbin_pt,bins_pt); h_PromptDs_Reschi2clScan_SysRel->Sumw2();


	// TH1D *h_PromptDs_CutScan_SysRel=new TH1D("h_PromptDs_CutScan_SysRel","h_PromptDs_CutScan_SysRel",nbin_pt,bins_pt); h_PromptDs_CutScan_SysRel->Sumw2();

	TLatex *tlt=new TLatex();
	// TF1 *f1_pol1=new TF1("f1_pol1","[0]+x*[1]");
	TF1 *f1_pol1=new TF1("f1_pol1","[0]+(x-[2])*[1]");
	f1_pol1->SetLineColor(2);
	f1_pol1->SetLineStyle(2);
	cout<<"reading file done "<<endl;

	TLine *tl=new TLine(0.065,1,0.22,1);
	tl->SetLineColor(2);
	tl->SetLineStyle(2);


	TLine *tl2=new TLine(0.00775,1,0.01125,1);
	tl2->SetLineColor(2);
	tl2->SetLineStyle(2);

	TH1D *h_PromptDs=(TH1D*)fin->Get(Form("h_PromptDs_%s",str_PbPb.Data())); // central value

	cout<<"reading one TH1D "<<endl;

	TH1D *h_PromptDs_DalphaMaxScan[nbin_pt];
	TH1D *h_PromptDs_DalphaMaxScan_ratio[nbin_pt];
	TH1D *h_PromptDs_DalphaMaxScan_Stat[nbin_pt];

	TGraphErrors *Gr_PromptDs_DalphaMaxScan_ratio[nbin_pt];


	TH1D *h_PromptDs_Dchi2clMinScan[nbin_pt];
	TH1D *h_PromptDs_Dchi2clMinScan_ratio[nbin_pt];
	TH1D *h_PromptDs_Dchi2clMinScan_Stat[nbin_pt];

	TH1D *h_PromptDs_Dchi2clMinScan_Yieldratio[nbin_pt];
	TH1D *h_PromptDs_Dchi2clMinScan_Effratio[nbin_pt];
	TH1D *h_NonPromptDs_Dchi2clMinScan_Effratio[nbin_pt];

	TH1D *h_PromptDs_Dchi2clMinScan_MCratio[nbin_pt];
	TH1D *h_NonPromptDs_Dchi2clMinScan_MCratio[nbin_pt];
	TH1D *h_All_Dchi2clMinScan_Yieldratio[nbin_pt];

	TGraphErrors *Gr_PromptDs_Dchi2clMinScan_ratio[nbin_pt];

	TH1D *h_PromptDs_DdlsMinScan[nbin_pt];
	TH1D *h_PromptDs_DdlsMinScan_ratio[nbin_pt];
	TH1D *h_PromptDs_DdlsMinScan_Stat[nbin_pt];

	TH1D *h_PromptDs_DdlsMinScan_Yieldratio[nbin_pt];
	TH1D *h_PromptDs_DdlsMinScan_Effratio[nbin_pt];
	TH1D *h_NonPromptDs_DdlsMinScan_Effratio[nbin_pt];
	TH1D *h_All_DdlsMinScan_Yieldratio[nbin_pt];

	TGraphErrors *Gr_PromptDs_DdlsMinScan_ratio[nbin_pt];

	TH1D *h_PromptDs_PhiMassScan[nbin_pt];
	TH1D *h_PromptDs_PhiMassScan_ratio[nbin_pt];
	TH1D *h_PromptDs_PhiMassScan_Stat[nbin_pt];

  TGraphErrors *Gr_PromptDs_PhiMassScan_ratio[nbin_pt];

	TH1D *h_PromptDs_Reschi2clScan[nbin_pt];
	TH1D *h_PromptDs_Reschi2clScan_ratio[nbin_pt];
	TH1D *h_PromptDs_Reschi2clScan_Stat[nbin_pt];

	TH1D *h_PromptDs_Reschi2clScan_Yieldratio[nbin_pt];
	TH1D *h_PromptDs_Reschi2clScan_Effratio[nbin_pt];
	TH1D *h_NonPromptDs_Reschi2clScan_Effratio[nbin_pt];

	TH1D *h_PromptDs_Reschi2clScan_MCratio[nbin_pt];
	TH1D *h_NonPromptDs_Reschi2clScan_MCratio[nbin_pt];
	TH1D *h_All_Reschi2clScan_Yieldratio[nbin_pt];


	TH1D *h_PromptDs_KaonPtScan[nbin_pt];
	TH1D *h_PromptDs_KaonPtScan_ratio[nbin_pt];
	TH1D *h_PromptDs_KaonPtScan_Stat[nbin_pt];

	TH1D *h_PromptDs_KaonPtScan_Yieldratio[nbin_pt];
	TH1D *h_PromptDs_KaonPtScan_Effratio[nbin_pt];
	TH1D *h_NonPromptDs_KaonPtScan_Effratio[nbin_pt];

	TH1D *h_PromptDs_KaonPtScan_MCratio[nbin_pt];
	TH1D *h_NonPromptDs_KaonPtScan_MCratio[nbin_pt];
	TH1D *h_All_KaonPtScan_Yieldratio[nbin_pt];

	TGraphErrors *Gr_PromptDs_KaonPtScan_ratio[nbin_pt];

	TH1D *h_PromptDs_PionPtScan[nbin_pt];
	TH1D *h_PromptDs_PionPtScan_ratio[nbin_pt];
	TH1D *h_PromptDs_PionPtScan_Stat[nbin_pt];

	TH1D *h_PromptDs_PionPtScan_Yieldratio[nbin_pt];
	TH1D *h_PromptDs_PionPtScan_Effratio[nbin_pt];
	TH1D *h_NonPromptDs_PionPtScan_Effratio[nbin_pt];

	TH1D *h_PromptDs_PionPtScan_MCratio[nbin_pt];
	TH1D *h_NonPromptDs_PionPtScan_MCratio[nbin_pt];
	TH1D *h_All_PionPtScan_Yieldratio[nbin_pt];

	TGraphErrors *Gr_PromptDs_PionPtScan_ratio[nbin_pt];

	TH1D *h_PromptDs_AllDauPtScan[nbin_pt];
	TH1D *h_PromptDs_AllDauPtScan_ratio[nbin_pt];
	TH1D *h_PromptDs_AllDauPtScan_Stat[nbin_pt];

	TH1D *h_PromptDs_AllDauPtScan_Yieldratio[nbin_pt];
	TH1D *h_PromptDs_AllDauPtScan_Effratio[nbin_pt];
	TH1D *h_NonPromptDs_AllDauPtScan_Effratio[nbin_pt];

	TH1D *h_PromptDs_AllDauPtScan_MCratio[nbin_pt];
	TH1D *h_NonPromptDs_AllDauPtScan_MCratio[nbin_pt];
	TH1D *h_All_AllDauPtScan_Yieldratio[nbin_pt];

	TGraphErrors *Gr_PromptDs_AllDauPtScan_ratio[nbin_pt];



	

	TLine *tl_Dalpha[nbin_pt];
	TLine *tl_Dchi2cl[nbin_pt];
	TLine *tl_Ddls[nbin_pt];
	TLine *tl_PhiMass[nbin_pt];
	TLine *tl_Reschi2cl[nbin_pt];

	TLine *tl_KaonPt[nbin_pt];
	TLine *tl_PionPt[nbin_pt];
	TLine *tl_AllDauPt[nbin_pt];



	for(int ibin_pt=startbin; ibin_pt<nbin_pt; ibin_pt++)
	// for(int ibin_pt=4; ibin_pt<5; ibin_pt++)
	{
		double DptLow=bins_pt[ibin_pt];
		double DptHigh=bins_pt[ibin_pt+1];

		double CSdNdpt=h_PromptDs->GetBinContent(ibin_pt+1);
		double CSdNdptErr=h_PromptDs->GetBinError(ibin_pt+1);
		double CSdNdptErrRel=CSdNdptErr/CSdNdpt;

		double cut_Dalpha=DalphaMax_bins[ibin_pt];
		double cut_Dchi2cl=Dchi2clMin_bins[ibin_pt];
		double cut_Ddls=DdlsMin_bins[ibin_pt];

		double cut_PhiMass=DtktkResmassCutWidth;
		double cut_Reschi2cl=Reschi2clCut;

		double cut_DauPt=DtrkPtMin_bins[ibin_pt];

		double Eff_Prompt_phi=h_Eff_Prompt_phikkpi->GetBinContent(ibin_pt+1);
		double Eff_Prompt_f0=h_Eff_Prompt_f0kkpi->GetBinContent(ibin_pt+1);
		double Eff_NonPrompt_phi=h_Eff_NonPrompt_phikkpi->GetBinContent(ibin_pt+1);
		double Eff_NonPrompt_f0=h_Eff_NonPrompt_f0kkpi->GetBinContent(ibin_pt+1);

		double Eff_Prompt_AllBR=h_Eff_Prompt_AllBR->GetBinContent(ibin_pt+1);
		double EffErr_Prompt_AllBR=h_Eff_Prompt_AllBR->GetBinError(ibin_pt+1);
		double EffErrRel_Prompt_AllBR=EffErr_Prompt_AllBR/Eff_Prompt_AllBR;


		double Eff_NonPrompt_AllBR=h_Eff_NonPrompt_AllBR->GetBinContent(ibin_pt+1);
		double EffErr_NonPrompt_AllBR=h_Eff_NonPrompt_AllBR->GetBinError(ibin_pt+1);
		double EffErrRel_NonPrompt_AllBR=EffErr_Prompt_AllBR/Eff_NonPrompt_AllBR;


		double Prompt_yield=CSdNdpt*2*LumiNevt*(BRphi*Eff_Prompt_phi+BRf0*Eff_Prompt_f0);
		double Prompt_yieldErr=Prompt_yield*CSdNdptErrRel;

	  double All_yield=h_RawFitYield->GetBinContent(ibin_pt+1);
		double All_yieldErr=h_RawFitYield->GetBinError(ibin_pt+1);


	//-- start DalphaMaxScan --//
	cout<<"start DalphaMaxScan Systematics "<<endl;

		h_PromptDs_DalphaMaxScan[ibin_pt]=(TH1D*)fin->Get(Form("h_PromptDs_DalphaMaxScan_pt%.0fto%.0f",DptLow,DptHigh));
		// h_PromptDs_DalphaMaxScan[ibin_pt]=(TH1D*)fin->Get(Form("h_PromptDs_DalphaMaxScan2_pt%.0fto%.0f",DptLow,DptHigh));
		cout<<"read DalphaMaxScan TH1D "<<endl;
		h_PromptDs_DalphaMaxScan_ratio[ibin_pt]=new TH1D(Form("h_PromptDs_DalphaMaxScan_ratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_DalphaMaxScan_ratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DalphaMaxScan, bins_DalphaMaxScan );
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]=new TH1D(Form("h_PromptDs_DalphaMaxScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_DalphaMaxScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),nbin_DalphaMaxScan, bins_DalphaMaxScan );

	// replace TH1D by TGraph 
		Gr_PromptDs_DalphaMaxScan_ratio[ibin_pt]=new TGraphErrors();



		cout<<"\n ibin_pt = "<<ibin_pt<<" , Ds pt = "<<DptLow<<" to "<<DptHigh<<endl;
		cout<<" cut_Dalpha = "<<cut_Dalpha<<" , cut_Dchi2cl = "<<cut_Dchi2cl<<" , cut_Ddls = "<<cut_Ddls<<endl;
		cout<<" CSdNdpt = "<<CSdNdpt<<" CSdNdptErr = "<<CSdNdptErr<<" , CSdNdptErrRel = "<<CSdNdptErrRel*100<<"%"<<endl;
		cout<<" \n DalphaMaxScan "<<endl;

		start_var=0;
		// if(isPbPb==3 &&ibin_pt<=2){ start_var=1;} // don't use first bin result, fitting failed
		
		double DalphaMaxScan_maxDiff=0;
		double DalphaMaxScan_maxDifferr=0;

		for(int ibin_var=start_var; ibin_var<nbin_DalphaMaxScan; ibin_var++){
			double CSdNdpt_var=h_PromptDs_DalphaMaxScan[ibin_pt]->GetBinContent(ibin_var+1);
			double CSdNdptErr_var=h_PromptDs_DalphaMaxScan[ibin_pt]->GetBinError(ibin_var+1);
			double CSdNdptErrRel_var=CSdNdptErr_var/CSdNdpt_var;	

			double cut_DalphaMaxScan=DalphaMaxScan_bins[ibin_var];

			double ratio=CSdNdpt_var/CSdNdpt;
			double ratioErrRel_temp = CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var - 2/CSdNdpt/CSdNdpt_var*CSdNdptErr*CSdNdptErr;
			// cout<<"ratioErrRel_temp = "<<ratioErrRel_temp<<endl;

/*
			if(cut_Dalpha > cut_DalphaMaxScan ){
			ratioErrRel_temp = CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var;
			}
			if(ratioErrRel_temp <0 ){
				cout<<"enter temp"<<endl;
			ratioErrRel_temp = CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var >0 ? CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var : CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var - 2/CSdNdpt/CSdNdpt_var*CSdNdptErr*CSdNdptErr;
			}
			if(ratioErrRel_temp <0 ){
				cout<<" --- warning , ratioErrRel_temp still < 0 !!! "<<endl;
			}
			cout<<"ratioErrRel_temp = "<<ratioErrRel_temp<<endl;
*/
			//double ratioErr = ratio*sqrt(ratioErrRel_temp);
			double ratioErr = sqrt(abs(CSdNdptErr_var*CSdNdptErr_var-CSdNdptErr*CSdNdptErr))/CSdNdpt;
			
			// special case for var_cut = default cut
			if(cut_Dalpha == cut_DalphaMaxScan){
				ratio=1;
				ratioErr=0;
			}





			// double ratioErr=ratio*sqrt( CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr*CSdNdptErr  );
//			if(cut_Dalpha < cut_DalphaMaxScan ){
	//			ratioErr=ratio*sqrt( CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var  );			
//			}

			cout<<"\n ibin_var = "<<ibin_var<<" , DalphaMaxScan var = "<<cut_DalphaMaxScan<<endl;
			cout<<"CSdNdpt_var = "<<CSdNdpt_var<<" , CSdNdptErr_var = "<<CSdNdptErr_var<<" , CSdNdptErrRel_var = "<<CSdNdptErrRel_var*100<<"%"<<endl;
			cout<<"ratio = "<<ratio<<" +- "<<ratioErr<<" relErr = "<<ratioErr/ratio*100<<"%"<<endl;

			h_PromptDs_DalphaMaxScan_ratio[ibin_pt]->SetBinContent(ibin_var+1, ratio);	
			h_PromptDs_DalphaMaxScan_ratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr);	

			h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetBinContent(ibin_var+1, 1);	
			h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetBinError(ibin_var+1, CSdNdptErrRel);	

			if(cut_Dalpha != cut_DalphaMaxScan){
			Gr_PromptDs_DalphaMaxScan_ratio[ibin_pt]->SetPoint(ibin_var,DalphaMaxScan_bins[ibin_var],ratio);
			Gr_PromptDs_DalphaMaxScan_ratio[ibin_pt]->SetPointError(ibin_var,0,ratioErr);
			}

			if(abs(ratio-1) >DalphaMaxScan_maxDiff){
				DalphaMaxScan_maxDiff = abs(ratio-1); 
				DalphaMaxScan_maxDifferr= ratioErr; 
			}

		} // end for ibin_var<nbin_DalphaMaxScan

		h_PromptDs_DalphaMaxScan_SysRel->SetBinContent(ibin_pt+1, DalphaMaxScan_maxDiff);

		fout->cd();
		h_PromptDs_DalphaMaxScan_ratio[ibin_pt]->Write("",TObject::kOverwrite);

		Gr_PromptDs_DalphaMaxScan_ratio[ibin_pt]->Write("",TObject::kOverwrite);
		
		// fit		

		// draw

		cav[c_count]= new TCanvas(Form("cav%i",c_count),Form("cav%i",c_count),800,800);
		cav[c_count]->cd();
//		Gr_PromptDs_DalphaMaxScan_ratio[ibin_pt]->SetMinimum(0);
	//	Gr_PromptDs_DalphaMaxScan_ratio[ibin_pt]->SetMaximum(2.5);
		//Gr_PromptDs_DalphaMaxScan_ratio[ibin_pt]->SetMarkerSize(0.25);
//		Gr_PromptDs_DalphaMaxScan_ratio[ibin_pt]->SetMarkerStyle(24);
	//	Gr_PromptDs_DalphaMaxScan_ratio[ibin_pt]->Draw("AP");
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetMinimum(RatioMin);
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->GetXaxis()->SetTitle("D_{S} pointing angle");
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->GetYaxis()->SetTitle("#sigma_{Varied cut}/#sigma_{Default cut}");
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetTitle("");
	
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetFillColor(kGray);
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetFillStyle(3001);
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetMarkerColor(kGray);
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetMarkerSize(0);
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetLineWidth(0);
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetLineColor(kGray);
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->Draw("SAME E2");

		Gr_PromptDs_DalphaMaxScan_ratio[ibin_pt]->SetMinimum(0);
		Gr_PromptDs_DalphaMaxScan_ratio[ibin_pt]->SetMaximum(2.5);
		Gr_PromptDs_DalphaMaxScan_ratio[ibin_pt]->SetMarkerSize(0.5);
		Gr_PromptDs_DalphaMaxScan_ratio[ibin_pt]->SetMarkerStyle(21);
		Gr_PromptDs_DalphaMaxScan_ratio[ibin_pt]->Draw("P");
	
		tl->Draw("same");

		shiftY=0;
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;
		if(isPbPb==3){
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,"centrality : 0-100%"); shiftY-=oneshift;
		}
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%% #pm %.1f%%",DalphaMaxScan_maxDiff*100, DalphaMaxScan_maxDifferr*100)); shiftY-=oneshift;

		tl_Dalpha[ibin_pt]=new TLine(cut_Dalpha,RatioMin,cut_Dalpha,RatioMax);
		tl_Dalpha[ibin_pt]->SetLineColor(4);
		tl_Dalpha[ibin_pt]->SetLineStyle(6);
		tl_Dalpha[ibin_pt]->Draw("same");
	


		// SavePlotDirs(cav[c_count],Form("%s_Sys_DalphaMaxScan_Graph_pt%.0fto%.0f",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );

		c_count++;




		cav[c_count]= new TCanvas(Form("cav%i",c_count),Form("cav%i",c_count),800,800);
		cav[c_count]->cd();

		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetMinimum(RatioMin);
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->GetXaxis()->SetTitle("D_{S} pointing angle");
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->GetYaxis()->SetTitle("#sigma_{Varied cut}/#sigma_{Default cut}");
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetTitle("");
	
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetFillColor(kGray);
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetFillStyle(3001);
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetMarkerColor(kGray);
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetMarkerSize(0);
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetLineWidth(0);
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->SetLineColor(kGray);
		h_PromptDs_DalphaMaxScan_Stat[ibin_pt]->Draw("SAME E2");

		h_PromptDs_DalphaMaxScan_ratio[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_DalphaMaxScan_ratio[ibin_pt]->SetMinimum(RatioMin);
		h_PromptDs_DalphaMaxScan_ratio[ibin_pt]->GetXaxis()->SetTitle("D_{S} pointing angle");
		h_PromptDs_DalphaMaxScan_ratio[ibin_pt]->SetTitle("");
//		h_PromptDs_DalphaMaxScan_ratio[ibin_pt]->GetYaxis()->SetTitle("");
		h_PromptDs_DalphaMaxScan_ratio[ibin_pt]->Draw("SAME");
		tl->Draw("same");


		shiftY=0;
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;
		if(isPbPb==3){
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,"centrality : 0-100%"); shiftY-=oneshift;
		}
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%% #pm %.1f%%",DalphaMaxScan_maxDiff*100, DalphaMaxScan_maxDifferr*100)); shiftY-=oneshift;

		tl_Dalpha[ibin_pt]=new TLine(cut_Dalpha,RatioMin,cut_Dalpha,RatioMax);
		tl_Dalpha[ibin_pt]->SetLineColor(4);
		tl_Dalpha[ibin_pt]->SetLineStyle(6);
		tl_Dalpha[ibin_pt]->Draw("same");
		

		SavePlotDirs(cav[c_count],Form("%s_Sys_DalphaMaxScan_pt%.0fto%.0f",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );

		c_count++;

		// special handle


	cout<<"end DalphaMaxScan Systematics "<<endl;
		//-- end DalphaMaxScan_bins --//



	//-- start Dchi2clMinScan --//
	cout<<"start Dchi2clMinScan Systematics "<<endl;

		h_PromptDs_Dchi2clMinScan[ibin_pt]=(TH1D*)fin->Get(Form("h_PromptDs_Dchi2clMinScan_pt%.0fto%.0f",DptLow,DptHigh));
		// h_PromptDs_Dchi2clMinScan[ibin_pt]=(TH1D*)fin->Get(Form("h_PromptDs_Dchi2clMinScan2_pt%.0fto%.0f",DptLow,DptHigh));
		cout<<"read Dchi2clMinScan TH1D "<<endl;
		h_PromptDs_Dchi2clMinScan_ratio[ibin_pt]=new TH1D(Form("h_PromptDs_Dchi2clMinScan_ratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_Dchi2clMinScan_ratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_Dchi2clMinScan, bins_Dchi2clMinScan );

		h_PromptDs_Dchi2clMinScan_Yieldratio[ibin_pt]=new TH1D(Form("h_PromptDs_Dchi2clMinScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_Dchi2clMinScan_Yiledratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_Dchi2clMinScan, bins_Dchi2clMinScan );
		h_PromptDs_Dchi2clMinScan_Effratio[ibin_pt]=new TH1D(Form("h_PromptDs_Dchi2clMinScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_Dchi2clMinScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_Dchi2clMinScan, bins_Dchi2clMinScan );
		h_NonPromptDs_Dchi2clMinScan_Effratio[ibin_pt]=new TH1D(Form("h_NonPromptDs_Dchi2clMinScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_NonPromptDs_Dchi2clMinScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_Dchi2clMinScan, bins_Dchi2clMinScan );
    Gr_PromptDs_Dchi2clMinScan_ratio[ibin_pt]=new TGraphErrors();

		h_PromptDs_Dchi2clMinScan_MCratio[ibin_pt]=new TH1D(Form("h_PromptDs_Dchi2clMinScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_Dchi2clMinScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_Dchi2clMinScan, bins_Dchi2clMinScan );
		h_NonPromptDs_Dchi2clMinScan_MCratio[ibin_pt]=new TH1D(Form("h_NonPromptDs_Dchi2clMinScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_NonPromptDs_Dchi2clMinScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_Dchi2clMinScan, bins_Dchi2clMinScan );
		h_All_Dchi2clMinScan_Yieldratio[ibin_pt]=new TH1D(Form("h_All_Dchi2clMinScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_All_Dchi2clMinScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_Dchi2clMinScan, bins_Dchi2clMinScan );

		// h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]=new TH1D(Form("h_PromptDs_Dchi2clMinScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_Dchi2clMinScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),nbin_Dchi2clMinScan, bins_Dchi2clMinScan );
		h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]=new TH1D(Form("h_PromptDs_Dchi2clMinScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_Dchi2clMinScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),1,0, bins_Dchi2clMinScan[nbin_Dchi2clMinScan] );

		cout<<"\n ibin_pt = "<<ibin_pt<<" , Ds pt = "<<DptLow<<" to "<<DptHigh<<endl;
		cout<<" cut_Dchi2cl = "<<cut_Dchi2cl<<" , cut_Dchi2cl = "<<cut_Dchi2cl<<" , cut_Ddls = "<<cut_Ddls<<endl;
		cout<<" CSdNdpt = "<<CSdNdpt<<" CSdNdptErr = "<<CSdNdptErr<<" , CSdNdptErrRel = "<<CSdNdptErrRel*100<<"%"<<endl;
		cout<<" \n Dchi2clMinScan "<<endl;

		start_var=0;
		if(isPbPb==3 &&ibin_pt<=2){ start_var=1;} // don't use first bin result, fitting failed


		// TH1D *h_eff_Dchi2clMinScan[nbin_Dchi2clMinScan];
		

		for(int ibin_var=start_var; ibin_var<nbin_Dchi2clMinScan; ibin_var++){

			double CSdNdpt_var=h_PromptDs_Dchi2clMinScan[ibin_pt]->GetBinContent(ibin_var+1);
			double CSdNdptErr_var=h_PromptDs_Dchi2clMinScan[ibin_pt]->GetBinError(ibin_var+1);
			double CSdNdptErrRel_var=CSdNdptErr_var/CSdNdpt_var;	

			double cut_Dchi2clMinScan=Dchi2clMinScan_bins[ibin_var];

			double ratio=CSdNdpt_var/CSdNdpt;
			double ratioErrRel_temp = CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var - 2/CSdNdpt/CSdNdpt_var*CSdNdptErr*CSdNdptErr;
			// cout<<"ratioErrRel_temp = "<<ratioErrRel_temp<<endl;
	// for single ratio 

		double Eff_Prompt_phi_var=h_Eff_Prompt_phikkpi_Dchi2clMinScan[ibin_var]->GetBinContent(ibin_pt+1);
		double Eff_Prompt_f0_var=h_Eff_Prompt_f0kkpi_Dchi2clMinScan[ibin_var]->GetBinContent(ibin_pt+1);
		double Eff_NonPrompt_phi_var=h_Eff_NonPrompt_phikkpi_Dchi2clMinScan[ibin_var]->GetBinContent(ibin_pt+1);
		double Eff_NonPrompt_f0_var=h_Eff_NonPrompt_f0kkpi_Dchi2clMinScan[ibin_var]->GetBinContent(ibin_pt+1);
	
		double Eff_Prompt_AllBR_var=h_Eff_Prompt_AllBR_Dchi2clMinScan[ibin_var]->GetBinContent(ibin_pt+1);
		double EffErr_Prompt_AllBR_var=h_Eff_Prompt_AllBR_Dchi2clMinScan[ibin_var]->GetBinError(ibin_pt+1);
		double EffErrRel_Prompt_AllBR_var=EffErr_Prompt_AllBR_var/Eff_Prompt_AllBR_var;

		double Eff_NonPrompt_AllBR_var=h_Eff_NonPrompt_AllBR_Dchi2clMinScan[ibin_var]->GetBinContent(ibin_pt+1);
		double EffErr_NonPrompt_AllBR_var=h_Eff_NonPrompt_AllBR_Dchi2clMinScan[ibin_var]->GetBinError(ibin_pt+1);
		double EffErrRel_NonPrompt_AllBR_var=EffErr_NonPrompt_AllBR_var/Eff_NonPrompt_AllBR_var;


		double Prompt_yield_var=CSdNdpt_var*2*LumiNevt*(BRphi*Eff_Prompt_phi_var+BRf0*Eff_Prompt_f0_var);
		double Prompt_yieldErr_var=Prompt_yield_var*CSdNdptErrRel_var;

		double ratio_yield=Prompt_yield_var/Prompt_yield;

		double ratio_Eff=Eff_Prompt_AllBR_var/Eff_Prompt_AllBR;
		double ratioErr_Eff=ratio_Eff*( sqrt(EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR+ EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var)  );


		double ratio_NPEff=Eff_NonPrompt_AllBR_var/Eff_NonPrompt_AllBR;
		double ratioErr_NPEff=ratio_NPEff*( sqrt(EffErrRel_NonPrompt_AllBR*EffErrRel_NonPrompt_AllBR+ EffErrRel_NonPrompt_AllBR_var*EffErrRel_NonPrompt_AllBR_var)  );


		double All_yield_var=h_RawFitYield_Dchi2clMinScan[ibin_var]->GetBinContent(ibin_pt+1);
		double ratio_All_yield=All_yield_var/All_yield;
		double ratio_All_yieldErr=0;

		// cout<<"Eff_Prompt_AllBR_var = "<<Eff_Prompt_AllBR_var<<" , Eff_Prompt_AllBR = "<<Eff_Prompt_AllBR<<" , ratio_Eff = "<<ratio_Eff<<endl;
	// return 1;
		double ratioErrRel_temp_Eff=EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR + EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var-2/Eff_Prompt_AllBR/EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR*EffErr_Prompt_AllBR;



//		double ratioErr_yield=ratio_yield*


/*
 // cout test
		cout<<"ibin_pt = "<<ibin_pt<<" , ibin_var = "<<ibin_var<<endl;
		cout<<"Prompt_yield = "<<Prompt_yield<<" , CSdNdpt = "<<CSdNdpt<<" , Eff_Prompt_phi = "<<Eff_Prompt_phi<<", Eff_Prompt_f0 = "<<Eff_Prompt_f0<<", Eff_Prompt_AllBR = "<< Eff_Prompt_AllBR<<endl;
		cout<<"Prompt_yield_var = "<<Prompt_yield_var<<" , CSdNdpt_var = "<<CSdNdpt_var<<" , Eff_Prompt_phi_var = "<<Eff_Prompt_phi_var<<", Eff_Prompt_f0_var = "<<Eff_Prompt_f0_var<<", Eff_Prompt_AllBR_var = "<< Eff_Prompt_AllBR_var<<endl;

		return 1;
*/
			
	// end single ratio

			if(cut_Dchi2cl > cut_Dchi2clMinScan ){
			ratioErrRel_temp = CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var;
			ratioErrRel_temp_Eff=EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR + EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var-2/Eff_Prompt_AllBR/EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR_var;
			}
			if(ratioErrRel_temp <0 ){
			ratioErrRel_temp = CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var >0 ? CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var : CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var - 2/CSdNdpt/CSdNdpt_var*CSdNdptErr*CSdNdptErr;
      ratioErrRel_temp_Eff= EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR + EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var-2/Eff_Prompt_AllBR/EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR_var >0 ? EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR + EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var-2/Eff_Prompt_AllBR/EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR_var : EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR + EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var-2/Eff_Prompt_AllBR/EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR*EffErr_Prompt_AllBR;

			}
			if(ratioErrRel_temp <0 ){
				cout<<" --- warning , ratioErrRel_temp still < 0 !!! "<<endl;
			}
			if(ratioErrRel_temp_Eff <0 ){
				// cout<<"temp 1 = "<<EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR + EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var-2/Eff_Prompt_AllBR/EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR_var<<endl;
				// cout<<"temp 2 = "<<EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR + EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var-2/Eff_Prompt_AllBR/EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR*EffErr_Prompt_AllBR<<endl;
			}


			cout<<"ratioErrRel_temp = "<<ratioErrRel_temp<<endl;
			double ratioErr = ratio*sqrt(ratioErrRel_temp);
			// double ratioErr_Eff=ratio_Eff*sqrt(ratioErrRel_temp_Eff);
			ratioErr = sqrt(abs(CSdNdptErr_var*CSdNdptErr_var-CSdNdptErr*CSdNdptErr))/CSdNdpt;	
			double ratioErr_yield=ratio_yield*sqrt(ratioErrRel_temp);
			ratio_All_yieldErr=ratio_All_yield*sqrt(ratioErrRel_temp);


	//// sigle ratio use full error //
			ratioErr_yield=Prompt_yieldErr_var/Prompt_yield;
			ratio_All_yieldErr=h_RawFitYield_Dchi2clMinScan[ibin_var]->GetBinError(ibin_pt+1)/All_yield;
			ratioErr_Eff=EffErr_Prompt_AllBR_var/Eff_Prompt_AllBR;
			ratioErr_NPEff=EffErr_NonPrompt_AllBR_var/Eff_NonPrompt_AllBR;

			cout<<"ratioErr_Eff = "<<ratioErr_Eff<<endl;

			// special case for var_cut = default cut
			if(cut_Dchi2cl == cut_Dchi2clMinScan){
				ratio=1;
				ratioErr=0;
			}

			// double ratioErr=ratio*sqrt( CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr*CSdNdptErr  );
//			if(cut_Dchi2cl < cut_Dchi2clMinScan ){
	//			ratioErr=ratio*sqrt( CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var  );			
//			}

			cout<<"\n ibin_var = "<<ibin_var<<" , Dchi2clMinScan var = "<<cut_Dchi2clMinScan<<endl;
			cout<<"CSdNdpt_var = "<<CSdNdpt_var<<" , CSdNdptErr_var = "<<CSdNdptErr_var<<" , CSdNdptErrRel_var = "<<CSdNdptErrRel_var*100<<"%"<<endl;
			cout<<"ratio = "<<ratio<<" +- "<<ratioErr<<" relErr = "<<ratioErr/ratio*100<<"%"<<endl;

			h_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->SetBinContent(ibin_var+1, ratio);	
			h_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr);	

			h_PromptDs_Dchi2clMinScan_Yieldratio[ibin_pt]->SetBinContent(ibin_var+1 , ratio_yield);
			h_PromptDs_Dchi2clMinScan_Yieldratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr_yield);
			h_PromptDs_Dchi2clMinScan_Effratio[ibin_pt]->SetBinContent(ibin_var+1, ratio_Eff);
			h_PromptDs_Dchi2clMinScan_Effratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr_Eff);

			h_NonPromptDs_Dchi2clMinScan_Effratio[ibin_pt]->SetBinContent(ibin_var+1, ratio_NPEff);
			h_NonPromptDs_Dchi2clMinScan_Effratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr_NPEff);
	
			h_All_Dchi2clMinScan_Yieldratio[ibin_pt]->SetBinContent(ibin_var+1, ratio_All_yield);
			h_All_Dchi2clMinScan_Yieldratio[ibin_pt]->SetBinError(ibin_var+1, ratio_All_yieldErr);

      if(cut_Dchi2cl != cut_Dchi2clMinScan){
      Gr_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->SetPoint(ibin_var,Dchi2clMinScan_bins[ibin_var],ratio);
      Gr_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->SetPointError(ibin_var,0,ratioErr);
      }


	
			// h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetBinContent(ibin_var+1, 1);	
			// h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetBinError(ibin_var+1, CSdNdptErrRel);	

		} // end for ibin_var<nbin_Dchi2clMinScan


			h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetBinContent(1, 1);	
			h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetBinError(1, CSdNdptErrRel);	

		fout->cd();
		h_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->Write("",TObject::kOverwrite);
		h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->Write("",TObject::kOverwrite);
		Gr_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->Write("",TObject::kOverwrite);
		
		// fit			
		// TF1 *f1_pol1=new TF1("f1_pol1","[0]+x*[1]");
		f1_pol1->SetRange(0,0.6);
		f1_pol1->SetParameter(0,1);
		// f1_pol1->FixParameter(0,1);
		f1_pol1->SetParameter(1,0);
		f1_pol1->SetParameter(2,cut_Dchi2cl);
		f1_pol1->FixParameter(2,cut_Dchi2cl);

		h_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->Fit("f1_pol1","QN0 W");
		h_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->Fit("f1_pol1","EMS0 ");
		h_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->Fit("f1_pol1","EMISN0");

    fitR=h_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->Fit("f1_pol1","EMIS0");
    // fitR=h_PromptDs_DdlsMinScan_ratio[ibin_pt]->Fit("f1_pol1","LEMS0","",2,8);
    // fitR=h_PromptDs_DdlsMinScan_ratio[ibin_pt]->Fit("f1_pol1","LEMS0I","",2,8);
    // fitR=h_PromptDs_DdlsMinScan_ratio[ibin_pt]->Fit("f1_pol1","EMS0","",2,8);
    fitR->GetConfidenceIntervals(1,1,1,x0,x0err,0.683,false);


    cav[c_count]= new TCanvas(Form("cav%i",c_count),Form("cav%i",c_count),800,800);
    cav[c_count]->cd();
    h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetMaximum(RatioMax);
    h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetMinimum(RatioMin);
    h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->GetXaxis()->SetTitle("D_{S} vertex probability");
    h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->GetYaxis()->SetTitle("#sigma_{Varied cut}/#sigma_{Default cut}");
    h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetTitle("");

    h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetFillColor(kGray);
    h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetFillStyle(3001);
    h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetMarkerColor(kGray);
    h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetMarkerSize(0);
    h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetLineWidth(0);
    h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetLineColor(kGray);
    h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->Draw("SAME E2");
    Gr_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->SetMarkerSize(0.5);
    Gr_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->SetMarkerStyle(21);
    Gr_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->Draw("P");

    f1_pol1->SetRange(0,0.6);
    f1_pol1->Draw("same");

    shiftY=0;
    tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;
    if(isPbPb==3){
    tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,"centrality : 0-100%"); shiftY-=oneshift;
    }
    tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%%",(abs(f1_pol1->GetParameter(0)-1))*100)); shiftY-=oneshift;

    tl_Dchi2cl[ibin_pt]=new TLine(cut_Dchi2cl,RatioMin,cut_Dchi2cl,RatioMax);
    tl_Dchi2cl[ibin_pt]->SetLineColor(4);
    tl_Dchi2cl[ibin_pt]->SetLineStyle(6);
    tl_Dchi2cl[ibin_pt]->Draw("same");


    // SavePlotDirs(cav[c_count],Form("%s_Sys_Dchi2clMinScan_Graph_pt%.0fto%.0f",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );

    c_count++;



		cav1[c_count]= new TCanvas(Form("cav1_%i",c_count),Form("cav1_%i",c_count),800,800);
		cav1[c_count]->cd();
		h_PromptDs_Dchi2clMinScan_Yieldratio[ibin_pt]->GetXaxis()->SetRangeUser(0,0.5);
		h_PromptDs_Dchi2clMinScan_Yieldratio[ibin_pt]->GetYaxis()->SetTitle("N_{vary cuts}/N_{default cuts}");
		h_PromptDs_Dchi2clMinScan_Yieldratio[ibin_pt]->GetXaxis()->SetTitle("D_{S} vertex probability");
		h_PromptDs_Dchi2clMinScan_Yieldratio[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_Dchi2clMinScan_Yieldratio[ibin_pt]->SetMinimum(0.0);
		h_PromptDs_Dchi2clMinScan_Yieldratio[ibin_pt]->SetLineColor(2);
		h_PromptDs_Dchi2clMinScan_Yieldratio[ibin_pt]->SetMarkerColor(2);
		h_PromptDs_Dchi2clMinScan_Yieldratio[ibin_pt]->Draw("SAME");
		h_PromptDs_Dchi2clMinScan_Effratio[ibin_pt]->SetLineColor(4);
		h_PromptDs_Dchi2clMinScan_Effratio[ibin_pt]->SetMarkerColor(4);
		h_PromptDs_Dchi2clMinScan_Effratio[ibin_pt]->Draw("SAME");
		h_All_Dchi2clMinScan_Yieldratio[ibin_pt]->SetLineColor(1);
		h_All_Dchi2clMinScan_Yieldratio[ibin_pt]->SetMarkerColor(1);
		h_All_Dchi2clMinScan_Yieldratio[ibin_pt]->Draw("SAME");
		h_NonPromptDs_Dchi2clMinScan_Effratio[ibin_pt]->SetLineColor(kMagenta);
		h_NonPromptDs_Dchi2clMinScan_Effratio[ibin_pt]->SetMarkerColor(kMagenta);
		h_NonPromptDs_Dchi2clMinScan_Effratio[ibin_pt]->Draw("SAME");


		TLegend *le_Dchi2clMinScan_SingleRatio=new TLegend(0.7,0.7,0.85,0.85);
		le_Dchi2clMinScan_SingleRatio->SetBorderSize(0);	
		le_Dchi2clMinScan_SingleRatio->AddEntry(h_All_Dchi2clMinScan_Yieldratio[ibin_pt],"Data","lp");
		le_Dchi2clMinScan_SingleRatio->AddEntry(h_PromptDs_Dchi2clMinScan_Yieldratio[ibin_pt],"Data Prompt","lp");
		le_Dchi2clMinScan_SingleRatio->AddEntry(h_PromptDs_Dchi2clMinScan_Effratio[ibin_pt],"Prompt MC","lp");
		le_Dchi2clMinScan_SingleRatio->AddEntry(h_NonPromptDs_Dchi2clMinScan_Effratio[ibin_pt],"NonPrompt MC","lp");
		le_Dchi2clMinScan_SingleRatio->Draw("same");

		shiftY=0;
	  tlt->DrawLatexNDC(textposx+0.0,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;

		TLine *tl_Dchi2cl_SigleRatio=new TLine(cut_Dchi2cl,RatioMin,cut_Dchi2cl,RatioMax);
		tl_Dchi2cl_SigleRatio->SetLineColor(4);
		tl_Dchi2cl_SigleRatio->SetLineStyle(6);
		tl_Dchi2cl_SigleRatio->Draw("same");
	

		SavePlotDirs(cav1[c_count],Form("%s_Sys_Dchi2clMinScan_pt%.0fto%.0f_sigleRatio",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );


		// draw
		cav[c_count]= new TCanvas(Form("cav%i",c_count),Form("cav%i",c_count),800,800);
		cav[c_count]->cd();

		h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetMinimum(RatioMin);
		h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->GetXaxis()->SetTitle("D_{S} vertex probability");
		h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->GetYaxis()->SetTitle("#sigma_{Varied cut}/#sigma_{Default cut}");
		h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetTitle("");
	
		h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetFillColor(kGray);
		h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetFillStyle(3001);
		h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetMarkerColor(kGray);
		h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetMarkerSize(0);
		h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetLineWidth(0);
		h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->SetLineColor(kGray);
		h_PromptDs_Dchi2clMinScan_Stat[ibin_pt]->Draw("SAME E2");


		h_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->GetXaxis()->SetRangeUser(0,0.6);
		h_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->SetMinimum(RatioMin);
		h_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->GetXaxis()->SetTitle("D_{S} vertex probability");
		h_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->SetTitle("");
//		h_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->GetYaxis()->SetTitle("");
		h_PromptDs_Dchi2clMinScan_ratio[ibin_pt]->Draw("SAME");
		f1_pol1->SetRange(0,0.6);
		f1_pol1->Draw("same");

		shiftY=0;
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;
		if(isPbPb==3){
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,"centrality : 0-100%"); shiftY-=oneshift;
		}
	  // tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%%",(abs(f1_pol1->GetParameter(0)-1))*100)); shiftY-=oneshift;
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%% #pm %.1f%%",(abs(f1_pol1->Eval(0)-1))*100, x0err[0]*100 )); shiftY-=oneshift;

		tl_Dchi2cl[ibin_pt]=new TLine(cut_Dchi2cl,RatioMin,cut_Dchi2cl,RatioMax);
		tl_Dchi2cl[ibin_pt]->SetLineColor(4);
		tl_Dchi2cl[ibin_pt]->SetLineStyle(6);
		tl_Dchi2cl[ibin_pt]->Draw("same");
	

		SavePlotDirs(cav[c_count],Form("%s_Sys_Dchi2clMinScan_pt%.0fto%.0f",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );

		c_count++;

		// special handle

		// h_PromptDs_Dchi2clMinScan_SysRel->SetBinContent(ibin_pt+1, abs(f1_pol1->GetParameter(0)-1));
		h_PromptDs_Dchi2clMinScan_SysRel->SetBinContent(ibin_pt+1, abs(f1_pol1->Eval(0)-1));

	cout<<"end Dchi2clMinScan Systematics "<<endl;
		//-- end Dchi2clMinScan_bins --//

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

	//-- start DdlsMinScan --//
	cout<<"start DdlsMinScan Systematics "<<endl;

    Gr_PromptDs_DdlsMinScan_ratio[ibin_pt]=new TGraphErrors();

		h_PromptDs_DdlsMinScan[ibin_pt]=(TH1D*)fin->Get(Form("h_PromptDs_DdlsMinScan_pt%.0fto%.0f",DptLow,DptHigh));
		// h_PromptDs_DdlsMinScan[ibin_pt]=(TH1D*)fin->Get(Form("h_PromptDs_DdlsMinScan2_pt%.0fto%.0f",DptLow,DptHigh));
		cout<<"read DdlsMinScan TH1D "<<endl;
		// h_PromptDs_DdlsMinScan_ratio[ibin_pt]=new TH1D(Form("h_PromptDs_DdlsMinScan_ratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_DdlsMinScan_ratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DdlsMinScan, bins_DdlsMinScan );
		h_PromptDs_DdlsMinScan_ratio[ibin_pt]=new TH1D(Form("h_PromptDs_DdlsMinScan_ratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_DdlsMinScan_ratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DdlsMinScan+1, bins_DdlsMinScan_new );
		// h_PromptDs_DdlsMinScan_Stat[ibin_pt]=new TH1D(Form("h_PromptDs_DdlsMinScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_DdlsMinScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),nbin_DdlsMinScan+1, bins_DdlsMinScan_new );
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]=new TH1D(Form("h_PromptDs_DdlsMinScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_DdlsMinScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),1,0, bins_DdlsMinScan_new[nbin_DdlsMinScan+1] );

    h_PromptDs_DdlsMinScan_Yieldratio[ibin_pt]=new TH1D(Form("h_PromptDs_DdlsMinScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_DdlsMinScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DdlsMinScan+1, bins_DdlsMinScan_new );
    h_PromptDs_DdlsMinScan_Effratio[ibin_pt]=new TH1D(Form("h_PromptDs_DdlsMinScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_DdlsMinScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DdlsMinScan+1, bins_DdlsMinScan_new );
    h_NonPromptDs_DdlsMinScan_Effratio[ibin_pt]=new TH1D(Form("h_NonPromptDs_DdlsMinScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_NonPromptDs_DdlsMinScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DdlsMinScan+1, bins_DdlsMinScan_new );
    h_All_DdlsMinScan_Yieldratio[ibin_pt]=new TH1D(Form("h_All_DdlsMinScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_All_DdlsMinScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DdlsMinScan+1, bins_DdlsMinScan_new );


		cout<<"\n ibin_pt = "<<ibin_pt<<" , Ds pt = "<<DptLow<<" to "<<DptHigh<<endl;
		cout<<" cut_Ddls = "<<cut_Ddls<<" , cut_Ddls = "<<cut_Ddls<<" , cut_Ddls = "<<cut_Ddls<<endl;
		cout<<" CSdNdpt = "<<CSdNdpt<<" CSdNdptErr = "<<CSdNdptErr<<" , CSdNdptErrRel = "<<CSdNdptErrRel*100<<"%"<<endl;
		cout<<" \n DdlsMinScan "<<endl;

		start_var=0;
		if(isPbPb==3 &&ibin_pt<=2){ start_var=1;} // don't use first bin result, fitting failed

		for(int ibin_var=start_var; ibin_var<nbin_DdlsMinScan; ibin_var++){
			double CSdNdpt_var=h_PromptDs_DdlsMinScan[ibin_pt]->GetBinContent(ibin_var+1);
			double CSdNdptErr_var=h_PromptDs_DdlsMinScan[ibin_pt]->GetBinError(ibin_var+1);
			double CSdNdptErrRel_var=CSdNdptErr_var/CSdNdpt_var;	

			double cut_DdlsMinScan=DdlsMinScan_bins[ibin_var];

			double ratio=CSdNdpt_var/CSdNdpt;
			double ratioErrRel_temp = CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var - 2/CSdNdpt/CSdNdpt_var*CSdNdptErr*CSdNdptErr;

    double Eff_Prompt_phi_var=h_Eff_Prompt_phikkpi_DdlsMinScan[ibin_var]->GetBinContent(ibin_pt+1);
    double Eff_Prompt_f0_var=h_Eff_Prompt_f0kkpi_DdlsMinScan[ibin_var]->GetBinContent(ibin_pt+1);
    double Eff_NonPrompt_phi_var=h_Eff_NonPrompt_phikkpi_DdlsMinScan[ibin_var]->GetBinContent(ibin_pt+1);
    double Eff_NonPrompt_f0_var=h_Eff_NonPrompt_f0kkpi_DdlsMinScan[ibin_var]->GetBinContent(ibin_pt+1);

    double Eff_Prompt_AllBR_var=h_Eff_Prompt_AllBR_DdlsMinScan[ibin_var]->GetBinContent(ibin_pt+1);
    double EffErr_Prompt_AllBR_var=h_Eff_Prompt_AllBR_DdlsMinScan[ibin_var]->GetBinError(ibin_pt+1);
    double EffErrRel_Prompt_AllBR_var=EffErr_Prompt_AllBR_var/Eff_Prompt_AllBR_var;

    double Eff_NonPrompt_AllBR_var=h_Eff_NonPrompt_AllBR_DdlsMinScan[ibin_var]->GetBinContent(ibin_pt+1);
    double EffErr_NonPrompt_AllBR_var=h_Eff_NonPrompt_AllBR_DdlsMinScan[ibin_var]->GetBinError(ibin_pt+1);
    double EffErrRel_NonPrompt_AllBR_var=EffErr_Prompt_AllBR_var/Eff_NonPrompt_AllBR_var;

    double Prompt_yield_var=CSdNdpt_var*2*LumiNevt*(BRphi*Eff_Prompt_phi_var+BRf0*Eff_Prompt_f0_var);
    double Prompt_yieldErr_var=Prompt_yield_var*CSdNdptErrRel_var;

    double ratio_yield=Prompt_yield_var/Prompt_yield;

    double ratio_Eff=Eff_Prompt_AllBR_var/Eff_Prompt_AllBR;
    double ratioErr_Eff=ratio_Eff*( sqrt(EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR+ EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var)  );

    double ratio_NPEff=Eff_NonPrompt_AllBR_var/Eff_NonPrompt_AllBR;
    double ratioErr_NPEff=ratio_NPEff*( sqrt(EffErrRel_NonPrompt_AllBR*EffErrRel_NonPrompt_AllBR+ EffErrRel_NonPrompt_AllBR_var*EffErrRel_NonPrompt_AllBR_var)  );

    double All_yield_var=h_RawFitYield_DdlsMinScan[ibin_var]->GetBinContent(ibin_pt+1);
    double ratio_All_yield=All_yield_var/All_yield;
    double ratio_All_yieldErr=0;



			// cout<<"ratioErrRel_temp = "<<ratioErrRel_temp<<endl;
			if(cut_Ddls > cut_DdlsMinScan ){
			ratioErrRel_temp = CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var;
			}
			if(ratioErrRel_temp <0 ){
			ratioErrRel_temp = CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var >0 ? CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var : CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var - 2/CSdNdpt/CSdNdpt_var*CSdNdptErr*CSdNdptErr;
			}
			if(ratioErrRel_temp <0 ){
				cout<<" --- warning , ratioErrRel_temp still < 0 !!! "<<endl;
			}
			cout<<"ratioErrRel_temp = "<<ratioErrRel_temp<<endl;
			double ratioErr = ratio*sqrt(ratioErrRel_temp);
      ratioErr = sqrt(abs(CSdNdptErr_var*CSdNdptErr_var-CSdNdptErr*CSdNdptErr))/CSdNdpt;
      double ratioErr_yield=ratio_yield*sqrt(ratioErrRel_temp);
			ratio_All_yieldErr=ratio_All_yield*sqrt(ratioErrRel_temp);


			// single ratio using full error //
			ratioErr_yield=Prompt_yieldErr_var/Prompt_yield;
			ratio_All_yieldErr=h_RawFitYield_DdlsMinScan[ibin_var]->GetBinError(ibin_pt+1)/All_yield;
			ratioErr_Eff=EffErr_Prompt_AllBR_var/Eff_Prompt_AllBR;
			ratioErr_NPEff=EffErr_NonPrompt_AllBR_var/Eff_NonPrompt_AllBR;



			// special case for var_cut = default cut
			if(cut_Ddls == cut_DdlsMinScan){
				ratio=1;
				ratioErr=0;
			}

			// double ratioErr=ratio*sqrt( CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr*CSdNdptErr  );
//			if(cut_Ddls < cut_DdlsMinScan ){
	//			ratioErr=ratio*sqrt( CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var  );			
//			}

			cout<<"\n ibin_var = "<<ibin_var<<" , DdlsMinScan var = "<<cut_DdlsMinScan<<endl;
			cout<<"CSdNdpt_var = "<<CSdNdpt_var<<" , CSdNdptErr_var = "<<CSdNdptErr_var<<" , CSdNdptErrRel_var = "<<CSdNdptErrRel_var*100<<"%"<<endl;
			cout<<"ratio = "<<ratio<<" +- "<<ratioErr<<" relErr = "<<ratioErr/ratio*100<<"%"<<endl;

			// h_PromptDs_DdlsMinScan_ratio[ibin_pt]->SetBinContent(ibin_var+1, ratio);	
			// h_PromptDs_DdlsMinScan_ratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr);	

			h_PromptDs_DdlsMinScan_ratio[ibin_pt]->SetBinContent(ibin_var+2, ratio);	
			h_PromptDs_DdlsMinScan_ratio[ibin_pt]->SetBinError(ibin_var+2, ratioErr);	

			h_PromptDs_DdlsMinScan_Yieldratio[ibin_pt]->SetBinContent(ibin_var+2, ratio_yield);	
			h_PromptDs_DdlsMinScan_Yieldratio[ibin_pt]->SetBinError(ibin_var+2, ratioErr_yield);	

			h_PromptDs_DdlsMinScan_Effratio[ibin_pt]->SetBinContent(ibin_var+2, ratio_Eff);	
			h_PromptDs_DdlsMinScan_Effratio[ibin_pt]->SetBinError(ibin_var+2, ratioErr_Eff);	

			h_NonPromptDs_DdlsMinScan_Effratio[ibin_pt]->SetBinContent(ibin_var+2, ratio_NPEff);	
			h_NonPromptDs_DdlsMinScan_Effratio[ibin_pt]->SetBinError(ibin_var+2, ratioErr_NPEff);	

			h_All_DdlsMinScan_Yieldratio[ibin_pt]->SetBinContent(ibin_var+2, ratio_All_yield);	
			h_All_DdlsMinScan_Yieldratio[ibin_pt]->SetBinError(ibin_var+2, ratio_All_yieldErr);	
			// h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetBinContent(ibin_var+2, 1);	
			// h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetBinError(ibin_var+2, CSdNdptErrRel);	

      if(cut_DdlsMinScan != cut_DdlsMinScan){
      Gr_PromptDs_DdlsMinScan_ratio[ibin_pt]->SetPoint(ibin_var,DdlsMinScan_bins[ibin_var],ratio);
      Gr_PromptDs_DdlsMinScan_ratio[ibin_pt]->SetPointError(ibin_var,0,ratioErr);
      }

		} // end for ibin_var<nbin_DdlsMinScan


			h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetBinContent(1, 1);	
			h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetBinError(1, CSdNdptErrRel);	

		fout->cd();
		h_PromptDs_DdlsMinScan_ratio[ibin_pt]->Write("",TObject::kOverwrite);
	  Gr_PromptDs_DdlsMinScan_ratio[ibin_pt]->Write("",TObject::kOverwrite);	
		// fit			
		// TF1 *f1_pol1=new TF1("f1_pol1","[0]+x*[1]");
		f1_pol1->SetRange(0,8);
		f1_pol1->SetParameter(0,1);
		// f1_pol1->FixParameter(0,1);
		f1_pol1->SetParameter(1,0);
		f1_pol1->SetParameter(2,cut_Ddls);
		f1_pol1->FixParameter(2,cut_Ddls);

		h_PromptDs_DdlsMinScan_ratio[ibin_pt]->Fit("f1_pol1","QN0");
		h_PromptDs_DdlsMinScan_ratio[ibin_pt]->Fit("f1_pol1","QN0");
		h_PromptDs_DdlsMinScan_ratio[ibin_pt]->Fit("f1_pol1","EMS0 ");
		h_PromptDs_DdlsMinScan_ratio[ibin_pt]->Fit("f1_pol1","EMS0");
		// fitR=h_PromptDs_DdlsMinScan_ratio[ibin_pt]->Fit("f1_pol1","QN0 WS");

    fitR=h_PromptDs_DdlsMinScan_ratio[ibin_pt]->Fit("f1_pol1","EMIS0");
    // fitR=h_PromptDs_DdlsMinScan_ratio[ibin_pt]->Fit("f1_pol1","LEMS0","",2,8);
    // fitR=h_PromptDs_DdlsMinScan_ratio[ibin_pt]->Fit("f1_pol1","LEMS0I","",2,8);
    // fitR=h_PromptDs_DdlsMinScan_ratio[ibin_pt]->Fit("f1_pol1","EMS0","",2,8);
    fitR->GetConfidenceIntervals(1,1,1,x0,x0err,0.683,false);


		// draw

		cav[c_count]= new TCanvas(Form("cav%i",c_count),Form("cav%i",c_count),800,800);
		cav[c_count]->cd();
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetMinimum(RatioMin);
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->GetXaxis()->SetTitle("D_{S} Decay Length Significance");
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->GetYaxis()->SetTitle("#sigma_{Varied cut}/#sigma_{Default cut}");
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetTitle("");
	
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetFillColor(kGray);
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetFillStyle(3001);
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetMarkerColor(kGray);
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetMarkerSize(0);
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetLineWidth(0);
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetLineColor(kGray);
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->Draw("SAME E2");
    Gr_PromptDs_DdlsMinScan_ratio[ibin_pt]->SetMarkerSize(0.5);
    Gr_PromptDs_DdlsMinScan_ratio[ibin_pt]->SetMarkerStyle(21);
    Gr_PromptDs_DdlsMinScan_ratio[ibin_pt]->Draw("P");


//		h_PromptDs_DdlsMinScan_ratio[ibin_pt]->GetYaxis()->SetTitle("");
		f1_pol1->SetRange(0,8);
		f1_pol1->Draw("same");

		shiftY=0;
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;
		if(isPbPb==3){
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,"centrality : 0-100%"); shiftY-=oneshift;
		}
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%%",(abs(f1_pol1->GetParameter(0)-1))*100)); shiftY-=oneshift;


		tl_Ddls[ibin_pt]=new TLine(cut_Ddls,RatioMin,cut_Ddls,RatioMax);
		tl_Ddls[ibin_pt]->SetLineColor(4);
		tl_Ddls[ibin_pt]->SetLineStyle(6);
		tl_Ddls[ibin_pt]->Draw("same");


		// SavePlotDirs(cav[c_count],Form("%s_Sys_DdlsMinScan_Graph_pt%.0fto%.0f",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );

		c_count++;



    cav1[c_count]= new TCanvas(Form("cav1_%i",c_count),Form("cav1_%i",c_count),800,800);
    cav1[c_count]->cd();
    h_PromptDs_DdlsMinScan_Yieldratio[ibin_pt]->GetXaxis()->SetRangeUser(0,7);
    h_PromptDs_DdlsMinScan_Yieldratio[ibin_pt]->GetYaxis()->SetTitle("N_{vary cuts}/N_{default cuts}");
    h_PromptDs_DdlsMinScan_Yieldratio[ibin_pt]->GetXaxis()->SetTitle("D_{S} Decay Length Significance");
    h_PromptDs_DdlsMinScan_Yieldratio[ibin_pt]->SetMaximum(RatioMax);
    h_PromptDs_DdlsMinScan_Yieldratio[ibin_pt]->SetMinimum(0.0);
    h_PromptDs_DdlsMinScan_Yieldratio[ibin_pt]->SetLineColor(2);
    h_PromptDs_DdlsMinScan_Yieldratio[ibin_pt]->SetMarkerColor(2);
    h_PromptDs_DdlsMinScan_Yieldratio[ibin_pt]->Draw("SAME");
    h_PromptDs_DdlsMinScan_Effratio[ibin_pt]->SetLineColor(4);
    h_PromptDs_DdlsMinScan_Effratio[ibin_pt]->SetMarkerColor(4);
    h_PromptDs_DdlsMinScan_Effratio[ibin_pt]->Draw("SAME");
    h_All_DdlsMinScan_Yieldratio[ibin_pt]->SetLineColor(1);
    h_All_DdlsMinScan_Yieldratio[ibin_pt]->SetMarkerColor(1);
    h_All_DdlsMinScan_Yieldratio[ibin_pt]->Draw("SAME");
    h_NonPromptDs_DdlsMinScan_Effratio[ibin_pt]->SetLineColor(kMagenta);
    h_NonPromptDs_DdlsMinScan_Effratio[ibin_pt]->SetMarkerColor(kMagenta);
    h_NonPromptDs_DdlsMinScan_Effratio[ibin_pt]->Draw("SAME");


    TLegend *le_DdlsMinScan_SingleRatio=new TLegend(0.7,0.7,0.85,0.85);
    le_DdlsMinScan_SingleRatio->SetBorderSize(0);
    le_DdlsMinScan_SingleRatio->AddEntry(h_All_DdlsMinScan_Yieldratio[ibin_pt],"Data","lp");
    le_DdlsMinScan_SingleRatio->AddEntry(h_PromptDs_DdlsMinScan_Yieldratio[ibin_pt],"Data Prompt","lp");
    le_DdlsMinScan_SingleRatio->AddEntry(h_PromptDs_DdlsMinScan_Effratio[ibin_pt],"Prompt MC","lp");
    le_DdlsMinScan_SingleRatio->AddEntry(h_NonPromptDs_DdlsMinScan_Effratio[ibin_pt],"NonPrompt MC","lp");
    le_DdlsMinScan_SingleRatio->Draw("same");

    shiftY=0;
    tlt->DrawLatexNDC(textposx+0.0,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;

    TLine *tl_Ddls_SigleRatio=new TLine(cut_Ddls,RatioMin,cut_Ddls,RatioMax);
    tl_Ddls_SigleRatio->SetLineColor(4);
    tl_Ddls_SigleRatio->SetLineStyle(6);
    tl_Ddls_SigleRatio->Draw("same");


    SavePlotDirs(cav1[c_count],Form("%s_Sys_DdlsMinScan_pt%.0fto%.0f_sigleRatio",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );




		cav[c_count]= new TCanvas(Form("cav%i",c_count),Form("cav%i",c_count),800,800);
		cav[c_count]->cd();
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetMinimum(RatioMin);
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->GetXaxis()->SetTitle("D_{S} Decay Length Significance");
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->GetYaxis()->SetTitle("#sigma_{Varied cut}/#sigma_{Default cut}");
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetTitle("");
	
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetFillColor(kGray);
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetFillStyle(3001);
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetMarkerColor(kGray);
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetMarkerSize(0);
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetLineWidth(0);
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->SetLineColor(kGray);
		h_PromptDs_DdlsMinScan_Stat[ibin_pt]->Draw("SAME E2");


		h_PromptDs_DdlsMinScan_ratio[ibin_pt]->GetXaxis()->SetRangeUser(0,8);
		h_PromptDs_DdlsMinScan_ratio[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_DdlsMinScan_ratio[ibin_pt]->SetMinimum(RatioMin);
		h_PromptDs_DdlsMinScan_ratio[ibin_pt]->GetXaxis()->SetTitle("D_{S} Decay Length Significance");
		h_PromptDs_DdlsMinScan_ratio[ibin_pt]->SetTitle("");
//		h_PromptDs_DdlsMinScan_ratio[ibin_pt]->GetYaxis()->SetTitle("");
		h_PromptDs_DdlsMinScan_ratio[ibin_pt]->Draw("SAME");
		f1_pol1->SetRange(0,8);
		f1_pol1->Draw("same");

		shiftY=0;
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;
		if(isPbPb==3){
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,"centrality : 0-100%"); shiftY-=oneshift;
		}
	  // tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%%",(abs(f1_pol1->GetParameter(0)-1))*100)); shiftY-=oneshift;
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%% #pm %.1f%%",(abs(f1_pol1->Eval(0)-1))*100, x0err[0]*100)); shiftY-=oneshift;


		tl_Ddls[ibin_pt]=new TLine(cut_Ddls,RatioMin,cut_Ddls,RatioMax);
		tl_Ddls[ibin_pt]->SetLineColor(4);
		tl_Ddls[ibin_pt]->SetLineStyle(6);
		tl_Ddls[ibin_pt]->Draw("same");


		SavePlotDirs(cav[c_count],Form("%s_Sys_DdlsMinScan_pt%.0fto%.0f",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );


	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Slope : %.3f #pm %.3f",f1_pol1->GetParameter(1), f1_pol1->GetParError(1)) ); shiftY-=oneshift;
		SavePlotDirs(cav[c_count],Form("%s_Sys_DdlsMinScan_pt%.0fto%.0f",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape","withSlope",Form("%s",str_PbPb.Data())} );

		c_count++;

		// special handle

		// h_PromptDs_DdlsMinScan_SysRel->SetBinContent(ibin_pt+1,abs(f1_pol1->GetParameter(0)-1));
		h_PromptDs_DdlsMinScan_SysRel->SetBinContent(ibin_pt+1,abs(f1_pol1->Eval(0)-1));



	cout<<"end DdlsMinScan Systematics "<<endl;
		//-- end DdlsMinScan_bins --//




	//-- start PhiMassScan --//
	cout<<"start PhiMassScan Systematics "<<endl;

		h_PromptDs_PhiMassScan[ibin_pt]=(TH1D*)fin->Get(Form("h_PromptDs_PhiMassScan_pt%.0fto%.0f",DptLow,DptHigh));
		// h_PromptDs_PhiMassScan[ibin_pt]=(TH1D*)fin->Get(Form("h_PromptDs_PhiMassScan2_pt%.0fto%.0f",DptLow,DptHigh));
		cout<<"read PhiMassScan TH1D "<<endl;
		h_PromptDs_PhiMassScan_ratio[ibin_pt]=new TH1D(Form("h_PromptDs_PhiMassScan_ratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_PhiMassScan_ratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_PhiMassScan, bins_PhiMassScan );
		h_PromptDs_PhiMassScan_Stat[ibin_pt]=new TH1D(Form("h_PromptDs_PhiMassScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_PhiMassScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),nbin_PhiMassScan, bins_PhiMassScan );

		cout<<"\n ibin_pt = "<<ibin_pt<<" , Ds pt = "<<DptLow<<" to "<<DptHigh<<endl;
		cout<<" cut_PhiMass = "<<cut_PhiMass<<" , cut_PhiMass = "<<cut_PhiMass<<" , cut_Ddls = "<<cut_Ddls<<endl;
		cout<<" CSdNdpt = "<<CSdNdpt<<" CSdNdptErr = "<<CSdNdptErr<<" , CSdNdptErrRel = "<<CSdNdptErrRel*100<<"%"<<endl;
		cout<<" \n PhiMassScan "<<endl;

		start_var=0;
//		if(isPbPb==3 &&ibin_pt<=2){ start_var=1;} // don't use first bin result, fitting failed

		//// --  start of PhiMass 
		double PhiMassScan_maxDiff=0; 
		double PhiMassScan_maxDiffErr=0; 

		for(int ibin_var=start_var; ibin_var<nbin_PhiMassScan; ibin_var++){
			double CSdNdpt_var=h_PromptDs_PhiMassScan[ibin_pt]->GetBinContent(ibin_var+1);
			double CSdNdptErr_var=h_PromptDs_PhiMassScan[ibin_pt]->GetBinError(ibin_var+1);
			double CSdNdptErrRel_var=CSdNdptErr_var/CSdNdpt_var;	

			double cut_PhiMassScan=PhiMassScan_bins[ibin_var];

			double ratio=CSdNdpt_var/CSdNdpt;
			double ratioErrRel_temp = CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var - 2/CSdNdpt/CSdNdpt_var*CSdNdptErr*CSdNdptErr;
			// cout<<"ratioErrRel_temp = "<<ratioErrRel_temp<<endl;
			if(cut_PhiMass > cut_PhiMassScan ){
			ratioErrRel_temp = CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var;
			}
			if(ratioErrRel_temp <0 ){
			ratioErrRel_temp = CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var >0 ? CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var : CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var - 2/CSdNdpt/CSdNdpt_var*CSdNdptErr*CSdNdptErr;
			}
			if(ratioErrRel_temp <0 ){
				cout<<" --- warning , ratioErrRel_temp still < 0 !!! "<<endl;
			}
			cout<<"ratioErrRel_temp = "<<ratioErrRel_temp<<endl;
			double ratioErr = ratio*sqrt(ratioErrRel_temp);
	    ratioErr = sqrt(abs(CSdNdptErr_var*CSdNdptErr_var-CSdNdptErr*CSdNdptErr))/CSdNdpt;


			// special case for var_cut = default cut
			cout<<"cut_PhiMass = "<<cut_PhiMass<<" , cut_PhiMassScan = "<<cut_PhiMassScan<<endl;
			if(cut_PhiMass == cut_PhiMassScan){
				// ratio=1;
				// ratioErr=0;
				cout<<"hello"<<endl; // somewhat not work..
			}

			// double ratioErr=ratio*sqrt( CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr*CSdNdptErr  );
//			if(cut_PhiMass < cut_PhiMassScan ){
	//			ratioErr=ratio*sqrt( CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var  );			
//			}

			cout<<"\n ibin_var = "<<ibin_var<<" , PhiMassScan var = "<<cut_PhiMassScan<<endl;
			cout<<"CSdNdpt_var = "<<CSdNdpt_var<<" , CSdNdptErr_var = "<<CSdNdptErr_var<<" , CSdNdptErrRel_var = "<<CSdNdptErrRel_var*100<<"%"<<endl;
			cout<<"ratio = "<<ratio<<" +- "<<ratioErr<<" relErr = "<<ratioErr/ratio*100<<"%"<<endl;

			h_PromptDs_PhiMassScan_ratio[ibin_pt]->SetBinContent(ibin_var+1, ratio);	
			h_PromptDs_PhiMassScan_ratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr);	
	
			h_PromptDs_PhiMassScan_Stat[ibin_pt]->SetBinContent(ibin_var+1, 1);	
			h_PromptDs_PhiMassScan_Stat[ibin_pt]->SetBinError(ibin_var+1, CSdNdptErrRel);	
      if(abs(ratio-1) >PhiMassScan_maxDiff){
				PhiMassScan_maxDiff = abs(ratio-1); 
				PhiMassScan_maxDiffErr=ratioErr;
			}

		} // end for ibin_var<nbin_PhiMassScan

		h_PromptDs_PhiMassScan_SysRel->SetBinContent(ibin_pt+1,PhiMassScan_maxDiff);
		fout->cd();
		h_PromptDs_PhiMassScan_ratio[ibin_pt]->Write("",TObject::kOverwrite);
		
		// fit			
		// TF1 *f1_pol1=new TF1("f1_pol1","[0]+x*[1]");
	/*
		f1_pol1->SetParameter(0,1);
		f1_pol1->SetParameter(1,0);
		h_PromptDs_PhiMassScan_ratio[ibin_pt]->Fit("f1_pol1","QN0 W");
		h_PromptDs_PhiMassScan_ratio[ibin_pt]->Fit("f1_pol1","QN0 W");
		h_PromptDs_PhiMassScan_ratio[ibin_pt]->Fit("f1_pol1","EMIS0 W");
		h_PromptDs_PhiMassScan_ratio[ibin_pt]->Fit("f1_pol1","EMIS0 W");
*/

		// draw
		cav[c_count]= new TCanvas(Form("cav%i",c_count),Form("cav%i",c_count),800,800);
		cav[c_count]->cd();

		h_PromptDs_PhiMassScan_Stat[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_PhiMassScan_Stat[ibin_pt]->SetMinimum(RatioMin);
		h_PromptDs_PhiMassScan_Stat[ibin_pt]->GetXaxis()->SetTitle("#phi Mass window");
		h_PromptDs_PhiMassScan_Stat[ibin_pt]->GetYaxis()->SetTitle("#sigma_{Varied cut}/#sigma_{Default cut}");
		h_PromptDs_PhiMassScan_Stat[ibin_pt]->SetTitle("");
	
		h_PromptDs_PhiMassScan_Stat[ibin_pt]->SetFillColor(kGray);
		h_PromptDs_PhiMassScan_Stat[ibin_pt]->SetFillStyle(3001);
		h_PromptDs_PhiMassScan_Stat[ibin_pt]->SetMarkerColor(kGray);
		h_PromptDs_PhiMassScan_Stat[ibin_pt]->SetMarkerSize(0);
		h_PromptDs_PhiMassScan_Stat[ibin_pt]->SetLineWidth(0);
		h_PromptDs_PhiMassScan_Stat[ibin_pt]->SetLineColor(kGray);
		h_PromptDs_PhiMassScan_Stat[ibin_pt]->Draw("SAME E2");



		h_PromptDs_PhiMassScan_ratio[ibin_pt]->GetXaxis()->SetRangeUser(8,12);
		h_PromptDs_PhiMassScan_ratio[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_PhiMassScan_ratio[ibin_pt]->SetMinimum(RatioMin);
		h_PromptDs_PhiMassScan_ratio[ibin_pt]->GetXaxis()->SetTitle("#phi Mass window");
		h_PromptDs_PhiMassScan_ratio[ibin_pt]->SetTitle("");
//		h_PromptDs_PhiMassScan_ratio[ibin_pt]->GetYaxis()->SetTitle("");
		h_PromptDs_PhiMassScan_ratio[ibin_pt]->Draw("SAME");
//	f1_pol1->SetRange(9,11);
//		f1_pol1->Draw("same");
		tl2->Draw("same");

		shiftY=0;
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;
		if(isPbPb==3){
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,"centrality : 0-100%"); shiftY-=oneshift;
		}
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%% #pm %.1f%%",(PhiMassScan_maxDiff*100), PhiMassScan_maxDiffErr*100)); shiftY-=oneshift;

		tl_PhiMass[ibin_pt]=new TLine(cut_PhiMass,RatioMin,cut_PhiMass,RatioMax);
		tl_PhiMass[ibin_pt]->SetLineColor(4);
		tl_PhiMass[ibin_pt]->SetLineStyle(6);
		tl_PhiMass[ibin_pt]->Draw("same");


		SavePlotDirs(cav[c_count],Form("%s_Sys_PhiMassScan_pt%.0fto%.0f",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );

		c_count++;

		// special handle

		// h_PromptDs_PhiMassScan_SysRel->SetBinContent(ibin_pt+1, abs(f1_pol1->GetParameter(0)-1));

	cout<<"end PhiMassScan Systematics "<<endl;
		//-- end PhiMassScan_bins --//



	if(doDauPtCut){

		//// --  start of KaonPt /// 
	cout<<"start KaonPtScan Systematics "<<endl;

		// h_PromptDs_KaonPtScan[ibin_pt]=(TH1D*)fin->Get(Form("h_PromptDs_KaonPtScan_pt%.0fto%.0f",DptLow,DptHigh));
		h_PromptDs_KaonPtScan[ibin_pt]=(TH1D*)fin->Get(Form("h_PromptDs_KaonPtScan2_pt%.0fto%.0f",DptLow,DptHigh));
		cout<<"read KaonPtScan TH1D "<<endl;
		h_PromptDs_KaonPtScan_ratio[ibin_pt]=new TH1D(Form("h_PromptDs_KaonPtScan_ratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_KaonPtScan_ratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );

		h_PromptDs_KaonPtScan_Yieldratio[ibin_pt]=new TH1D(Form("h_PromptDs_KaonPtScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_KaonPtScan_Yiledratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );
		h_PromptDs_KaonPtScan_Effratio[ibin_pt]=new TH1D(Form("h_PromptDs_KaonPtScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_KaonPtScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );
		h_NonPromptDs_KaonPtScan_Effratio[ibin_pt]=new TH1D(Form("h_NonPromptDs_KaonPtScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_NonPromptDs_KaonPtScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );
    Gr_PromptDs_KaonPtScan_ratio[ibin_pt]=new TGraphErrors();

		h_PromptDs_KaonPtScan_MCratio[ibin_pt]=new TH1D(Form("h_PromptDs_KaonPtScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_KaonPtScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );
		h_NonPromptDs_KaonPtScan_MCratio[ibin_pt]=new TH1D(Form("h_NonPromptDs_KaonPtScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_NonPromptDs_KaonPtScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );
		h_All_KaonPtScan_Yieldratio[ibin_pt]=new TH1D(Form("h_All_KaonPtScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_All_KaonPtScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );

		// h_PromptDs_KaonPtScan_Stat[ibin_pt]=new TH1D(Form("h_PromptDs_KaonPtScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_KaonPtScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );
		h_PromptDs_KaonPtScan_Stat[ibin_pt]=new TH1D(Form("h_PromptDs_KaonPtScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_KaonPtScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),1,0, bins_DauPtScan[nbin_DauPtScan] );


	cout<<__LINE__<<endl;
		// TH1D *h_eff_KaonPtScan[nbin_DauPtScan];
		
		for(int ibin_var=start_var; ibin_var<nbin_DauPtScan; ibin_var++){

			double CSdNdpt_var=h_PromptDs_KaonPtScan[ibin_pt]->GetBinContent(ibin_var+1);
			double CSdNdptErr_var=h_PromptDs_KaonPtScan[ibin_pt]->GetBinError(ibin_var+1);
			double CSdNdptErrRel_var=CSdNdptErr_var/CSdNdpt_var;	

			double cut_KaonPtScan=DauPtScan_bins[ibin_var];

			double ratio=CSdNdpt_var/CSdNdpt;
			double ratioErrRel_temp = CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var - 2/CSdNdpt/CSdNdpt_var*CSdNdptErr*CSdNdptErr;
			// cout<<"ratioErrRel_temp = "<<ratioErrRel_temp<<endl;
	// for single ratio 

		double Eff_Prompt_phi_var=h_Eff_Prompt_phikkpi_KaonPtScan[ibin_var]->GetBinContent(ibin_pt+1);
		double Eff_Prompt_f0_var=h_Eff_Prompt_f0kkpi_KaonPtScan[ibin_var]->GetBinContent(ibin_pt+1);
		double Eff_NonPrompt_phi_var=h_Eff_NonPrompt_phikkpi_KaonPtScan[ibin_var]->GetBinContent(ibin_pt+1);
		double Eff_NonPrompt_f0_var=h_Eff_NonPrompt_f0kkpi_KaonPtScan[ibin_var]->GetBinContent(ibin_pt+1);
	
		double Eff_Prompt_AllBR_var=h_Eff_Prompt_AllBR_KaonPtScan[ibin_var]->GetBinContent(ibin_pt+1);
		double EffErr_Prompt_AllBR_var=h_Eff_Prompt_AllBR_KaonPtScan[ibin_var]->GetBinError(ibin_pt+1);
		double EffErrRel_Prompt_AllBR_var=EffErr_Prompt_AllBR_var/Eff_Prompt_AllBR_var;

		double Eff_NonPrompt_AllBR_var=h_Eff_NonPrompt_AllBR_KaonPtScan[ibin_var]->GetBinContent(ibin_pt+1);
		double EffErr_NonPrompt_AllBR_var=h_Eff_NonPrompt_AllBR_KaonPtScan[ibin_var]->GetBinError(ibin_pt+1);
		double EffErrRel_NonPrompt_AllBR_var=EffErr_NonPrompt_AllBR_var/Eff_NonPrompt_AllBR_var;

	cout<<__LINE__<<endl;

		double Prompt_yield_var=CSdNdpt_var*2*LumiNevt*(BRphi*Eff_Prompt_phi_var+BRf0*Eff_Prompt_f0_var);
		double Prompt_yieldErr_var=Prompt_yield_var*CSdNdptErrRel_var;

		double ratio_yield=Prompt_yield_var/Prompt_yield;

		double ratio_Eff=Eff_Prompt_AllBR_var/Eff_Prompt_AllBR;
		double ratioErr_Eff=ratio_Eff*( sqrt(EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR+ EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var)  );


		double ratio_NPEff=Eff_NonPrompt_AllBR_var/Eff_NonPrompt_AllBR;
		double ratioErr_NPEff=ratio_NPEff*( sqrt(EffErrRel_NonPrompt_AllBR*EffErrRel_NonPrompt_AllBR+ EffErrRel_NonPrompt_AllBR_var*EffErrRel_NonPrompt_AllBR_var)  );


		double All_yield_var=h_RawFitYield_KaonPtScan[ibin_var]->GetBinContent(ibin_pt+1);
		double ratio_All_yield=All_yield_var/All_yield;
		double ratio_All_yieldErr=0;

	cout<<__LINE__<<endl;
		double ratioErrRel_temp_Eff=EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR + EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var-2/Eff_Prompt_AllBR/EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR*EffErr_Prompt_AllBR;


//		double ratioErr_yield=ratio_yield*
			
	// end single ratio


	cout<<__LINE__<<endl;

			// double ratioErr_Eff=ratio_Eff*sqrt(ratioErrRel_temp_Eff);
			double ratioErr = sqrt(abs(CSdNdptErr_var*CSdNdptErr_var-CSdNdptErr*CSdNdptErr))/CSdNdpt;	
			double ratioErr_yield=ratio_yield*sqrt(ratioErrRel_temp);
			ratio_All_yieldErr=ratio_All_yield*sqrt(ratioErrRel_temp);

	//// sigle ratio use full error //
			ratioErr_yield=Prompt_yieldErr_var/Prompt_yield;
			ratio_All_yieldErr=h_RawFitYield_KaonPtScan[ibin_var]->GetBinError(ibin_pt+1)/All_yield;
			ratioErr_Eff=EffErr_Prompt_AllBR_var/Eff_Prompt_AllBR;
			ratioErr_NPEff=EffErr_NonPrompt_AllBR_var/Eff_NonPrompt_AllBR;

			cout<<"ratioErr_Eff = "<<ratioErr_Eff<<endl;

			// special case for var_cut = default cut
			if(cut_DauPt == cut_KaonPtScan){
				ratio=1;
				ratioErr=0;
			}

			cout<<"\n ibin_var = "<<ibin_var<<" , KaonPtScan var = "<<cut_KaonPtScan<<endl;
			cout<<"CSdNdpt_var = "<<CSdNdpt_var<<" , CSdNdptErr_var = "<<CSdNdptErr_var<<" , CSdNdptErrRel_var = "<<CSdNdptErrRel_var*100<<"%"<<endl;
			cout<<"ratio = "<<ratio<<" +- "<<ratioErr<<" relErr = "<<ratioErr/ratio*100<<"%"<<endl;

			h_PromptDs_KaonPtScan_ratio[ibin_pt]->SetBinContent(ibin_var+1, ratio);	
			h_PromptDs_KaonPtScan_ratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr);	

			h_PromptDs_KaonPtScan_Yieldratio[ibin_pt]->SetBinContent(ibin_var+1 , ratio_yield);
			h_PromptDs_KaonPtScan_Yieldratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr_yield);
			h_PromptDs_KaonPtScan_Effratio[ibin_pt]->SetBinContent(ibin_var+1, ratio_Eff);
			h_PromptDs_KaonPtScan_Effratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr_Eff);

			h_NonPromptDs_KaonPtScan_Effratio[ibin_pt]->SetBinContent(ibin_var+1, ratio_NPEff);
			h_NonPromptDs_KaonPtScan_Effratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr_NPEff);
	
			h_All_KaonPtScan_Yieldratio[ibin_pt]->SetBinContent(ibin_var+1, ratio_All_yield);
			h_All_KaonPtScan_Yieldratio[ibin_pt]->SetBinError(ibin_var+1, ratio_All_yieldErr);

      if(cut_DauPt != cut_KaonPtScan){
      Gr_PromptDs_KaonPtScan_ratio[ibin_pt]->SetPoint(ibin_var,DauPtScan_bins[ibin_var],ratio);
      Gr_PromptDs_KaonPtScan_ratio[ibin_pt]->SetPointError(ibin_var,0,ratioErr);
      }


		} // end for ibin_var<nbin_DauPtScan


			h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetBinContent(1, 1);	
			h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetBinError(1, CSdNdptErrRel);	

		fout->cd();
		h_PromptDs_KaonPtScan_ratio[ibin_pt]->Write("",TObject::kOverwrite);
		Gr_PromptDs_KaonPtScan_ratio[ibin_pt]->Write("",TObject::kOverwrite);
		
		// TF1 *f1_pol1=new TF1("f1_pol1","[0]+x*[1]");
		f1_pol1->SetParameter(0,1);
		// f1_pol1->FixParameter(0,1);
		f1_pol1->SetParameter(1,0);
		f1_pol1->SetParameter(2,cut_DauPt);
		f1_pol1->FixParameter(2,cut_DauPt);

		h_PromptDs_KaonPtScan_ratio[ibin_pt]->Fit("f1_pol1","QN0 W");
		h_PromptDs_KaonPtScan_ratio[ibin_pt]->Fit("f1_pol1","QN0 W");
		h_PromptDs_KaonPtScan_ratio[ibin_pt]->Fit("f1_pol1","EMIS0 W");
		h_PromptDs_KaonPtScan_ratio[ibin_pt]->Fit("f1_pol1","EMIS0 W");


    cav[c_count]= new TCanvas(Form("cav%i",c_count),Form("cav%i",c_count),800,800);
    cav[c_count]->cd();
    h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetMaximum(RatioMax);
    h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetMinimum(RatioMin);
    h_PromptDs_KaonPtScan_Stat[ibin_pt]->GetXaxis()->SetTitle("Kaon p_{T}");
    h_PromptDs_KaonPtScan_Stat[ibin_pt]->GetYaxis()->SetTitle("#sigma_{Varied cut}/#sigma_{Default cut}");
    h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetTitle("");

    h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetFillColor(kGray);
    h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetFillStyle(3001);
    h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetMarkerColor(kGray);
    h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetMarkerSize(0);
    h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetLineWidth(0);
    h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetLineColor(kGray);
    h_PromptDs_KaonPtScan_Stat[ibin_pt]->Draw("SAME E2");
    Gr_PromptDs_KaonPtScan_ratio[ibin_pt]->SetMarkerSize(0.5);
    Gr_PromptDs_KaonPtScan_ratio[ibin_pt]->SetMarkerStyle(21);
    Gr_PromptDs_KaonPtScan_ratio[ibin_pt]->Draw("P");

    f1_pol1->SetRange(0,3.75);
    f1_pol1->Draw("same");

    shiftY=0;
    tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;
    if(isPbPb==3){
    tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,"centrality : 0-100%"); shiftY-=oneshift;
    }
    tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%%",(abs(f1_pol1->GetParameter(0)-1))*100)); shiftY-=oneshift;

    tl_KaonPt[ibin_pt]=new TLine(cut_DauPt,RatioMin,cut_DauPt,RatioMax);
    tl_KaonPt[ibin_pt]->SetLineColor(4);
    tl_KaonPt[ibin_pt]->SetLineStyle(6);
    tl_KaonPt[ibin_pt]->Draw("same");


    // SavePlotDirs(cav[c_count],Form("%s_Sys_KaonPtScan_Graph_pt%.0fto%.0f",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );

    c_count++;



		cav1[c_count]= new TCanvas(Form("cav1_%i",c_count),Form("cav1_%i",c_count),800,800);
		cav1[c_count]->cd();
		h_PromptDs_KaonPtScan_Yieldratio[ibin_pt]->GetXaxis()->SetRangeUser(0,0.5);
		h_PromptDs_KaonPtScan_Yieldratio[ibin_pt]->GetYaxis()->SetTitle("N_{vary cuts}/N_{default cuts}");
		h_PromptDs_KaonPtScan_Yieldratio[ibin_pt]->GetXaxis()->SetTitle("Kaon p_{T}");
		h_PromptDs_KaonPtScan_Yieldratio[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_KaonPtScan_Yieldratio[ibin_pt]->SetMinimum(0.0);
		h_PromptDs_KaonPtScan_Yieldratio[ibin_pt]->SetLineColor(2);
		h_PromptDs_KaonPtScan_Yieldratio[ibin_pt]->SetMarkerColor(2);
		h_PromptDs_KaonPtScan_Yieldratio[ibin_pt]->Draw("SAME");
		h_PromptDs_KaonPtScan_Effratio[ibin_pt]->SetLineColor(4);
		h_PromptDs_KaonPtScan_Effratio[ibin_pt]->SetMarkerColor(4);
		h_PromptDs_KaonPtScan_Effratio[ibin_pt]->Draw("SAME");
		h_All_KaonPtScan_Yieldratio[ibin_pt]->SetLineColor(1);
		h_All_KaonPtScan_Yieldratio[ibin_pt]->SetMarkerColor(1);
		h_All_KaonPtScan_Yieldratio[ibin_pt]->Draw("SAME");
		h_NonPromptDs_KaonPtScan_Effratio[ibin_pt]->SetLineColor(kMagenta);
		h_NonPromptDs_KaonPtScan_Effratio[ibin_pt]->SetMarkerColor(kMagenta);
		h_NonPromptDs_KaonPtScan_Effratio[ibin_pt]->Draw("SAME");


		TLegend *le_KaonPtScan_SingleRatio=new TLegend(0.7,0.7,0.85,0.85);
		le_KaonPtScan_SingleRatio->SetBorderSize(0);	
		le_KaonPtScan_SingleRatio->AddEntry(h_All_KaonPtScan_Yieldratio[ibin_pt],"Data","lp");
		le_KaonPtScan_SingleRatio->AddEntry(h_PromptDs_KaonPtScan_Yieldratio[ibin_pt],"Data Prompt","lp");
		le_KaonPtScan_SingleRatio->AddEntry(h_PromptDs_KaonPtScan_Effratio[ibin_pt],"Prompt MC","lp");
		le_KaonPtScan_SingleRatio->AddEntry(h_NonPromptDs_KaonPtScan_Effratio[ibin_pt],"NonPrompt MC","lp");
		le_KaonPtScan_SingleRatio->Draw("same");

		shiftY=0;
	  tlt->DrawLatexNDC(textposx+0.0,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;

		TLine *tl_KaonPt_SigleRatio=new TLine(cut_DauPt,RatioMin,cut_DauPt,RatioMax);
		tl_KaonPt_SigleRatio->SetLineColor(4);
		tl_KaonPt_SigleRatio->SetLineStyle(6);
		tl_KaonPt_SigleRatio->Draw("same");
	

		SavePlotDirs(cav1[c_count],Form("%s_Sys_KaonPtScan_pt%.0fto%.0f_sigleRatio",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );


		// draw
		cav[c_count]= new TCanvas(Form("cav%i",c_count),Form("cav%i",c_count),800,800);
		cav[c_count]->cd();

		h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetMinimum(RatioMin);
		h_PromptDs_KaonPtScan_Stat[ibin_pt]->GetXaxis()->SetTitle("Kaon p_{T}");
		h_PromptDs_KaonPtScan_Stat[ibin_pt]->GetYaxis()->SetTitle("#sigma_{Varied cut}/#sigma_{Default cut}");
		h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetTitle("");
	
		h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetFillColor(kGray);
		h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetFillStyle(3001);
		h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetMarkerColor(kGray);
		h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetMarkerSize(0);
		h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetLineWidth(0);
		h_PromptDs_KaonPtScan_Stat[ibin_pt]->SetLineColor(kGray);
		h_PromptDs_KaonPtScan_Stat[ibin_pt]->Draw("SAME E2");


		h_PromptDs_KaonPtScan_ratio[ibin_pt]->GetXaxis()->SetRangeUser(0,0.6);
		h_PromptDs_KaonPtScan_ratio[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_KaonPtScan_ratio[ibin_pt]->SetMinimum(RatioMin);
		h_PromptDs_KaonPtScan_ratio[ibin_pt]->GetXaxis()->SetTitle("Kaon p_{T}");
		h_PromptDs_KaonPtScan_ratio[ibin_pt]->SetTitle("");
//		h_PromptDs_KaonPtScan_ratio[ibin_pt]->GetYaxis()->SetTitle("");
		h_PromptDs_KaonPtScan_ratio[ibin_pt]->Draw("SAME");
		f1_pol1->SetRange(0,3.75);
		f1_pol1->Draw("same");

		shiftY=0;
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;
		if(isPbPb==3){
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,"centrality : 0-100%"); shiftY-=oneshift;
		}
	  // tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%%",(abs(f1_pol1->GetParameter(0)-1))*100)); shiftY-=oneshift;
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%%",(abs(f1_pol1->Eval(0)-1))*100)); shiftY-=oneshift;

		tl_KaonPt[ibin_pt]=new TLine(cut_DauPt,RatioMin,cut_DauPt,RatioMax);
		tl_KaonPt[ibin_pt]->SetLineColor(4);
		tl_KaonPt[ibin_pt]->SetLineStyle(6);
		tl_KaonPt[ibin_pt]->Draw("same");
	

		SavePlotDirs(cav[c_count],Form("%s_Sys_KaonPtScan_pt%.0fto%.0f",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );

		c_count++;

		// special handle

		// h_PromptDs_KaonPtScan_SysRel->SetBinContent(ibin_pt+1, abs(f1_pol1->GetParameter(0)-1));
		// h_PromptDs_KaonPtScan_SysRel->SetBinContent(ibin_pt+1, abs(f1_pol1->Eval(0)-1));

	cout<<"end KaonPtScan Systematics "<<endl;
		//-- end KaonPtScan_bins --//


		//// --  start of PionPt /// 
	cout<<"start PionPtScan Systematics "<<endl;

		// h_PromptDs_PionPtScan[ibin_pt]=(TH1D*)fin->Get(Form("h_PromptDs_PionPtScan_pt%.0fto%.0f",DptLow,DptHigh));
		h_PromptDs_PionPtScan[ibin_pt]=(TH1D*)fin->Get(Form("h_PromptDs_PionPtScan2_pt%.0fto%.0f",DptLow,DptHigh));
		cout<<"read PionPtScan TH1D "<<endl;
		h_PromptDs_PionPtScan_ratio[ibin_pt]=new TH1D(Form("h_PromptDs_PionPtScan_ratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_PionPtScan_ratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );

		h_PromptDs_PionPtScan_Yieldratio[ibin_pt]=new TH1D(Form("h_PromptDs_PionPtScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_PionPtScan_Yiledratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );
		h_PromptDs_PionPtScan_Effratio[ibin_pt]=new TH1D(Form("h_PromptDs_PionPtScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_PionPtScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );
		h_NonPromptDs_PionPtScan_Effratio[ibin_pt]=new TH1D(Form("h_NonPromptDs_PionPtScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_NonPromptDs_PionPtScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );
    Gr_PromptDs_PionPtScan_ratio[ibin_pt]=new TGraphErrors();

		h_PromptDs_PionPtScan_MCratio[ibin_pt]=new TH1D(Form("h_PromptDs_PionPtScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_PionPtScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );
		h_NonPromptDs_PionPtScan_MCratio[ibin_pt]=new TH1D(Form("h_NonPromptDs_PionPtScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_NonPromptDs_PionPtScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );
		h_All_PionPtScan_Yieldratio[ibin_pt]=new TH1D(Form("h_All_PionPtScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_All_PionPtScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );

		// h_PromptDs_PionPtScan_Stat[ibin_pt]=new TH1D(Form("h_PromptDs_PionPtScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_PionPtScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );
		h_PromptDs_PionPtScan_Stat[ibin_pt]=new TH1D(Form("h_PromptDs_PionPtScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_PionPtScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),1,0, bins_DauPtScan[nbin_DauPtScan] );


	cout<<__LINE__<<endl;
		// TH1D *h_eff_PionPtScan[nbin_DauPtScan];
		
		for(int ibin_var=start_var; ibin_var<nbin_DauPtScan; ibin_var++){

			double CSdNdpt_var=h_PromptDs_PionPtScan[ibin_pt]->GetBinContent(ibin_var+1);
			double CSdNdptErr_var=h_PromptDs_PionPtScan[ibin_pt]->GetBinError(ibin_var+1);
			double CSdNdptErrRel_var=CSdNdptErr_var/CSdNdpt_var;	

			double cut_PionPtScan=DauPtScan_bins[ibin_var];

			double ratio=CSdNdpt_var/CSdNdpt;
			double ratioErrRel_temp = CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var - 2/CSdNdpt/CSdNdpt_var*CSdNdptErr*CSdNdptErr;
			// cout<<"ratioErrRel_temp = "<<ratioErrRel_temp<<endl;
	// for single ratio 

		double Eff_Prompt_phi_var=h_Eff_Prompt_phikkpi_PionPtScan[ibin_var]->GetBinContent(ibin_pt+1);
		double Eff_Prompt_f0_var=h_Eff_Prompt_f0kkpi_PionPtScan[ibin_var]->GetBinContent(ibin_pt+1);
		double Eff_NonPrompt_phi_var=h_Eff_NonPrompt_phikkpi_PionPtScan[ibin_var]->GetBinContent(ibin_pt+1);
		double Eff_NonPrompt_f0_var=h_Eff_NonPrompt_f0kkpi_PionPtScan[ibin_var]->GetBinContent(ibin_pt+1);
	
		double Eff_Prompt_AllBR_var=h_Eff_Prompt_AllBR_PionPtScan[ibin_var]->GetBinContent(ibin_pt+1);
		double EffErr_Prompt_AllBR_var=h_Eff_Prompt_AllBR_PionPtScan[ibin_var]->GetBinError(ibin_pt+1);
		double EffErrRel_Prompt_AllBR_var=EffErr_Prompt_AllBR_var/Eff_Prompt_AllBR_var;

		double Eff_NonPrompt_AllBR_var=h_Eff_NonPrompt_AllBR_PionPtScan[ibin_var]->GetBinContent(ibin_pt+1);
		double EffErr_NonPrompt_AllBR_var=h_Eff_NonPrompt_AllBR_PionPtScan[ibin_var]->GetBinError(ibin_pt+1);
		double EffErrRel_NonPrompt_AllBR_var=EffErr_NonPrompt_AllBR_var/Eff_NonPrompt_AllBR_var;

	cout<<__LINE__<<endl;

		double Prompt_yield_var=CSdNdpt_var*2*LumiNevt*(BRphi*Eff_Prompt_phi_var+BRf0*Eff_Prompt_f0_var);
		double Prompt_yieldErr_var=Prompt_yield_var*CSdNdptErrRel_var;

		double ratio_yield=Prompt_yield_var/Prompt_yield;

		double ratio_Eff=Eff_Prompt_AllBR_var/Eff_Prompt_AllBR;
		double ratioErr_Eff=ratio_Eff*( sqrt(EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR+ EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var)  );


		double ratio_NPEff=Eff_NonPrompt_AllBR_var/Eff_NonPrompt_AllBR;
		double ratioErr_NPEff=ratio_NPEff*( sqrt(EffErrRel_NonPrompt_AllBR*EffErrRel_NonPrompt_AllBR+ EffErrRel_NonPrompt_AllBR_var*EffErrRel_NonPrompt_AllBR_var)  );


		double All_yield_var=h_RawFitYield_PionPtScan[ibin_var]->GetBinContent(ibin_pt+1);
		double ratio_All_yield=All_yield_var/All_yield;
		double ratio_All_yieldErr=0;

	cout<<__LINE__<<endl;
		double ratioErrRel_temp_Eff=EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR + EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var-2/Eff_Prompt_AllBR/EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR*EffErr_Prompt_AllBR;


//		double ratioErr_yield=ratio_yield*
			
	// end single ratio


	cout<<__LINE__<<endl;

			// double ratioErr_Eff=ratio_Eff*sqrt(ratioErrRel_temp_Eff);
			double ratioErr = sqrt(abs(CSdNdptErr_var*CSdNdptErr_var-CSdNdptErr*CSdNdptErr))/CSdNdpt;	
			double ratioErr_yield=ratio_yield*sqrt(ratioErrRel_temp);
			ratio_All_yieldErr=ratio_All_yield*sqrt(ratioErrRel_temp);

	//// sigle ratio use full error //
			ratioErr_yield=Prompt_yieldErr_var/Prompt_yield;
			ratio_All_yieldErr=h_RawFitYield_PionPtScan[ibin_var]->GetBinError(ibin_pt+1)/All_yield;
			ratioErr_Eff=EffErr_Prompt_AllBR_var/Eff_Prompt_AllBR;
			ratioErr_NPEff=EffErr_NonPrompt_AllBR_var/Eff_NonPrompt_AllBR;

			cout<<"ratioErr_Eff = "<<ratioErr_Eff<<endl;

			// special case for var_cut = default cut
			if(cut_DauPt == cut_PionPtScan){
				ratio=1;
				ratioErr=0;
			}

			cout<<"\n ibin_var = "<<ibin_var<<" , PionPtScan var = "<<cut_PionPtScan<<endl;
			cout<<"CSdNdpt_var = "<<CSdNdpt_var<<" , CSdNdptErr_var = "<<CSdNdptErr_var<<" , CSdNdptErrRel_var = "<<CSdNdptErrRel_var*100<<"%"<<endl;
			cout<<"ratio = "<<ratio<<" +- "<<ratioErr<<" relErr = "<<ratioErr/ratio*100<<"%"<<endl;

			h_PromptDs_PionPtScan_ratio[ibin_pt]->SetBinContent(ibin_var+1, ratio);	
			h_PromptDs_PionPtScan_ratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr);	

			h_PromptDs_PionPtScan_Yieldratio[ibin_pt]->SetBinContent(ibin_var+1 , ratio_yield);
			h_PromptDs_PionPtScan_Yieldratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr_yield);
			h_PromptDs_PionPtScan_Effratio[ibin_pt]->SetBinContent(ibin_var+1, ratio_Eff);
			h_PromptDs_PionPtScan_Effratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr_Eff);

			h_NonPromptDs_PionPtScan_Effratio[ibin_pt]->SetBinContent(ibin_var+1, ratio_NPEff);
			h_NonPromptDs_PionPtScan_Effratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr_NPEff);
	
			h_All_PionPtScan_Yieldratio[ibin_pt]->SetBinContent(ibin_var+1, ratio_All_yield);
			h_All_PionPtScan_Yieldratio[ibin_pt]->SetBinError(ibin_var+1, ratio_All_yieldErr);

      if(cut_DauPt != cut_PionPtScan){
      Gr_PromptDs_PionPtScan_ratio[ibin_pt]->SetPoint(ibin_var,DauPtScan_bins[ibin_var],ratio);
      Gr_PromptDs_PionPtScan_ratio[ibin_pt]->SetPointError(ibin_var,0,ratioErr);
      }


		} // end for ibin_var<nbin_DauPtScan


			h_PromptDs_PionPtScan_Stat[ibin_pt]->SetBinContent(1, 1);	
			h_PromptDs_PionPtScan_Stat[ibin_pt]->SetBinError(1, CSdNdptErrRel);	

		fout->cd();
		h_PromptDs_PionPtScan_ratio[ibin_pt]->Write("",TObject::kOverwrite);
		Gr_PromptDs_PionPtScan_ratio[ibin_pt]->Write("",TObject::kOverwrite);
		
		// TF1 *f1_pol1=new TF1("f1_pol1","[0]+x*[1]");
		f1_pol1->SetParameter(0,1);
		// f1_pol1->FixParameter(0,1);
		f1_pol1->SetParameter(1,0);
		f1_pol1->SetParameter(2,cut_DauPt);
		f1_pol1->FixParameter(2,cut_DauPt);

		h_PromptDs_PionPtScan_ratio[ibin_pt]->Fit("f1_pol1","QN0 W");
		h_PromptDs_PionPtScan_ratio[ibin_pt]->Fit("f1_pol1","QN0 W");
		h_PromptDs_PionPtScan_ratio[ibin_pt]->Fit("f1_pol1","EMIS0 W");
		h_PromptDs_PionPtScan_ratio[ibin_pt]->Fit("f1_pol1","EMIS0 W");


    cav[c_count]= new TCanvas(Form("cav%i",c_count),Form("cav%i",c_count),800,800);
    cav[c_count]->cd();
    h_PromptDs_PionPtScan_Stat[ibin_pt]->SetMaximum(RatioMax);
    h_PromptDs_PionPtScan_Stat[ibin_pt]->SetMinimum(RatioMin);
    h_PromptDs_PionPtScan_Stat[ibin_pt]->GetXaxis()->SetTitle("Pion p_{T}");
    h_PromptDs_PionPtScan_Stat[ibin_pt]->GetYaxis()->SetTitle("#sigma_{Varied cut}/#sigma_{Default cut}");
    h_PromptDs_PionPtScan_Stat[ibin_pt]->SetTitle("");

    h_PromptDs_PionPtScan_Stat[ibin_pt]->SetFillColor(kGray);
    h_PromptDs_PionPtScan_Stat[ibin_pt]->SetFillStyle(3001);
    h_PromptDs_PionPtScan_Stat[ibin_pt]->SetMarkerColor(kGray);
    h_PromptDs_PionPtScan_Stat[ibin_pt]->SetMarkerSize(0);
    h_PromptDs_PionPtScan_Stat[ibin_pt]->SetLineWidth(0);
    h_PromptDs_PionPtScan_Stat[ibin_pt]->SetLineColor(kGray);
    h_PromptDs_PionPtScan_Stat[ibin_pt]->Draw("SAME E2");
    Gr_PromptDs_PionPtScan_ratio[ibin_pt]->SetMarkerSize(0.5);
    Gr_PromptDs_PionPtScan_ratio[ibin_pt]->SetMarkerStyle(21);
    Gr_PromptDs_PionPtScan_ratio[ibin_pt]->Draw("P");

    f1_pol1->SetRange(0,3.75);
    f1_pol1->Draw("same");

    shiftY=0;
    tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;
    if(isPbPb==3){
    tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,"centrality : 0-100%"); shiftY-=oneshift;
    }
    tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%%",(abs(f1_pol1->GetParameter(0)-1))*100)); shiftY-=oneshift;

    tl_PionPt[ibin_pt]=new TLine(cut_DauPt,RatioMin,cut_DauPt,RatioMax);
    tl_PionPt[ibin_pt]->SetLineColor(4);
    tl_PionPt[ibin_pt]->SetLineStyle(6);
    tl_PionPt[ibin_pt]->Draw("same");


    // SavePlotDirs(cav[c_count],Form("%s_Sys_PionPtScan_Graph_pt%.0fto%.0f",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );

    c_count++;



		cav1[c_count]= new TCanvas(Form("cav1_%i",c_count),Form("cav1_%i",c_count),800,800);
		cav1[c_count]->cd();
		h_PromptDs_PionPtScan_Yieldratio[ibin_pt]->GetXaxis()->SetRangeUser(0,0.5);
		h_PromptDs_PionPtScan_Yieldratio[ibin_pt]->GetYaxis()->SetTitle("N_{vary cuts}/N_{default cuts}");
		h_PromptDs_PionPtScan_Yieldratio[ibin_pt]->GetXaxis()->SetTitle("Pion p_{T}");
		h_PromptDs_PionPtScan_Yieldratio[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_PionPtScan_Yieldratio[ibin_pt]->SetMinimum(0.0);
		h_PromptDs_PionPtScan_Yieldratio[ibin_pt]->SetLineColor(2);
		h_PromptDs_PionPtScan_Yieldratio[ibin_pt]->SetMarkerColor(2);
		h_PromptDs_PionPtScan_Yieldratio[ibin_pt]->Draw("SAME");
		h_PromptDs_PionPtScan_Effratio[ibin_pt]->SetLineColor(4);
		h_PromptDs_PionPtScan_Effratio[ibin_pt]->SetMarkerColor(4);
		h_PromptDs_PionPtScan_Effratio[ibin_pt]->Draw("SAME");
		h_All_PionPtScan_Yieldratio[ibin_pt]->SetLineColor(1);
		h_All_PionPtScan_Yieldratio[ibin_pt]->SetMarkerColor(1);
		h_All_PionPtScan_Yieldratio[ibin_pt]->Draw("SAME");
		h_NonPromptDs_PionPtScan_Effratio[ibin_pt]->SetLineColor(kMagenta);
		h_NonPromptDs_PionPtScan_Effratio[ibin_pt]->SetMarkerColor(kMagenta);
		h_NonPromptDs_PionPtScan_Effratio[ibin_pt]->Draw("SAME");


		TLegend *le_PionPtScan_SingleRatio=new TLegend(0.7,0.7,0.85,0.85);
		le_PionPtScan_SingleRatio->SetBorderSize(0);	
		le_PionPtScan_SingleRatio->AddEntry(h_All_PionPtScan_Yieldratio[ibin_pt],"Data","lp");
		le_PionPtScan_SingleRatio->AddEntry(h_PromptDs_PionPtScan_Yieldratio[ibin_pt],"Data Prompt","lp");
		le_PionPtScan_SingleRatio->AddEntry(h_PromptDs_PionPtScan_Effratio[ibin_pt],"Prompt MC","lp");
		le_PionPtScan_SingleRatio->AddEntry(h_NonPromptDs_PionPtScan_Effratio[ibin_pt],"NonPrompt MC","lp");
		le_PionPtScan_SingleRatio->Draw("same");

		shiftY=0;
	  tlt->DrawLatexNDC(textposx+0.0,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;

		TLine *tl_PionPt_SigleRatio=new TLine(cut_DauPt,RatioMin,cut_DauPt,RatioMax);
		tl_PionPt_SigleRatio->SetLineColor(4);
		tl_PionPt_SigleRatio->SetLineStyle(6);
		tl_PionPt_SigleRatio->Draw("same");
	

		SavePlotDirs(cav1[c_count],Form("%s_Sys_PionPtScan_pt%.0fto%.0f_sigleRatio",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );


		// draw
		cav[c_count]= new TCanvas(Form("cav%i",c_count),Form("cav%i",c_count),800,800);
		cav[c_count]->cd();

		h_PromptDs_PionPtScan_Stat[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_PionPtScan_Stat[ibin_pt]->SetMinimum(RatioMin);
		h_PromptDs_PionPtScan_Stat[ibin_pt]->GetXaxis()->SetTitle("Pion p_{T}");
		h_PromptDs_PionPtScan_Stat[ibin_pt]->GetYaxis()->SetTitle("#sigma_{Varied cut}/#sigma_{Default cut}");
		h_PromptDs_PionPtScan_Stat[ibin_pt]->SetTitle("");
	
		h_PromptDs_PionPtScan_Stat[ibin_pt]->SetFillColor(kGray);
		h_PromptDs_PionPtScan_Stat[ibin_pt]->SetFillStyle(3001);
		h_PromptDs_PionPtScan_Stat[ibin_pt]->SetMarkerColor(kGray);
		h_PromptDs_PionPtScan_Stat[ibin_pt]->SetMarkerSize(0);
		h_PromptDs_PionPtScan_Stat[ibin_pt]->SetLineWidth(0);
		h_PromptDs_PionPtScan_Stat[ibin_pt]->SetLineColor(kGray);
		h_PromptDs_PionPtScan_Stat[ibin_pt]->Draw("SAME E2");


		h_PromptDs_PionPtScan_ratio[ibin_pt]->GetXaxis()->SetRangeUser(0,0.6);
		h_PromptDs_PionPtScan_ratio[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_PionPtScan_ratio[ibin_pt]->SetMinimum(RatioMin);
		h_PromptDs_PionPtScan_ratio[ibin_pt]->GetXaxis()->SetTitle("Pion p_{T}");
		h_PromptDs_PionPtScan_ratio[ibin_pt]->SetTitle("");
//		h_PromptDs_PionPtScan_ratio[ibin_pt]->GetYaxis()->SetTitle("");
		h_PromptDs_PionPtScan_ratio[ibin_pt]->Draw("SAME");
		f1_pol1->SetRange(0,3.75);
		f1_pol1->Draw("same");

		shiftY=0;
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;
		if(isPbPb==3){
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,"centrality : 0-100%"); shiftY-=oneshift;
		}
	  // tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%%",(abs(f1_pol1->GetParameter(0)-1))*100)); shiftY-=oneshift;
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%%",(abs(f1_pol1->Eval(0)-1))*100)); shiftY-=oneshift;

		tl_PionPt[ibin_pt]=new TLine(cut_DauPt,RatioMin,cut_DauPt,RatioMax);
		tl_PionPt[ibin_pt]->SetLineColor(4);
		tl_PionPt[ibin_pt]->SetLineStyle(6);
		tl_PionPt[ibin_pt]->Draw("same");
	

		SavePlotDirs(cav[c_count],Form("%s_Sys_PionPtScan_pt%.0fto%.0f",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );

		c_count++;

		// special handle

		// h_PromptDs_PionPtScan_SysRel->SetBinContent(ibin_pt+1, abs(f1_pol1->GetParameter(0)-1));
		// h_PromptDs_PionPtScan_SysRel->SetBinContent(ibin_pt+1, abs(f1_pol1->Eval(0)-1));

	cout<<"end PionPtScan Systematics "<<endl;
		//-- end PionPtScan_bins --//




		//// --  start of AllDauPt /// 
	cout<<"start AllDauPtScan Systematics "<<endl;

		// h_PromptDs_AllDauPtScan[ibin_pt]=(TH1D*)fin->Get(Form("h_PromptDs_AllDauPtScan_pt%.0fto%.0f",DptLow,DptHigh));
		h_PromptDs_AllDauPtScan[ibin_pt]=(TH1D*)fin->Get(Form("h_PromptDs_AllDauPtScan2_pt%.0fto%.0f",DptLow,DptHigh));
		cout<<"read AllDauPtScan TH1D "<<endl;
		h_PromptDs_AllDauPtScan_ratio[ibin_pt]=new TH1D(Form("h_PromptDs_AllDauPtScan_ratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_AllDauPtScan_ratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );

		h_PromptDs_AllDauPtScan_Yieldratio[ibin_pt]=new TH1D(Form("h_PromptDs_AllDauPtScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_AllDauPtScan_Yiledratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );
		h_PromptDs_AllDauPtScan_Effratio[ibin_pt]=new TH1D(Form("h_PromptDs_AllDauPtScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_AllDauPtScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );
		h_NonPromptDs_AllDauPtScan_Effratio[ibin_pt]=new TH1D(Form("h_NonPromptDs_AllDauPtScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_NonPromptDs_AllDauPtScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );
    Gr_PromptDs_AllDauPtScan_ratio[ibin_pt]=new TGraphErrors();

		h_PromptDs_AllDauPtScan_MCratio[ibin_pt]=new TH1D(Form("h_PromptDs_AllDauPtScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_AllDauPtScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );
		h_NonPromptDs_AllDauPtScan_MCratio[ibin_pt]=new TH1D(Form("h_NonPromptDs_AllDauPtScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_NonPromptDs_AllDauPtScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );
		h_All_AllDauPtScan_Yieldratio[ibin_pt]=new TH1D(Form("h_All_AllDauPtScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_All_AllDauPtScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );

		// h_PromptDs_AllDauPtScan_Stat[ibin_pt]=new TH1D(Form("h_PromptDs_AllDauPtScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_AllDauPtScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),nbin_DauPtScan, bins_DauPtScan );
		h_PromptDs_AllDauPtScan_Stat[ibin_pt]=new TH1D(Form("h_PromptDs_AllDauPtScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_AllDauPtScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),1,0, bins_DauPtScan[nbin_DauPtScan] );


	cout<<__LINE__<<endl;
		// TH1D *h_eff_AllDauPtScan[nbin_DauPtScan];
		
		for(int ibin_var=start_var; ibin_var<nbin_DauPtScan; ibin_var++){

			double CSdNdpt_var=h_PromptDs_AllDauPtScan[ibin_pt]->GetBinContent(ibin_var+1);
			double CSdNdptErr_var=h_PromptDs_AllDauPtScan[ibin_pt]->GetBinError(ibin_var+1);
			double CSdNdptErrRel_var=CSdNdptErr_var/CSdNdpt_var;	

			double cut_AllDauPtScan=DauPtScan_bins[ibin_var];

			double ratio=CSdNdpt_var/CSdNdpt;
			double ratioErrRel_temp = CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var - 2/CSdNdpt/CSdNdpt_var*CSdNdptErr*CSdNdptErr;
			// cout<<"ratioErrRel_temp = "<<ratioErrRel_temp<<endl;
	// for single ratio 

		double Eff_Prompt_phi_var=h_Eff_Prompt_phikkpi_AllDauPtScan[ibin_var]->GetBinContent(ibin_pt+1);
		double Eff_Prompt_f0_var=h_Eff_Prompt_f0kkpi_AllDauPtScan[ibin_var]->GetBinContent(ibin_pt+1);
		double Eff_NonPrompt_phi_var=h_Eff_NonPrompt_phikkpi_AllDauPtScan[ibin_var]->GetBinContent(ibin_pt+1);
		double Eff_NonPrompt_f0_var=h_Eff_NonPrompt_f0kkpi_AllDauPtScan[ibin_var]->GetBinContent(ibin_pt+1);
	
		double Eff_Prompt_AllBR_var=h_Eff_Prompt_AllBR_AllDauPtScan[ibin_var]->GetBinContent(ibin_pt+1);
		double EffErr_Prompt_AllBR_var=h_Eff_Prompt_AllBR_AllDauPtScan[ibin_var]->GetBinError(ibin_pt+1);
		double EffErrRel_Prompt_AllBR_var=EffErr_Prompt_AllBR_var/Eff_Prompt_AllBR_var;

		double Eff_NonPrompt_AllBR_var=h_Eff_NonPrompt_AllBR_AllDauPtScan[ibin_var]->GetBinContent(ibin_pt+1);
		double EffErr_NonPrompt_AllBR_var=h_Eff_NonPrompt_AllBR_AllDauPtScan[ibin_var]->GetBinError(ibin_pt+1);
		double EffErrRel_NonPrompt_AllBR_var=EffErr_NonPrompt_AllBR_var/Eff_NonPrompt_AllBR_var;

	cout<<__LINE__<<endl;

		double Prompt_yield_var=CSdNdpt_var*2*LumiNevt*(BRphi*Eff_Prompt_phi_var+BRf0*Eff_Prompt_f0_var);
		double Prompt_yieldErr_var=Prompt_yield_var*CSdNdptErrRel_var;

		double ratio_yield=Prompt_yield_var/Prompt_yield;

		double ratio_Eff=Eff_Prompt_AllBR_var/Eff_Prompt_AllBR;
		double ratioErr_Eff=ratio_Eff*( sqrt(EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR+ EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var)  );


		double ratio_NPEff=Eff_NonPrompt_AllBR_var/Eff_NonPrompt_AllBR;
		double ratioErr_NPEff=ratio_NPEff*( sqrt(EffErrRel_NonPrompt_AllBR*EffErrRel_NonPrompt_AllBR+ EffErrRel_NonPrompt_AllBR_var*EffErrRel_NonPrompt_AllBR_var)  );


		double All_yield_var=h_RawFitYield_AllDauPtScan[ibin_var]->GetBinContent(ibin_pt+1);
		double ratio_All_yield=All_yield_var/All_yield;
		double ratio_All_yieldErr=0;

	cout<<__LINE__<<endl;
		double ratioErrRel_temp_Eff=EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR + EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var-2/Eff_Prompt_AllBR/EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR*EffErr_Prompt_AllBR;


//		double ratioErr_yield=ratio_yield*
			
	// end single ratio


	cout<<__LINE__<<endl;

			// double ratioErr_Eff=ratio_Eff*sqrt(ratioErrRel_temp_Eff);
			double ratioErr = sqrt(abs(CSdNdptErr_var*CSdNdptErr_var-CSdNdptErr*CSdNdptErr))/CSdNdpt;	
			double ratioErr_yield=ratio_yield*sqrt(ratioErrRel_temp);
			ratio_All_yieldErr=ratio_All_yield*sqrt(ratioErrRel_temp);

	//// sigle ratio use full error //
			ratioErr_yield=Prompt_yieldErr_var/Prompt_yield;
			ratio_All_yieldErr=h_RawFitYield_AllDauPtScan[ibin_var]->GetBinError(ibin_pt+1)/All_yield;
			ratioErr_Eff=EffErr_Prompt_AllBR_var/Eff_Prompt_AllBR;
			ratioErr_NPEff=EffErr_NonPrompt_AllBR_var/Eff_NonPrompt_AllBR;

			cout<<"ratioErr_Eff = "<<ratioErr_Eff<<endl;

			// special case for var_cut = default cut
			if(cut_DauPt == cut_AllDauPtScan){
				ratio=1;
				ratioErr=0;
			}

			cout<<"\n ibin_var = "<<ibin_var<<" , AllDauPtScan var = "<<cut_AllDauPtScan<<endl;
			cout<<"CSdNdpt_var = "<<CSdNdpt_var<<" , CSdNdptErr_var = "<<CSdNdptErr_var<<" , CSdNdptErrRel_var = "<<CSdNdptErrRel_var*100<<"%"<<endl;
			cout<<"ratio = "<<ratio<<" +- "<<ratioErr<<" relErr = "<<ratioErr/ratio*100<<"%"<<endl;

			h_PromptDs_AllDauPtScan_ratio[ibin_pt]->SetBinContent(ibin_var+1, ratio);	
			h_PromptDs_AllDauPtScan_ratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr);	

			h_PromptDs_AllDauPtScan_Yieldratio[ibin_pt]->SetBinContent(ibin_var+1 , ratio_yield);
			h_PromptDs_AllDauPtScan_Yieldratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr_yield);
			h_PromptDs_AllDauPtScan_Effratio[ibin_pt]->SetBinContent(ibin_var+1, ratio_Eff);
			h_PromptDs_AllDauPtScan_Effratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr_Eff);

			h_NonPromptDs_AllDauPtScan_Effratio[ibin_pt]->SetBinContent(ibin_var+1, ratio_NPEff);
			h_NonPromptDs_AllDauPtScan_Effratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr_NPEff);
	
			h_All_AllDauPtScan_Yieldratio[ibin_pt]->SetBinContent(ibin_var+1, ratio_All_yield);
			h_All_AllDauPtScan_Yieldratio[ibin_pt]->SetBinError(ibin_var+1, ratio_All_yieldErr);

      if(cut_DauPt != cut_AllDauPtScan){
      Gr_PromptDs_AllDauPtScan_ratio[ibin_pt]->SetPoint(ibin_var,DauPtScan_bins[ibin_var],ratio);
      Gr_PromptDs_AllDauPtScan_ratio[ibin_pt]->SetPointError(ibin_var,0,ratioErr);
      }


		} // end for ibin_var<nbin_DauPtScan


			h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetBinContent(1, 1);	
			h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetBinError(1, CSdNdptErrRel);	

		fout->cd();
		h_PromptDs_AllDauPtScan_ratio[ibin_pt]->Write("",TObject::kOverwrite);
		Gr_PromptDs_AllDauPtScan_ratio[ibin_pt]->Write("",TObject::kOverwrite);
		
		// TF1 *f1_pol1=new TF1("f1_pol1","[0]+x*[1]");
		f1_pol1->SetParameter(0,1);
		// f1_pol1->FixParameter(0,1);
		f1_pol1->SetParameter(1,0);
		f1_pol1->SetParameter(2,cut_DauPt);
		f1_pol1->FixParameter(2,cut_DauPt);

		h_PromptDs_AllDauPtScan_ratio[ibin_pt]->Fit("f1_pol1","QN0 W");
		h_PromptDs_AllDauPtScan_ratio[ibin_pt]->Fit("f1_pol1","QN0 W");
		h_PromptDs_AllDauPtScan_ratio[ibin_pt]->Fit("f1_pol1","EMIS0 W");
		h_PromptDs_AllDauPtScan_ratio[ibin_pt]->Fit("f1_pol1","EMIS0 W");


    cav[c_count]= new TCanvas(Form("cav%i",c_count),Form("cav%i",c_count),800,800);
    cav[c_count]->cd();
    h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetMaximum(RatioMax);
    h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetMinimum(RatioMin);
    h_PromptDs_AllDauPtScan_Stat[ibin_pt]->GetXaxis()->SetTitle("AllDau p_{T}");
    h_PromptDs_AllDauPtScan_Stat[ibin_pt]->GetYaxis()->SetTitle("#sigma_{Varied cut}/#sigma_{Default cut}");
    h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetTitle("");

    h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetFillColor(kGray);
    h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetFillStyle(3001);
    h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetMarkerColor(kGray);
    h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetMarkerSize(0);
    h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetLineWidth(0);
    h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetLineColor(kGray);
    h_PromptDs_AllDauPtScan_Stat[ibin_pt]->Draw("SAME E2");
    Gr_PromptDs_AllDauPtScan_ratio[ibin_pt]->SetMarkerSize(0.5);
    Gr_PromptDs_AllDauPtScan_ratio[ibin_pt]->SetMarkerStyle(21);
    Gr_PromptDs_AllDauPtScan_ratio[ibin_pt]->Draw("P");

    f1_pol1->SetRange(0,3.75);
    f1_pol1->Draw("same");

    shiftY=0;
    tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;
    if(isPbPb==3){
    tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,"centrality : 0-100%"); shiftY-=oneshift;
    }
    tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%%",(abs(f1_pol1->GetParameter(0)-1))*100)); shiftY-=oneshift;

    tl_AllDauPt[ibin_pt]=new TLine(cut_DauPt,RatioMin,cut_DauPt,RatioMax);
    tl_AllDauPt[ibin_pt]->SetLineColor(4);
    tl_AllDauPt[ibin_pt]->SetLineStyle(6);
    tl_AllDauPt[ibin_pt]->Draw("same");


    // SavePlotDirs(cav[c_count],Form("%s_Sys_AllDauPtScan_Graph_pt%.0fto%.0f",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );

    c_count++;



		cav1[c_count]= new TCanvas(Form("cav1_%i",c_count),Form("cav1_%i",c_count),800,800);
		cav1[c_count]->cd();
		h_PromptDs_AllDauPtScan_Yieldratio[ibin_pt]->GetXaxis()->SetRangeUser(0,0.5);
		h_PromptDs_AllDauPtScan_Yieldratio[ibin_pt]->GetYaxis()->SetTitle("N_{vary cuts}/N_{default cuts}");
		h_PromptDs_AllDauPtScan_Yieldratio[ibin_pt]->GetXaxis()->SetTitle("AllDau p_{T}");
		h_PromptDs_AllDauPtScan_Yieldratio[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_AllDauPtScan_Yieldratio[ibin_pt]->SetMinimum(0.0);
		h_PromptDs_AllDauPtScan_Yieldratio[ibin_pt]->SetLineColor(2);
		h_PromptDs_AllDauPtScan_Yieldratio[ibin_pt]->SetMarkerColor(2);
		h_PromptDs_AllDauPtScan_Yieldratio[ibin_pt]->Draw("SAME");
		h_PromptDs_AllDauPtScan_Effratio[ibin_pt]->SetLineColor(4);
		h_PromptDs_AllDauPtScan_Effratio[ibin_pt]->SetMarkerColor(4);
		h_PromptDs_AllDauPtScan_Effratio[ibin_pt]->Draw("SAME");
		h_All_AllDauPtScan_Yieldratio[ibin_pt]->SetLineColor(1);
		h_All_AllDauPtScan_Yieldratio[ibin_pt]->SetMarkerColor(1);
		h_All_AllDauPtScan_Yieldratio[ibin_pt]->Draw("SAME");
		h_NonPromptDs_AllDauPtScan_Effratio[ibin_pt]->SetLineColor(kMagenta);
		h_NonPromptDs_AllDauPtScan_Effratio[ibin_pt]->SetMarkerColor(kMagenta);
		h_NonPromptDs_AllDauPtScan_Effratio[ibin_pt]->Draw("SAME");


		TLegend *le_AllDauPtScan_SingleRatio=new TLegend(0.7,0.7,0.85,0.85);
		le_AllDauPtScan_SingleRatio->SetBorderSize(0);	
		le_AllDauPtScan_SingleRatio->AddEntry(h_All_AllDauPtScan_Yieldratio[ibin_pt],"Data","lp");
		le_AllDauPtScan_SingleRatio->AddEntry(h_PromptDs_AllDauPtScan_Yieldratio[ibin_pt],"Data Prompt","lp");
		le_AllDauPtScan_SingleRatio->AddEntry(h_PromptDs_AllDauPtScan_Effratio[ibin_pt],"Prompt MC","lp");
		le_AllDauPtScan_SingleRatio->AddEntry(h_NonPromptDs_AllDauPtScan_Effratio[ibin_pt],"NonPrompt MC","lp");
		le_AllDauPtScan_SingleRatio->Draw("same");

		shiftY=0;
	  tlt->DrawLatexNDC(textposx+0.0,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;

		TLine *tl_AllDauPt_SigleRatio=new TLine(cut_DauPt,RatioMin,cut_DauPt,RatioMax);
		tl_AllDauPt_SigleRatio->SetLineColor(4);
		tl_AllDauPt_SigleRatio->SetLineStyle(6);
		tl_AllDauPt_SigleRatio->Draw("same");
	

		SavePlotDirs(cav1[c_count],Form("%s_Sys_AllDauPtScan_pt%.0fto%.0f_sigleRatio",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );


		// draw
		cav[c_count]= new TCanvas(Form("cav%i",c_count),Form("cav%i",c_count),800,800);
		cav[c_count]->cd();

		h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetMinimum(RatioMin);
		h_PromptDs_AllDauPtScan_Stat[ibin_pt]->GetXaxis()->SetTitle("AllDau p_{T}");
		h_PromptDs_AllDauPtScan_Stat[ibin_pt]->GetYaxis()->SetTitle("#sigma_{Varied cut}/#sigma_{Default cut}");
		h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetTitle("");
	
		h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetFillColor(kGray);
		h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetFillStyle(3001);
		h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetMarkerColor(kGray);
		h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetMarkerSize(0);
		h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetLineWidth(0);
		h_PromptDs_AllDauPtScan_Stat[ibin_pt]->SetLineColor(kGray);
		h_PromptDs_AllDauPtScan_Stat[ibin_pt]->Draw("SAME E2");


		h_PromptDs_AllDauPtScan_ratio[ibin_pt]->GetXaxis()->SetRangeUser(0,0.6);
		h_PromptDs_AllDauPtScan_ratio[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_AllDauPtScan_ratio[ibin_pt]->SetMinimum(RatioMin);
		h_PromptDs_AllDauPtScan_ratio[ibin_pt]->GetXaxis()->SetTitle("AllDau p_{T}");
		h_PromptDs_AllDauPtScan_ratio[ibin_pt]->SetTitle("");
//		h_PromptDs_AllDauPtScan_ratio[ibin_pt]->GetYaxis()->SetTitle("");
		h_PromptDs_AllDauPtScan_ratio[ibin_pt]->Draw("SAME");
		f1_pol1->SetRange(0,3.75);
		f1_pol1->Draw("same");

		shiftY=0;
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;
		if(isPbPb==3){
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,"centrality : 0-100%"); shiftY-=oneshift;
		}
	  // tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%%",(abs(f1_pol1->GetParameter(0)-1))*100)); shiftY-=oneshift;
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%%",(abs(f1_pol1->Eval(0)-1))*100)); shiftY-=oneshift;

		tl_AllDauPt[ibin_pt]=new TLine(cut_DauPt,RatioMin,cut_DauPt,RatioMax);
		tl_AllDauPt[ibin_pt]->SetLineColor(4);
		tl_AllDauPt[ibin_pt]->SetLineStyle(6);
		tl_AllDauPt[ibin_pt]->Draw("same");
	

		SavePlotDirs(cav[c_count],Form("%s_Sys_AllDauPtScan_pt%.0fto%.0f",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );

		c_count++;

		// special handle

		// h_PromptDs_AllDauPtScan_SysRel->SetBinContent(ibin_pt+1, abs(f1_pol1->GetParameter(0)-1));
		// h_PromptDs_AllDauPtScan_SysRel->SetBinContent(ibin_pt+1, abs(f1_pol1->Eval(0)-1));

	cout<<"end AllDauPtScan Systematics "<<endl;
		//-- end AllDauPtScan_bins --//


		} // end doDauPtCut


/*
	//-- start Reschi2clScan --//
	cout<<"start Reschi2clScan Systematics "<<endl;

		h_PromptDs_Reschi2clScan[ibin_pt]=(TH1D*)fin->Get(Form("h_PromptDs_Reschi2clScan_pt%.0fto%.0f",DptLow,DptHigh));
		// h_PromptDs_Reschi2clScan[ibin_pt]=(TH1D*)fin->Get(Form("h_PromptDs_Reschi2clScan2_pt%.0fto%.0f",DptLow,DptHigh));
		cout<<"read Reschi2clScan TH1D "<<endl;
		h_PromptDs_Reschi2clScan_ratio[ibin_pt]=new TH1D(Form("h_PromptDs_Reschi2clScan_ratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_Reschi2clScan_ratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_Reschi2clScan, bins_Reschi2clScan );

		h_PromptDs_Reschi2clScan_Yieldratio[ibin_pt]=new TH1D(Form("h_PromptDs_Reschi2clScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_Reschi2clScan_Yiledratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_Reschi2clScan, bins_Reschi2clScan );
		h_PromptDs_Reschi2clScan_Effratio[ibin_pt]=new TH1D(Form("h_PromptDs_Reschi2clScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_Reschi2clScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_Reschi2clScan, bins_Reschi2clScan );
		h_NonPromptDs_Reschi2clScan_Effratio[ibin_pt]=new TH1D(Form("h_NonPromptDs_Reschi2clScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_NonPromptDs_Reschi2clScan_Effratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_Reschi2clScan, bins_Reschi2clScan );


		h_PromptDs_Reschi2clScan_MCratio[ibin_pt]=new TH1D(Form("h_PromptDs_Reschi2clScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_Reschi2clScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_Reschi2clScan, bins_Reschi2clScan );
		h_NonPromptDs_Reschi2clScan_MCratio[ibin_pt]=new TH1D(Form("h_NonPromptDs_Reschi2clScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_NonPromptDs_Reschi2clScan_MCratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_Reschi2clScan, bins_Reschi2clScan );
		h_All_Reschi2clScan_Yieldratio[ibin_pt]=new TH1D(Form("h_All_Reschi2clScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),Form("h_All_Reschi2clScan_Yieldratio_pt%.0fto%.0f",DptLow,DptHigh),nbin_Reschi2clScan, bins_Reschi2clScan );


		// h_PromptDs_Reschi2clScan_Stat[ibin_pt]=new TH1D(Form("h_PromptDs_Reschi2clScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_Reschi2clScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),nbin_Reschi2clScan, bins_Reschi2clScan );
		h_PromptDs_Reschi2clScan_Stat[ibin_pt]=new TH1D(Form("h_PromptDs_Reschi2clScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),Form("h_PromptDs_Reschi2clScan_Stat_pt%.0fto%.0f",DptLow,DptHigh),1,0, bins_Reschi2clScan[nbin_Reschi2clScan] );

		cout<<"\n ibin_pt = "<<ibin_pt<<" , Ds pt = "<<DptLow<<" to "<<DptHigh<<endl;
		// cout<<" cut_Reschi2cl = "<<cut_Reschi2cl<<" , cut_Reschi2cl = "<<cut_Reschi2cl<<" , cut_Ddls = "<<cut_Ddls<<endl;
		cout<<" CSdNdpt = "<<CSdNdpt<<" CSdNdptErr = "<<CSdNdptErr<<" , CSdNdptErrRel = "<<CSdNdptErrRel*100<<"%"<<endl;
		cout<<" \n Reschi2clScan "<<endl;

		start_var=0;
		if(isPbPb==3 &&ibin_pt<=2){ start_var=0;} // don't use first bin result, fitting failed

		//// --  start of Reschi2cl 

		// TH1D *h_eff_Reschi2clScan[nbin_Reschi2clScan];
		

		for(int ibin_var=start_var; ibin_var<nbin_Reschi2clScan; ibin_var++){

			double CSdNdpt_var=h_PromptDs_Reschi2clScan[ibin_pt]->GetBinContent(ibin_var+1);
			double CSdNdptErr_var=h_PromptDs_Reschi2clScan[ibin_pt]->GetBinError(ibin_var+1);
			double CSdNdptErrRel_var=CSdNdptErr_var/CSdNdpt_var;	

			double cut_Reschi2clScan=Reschi2clScan_bins[ibin_var];

			double ratio=CSdNdpt_var/CSdNdpt;
			double ratioErrRel_temp = CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var - 2/CSdNdpt/CSdNdpt_var*CSdNdptErr*CSdNdptErr;
			// cout<<"ratioErrRel_temp = "<<ratioErrRel_temp<<endl;
	// for single ratio 

		double Eff_Prompt_phi_var=h_Eff_Prompt_phikkpi_Reschi2clScan[ibin_var]->GetBinContent(ibin_pt+1);
		double Eff_Prompt_f0_var=h_Eff_Prompt_f0kkpi_Reschi2clScan[ibin_var]->GetBinContent(ibin_pt+1);
		double Eff_NonPrompt_phi_var=h_Eff_NonPrompt_phikkpi_Reschi2clScan[ibin_var]->GetBinContent(ibin_pt+1);
		double Eff_NonPrompt_f0_var=h_Eff_NonPrompt_f0kkpi_Reschi2clScan[ibin_var]->GetBinContent(ibin_pt+1);
	
		double Eff_Prompt_AllBR_var=h_Eff_Prompt_AllBR_Reschi2clScan[ibin_var]->GetBinContent(ibin_pt+1);
		double EffErr_Prompt_AllBR_var=h_Eff_Prompt_AllBR_Reschi2clScan[ibin_var]->GetBinError(ibin_pt+1);
		double EffErrRel_Prompt_AllBR_var=EffErr_Prompt_AllBR_var/Eff_Prompt_AllBR_var;

		double Eff_NonPrompt_AllBR_var=h_Eff_NonPrompt_AllBR_Reschi2clScan[ibin_var]->GetBinContent(ibin_pt+1);
		double EffErr_NonPrompt_AllBR_var=h_Eff_NonPrompt_AllBR_Reschi2clScan[ibin_var]->GetBinError(ibin_pt+1);
		double EffErrRel_NonPrompt_AllBR_var=EffErr_NonPrompt_AllBR_var/Eff_NonPrompt_AllBR_var;


		double Prompt_yield_var=CSdNdpt_var*2*LumiNevt*(BRphi*Eff_Prompt_phi_var+BRf0*Eff_Prompt_f0_var);
		double Prompt_yieldErr_var=Prompt_yield_var*CSdNdptErrRel_var;

		double ratio_yield=Prompt_yield_var/Prompt_yield;

		double ratio_Eff=Eff_Prompt_AllBR_var/Eff_Prompt_AllBR;
		double ratioErr_Eff=ratio_Eff*( sqrt(EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR+ EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var)  );


		double ratio_NPEff=Eff_NonPrompt_AllBR_var/Eff_NonPrompt_AllBR;
		double ratioErr_NPEff=ratio_NPEff*( sqrt(EffErrRel_NonPrompt_AllBR*EffErrRel_NonPrompt_AllBR+ EffErrRel_NonPrompt_AllBR_var*EffErrRel_NonPrompt_AllBR_var)  );


		double All_yield_var=h_RawFitYield_Reschi2clScan[ibin_var]->GetBinContent(ibin_pt+1);
		double ratio_All_yield=All_yield_var/All_yield;
		double ratio_All_yieldErr=0;

		// cout<<"Eff_Prompt_AllBR_var = "<<Eff_Prompt_AllBR_var<<" , Eff_Prompt_AllBR = "<<Eff_Prompt_AllBR<<" , ratio_Eff = "<<ratio_Eff<<endl;
	// return 1;
		double ratioErrRel_temp_Eff=EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR + EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var-2/Eff_Prompt_AllBR/EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR*EffErr_Prompt_AllBR;



//		double ratioErr_yield=ratio_yield*


 // cout test
		// cout<<"ibin_pt = "<<ibin_pt<<" , ibin_var = "<<ibin_var<<endl;
		// cout<<"Prompt_yield = "<<Prompt_yield<<" , CSdNdpt = "<<CSdNdpt<<" , Eff_Prompt_phi = "<<Eff_Prompt_phi<<", Eff_Prompt_f0 = "<<Eff_Prompt_f0<<", Eff_Prompt_AllBR = "<< Eff_Prompt_AllBR<<endl;
		// cout<<"Prompt_yield_var = "<<Prompt_yield_var<<" , CSdNdpt_var = "<<CSdNdpt_var<<" , Eff_Prompt_phi_var = "<<Eff_Prompt_phi_var<<", Eff_Prompt_f0_var = "<<Eff_Prompt_f0_var<<", Eff_Prompt_AllBR_var = "<< Eff_Prompt_AllBR_var<<endl;

		// return 1;
			
	// end single ratio

			if(cut_Reschi2cl > cut_Reschi2clScan ){
			ratioErrRel_temp = CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var;
			ratioErrRel_temp_Eff=EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR + EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var-2/Eff_Prompt_AllBR/EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR_var;
			}
			if(ratioErrRel_temp <0 ){
			ratioErrRel_temp = CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var >0 ? CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var : CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var - 2/CSdNdpt/CSdNdpt_var*CSdNdptErr*CSdNdptErr;
      ratioErrRel_temp_Eff= EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR + EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var-2/Eff_Prompt_AllBR/EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR_var >0 ? EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR + EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var-2/Eff_Prompt_AllBR/EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR_var : EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR + EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var-2/Eff_Prompt_AllBR/EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR*EffErr_Prompt_AllBR;

			}
			if(ratioErrRel_temp <0 ){
				cout<<" --- warning , ratioErrRel_temp still < 0 !!! "<<endl;
			}
			if(ratioErrRel_temp_Eff <0 ){
				// cout<<"temp 1 = "<<EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR + EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var-2/Eff_Prompt_AllBR/EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR_var<<endl;
				// cout<<"temp 2 = "<<EffErrRel_Prompt_AllBR*EffErrRel_Prompt_AllBR + EffErrRel_Prompt_AllBR_var*EffErrRel_Prompt_AllBR_var-2/Eff_Prompt_AllBR/EffErr_Prompt_AllBR_var*EffErr_Prompt_AllBR*EffErr_Prompt_AllBR<<endl;
			}


			cout<<"ratioErrRel_temp = "<<ratioErrRel_temp<<endl;
			double ratioErr = ratio*sqrt(ratioErrRel_temp);
			// double ratioErr_Eff=ratio_Eff*sqrt(ratioErrRel_temp_Eff);
      ratioErr = sqrt(abs(CSdNdptErr_var*CSdNdptErr_var-CSdNdptErr*CSdNdptErr))/CSdNdpt;
			double ratioErr_yield=ratio_yield*sqrt(ratioErrRel_temp);
			ratio_All_yieldErr=ratio_All_yield*sqrt(ratioErrRel_temp);

			cout<<"ratioErr_Eff = "<<ratioErr_Eff<<endl;

			// special case for var_cut = default cut
			if(cut_Reschi2cl == cut_Reschi2clScan){
				ratio=1;
				ratioErr=0;
			}

			// double ratioErr=ratio*sqrt( CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr*CSdNdptErr  );
//			if(cut_Reschi2cl < cut_Reschi2clScan ){
	//			ratioErr=ratio*sqrt( CSdNdptErrRel*CSdNdptErrRel + CSdNdptErrRel_var*CSdNdptErrRel_var -2/CSdNdpt/CSdNdpt_var*CSdNdptErr_var*CSdNdptErr_var  );			
//			}

			cout<<"\n ibin_var = "<<ibin_var<<" , Reschi2clScan var = "<<cut_Reschi2clScan<<endl;
			cout<<"CSdNdpt_var = "<<CSdNdpt_var<<" , CSdNdptErr_var = "<<CSdNdptErr_var<<" , CSdNdptErrRel_var = "<<CSdNdptErrRel_var*100<<"%"<<endl;
			cout<<"ratio = "<<ratio<<" +- "<<ratioErr<<" relErr = "<<ratioErr/ratio*100<<"%"<<endl;

			h_PromptDs_Reschi2clScan_ratio[ibin_pt]->SetBinContent(ibin_var+1, ratio);	
			h_PromptDs_Reschi2clScan_ratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr);	

			h_PromptDs_Reschi2clScan_Yieldratio[ibin_pt]->SetBinContent(ibin_var+1 , ratio_yield);
			h_PromptDs_Reschi2clScan_Yieldratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr_yield);
			h_PromptDs_Reschi2clScan_Effratio[ibin_pt]->SetBinContent(ibin_var+1, ratio_Eff);
			h_PromptDs_Reschi2clScan_Effratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr_Eff);

			h_NonPromptDs_Reschi2clScan_Effratio[ibin_pt]->SetBinContent(ibin_var+1, ratio_NPEff);
			h_NonPromptDs_Reschi2clScan_Effratio[ibin_pt]->SetBinError(ibin_var+1, ratioErr_NPEff);
	
			h_All_Reschi2clScan_Yieldratio[ibin_pt]->SetBinContent(ibin_var+1, ratio_All_yield);
			h_All_Reschi2clScan_Yieldratio[ibin_pt]->SetBinError(ibin_var+1, ratio_All_yieldErr);

	
			// h_PromptDs_Reschi2clScan_Stat[ibin_pt]->SetBinContent(ibin_var+1, 1);	
			// h_PromptDs_Reschi2clScan_Stat[ibin_pt]->SetBinError(ibin_var+1, CSdNdptErrRel);	
		} // end for ibin_var<nbin_Reschi2clScan


			h_PromptDs_Reschi2clScan_Stat[ibin_pt]->SetBinContent(1, 1);	
			h_PromptDs_Reschi2clScan_Stat[ibin_pt]->SetBinError(1, CSdNdptErrRel);	

		fout->cd();
		h_PromptDs_Reschi2clScan_ratio[ibin_pt]->Write("",TObject::kOverwrite);
		h_PromptDs_Reschi2clScan_Stat[ibin_pt]->Write("",TObject::kOverwrite);		

		// fit			
		// TF1 *f1_pol1=new TF1("f1_pol1","[0]+x*[1]");
		f1_pol1->SetParameter(0,1);
		f1_pol1->SetParameter(1,0);
		h_PromptDs_Reschi2clScan_ratio[ibin_pt]->Fit("f1_pol1","QN0 W");
		h_PromptDs_Reschi2clScan_ratio[ibin_pt]->Fit("f1_pol1","QN0 W");
		h_PromptDs_Reschi2clScan_ratio[ibin_pt]->Fit("f1_pol1","EMIS0 W");
		h_PromptDs_Reschi2clScan_ratio[ibin_pt]->Fit("f1_pol1","EMIS0 W");

		cav1[c_count]= new TCanvas(Form("cav1_%i",c_count),Form("cav1_%i",c_count),800,800);
		cav1[c_count]->cd();
		h_PromptDs_Reschi2clScan_Yieldratio[ibin_pt]->GetXaxis()->SetRangeUser(0,0.6);
		h_PromptDs_Reschi2clScan_Yieldratio[ibin_pt]->GetYaxis()->SetTitle("N_{vary cuts}/N_{default cuts}");
		h_PromptDs_Reschi2clScan_Yieldratio[ibin_pt]->GetXaxis()->SetTitle("Res. vertex probability");
		h_PromptDs_Reschi2clScan_Yieldratio[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_Reschi2clScan_Yieldratio[ibin_pt]->SetMinimum(0.0);
		h_PromptDs_Reschi2clScan_Yieldratio[ibin_pt]->SetLineColor(2);
		h_PromptDs_Reschi2clScan_Yieldratio[ibin_pt]->SetMarkerColor(2);
		h_PromptDs_Reschi2clScan_Yieldratio[ibin_pt]->Draw("SAME");
		h_PromptDs_Reschi2clScan_Effratio[ibin_pt]->SetLineColor(4);
		h_PromptDs_Reschi2clScan_Effratio[ibin_pt]->SetMarkerColor(4);
		h_PromptDs_Reschi2clScan_Effratio[ibin_pt]->Draw("SAME");
		h_All_Reschi2clScan_Yieldratio[ibin_pt]->SetLineColor(1);
		h_All_Reschi2clScan_Yieldratio[ibin_pt]->SetMarkerColor(1);
		h_All_Reschi2clScan_Yieldratio[ibin_pt]->Draw("SAME");
		h_NonPromptDs_Reschi2clScan_Effratio[ibin_pt]->SetLineColor(kMagenta);
		h_NonPromptDs_Reschi2clScan_Effratio[ibin_pt]->SetMarkerColor(kMagenta);
		h_NonPromptDs_Reschi2clScan_Effratio[ibin_pt]->Draw("SAME");


		TLegend *le_Reschi2clScan_SingleRatio=new TLegend(0.7,0.7,0.85,0.85);
		le_Reschi2clScan_SingleRatio->SetBorderSize(0);	
		le_Reschi2clScan_SingleRatio->AddEntry(h_All_Reschi2clScan_Yieldratio[ibin_pt],"Data","lp");
		le_Reschi2clScan_SingleRatio->AddEntry(h_PromptDs_Reschi2clScan_Yieldratio[ibin_pt],"Data Prompt","lp");
		le_Reschi2clScan_SingleRatio->AddEntry(h_PromptDs_Reschi2clScan_Effratio[ibin_pt],"Prompt MC","lp");
		le_Reschi2clScan_SingleRatio->AddEntry(h_NonPromptDs_Reschi2clScan_Effratio[ibin_pt],"NonPrompt MC","lp");
		le_Reschi2clScan_SingleRatio->Draw("same");

		shiftY=0;
	  tlt->DrawLatexNDC(textposx+0.0,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;

		TLine *tl_Reschi2cl_SigleRatio=new TLine(cut_Reschi2cl,RatioMin,cut_Reschi2cl,RatioMax);
		tl_Reschi2cl_SigleRatio->SetLineColor(4);
		tl_Reschi2cl_SigleRatio->SetLineStyle(6);
		tl_Reschi2cl_SigleRatio->Draw("same");
	

		SavePlotDirs(cav1[c_count],Form("%s_Sys_Reschi2clScan_pt%.0fto%.0f_sigleRatio",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );


		// draw
		cav[c_count]= new TCanvas(Form("cav%i",c_count),Form("cav%i",c_count),800,800);
		cav[c_count]->cd();

		h_PromptDs_Reschi2clScan_Stat[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_Reschi2clScan_Stat[ibin_pt]->SetMinimum(RatioMin);
		h_PromptDs_Reschi2clScan_Stat[ibin_pt]->GetXaxis()->SetTitle("Res. vertex probability");
		h_PromptDs_Reschi2clScan_Stat[ibin_pt]->GetYaxis()->SetTitle("#sigma_{Varied cut}/#sigma_{Default cut}");
		h_PromptDs_Reschi2clScan_Stat[ibin_pt]->SetTitle("");
	
		h_PromptDs_Reschi2clScan_Stat[ibin_pt]->SetFillColor(kGray);
		h_PromptDs_Reschi2clScan_Stat[ibin_pt]->SetFillStyle(3001);
		h_PromptDs_Reschi2clScan_Stat[ibin_pt]->SetMarkerColor(kGray);
		h_PromptDs_Reschi2clScan_Stat[ibin_pt]->SetMarkerSize(0);
		h_PromptDs_Reschi2clScan_Stat[ibin_pt]->SetLineWidth(0);
		h_PromptDs_Reschi2clScan_Stat[ibin_pt]->SetLineColor(kGray);
		h_PromptDs_Reschi2clScan_Stat[ibin_pt]->Draw("SAME E2");


		h_PromptDs_Reschi2clScan_ratio[ibin_pt]->GetXaxis()->SetRangeUser(0,0.6);
		h_PromptDs_Reschi2clScan_ratio[ibin_pt]->SetMaximum(RatioMax);
		h_PromptDs_Reschi2clScan_ratio[ibin_pt]->SetMinimum(RatioMin);
		h_PromptDs_Reschi2clScan_ratio[ibin_pt]->GetXaxis()->SetTitle("D_{S} vertex probability");
		h_PromptDs_Reschi2clScan_ratio[ibin_pt]->SetTitle("");
//		h_PromptDs_Reschi2clScan_ratio[ibin_pt]->GetYaxis()->SetTitle("");
		h_PromptDs_Reschi2clScan_ratio[ibin_pt]->Draw("SAME");
		f1_pol1->SetRange(0,0.6);
		f1_pol1->Draw("same");

		shiftY=0;
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s : %.0f < p_{T} < %.0f GeV",str_PbPbtext.Data(),DptLow,DptHigh)); shiftY-=oneshift;
		if(isPbPb==3){
	  tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,"centrality : 0-100%"); shiftY-=oneshift;
		}
	  // tlt->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Systematics : %.1f%%",(abs(f1_pol1->GetParameter(0)-1))*100)); shiftY-=oneshift;

		tl_Reschi2cl[ibin_pt]=new TLine(cut_Reschi2cl,RatioMin,cut_Reschi2cl,RatioMax);
		tl_Reschi2cl[ibin_pt]->SetLineColor(4);
		tl_Reschi2cl[ibin_pt]->SetLineStyle(6);
		tl_Reschi2cl[ibin_pt]->Draw("same");
	

		SavePlotDirs(cav[c_count],Form("%s_Sys_Reschi2clScan_pt%.0fto%.0f",str_PbPb.Data(),DptLow,DptHigh),{"Systematics","CutScan_FixShape",Form("%s",str_PbPb.Data())} );

		c_count++;

		// special handle

		h_PromptDs_Reschi2clScan_SysRel->SetBinContent(ibin_pt+1, abs(f1_pol1->GetParameter(0)-1));

	cout<<"end Reschi2clScan Systematics "<<endl;
		//-- end Reschi2clScan_bins --//
*/



	}// end ibin_pt<nbin_pt


	fout->cd();
	h_PromptDs_DalphaMaxScan_SysRel->Write("",TObject::kOverwrite);
	h_PromptDs_Dchi2clMinScan_SysRel->Write("",TObject::kOverwrite);
	h_PromptDs_DdlsMinScan_SysRel->Write("",TObject::kOverwrite);
	h_PromptDs_PhiMassScan_SysRel->Write("",TObject::kOverwrite);
	// h_PromptDs_Reschi2clScan_SysRel->Write("",TObject::kOverwrite);

	// print out result
		cout<<"\n\n --- "<<str_PbPb<<endl;
	for(int ibin_pt=0; ibin_pt<nbin_pt; ibin_pt++){
		h_PromptDs_CutScanAll_SysRel->SetBinContent(ibin_pt+1, sqrt( h_PromptDs_DalphaMaxScan_SysRel->GetBinContent(ibin_pt+1)*h_PromptDs_DalphaMaxScan_SysRel->GetBinContent(ibin_pt+1)+ h_PromptDs_Dchi2clMinScan_SysRel->GetBinContent(ibin_pt+1)*h_PromptDs_Dchi2clMinScan_SysRel->GetBinContent(ibin_pt+1)+ h_PromptDs_DdlsMinScan_SysRel->GetBinContent(ibin_pt+1)*h_PromptDs_DdlsMinScan_SysRel->GetBinContent(ibin_pt+1) + h_PromptDs_PhiMassScan_SysRel->GetBinContent(ibin_pt+1)*h_PromptDs_PhiMassScan_SysRel->GetBinContent(ibin_pt+1) ) );

		cout<<"\n Ds pT : "<<bins_pt[ibin_pt]<<" - "<<bins_pt[ibin_pt+1]<<" GeV"<<endl;
		cout<<setw(18)<<"Dalpha Sys = "<<setw(10)<<setprecision(1)<<std::fixed<<h_PromptDs_DalphaMaxScan_SysRel->GetBinContent(ibin_pt+1)*100<<" %"<<endl;
		cout<<setw(18)<<"Dchi2cl Sys = "<<setw(10)<<setprecision(1)<<std::fixed<<h_PromptDs_Dchi2clMinScan_SysRel->GetBinContent(ibin_pt+1)*100<<" %"<<endl;
		cout<<setw(18)<<"Ddls Sys = "<<setw(10)<<setprecision(1)<<std::fixed<<h_PromptDs_DdlsMinScan_SysRel->GetBinContent(ibin_pt+1)*100<<" %"<<endl;
		cout<<setw(18)<<"PhiMass Sys = "<<setw(10)<<setprecision(1)<<std::fixed<<h_PromptDs_PhiMassScan_SysRel->GetBinContent(ibin_pt+1)*100<<" %"<<endl;
		
		cout<<setw(18)<<"Total Cut Sys = "<<setw(10)<<setprecision(1)<<std::fixed<<h_PromptDs_CutScanAll_SysRel->GetBinContent(ibin_pt+1)*100<<" %"<<endl;

	}

	cout<<"finish "<<str_PbPb<<endl;

	fout->cd();
	h_PromptDs_CutScanAll_SysRel->Write("",TObject::kOverwrite);

	return 0;
}
