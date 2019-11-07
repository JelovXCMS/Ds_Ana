#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"
#include <TGraphErrors.h>
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

#include <RooFit.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooHistPdf.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooGlobalFunc.h>
#include <RooGaussian.h>
#include <RooPolynomial.h>
#include <RooGenericPdf.h>
#include <RooFormulaVar.h>
#include <RooChebychev.h>
#include <RooExponential.h>
#include <RooFitResult.h>
#include <RooWorkspace.h>
#include <RooMCStudy.h>
#include "RooHist.h"
#include "RooConstVar.h"
#include "RooMsgService.h"


#include "RooStats/SPlot.h"
#include "RooConstVar.h"

#include "varCompare_para.h"

#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooAbsPdf.h"
#include <iomanip>

using namespace std;

int MergeWeight(){

	double nbin_pt_pp[]={2,4,6,8,10,40};
	double nbin_pt_PbPb[]={6,8,10,40};

	int nbins_pp=sizeof(nbin_pt_pp)/sizeof(nbin_pt_pp[0])-1;
	int nbins_PbPb=sizeof(nbin_pt_PbPb)/sizeof(nbin_pt_PbPb[0])-1;

	TFile *f_pp[nbins_pp];
	TFile *f_PbPb[nbins_PbPb];

	TFile *fout=TFile::Open("MergeWeight_Dalpha.root","recreate");

	TGraphErrors *gr_pp[nbins_pp];
	TH1D *h_MCP_WtScale_pp[nbins_pp];
	TH1D *h_MCNP_WtScale_pp[nbins_pp];

	TGraph *gr_pp_fit[nbins_pp];
	TH1D *h_MCP_FitWtScale_pp[nbins_pp];
	TH1D *h_MCNP_FitWtScale_pp[nbins_pp];

	for(int i=0; i<nbins_pp;i++){
		f_pp[i]=TFile::Open(Form("./Compare_DataMC_out/pp_Dpt%.0fto%.0f.root",nbin_pt_pp[i]*100,nbin_pt_pp[i+1]*100));
		gr_pp[i]=(TGraphErrors*)f_pp[i]->Get("gr_weight");
		gr_pp[i]->Draw();
		h_MCP_WtScale_pp[i]=(TH1D*)f_pp[i]->Get("h_MCP_WtScale");
		h_MCNP_WtScale_pp[i]=(TH1D*)f_pp[i]->Get("h_MCNP_WtScale");
	
		fout->cd();
		gr_pp[i]->SetName(Form("gr_weight_pp_Dpt%.0fto%.0f",nbin_pt_pp[i],nbin_pt_pp[i+1]));	
		gr_pp[i]->Write();
		h_MCP_WtScale_pp[i]->Write(Form("h_MCP_WtScale_pp_Dpt%.0fto%.0f",nbin_pt_pp[i],nbin_pt_pp[i+1]));
		h_MCNP_WtScale_pp[i]->Write(Form("h_MCNP_WtScale_pp_Dpt%.0fto%.0f",nbin_pt_pp[i],nbin_pt_pp[i+1]));


		gr_pp_fit[i]=(TGraphErrors*)f_pp[i]->Get("gr_Fitweight");
		gr_pp_fit[i]->Draw();
		h_MCP_FitWtScale_pp[i]=(TH1D*)f_pp[i]->Get("h_MCP_FitWtScale");
		h_MCNP_FitWtScale_pp[i]=(TH1D*)f_pp[i]->Get("h_MCNP_FitWtScale");
	
		gr_pp_fit[i]->SetName(Form("gr_Fitweight_pp_Dpt%.0fto%.0f",nbin_pt_pp[i],nbin_pt_pp[i+1]));	
		gr_pp_fit[i]->Write();
		h_MCP_FitWtScale_pp[i]->Write(Form("h_MCP_FitWtScale_pp_Dpt%.0fto%.0f",nbin_pt_pp[i],nbin_pt_pp[i+1]));
		h_MCNP_FitWtScale_pp[i]->Write(Form("h_MCNP_FitWtScale_pp_Dpt%.0fto%.0f",nbin_pt_pp[i],nbin_pt_pp[i+1]));


	
	}



	TGraphErrors *gr_PbPb[nbins_PbPb];
	TH1D *h_MCP_WtScale_PbPb[nbins_PbPb];
	TH1D *h_MCNP_WtScale_PbPb[nbins_PbPb];


	TGraph *gr_PbPb_fit[nbins_PbPb];
	TH1D *h_MCP_FitWtScale_PbPb[nbins_PbPb];
	TH1D *h_MCNP_FitWtScale_PbPb[nbins_PbPb];



	for(int i=0; i<nbins_PbPb;i++){
		f_PbPb[i]=TFile::Open(Form("./Compare_DataMC_out/PbPb_Dpt%.0fto%.0f.root",nbin_pt_PbPb[i]*100,nbin_pt_PbPb[i+1]*100));
		gr_PbPb[i]=(TGraphErrors*)f_PbPb[i]->Get("gr_weight");
		gr_PbPb[i]->Draw();
		h_MCP_WtScale_PbPb[i]=(TH1D*)f_PbPb[i]->Get("h_MCP_WtScale");
		h_MCNP_WtScale_PbPb[i]=(TH1D*)f_PbPb[i]->Get("h_MCNP_WtScale");
	
		fout->cd();
		gr_PbPb[i]->SetName(Form("gr_weight_PbPb_Dpt%.0fto%.0f",nbin_pt_PbPb[i],nbin_pt_PbPb[i+1]));	
		gr_PbPb[i]->Write();
		h_MCP_WtScale_PbPb[i]->Write(Form("h_MCP_WtScale_PbPb_Dpt%.0fto%.0f",nbin_pt_PbPb[i],nbin_pt_PbPb[i+1]));
		h_MCNP_WtScale_PbPb[i]->Write(Form("h_MCNP_WtScale_PbPb_Dpt%.0fto%.0f",nbin_pt_PbPb[i],nbin_pt_PbPb[i+1]));


		gr_PbPb_fit[i]=(TGraphErrors*)f_PbPb[i]->Get("gr_Fitweight");
		gr_PbPb_fit[i]->Draw();
		h_MCP_FitWtScale_PbPb[i]=(TH1D*)f_PbPb[i]->Get("h_MCP_FitWtScale");
		h_MCNP_FitWtScale_PbPb[i]=(TH1D*)f_PbPb[i]->Get("h_MCNP_FitWtScale");
	
		gr_PbPb_fit[i]->SetName(Form("gr_Fitweight_PbPb_Dpt%.0fto%.0f",nbin_pt_PbPb[i],nbin_pt_PbPb[i+1]));	
		gr_PbPb_fit[i]->Write();
		h_MCP_FitWtScale_PbPb[i]->Write(Form("h_MCP_FitWtScale_PbPb_Dpt%.0fto%.0f",nbin_pt_PbPb[i],nbin_pt_PbPb[i+1]));
		h_MCNP_FitWtScale_PbPb[i]->Write(Form("h_MCNP_FitWtScale_PbPb_Dpt%.0fto%.0f",nbin_pt_PbPb[i],nbin_pt_PbPb[i+1]));


	
	}


// #<{(|
// special treatment for PbPb 6-8 
  // TString mcName_Prompt="./rootF/ppMC_phiPrompt_fitFile.root";
  // TString mcName_NonPrompt="./rootF/ppMC_phiNonPrompt_fitFile.root";
  TString mcName_Prompt="../sPlot_varCompare_v2/rootF/PbPb3MC_phiPrompt_fitFile.root";
  TString mcName_NonPrompt="../sPlot_varCompare_v2/rootF/PbPb3MC_phiNonPrompt_fitFile.root";

  double Dalpha_cut=0.2;
  double Dchi2cl_cut=0.3;
  double Ddls_cut=5.0;
  double Dpt_Low=6;
  double Dpt_High=8;

	if(Dpt_High<=10){
		Dchi2cl_cut=0.3;
    Ddls_cut=5.0;
  }

  TFile *f_mc_Prompt=TFile::Open(mcName_Prompt.Data());
  TFile *f_mc_NonPrompt=TFile::Open(mcName_NonPrompt.Data());

  TTree *t_Ds_MCPrompt=(TTree*)f_mc_Prompt->Get(Form("t_fit"));
  TTree *t_Ds_MCNonPrompt=(TTree*)f_mc_NonPrompt->Get(Form("t_fit"));


  Float_t Dpt;
  Float_t Dalpha;
  Float_t Dchi2cl;
  Float_t Ddls;
  Float_t Dmass;
  Float_t TotalWeight;

	int nbin=200;
	double binLow=0.0;
	double binHigh=0.2;

  t_Ds_MCPrompt->SetBranchAddress("Dpt",&Dpt);
  t_Ds_MCPrompt->SetBranchAddress("Dalpha",&Dalpha);
  t_Ds_MCPrompt->SetBranchAddress("Dchi2cl",&Dchi2cl);
  t_Ds_MCPrompt->SetBranchAddress("Ddls",&Ddls);
  t_Ds_MCPrompt->SetBranchAddress("Dmass",&Dmass);
  t_Ds_MCPrompt->SetBranchAddress("TotalWeight",&TotalWeight);


  TH1D *h_MCP_Dalpha_Ori=new TH1D("h_MCP_Dalpha_Ori","h_MCP_Dalpha_Ori",nbin,binLow,binHigh); h_MCP_Dalpha_Ori->Sumw2();
  TH1D *h_MCP_Dalpha_Wet=new TH1D("h_MCP_Dalpha_Wet","h_MCP_Dalpha_Wet",nbin,binLow,binHigh); h_MCP_Dalpha_Wet->Sumw2();
  TH1D *h_MCP_Dalpha_FitWet=new TH1D("h_MCP_Dalpha_FitWet","h_MCP_Dalpha_FitWet",nbin,binLow,binHigh); h_MCP_Dalpha_FitWet->Sumw2();

  long int nEntries=t_Ds_MCPrompt->GetEntries();
  for(int i =0; i<nEntries; i++){
    t_Ds_MCPrompt->GetEntry(i);
    if(Dpt> Dpt_Low && Dpt < Dpt_High && Dalpha< Dalpha_cut && Dchi2cl > Dchi2cl_cut && Ddls> Ddls_cut && Dpt > Dpt_Low && Dpt < Dpt_High){
      h_MCP_Dalpha_Ori->Fill(Dalpha,TotalWeight);
      double alpha_weight= gr_PbPb[1]->Eval(Dalpha);
      if(alpha_weight<0) {alpha_weight=0;}
      h_MCP_Dalpha_Wet->Fill(Dalpha,TotalWeight*alpha_weight);
      h_MCP_Dalpha_FitWet->Fill(Dalpha,TotalWeight*gr_PbPb_fit[1]->Eval(Dalpha));

    }// end if cut

  } // end for  i<nEntries


    cout<<"Ori integral = "<<h_MCP_Dalpha_Ori->Integral()<<endl;
    cout<<"Wet integral = "<<h_MCP_Dalpha_Wet->Integral()<<endl;
    cout<<"Wet integral = "<<h_MCP_Dalpha_FitWet->Integral()<<endl;
  TH1D *h_MCP_WtScale=new TH1D("h_MCP_WtScale","h_MCP_WtScale",1,0,1);
  h_MCP_WtScale->SetBinContent(1,h_MCP_Dalpha_Ori->Integral()/h_MCP_Dalpha_Wet->Integral());
// |)}>#
  TH1D *h_MCP_FitWtScale=new TH1D("h_MCP_FitWtScale","h_MCP_FitWtScale",1,0,1);
  h_MCP_FitWtScale->SetBinContent(1,h_MCP_Dalpha_Ori->Integral()/h_MCP_Dalpha_FitWet->Integral());

	h_MCP_FitWtScale->Draw();

	fout->cd();
	h_MCP_FitWtScale->Write("h_MCP_FitWtScale_PbPb_Dpt6to8",TObject::kOverwrite);

/// nonprompt

  t_Ds_MCNonPrompt->SetBranchAddress("Dpt",&Dpt);
  t_Ds_MCNonPrompt->SetBranchAddress("Dalpha",&Dalpha);
  t_Ds_MCNonPrompt->SetBranchAddress("Dchi2cl",&Dchi2cl);
  t_Ds_MCNonPrompt->SetBranchAddress("Ddls",&Ddls);
  t_Ds_MCNonPrompt->SetBranchAddress("Dmass",&Dmass);
  t_Ds_MCNonPrompt->SetBranchAddress("TotalWeight",&TotalWeight);


  TH1D *h_MCNP_Dalpha_Ori=new TH1D("h_MCNP_Dalpha_Ori","h_MCNP_Dalpha_Ori",nbin,binLow,binHigh); h_MCNP_Dalpha_Ori->Sumw2();
  TH1D *h_MCNP_Dalpha_Wet=new TH1D("h_MCNP_Dalpha_Wet","h_MCNP_Dalpha_Wet",nbin,binLow,binHigh); h_MCNP_Dalpha_Wet->Sumw2();
  TH1D *h_MCNP_Dalpha_FitWet=new TH1D("h_MCNP_Dalpha_FitWet","h_MCNP_Dalpha_FitWet",nbin,binLow,binHigh); h_MCNP_Dalpha_FitWet->Sumw2();

  nEntries=t_Ds_MCNonPrompt->GetEntries();
  for(int i =0; i<nEntries; i++){
    t_Ds_MCNonPrompt->GetEntry(i);
    if(Dpt> Dpt_Low && Dpt < Dpt_High && Dalpha< Dalpha_cut && Dchi2cl > Dchi2cl_cut && Ddls> Ddls_cut && Dpt > Dpt_Low && Dpt < Dpt_High){
      h_MCNP_Dalpha_Ori->Fill(Dalpha,TotalWeight);
      double alpha_weight= gr_PbPb[1]->Eval(Dalpha);
      if(alpha_weight<0) {alpha_weight=0;}
      h_MCNP_Dalpha_Wet->Fill(Dalpha,TotalWeight*alpha_weight);
      h_MCNP_Dalpha_FitWet->Fill(Dalpha,TotalWeight*gr_PbPb_fit[1]->Eval(Dalpha));

    }// end if cut

  } // end for  i<nEntries


    cout<<"Ori integral = "<<h_MCNP_Dalpha_Ori->Integral()<<endl;
    cout<<"Wet integral = "<<h_MCNP_Dalpha_Wet->Integral()<<endl;
    cout<<"Wet integral = "<<h_MCNP_Dalpha_FitWet->Integral()<<endl;
  TH1D *h_MCNP_WtScale=new TH1D("h_MCNP_WtScale","h_MCNP_WtScale",1,0,1);
  h_MCNP_WtScale->SetBinContent(1,h_MCNP_Dalpha_Ori->Integral()/h_MCNP_Dalpha_Wet->Integral());
// |)}>#
  TH1D *h_MCNP_FitWtScale=new TH1D("h_MCNP_FitWtScale","h_MCNP_FitWtScale",1,0,1);
  h_MCNP_FitWtScale->SetBinContent(1,h_MCNP_Dalpha_Ori->Integral()/h_MCNP_Dalpha_FitWet->Integral());

	h_MCNP_FitWtScale->Draw();

	fout->cd();
	h_MCNP_FitWtScale->Write("h_MCNP_FitWtScale_PbPb_Dpt6to8",TObject::kOverwrite);


		gr_PbPb_fit[1]->SetName(Form("gr_Fitweight_PbPb_Dpt6to8"));	
		gr_PbPb_fit[1]->Write("",TObject::kOverwrite);

	return 0;


}
