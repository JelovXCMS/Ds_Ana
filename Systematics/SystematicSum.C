#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"


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

// using namespace RooFit;
using namespace std;


int SystematicSum(int isPbPb=0){

  double tex_upperY=0.95;

  TLatex* texCmsPre = new TLatex(0.15,tex_upperY, "#scale[1.25]{CMS} Preliminary");
  texCmsPre->SetNDC();
  texCmsPre->SetTextAlign(12);
  texCmsPre->SetTextSize(0.035);
  texCmsPre->SetTextFont(42);

  TLatex* texCmsSim = new TLatex(0.15,tex_upperY, "#scale[1.25]{CMS} Simulations");
  texCmsSim->SetNDC();
  texCmsSim->SetTextAlign(12);
  texCmsSim->SetTextSize(0.035);
  texCmsSim->SetTextFont(42);


  // TLatex* texColPbPb = new TLatex(0.88,tex_upperY, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  TLatex* texColPbPb = new TLatex(0.95,tex_upperY, "44 #mub^{-1} (5.02 TeV PbPb)");
  texColPbPb->SetNDC();
  texColPbPb->SetTextAlign(32);
  texColPbPb->SetTextSize(0.035);
  texColPbPb->SetTextFont(42);

  // TLatex* texColpp = new TLatex(0.88,tex_upperY, "pp #sqrt{s} = 5.02 TeV");
  TLatex* texColpp = new TLatex(0.95,tex_upperY, "38 nb^{-1} (5.02 TeV pp )");
  texColpp->SetNDC();
  texColpp->SetTextAlign(32);
  texColpp->SetTextSize(0.035);
  texColpp->SetTextFont(42);


	int startbin=0;

	  if(isPbPb==3){startbin=2;}
  InitStyle();
  initParameter();
  setTDRStyle();
	

  gStyle->SetOptStat(0);
  TString str_PbPb="pp";
  TString str_PbPbtext="pp";
  int nbin_pt=nbin_pt_pp;
  double *bins_pt=bins_pt_pp;


  if(isPbPb==3){
    str_PbPb="PbPb3";
    str_PbPbtext="PbPb";
    nbin_pt=nbin_pt_PbPb3;
		bins_pt=bins_pt_PbPb3;
	}

	cout<<"check 1"<<endl;

	TFile *fout=TFile::Open(Form("./output/SysSum_%s.root",str_PbPb.Data()),"RECREATE");
	fout->cd();

	TH1D *h_Tracking_SysRel=new TH1D("h_Tracking_SysRel","h_Tracking_SysRel",nbin_pt,bins_pt);
	TH1D *h_DsBR_SysRel=new TH1D("h_DsBR_SysRel","h_DsBR_SysRel",nbin_pt,bins_pt);
	TH1D *h_Lumi_SysRel=new TH1D("h_Lumi_SysRel","h_Lumi_SysRel",nbin_pt,bins_pt);
	TH1D *h_f0_SysRel=new TH1D("h_f0_SysRel","h_f0_SysRel",nbin_pt,bins_pt);
	TH1D *h_PhiRatio_SysRel=new TH1D("h_PhiRatio_SysRel","h_PhiRatio_SysRel",nbin_pt,bins_pt);
	TH1D *h_Total_SysRel=new TH1D("h_Total_SysRel","h_Total_SysRel",nbin_pt,bins_pt);


	TH1D *h_D0_SysRel=new TH1D("h_D0_SysRel","h_D0_SysRel",nbin_pt,bins_pt);
	TH1D *h_DsOverD0_SysRel=new TH1D("h_DsOverD0_SysRel","h_DsOverD0_SysRel",nbin_pt,bins_pt);
	TH1D *h_TrkOne_SysRel=new TH1D("h_TrkOne_SysRel","h_TrkOne_SysRel",nbin_pt,bins_pt);
	TH1D *h_partial_forRaa_SysRel=new TH1D("h_partial_forRaa_SysRel","h_partial_forRaa_SysRel",nbin_pt,bins_pt);

	cout<<"check 2"<<endl;

	TFile *fBtoDs=TFile::Open(Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/Systematics/output/BtoDsSys_%s.root",str_PbPb.Data()));
	TH1D *h_PromptDs_BtoDs_SysRel=(TH1D*)fBtoDs->Get("h_PromptDs_BtoDs_SysRel");
	TH1D *h_PromptDs_BtoDs_SysRel_ForRaa=(TH1D*)fBtoDs->Get("h_PromptDs_BtoDs_SysRel_ForRaa");
	TH1D *h_PromptDs_BtoDs_SysRel_ForRaaCancel=(TH1D*)fBtoDs->Get("h_PromptDs_BtoDs_SysRel_ForRaaCancel");

	cout<<"check 3"<<endl;
	// TFile *fCutScan=TFile::Open(Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/Systematics/output/CutScanSys_%s.root",str_PbPb.Data()));
	// TFile *fCutScan_old=TFile::Open(Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/Systematics/output/CutScanSys_FixShape_%s.root",str_PbPb.Data()));
	TFile *fCutScan=TFile::Open(Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/ARC_HW/FitCutScan/RatioOut_BestScale/CutScanSys_FixShape_%s.root",str_PbPb.Data()));
	TH1D *h_PromptDs_CutScanAll_SysRel=(TH1D*)fCutScan->Get("h_PromptDs_CutScanAll_SysRel");
	// TH1D *h_PromptDs_CutScanAll_SysRel=(TH1D*)fCutScan_old->Get("h_PromptDs_CutScanAll_SysRel");

	// this broken for unknow reason
	TH1D *h_DdlsScale_SysRel=(TH1D*)fCutScan->Get("h_DdlsScale_SysRel");
	if(!isPbPb){
		h_DdlsScale_SysRel->SetBinContent(1,0.049);
		h_DdlsScale_SysRel->SetBinContent(2,0.024);
		h_DdlsScale_SysRel->SetBinContent(3,0.018);
		h_DdlsScale_SysRel->SetBinContent(4,0.028);
		h_DdlsScale_SysRel->SetBinContent(5,0.015);
		h_DdlsScale_SysRel->SetBinContent(6,0.017);
		h_DdlsScale_SysRel->SetBinContent(7,0.011);
		h_DdlsScale_SysRel->SetBinContent(8,0.024);
	}else{
		h_DdlsScale_SysRel->SetBinContent(3,0.039);
		h_DdlsScale_SysRel->SetBinContent(4,0.069);
		h_DdlsScale_SysRel->SetBinContent(5,0.026);
		h_DdlsScale_SysRel->SetBinContent(6,0.012);
	}


	// h_PromptDs_CutScanAll_SysRel->Draw();	

	// return 1;

	cout<<"check 4"<<endl;
	TFile *fMCShape=TFile::Open(Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/Systematics/output/MCShapeSys_%s.root",str_PbPb.Data()));
	TH1D *h_PromptDs_f0Eff_SysRel=(TH1D*)fMCShape->Get("h_PromptDs_f0Eff_SysRel");
	TH1D *h_PromptDs_PhiRatio_SysRel=(TH1D*)fMCShape->Get("h_PromptDs_PhiRatio_SysRel");
	// TH1D *h_PromptDs_f0Eff_SysRel=new TH1D("h_PromptDs_f0Eff_SysRel","h_PromptDs_f0Eff_SysRel",nbin_pt,bins_pt);
	TH1D *h_PromptDs_MCShape_SysRel=(TH1D*)fMCShape->Get("h_PromptDs_MCShape_SysRel");

	cout<<"check 5"<<endl;
	TFile *fpdf=TFile::Open(Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/SignalFit/output%s/RawFitYield_FixShapeTrkPtScan_%s.root",s_CutSet.Data(),str_PbPb.Data()));
	TH1D *h_RawRooFitYield_pdfVar_RelErr=(TH1D*)fpdf->Get("h_RawRooFitYield_pdfVar_RelErr");

	cout<<"check 6"<<endl;

	if(isPbPb==0){
		h_D0_SysRel->SetBinContent(1,0.1438);
		h_D0_SysRel->SetBinContent(2,0.1347);
		h_D0_SysRel->SetBinContent(3,0.1144);
		h_D0_SysRel->SetBinContent(4,0.1130);
		h_D0_SysRel->SetBinContent(5,0.1116);
		h_D0_SysRel->SetBinContent(6,0.1089);
		h_D0_SysRel->SetBinContent(7,0.1147);
		h_D0_SysRel->SetBinContent(8,0.1095);
	}else{
		h_D0_SysRel->SetBinContent(3,0.1191);
		h_D0_SysRel->SetBinContent(4,0.1116);
		h_D0_SysRel->SetBinContent(5,0.1254);
		h_D0_SysRel->SetBinContent(6,0.1422);
	}


	double Tracking_Sys=0.12;
	double Tracking_one=0.04;
	if(isPbPb==3){
		Tracking_Sys=0.15;
		Tracking_one=0.05;
	}
	double Total_Sys=0;
	double DsOverD0_Sys=0;
	double partialRaa_Sys=0;

  // double PhiRatioErr_pp=0.014/0.954;
  // double PhiRatioErr_PbPb=0.009/0.896;

  double PhiRatioErr=PhiRatioErr_pp;
  if(isPbPb){
    PhiRatioErr=PhiRatioErr_PbPb;
  }


	for(int i=startbin; i<nbin_pt; i++){
		h_Tracking_SysRel->SetBinContent(i+1,Tracking_Sys);
		h_DsBR_SysRel->SetBinContent(i+1,0.035);
		h_Lumi_SysRel->SetBinContent(i+1,0.023);
		h_f0_SysRel->SetBinContent(i+1,0.13);
		h_PhiRatio_SysRel->SetBinContent(i+1,PhiRatioErr);
		
		// Total_Sys=sqrt( pow(h_Tracking_SysRel->GetBinContent(i+1),2) + pow(h_DsBR_SysRel->GetBinContent(i+1),2) + pow(h_f0_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_BtoDs_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_CutScanAll_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_MCShape_SysRel->GetBinContent(i+1),2) + pow(h_RawRooFitYield_pdfVar_RelErr->GetBinContent(i+1),2));
		// Total_Sys=sqrt( pow(h_Tracking_SysRel->GetBinContent(i+1),2) + pow(h_DsBR_SysRel->GetBinContent(i+1),2) + pow(h_PhiRatio_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_BtoDs_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_CutScanAll_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_MCShape_SysRel->GetBinContent(i+1),2) + pow(h_RawRooFitYield_pdfVar_RelErr->GetBinContent(i+1),2));
		// Total_Sys=sqrt( pow(h_Tracking_SysRel->GetBinContent(i+1),2) + pow(h_DsBR_SysRel->GetBinContent(i+1),2) + pow(h_PhiRatio_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_BtoDs_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_CutScanAll_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_MCShape_SysRel->GetBinContent(i+1),2) + pow(h_RawRooFitYield_pdfVar_RelErr->GetBinContent(i+1),2) + pow(h_DdlsScale_SysRel->GetBinContent(i+1),2));

		// remove PhiRatio_sys from sys 
		Total_Sys=sqrt( pow(h_Tracking_SysRel->GetBinContent(i+1),2) + pow(h_DsBR_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_BtoDs_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_CutScanAll_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_MCShape_SysRel->GetBinContent(i+1),2) + pow(h_RawRooFitYield_pdfVar_RelErr->GetBinContent(i+1),2) + pow(h_DdlsScale_SysRel->GetBinContent(i+1),2));
		h_Total_SysRel->SetBinContent(i+1, Total_Sys);

		cout<<"pt bin : "<<i<<endl;
		cout<<"Total_SysRel = "<<Total_Sys<<endl;

		// DsOverD0_Sys=sqrt( pow(Tracking_one,2) + pow(h_DsBR_SysRel->GetBinContent(i+1),2) + pow(h_f0_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_BtoDs_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_CutScanAll_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_MCShape_SysRel->GetBinContent(i+1),2) + pow(h_RawRooFitYield_pdfVar_RelErr->GetBinContent(i+1),2) +pow(h_D0_SysRel->GetBinContent(i+1),2) );
		// DsOverD0_Sys=sqrt( pow(Tracking_one,2) + pow(h_DsBR_SysRel->GetBinContent(i+1),2) + pow(h_PhiRatio_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_BtoDs_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_CutScanAll_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_MCShape_SysRel->GetBinContent(i+1),2) + pow(h_RawRooFitYield_pdfVar_RelErr->GetBinContent(i+1),2) +pow(h_D0_SysRel->GetBinContent(i+1),2) );
		// DsOverD0_Sys=sqrt( pow(Tracking_one,2) + pow(h_DsBR_SysRel->GetBinContent(i+1),2) + pow(h_PhiRatio_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_BtoDs_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_CutScanAll_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_MCShape_SysRel->GetBinContent(i+1),2) + pow(h_RawRooFitYield_pdfVar_RelErr->GetBinContent(i+1),2) +pow(h_D0_SysRel->GetBinContent(i+1),2) + pow(h_DdlsScale_SysRel->GetBinContent(i+1),2));
	
		// remove PhiRatioErr from sys
		DsOverD0_Sys=sqrt( pow(Tracking_one,2) + pow(h_DsBR_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_BtoDs_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_CutScanAll_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_MCShape_SysRel->GetBinContent(i+1),2) + pow(h_RawRooFitYield_pdfVar_RelErr->GetBinContent(i+1),2) +pow(h_D0_SysRel->GetBinContent(i+1),2) + pow(h_DdlsScale_SysRel->GetBinContent(i+1),2));

		h_DsOverD0_SysRel->SetBinContent(i+1,DsOverD0_Sys);
		h_TrkOne_SysRel->SetBinContent(i+1,Tracking_one);

		cout<<"DsOverD0_Sys = "<<DsOverD0_Sys<<endl;

		// partialRaa_Sys=sqrt( pow(h_PromptDs_CutScanAll_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_MCShape_SysRel->GetBinContent(i+1),2) + pow(h_RawRooFitYield_pdfVar_RelErr->GetBinContent(i+1),2));
		partialRaa_Sys=sqrt( pow(h_PromptDs_CutScanAll_SysRel->GetBinContent(i+1),2) + pow(h_PromptDs_MCShape_SysRel->GetBinContent(i+1),2) + pow(h_RawRooFitYield_pdfVar_RelErr->GetBinContent(i+1),2) + pow(h_DdlsScale_SysRel->GetBinContent(i+1),2) );

		h_partial_forRaa_SysRel->SetBinContent(i+1,partialRaa_Sys);
		cout<<"partialRaa_Sys = "<<partialRaa_Sys<<endl;

	}



	double rangehi=40;
	double rangelo=6;

	if(isPbPb==3){
		h_Total_SysRel->GetXaxis()->SetRangeUser(rangelo,rangehi);
		h_DsOverD0_SysRel->GetXaxis()->SetRangeUser(rangelo,rangehi);
//		h_PromptDs_BtoDs_SysRel->GetXaxis()->SetRangeUser(rangelo,rangehi);
//		h_PromptDs_CutScanAll_SysRel->GetXaxis()->SetRangeUser(rangelo,rangehi);
//		h_PromptDs_MCShape_SysRel->GetXaxis()->SetRangeUser(rangelo,rangehi);
//		h_RawRooFitYield_pdfVar_RelErr->GetXaxis()->SetRangeUser(rangelo,rangehi);
//		h_Tracking_SysRel->GetXaxis()->SetRangeUser(rangelo,rangehi);
//		h_DsBR_SysRel->GetXaxis()->SetRangeUser(rangelo,rangehi);
	}

	


	TCanvas *c_Sys=new TCanvas("c_Sys","c_Sys",c_wtopx,c_wtopy,c_W,c_H);
	c_Sys->cd();
  c_Sys->SetFillColor(0);
  c_Sys->SetBorderMode(0);
  c_Sys->SetFrameFillStyle(0);
  c_Sys->SetFrameBorderMode(0);
  c_Sys->SetLeftMargin( c_Lmg );
  c_Sys->SetRightMargin( c_Rmg );
  c_Sys->SetTopMargin( c_Tmg );
  c_Sys->SetBottomMargin( c_Bmg );
  c_Sys->SetTickx(0);
  c_Sys->SetTicky(0);

	gPad->SetLogx();

	h_Total_SysRel->SetMaximum(0.7);
	h_Total_SysRel->SetMinimum(0);
	h_Total_SysRel->SetFillColor(kGray);
	h_Total_SysRel->SetFillStyle(1001);
	h_Total_SysRel->SetTitle("");
	h_Total_SysRel->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h_Total_SysRel->GetXaxis()->CenterTitle();
	h_Total_SysRel->GetYaxis()->SetTitle("uncertainty");
	h_Total_SysRel->GetYaxis()->CenterTitle();

	h_Total_SysRel->Draw("SAME E2");
	
	h_PromptDs_BtoDs_SysRel->SetLineColor(kRed+2);
	h_PromptDs_BtoDs_SysRel->Draw("SAME");
	h_PromptDs_CutScanAll_SysRel->SetLineColor(kOrange+2);
	h_PromptDs_CutScanAll_SysRel->Draw("SAME");
  cout<<"check 7"<<endl;

	h_DdlsScale_SysRel->SetLineColor(kOrange-6);
	h_DdlsScale_SysRel->Draw("same");

	h_PromptDs_MCShape_SysRel->SetLineColor(kGreen+4);
	h_PromptDs_MCShape_SysRel->Draw("SAME");
	h_RawRooFitYield_pdfVar_RelErr->SetLineColor(kBlue+1);
	h_RawRooFitYield_pdfVar_RelErr->Draw("SAME");

  cout<<"check 8"<<endl;

	h_Tracking_SysRel->SetLineColor(kViolet+2);
	h_Tracking_SysRel->Draw("SAME");
	h_DsBR_SysRel->SetLineColor(kCyan+2);
	h_DsBR_SysRel->Draw("SAME");

	h_f0_SysRel->SetLineColor(kYellow+3);
	// h_f0_SysRel->Draw("SAME");

	// h_PhiRatio_SysRel->SetLineColor(kYellow+3);
	// h_PhiRatio_SysRel->Draw("SAME");

  cout<<"check 9"<<endl;

	TLegend *le=new TLegend(0.4,0.56,0.75,0.88);
	le->SetBorderSize(0);
	le->AddEntry(h_Tracking_SysRel,"Tracking Efficiency","l");
	le->AddEntry(h_PromptDs_CutScanAll_SysRel,"Selection efficiency","l");
	le->AddEntry(h_DdlsScale_SysRel,"Decay Length Error Scale","l");
	le->AddEntry(h_RawRooFitYield_pdfVar_RelErr,"Signal extraction","l");
	le->AddEntry(h_PromptDs_MCShape_SysRel,"MC p_{T} Shape","l");
	// le->AddEntry(h_f0_SysRel,"MC f0 Mass shape","l");
	// le->AddEntry(h_PhiRatio_SysRel,"#phi#pi Channel Ratio","l");
	le->AddEntry(h_PromptDs_BtoDs_SysRel,"Non-Prompt D_{S}","l");
	le->AddEntry(h_DsBR_SysRel,"D_{S} Branching ratio","l");
	le->AddEntry(h_Total_SysRel,"Total cross section","lf");
	le->Draw("same");

//	TLatex tl=new TLatex();
//	tl->Draw()
	texCmsPre->Draw("SAME");
	if(str_PbPb=="pp"){
	texColpp->Draw("SAME");
	}else{
	texColPbPb->Draw("SAME");
	}

	SavePlotDirs(c_Sys,Form("SystSum_%s",str_PbPb.Data()),{"Systematics"});

	// h_Lumi_SysRel->Draw("SAME");
	// h_f0_SysRel->Draw("SAME");

	TCanvas *c_SysD0Ds=new TCanvas("c_SysD0Ds","c_SysD0Ds",c_wtopx,c_wtopy,c_W,c_H);
	c_SysD0Ds->cd();
	SetCanvas(c_SysD0Ds);
	gPad->SetLogx();

	h_DsOverD0_SysRel->SetMaximum(0.7);
	h_DsOverD0_SysRel->SetMinimum(0);
	h_DsOverD0_SysRel->SetFillColor(kGray);
	h_DsOverD0_SysRel->SetFillStyle(1001);
	h_DsOverD0_SysRel->SetTitle("");
	h_DsOverD0_SysRel->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h_DsOverD0_SysRel->GetXaxis()->CenterTitle();
	h_DsOverD0_SysRel->GetYaxis()->SetTitle("uncertainty");
	h_DsOverD0_SysRel->GetYaxis()->CenterTitle();
	
	h_DsOverD0_SysRel->Draw("SAME E2");

	h_PromptDs_BtoDs_SysRel->Draw("SAME");
	h_PromptDs_CutScanAll_SysRel->Draw("SAME");
	h_DdlsScale_SysRel->Draw("SAME");
	h_PromptDs_MCShape_SysRel->Draw("SAME");
	h_RawRooFitYield_pdfVar_RelErr->Draw("SAME");
	h_DsBR_SysRel->Draw("SAME");
	// h_f0_SysRel->Draw("SAME");
	// h_PhiRatio_SysRel->Draw("SAME");

	h_TrkOne_SysRel->SetLineColor(kViolet+2);
	h_TrkOne_SysRel->Draw("SAME");
	h_D0_SysRel->SetLineColor(kMagenta+1);
	h_D0_SysRel->Draw("SAME");

	texCmsPre->Draw("SAME");
	if(str_PbPb=="pp"){
	texColpp->Draw("SAME");
	}else{
	texColPbPb->Draw("SAME");
	}

	TLegend *le_DsD0=new TLegend(0.4,0.56,0.75,0.88);
	le_DsD0->SetBorderSize(0);
	le_DsD0->AddEntry(h_TrkOne_SysRel,"Tracking Efficiency","l");
	le_DsD0->AddEntry(h_PromptDs_CutScanAll_SysRel,"Selection efficiency","l");
	le_DsD0->AddEntry(h_DdlsScale_SysRel,"Decay Length Error Scale","l");
	le_DsD0->AddEntry(h_RawRooFitYield_pdfVar_RelErr,"Signal extraction","l");
	le_DsD0->AddEntry(h_PromptDs_MCShape_SysRel,"MC p_{T} Shape","l");
	// le_DsD0->AddEntry(h_f0_SysRel,"MC f0 Mass shape","l");
	// le_DsD0->AddEntry(h_PhiRatio_SysRel,"#phi#pi Channel Ratio","l");
	le_DsD0->AddEntry(h_PromptDs_BtoDs_SysRel,"Non-Prompt D_{S}","l");
	le_DsD0->AddEntry(h_DsBR_SysRel,"D_{S} Branching ratio","l");
	le_DsD0->AddEntry(h_D0_SysRel,"D^{0}","l");
	le_DsD0->AddEntry(h_DsOverD0_SysRel,"Total D_{S}/D^{0}","lf");
	le_DsD0->Draw("same");

	SavePlotDirs(c_SysD0Ds,Form("SystSum_D0Ds_%s",str_PbPb.Data()),{"Systematics"});


	for(int i=startbin; i<nbin_pt; i++){
		cout<<"\n Cross section "<<str_PbPb<<" , bin : "<<bins_pt[i]<<" - "<<bins_pt[i+1]<<endl;
		cout<<setw(18)<<"Selection eff : "<<setw(15)<<setprecision(1)<<std::fixed<<h_PromptDs_CutScanAll_SysRel->GetBinContent(i+1)*100<<endl;
		cout<<setw(18)<<"Ddl Err Scale : "<<setw(15)<<setprecision(1)<<std::fixed<<h_DdlsScale_SysRel->GetBinContent(i+1)*100<<endl;
		cout<<setw(18)<<"SignalExt eff : "<<setw(15)<<setprecision(1)<<std::fixed<<h_RawRooFitYield_pdfVar_RelErr->GetBinContent(i+1)*100<<endl;
		cout<<setw(18)<<"MC pt : "<<setw(15)<<setprecision(1)<<std::fixed<<h_PromptDs_MCShape_SysRel->GetBinContent(i+1)*100<<endl;
		// cout<<setw(18)<<"MC f0 : "<<setw(15)<<h_PromptDs_f0Eff_SysRel->GetBinContent(i+1)*100<<endl;
		cout<<setw(18)<<"MC nonprompt : "<<setw(15)<<setprecision(1)<<std::fixed<<h_PromptDs_BtoDs_SysRel->GetBinContent(i+1)*100<<endl;
		cout<<setw(18)<<"Total : "<<setw(15)<<setprecision(1)<<std::fixed<<h_Total_SysRel->GetBinContent(i+1)*100<<endl;

	}


	for(int i=startbin; i<nbin_pt; i++){
		cout<<"\n DsD0 "<<str_PbPb<<" , bin : "<<bins_pt[i]<<" - "<<bins_pt[i+1]<<endl;

		cout<<setw(18)<<"Selection eff : "<<setw(15)<<setprecision(1)<<h_PromptDs_CutScanAll_SysRel->GetBinContent(i+1)*100<<endl;
		cout<<setw(18)<<"Ddl Err Scale  : "<<setw(15)<<setprecision(1)<<h_DdlsScale_SysRel->GetBinContent(i+1)*100<<endl;
		cout<<setw(18)<<"SignalExt eff : "<<setw(15)<<setprecision(1)<<h_RawRooFitYield_pdfVar_RelErr->GetBinContent(i+1)*100<<endl;
		cout<<setw(18)<<"MC pt : "<<setw(15)<<setprecision(1)<<h_PromptDs_MCShape_SysRel->GetBinContent(i+1)*100<<endl;
		// cout<<setw(18)<<"MC f0 : "<<setw(15)<<setprecision(1)<<h_PromptDs_f0Eff_SysRel->GetBinContent(i+1)*100<<endl;
		cout<<setw(18)<<"MC nonprompt : "<<setw(15)<<setprecision(1)<<h_PromptDs_BtoDs_SysRel->GetBinContent(i+1)*100<<endl;
		cout<<setw(18)<<"Total : "<<setw(15)<<setprecision(1)<<h_DsOverD0_SysRel->GetBinContent(i+1)*100<<endl;

	}

		// cout<<"\n DsD0 "<<str_PbPb<<" , bin : "<<bins_pt[i]<<" - "<<bins_pt[i+1]<<endl;

		cout<<"\n\n --------------------\nCross section "<<str_PbPb<<endl;
		cout<<"Selection efficiency ";
	for(int i=startbin; i<nbin_pt; i++){
		cout<<" & "<<setprecision(1)<<setw(4)<<h_PromptDs_CutScanAll_SysRel->GetBinContent(i+1)*100;

	}
		cout<<" \\\\ \\hline"<<endl;

		cout<<"Signal extraction ";
	for(int i=startbin; i<nbin_pt; i++){
		cout<<" & "<<setprecision(1)<<setw(4)<<h_RawRooFitYield_pdfVar_RelErr->GetBinContent(i+1)*100;

	}
		cout<<" \\\\ \\hline"<<endl;

		cout<<"MC $p_{T}$ shape ";
	for(int i=startbin; i<nbin_pt; i++){
		cout<<" & "<<setprecision(1)<<setw(4)<<h_PromptDs_MCShape_SysRel->GetBinContent(i+1)*100;

	}
		cout<<" \\\\ \\hline"<<endl;

		cout<<"MC Decay Length Tune ";
	for(int i=startbin; i<nbin_pt; i++){
		cout<<" & "<<setprecision(1)<<setw(4)<<h_DdlsScale_SysRel->GetBinContent(i+1)*100;

	}
		cout<<" \\\\ \\hline"<<endl;

	int ncol=8;
	if(isPbPb){ncol=4;}

		cout<<"$\\Ds$ $\\to$ $\\phi$ $\\pi^{\\pm}$  ratio ";
		cout<<" & \\multicolumn{"<<ncol<<"}{c|}{"<<setprecision(1)<<PhiRatioErr*100<<"}";
		cout<<" \\\\ \\hline"<<endl;

		cout<<"Non-prompt \\Ds";
	for(int i=startbin; i<nbin_pt; i++){
		cout<<" & "<<setprecision(1)<<setw(4)<<h_PromptDs_BtoDs_SysRel->GetBinContent(i+1)*100;

	}
		cout<<" \\\\ \\hline"<<endl;

		cout<<"Branching ratio";
		cout<<" & \\multicolumn{"<<ncol<<"}{c|}{"<<setprecision(1)<<3.5<<"}";
		cout<<" \\\\ \\hline"<<endl;

		cout<<"Total bin by bin ";
	for(int i=startbin; i<nbin_pt; i++){
		cout<<" & "<<setprecision(1)<<setw(4)<<h_Total_SysRel->GetBinContent(i+1)*100;

	}
		cout<<" \\\\ \\hline"<<endl;

		cout<<"\n\n------------\n  DsOverD0 total"<<endl;
		cout<<"Total bin by bin ";
	for(int i=startbin; i<nbin_pt; i++){
		cout<<" & "<<setprecision(1)<<setw(4)<<h_DsOverD0_SysRel->GetBinContent(i+1)*100;

	}
		cout<<" \\\\ \\hline"<<endl;









	fout->cd();
	h_Total_SysRel->Write();
	h_partial_forRaa_SysRel->Write();
	h_DsOverD0_SysRel->Write();
	h_PromptDs_BtoDs_SysRel->Write();
	h_PromptDs_BtoDs_SysRel_ForRaa->Write();
	h_PromptDs_BtoDs_SysRel_ForRaaCancel->Write();
	h_PromptDs_CutScanAll_SysRel->Write();
	h_DdlsScale_SysRel->Write();
	h_RawRooFitYield_pdfVar_RelErr->Write();
	h_PromptDs_MCShape_SysRel->Write();
	// fout->Close();


	return 0;

}
