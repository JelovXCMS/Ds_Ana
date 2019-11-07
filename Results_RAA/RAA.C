
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting_simple.h"
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
#include <TGraph.h>
#include <TGraphAsymmErrors.h>


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
#include <RooFitResult.h>

#include "TFitter.h"
#include "TFitResult.h"

using namespace RooFit;
using namespace std;

// Bool_t verbose=true;
Int_t StartBin=0;



void RAA(int startbin=1){

	// gStyle->SetOptStat(0);
	InitStyle();
	// c_Tmg=0.15;
	// c_Bmg=0.15;
	setTDRStyle();
  double tex_upperY=0.95;

  // TLatex* texCmsPre = new TLatex(0.15,tex_upperY, "#scale[1.25]{CMS} Preliminary");
  TLatex* texCmsPre = new TLatex(0.15,tex_upperY, "#scale[1.25]{#bf{CMS}}");
  texCmsPre->SetNDC();
  texCmsPre->SetTextAlign(12);
  texCmsPre->SetTextSize(0.045);
  texCmsPre->SetTextFont(42);

  TLatex* texCms2 = new TLatex(0.20,tex_upperY-0.08, "#scale[1.25]{#bf{CMS}}");
  texCms2->SetNDC();
  texCms2->SetTextAlign(12);
  texCms2->SetTextSize(0.045);
  texCms2->SetTextFont(42);

  TLatex* texPre2 = new TLatex(0.20,tex_upperY-0.16, "Preliminary");
  texPre2->SetNDC();
  texPre2->SetTextAlign(12);
  texPre2->SetTextSize(0.045);
  texPre2->SetTextFont(42);



  TLatex* texCmsSim = new TLatex(0.15,tex_upperY, "#scale[1.25]{#bf{CMS}} Simulations");
  texCmsSim->SetNDC();
  texCmsSim->SetTextAlign(12);
  texCmsSim->SetTextSize(0.045);
  texCmsSim->SetTextFont(42);


  TLatex* texColPbPb = new TLatex(0.95,tex_upperY, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  texColPbPb->SetNDC();
  texColPbPb->SetTextAlign(32);
  texColPbPb->SetTextSize(0.045);
  texColPbPb->SetTextFont(42);

  TLatex* texColpp = new TLatex(0.95,tex_upperY, "pp #sqrt{s} = 5.02 TeV");
  texColpp->SetNDC();
  texColpp->SetTextAlign(32);
  texColpp->SetTextSize(0.045);
  texColpp->SetTextFont(42);

	// TLatex* texppPbPb = new TLatex(0.95, tex_upperY, "27.4 pb^{-1} (5.02 TeV pp) + 530 #mub^{-1} (5.02 TeV PbPb)");
	TLatex* texppPbPb = new TLatex(0.95, tex_upperY, "pp 38 nb^{-1}, PbPb 44 #mub^{-1} (5.02 TeV)");
  texppPbPb->SetNDC();
  texppPbPb->SetTextAlign(32);
  texppPbPb->SetTextSize(0.045);
  texppPbPb->SetTextFont(42);

	TLatex *textemp=new TLatex();

	// get Theory prediction

	TFile *f_theo=TFile::Open("./output/TheoryPredict.root");
	TGraph* gr_TAMU_DsD0_pp=(TGraph*)f_theo->Get("gr_TAMU_DsD0_pp");
	TGraph* gr_TAMU_DsD0_PbPb020_noHadDif=(TGraph*)f_theo->Get("gr_TAMU_DsD0_PbPb020_noHadDif");
	TGraph* gr_TAMU_DsD0_PbPb020_withHadDif=(TGraph*)f_theo->Get("gr_TAMU_DsD0_PbPb020_withHadDif");
	TGraph* gr_CuJet_Raa_alpha1cm0p28=(TGraph*)f_theo->Get("gr_CuJet_Raa_alpha1cm0p28");
	TGraph* gr_CuJet_Raa_alpha0p8cm0p22=(TGraph*)f_theo->Get("gr_CuJet_Raa_alpha0p8cm0p22");
	TGraph* gr_PHSD_Raa=(TGraph*)f_theo->Get("gr_PHSD_Raa");
	TGraph* gr_PHSD_DsD0_pp=(TGraph*)f_theo->Get("gr_PHSD_DsD0_pp");
	TGraph* gr_PHSD_DsD0_PbPb=(TGraph*)f_theo->Get("gr_PHSD_DsD0_PbPb");
	TGraph* gr_PHSD_DsD0_PbPbpp=(TGraph*)f_theo->Get("gr_PHSD_DsD0_PbPbpp");
	// gr_TAMU_DsD0_pp->Draw("l");


	// StartBin=startbin;

	TFile *fout=TFile::Open("output/RAA.root","RECREATE");

	fout->cd();

	TFile *fin_pp=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/CrossSection_dNdpt/output/PromptDsCrossSectionPP_FixShape.root","OPEN");

	TFile *fin_PbPb3=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/CrossSection_dNdpt/output/PromptDsdNdptPbPb3_FixShape.root","OPEN");

	TH1D *h_PromptDs_CS_pp=(TH1D*)fin_pp->Get("h_PromptDs_pp");
	TH1D *h_PromptDs_dNdpt_PbPb3=(TH1D*)fin_PbPb3->Get("h_PromptDs_PbPb3");

 // adding phi channel statistic error to CS
	double phiRatio_sta_pp=0.002/0.96;
	double phiRatio_sta_PbPb=0.007/0.901;

	for(int i=1; i<=h_PromptDs_CS_pp->GetNbinsX(); i++){
		double ori_sta_rel=h_PromptDs_CS_pp->GetBinError(i)/h_PromptDs_CS_pp->GetBinContent(i);
		double new_sta_rel=sqrt(ori_sta_rel*ori_sta_rel+phiRatio_sta_pp*phiRatio_sta_pp);
		double new_sta_abs=new_sta_rel*h_PromptDs_CS_pp->GetBinContent(i);
		h_PromptDs_CS_pp->SetBinError(i,new_sta_abs);
		cout<<"bin : "<<i<<" , ori_stat = "<<ori_sta_rel<<" , new stat = "<<new_sta_rel<<endl;
	}

	for(int i=3; i<=h_PromptDs_dNdpt_PbPb3->GetNbinsX(); i++){
		double ori_sta_rel=h_PromptDs_dNdpt_PbPb3->GetBinError(i)/h_PromptDs_dNdpt_PbPb3->GetBinContent(i);
		double new_sta_rel=sqrt(ori_sta_rel*ori_sta_rel+phiRatio_sta_pp*phiRatio_sta_pp);
		double new_sta_abs=new_sta_rel*h_PromptDs_dNdpt_PbPb3->GetBinContent(i);
		h_PromptDs_dNdpt_PbPb3->SetBinError(i,new_sta_abs);
		cout<<"bin : "<<i<<" , ori_stat = "<<ori_sta_rel<<" , new stat = "<<new_sta_rel<<endl;
	}

		// h_PromptDs_dNdpt_PbPb3->Draw();


	// return ;


	// TH1D *h_PromptDs_CS_pp_Sys=(TH1D*)h_PromptDs_CS_pp->Clone("h_PromptDs_CS_pp_Sys");
	TH1D *h_PromptDs_CS_pp_Sys=new TH1D("h_PromptDs_CS_pp_Sys","h_PromptDs_CS_pp_Sys",nbin_pt_pp,bins_pt_pp); h_PromptDs_CS_pp_Sys->Sumw2();
	TH1D *h_PromptDs_CS_pp_Sta=new TH1D("h_PromptDs_CS_pp_Sta","h_PromptDs_CS_pp_Sta",nbin_pt_pp,bins_pt_pp); h_PromptDs_CS_pp_Sta->Sumw2(); // what is this ??
//	TH1D *h_PromptDs_dNdpt_PbPb3_Sys=(TH1D*)h_PromptDs_dNdpt_PbPb3->Clone("h_PromptDs_dNdpt_PbPb3_Sys");
	TH1D *h_PromptDs_dNdpt_PbPb3_Sys=new TH1D("h_PromptDs_dNdpt_PbPb3_Sys","h_PromptDs_dNdpt_PbPb3_Sys",nbin_pt_PbPb3,bins_pt_PbPb3); h_PromptDs_dNdpt_PbPb3_Sys->Sumw2();
	TH1D *h_PromptDs_dNdpt_PbPb3_Sta=new TH1D("h_PromptDs_dNdpt_PbPb3_Sta","h_PromptDs_dNdpt_PbPb3_Sta",nbin_pt_PbPb3,bins_pt_PbPb3); h_PromptDs_dNdpt_PbPb3_Sta->Sumw2();

	TH1D *h_RAA=new TH1D("h_RAA","h_RAA",nbin_pt_PbPb3,bins_pt_PbPb3); h_RAA->Sumw2();
	TH1D *h_RAA_Sta=new TH1D("h_RAA_Sta","h_RAA_Sta",nbin_pt_PbPb3,bins_pt_PbPb3); h_RAA_Sta->Sumw2();

	// for Ds/D0 ratio

  TFile *f_NonPromptDs=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BtoDs_Results/output/BtoDs.root");
  TH1D *hPromptDCrossSectionPP_AnaBin=(TH1D*)f_NonPromptDs->Get("hPromptDCrossSectionPP_AnaBin");
  TH1D *hPromptDdNdPtPbPb_AnaBin=(TH1D*)f_NonPromptDs->Get("hPromptDdNdPtPbPb_AnaBin");

	TH1D *h_DsoverD0_pp_Sta=new TH1D ("h_DsoverD0_pp_Sta","h_DsoverD0_pp_Sta",nbin_pt_pp,bins_pt_pp);
	TH1D *h_DsoverD0_pp_Sys=new TH1D ("h_DsoverD0_pp_Sys","h_DsoverD0_pp_Sys",nbin_pt_pp,bins_pt_pp);

	TH1D *h_DsoverD0_PbPb_Sta=new TH1D ("h_DsoverD0_PbPb_Sta","h_DsoverD0_PbPb_Sta",nbin_pt_PbPb3,bins_pt_PbPb3);
	TH1D *h_DsoverD0_PbPb_Sys=new TH1D ("h_DsoverD0_PbPb_Sys","h_DsoverD0_PbPb_Sys",nbin_pt_PbPb3,bins_pt_PbPb3);

	// TH1D *h_DsoverD0_DRatio_Sta=new TH1D ("h_DsoverD0_DRatio_Sta","h_DsoverD0_DRatio_Sta",nbin_pt_PbPb3,bins_pt_PbPb3);
	// TH1D *h_DsoverD0_DRatio_Sys=new TH1D ("h_DsoverD0_DRatio_Sys","h_DsoverD0_DRatio_Sys",nbin_pt_PbPb3,bins_pt_PbPb3);
	TH1D *h_DsoverD0_DRatio_Sta=new TH1D ("h_DsoverD0_DRatio_Sta","h_DsoverD0_DRatio_Sta",nbin_pt_pp,bins_pt_pp);
	TH1D *h_DsoverD0_DRatio_Sys=new TH1D ("h_DsoverD0_DRatio_Sys","h_DsoverD0_DRatio_Sys",nbin_pt_pp,bins_pt_pp);

	TFile *f_Sys=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/Systematics/output/SysSumRaa.root");
	TH1D *h_Total_SysRel_pp=(TH1D*)f_Sys->Get("h_Total_SysRel_pp");
	TH1D *h_DsOverD0_SysRel_pp=(TH1D*)f_Sys->Get("h_DsOverD0_SysRel_pp");
	TH1D *h_Total_SysRel_PbPb3=(TH1D*)f_Sys->Get("h_Total_SysRel_PbPb3");
	TH1D *h_DsOverD0_SysRel_PbPb3=(TH1D*)f_Sys->Get("h_DsOverD0_SysRel_PbPb3");
	TH1D *h_Raa_SysRel=(TH1D*)f_Sys->Get("h_Raa_SysRel");

	Double_t SysRel_pp[nbin_pt_pp]={0.368,0.2701,0.2089,0.1885,0.1726,0.1866,0.2053,0.3813};
	Double_t SysRel_PbPb[nbin_pt_PbPb3]={0,0,0.386,0.398,0.4195,0.2505};
	Double_t SysRel_BR=0.035;

	Double_t SysRel_Raa[nbin_pt_PbPb3]={0,0,0.3875,0.4254,0.44746,0.44131};

	// Double_t DsD0SysRel_pp[nbin_pt_pp]={0.411,0.347,0.293,0.279,0.268,0.276,0.291,0.432};
	// Double_t DsD0SysRel_PbPb[nbin_pt_PbPb3]={0,0,0.423,0.451,0.474,0.340};

	Double_t DsD0SysRel_pp[nbin_pt_pp]={0.381,0.311,0.249,0.232,0.218,0.228,0.246,0.403};
	Double_t DsD0SysRel_PbPb[nbin_pt_PbPb3]={0,0,0.380,0.411,0.436,0.285};

	Double_t DsD0SysRel_Temp=0;
	Double_t DsD0Ratio_Temp=0;

	Double_t SysRel_Temp=0;

	// Double_t D0Sys_ForDsD0DRatio_pp[nbin_pt_PbPb3]={};
	Double_t DsD0_DoubleRatio_Temp=0;
	Double_t D0Sys_ForDsD0DRatio[nbin_pt_PbPb3]={0,0,0.1277106,0.1182328,0.1361837,0.1479392};



	for(int i=0; i<nbin_pt_pp; i++){

		SysRel_pp[i]=h_Total_SysRel_pp->GetBinContent(i+1);
		DsD0SysRel_pp[i]=h_DsOverD0_SysRel_pp->GetBinContent(i+1);
		// SysRel_Temp=sqrt(SysRel_pp[i]*SysRel_pp[i]+SysRel_BR*SysRel_BR);
		SysRel_Temp=sqrt(SysRel_pp[i]*SysRel_pp[i]);
		h_PromptDs_CS_pp_Sys->SetBinError(i+1, h_PromptDs_CS_pp->GetBinContent(i+1)*SysRel_Temp );
		h_PromptDs_CS_pp_Sys->SetBinContent(i+1, h_PromptDs_CS_pp->GetBinContent(i+1));

		h_PromptDs_CS_pp_Sta->SetBinError(i+1, h_PromptDs_CS_pp->GetBinError(i+1) );
		h_PromptDs_CS_pp_Sta->SetBinContent(i+1, h_PromptDs_CS_pp->GetBinContent(i+1));

		cout<<"bin i = "<<i<<" , SysRel_Temp = "<<SysRel_Temp<<endl;

		DsD0Ratio_Temp= h_PromptDs_CS_pp->GetBinContent(i+1)/hPromptDCrossSectionPP_AnaBin->GetBinContent(i+1);
		h_DsoverD0_pp_Sta->SetBinContent(i+1,DsD0Ratio_Temp);
		h_DsoverD0_pp_Sta->SetBinError(i+1, DsD0Ratio_Temp * sqrt(   hPromptDCrossSectionPP_AnaBin->GetBinError(i+1)/hPromptDCrossSectionPP_AnaBin->GetBinContent(i+1)*hPromptDCrossSectionPP_AnaBin->GetBinError(i+1)/hPromptDCrossSectionPP_AnaBin->GetBinContent(i+1) + h_PromptDs_CS_pp->GetBinError(i+1)/h_PromptDs_CS_pp->GetBinContent(i+1)*h_PromptDs_CS_pp->GetBinError(i+1)/h_PromptDs_CS_pp->GetBinContent(i+1)  ));
		h_DsoverD0_pp_Sys->SetBinContent(i+1,DsD0Ratio_Temp);
		h_DsoverD0_pp_Sys->SetBinError(i+1,DsD0Ratio_Temp*DsD0SysRel_pp[i]);

	}

	for(int i=0; i<nbin_pt_PbPb3; i++){
		SysRel_PbPb[i]=h_Total_SysRel_PbPb3->GetBinContent(i+1);
		DsD0SysRel_PbPb[i]=h_DsOverD0_SysRel_PbPb3->GetBinContent(i+1);
		// SysRel_Temp=sqrt(SysRel_PbPb[i]*SysRel_PbPb[i]+SysRel_BR*SysRel_BR);
		SysRel_Temp=sqrt(SysRel_PbPb[i]*SysRel_PbPb[i]);
		SysRel_Raa[i]=h_Raa_SysRel->GetBinContent(i+1);
		if(i>=2){
		h_PromptDs_dNdpt_PbPb3_Sys->SetBinContent(i+1, h_PromptDs_dNdpt_PbPb3->GetBinContent(i+1)/TAA0to100);
		h_PromptDs_dNdpt_PbPb3_Sys->SetBinError(i+1, h_PromptDs_dNdpt_PbPb3->GetBinContent(i+1)*SysRel_Temp/TAA0to100);

		h_PromptDs_dNdpt_PbPb3_Sta->SetBinContent(i+1, h_PromptDs_dNdpt_PbPb3->GetBinContent(i+1)/TAA0to100);
		h_PromptDs_dNdpt_PbPb3_Sta->SetBinError(i+1, h_PromptDs_dNdpt_PbPb3->GetBinError(i+1)/TAA0to100);

		h_RAA->SetBinContent(i+1, h_PromptDs_dNdpt_PbPb3_Sys->GetBinContent(i+1)/ h_PromptDs_CS_pp_Sys->GetBinContent(i+3));
		h_RAA->SetBinError(i+1, h_RAA->GetBinContent(i+1)*SysRel_Raa[i]);	

		h_RAA_Sta->SetBinContent(i+1, h_PromptDs_dNdpt_PbPb3_Sys->GetBinContent(i+1)/ h_PromptDs_CS_pp_Sys->GetBinContent(i+3));
		h_RAA_Sta->SetBinError(i+1, h_PromptDs_dNdpt_PbPb3_Sys->GetBinContent(i+1)/ h_PromptDs_CS_pp_Sys->GetBinContent(i+3) * sqrt( h_PromptDs_CS_pp_Sta->GetBinError(i+3)/h_PromptDs_CS_pp_Sta->GetBinContent(i+3)*h_PromptDs_CS_pp_Sta->GetBinError(i+3)/h_PromptDs_CS_pp_Sta->GetBinContent(i+3) + h_PromptDs_dNdpt_PbPb3_Sta->GetBinError(i+1)/h_PromptDs_dNdpt_PbPb3_Sta->GetBinContent(i+1)* h_PromptDs_dNdpt_PbPb3_Sta->GetBinError(i+1)/h_PromptDs_dNdpt_PbPb3_Sta->GetBinContent(i+1)   )    );	

		// for DsoverD0		
		DsD0Ratio_Temp= h_PromptDs_dNdpt_PbPb3->GetBinContent(i+1)/hPromptDdNdPtPbPb_AnaBin->GetBinContent(i+1);
		h_DsoverD0_PbPb_Sta->SetBinContent(i+1,DsD0Ratio_Temp);
		h_DsoverD0_PbPb_Sta->SetBinError(i+1, DsD0Ratio_Temp * sqrt(   hPromptDdNdPtPbPb_AnaBin->GetBinError(i+1)/hPromptDdNdPtPbPb_AnaBin->GetBinContent(i+1)*hPromptDdNdPtPbPb_AnaBin->GetBinError(i+1)/hPromptDdNdPtPbPb_AnaBin->GetBinContent(i+1) + h_PromptDs_dNdpt_PbPb3->GetBinError(i+1)/h_PromptDs_dNdpt_PbPb3->GetBinContent(i+1)*h_PromptDs_dNdpt_PbPb3->GetBinError(i+1)/h_PromptDs_dNdpt_PbPb3->GetBinContent(i+1)  ));
		h_DsoverD0_PbPb_Sys->SetBinContent(i+1,DsD0Ratio_Temp);
		h_DsoverD0_PbPb_Sys->SetBinError(i+1,DsD0Ratio_Temp*DsD0SysRel_PbPb[i]);

		// DsoverD0 pp/PbPb Double Ratio
		DsD0_DoubleRatio_Temp=h_DsoverD0_PbPb_Sys->GetBinContent(i+1)/h_DsoverD0_pp_Sys->GetBinContent(i+1+2);	
		h_DsoverD0_DRatio_Sta->SetBinContent(i+1+2,DsD0_DoubleRatio_Temp);
		h_DsoverD0_DRatio_Sta->SetBinError(i+1+2,DsD0_DoubleRatio_Temp*sqrt( pow(h_DsoverD0_PbPb_Sta->GetBinError(i+1)/h_DsoverD0_PbPb_Sta->GetBinContent(i+1),2)+ pow( h_DsoverD0_pp_Sta->GetBinError(i+3)/h_DsoverD0_pp_Sta->GetBinContent(i+3),2 )   ) );

		h_DsoverD0_DRatio_Sys->SetBinContent(i+1+2,DsD0_DoubleRatio_Temp);
		h_DsoverD0_DRatio_Sys->SetBinError(i+1+2,DsD0_DoubleRatio_Temp*sqrt( pow(D0Sys_ForDsD0DRatio[i],2)+ pow(SysRel_Raa[i],2)-pow(0.18,2)+pow(0.06,2)-2*pow(0.03,2)  ));
// 0.18 for tracking , 0.06 for a single track,  0.03 for luminosity (cancel in Ds/D0)

		}
		cout<<"bin i = "<<i<<" , SysRel_Temp = "<<SysRel_Temp<<endl;
	}

	// load PythiaRef
	TFile *f_pythia=TFile::Open("output/PythiaRef.root","READ");
	TH1D *h_Genpt_D0=(TH1D*)f_pythia->Get("h_Genpt_D0");
	TH1D *h_Genpt_Ds=(TH1D*)f_pythia->Get("h_Genpt_Ds");
	TH1D *h_DsoverD0_Pythia=(TH1D*)f_pythia->Get("h_DsOverD0_Pythia");

	TGraph *gr_DsD0_Pythia=(TGraph*)f_pythia->Get("gr_DsD0");
	TGraph *gr_Genpt_Ds_Pythia=(TGraph*)f_pythia->Get("gr_Genpt_Ds");

	// load RAA models 
	
	// PHSD this is old 
	/*
	Double_t RAA_PHSD_pt[14]={0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,6,10,14,18,22,26};
	Double_t RAA_PHSD_val[14]={0.604,0.862,1.07,1.24,1.27,1.12,0.969,0.725,0.407,0.23,0.24,0.208,0.398,0.458};
	TGraph *gr_RAA_PHSD=new TGraph(14,RAA_PHSD_pt,RAA_PHSD_val);
*/


	gStyle->SetOptStat(0);
	TCanvas *c_DsOverD0=new TCanvas("c_DsOverD0","c_DsOverD0",c_wtopx,c_wtopy,c_W,c_H);
	c_DsOverD0->cd();
  c_DsOverD0->SetFillColor(0);
  c_DsOverD0->SetBorderMode(0);
  c_DsOverD0->SetFrameFillStyle(0);
  c_DsOverD0->SetFrameBorderMode(0);
  c_DsOverD0->SetLeftMargin( c_Lmg );
  c_DsOverD0->SetRightMargin( c_Rmg );
  c_DsOverD0->SetTopMargin( c_Tmg );
  c_DsOverD0->SetBottomMargin( c_Bmg );
  c_DsOverD0->SetTickx(0);
  c_DsOverD0->SetTicky(0);

// gStyle->SetOptStat(0);
	// gPad->SetLogy();
	// gStyle->SetOptStat(0);


	// h_DsoverD0_pp_Sys->SetLineColor(2);
	h_DsoverD0_pp_Sys->SetMinimum(0);
	h_DsoverD0_pp_Sys->SetMaximum(1);
	h_DsoverD0_pp_Sys->SetTitle("");
	h_DsoverD0_pp_Sys->GetYaxis()->SetTitle("D_{S}^{+}/D^{0} Ratio");
	h_DsoverD0_pp_Sys->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h_DsoverD0_pp_Sys->GetXaxis()->CenterTitle();
	h_DsoverD0_pp_Sys->GetYaxis()->CenterTitle();
	h_DsoverD0_pp_Sys->SetFillColorAlpha(kRed-4,0.4);
	// h_DsoverD0_pp_Sys->SetFillColor(kRed-4);
	h_DsoverD0_pp_Sys->SetLineColor(kRed);
	h_DsoverD0_pp_Sys->SetMarkerColor(kRed);
	// h_DsoverD0_pp_Sys->SetFillStyle(3005);
	h_DsoverD0_pp_Sys->Draw("SAME E2");
	h_DsoverD0_pp_Sta->SetLineColor(kRed);
	h_DsoverD0_pp_Sta->SetMarkerColor(kRed);
	h_DsoverD0_pp_Sta->Draw("SAME");

	h_DsoverD0_PbPb_Sys->SetLineColor(kBlue);
	h_DsoverD0_PbPb_Sys->SetMarkerColor(kBlue);
	// h_DsoverD0_PbPb_Sys->SetFillColor(kBlue-4);
	h_DsoverD0_PbPb_Sys->SetFillColorAlpha(kBlue-8,0.4);
	// h_DsoverD0_PbPb_Sys->SetFillStyle(3004);
	h_DsoverD0_PbPb_Sys->Draw("SAME E2");

	h_DsoverD0_PbPb_Sta->SetLineColor(kBlue);
	h_DsoverD0_PbPb_Sta->SetMarkerColor(kBlue);
	h_DsoverD0_PbPb_Sta->Draw("SAME");

	// h_DsoverD0_Pythia->SetLineColor(1);
	// h_DsoverD0_Pythia->SetMarkerColor(1);
	// h_DsoverD0_Pythia->Draw("SAME");

	gr_DsD0_Pythia->SetLineColor(1);
	gr_DsD0_Pythia->SetLineWidth(3);
	gr_DsD0_Pythia->Draw("l");

	gr_TAMU_DsD0_pp->SetLineColor(kGreen+2);
	gr_TAMU_DsD0_pp->SetLineWidth(3);
	gr_TAMU_DsD0_pp->Draw("l");

	gr_PHSD_DsD0_pp->SetLineColor(kViolet+2);
	gr_PHSD_DsD0_pp->SetLineWidth(3);
	gr_PHSD_DsD0_pp->Draw("l");
	gr_PHSD_DsD0_PbPb->SetLineColor(kMagenta+1);
	gr_PHSD_DsD0_PbPb->SetLineWidth(3);
	gr_PHSD_DsD0_PbPb->Draw("l");

	gr_PHSD_DsD0_PbPbpp->SetLineColor(kMagenta+1);
	gr_PHSD_DsD0_PbPbpp->SetLineWidth(3);


	TLegend *le_DsD0= new TLegend(0.62,0.58,0.85,0.85,NULL,"brNDC");
	le_DsD0->SetBorderSize(0);
	le_DsD0->SetTextSize(0.04);
	le_DsD0->AddEntry(h_DsoverD0_pp_Sys,"pp","flp");
	le_DsD0->AddEntry(h_DsoverD0_PbPb_Sys,"PbPb 0-100\%","flp");
	// le_DsD0->AddEntry(h_DsoverD0_Pythia,"PYTHIA 8","lp");
	le_DsD0->AddEntry(gr_DsD0_Pythia,"PYTHIA 8","l");
	le_DsD0->Draw("SAME");

	TLegend *le_DsD02= new TLegend(0.20,0.58,0.45,0.85,NULL,"brNDC");
	le_DsD02->SetBorderSize(0);
	le_DsD02->SetTextSize(0.04);
	le_DsD02->AddEntry(gr_TAMU_DsD0_pp,"TAMU pp","l");
	le_DsD02->AddEntry(gr_PHSD_DsD0_pp,"PHSD pp","l");
	le_DsD02->AddEntry(gr_PHSD_DsD0_PbPb,"PHSD PbPb 0-80\%","l");
	le_DsD02->Draw("SAME");



	texppPbPb->Draw("SAME");
	texCmsPre->Draw("SAME");
	
  // gStyle->SetOptStat(0);
	gPad->Modified();
	gPad->Update();


	SavePlotDirs(c_DsOverD0,"DsOverD0",{"Result"});

	// return;

	TCanvas *c_DsD0_DRatio= new TCanvas("c_DsD0_DRatio","c_DsD0_DRatio",c_wtopx,c_wtopy,c_W,c_H);
	c_DsD0_DRatio->cd();
	c_DsD0_DRatio->SetFillColor(0);
  c_DsD0_DRatio->SetBorderMode(0);
  c_DsD0_DRatio->SetFrameFillStyle(0);
  c_DsD0_DRatio->SetFrameBorderMode(0);
  c_DsD0_DRatio->SetLeftMargin( c_Lmg );
  c_DsD0_DRatio->SetRightMargin( c_Rmg );
  c_DsD0_DRatio->SetTopMargin( c_Tmg );
  c_DsD0_DRatio->SetBottomMargin( c_Bmg );
  c_DsD0_DRatio->SetTickx(0);
  c_DsD0_DRatio->SetTicky(0);
	h_DsoverD0_DRatio_Sys->SetTitle("");
	h_DsoverD0_DRatio_Sys->SetMaximum(2);
	h_DsoverD0_DRatio_Sys->SetMinimum(0.4);
	h_DsoverD0_DRatio_Sys->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h_DsoverD0_DRatio_Sys->GetYaxis()->SetTitle("#frac{D_{S}^{+}/D^{0} (PbPb)}{ D_{S}^{+}/D^{0} (pp) }");
	h_DsoverD0_DRatio_Sys->GetXaxis()->CenterTitle();
	h_DsoverD0_DRatio_Sys->GetYaxis()->CenterTitle();
	h_DsoverD0_DRatio_Sys->SetLineColor(kBlue);
	h_DsoverD0_DRatio_Sys->SetMarkerColor(kBlue);
	// h_DsoverD0_DRatio_Sys->SetFillColor(kBlue-9);
	h_DsoverD0_DRatio_Sys->SetFillColorAlpha(kBlue-8,0.4);
	// h_DsoverD0_DRatio_Sys->SetFillStyle(3003);
	// h_DsoverD0_DRatio_Sys->SetFillStyle(4095);
	h_DsoverD0_DRatio_Sys->Draw("SAME E2");
	// h_DsoverD0_DRatio_Sys->Draw("2SAME");
	h_DsoverD0_DRatio_Sta->SetLineColor(kBlue);
	h_DsoverD0_DRatio_Sta->SetMarkerColor(kBlue);
	h_DsoverD0_DRatio_Sta->Draw("SAME");
	// h_DsoverD0_DRatio_Sys->Draw("SAME E2");
	gr_PHSD_DsD0_PbPbpp->Draw("l");


  SavePlotDirs(c_DsD0_DRatio,"DsD0_DRatio",{"Result"});


	TCanvas *c_DsOverD0_withDRatio=new TCanvas("c_DsOverD0_withDRatio","c_DsOverD0_withDRatio",c_wtopx,c_wtopy,c_W,800);
	SetCanvas(c_DsOverD0_withDRatio);

//	TCanvas *c_CSdNdpt_ratio=new TCanvas("c_CSdNdpt_ratio","c_CSdNdpt_ratio",c_wtopx,c_wtopy,c_W,800);
/*	c_CSdNdpt_ratio->cd();
  c_CSdNdpt_ratio->SetFillColor(0);
  c_CSdNdpt_ratio->SetBorderMode(0);
  c_CSdNdpt_ratio->SetFrameFillStyle(0);
  c_CSdNdpt_ratio->SetFrameBorderMode(0);
  c_CSdNdpt_ratio->SetLeftMargin( c_Lmg );
  c_CSdNdpt_ratio->SetRightMargin( c_Rmg );
  c_CSdNdpt_ratio->SetTopMargin( c_Tmg );
  c_CSdNdpt_ratio->SetBottomMargin( c_Bmg );
  c_CSdNdpt_ratio->SetTickx(0);
  c_CSdNdpt_ratio->SetTicky(0);
*/
  TPad *pad1a = new TPad("pad1a","top pad",0.0,0.25,1.0,1.0);
  pad1a->SetTopMargin(0.08);
  pad1a->SetBottomMargin(0.0);
  pad1a->SetRightMargin(c_Rmg);
  pad1a->SetLeftMargin(c_Lmg);
  pad1a->Draw();
  TPad *pad2a = new TPad("pad2a","bottom pad",0.0,0.00,1.0,0.25);
  pad2a->SetTopMargin(0.0);
  pad2a->SetBottomMargin(0.30);
  pad2a->SetRightMargin(c_Rmg);
  pad2a->SetLeftMargin(c_Lmg);
  pad2a->Draw();

  pad1a->cd();
	// h_DsoverD0_pp_Sys->Draw("SAME E2");
	// h_DsoverD0_pp_Sta->SetLineColor(kRed);
	// h_DsoverD0_pp_Sta->SetMarkerColor(kRed);
	h_DsoverD0_pp_Sys->SetMinimum(0.02);
	h_DsoverD0_pp_Sys->GetYaxis()->SetTitleOffset(1.1);
	// h_DsoverD0_pp_Sys->SetMaximum(1);
	// h_DsoverD0_pp_Sys->SetTitle("");
	// h_DsoverD0_pp_Sys->GetYaxis()->SetTitle("D_{S}/D^{0} Ratio");
	// h_DsoverD0_pp_Sys->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	// h_DsoverD0_pp_Sys->GetXaxis()->CenterTitle();
	// h_DsoverD0_pp_Sys->GetYaxis()->CenterTitle();
	// h_DsoverD0_pp_Sys->SetFillColor(kRed-4);
	// h_DsoverD0_pp_Sys->SetLineColor(kRed);
	// h_DsoverD0_pp_Sys->SetMarkerColor(kRed);
	// h_DsoverD0_pp_Sys->SetFillStyle(3005);
	h_DsoverD0_pp_Sys->Draw("SAME E2");
	
	h_DsoverD0_pp_Sta->Draw("SAME");

	// h_DsoverD0_PbPb_Sys->SetLineColor(kBlue);
	// h_DsoverD0_PbPb_Sys->SetMarkerColor(kBlue);
	// h_DsoverD0_PbPb_Sys->SetFillColor(kBlue-4);
	// h_DsoverD0_PbPb_Sys->SetFillStyle(3004);
	h_DsoverD0_PbPb_Sys->Draw("SAME E2");

	// h_DsoverD0_PbPb_Sta->SetLineColor(kBlue);
	// h_DsoverD0_PbPb_Sta->SetMarkerColor(kBlue);
	h_DsoverD0_PbPb_Sta->Draw("SAME");

	// h_DsoverD0_Pythia->SetLineColor(1);
	// h_DsoverD0_Pythia->SetMarkerColor(1);
	// h_DsoverD0_Pythia->Draw("SAME");
	gr_DsD0_Pythia->Draw("l");

	// gr_TAMU_DsD0_pp->SetLineColor(kGreen+2);
	gr_TAMU_DsD0_pp->Draw("l");
	// gr_PHSD_DsD0_pp->SetLineColor(kMagenta+1);
	gr_PHSD_DsD0_pp->Draw("l");
	// gr_PHSD_DsD0_PbPb->SetLineColor(kViolet+2);
	gr_PHSD_DsD0_PbPb->Draw("l");



	le_DsD0->Draw("SAME");
	le_DsD02->Draw("SAME");
	texppPbPb->Draw("SAME");
	texCmsPre->Draw("SAME");
	
  // gStyle->SetOptStat(0);
	gPad->Modified();
	gPad->Update();

	pad2a->cd();
  h_DsoverD0_DRatio_Sys->GetXaxis()->SetTitleSize(0.14);
  h_DsoverD0_DRatio_Sys->GetXaxis()->SetLabelSize(0.12);
  h_DsoverD0_DRatio_Sys->GetYaxis()->SetTitleSize(0.1);
  h_DsoverD0_DRatio_Sys->GetYaxis()->SetLabelSize(0.08);
  h_DsoverD0_DRatio_Sys->GetYaxis()->SetTitleOffset(0.58);

	h_DsoverD0_DRatio_Sys->SetMaximum(1.9);
	h_DsoverD0_DRatio_Sys->SetMinimum(0.35);
	h_DsoverD0_DRatio_Sys->Draw("SAME E2");
	h_DsoverD0_DRatio_Sta->SetLineColor(kBlue);
	h_DsoverD0_DRatio_Sta->SetMarkerColor(kBlue);
	h_DsoverD0_DRatio_Sta->Draw("SAME");
	// gr_PHSD_DsD0_PbPbpp->SetLineWidth(3);
	// gr_PHSD_DsD0_PbPbpp->SetLineColor(kMagenta+1);
	gr_PHSD_DsD0_PbPbpp->Draw("l");

	SavePlotDirs(c_DsOverD0_withDRatio,"DsOverD0_withDRatio",{"Result"},"cf");


	// return;


	// TCanvas *c_RAA=new TCanvas("c_RAA","c_RAA",c_wtopx,c_wtopy,c_W,c_H);
	TCanvas *c_RAA=new TCanvas("c_RAA","c_RAA",c_wtopx,c_wtopy,c_W,800); // square
	c_RAA->cd();
  c_RAA->SetFillColor(0);
  c_RAA->SetBorderMode(0);
  c_RAA->SetFrameFillStyle(0);
  c_RAA->SetFrameBorderMode(0);
  c_RAA->SetLeftMargin( c_Lmg );
  c_RAA->SetRightMargin( c_Rmg );
  c_RAA->SetTopMargin( c_Tmg );
  c_RAA->SetBottomMargin( c_Bmg );
  c_RAA->SetTickx(0);
  c_RAA->SetTicky(0);

	h_RAA->SetTitle("");
	h_RAA->SetMaximum(1.35);
	h_RAA->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h_RAA->GetXaxis()->SetTitleSize(0.045);
	h_RAA->GetXaxis()->SetTitleOffset(1.2);
	h_RAA->GetYaxis()->SetTitleOffset(1.1);
	h_RAA->GetYaxis()->SetTitle("R_{AA}");
	h_RAA->GetXaxis()->CenterTitle();
	h_RAA->GetYaxis()->CenterTitle();
	h_RAA->SetLineColor(kBlue);
	h_RAA->SetMarkerColor(kBlue);
	h_RAA->SetFillColorAlpha(kBlue-9,0.4);
//	h_RAA->SetFillStyle(3004);
	h_RAA->Draw("SAME E2");
	h_RAA_Sta->SetLineColor(kBlue);
	h_RAA_Sta->SetMarkerColor(kBlue);
	h_RAA_Sta->Draw("SAME");

  texppPbPb->Draw("SAME");
	// texCmsPre->Draw("SAME");
	texCms2->Draw("SAME");
	// texPre2->Draw("SAME");

	// gr_RAA_PHSD->Draw("CP");
	gr_PHSD_Raa->SetLineWidth(3);
	gr_PHSD_Raa->SetLineColor(kMagenta+1);
	gr_PHSD_Raa->Draw("l");

  TLegend *le_RAA= new TLegend(0.65,0.50,0.85,0.7,NULL,"brNDC");
  le_RAA->SetBorderSize(0);
	le_RAA->SetTextSize(0.03);
  le_RAA->AddEntry(h_RAA,"R_{AA} 0-100%","flp");
	// le_RAA->AddEntry(gr_RAA_PHSD,"PHSD 0-80%","lp");
	le_RAA->AddEntry(gr_PHSD_Raa,"PHSD 0-80%","l");
	le_RAA->Draw("same");

	TGraphAsymmErrors *gr_RAA_GloSys=new TGraphAsymmErrors();
	gr_RAA_GloSys->SetPoint(0,5,1);
	gr_RAA_GloSys->SetPointError(0,1,1,0.04566,0.04139);
	gr_RAA_GloSys->SetFillColorAlpha(kGray+2,0.4);
	gr_RAA_GloSys->Draw("F2");

	TLine *l_Raa1=new TLine(4,1,40,1);
	l_Raa1->SetLineStyle(7);
	l_Raa1->Draw();

	textemp->DrawLatexNDC(0.68,0.85,"D_{S}^{+}+D_{S}^{-}");

	
	SavePlotDirs(c_RAA,"RAA",{"Result"},"cf");

	// return;


	TCanvas *c_CSdNdpt=new TCanvas("c_CSdNdpt","c_CSdNdpt",c_wtopx,c_wtopy,c_W,c_H);
	c_CSdNdpt->cd();
  c_CSdNdpt->SetFillColor(0);
  c_CSdNdpt->SetBorderMode(0);
  c_CSdNdpt->SetFrameFillStyle(0);
  c_CSdNdpt->SetFrameBorderMode(0);
  c_CSdNdpt->SetLeftMargin( c_Lmg );
  c_CSdNdpt->SetRightMargin( c_Rmg );
  c_CSdNdpt->SetTopMargin( c_Tmg );
  c_CSdNdpt->SetBottomMargin( c_Bmg );
  c_CSdNdpt->SetTickx(0);
  c_CSdNdpt->SetTicky(0);

	// gStyle->SetOptStat(0);
	gPad->SetLogy();
	// gStyle->SetOptStat(0);

	double CSMax=h_PromptDs_CS_pp_Sys->GetMaximum();
	if(h_Genpt_Ds->GetMaximum()>CSMax){
		CSMax=h_Genpt_Ds->GetMaximum();
	}


	// h_PromptDs_CS_pp_Sys->SetLineColor(2);
	h_PromptDs_CS_pp_Sys->SetMaximum(CSMax*3);
	h_PromptDs_CS_pp_Sys->SetMinimum(100);
	h_PromptDs_CS_pp_Sys->SetTitle("");
	h_PromptDs_CS_pp_Sys->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{T}} (or #frac{1}{T_{AA}}#frac{dN}{dp_{T}}) (pb GeV^{-1})");
	h_PromptDs_CS_pp_Sys->GetYaxis()->CenterTitle();
	h_PromptDs_CS_pp_Sys->GetYaxis()->SetTitleOffset(1.45);
  h_PromptDs_CS_pp_Sys->GetYaxis()->SetTitleSize(0.045);
	h_PromptDs_CS_pp_Sys->GetXaxis()->CenterTitle();
	h_PromptDs_CS_pp_Sys->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h_PromptDs_CS_pp_Sys->SetFillColorAlpha(kRed-9,0.4);
	h_PromptDs_CS_pp_Sys->SetLineColor(kRed);
	h_PromptDs_CS_pp_Sys->SetMarkerColor(kRed);
	// h_PromptDs_CS_pp_Sys->SetFillStyle(3005);
	h_PromptDs_CS_pp_Sys->Draw("SAME E2");
	h_PromptDs_CS_pp_Sta->SetLineColor(kRed);
	h_PromptDs_CS_pp_Sta->SetMarkerColor(kRed);
	h_PromptDs_CS_pp_Sta->Draw("SAME");

	h_PromptDs_dNdpt_PbPb3_Sys->SetLineColor(kBlue);
	h_PromptDs_dNdpt_PbPb3_Sys->SetMarkerColor(kBlue);
	h_PromptDs_dNdpt_PbPb3_Sys->SetFillColorAlpha(kBlue-9,0.4);
	// h_PromptDs_dNdpt_PbPb3_Sys->SetFillStyle(3004);
	h_PromptDs_dNdpt_PbPb3_Sys->Draw("SAME E2");

	h_PromptDs_dNdpt_PbPb3_Sta->SetLineColor(kBlue);
	h_PromptDs_dNdpt_PbPb3_Sta->SetMarkerColor(kBlue);
	h_PromptDs_dNdpt_PbPb3_Sta->Draw("SAME");

	// h_Genpt_Ds->SetLineColor(1);
	// h_Genpt_Ds->SetMarkerColor(1);
	// h_Genpt_Ds->Draw("SAME");

	gr_Genpt_Ds_Pythia->SetLineColor(1);
	gr_Genpt_Ds_Pythia->SetLineWidth(3);
	gr_Genpt_Ds_Pythia->Draw("l");


	TLegend *le= new TLegend(0.65,0.53,0.82,0.77,NULL,"brNDC");
	le->SetBorderSize(0);
	le->SetTextSize(0.035);
	le->AddEntry(h_PromptDs_CS_pp_Sys,"pp","flp");
	le->AddEntry(h_PromptDs_dNdpt_PbPb3_Sys,"PbPb 0-100%","flp");
	// le->AddEntry(h_Genpt_Ds,"PYTHIA 8","lp");
	le->AddEntry(gr_Genpt_Ds_Pythia,"PYTHIA 8","l");
	le->Draw("SAME");
  texppPbPb->Draw("SAME");
	texCmsPre->Draw("SAME");
	
	textemp->DrawLatexNDC(0.68,0.85,"#frac{D_{S}^{+}+D_{S}^{-}}{2}");
  // gStyle->SetOptStat(0);
	gPad->Modified();
	gPad->Update();

	SavePlotDirs(c_CSdNdpt,"CSdNdpt",{"Result"},"cf");

	TCanvas *c_CSdNdpt_ratio=new TCanvas("c_CSdNdpt_ratio","c_CSdNdpt_ratio",c_wtopx,c_wtopy,c_W,800);
	c_CSdNdpt_ratio->cd();
  c_CSdNdpt_ratio->SetFillColor(0);
  c_CSdNdpt_ratio->SetBorderMode(0);
  c_CSdNdpt_ratio->SetFrameFillStyle(0);
  c_CSdNdpt_ratio->SetFrameBorderMode(0);
  c_CSdNdpt_ratio->SetLeftMargin( c_Lmg );
  c_CSdNdpt_ratio->SetRightMargin( c_Rmg );
  c_CSdNdpt_ratio->SetTopMargin( c_Tmg );
  c_CSdNdpt_ratio->SetBottomMargin( c_Bmg );
  c_CSdNdpt_ratio->SetTickx(0);
  c_CSdNdpt_ratio->SetTicky(0);

  TPad *pad1 = new TPad("pad1","top pad",0.0,0.25,1.0,1.0);
  pad1->SetTopMargin(0.08);
  pad1->SetBottomMargin(0.0);
  pad1->SetRightMargin(c_Rmg);
  pad1->SetLeftMargin(c_Lmg);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","bottom pad",0.0,0.00,1.0,0.25);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.30);
  pad2->SetRightMargin(c_Rmg);
  pad2->SetLeftMargin(c_Lmg);
  pad2->Draw();

  pad1->cd();
	h_PromptDs_CS_pp_Sys->Draw("SAME E2");
	h_PromptDs_CS_pp_Sta->Draw("SAME");

	h_PromptDs_dNdpt_PbPb3_Sys->Draw("SAME E2");

	h_PromptDs_dNdpt_PbPb3_Sta->Draw("SAME");

	// h_Genpt_Ds->Draw("SAME");
	gr_Genpt_Ds_Pythia->Draw("l");
	le->Draw("SAME");
  texppPbPb->Draw("SAME");
	texCmsPre->Draw("SAME");
	
	textemp->DrawLatexNDC(0.68,0.8,"#frac{D_{S}^{+}+D_{S}^{-}}{2}");
	// textemp->DrawLatexNDC(0.68,0.85,"D_{S}^{#pm}");
  // gStyle->SetOptStat(0);
	gPad->Modified();
	gPad->Update();

	// gStyle->SetOptStat(0);
	gPad->SetLogy();


	pad2->cd();

	TH1D *h_Genpt_Ds_Ratio_pp_Sys=(TH1D*)h_PromptDs_CS_pp_Sys->Clone("h_Genpt_Ds_Ratio_pp_Sys");
	TH1D *h_Genpt_Ds_Ratio_pp_Sta=(TH1D*)h_PromptDs_CS_pp_Sta->Clone("h_Genpt_Ds_Ratio_pp_Sta");

	TH1D *h_Genpt_Ds_Ratio_PbPb3_Sys=(TH1D*)h_PromptDs_dNdpt_PbPb3_Sys->Clone("h_Genpt_Ds_Ratio_PbPb3_Sys");
	TH1D *h_Genpt_Ds_Ratio_PbPb3_Sta=(TH1D*)h_PromptDs_dNdpt_PbPb3_Sta->Clone("h_Genpt_Ds_Ratio_PbPb3_Sta");

	TH1D *h_Genpt_Ds_temp_pp=new TH1D("h_Genpt_Ds_temp_pp","h_Genpt_Ds_temp_pp",nbin_pt_pp,bins_pt_pp); h_Genpt_Ds_temp_pp->Sumw2();
	TH1D *h_Genpt_Ds_temp_PbPb3=new TH1D("h_Genpt_Ds_temp_PbPb3","h_Genpt_Ds_temp_PbPb3",nbin_pt_PbPb3,bins_pt_PbPb3);
	h_Genpt_Ds_temp_PbPb3->Sumw2();

	for(int i=0; i<nbin_pt_pp; i++){
		h_Genpt_Ds_temp_pp->SetBinContent(i+1, h_Genpt_Ds->GetBinContent(h_Genpt_Ds->FindBin(h_Genpt_Ds_temp_pp->GetBinCenter(i+1))));

		cout<<"bin center = "<<h_Genpt_Ds_temp_pp->GetBinCenter(i+1)<<" , bincontent = "<<h_Genpt_Ds->GetBinContent(h_Genpt_Ds->FindBin(h_Genpt_Ds_temp_pp->GetBinCenter(i+1))) <<endl;

	}
	for(int i=2; i<nbin_pt_PbPb3;i++){
		h_Genpt_Ds_temp_PbPb3->SetBinContent(i+1, h_Genpt_Ds->GetBinContent(h_Genpt_Ds->FindBin(h_Genpt_Ds_temp_PbPb3->GetBinCenter(i+1))));
	}

	h_Genpt_Ds_Ratio_pp_Sys->Divide(h_Genpt_Ds_temp_pp);
	h_Genpt_Ds_Ratio_pp_Sta->Divide(h_Genpt_Ds_temp_pp);

	for(int i =0; i<nbin_pt_pp; i++){
		cout<<"bin = "<<i<<" , h_Genpt_Ds_Ratio_pp_Sys->value "<<h_Genpt_Ds_Ratio_pp_Sys->GetBinContent(i+1)<<" , h_PromptDs_CS_pp_Sys-val = "<<h_PromptDs_CS_pp_Sys->GetBinContent(i+1)<<endl;
	}

	h_Genpt_Ds_Ratio_PbPb3_Sys->Divide(h_Genpt_Ds_temp_PbPb3);
	h_Genpt_Ds_Ratio_PbPb3_Sta->Divide(h_Genpt_Ds_temp_PbPb3);

	// h_Genpt_Ds_temp_pp->Draw();

	h_Genpt_Ds_Ratio_pp_Sys->SetMaximum(1.95);
	h_Genpt_Ds_Ratio_pp_Sys->SetMinimum(0);
	h_Genpt_Ds_Ratio_pp_Sys->GetYaxis()->SetTitle("#frac{Data}{PYTHIA}");
	h_Genpt_Ds_Ratio_pp_Sys->GetXaxis()->SetTitleSize(0.13);
	h_Genpt_Ds_Ratio_pp_Sys->GetXaxis()->SetLabelSize(0.09);
	h_Genpt_Ds_Ratio_pp_Sys->GetYaxis()->SetTitleOffset(0.45);
	h_Genpt_Ds_Ratio_pp_Sys->GetYaxis()->SetTitleSize(0.1);
	h_Genpt_Ds_Ratio_pp_Sys->GetYaxis()->SetLabelSize(0.08);
	h_Genpt_Ds_Ratio_pp_Sys->GetYaxis()->SetTitleOffset(0.58);
	h_Genpt_Ds_Ratio_pp_Sys->Draw("SAMEE2");
	h_Genpt_Ds_Ratio_pp_Sta->Draw("SAME");	
	h_Genpt_Ds_Ratio_PbPb3_Sys->Draw("SAMEE2");
	h_Genpt_Ds_Ratio_PbPb3_Sta->Draw("SAME");
	// h_Genpt_Ds_Ratio_PbPb3_Sys->Draw("SAME");


	SavePlotDirs(c_CSdNdpt_ratio,"CSdNdpt_ratio",{"Result"},"cf");


// print out Raa value

	for (int i=1; i<=h_RAA->GetNbinsX(); i++){

		cout<<"i = "<<i<<" ,h_RAA = "<< h_RAA->GetBinContent(i)<<" , stat = "<<h_RAA_Sta->GetBinError(i)/h_RAA->GetBinContent(i)*100 <<" %, syst = "<<h_RAA->GetBinError(i)/h_RAA->GetBinContent(i)*100<<" %"<<endl;

	}


	// plot Alice comparison

	TFile *f_Alice=TFile::Open("./output/Alice_result.root","read");

	TGraphAsymmErrors *gr_DsD0_0to10=(TGraphAsymmErrors*)f_Alice->Get("gr_DsD0_0to10");
	TGraphAsymmErrors *gr_DsD0_pp=(TGraphAsymmErrors*)f_Alice->Get("gr_DsD0_pp");
	TGraphAsymmErrors *gr_Raa=(TGraphAsymmErrors*)f_Alice->Get("gr_Raa");

	TCanvas *c_DsOverD0_Alice=new TCanvas("c_DsOverD0_Alice","c_DsOverD0_Alice",c_wtopx,c_wtopy,c_W,c_H);
	c_DsOverD0_Alice->cd();
	SetCanvas(c_DsOverD0_Alice);

// gStyle->SetOptStat(0);
	// gPad->SetLogy();
	// gStyle->SetOptStat(0);

	// h_DsoverD0_pp_Sys->SetLineColor(2);
	h_DsoverD0_pp_Sys->SetMinimum(0);
	h_DsoverD0_pp_Sys->SetMaximum(1);
	h_DsoverD0_pp_Sys->SetTitle("");
	h_DsoverD0_pp_Sys->GetYaxis()->SetTitle("D_{S}/D^{0} Ratio");
	h_DsoverD0_pp_Sys->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h_DsoverD0_pp_Sys->GetXaxis()->CenterTitle();
	h_DsoverD0_pp_Sys->GetYaxis()->CenterTitle();
	h_DsoverD0_pp_Sys->SetFillColor(kRed-4);
	h_DsoverD0_pp_Sys->SetLineColor(kRed);
	h_DsoverD0_pp_Sys->SetMarkerColor(kRed);
	h_DsoverD0_pp_Sys->SetFillStyle(3005);
	// h_DsoverD0_pp_Sys->SetFillStyle(4020);
	h_DsoverD0_pp_Sys->Draw("SAME E2");
	h_DsoverD0_pp_Sta->SetLineColor(kRed);
	h_DsoverD0_pp_Sta->SetMarkerColor(kRed);
	h_DsoverD0_pp_Sta->Draw("SAME");

	h_DsoverD0_PbPb_Sys->SetLineColor(kBlue);
	h_DsoverD0_PbPb_Sys->SetMarkerColor(kBlue);
	h_DsoverD0_PbPb_Sys->SetFillColor(kBlue-4);
	h_DsoverD0_PbPb_Sys->SetFillStyle(3004);
	// h_DsoverD0_PbPb_Sys->SetFillStyle(4020);
	h_DsoverD0_PbPb_Sys->Draw("SAME E2");

	h_DsoverD0_PbPb_Sta->SetLineColor(kBlue);
	h_DsoverD0_PbPb_Sta->SetMarkerColor(kBlue);
	h_DsoverD0_PbPb_Sta->Draw("SAME");

	// h_DsoverD0_Pythia->SetLineColor(1);
	// h_DsoverD0_Pythia->SetMarkerColor(1);
	// h_DsoverD0_Pythia->Draw("SAME");

	gr_DsD0_Pythia->Draw("l");

/*
	gr_DsD0_0to10->SetFillColor(kGreen-3);
	gr_DsD0_0to10->SetFillStyle(3005);
	gr_DsD0_0to10->SetLineColor(kGreen+2);
	gr_DsD0_0to10->SetMarkerColor(kGreen+2);
	gr_DsD0_0to10->SetMarkerStyle(28);	
*/
	gr_DsD0_0to10->SetFillColor(13);
	gr_DsD0_0to10->SetFillStyle(3003);
	gr_DsD0_0to10->SetLineColor(1);
	gr_DsD0_0to10->SetMarkerColor(1);
	gr_DsD0_0to10->SetMarkerStyle(24);	

	gr_DsD0_0to10->Draw("p2same");

/*
	gr_DsD0_pp->SetFillColor(kMagenta-6);
	gr_DsD0_pp->SetFillStyle(3004);
	gr_DsD0_pp->SetLineColor(kMagenta+3);
	gr_DsD0_pp->SetMarkerColor(kMagenta+3);
	gr_DsD0_pp->SetMarkerStyle(28);	
*/
	gr_DsD0_pp->SetFillColor(14);
	gr_DsD0_pp->SetFillStyle(3003);
	gr_DsD0_pp->SetLineColor(1);
	gr_DsD0_pp->SetMarkerColor(1);
	gr_DsD0_pp->SetMarkerStyle(32);	


	gr_DsD0_pp->Draw("p2same");


	TLegend *le_DsD0Alice= new TLegend(0.65,0.62,0.85,0.85,NULL,"brNDC");
	le_DsD0Alice->SetBorderSize(0);
	le_DsD0Alice->AddEntry(h_DsoverD0_pp_Sys,"pp","flp");
	le_DsD0Alice->AddEntry(h_DsoverD0_PbPb_Sys,"PbPb 0-100%","flp");
	// le_DsD0Alice->AddEntry(h_DsoverD0_Pythia,"PYTHIA 8","lp");
	le_DsD0Alice->AddEntry(gr_DsD0_Pythia,"PYTHIA 8","l");
	le_DsD0Alice->Draw("SAME");
	texppPbPb->Draw("SAME");
	texCmsPre->Draw("SAME");

	TLegend *le_DsD0Alice2=new TLegend(0.30,0.62,0.50,0.85);
	le_DsD0Alice2->SetBorderSize(0);
	le_DsD0Alice2->AddEntry(gr_DsD0_0to10,"ALICE PbPb 0-10%","flp");
	le_DsD0Alice2->AddEntry(gr_DsD0_pp,"ALICE pp 7TeV","flp");
	le_DsD0Alice2->Draw("SAME");

	
  // gStyle->SetOptStat(0);
	gPad->Modified();
	gPad->Update();


	SavePlotDirs(c_DsOverD0_Alice,"DsOverD0_Alice",{"Result"});


	TCanvas *c_DsOverD0_withDRatio_Alice=new TCanvas("c_DsOverD0_withDRatio_Alice","c_DsOverD0_withDRatio_Alice",c_wtopx,c_wtopy,c_W,800);
	SetCanvas(c_DsOverD0_withDRatio_Alice);

//	TCanvas *c_CSdNdpt_ratio=new TCanvas("c_CSdNdpt_ratio","c_CSdNdpt_ratio",c_wtopx,c_wtopy,c_W,800);
/*	c_CSdNdpt_ratio->cd();
  c_CSdNdpt_ratio->SetFillColor(0);
  c_CSdNdpt_ratio->SetBorderMode(0);
  c_CSdNdpt_ratio->SetFrameFillStyle(0);
  c_CSdNdpt_ratio->SetFrameBorderMode(0);
  c_CSdNdpt_ratio->SetLeftMargin( c_Lmg );
  c_CSdNdpt_ratio->SetRightMargin( c_Rmg );
  c_CSdNdpt_ratio->SetTopMargin( c_Tmg );
  c_CSdNdpt_ratio->SetBottomMargin( c_Bmg );
  c_CSdNdpt_ratio->SetTickx(0);
  c_CSdNdpt_ratio->SetTicky(0);
*/
  TPad *pad1b = new TPad("pad1b","top pad",0.0,0.25,1.0,1.0);
  pad1b->SetTopMargin(0.08);
  pad1b->SetBottomMargin(0.0);
  pad1b->SetRightMargin(c_Rmg);
  pad1b->SetLeftMargin(c_Lmg);
  pad1b->Draw();
  TPad *pad2b = new TPad("pad2b","bottom pad",0.0,0.00,1.0,0.25);
  pad2b->SetTopMargin(0.0);
  pad2b->SetBottomMargin(0.30);
  pad2b->SetRightMargin(c_Rmg);
  pad2b->SetLeftMargin(c_Lmg);
  pad2b->Draw();

  pad1b->cd();
	// h_DsoverD0_pp_Sys->Draw("SAME E2");
	// h_DsoverD0_pp_Sta->SetLineColor(kRed);
	// h_DsoverD0_pp_Sta->SetMarkerColor(kRed);
	h_DsoverD0_pp_Sys->SetMinimum(0.02);
	// h_DsoverD0_pp_Sys->SetMaximum(1);
	// h_DsoverD0_pp_Sys->SetTitle("");
	// h_DsoverD0_pp_Sys->GetYaxis()->SetTitle("D_{S}/D^{0} Ratio");
	// h_DsoverD0_pp_Sys->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	// h_DsoverD0_pp_Sys->GetXaxis()->CenterTitle();
	// h_DsoverD0_pp_Sys->GetYaxis()->CenterTitle();
	// h_DsoverD0_pp_Sys->SetFillColor(kRed-4);
	// h_DsoverD0_pp_Sys->SetLineColor(kRed);
	// h_DsoverD0_pp_Sys->SetMarkerColor(kRed);
	// h_DsoverD0_pp_Sys->SetFillStyle(3005);
	h_DsoverD0_pp_Sys->Draw("SAME E2");
	
	h_DsoverD0_pp_Sta->Draw("SAME");

	// h_DsoverD0_PbPb_Sys->SetLineColor(kBlue);
	// h_DsoverD0_PbPb_Sys->SetMarkerColor(kBlue);
	// h_DsoverD0_PbPb_Sys->SetFillColor(kBlue-4);
	// h_DsoverD0_PbPb_Sys->SetFillStyle(3004);
	h_DsoverD0_PbPb_Sys->Draw("SAME E2");

	// h_DsoverD0_PbPb_Sta->SetLineColor(kBlue);
	// h_DsoverD0_PbPb_Sta->SetMarkerColor(kBlue);
	h_DsoverD0_PbPb_Sta->Draw("SAME");

	// h_DsoverD0_Pythia->SetLineColor(1);
	// h_DsoverD0_Pythia->SetMarkerColor(1);
	// h_DsoverD0_Pythia->Draw("SAME");
	gr_DsD0_Pythia->Draw("l");

	le_DsD0->Draw("SAME");
	texppPbPb->Draw("SAME");
	texCmsPre->Draw("SAME");

	gr_DsD0_0to10->Draw("p2same");

	gr_DsD0_pp->Draw("p2same");

	// le_DsD0Alice2->Draw("SAME");
/*
	TLegend *le_DsD0Alice= new TLegend(0.65,0.62,0.85,0.85,NULL,"brNDC");
	le_DsD0Alice->SetBorderSize(0);
	le_DsD0Alice->AddEntry(h_DsoverD0_pp_Sys,"pp","flp");
	le_DsD0Alice->AddEntry(h_DsoverD0_PbPb_Sys,"PbPb 0-100%","flp");
	le_DsD0Alice->AddEntry(h_DsoverD0_Pythia,"PYTHIA 8","lp");
*/
//	le_DsD0Alice->Draw("SAME");
	// texppPbPb->Draw("SAME");
	// texCmsPre->Draw("SAME");

	TLegend *le_DsD0Alice3=new TLegend(0.30,0.62,0.50,0.85,NULL,"brNDC");
	le_DsD0Alice3->SetBorderSize(0);
	le_DsD0Alice3->AddEntry(gr_DsD0_0to10,"ALICE PbPb 0-10%","flp");
	le_DsD0Alice3->AddEntry(gr_DsD0_pp,"ALICE pp 7TeV","flp");
	le_DsD0Alice3->Draw("SAME");


	
  // gStyle->SetOptStat(0);
	gPad->Modified();
	gPad->Update();

	pad2b->cd();
  h_DsoverD0_DRatio_Sys->GetXaxis()->SetTitleSize(0.15);
  h_DsoverD0_DRatio_Sys->GetXaxis()->SetLabelSize(0.12);
  h_DsoverD0_DRatio_Sys->GetYaxis()->SetTitleSize(0.1);
  h_DsoverD0_DRatio_Sys->GetYaxis()->SetLabelSize(0.08);
  h_DsoverD0_DRatio_Sys->GetYaxis()->SetTitleOffset(0.58);

	h_DsoverD0_DRatio_Sys->Draw("SAME E2");
	h_DsoverD0_DRatio_Sta->SetLineColor(kBlue);
	h_DsoverD0_DRatio_Sta->SetMarkerColor(kBlue);
	h_DsoverD0_DRatio_Sta->Draw("SAME");

	SavePlotDirs(c_DsOverD0_withDRatio_Alice,"DsOverD0_withDRatioAlice",{"Result"});






	TCanvas *c_RAA_Alice=new TCanvas("c_RAA_Alice","c_RAA_Alice",c_wtopx,c_wtopy,c_W,c_H);
	c_RAA_Alice->cd();
  c_RAA_Alice->SetFillColor(0);
  c_RAA_Alice->SetBorderMode(0);
  c_RAA_Alice->SetFrameFillStyle(0);
  c_RAA_Alice->SetFrameBorderMode(0);
  c_RAA_Alice->SetLeftMargin( c_Lmg );
  c_RAA_Alice->SetRightMargin( c_Rmg );
  c_RAA_Alice->SetTopMargin( c_Tmg );
  c_RAA_Alice->SetBottomMargin( c_Bmg );
  c_RAA_Alice->SetTickx(0);
  c_RAA_Alice->SetTicky(0);


	h_RAA->SetTitle("");
	h_RAA->SetMaximum(1);
	h_RAA->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h_RAA->GetYaxis()->SetTitle("R_{AA}");
	h_RAA->GetXaxis()->CenterTitle();
	h_RAA->GetYaxis()->CenterTitle();
	h_RAA->SetLineColor(kBlue);
	h_RAA->SetMarkerColor(kBlue);
	h_RAA->SetFillColor(kBlue-9);
	h_RAA->SetFillStyle(3004);
	h_RAA->Draw("SAME E2");
	h_RAA_Sta->SetLineColor(kBlue);
	h_RAA_Sta->SetMarkerColor(kBlue);
	h_RAA_Sta->Draw("SAME");
  texppPbPb->Draw("SAME");
	texCmsPre->Draw("SAME");
	// gr_RAA_PHSD->Draw("CP");
	gr_PHSD_Raa->Draw("l");
/*
	gr_Raa->SetFillColor(kMagenta-6);
	gr_Raa->SetFillStyle(3004);
	gr_Raa->SetLineColor(kMagenta+3);
	gr_Raa->SetMarkerColor(kMagenta+3);
*/
	gr_Raa->SetFillColor(13);
	gr_Raa->SetFillStyle(3003);
	gr_Raa->SetLineColor(1);
	gr_Raa->SetMarkerColor(1);

	gr_Raa->SetMarkerStyle(24);	

	gr_Raa->Draw("p2same");


  TLegend *le_RAA_Alice= new TLegend(0.65,0.66,0.82,0.82,NULL,"brNDC");
  le_RAA_Alice->SetBorderSize(0);
  le_RAA_Alice->AddEntry(h_RAA,"R_{AA} 0-100%","flp");
	// le_RAA_Alice->AddEntry(gr_RAA_PHSD,"PHSD 0-80%","lp");
	le_RAA_Alice->AddEntry(gr_PHSD_Raa,"PHSD 0-80%","lp");
	le_RAA_Alice->AddEntry(gr_Raa,"ALICE R_{AA} 0-10%","flp");
	le_RAA_Alice->Draw("same");

	textemp->DrawLatexNDC(0.68,0.85,"D_{S}^{#pm}");

	
	SavePlotDirs(c_RAA_Alice,"RAA_Alice",{"Result"});


	// write histogram/graph into output

	fout->cd();
	h_RAA->Write();
	h_RAA_Sta->Write();
	
	h_DsoverD0_pp_Sta->Write();
	h_DsoverD0_pp_Sys->Write();
	
	h_DsoverD0_PbPb_Sys->Write();
	h_DsoverD0_PbPb_Sta->Write();

	h_PromptDs_CS_pp->Write();
	h_PromptDs_CS_pp_Sys->Write();
	h_PromptDs_CS_pp_Sta->Write();

	h_PromptDs_dNdpt_PbPb3->Write();
	h_PromptDs_dNdpt_PbPb3_Sys->Write();
	h_PromptDs_dNdpt_PbPb3_Sta->Write();

	gr_Raa->Write();
	gr_DsD0_pp->Write();
	gr_DsD0_0to10->Write();
	gr_RAA_GloSys->Write();
	// gr_RAA_PHSD->Write();
	// gr_DsD0->Write();
	gr_TAMU_DsD0_pp->Write();

/*

	TH1D *h_RAA_statonly=new TH1D("h_RAA_statonly","h_RAA_statonly",nbin_pt_PbPb3,bins_pt_PbPb3);	h_RAA_statonly->Sumw2();

	for(int i=StartBin; i<nbin_pt_PbPb3; i++){
		double ibin_pp=h_PromptDs_CS_pp->FindBin(h_PromptDs_dNdpt_PbPb3->GetBinCenter(i+1));
		double dNdpt=h_PromptDs_dNdpt_PbPb3->GetBinContent(i+1);		
		double dNdptErr=h_PromptDs_dNdpt_PbPb3->GetBinError(i+1);	
		double CS=h_PromptDs_CS_pp->GetBinContent(ibin_pp);
		double CSErr=h_PromptDs_CS_pp->GetBinError(ibin_pp);

		double RAA=1/TAA0to100*dNdpt/CS;
		double RAAErr=1/TAA0to100*ErrorPro_AoverB(dNdptErr,CSErr,dNdpt,CS);

		h_RAA_statonly->SetBinContent(i+1,RAA);
		h_RAA_statonly->SetBinError(i+1,RAAErr);

		if(verbose){
			cout<<"\nbin i ="<<i<<" , Dpt : "<<bins_pt_PbPb3[i]<<" to "<<bins_pt_PbPb3[i+1]<<endl;
			cout<<"RAA = "<<RAA<<" +- "<<RAAErr<<endl;
		}

	} // end for

	fout->cd();
	h_RAA_statonly->Write("",TObject::kOverwrite);

// plot?
	TCanvas *c_RAA_statonly=new TCanvas("c_RAA_statonly","c_RAA_statonly",800,800);
	c_RAA_statonly->cd();
	h_RAA_statonly->Draw();
	gPad->SetLogx();
	gPad->BuildLegend();
	SavePlotDirs(c_RAA_statonly,"RAA_statonly",{"Result"});
*/

}
