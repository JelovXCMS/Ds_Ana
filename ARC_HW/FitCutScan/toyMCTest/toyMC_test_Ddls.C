#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting_simple.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/CMS_lumi.C"

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

#include "TFitter.h"
#include "TFitResult.h"
#include "TRandom3.h"
using namespace RooFit;
using namespace std;

int toyMC_test_Ddls(){

	TRandom3 Rdm;
	Rdm.SetSeed(0);


	TFile *fout=new TFile("test.root","recreate");
	fout->cd();
 
	TTree *nt=new TTree("nt","nt");


  Float_t Ddl_True;
  Float_t Ddl_FromErr2;
  Float_t Ddls_FromErr2;
  Float_t Ddl_FromErr3;
  Float_t Ddls_FromErr3;

	Float_t Ddl_True_exp;
	Float_t Ddl_exp_res1;
	Float_t Ddl_exp_res3;

  Float_t Ddl_FromErr2Gaus1;
  Float_t Ddls_FromErr2Gaus1;
  Float_t Ddl_FromErr3Gaus1;
  Float_t Ddls_FromErr3Gaus1;

  Float_t Ddl_FromErr1;
  Float_t Ddl_FromErr1_smear1;
  Float_t Ddls_FromErr1;
  
  // Float_t Ddl_FromErr2;
  Float_t Ddls_2to3;
  Float_t Ddls_2to3Gaus1;
  Float_t Ddls_1to3;
  Float_t Ddls_1to3_noSm;

	nt->Branch("Ddl_True",&Ddl_True);
	nt->Branch("Ddl_True_exp",&Ddl_True_exp);
	nt->Branch("Ddl_exp_res1",&Ddl_exp_res1);
	nt->Branch("Ddl_exp_res3",&Ddl_exp_res3);



  nt->Branch("Ddl_FromErr1",&Ddl_FromErr1);
  nt->Branch("Ddl_FromErr1_smear1",&Ddl_FromErr1_smear1);
  nt->Branch("Ddls_FromErr1",&Ddls_FromErr1);
  nt->Branch("Ddl_FromErr2",&Ddl_FromErr2);
  nt->Branch("Ddls_FromErr2",&Ddls_FromErr2);
  nt->Branch("Ddl_FromErr3",&Ddl_FromErr3);
  nt->Branch("Ddls_FromErr3",&Ddls_FromErr3);
  nt->Branch("Ddls_1to3",&Ddls_1to3);
  nt->Branch("Ddls_1to3_noSm",&Ddls_1to3_noSm);
  nt->Branch("Ddls_2to3",&Ddls_2to3);

	int nEntry=5000;

	for(int i =0; i<nEntry; i++){
		Ddl_True=10;
	
		Ddl_True_exp=Rdm.Exp(10);
		Ddl_exp_res1=Ddl_True+Rdm.Gaus(0,1);
		Ddl_exp_res3=Ddl_True+Rdm.Gaus(0,3);

	
		// Ddl_FromErr1=(Ddl_True+Rdm.Gaus(0,1))/Rdm.Gaus(1,0.1);
		Ddl_FromErr1=(Ddl_True+Rdm.Gaus(0,1));
		Ddl_FromErr1_smear1=Rdm.Gaus(Ddl_FromErr1,2);
		Ddl_FromErr2=Ddl_True+Rdm.Gaus(0,2);
		Ddl_FromErr3=Ddl_True+Rdm.Gaus(0,3);
		Ddls_FromErr1=Ddl_FromErr1/1;
		Ddls_FromErr2=Ddl_FromErr2/2;
		Ddls_FromErr3=Ddl_FromErr3/3;
		// Ddls_FromErr1=Ddl_FromErr1/Rdm.Gaus(1,0.2);
		// Ddls_FromErr2=Ddl_FromErr2/Rdm.Gaus(2,0.4);
		// Ddls_FromErr3=Ddl_FromErr3/Rdm.Gaus(3,0.6);

		Ddls_1to3=Ddls_FromErr1*1/3;
		Ddls_1to3_noSm=Ddls_1to3;
		Ddls_1to3=Rdm.Gaus(Ddls_1to3,sqrt(1.0*1.0-1.0/3.0*1.0/3.0));

		Ddls_2to3=Ddls_FromErr2*2/3;
		Ddls_2to3=Rdm.Gaus(Ddls_2to3,sqrt(1.0*1.0-2.0/3.0*2.0/3.0));

		nt->Fill();

	}


	nt->Draw("Ddl_True_exp");

	return 1;

	TH1D *h_Ddl_Err1=new TH1D("h_Ddl_Err1","",200,0,20); h_Ddl_Err1->Sumw2();
	TH1D *h_Ddl_Err1_Smear1=new TH1D("h_Ddl_Err1_Smear1","",200,0,20); h_Ddl_Err1_Smear1->Sumw2();
	TH1D *h_Ddl_Err2=new TH1D("h_Ddl_Err2","",200,0,20); h_Ddl_Err2->Sumw2();
	TH1D *h_Ddl_Err3=new TH1D("h_Ddl_Err3","",200,0,20); h_Ddl_Err3->Sumw2();

	nt->Project("h_Ddl_Err1","Ddl_FromErr1");
	nt->Project("h_Ddl_Err1_Smear1","Ddl_FromErr1_smear1");
	nt->Project("h_Ddl_Err2","Ddl_FromErr2");
	nt->Project("h_Ddl_Err3","Ddl_FromErr3");

	TF1 *f1_Gaus= new TF1("f1_Gaus","[0]*TMath::Gaus(x,[1],[2]) /( sqrt(2*TMath::Pi())*[2] )");
	f1_Gaus->SetParameter(0,nEntry);
	f1_Gaus->SetParameter(1,10);
	f1_Gaus->SetParameter(2,1);
	f1_Gaus->SetLineColor(2);
	f1_Gaus->SetRange(0,20);

	h_Ddl_Err1->Fit("f1_Gaus","MIN0");
	h_Ddl_Err1->Fit("f1_Gaus","MIN0");
	h_Ddl_Err1->Fit("f1_Gaus","LMIN0");
	h_Ddl_Err1->Fit("f1_Gaus","LMIN0");
	h_Ddl_Err1->Fit("f1_Gaus","LMIN0");
	h_Ddl_Err1->Fit("f1_Gaus","LMIN0");
	h_Ddl_Err1->Fit("f1_Gaus","LMIN0");

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1111);
	TCanvas *c_Ddl=new TCanvas("c_Ddl","c_Ddl",600,600);
	c_Ddl->cd();
	h_Ddl_Err1->Draw();
	f1_Gaus->Draw("same");

	// return 1;	


	TH1D *h_Ddls_Err1=new TH1D("h_Ddls_Err1","",100,0,10); h_Ddls_Err1->Sumw2(); 
	TH1D *h_Ddls_Err2=new TH1D("h_Ddls_Err2","",100,0,10); h_Ddls_Err2->Sumw2(); 
	TH1D *h_Ddls_Err3=new TH1D("h_Ddls_Err3","",100,0,10); h_Ddls_Err3->Sumw2(); 
	TH1D *h_Ddls_2to3=new TH1D("h_Ddls_2to3","",100,0,10); h_Ddls_2to3->Sumw2(); 
	TH1D *h_Ddls_1to3=new TH1D("h_Ddls_1to3","",100,0,10); h_Ddls_1to3->Sumw2(); 
	TH1D *h_Ddls_1to3_noSm=new TH1D("h_Ddls_1to3_noSm","",100,0,10); h_Ddls_1to3_noSm->Sumw2(); 

	nt->Project("h_Ddls_Err1","Ddls_FromErr1");
	nt->Project("h_Ddls_Err2","Ddls_FromErr2");
	nt->Project("h_Ddls_Err3","Ddls_FromErr3");
	nt->Project("h_Ddls_2to3","Ddls_2to3");
	nt->Project("h_Ddls_1to3","Ddls_1to3");
	nt->Project("h_Ddls_1to3_noSm","Ddls_1to3_noSm");
	h_Ddls_Err3->SetLineColor(2);
	h_Ddls_Err3->SetMarkerColor(2);
	h_Ddls_2to3->SetLineColor(4);
	h_Ddls_2to3->SetMarkerColor(4);
	h_Ddls_1to3->SetLineColor(kGreen+2);
	h_Ddls_1to3->SetMarkerColor(kGreen+2);

	// f1_Gaus->SetParameter(1,3.5);
	f1_Gaus->SetParameter(1,10);
	TCanvas *c_test=new TCanvas("c_test","c_test",600,600);
	c_test->cd();
	h_Ddls_Err3->Draw();
	h_Ddls_1to3_noSm->Draw("same");

	h_Ddls_Err3->Fit("f1_Gaus","MIN0");
	h_Ddls_Err3->Fit("f1_Gaus","LMIN0");
	h_Ddls_Err3->Fit("f1_Gaus","LMIN0");


	// h_Ddls_2to3->Draw("same");

	// f1_Gaus->SetParameter(1,3.5);

	// h_Ddls_2to3->Fit("f1_Gaus","MIN0");
	// h_Ddls_2to3->Fit("f1_Gaus","LMIN0");
	// h_Ddls_2to3->Fit("f1_Gaus","LMIN0");
	// f1_Gaus->Draw("same");

	h_Ddls_1to3->Draw("same");


	nt->Write();	

	h_Ddls_Err2->Write();
	h_Ddls_Err3->Write();
	h_Ddls_2to3->Write();
	h_Ddls_1to3->Write();


	return 0;

}



