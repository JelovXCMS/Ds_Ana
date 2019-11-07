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
#include "RooGoF.C"

#include "TFitter.h"
#include "TFitResult.h"
#include "TRandom3.h"
using namespace RooFit;
using namespace std;

// double textposx=0.2;
// double textposy=0.77;

double shiftY=0.0;
double shiftX=0.35;
double oneshift=0.075;

  double DsMassMCMean=1.96834;
  double DsMass3SigWidth=0.03;


int smearTest(int isPbPb=0,int PNP=1,int nSmear=10){
	
	TRandom3 Rdm;
	Rdm.SetSeed(0);

	TString s_PbPb3="pp";
	if(isPbPb){
		s_PbPb3="PbPb3";
	}	
	TString s_PNP="Prompt";
	if(PNP){
		s_PNP="NonPrompt";
	}


  TString mcName=Form("/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_cutScan/%sMC_phi%s_fitFile.root",s_PbPb3.Data(),s_PNP.Data());
	TFile *f_mcp=new TFile(mcName.Data(),"read");
	TString mcName_New=Form("/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_cutScan/%sMC_phi%s_fitFile_Smear.root",s_PbPb3.Data(),s_PNP.Data());

	TFile *f_mcp_new=new TFile(mcName_New.Data(),"recreate");

	TTree *t_mcp=(TTree*)f_mcp->Get("t_fit");
	// t_mcp->Draw("Dpt");

  Float_t Dmass;
  Float_t Dpt;
  Float_t Ddls;
  Float_t Ddl;
  Float_t DdlErr;
  Float_t Dalpha;
  Float_t Dchi2cl;
	Float_t DtktkResmass;
	Float_t TotalWeight;	
  Float_t DlxyBS;
  Float_t DlxyBSErr;
  Float_t DlxyBSs;
	Float_t Dalpha_BS_2D;

	t_mcp->SetBranchAddress("Dmass",&Dmass);
	t_mcp->SetBranchAddress("Ddls",&Ddls);
	t_mcp->SetBranchAddress("Dpt",&Dpt);
	t_mcp->SetBranchAddress("DdlErr",&DdlErr);
	t_mcp->SetBranchAddress("Ddl",&Ddl);
	t_mcp->SetBranchAddress("DtktkResmass",&DtktkResmass);
	t_mcp->SetBranchAddress("Dalpha",&Dalpha);
	t_mcp->SetBranchAddress("Dchi2cl",&Dchi2cl);
	t_mcp->SetBranchAddress("TotalWeight",&TotalWeight);

	t_mcp->SetBranchAddress("Dalpha_BS_2D",&Dalpha_BS_2D);
	t_mcp->SetBranchAddress("DlxyBS",&DlxyBS);
	t_mcp->SetBranchAddress("DlxyBSErr",&DlxyBSErr);
	t_mcp->SetBranchAddress("DlxyBSs",&DlxyBSs);


	// Float_t Ddls_new;
	// Float_t Ddls_1p;
	// Float_t Ddls_5p;
	// Float_t Ddls_10p;
	// Float_t Ddls_20p;
	Float_t SMrWt=(float)1.0/nSmear;

	Float_t SmearRelFactorArr[]={0.01,0.05,0.1,0.2,0.3};
	Float_t SmearAbsFactorArr[]={0.0001,0.001,0.02,0.05,0.1};
	// Float_t SmearFactorArr[]={0.01,0.05};
	const	int nSmrRelF=sizeof(SmearRelFactorArr)/sizeof(SmearRelFactorArr[0]);
	const	int nSmrAbsF=sizeof(SmearAbsFactorArr)/sizeof(SmearAbsFactorArr[0]);
	Float_t Ddls_SmrRelF[nSmrRelF];
	Float_t Ddls_SmrAbsF[nSmrAbsF];

	f_mcp_new->cd();
	TTree *t_new=new TTree("t_fit_smear","t_fit_smear");
	t_new->Branch("Dmass",&Dmass);
	t_new->Branch("Ddls",&Ddls);
	t_new->Branch("Ddl",&Ddl);
	t_new->Branch("DdlErr",&DdlErr);
	t_new->Branch("Dalpha",&Dalpha);
	t_new->Branch("Dchi2cl",&Dchi2cl);
	t_new->Branch("DtktkResmass",&DtktkResmass);
	t_new->Branch("TotalWeight",&TotalWeight);
	// t_new->Branch("Ddls_new",&Ddls_new);
	// t_new->Branch("Ddls_1p",&Ddls_1p);
	// t_new->Branch("Ddls_5p",&Ddls_5p);
	// t_new->Branch("Ddls_10p",&Ddls_10p);
	// t_new->Branch("Ddls_20p",&Ddls_20p);
	t_new->Branch("SMrWt",&SMrWt);
	for(int i=0; i<nSmrRelF; i++){
		t_new->Branch(Form("Ddls_Rel%.0fem5",SmearRelFactorArr[i]*1e5),&Ddls_SmrRelF[i]);
	}
	for(int i=0; i<nSmrAbsF; i++){
		t_new->Branch(Form("Ddls_Abs%.0fem5",SmearAbsFactorArr[i]*1e5),&Ddls_SmrAbsF[i]);
	}


	Long64_t nentries=t_mcp->GetEntries();
	double smF=0.1;

		

	for(Long64_t i=0;i<nentries; i++)
	// for(Long64_t i=0;i<500; i++)
	{
		
		if(i%200000==0) {cout<<setw(10)<<i<<" / "<<nentries<<endl;}
		t_mcp->GetEntry(i);
		for(int j=0; j<nSmear;j++){


			for(int k=0;k<nSmrRelF;k++){
		   Ddls_SmrRelF[k]=Ddl / Rdm.Gaus(DdlErr,DdlErr*SmearRelFactorArr[k]);
			 while(Ddls_SmrRelF[k]<0){
		       Ddls_SmrRelF[k]=Ddl / Rdm.Gaus(DdlErr,DdlErr*SmearRelFactorArr[k]);
			 }			 
			}
			for(int k=0;k<nSmrAbsF;k++){
		   Ddls_SmrAbsF[k]=Ddl / Rdm.Gaus(DdlErr,DdlErr*SmearAbsFactorArr[k]);
			 while(Ddls_SmrAbsF[k]<0){
		       Ddls_SmrAbsF[k]=Ddl / Rdm.Gaus(DdlErr,DdlErr*SmearAbsFactorArr[k]);
			 }			 
			}

		// Ddls_new=Rdm.Gaus(Ddls,Ddls*smF);
		// Ddls_1p=Rdm.Gaus(Ddls,Ddls*0.01);
		// Ddls_5p=Rdm.Gaus(Ddls,Ddls*0.05);
		// Ddls_10p=Rdm.Gaus(Ddls,Ddls*0.10);
		// Ddls_20p=Rdm.Gaus(Ddls,Ddls*0.20);
		
		// cout<<"Ddls = "<<Ddls<<" , Ddls_new = "<<Ddls_new<<endl;
		
		t_new->Fill();
		}
		
	}

	f_mcp_new->cd();
	t_new->Write("",TObject::kOverwrite);

/*
	TH1D *h_Ddls_ori=new TH1D("h_Ddls_ori","h_Ddls_ori",200,0,20); h_Ddls_ori->Sumw2();
	TH1D *h_Ddls_1p=new TH1D("h_Ddls_1p","h_Ddls_1p",200,0,20); h_Ddls_1p->Sumw2();
	TH1D *h_Ddls_10p=new TH1D("h_Ddls_10p","h_Ddls_10p",200,0,20); h_Ddls_10p->Sumw2();

	t_new->Project("h_Ddls_ori","Ddls","Ddls>2.5");
	t_new->Project("h_Ddls_1p","Ddls_1p","Ddls_1p>2.5");
	t_new->Project("h_Ddls_10p","Ddls_10p","Ddls_10p>2.5");


	h_Ddls_ori->Draw();
	h_Ddls_1p->SetLineColor(2);
	h_Ddls_1p->Draw("same");
	h_Ddls_10p->SetLineColor(4);
	h_Ddls_10p->Draw("same");

	cout<<"h_Ddls_ori "<<h_Ddls_ori->GetMean()<<endl;
	cout<<"h_Ddls_1p "<<h_Ddls_1p->GetMean()<<endl;
	cout<<"h_Ddls_10p "<<h_Ddls_10p->GetMean()<<endl;
*/

	cout<<"fout : "<<mcName_New<<endl;

	return 1;
}
