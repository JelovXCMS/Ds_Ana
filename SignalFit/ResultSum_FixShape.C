// Summarize all fit result


#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
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

void ResultSum_FixShape(Int_t isPbPb=0, Int_t doFitFunVar=1, Int_t doCutScanSet=0){

  InitStyle();
  initParameter();

  int nbin_pt=nbin_pt_pp;
  double *bins_pt=bins_pt_pp;

  double *Dchi2clMinScan_bins=Dchi2clMinScan_bins_pp;
  double *DalphaMaxScan_bins=DalphaMaxScan_bins_pp;
  double *DdlsMinScan_bins=DdlsMinScan_bins_pp;

  // TString outName=Form("output/FitResult_pp_pt%.0fto%.0f.root",DptLow,DptHigh);
  TString inName=Form("output%s/FitResult_FixShapeTrkPtScan_pp",s_CutSet.Data());
  TString str_PbPb="pp";

  if(isPbPb==3){
    nbin_pt=nbin_pt_PbPb3;
    bins_pt=bins_pt_PbPb3;

    inName=Form("output%s/FitResult_FixShapeTrkPtScan_PbPb3",s_CutSet.Data());
    str_PbPb="PbPb3";

    Dchi2clMinScan_bins=Dchi2clMinScan_bins_PbPb3;
    DalphaMaxScan_bins=DalphaMaxScan_bins_PbPb3;
    DdlsMinScan_bins=DdlsMinScan_bins_PbPb3;
		
  }

	TFile *fin[nbin_pt];
	TFile *fout=new TFile(Form("output%s/RawFitYield_FixShapeTrkPtScan_%s.root",s_CutSet.Data(),str_PbPb.Data()),"RECREATE");

	TH1F *h_RawBinFitYield_ptbins[nbin_pt];
	TH1F *h_RawRooFitYield_ptbins[nbin_pt];

	TH1F *h_RawBinFitYield = new TH1F("h_RawBinFitYield","h_RawBinFitYield",nbin_pt,bins_pt);
	TH1F *h_RawRooFitYield = new TH1F("h_RawRooFitYield","h_RawRooFitYield",nbin_pt,bins_pt);

	TH1F *h_RawRooFitYield_sigVarP_ptbins[nbin_pt];
	TH1F *h_RawRooFitYield_sigVarN_ptbins[nbin_pt];
	TH1F *h_RawRooFitYield_sigVar_ptbins[nbin_pt];
	TH1F *h_RawRooFitYield_sigVarBkgFix_ptbins[nbin_pt];

	TH1F *h_RawRooFitYield_bkgVar_ptbins[nbin_pt];
	TH1F *h_RawRooFitYield_bkgVarFixSig_ptbins[nbin_pt];
	TH1F *h_RawRooFitYield_bkgVarExp_ptbins[nbin_pt];

	TH1F *h_RawRooFitYield_sigVarP = new TH1F("h_RawRooFitYield_sigVarP","h_RawRooFitYield_sigVarP",nbin_pt,bins_pt);
	TH1F *h_RawRooFitYield_sigVarN = new TH1F("h_RawRooFitYield_sigVarN","h_RawRooFitYield_sigVarN",nbin_pt,bins_pt);
	TH1F *h_RawRooFitYield_sigVar = new TH1F("h_RawRooFitYield_sigVar","h_RawRooFitYield_sigVar",nbin_pt,bins_pt);
	TH1F *h_RawRooFitYield_sigVarBkgFix = new TH1F("h_RawRooFitYield_sigVarBkgFix","h_RawRooFitYield_sigVarBkgFix",nbin_pt,bins_pt);

	TH1F *h_RawRooFitYield_bkgVar = new TH1F("h_RawRooFitYield_bkgVar","h_RawRooFitYield_bkgVar",nbin_pt,bins_pt);
	TH1F *h_RawRooFitYield_bkgVarFixSig = new TH1F("h_RawRooFitYield_bkgVarFixSig","h_RawRooFitYield_bkgVarFixSig",nbin_pt,bins_pt);
	TH1F *h_RawRooFitYield_bkgVarExp = new TH1F("h_RawRooFitYield_bkgVarExp","h_RawRooFitYield_bkgVarExp",nbin_pt,bins_pt);


	TH1F *h_RawRooFitYield_sigVar_RelErr= new TH1F("h_RawRooFitYield_sigVar_RelErr","h_RawRooFitYield_sigVar_RelErr",nbin_pt,bins_pt); h_RawRooFitYield_sigVar_RelErr->Sumw2();
	TH1F *h_RawRooFitYield_sigVarBkgFix_RelErr= new TH1F("h_RawRooFitYield_sigVarBkgFix_RelErr","h_RawRooFitYield_sigVarBkgFix_RelErr",nbin_pt,bins_pt); h_RawRooFitYield_sigVarBkgFix_RelErr->Sumw2();
	TH1F *h_RawRooFitYield_sigVar2_RelErr= new TH1F("h_RawRooFitYield_sigVar2_RelErr","h_RawRooFitYield_sigVar2_RelErr",nbin_pt,bins_pt); h_RawRooFitYield_sigVar2_RelErr->Sumw2();

	TH1F *h_RawRooFitYield_bkgVar_RelErr= new TH1F("h_RawRooFitYield_bkgVar_RelErr","h_RawRooFitYield_bkgVar_RelErr",nbin_pt,bins_pt); h_RawRooFitYield_bkgVar_RelErr->Sumw2();
	TH1F *h_RawRooFitYield_bkgVarFixSig_RelErr= new TH1F("h_RawRooFitYield_bkgVarFixSig_RelErr","h_RawRooFitYield_bkgVarFixSig_RelErr",nbin_pt,bins_pt); h_RawRooFitYield_bkgVarFixSig_RelErr->Sumw2();
	TH1F *h_RawRooFitYield_bkgVarExp_RelErr= new TH1F("h_RawRooFitYield_bkgVarExp_RelErr","h_RawRooFitYield_bkgVarExp_RelErr",nbin_pt,bins_pt); h_RawRooFitYield_bkgVarExp_RelErr->Sumw2();
	TH1F *h_RawRooFitYield_pdfVar_RelErr= new TH1F("h_RawRooFitYield_pdfVar_RelErr","h_RawRooFitYield_pdfVar_RelErr",nbin_pt,bins_pt); h_RawRooFitYield_pdfVar_RelErr->Sumw2();

	// doCutScanSet
	TH1F *h_RawRooFitYield_DalphaMaxScan_ptbins[nbin_pt][nbin_DalphaMaxScan];
	TH1F *h_RawRooFitYield_Dchi2clMinScan_ptbins[nbin_pt][nbin_Dchi2clMinScan];
	TH1F *h_RawRooFitYield_DdlsMinScan_ptbins[nbin_pt][nbin_DdlsMinScan];
	TH1F *h_RawRooFitYield_PhiMassScan_ptbins[nbin_pt][nbin_PhiMassScan];
	TH1F *h_RawRooFitYield_Reschi2clScan_ptbins[nbin_pt][nbin_Reschi2clScan];

	TH1F *h_RawRooFitYield_DalphaMaxScan[nbin_DalphaMaxScan];
	TH1F *h_RawRooFitYield_Dchi2clMinScan[nbin_Dchi2clMinScan];
  TH1F *h_RawRooFitYield_DdlsMinScan[nbin_DdlsMinScan];
  TH1F *h_RawRooFitYield_PhiMassScan[nbin_PhiMassScan];
  TH1F *h_RawRooFitYield_Reschi2clScan[nbin_Reschi2clScan];


	TH1F *h_RawRooFitYield_KaonPtScan_ptbins[nbin_pt][nbin_DauPtScan];
  TH1F *h_RawRooFitYield_KaonPtScan[nbin_DauPtScan];

	TH1F *h_RawRooFitYield_PionPtScan_ptbins[nbin_pt][nbin_DauPtScan];
  TH1F *h_RawRooFitYield_PionPtScan[nbin_DauPtScan];

	TH1F *h_RawRooFitYield_AllDauPtScan_ptbins[nbin_pt][nbin_DauPtScan];
  TH1F *h_RawRooFitYield_AllDauPtScan[nbin_DauPtScan];





	
	for(int i=0; i<nbin_DalphaMaxScan; i++){
		h_RawRooFitYield_DalphaMaxScan[i]=new TH1F(Form("h_RawRooFitYield_DalphaMaxScan_%i",i) , Form("h_RawRooFitYield_DalphaMaxScan_%i",i), nbin_pt,bins_pt );
	}
	for(int i=0; i<nbin_Dchi2clMinScan; i++){
		h_RawRooFitYield_Dchi2clMinScan[i]=new TH1F(Form("h_RawRooFitYield_Dchi2clMinScan_%i",i) , Form("h_RawRooFitYield_Dchi2clMinScan_%i",i), nbin_pt,bins_pt );
	}
	for(int i=0; i<nbin_DdlsMinScan; i++){
		h_RawRooFitYield_DdlsMinScan[i]=new TH1F(Form("h_RawRooFitYield_DdlsMinScan_%i",i) , Form("h_RawRooFitYield_DdlsMinScan_%i",i), nbin_pt,bins_pt );
	}
	for(int i=0; i<nbin_PhiMassScan; i++){
		h_RawRooFitYield_PhiMassScan[i]=new TH1F(Form("h_RawRooFitYield_PhiMassScan_%i",i) , Form("h_RawRooFitYield_PhiMassScan_%i",i), nbin_pt,bins_pt );
	}
	for(int i=0; i<nbin_Reschi2clScan; i++){
		h_RawRooFitYield_Reschi2clScan[i]=new TH1F(Form("h_RawRooFitYield_Reschi2clScan_%i",i) , Form("h_RawRooFitYield_Reschi2clScan_%i",i), nbin_pt,bins_pt );
	}

	for(int i=0; i<nbin_DauPtScan; i++){
		h_RawRooFitYield_KaonPtScan[i]=new TH1F(Form("h_RawRooFitYield_KaonPtScan_%i",i) , Form("h_RawRooFitYield_KaonPtScan_%i",i), nbin_pt,bins_pt );
	}
	for(int i=0; i<nbin_DauPtScan; i++){
		h_RawRooFitYield_PionPtScan[i]=new TH1F(Form("h_RawRooFitYield_PionPtScan_%i",i) , Form("h_RawRooFitYield_PionPtScan_%i",i), nbin_pt,bins_pt );
	}
	for(int i=0; i<nbin_DauPtScan; i++){
		h_RawRooFitYield_AllDauPtScan[i]=new TH1F(Form("h_RawRooFitYield_AllDauPtScan_%i",i) , Form("h_RawRooFitYield_AllDauPtScan_%i",i), nbin_pt,bins_pt );
	}




	Int_t firstbin=0;
	Int_t lastbin=nbin_pt;
	
	double DptLow=0;
	double DptHigh=0;

	if(isPbPb==3){firstbin=2;} // temp, start from 6-8 bin
//	if(isPbPb==3){firstbin=4;} // temp, start from 5-6 bin
	if(isPbPb==0){firstbin=0;} // temp, start from 5-6 bin


	for(Int_t ibin_pt=firstbin; ibin_pt<lastbin; ibin_pt ++){
		DptLow=bins_pt[ibin_pt];
		DptHigh=bins_pt[ibin_pt+1];
		cout<<"\n\nbin : "<<ibin_pt<<" Dpt : "<<DptLow<<" to "<<DptHigh<<endl;


		if(!TFile::Open(Form("%s_pt%.0fto%.0f.root",inName.Data(),DptLow,DptHigh))){
			cout<<Form("file %s_pt%.0fto%.0f.root not exist, skip to next ",inName.Data(),DptLow,DptHigh)<<endl;		
			continue;	
		}

		cout<<"check 1 "<<endl;
	
		fin[ibin_pt]=TFile::Open(Form("%s_pt%.0fto%.0f.root",inName.Data(),DptLow,DptHigh));
		fout->cd();

		cout<<"check 2"<<endl;

		h_RawBinFitYield_ptbins[ibin_pt]=(TH1F*)fin[ibin_pt]->Get("h_RawBinFitYield");
		h_RawBinFitYield->SetBinContent(ibin_pt+1,h_RawBinFitYield_ptbins[ibin_pt]->GetBinContent(1));
		h_RawBinFitYield->SetBinError(ibin_pt+1,h_RawBinFitYield_ptbins[ibin_pt]->GetBinError(1));

		h_RawRooFitYield_ptbins[ibin_pt]=(TH1F*)fin[ibin_pt]->Get("h_RawRooFitYield_FixShape");
		h_RawRooFitYield->SetBinContent(ibin_pt+1,h_RawRooFitYield_ptbins[ibin_pt]->GetBinContent(1));
		h_RawRooFitYield->SetBinError(ibin_pt+1,h_RawRooFitYield_ptbins[ibin_pt]->GetBinError(1));

		double YieldNorm=h_RawRooFitYield_ptbins[ibin_pt]->GetBinContent(1);

		if(doFitFunVar==1){
			
		h_RawRooFitYield_sigVarP_ptbins[ibin_pt]=(TH1F*)fin[ibin_pt]->Get("h_RawRooFitYield_sigVarP");
		h_RawRooFitYield_sigVarP->SetBinContent(ibin_pt+1,h_RawRooFitYield_sigVarP_ptbins[ibin_pt]->GetBinContent(1));
		h_RawRooFitYield_sigVarP->SetBinError(ibin_pt+1,h_RawRooFitYield_sigVarP_ptbins[ibin_pt]->GetBinError(1));

		h_RawRooFitYield_sigVarN_ptbins[ibin_pt]=(TH1F*)fin[ibin_pt]->Get("h_RawRooFitYield_sigVarN");
		h_RawRooFitYield_sigVarN->SetBinContent(ibin_pt+1,h_RawRooFitYield_sigVarN_ptbins[ibin_pt]->GetBinContent(1));
		h_RawRooFitYield_sigVarN->SetBinError(ibin_pt+1,h_RawRooFitYield_sigVarN_ptbins[ibin_pt]->GetBinError(1));

		h_RawRooFitYield_sigVar_ptbins[ibin_pt]=(TH1F*)fin[ibin_pt]->Get("h_RawRooFitYield_sigVar");
		h_RawRooFitYield_sigVar->SetBinContent(ibin_pt+1,h_RawRooFitYield_sigVar_ptbins[ibin_pt]->GetBinContent(1));
		h_RawRooFitYield_sigVar->SetBinError(ibin_pt+1,h_RawRooFitYield_sigVar_ptbins[ibin_pt]->GetBinError(1));

		h_RawRooFitYield_sigVarBkgFix_ptbins[ibin_pt]=(TH1F*)fin[ibin_pt]->Get("h_RawRooFitYield_sigVarBkgFix");
		h_RawRooFitYield_sigVarBkgFix->SetBinContent(ibin_pt+1,h_RawRooFitYield_sigVarBkgFix_ptbins[ibin_pt]->GetBinContent(1));
		h_RawRooFitYield_sigVarBkgFix->SetBinError(ibin_pt+1,h_RawRooFitYield_sigVarBkgFix_ptbins[ibin_pt]->GetBinError(1));


		h_RawRooFitYield_bkgVar_ptbins[ibin_pt]=(TH1F*)fin[ibin_pt]->Get("h_RawRooFitYield_bkgVar");
		h_RawRooFitYield_bkgVar->SetBinContent(ibin_pt+1,h_RawRooFitYield_bkgVar_ptbins[ibin_pt]->GetBinContent(1));
		h_RawRooFitYield_bkgVar->SetBinError(ibin_pt+1,h_RawRooFitYield_bkgVar_ptbins[ibin_pt]->GetBinError(1));

		h_RawRooFitYield_bkgVarFixSig_ptbins[ibin_pt]=(TH1F*)fin[ibin_pt]->Get("h_RawRooFitYield_bkgVarFixSig");
		h_RawRooFitYield_bkgVarFixSig->SetBinContent(ibin_pt+1,h_RawRooFitYield_bkgVarFixSig_ptbins[ibin_pt]->GetBinContent(1));
		h_RawRooFitYield_bkgVarFixSig->SetBinError(ibin_pt+1,h_RawRooFitYield_bkgVarFixSig_ptbins[ibin_pt]->GetBinError(1));


		h_RawRooFitYield_bkgVarExp_ptbins[ibin_pt]=(TH1F*)fin[ibin_pt]->Get("h_RawRooFitYield_bkgVarExp");
		h_RawRooFitYield_bkgVarExp->SetBinContent(ibin_pt+1,h_RawRooFitYield_bkgVarExp_ptbins[ibin_pt]->GetBinContent(1));
		h_RawRooFitYield_bkgVarExp->SetBinError(ibin_pt+1,h_RawRooFitYield_bkgVarExp_ptbins[ibin_pt]->GetBinError(1));

		double Yield_sigVarP=h_RawRooFitYield_sigVarP_ptbins[ibin_pt]->GetBinContent(1);
		double Yield_sigVarN=h_RawRooFitYield_sigVarN_ptbins[ibin_pt]->GetBinContent(1);
		double Yield_sigVar2=h_RawRooFitYield_sigVar_ptbins[ibin_pt]->GetBinContent(1);
		double Yield_sigVarBkgFix=h_RawRooFitYield_sigVarBkgFix_ptbins[ibin_pt]->GetBinContent(1);

		double Yield_bkgVar=h_RawRooFitYield_bkgVar_ptbins[ibin_pt]->GetBinContent(1);
		double Yield_bkgVarFixSig=h_RawRooFitYield_bkgVarFixSig_ptbins[ibin_pt]->GetBinContent(1);
		double Yield_bkgVarExp=h_RawRooFitYield_bkgVarExp_ptbins[ibin_pt]->GetBinContent(1);

		double Yield_sigVar=abs(Yield_sigVarP-YieldNorm) > abs(Yield_sigVarN-YieldNorm) ? abs(Yield_sigVarP-YieldNorm) : abs(Yield_sigVarN-YieldNorm);
		double Yield_sigVar_RelErr=Yield_sigVar/YieldNorm;
		double Yield_sigVarBkgFix_RelErr=abs(Yield_sigVarBkgFix-YieldNorm)/YieldNorm;
		double Yield_sigVar2_RelErr=abs(Yield_sigVar2-YieldNorm)/YieldNorm;

		double Yield_bkgVar_RelErr=abs(Yield_bkgVar-YieldNorm)/YieldNorm;
		double Yield_bkgVarFixSig_RelErr=abs(Yield_bkgVarFixSig-YieldNorm)/YieldNorm;
		double Yield_bkgVarExp_RelErr=abs(Yield_bkgVarExp-YieldNorm)/YieldNorm;
		// double Yield_bkgVarSum_RelErr=Yield_bkgVar_RelErr > Yield_bkgVarExp_RelErr ? Yield_bkgVar_RelErr : Yield_bkgVarExp_RelErr;
		// double Yield_bkgVarSum_RelErr=Yield_bkgVar_RelErr;
		double Yield_bkgVarSum_RelErr=Yield_bkgVarFixSig_RelErr;
		double Yield_pdfVar_RelErr=sqrt(Yield_bkgVarSum_RelErr*Yield_bkgVarSum_RelErr + Yield_sigVarBkgFix_RelErr*Yield_sigVarBkgFix_RelErr );

		h_RawRooFitYield_sigVar_RelErr->SetBinContent(ibin_pt+1, Yield_sigVar_RelErr);
		h_RawRooFitYield_sigVarBkgFix_RelErr->SetBinContent(ibin_pt+1, Yield_sigVarBkgFix_RelErr);
		h_RawRooFitYield_sigVar2_RelErr->SetBinContent(ibin_pt+1, Yield_sigVar2_RelErr);

		h_RawRooFitYield_bkgVar_RelErr->SetBinContent(ibin_pt+1, Yield_bkgVar_RelErr);
		h_RawRooFitYield_bkgVarFixSig_RelErr->SetBinContent(ibin_pt+1, Yield_bkgVarFixSig_RelErr);
		h_RawRooFitYield_bkgVarExp_RelErr->SetBinContent(ibin_pt+1, Yield_bkgVarExp_RelErr);

		h_RawRooFitYield_pdfVar_RelErr->SetBinContent(ibin_pt+1, Yield_pdfVar_RelErr );
		// h_RawRooFitYield_pdfVar_RelErr->SetBinContent(ibin_pt+1, sqrt(Yield_bkgVar_RelErr*Yield_bkgVar_RelErr + Yield_sigVar2_RelErr*Yield_sigVar2_RelErr ));

		cout<<"YieldNorm = "<<YieldNorm<<endl;
//		cout<<"Yield_sigVar_RelErr (Not used) = "<<Yield_sigVar_RelErr*100<<"% , Yield_sigVarP = "<<Yield_sigVarP<<" , Yield_sigVarN = "<<Yield_sigVarN<<endl;
		cout<<"Yield_sigVar2_RelErr (Not use this, single gaussian )= "<<Yield_sigVar2_RelErr*100<<"% , Yield_sigVar2 = "<<Yield_sigVar2<<endl;
		cout<<"Yield_sigVarBkgFix_RelErr (use this, single gaussian )= "<<Yield_sigVarBkgFix_RelErr*100<<"% , Yield_sigVarBkgFix = "<<Yield_sigVarBkgFix<<endl;
		cout<<"Yield_bkgVar_RelErr = "<<Yield_bkgVar_RelErr*100<<"% , Yield_bkgVar = "<<Yield_bkgVar<<endl;
		cout<<"Yield_bkgVarFixSig_RelErr = "<<Yield_bkgVarFixSig_RelErr*100<<"% , Yield_bkgVarFixSig = "<<Yield_bkgVarFixSig<<endl;
//		cout<<"Yield_bkgVarExp_RelErr = "<<Yield_bkgVarExp_RelErr*100<<"% , Yield_bkgVarExp = "<<Yield_bkgVarExp<<endl;
//		cout<<"Yield_bkgVarSum_RelErr = "<<Yield_bkgVarSum_RelErr*100<<"% "<<endl;
		cout<<"Yield_pdfVar_RelErr Sum = "<<Yield_pdfVar_RelErr*100<<"% "<<endl;


		// consider directly write a latex format output here

		} // end doFitFunVar

		cout<<"check 3 "<<endl;

		//// -- doCutScanSet --////
		if(doCutScanSet==1){
			cout<<"doCutScan"<<endl;
		// DalphaMaxScan 
		for(int ibin=0; ibin<nbin_DalphaMaxScan;ibin++){
			h_RawRooFitYield_DalphaMaxScan_ptbins[ibin_pt][ibin]=(TH1F*)fin[ibin_pt]->Get(Form("h_RawRooFitYield_DalphaMaxScan%.0f",DalphaMaxScan_bins[ibin]*100));
			h_RawRooFitYield_DalphaMaxScan[ibin]->SetBinContent(ibin_pt+1, h_RawRooFitYield_DalphaMaxScan_ptbins[ibin_pt][ibin]->GetBinContent(1));
			h_RawRooFitYield_DalphaMaxScan[ibin]->SetBinError(ibin_pt+1, h_RawRooFitYield_DalphaMaxScan_ptbins[ibin_pt][ibin]->GetBinError(1));
		}

		// DdlsMinScan 
		for(int ibin=0; ibin<nbin_DdlsMinScan;ibin++){
			h_RawRooFitYield_DdlsMinScan_ptbins[ibin_pt][ibin]=(TH1F*)fin[ibin_pt]->Get(Form("h_RawRooFitYield_DdlsMinScan%.0f",DdlsMinScan_bins[ibin]*100));
			h_RawRooFitYield_DdlsMinScan[ibin]->SetBinContent(ibin_pt+1, h_RawRooFitYield_DdlsMinScan_ptbins[ibin_pt][ibin]->GetBinContent(1));
			h_RawRooFitYield_DdlsMinScan[ibin]->SetBinError(ibin_pt+1, h_RawRooFitYield_DdlsMinScan_ptbins[ibin_pt][ibin]->GetBinError(1));
		}

		// Dchi2clMinScan 
		for(int ibin=0; ibin<nbin_Dchi2clMinScan;ibin++){
			h_RawRooFitYield_Dchi2clMinScan_ptbins[ibin_pt][ibin]=(TH1F*)fin[ibin_pt]->Get(Form("h_RawRooFitYield_Dchi2clMinScan%.0f",Dchi2clMinScan_bins[ibin]*100));
			h_RawRooFitYield_Dchi2clMinScan[ibin]->SetBinContent(ibin_pt+1, h_RawRooFitYield_Dchi2clMinScan_ptbins[ibin_pt][ibin]->GetBinContent(1));
			h_RawRooFitYield_Dchi2clMinScan[ibin]->SetBinError(ibin_pt+1, h_RawRooFitYield_Dchi2clMinScan_ptbins[ibin_pt][ibin]->GetBinError(1));
		}
		// PhiMassScan 
		for(int ibin=0; ibin<nbin_PhiMassScan;ibin++){
			h_RawRooFitYield_PhiMassScan_ptbins[ibin_pt][ibin]=(TH1F*)fin[ibin_pt]->Get(Form("h_RawRooFitYield_PhiMassScan_%i",ibin));
			h_RawRooFitYield_PhiMassScan[ibin]->SetBinContent(ibin_pt+1, h_RawRooFitYield_PhiMassScan_ptbins[ibin_pt][ibin]->GetBinContent(1));
			h_RawRooFitYield_PhiMassScan[ibin]->SetBinError(ibin_pt+1, h_RawRooFitYield_PhiMassScan_ptbins[ibin_pt][ibin]->GetBinError(1));
		}

		cout<<"finished PhimassScan"<<endl;
/*
		// KaonPtScan 
		for(int ibin=0; ibin<nbin_DauPtScan;ibin++){
			h_RawRooFitYield_KaonPtScan_ptbins[ibin_pt][ibin]=(TH1F*)fin[ibin_pt]->Get(Form("h_RawRooFitYield_KaonPtScan_%i",ibin));
			h_RawRooFitYield_KaonPtScan[ibin]->SetBinContent(ibin_pt+1, h_RawRooFitYield_KaonPtScan_ptbins[ibin_pt][ibin]->GetBinContent(1));
			h_RawRooFitYield_KaonPtScan[ibin]->SetBinError(ibin_pt+1, h_RawRooFitYield_KaonPtScan_ptbins[ibin_pt][ibin]->GetBinError(1));
		}

		cout<<"finished KaonScan"<<endl;
		// PionPtScan 
		for(int ibin=0; ibin<nbin_DauPtScan;ibin++){
			h_RawRooFitYield_PionPtScan_ptbins[ibin_pt][ibin]=(TH1F*)fin[ibin_pt]->Get(Form("h_RawRooFitYield_PionPtScan_%i",ibin));
			h_RawRooFitYield_PionPtScan[ibin]->SetBinContent(ibin_pt+1, h_RawRooFitYield_PionPtScan_ptbins[ibin_pt][ibin]->GetBinContent(1));
			h_RawRooFitYield_PionPtScan[ibin]->SetBinError(ibin_pt+1, h_RawRooFitYield_PionPtScan_ptbins[ibin_pt][ibin]->GetBinError(1));
		}

		cout<<"finished PionScan"<<endl;

		// AllDauPtScan 
		for(int ibin=0; ibin<nbin_DauPtScan;ibin++){
			h_RawRooFitYield_AllDauPtScan_ptbins[ibin_pt][ibin]=(TH1F*)fin[ibin_pt]->Get(Form("h_RawRooFitYield_AllDauPtScan_%i",ibin));
			h_RawRooFitYield_AllDauPtScan[ibin]->SetBinContent(ibin_pt+1, h_RawRooFitYield_AllDauPtScan_ptbins[ibin_pt][ibin]->GetBinContent(1));
			h_RawRooFitYield_AllDauPtScan[ibin]->SetBinError(ibin_pt+1, h_RawRooFitYield_AllDauPtScan_ptbins[ibin_pt][ibin]->GetBinError(1));
		}

		cout<<"finished AllDauPtScan"<<endl;



		// Reschi2clScan 
		for(int ibin=0; ibin<nbin_Reschi2clScan;ibin++){
			h_RawRooFitYield_Reschi2clScan_ptbins[ibin_pt][ibin]=(TH1F*)fin[ibin_pt]->Get(Form("h_RawRooFitYield_Reschi2clScan_%i",ibin));
			h_RawRooFitYield_Reschi2clScan[ibin]->SetBinContent(ibin_pt+1, h_RawRooFitYield_Reschi2clScan_ptbins[ibin_pt][ibin]->GetBinContent(1));
			h_RawRooFitYield_Reschi2clScan[ibin]->SetBinError(ibin_pt+1, h_RawRooFitYield_Reschi2clScan_ptbins[ibin_pt][ibin]->GetBinError(1));
		}
*/


	
		} //// -- end doCutScanSet


	}

	fout->cd();
	h_RawBinFitYield->Write("",TObject::kOverwrite);
	h_RawRooFitYield->Write("",TObject::kOverwrite);

	if(doFitFunVar){
	h_RawRooFitYield_sigVarP->Write("",TObject::kOverwrite);
	h_RawRooFitYield_sigVarN->Write("",TObject::kOverwrite);
	h_RawRooFitYield_sigVar->Write("",TObject::kOverwrite);
	h_RawRooFitYield_sigVarBkgFix->Write("",TObject::kOverwrite);

	h_RawRooFitYield_bkgVar->Write("",TObject::kOverwrite);
	h_RawRooFitYield_bkgVarFixSig->Write("",TObject::kOverwrite);
	h_RawRooFitYield_bkgVarExp->Write("",TObject::kOverwrite);

  h_RawRooFitYield_sigVar_RelErr->Write("",TObject::kOverwrite);
  h_RawRooFitYield_sigVarBkgFix_RelErr->Write("",TObject::kOverwrite);
  h_RawRooFitYield_sigVar2_RelErr->Write("",TObject::kOverwrite);

  h_RawRooFitYield_bkgVarExp_RelErr->Write("",TObject::kOverwrite);
  h_RawRooFitYield_bkgVarFixSig_RelErr->Write("",TObject::kOverwrite);
  h_RawRooFitYield_bkgVar_RelErr->Write("",TObject::kOverwrite);

  h_RawRooFitYield_pdfVar_RelErr->Write("",TObject::kOverwrite);
	}

	if(doCutScanSet){
		for(int ibin=0; ibin<nbin_DalphaMaxScan;ibin++){
			h_RawRooFitYield_DalphaMaxScan[ibin]->Write("",TObject::kOverwrite);
		}

		for(int ibin=0; ibin<nbin_Dchi2clMinScan;ibin++){
			h_RawRooFitYield_Dchi2clMinScan[ibin]->Write("",TObject::kOverwrite);
		}
		for(int ibin=0; ibin<nbin_DdlsMinScan;ibin++){
			h_RawRooFitYield_DdlsMinScan[ibin]->Write("",TObject::kOverwrite);
		}
		for(int ibin=0; ibin<nbin_PhiMassScan;ibin++){
			h_RawRooFitYield_PhiMassScan[ibin]->Write("",TObject::kOverwrite);
		}
/*
		for(int ibin=0; ibin<nbin_Reschi2clScan;ibin++){
			h_RawRooFitYield_Reschi2clScan[ibin]->Write("",TObject::kOverwrite);
		}

		for(int ibin=0; ibin<nbin_DauPtScan;ibin++){
			h_RawRooFitYield_KaonPtScan[ibin]->Write("",TObject::kOverwrite);
		}
		for(int ibin=0; ibin<nbin_DauPtScan;ibin++){
			h_RawRooFitYield_PionPtScan[ibin]->Write("",TObject::kOverwrite);
		}
		for(int ibin=0; ibin<nbin_DauPtScan;ibin++){
			h_RawRooFitYield_AllDauPtScan[ibin]->Write("",TObject::kOverwrite);
		}
*/
	}




}
