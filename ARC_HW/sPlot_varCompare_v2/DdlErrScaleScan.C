#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"

#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting_simple.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"

#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"

#include <TObjArray.h>
#include <TFractionFitter.h>
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
#include "RooWorkspace.h"
#include "RooConstVar.h"

#include "varCompare_para.h"

using namespace RooStats;


using namespace RooFit;
using namespace std;
double DsDataFitRangeLow =1.91;
double DsDataFitRangeHigh = 2.11;
Int_t  nbin_DmassDraw=50;

// double textposx=0.2;
// double textposy=0.77;

double shiftY=0.0;
double shiftX=0.45;
double oneshift=0.075;

TCanvas *c_binfit_MC[500];
TCanvas *c_binfit_Data[500];
TCanvas *c_roofit_MC[500] ;
TCanvas *c_roofit_Data[500];
TCanvas *c_roofit_Data_pull[500];
TCanvas *c_roofit_Data_withpull[500];

TCanvas *c_roofit_MCstudy[500];
TCanvas *c_roofit_GoF[500];

TCanvas *c_roofit_Data_Dca[500][30];

int count_c_binfit=0;
int count_c_rootfit=0;

double FloatWidthErr_Norm=0;
double FloatWidthVal_Norm=0;
double FloatWidth_DsMass_Norm=1.965;

double DmassSideBand1Low=1.91;
double DmassSideBand1High=1.93;
double DmassSideBand2Low=2.01;
double DmassSideBand2High=2.03;
double Dmass2SigLow=1.949;
double Dmass2SigHigh=1.989;

Float_t SmearRelFactorArr[]={0.01,0.05,0.1,0.2,0.3};
Float_t SmearAbsFactorArr[]={0.0001,0.001,0.005,0.02,0.03,0.05,0.08,0.1,0.15};
// Float_t ScaleErrFactorArr[]={0.7,0.8,0.85,0.9,0.95,0.98,1,1.02,1.05,1.1,1.15,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9};
// Float_t SmearFactorArr[]={0.01,0.05};
const int nSmrRelF=sizeof(SmearRelFactorArr)/sizeof(SmearRelFactorArr[0]);
const int nSmrAbsF=sizeof(SmearAbsFactorArr)/sizeof(SmearAbsFactorArr[0]);
// const int nSclErrF=sizeof(ScaleErrFactorArr)/sizeof(ScaleErrFactorArr[0]);
Float_t Ddls_SmrRelF[nSmrRelF];
Float_t DdlErr_SmrRelF[nSmrRelF];
Float_t Ddls_SmrAbsF[nSmrAbsF];
Float_t DdlErr_SmrAbsF[nSmrAbsF];
// Float_t Ddls_SclErrF[nSclErrF];
// Float_t Ddl_SclErrF[nSclErrF];
// Float_t DdlErr_SclErrF[nSclErrF];

TH1D *hDsDlsDataBkg;
TH1D *hDsDlsMCPSignal;
TH1D *hDsDlsMCNPSignal;

void normalize(TH1D* h)
{
	h->Sumw2();
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		Float_t val=h->GetBinContent(i);
		Float_t valErr=h->GetBinError(i);
		h->SetBinContent(i,val/h->GetBinWidth(i));
		h->SetBinError(i,valErr/h->GetBinWidth(i));
	}
	h->Scale(1./h->Integral(0,100,"width"));
}

void fixNegBin(TH1D* h){
	for(int i=1;i<=h->GetNbinsX();i++){
		if(h->GetBinContent(i)<0 ){
			h->SetBinContent(i,0);
			h->SetBinError(i,0);
		}
	}
}

Double_t funMix(Double_t* x_, Double_t* para)
{
	float x = x_[0];
	float APrompt = para[0];
	// float ANonPrompt = para[1];
	float ANonPrompt = 1-APrompt;
	float promptYield = 0;
	float nonPromptYield = 0;

	promptYield = hDsDlsMCPSignal->GetBinContent(hDsDlsMCPSignal->GetXaxis()->FindBin(x));
	nonPromptYield = hDsDlsMCNPSignal->GetBinContent(hDsDlsMCNPSignal->GetXaxis()->FindBin(x));

	return APrompt*promptYield+ANonPrompt*nonPromptYield;
}

Double_t funMixWBkg(Double_t* x_, Double_t* para)
{
	float x = x_[0];
	float ABkg = para[0];
	float APrompt = para[1];
	// float ANonPrompt = para[1];
	float ANonPrompt = 1-APrompt;
	float promptYield = 0;
	float nonPromptYield = 0;
	float BkgYield =0;

	BkgYield =hDsDlsDataBkg->GetBinContent(hDsDlsDataBkg->GetXaxis()->FindBin(x));
	promptYield = hDsDlsMCPSignal->GetBinContent(hDsDlsMCPSignal->GetXaxis()->FindBin(x));
	nonPromptYield = hDsDlsMCNPSignal->GetBinContent(hDsDlsMCNPSignal->GetXaxis()->FindBin(x));

	return ABkg*BkgYield+ (1-ABkg)*(APrompt*promptYield+ANonPrompt*nonPromptYield);
}


Double_t funSigDdls(Double_t* x_, Double_t* para)
{
	float x = x_[0];
	float ABkg = para[0];
	float APrompt = para[1];
	// float ANonPrompt = para[1];
	float ANonPrompt = 1-APrompt;
	float promptYield = 0;
	float nonPromptYield = 0;
	float BkgYield =0;

	// BkgYield =hDsDlsDataBkg->GetBinContent(hDsDlsDataBkg->GetXaxis()->FindBin(x));
	promptYield = hDsDlsMCPSignal->GetBinContent(hDsDlsMCPSignal->GetXaxis()->FindBin(x));
	nonPromptYield = hDsDlsMCNPSignal->GetBinContent(hDsDlsMCNPSignal->GetXaxis()->FindBin(x));

	return (1-ABkg)*(APrompt*promptYield+ANonPrompt*nonPromptYield);
}





Double_t funNonPrompt(Double_t* x_, Double_t* para)
{
	float x = x_[0];
	float APrompt = para[0];
	float ANonPrompt = 1-APrompt;
	float nonPromptYield = 0;
	nonPromptYield = hDsDlsMCNPSignal->GetBinContent(hDsDlsMCNPSignal->GetXaxis()->FindBin(x));
	return ANonPrompt*nonPromptYield;
}


int smearTree(TFile *,TFile *,int,Float_t*, int , double ,double);

TH1D *ProjectData(TString var_compare, TString DataCuts, TTree *t_Ds_Data, int nbin_var, double bins_var_Low, double bins_var_High){

	// TH1D *htemp= 
	TH1D *h_temp= new TH1D(Form("h_%s",var_compare.Data()),Form(";%s",var_compare.Data()) , nbin_var, bins_var_Low,bins_var_High); h_temp->Sumw2();
	t_Ds_Data->Project(Form("h_%s", var_compare.Data()) , var_compare.Data() , Form("%s", DataCuts.Data()) );
	h_temp->Scale(1/h_temp->Integral());
	return h_temp;

}

int GetDisDat(double sigFrac, TString var_compare, TString DataCuts, TH1D *h_var_Data2sig, TH1D *h_var_sideband, TH1D *h_var_Data2sigMix, TTree *t_Ds_Data, int nbin_var, double bins_var_Low, double bins_var_High, TFile *ftemp){

	ftemp->cd();

	h_var_sideband= new TH1D(Form("h_%s_sideband",var_compare.Data()),Form("h_%s_sideband",var_compare.Data()),nbin_var, bins_var_Low,bins_var_High); h_var_sideband->Sumw2();
	t_Ds_Data->Project(Form("h_%s_sideband",var_compare.Data()),var_compare.Data(),Form("( (Dmass>%f && Dmass <%f ) || ( Dmass >%f && Dmass<%f) ) && %s",DmassSideBand1Low, DmassSideBand1High, DmassSideBand2Low, DmassSideBand2High, DataCuts.Data()));
	h_var_sideband->Scale(1/h_var_sideband->Integral());

	// h_var_sideband->Draw();

	// h_var_Data2sig= new TH1D(Form("h_%s_Data2sig",var_compare.Data()),Form("h_%s_Data2sig",var_compare.Data()),nbin_var, bins_var_Low,bins_var_High); h_var_Data2sig->Sumw2();
	t_Ds_Data->Project(Form("h_%s_Data2sig",var_compare.Data()),var_compare.Data(),Form("(Dmass>%f && Dmass<%f) && %s",Dmass2SigLow,Dmass2SigHigh, DataCuts.Data()));
	h_var_Data2sig->Scale(1/h_var_Data2sig->Integral());


	h_var_Data2sigMix=(TH1D*)h_var_Data2sig->Clone(Form("h_%s_Data2sigMix",var_compare.Data()));

	h_var_Data2sig->Add(h_var_sideband,-(1-sigFrac));
	h_var_Data2sig->Scale(1/h_var_Data2sig->Integral());
	h_var_Data2sig->SetTitle("");
	// h_var_Data2sig->Rebin(nRebin);
	h_var_Data2sig->GetXaxis()->SetRangeUser(0,30);
	h_var_Data2sig->Draw();

	h_var_Data2sig->Write();


	return 0;

}


int DdlErrScaleScan(int isPbPb=0, TString var_compare="DdlErr", TString var_cut="Dpt" , double var_cutLow=10, double var_cutHigh=20, double bins_var_Low=0, double bins_var_High=1, double bins_var_DrawLow=0, double bins_var_DrawHigh=0.1 , int nRebin=1, TString MC_DefaultWeight="TotalWeight",int doSmear=0, int nSmear=10, int doCreateSMTree=1){

	// setTDRStyle();
	// InitStyle();


	TString GenWt="";
	if(MC_DefaultWeight=="D0Weight"){
		GenWt="_D0Weight";
	}

	gSystem->Exec("mkdir -p fitout");
	gSystem->Exec("mkdir -p plots/fit");
	gSystem->Exec(Form("mkdir -p plots/%s",var_compare.Data()));

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1111);
	TString Str_PbPb="pp";
	if(isPbPb){Str_PbPb="PbPb";}

	TLatex *tltx=new TLatex();
	TString dataName="./rootF/pp_fitFile.root";
	TString mcName_Prompt="./rootF/ppMC_phiPrompt_fitFile.root";
	TString mcName_NonPrompt="./rootF/ppMC_phiNonPrompt_fitFile.root";
	TString mcName_NonPrompt_f0="./rootF/ppMC_f0NonPrompt_fitFile.root";

	double Dpt_Low=Dpt_Low_pp;
	double Dpt_High=Dpt_Hight_pp;

	double Dalpha_cut=Dalpha_cut_pp;
	double Dchi2cl_cut=Dchi2cl_cut_pp;
	double Ddls_cut=Ddls_cut_pp;

	double *bins_var=bins_DdxyzErr_pp;
	int nbin_var=nbin_DdxyzErr_pp;
	nbin_var=1000;
	nbin_var=nbin_var/nRebin;


	if(var_compare=="Dchi2cl"){
		Dchi2cl_cut=0.02;
	}
	TString s_PbPb="pp";

	if(isPbPb){
		Dpt_Low=Dpt_Low_PbPb;
		Dpt_High=Dpt_Hight_PbPb;

		dataName="./rootF/PbPb3_fitFile.root";
		mcName_Prompt="./rootF/PbPb3MC_phiPrompt_fitFile.root";
		mcName_NonPrompt="./rootF/PbPb3MC_phiNonPrompt_fitFile.root";
		mcName_NonPrompt_f0="./rootF/PbPb3MC_f0NonPrompt_fitFile.root";

		Dalpha_cut=Dalpha_cut_PbPb3;
		Dchi2cl_cut=Dchi2cl_cut_PbPb3;
		Ddls_cut=Ddls_cut_PbPb3;
		s_PbPb="PbPb";

		if(var_compare=="Dchi2cl"){
			Dchi2cl_cut=0.05;
		}

	}

	if(var_cut=="Dpt"){
		Dpt_Low=var_cutLow;
		Dpt_High=var_cutHigh;
	}


		int rebinN=4;
		if(var_compare=="Dchi2cl"){
			rebinN=20;
		} 
		if(Dpt_Low>=20){
			rebinN=8;
		} 



	TString DataCuts=Form("Dpt>%f && Dpt<%f && Dalpha < %f && Dchi2cl > %f && Ddls >%f && %s > %f && %s <%f", Dpt_Low, Dpt_High, Dalpha_cut, Dchi2cl_cut, Ddls_cut, var_cut.Data(), var_cutLow,var_cut.Data(), var_cutHigh);

	TString DataCuts_2sig=Form("(Dmass>%f && Dmass<%f) && Dpt>%f && Dpt<%f && Dalpha < %f && Dchi2cl > %f && Ddls >%f && %s > %f && %s <%f",Dmass2SigLow,Dmass2SigHigh , Dpt_Low, Dpt_High, Dalpha_cut, Dchi2cl_cut, Ddls_cut, var_cut.Data(), var_cutLow,var_cut.Data(), var_cutHigh);

	TString DataCuts_sideband=Form("( (Dmass>%f && Dmass <%f ) || ( Dmass >%f && Dmass<%f) ) && Dpt>%f && Dpt<%f && Dalpha < %f && Dchi2cl > %f && Ddls >%f && %s > %f && %s <%f",DmassSideBand1Low, DmassSideBand1High, DmassSideBand2Low, DmassSideBand2High , Dpt_Low, Dpt_High, Dalpha_cut, Dchi2cl_cut, Ddls_cut, var_cut.Data(), var_cutLow,var_cut.Data(), var_cutHigh);

	TFile *f_data=TFile::Open(dataName.Data());
	TFile *f_mc_Prompt=TFile::Open(mcName_Prompt.Data());
	TFile *f_mc_NonPrompt=TFile::Open(mcName_NonPrompt.Data());
	TTree *t_Ds_MCPrompt=(TTree*)f_mc_Prompt->Get(Form("t_fit"));
	TTree *t_Ds_MCNonPrompt=(TTree*)f_mc_NonPrompt->Get(Form("t_fit"));
	TTree *t_Ds_Data=(TTree*)f_data->Get(Form("t_fit"));



	TFile *fin=TFile::Open(Form("./fitout/%s_%s%.0fto%.0f.root",Str_PbPb.Data(), var_cut.Data(), var_cutLow*100, var_cutHigh*100) ,"read");
	gSystem->Exec("mkdir -p ./DdlErrScaleOut");
	TFile *fout=TFile::Open(Form("./DdlErrScaleOut/%s_%s%.0fto%.0f_DdlErrScale.root",Str_PbPb.Data(), var_cut.Data(), var_cutLow, var_cutHigh) ,"recreate");
	TH1D *h_DataRawYield=(TH1D*)fin->Get("h_DataRawYield");
	TH1D *h_Data2sig_SigFrac=(TH1D*)fin->Get("h_Data2sig_SigFrac");	

	TFile *ftemp=new TFile("ftemp.root","recreate");

	double sigFrac=h_Data2sig_SigFrac->GetBinContent(1);

	var_compare="Ddl";
	bins_var_Low=0;
	bins_var_High=10;
	bins_var_DrawLow=0;
	bins_var_DrawHigh=1;

/*
	 var_compare="DdlErr";
	 bins_var_Low=0;
	 bins_var_High=1;
	 bins_var_DrawLow=0;
	 bins_var_DrawHigh=0.1;

	if(Dpt_High<=10){
		bins_var_DrawHigh=0.05;
	}
*/
	// TH1D *h_var_sideband;
	// TH1D *h_var_Data2sig;
	// TH1D *h_var_Data2sigMix;
	// GetDisDat(sigFrac, var_compare, DataCuts, h_var_Data2sig, h_var_sideband, h_var_Data2sigMix, t_Ds_Data, nbin_var, bins_var_Low, bins_var_High,ftemp); // somehow not work...

	TH1D *h_DdlErr_sideband=ProjectData(var_compare, DataCuts_sideband, t_Ds_Data, nbin_var, bins_var_Low, bins_var_High);
	TH1D *h_DdlErr_Data2sigMix=ProjectData(var_compare, DataCuts_2sig, t_Ds_Data, nbin_var, bins_var_Low, bins_var_High);

	TH1D *h_DdlErr_DataSig=(TH1D*)h_DdlErr_Data2sigMix->Clone("h_DdlErr_DataSig");
	// h_DdlErr_DataSig->Scale(1/h_DdlErr_DataSig->Integral());
	h_DdlErr_DataSig->Add(h_DdlErr_sideband,-(1-sigFrac));
	fixNegBin(h_DdlErr_DataSig);
	h_DdlErr_DataSig->Scale(1/h_DdlErr_DataSig->Integral());

	/*
		 h_DdlErr_Data2sigMix->Rebin(8);
		 h_DdlErr_Data2sigMix->GetXaxis()->SetRangeUser(bins_var_DrawLow,bins_var_DrawHigh);
		 h_DdlErr_Data2sigMix->Draw();

		 h_DdlErr_sideband->SetLineColor(2);
		 h_DdlErr_sideband->Rebin(8);
	// h_DdlErr_sideband->GetXaxis()->SetRangeUser(bins_var_DrawLow,bins_var_DrawHigh);
	h_DdlErr_sideband->Draw("same");

	h_DdlErr_DataSig->SetLineColor(1);
	h_DdlErr_DataSig->Rebin(8);
	h_DdlErr_DataSig->GetXaxis()->SetRangeUser(bins_var_DrawLow,bins_var_DrawHigh);
	h_DdlErr_DataSig->Draw("same");
	*/

	//	return 1;
	/*
	*/


	TFile *f_NonPromptDs=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BtoDs_Results/output/BtoDs.root");
	// TH1D *hBtoDsCrossSectionPP_AnaBin=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin");
	TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight");

	// TH1D *hBtoDsdNdPtPbPb_AnaBin=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin");
	TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight");

	TH1D *hBtoDs_AnaBin_pythiaWeight=hBtoDsCrossSectionPP_AnaBin_pythiaWeight;  // this is differential cross section
	if(isPbPb){ hBtoDs_AnaBin_pythiaWeight=hBtoDsdNdPtPbPb_AnaBin_pythiaWeight; }

	MutiplyBinWidth(hBtoDs_AnaBin_pythiaWeight);
	double CS_Integral=hBtoDs_AnaBin_pythiaWeight->Integral(hBtoDs_AnaBin_pythiaWeight->GetXaxis()->FindBin(Dpt_Low+0.00001), hBtoDs_AnaBin_pythiaWeight->GetXaxis()->FindBin(Dpt_High-0.00001) );

	cout<<"CS_Integral = "<<CS_Integral<<endl;


	TString mcName_Prompt_new=Form("/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_sPlot/%sMC_phiPrompt_fitFile_smear.root",s_PbPb.Data());
	TString mcName_NonPrompt_new=Form("/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_sPlot/%sMC_phiNonPrompt_fitFile_smear.root",s_PbPb.Data());
	////////////////////////
	/// create smear Tree //
	////////////////////////

	TFile *f_mcp_new=NULL;
	TFile *f_mcnp_new=NULL;

	// double step=0.01;
	double step=0.1;
	double scaleHi=1.30;
	double scaleLow=1.0;

	if(Dpt_High<=3){
		scaleHi=1.20;
		scaleLow=0.90;
	}
	if( isPbPb==0 && Dpt_High>20){
		scaleHi=1.38;
		scaleLow=1.15;
	}

	int nSclErrF=(int)((scaleHi-scaleLow)/step);
	Float_t ScaleErrFactorArr[nSclErrF];
	Float_t ScaleErrChi2[nSclErrF];
	double Best_Scale_PromptEff=0;
	double Best_Scale_NonPromptEff=0;

	for(int i=0; i<nSclErrF;i++){
		ScaleErrFactorArr[i]=scaleLow+i*step;
	}
	// Float_t ScaleErrFactorArr[]={1.2,1.3,1.4};
	// Float_t ScaleErrFactorArr[]={0.7,0.8,0.85,0.9,0.95,0.98,1,1.02,1.05,1.1,1.15,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9};
	// const int nSclErrF=sizeof(ScaleErrFactorArr)/sizeof(ScaleErrFactorArr[0]);

	// doCreateSMTree=1;
	if(doCreateSMTree){
		double DptLow_tune=20;
		double DptHigh_tune=40;
		DptLow_tune=Dpt_Low;
		DptHigh_tune=Dpt_High;

		if(isPbPb){
			DptLow_tune=6;
		}
		f_mcp_new=new TFile(mcName_Prompt_new.Data(),"recreate");
		smearTree(f_mcp_new,f_mc_Prompt,nSclErrF, ScaleErrFactorArr,nSmear,DptLow_tune,DptHigh_tune);
		f_mcnp_new=new TFile(mcName_NonPrompt_new.Data(),"recreate");
		smearTree(f_mcnp_new,f_mc_NonPrompt,nSclErrF,ScaleErrFactorArr,nSmear,DptLow_tune,DptHigh_tune);

	}else{
		f_mcp_new=new TFile(mcName_Prompt_new.Data(),"read");
		f_mcnp_new=new TFile(mcName_NonPrompt_new.Data(),"read");
	}
	TTree *t_mcp_new=(TTree*)f_mcp_new->Get("t_fit");
	TTree *t_mcnp_new=(TTree*)f_mcnp_new->Get("t_fit");

	// t_mcp_new->Draw("Dpt");

	double Best_Scale=0;
	double Best_Scale_index=0;
	double Best_Scale_Prob=0;
	double Best_Scale_Chi2=0;
	double NonPrompt_phi_Eff=0;
	double Prompt_phi_Eff=0;


	for(int m=0;m<nSclErrF;m++)
	// for(int m=0;m<0;m++)
	{

		// cout<<__LINE__<<endl;
		TString s_Ddls=Form("Ddls_Scl%.0fem3",ScaleErrFactorArr[m]*1e3);
		TString s_Ddls_text=Form("Ddls Err. Scale : %.3f",ScaleErrFactorArr[m]);
		TString s_SMrWt="SMrWt";

		TString s_var_inScan=var_compare.Data();
		if(var_compare=="Ddls"){
			s_var_inScan=Form("Ddls_Scl%.0fem3",ScaleErrFactorArr[m]*1e3);
		}else if (var_compare=="DdlErr"){
			s_var_inScan=Form("DdlErr_Scl%.0fem3",ScaleErrFactorArr[m]*1e3);
		} else if (var_compare=="Ddl"){
			s_var_inScan=Form("Ddl_Scl%.0fem3",ScaleErrFactorArr[m]*1e3);
		}


		// cout<<__LINE__<<endl;

		cout<<"m = "<<m<<" s_Ddls_text ="<<s_Ddls_text<<endl;

		TString DataCuts_Smear=Form("Dpt>%f && Dpt<%f && Dalpha < %f && Dchi2cl > %f && %s >%f && %s > %f && %s <%f", Dpt_Low, Dpt_High, Dalpha_cut, Dchi2cl_cut, s_Ddls.Data(),Ddls_cut, var_cut.Data(), var_cutLow,var_cut.Data(), var_cutHigh);
		// build prompt /nonprompt MC historgram 
		TH1D *h_var_PromptMC_scl =new TH1D(Form("h_%s_PromptMC_scl",var_compare.Data()), Form("h_%s_PromptMC_scl",var_compare.Data()), nbin_var, bins_var_Low,bins_var_High); h_var_PromptMC_scl->Sumw2();
		t_mcp_new->Project(Form("h_%s_PromptMC_scl",var_compare.Data()),s_var_inScan.Data(),Form("(Dmass>%f && Dmass<%f && %s)*%s",Dmass2SigLow, Dmass2SigHigh, DataCuts_Smear.Data(),MC_DefaultWeight.Data() ));

		h_var_PromptMC_scl->Scale(1/h_var_PromptMC_scl->Integral());

		TH1D *h_var_NonPromptMC_scl =new TH1D(Form("h_%s_NonPromptMC_scl",var_compare.Data()), Form("h_%s_NonPromptMC_scl",var_compare.Data()), nbin_var, bins_var_Low,bins_var_High); h_var_NonPromptMC_scl->Sumw2();
		t_mcnp_new->Project(Form("h_%s_NonPromptMC_scl",var_compare.Data()),s_var_inScan.Data(),Form("(Dmass>%f && Dmass<%f && %s)*%s", Dmass2SigLow, Dmass2SigHigh, DataCuts_Smear.Data(),MC_DefaultWeight.Data()));

		h_var_NonPromptMC_scl->Scale(1/h_var_NonPromptMC_scl->Integral());


		h_var_PromptMC_scl->SetLineColor(2);
		h_var_NonPromptMC_scl->Draw("same");
		h_var_NonPromptMC_scl->SetLineColor(4);
		h_var_PromptMC_scl->Draw("same");


		// cout<<__LINE__<<endl;
		// prompt fraction cal
		double fr_withcut_prompt=0.85; // just for get quick result
		double fr_prompt=0.85; // just for get quick result

		// extract fr_prompt from fit and BtoDs estimation
		bool use_fr_fromBtoD=true;

		if(use_fr_fromBtoD){
			double N_yield=h_DataRawYield->GetBinContent(1);
			double LumiNevt=LumiSum;
			if(isPbPb){LumiNevt=NevtPbPb3;}

			TH1D *h_eff_nonprompt_phi=new TH1D("h_eff_nonprompt_phi","h_eff_nonprompt_phi",40,0,40);  h_eff_nonprompt_phi->Sumw2();
			t_mcnp_new->Project("h_eff_nonprompt_phi","Dpt",(TCut)Form("%s*%s*(%s)",MC_DefaultWeight.Data(),s_SMrWt.Data(),DataCuts_Smear.Data()));

			double NonPrompt_phi_reco=h_eff_nonprompt_phi->Integral(h_eff_nonprompt_phi->FindBin(Dpt_Low+0.00001), h_eff_nonprompt_phi->FindBin(Dpt_High-0.00001));
			TH1D *hGen_pt_nonprompt_phi=(TH1D*)f_mc_NonPrompt->Get(Form("hGen_pt%s",GenWt.Data() ) );
			double NonPrompt_phi_Gen=hGen_pt_nonprompt_phi->Integral(hGen_pt_nonprompt_phi->FindBin(Dpt_Low+0.00001),hGen_pt_nonprompt_phi->FindBin(Dpt_High-0.0001) );

			NonPrompt_phi_Eff=NonPrompt_phi_reco/NonPrompt_phi_Gen;

			cout<<"NonPrompt_phi_Eff = "<<NonPrompt_phi_Eff<<endl;


			TH1D *h_eff_prompt_phi=new TH1D("h_eff_prompt_phi","h_eff_prompt_phi",40,0,40);  h_eff_prompt_phi->Sumw2();
			t_mcp_new->Project("h_eff_prompt_phi","Dpt",(TCut)Form("%s*%s*(%s)",MC_DefaultWeight.Data(),s_SMrWt.Data(),DataCuts_Smear.Data()));

			double Prompt_phi_reco=h_eff_prompt_phi->Integral(h_eff_prompt_phi->FindBin(Dpt_Low+0.00001), h_eff_prompt_phi->FindBin(Dpt_High-0.00001));
			TH1D *hGen_pt_prompt_phi=(TH1D*)f_mc_Prompt->Get(Form("hGen_pt%s",GenWt.Data() ) );
			double Prompt_phi_Gen=hGen_pt_prompt_phi->Integral(hGen_pt_prompt_phi->FindBin(Dpt_Low+0.00001),hGen_pt_prompt_phi->FindBin(Dpt_High-0.0001) );

			Prompt_phi_Eff=Prompt_phi_reco/Prompt_phi_Gen;

			cout<<"Prompt_phi_Eff = "<<Prompt_phi_Eff<<endl;




			//    cout<<"NonPrompt_f0_Eff = "<<NonPrompt_f0_Eff<<endl;
			double phiFr=0.95 ;
			if(isPbPb){
				phiFr=0.89;
			}

			double N_NonPrompt_yield=CS_Integral*(2*LumiNevt)*( BRphi*NonPrompt_phi_Eff)/phiFr;
			// double N_NonPrompt_yield=CS_Integral*(2*LumiNevt)*( BRphi*NonPrompt_phi_Eff+ BRf0*NonPrompt_f0_Eff );

			fr_withcut_prompt=1-(N_NonPrompt_yield/N_yield);
			fr_prompt=(fr_withcut_prompt/Prompt_phi_Eff)/(fr_withcut_prompt/Prompt_phi_Eff + (1-fr_withcut_prompt)/NonPrompt_phi_Eff );


			cout<<"fr_withcut_prompt = "<<fr_withcut_prompt<<" , N_NonPrompt_yield = "<<N_NonPrompt_yield<<" , N_yield = "<<N_yield<<endl;
			cout<<"fr_prompt = "<<fr_prompt<<endl;

			delete h_eff_nonprompt_phi;
			delete h_eff_prompt_phi;

		} // end fr_prompt calculation

		TH1D *h_var_MixMC= new TH1D(Form("h_%s_MixMC",var_compare.Data()),"h_var_MixMC",nbin_var, bins_var_Low,bins_var_High); h_var_MixMC->Sumw2();

		h_var_MixMC->Add(h_var_PromptMC_scl,h_var_NonPromptMC_scl,fr_withcut_prompt,1-fr_withcut_prompt);

		//      h_var_MixMC->Draw("same");
		// cout<<__LINE__<<endl;


		TCanvas *c_var_compare= new TCanvas("c_var_compare","c_var_compare",800,800);
		c_var_compare->cd();

		TPad *pad1 = new TPad("pad1","top pad", 0.0,0.25,1.0,1.0);
		pad1->SetBottomMargin(0.0);
		pad1->Draw();
		TPad *pad2 = new TPad("pad2","bottom pad",0.0,0.0,1.0,0.25);
		pad2->SetTopMargin(0.0);
		pad2->SetBottomMargin(0.30);
		pad2->Draw();

		pad1->cd();

		// cout<<__LINE__<<endl;

		TH1D *h_var_Data2sig_temp=(TH1D*)h_DdlErr_DataSig->Clone("h_var_Data2sig_temp");

		h_var_Data2sig_temp->Rebin(rebinN);
		h_var_PromptMC_scl->Rebin(rebinN);
		h_var_NonPromptMC_scl->Rebin(rebinN);
		h_var_MixMC->Rebin(rebinN);

		// cout<<__LINE__<<endl;

		h_var_Data2sig_temp->GetXaxis()->SetRangeUser(bins_var_DrawLow,bins_var_DrawHigh); // normalized before setRange   
		h_var_Data2sig_temp->GetXaxis()->SetTitle(var_compare.Data());
		h_var_Data2sig_temp->SetTitle("");
		h_var_Data2sig_temp->SetLineColor(1);
		h_var_Data2sig_temp->SetMarkerColor(1);
		h_var_Data2sig_temp->SetMarkerStyle(22);
		h_var_Data2sig_temp->Draw();

		// h_DdlErr_Data2sigMix->SetLineColor(2);
		// h_DdlErr_Data2sigMix->Draw("same");      

		// return 1;

		h_var_PromptMC_scl->GetXaxis()->SetRangeUser(bins_var_DrawLow,bins_var_DrawHigh);
		h_var_PromptMC_scl->SetLineColor(2);
		h_var_PromptMC_scl->SetMarkerColor(2);
		h_var_PromptMC_scl->SetMarkerStyle(26);
		h_var_PromptMC_scl->Draw("same");
		h_var_NonPromptMC_scl->SetLineColor(4); 
		h_var_NonPromptMC_scl->SetMarkerColor(4);
		h_var_NonPromptMC_scl->SetMarkerStyle(26);
		h_var_NonPromptMC_scl->Draw("same");
		h_var_MixMC->SetLineColor(1);
		h_var_MixMC->SetMarkerColor(1);
		h_var_MixMC->SetMarkerStyle(26);
		h_var_MixMC->Draw("same");


		// cout<<__LINE__<<endl;

		TLegend *le_var = new TLegend(0.65,0.37,0.88,0.62);
		le_var->SetBorderSize(0);
		le_var->AddEntry(h_var_Data2sig_temp,"Data","pl");
		le_var->AddEntry(h_var_PromptMC_scl,"MC Prompt D_{S}","pl");
		le_var->AddEntry(h_var_NonPromptMC_scl,"MC NonPrompt D_{S}","pl");
		le_var->AddEntry(h_var_MixMC,"MC Mix D_{S}","pl");
		le_var->Draw("same");

		shiftY=0.035;
		tltx->SetTextSize(0.03);
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Ds %s ",Str_PbPb.Data(), var_compare.Data())); shiftY-=oneshift;
		// tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Ds %s %.0f<p_{T}<%.0f",Str_PbPb.Data(), var_compare.Data(),Dpt_Low,Dpt_High)); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Var : %s %.2f to %.2f ", var_cut.Data(), var_cutLow, var_cutHigh)); shiftY-=oneshift;
		// tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Prompt fraction (CS) : %.2f ", fr_prompt)); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("DdlErr Scale : %.2f ", ScaleErrFactorArr[m] )); shiftY-=oneshift;

		pad2->cd();

		TH1D *h_ratio=(TH1D*)h_var_Data2sig_temp->Clone("h_ratio");
		h_ratio->Divide(h_var_MixMC);
		h_ratio->GetYaxis()->SetTitle("Data/MC");
		h_ratio->GetYaxis()->CenterTitle();
		h_ratio->GetYaxis()->SetLabelSize(0.05);
		h_ratio->GetYaxis()->SetTitleSize(0.08);
		h_ratio->GetYaxis()->SetTitleOffset(0.38);
		h_ratio->GetXaxis()->SetLabelSize(0.12);
		h_ratio->GetXaxis()->SetTitleSize(0.12);
		h_ratio->SetMaximum(2.9);
		h_ratio->SetMinimum(0.1);
		h_ratio->Draw();

		TF1 *f1_line=new TF1("f1_line","[0]");
		f1_line->SetParameter(0,1);
		f1_line->FixParameter(0,1);
		f1_line->SetRange(bins_var_DrawLow,bins_var_DrawHigh);
		h_ratio->Fit("f1_line","q");
		f1_line->Draw("same");

		// cout<<__LINE__<<endl;

		gSystem->Exec(Form("mkdir -p plots_smear/%s",var_compare.Data()));
		c_var_compare->SaveAs(Form("./plots_smear/%s/%s_%s%.0fto%.0f_DdlErrScale%.0fem3.png", var_compare.Data(),Str_PbPb.Data(), var_cut.Data(),var_cutLow*100,var_cutHigh*100, ScaleErrFactorArr[m]*1e3));
		c_var_compare->SaveAs(Form("./plots_smear/%s/C/%s_%s%.0fto%.0f_DdlErrScale%.0fem3.C", var_compare.Data(),Str_PbPb.Data(), var_cut.Data(),var_cutLow*100,var_cutHigh*100, ScaleErrFactorArr[m]*1e3));

		TCanvas *c_var_compareRatio=new TCanvas("c_var_compareRatio","",600,600);
		c_var_compareRatio->cd();

		h_ratio->SetMaximum(3);
		h_ratio->SetMinimum(-1);
		h_ratio->GetYaxis()->SetLabelSize(0.03);
		h_ratio->GetYaxis()->SetTitleSize(0.04);
		h_ratio->GetYaxis()->SetTitleOffset(1.0);
		h_ratio->GetXaxis()->SetTitleOffset(1.0);
		h_ratio->GetXaxis()->SetLabelSize(0.04);
		h_ratio->GetXaxis()->SetTitleSize(0.04);

		h_ratio->Draw();
		f1_line->Draw("same");

		cout<<"DdlErr Scale : "<<ScaleErrFactorArr[m]<<" ,GetNDF() = "<<f1_line->GetNDF()<<" GetChisquare() = "<<f1_line->GetChisquare()<<" prob = "<<TMath::Prob(f1_line->GetChisquare(),f1_line->GetNDF())<<endl;

		// c_var_compareRatio->SaveAs(Form("./plots_smear/%s/%s_Ratio_%s%.0fto%.0f_DdlErrScale%.0fem3.png", var_compare.Data(),Str_PbPb.Data(), var_cut.Data(),var_cutLow*100,var_cutHigh*100, ScaleErrFactorArr[m]*1e3));
		 ScaleErrChi2[m]=f1_line->GetChisquare();

		if(Best_Scale_Prob < TMath::Prob(f1_line->GetChisquare(),f1_line->GetNDF()) ){
			Best_Scale_Prob=TMath::Prob(f1_line->GetChisquare(),f1_line->GetNDF());
			Best_Scale_Chi2=f1_line->GetChisquare();
			Best_Scale=ScaleErrFactorArr[m];
			Best_Scale_index=m;
			Best_Scale_PromptEff=Prompt_phi_Eff;
			Best_Scale_NonPromptEff=NonPrompt_phi_Eff;

		}


		// cout<<__LINE__<<endl;
		// return 1;

		// delete f1_line;
		delete le_var;
		delete h_ratio;
		delete h_var_MixMC;
		delete h_var_Data2sig_temp;
		delete h_var_PromptMC_scl;
		delete h_var_NonPromptMC_scl;
		// delete c_var_compare;
		delete c_var_compareRatio;

	} // end for m<nSclErrF

	cout<<"Best_Scale = "<<Best_Scale<<" , Best Scale Prob = "<<Best_Scale_Prob<<" , Best_Scale_Chi2 = " <<Best_Scale_Chi2<<" index = "<<Best_Scale_index<<endl;
	cout<<"Best_Scale Prompt Eff = "<<Best_Scale_PromptEff<<" , NonPrompt Eff = "<<Best_Scale_NonPromptEff<<endl;	

	// loop for chi2+1 scale
	int chi2InRange=0;
	double Best_Scale_Low=0;
	double Best_Scale_High=0;

	TGraph *gr_scale_chi2=new TGraph(nSclErrF,ScaleErrFactorArr,ScaleErrChi2);
	gr_scale_chi2->SetName("gr_scale_chi2");
	gr_scale_chi2->GetXaxis()->SetTitle("DdlErr. Scale");
	gr_scale_chi2->GetYaxis()->SetTitle("#Chi2^{2}");
	fout->cd();
	gr_scale_chi2->Write("",TObject::kOverwrite);

	for(int m=0; m<nSclErrF; m++){
		if(ScaleErrChi2[m]<Best_Scale_Chi2+1 && chi2InRange==0){
				Best_Scale_Low=ScaleErrFactorArr[m-1];
				chi2InRange=1;
		}
		if(ScaleErrChi2[m]>Best_Scale_Chi2+1 && chi2InRange==1){
				Best_Scale_High=ScaleErrFactorArr[m];
				break;
		}
	}

	cout<<"Best_Scale_Low = "<<Best_Scale_Low<<" , Best_Scale_High = "<<Best_Scale_High<<endl;

	// return 1;

	TH1D *h_Best_Scale=new TH1D("h_Best_Scale","",1,0,1);
	h_Best_Scale->SetBinContent(1,Best_Scale);

	TGraphAsymmErrors *gr_Best_Scale=new TGraphAsymmErrors();
	gr_Best_Scale->SetName("gr_Best_Scale");
	gr_Best_Scale->SetPoint(0,(Dpt_Low+Dpt_High)/2,Best_Scale);
	gr_Best_Scale->SetPointError(0 , (Dpt_Low+Dpt_High)/2-Dpt_Low , (Dpt_Low+Dpt_High)/2-Dpt_Low , Best_Scale-Best_Scale_Low , Best_Scale_High-Best_Scale);


	fout->cd();
	gr_Best_Scale->Write("",TObject::kOverwrite);
	h_Best_Scale->Write("",TObject::kOverwrite);

	// end DdlErrScale Scan part


	// 

	// Best_Scale_index=49;
	int m=Best_Scale_index;
	// start Ddls PNP fraction part

	var_compare="Ddls";
	bins_var_Low=0;
	bins_var_High=200;
	bins_var_DrawLow=0;
	bins_var_DrawHigh=30;
	if(Dpt_High<=6){
		bins_var_DrawHigh=20;
	}

	TH1D *h_Ddls_sideband=ProjectData(var_compare, DataCuts_sideband, t_Ds_Data, nbin_var, bins_var_Low, bins_var_High);
	TH1D *h_Ddls_Data2sigMix=ProjectData(var_compare, DataCuts_2sig, t_Ds_Data, nbin_var, bins_var_Low, bins_var_High);

	TH1D *h_Ddls_DataSig=(TH1D*)h_Ddls_Data2sigMix->Clone("h_Ddls_DataSig");
	// h_Ddls_DataSig->Scale(1/h_Ddls_DataSig->Integral());
	h_Ddls_DataSig->Add(h_Ddls_sideband,-(1-sigFrac));
	h_Ddls_DataSig->Scale(1/h_Ddls_DataSig->Integral());

	h_Ddls_DataSig->Rebin(rebinN);
	h_Ddls_DataSig->GetXaxis()->SetRangeUser(bins_var_DrawLow,bins_var_DrawHigh);
	h_Ddls_DataSig->SetLineColor(1);
	h_Ddls_DataSig->Draw();

	h_Ddls_sideband->Rebin(rebinN);
	h_Ddls_sideband->GetXaxis()->SetRangeUser(bins_var_DrawLow,bins_var_DrawHigh);
	h_Ddls_Data2sigMix->Rebin(rebinN);
	h_Ddls_Data2sigMix->GetXaxis()->SetRangeUser(bins_var_DrawLow,bins_var_DrawHigh);
	



	TString s_Ddls=Form("Ddls_Scl%.0fem3",ScaleErrFactorArr[m]*1e3);
	TString s_Ddls_text=Form("Ddls Err. Scale : %.2f",ScaleErrFactorArr[m]);
	TString s_SMrWt="SMrWt";

	TString s_var_inScan=var_compare.Data();
	if(var_compare=="Ddls"){
		s_var_inScan=Form("Ddls_Scl%.0fem3",ScaleErrFactorArr[m]*1e3);
	}else if (var_compare=="DdlErr"){
		s_var_inScan=Form("DdlErr_Scl%.0fem3",ScaleErrFactorArr[m]*1e3);
	} else if (var_compare=="Ddl"){
		s_var_inScan=Form("Ddl_Scl%.0fem3",ScaleErrFactorArr[m]*1e3);
	}


	// cout<<__LINE__<<endl;

	cout<<"m = "<<m<<" s_Ddls_text ="<<s_Ddls_text<<endl;

	TString DataCuts_Smear=Form("Dpt>%f && Dpt<%f && Dalpha < %f && Dchi2cl > %f && %s >%f && %s > %f && %s <%f", Dpt_Low, Dpt_High, Dalpha_cut, Dchi2cl_cut, s_Ddls.Data(),Ddls_cut, var_cut.Data(), var_cutLow,var_cut.Data(), var_cutHigh);
	// build prompt /nonprompt MC historgram 
	TH1D *h_var_PromptMC_scl =new TH1D(Form("h_%s_PromptMC_scl",var_compare.Data()), Form("h_%s_PromptMC_scl",var_compare.Data()), nbin_var, bins_var_Low,bins_var_High); h_var_PromptMC_scl->Sumw2();
	t_mcp_new->Project(Form("h_%s_PromptMC_scl",var_compare.Data()),s_var_inScan.Data(),Form("(Dmass>%f && Dmass<%f && %s)*%s",Dmass2SigLow, Dmass2SigHigh, DataCuts_Smear.Data(),MC_DefaultWeight.Data() ));

	h_var_PromptMC_scl->Scale(1/h_var_PromptMC_scl->Integral());

	TH1D *h_var_NonPromptMC_scl =new TH1D(Form("h_%s_NonPromptMC_scl",var_compare.Data()), Form("h_%s_NonPromptMC_scl",var_compare.Data()), nbin_var, bins_var_Low,bins_var_High); h_var_NonPromptMC_scl->Sumw2();
	t_mcnp_new->Project(Form("h_%s_NonPromptMC_scl",var_compare.Data()),s_var_inScan.Data(),Form("(Dmass>%f && Dmass<%f && %s)*%s", Dmass2SigLow, Dmass2SigHigh, DataCuts_Smear.Data(),MC_DefaultWeight.Data()));

	h_var_NonPromptMC_scl->Scale(1/h_var_NonPromptMC_scl->Integral());


	h_var_PromptMC_scl->SetLineColor(2);
	h_var_NonPromptMC_scl->Rebin(rebinN);
	// h_var_NonPromptMC_scl->Draw("same");
	h_var_NonPromptMC_scl->SetLineColor(4);
	h_var_PromptMC_scl->Rebin(rebinN);
	h_var_PromptMC_scl->Draw("same");

	hDsDlsDataBkg=(TH1D*)h_Ddls_sideband->Clone("hDsDlsDataBkg");
	hDsDlsMCPSignal=(TH1D*)h_var_PromptMC_scl->Clone("hDsDlsMCPSignal");
	hDsDlsMCNPSignal=(TH1D*)h_var_NonPromptMC_scl->Clone("hDsDlsMCNPSignal");

	hDsDlsMCPSignal->Scale(1/hDsDlsMCPSignal->Integral(hDsDlsMCPSignal->FindBin(bins_var_DrawLow),hDsDlsMCPSignal->FindBin(bins_var_DrawHigh)));
	hDsDlsMCNPSignal->Scale(1/hDsDlsMCNPSignal->Integral(hDsDlsMCNPSignal->FindBin(bins_var_DrawLow),hDsDlsMCNPSignal->FindBin(bins_var_DrawHigh)));


	// normalize(hDsDlsMCPSignal); // this consider unequal binwidth
	// normalize(hDsDlsMCNPSignal);

	// hDsDlsMCPSignal->SetLineColor(4);
	// hDsDlsMCPSignal->Draw("same");

	cout<<"hDsDlsMCPSignal Integral()"<<hDsDlsMCPSignal->Integral(hDsDlsMCPSignal->FindBin(0),hDsDlsMCPSignal->FindBin(30))<<endl;



	TF1* fMix = new TF1("fMix",&funMix, bins_var_DrawLow, bins_var_DrawHigh, 1);
	// fMix->SetParameter(0,0.9);
	fMix->SetParLimits(0,0,1.0);
	// fMix->FixParameter(0,0.9);
	// fMix->SetParLimits(1,0,2*integralTotalYield);
	fMix->SetLineColor(2);
	// fMix->SetFillColor(kRed-9);
	fMix->SetFillStyle(1001);

	// float fitRangeL = 0.;
	// float fitRangeH = 0.08;

	// hD0DcaData->GetXaxis()->SetRangeUser(0.,0.07);
	// hD0DcaData->SetMaximum(hD0DcaData->GetMaximum()*40);

	TCanvas *c_Ddls_PfrFit=new TCanvas("c_Ddls_PfrFit","",800,800);
	c_Ddls_PfrFit->cd();

	// h->Draw("p");
	h_Ddls_DataSig->Scale(1/h_Ddls_DataSig->Integral(h_Ddls_DataSig->FindBin(bins_var_DrawLow), h_Ddls_DataSig->FindBin(bins_var_DrawHigh) ));

	TH1D *h_Ddls_DataSig_temp=(TH1D*)h_Ddls_DataSig->Clone("h_Ddls_DataSig_temp");
	h_Ddls_DataSig_temp->Draw();

	h_Ddls_DataSig_temp->Fit("fMix","SNQ0","",bins_var_DrawLow,bins_var_DrawHigh);
	h_Ddls_DataSig_temp->Fit("fMix","SNQ0","",bins_var_DrawLow,bins_var_DrawHigh);
	h_Ddls_DataSig_temp->Fit("fMix","SNQ0","",bins_var_DrawLow,bins_var_DrawHigh);
	h_Ddls_DataSig_temp->Fit("fMix","SNQ0","",bins_var_DrawLow,bins_var_DrawHigh);
	h_Ddls_DataSig_temp->Fit("fMix","SNQ0","",bins_var_DrawLow,bins_var_DrawHigh);
	h_Ddls_DataSig_temp->Fit("fMix","","",bins_var_DrawLow,bins_var_DrawHigh);


	fMix->Draw("same");


	TF1* fNP = new TF1("fNP",&funNonPrompt, bins_var_DrawLow, bins_var_DrawHigh, 1);
	fNP->SetParameter(0,fMix->GetParameter(0));
	fNP->SetRange(bins_var_DrawLow,bins_var_DrawHigh);
	fNP->SetLineColor(4);
	fNP->SetFillStyle(1001);
	fNP->SetFillColor(kBlue-9);
	fNP->SetNpx(10000);
	fNP->Draw("same");

	double Pfr_withcut=fMix->GetParameter(0);
	double Pfr_cs=(Pfr_withcut/Best_Scale_PromptEff) / (Pfr_withcut/Best_Scale_PromptEff+ (1-Pfr_withcut)/Best_Scale_NonPromptEff );
	double Pfr_csErr=fMix->GetParError(0)/Pfr_withcut*Pfr_cs;

	gStyle->SetOptFit(0);

		shiftY=0.030;
		tltx->SetTextSize(0.027);
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Ds %s ",Str_PbPb.Data(), var_compare.Data())); shiftY-=oneshift;
		// tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Ds %s %.0f<p_{T}<%.0f",Str_PbPb.Data(), var_compare.Data(),Dpt_Low,Dpt_High)); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Var : %s %.2f to %.2f ", var_cut.Data(), var_cutLow, var_cutHigh)); shiftY-=oneshift;
		// tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Prompt fraction : %.2f #pm %.2f", fMix->GetParameter(0),fMix->GetParError(0))); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Prompt fraction (CS) : %.2f #pm %.2f", Pfr_cs,Pfr_csErr ) ); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("DdlErr Scale : %.2f ", ScaleErrFactorArr[m] )); shiftY-=oneshift;


	c_Ddls_PfrFit->SaveAs(Form("./plots_smear/%s/%s_pfrFit_%s%.0fto%.0f_DdlErrScale%.0fem3.png", var_compare.Data(),Str_PbPb.Data(), var_cut.Data(),var_cutLow*100,var_cutHigh*100, ScaleErrFactorArr[m]*1e3));

	TGraphErrors *gr_Fit_PNP=new TGraphErrors();
	gr_Fit_PNP->SetName("gr_Fit_PNP");
	gr_Fit_PNP->SetPoint(0,(Dpt_Low+Dpt_High)/2, Pfr_cs);
	gr_Fit_PNP->SetPointError(0,(Dpt_Low+Dpt_High)/2-Dpt_Low, Pfr_csErr);

	fout->cd();
	gr_Fit_PNP->Write("",TObject::kOverwrite);

	// include Bkg fit

	TF1* fMixWBkg = new TF1("fMixWBkg",&funMixWBkg, bins_var_DrawLow, bins_var_DrawHigh, 2);
	fMixWBkg->SetParameter(0,0.5);
	fMixWBkg->SetParameter(1,0.9);
	fMixWBkg->SetParLimits(0,0,1.0);
	fMixWBkg->SetParLimits(1,0,1.0);
	fMixWBkg->SetLineColor(2);
	// fMixWBkg->SetFillColor(kRed-9);
	fMixWBkg->SetFillStyle(1001);


	TCanvas *c_Ddls_PfrFitWBkg=new TCanvas("c_Ddls_PfrFitWBkg","",800,800);
	c_Ddls_PfrFitWBkg->cd();

	// h->Draw("p");
	h_Ddls_Data2sigMix->Scale(1/h_Ddls_Data2sigMix->Integral(h_Ddls_Data2sigMix->FindBin(bins_var_DrawLow), h_Ddls_Data2sigMix->FindBin(bins_var_DrawHigh) ));
	h_Ddls_sideband->Scale(1/h_Ddls_sideband->Integral(h_Ddls_sideband->FindBin(bins_var_DrawLow), h_Ddls_sideband->FindBin(bins_var_DrawHigh) ));

	TH1D *h_Ddls_Data2SigMix_temp=(TH1D*)h_Ddls_Data2sigMix->Clone("h_Ddls_Data2SigMix_temp");
	h_Ddls_Data2SigMix_temp->Draw();

	h_Ddls_Data2SigMix_temp->Fit("fMixWBkg","SNQ0","",bins_var_DrawLow,bins_var_DrawHigh);
	h_Ddls_Data2SigMix_temp->Fit("fMixWBkg","SNQ0","",bins_var_DrawLow,bins_var_DrawHigh);
	h_Ddls_Data2SigMix_temp->Fit("fMixWBkg","SNQ0","",bins_var_DrawLow,bins_var_DrawHigh);
	h_Ddls_Data2SigMix_temp->Fit("fMixWBkg","SNQ0","",bins_var_DrawLow,bins_var_DrawHigh);
	h_Ddls_Data2SigMix_temp->Fit("fMixWBkg","SNQ0","",bins_var_DrawLow,bins_var_DrawHigh);
	h_Ddls_Data2SigMix_temp->Fit("fMixWBkg","","",bins_var_DrawLow,bins_var_DrawHigh);

	fMixWBkg->Draw("same");

	TF1* fSigDdls = new TF1("fSigDdls",&funSigDdls, bins_var_DrawLow, bins_var_DrawHigh, 1);
	fSigDdls->SetParameter(0,fMixWBkg->GetParameter(0));
	fSigDdls->SetParameter(1,fMixWBkg->GetParameter(1));
	fSigDdls->SetRange(bins_var_DrawLow,bins_var_DrawHigh);
	fSigDdls->SetLineColor(4);
	fSigDdls->SetFillStyle(1001);
	fSigDdls->SetFillColor(kBlue-9);
	fSigDdls->SetNpx(10000);
	fSigDdls->Draw("same");


		shiftY=0.030;
		tltx->SetTextSize(0.027);
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Ds %s ",Str_PbPb.Data(), var_compare.Data())); shiftY-=oneshift;
		// tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Ds %s %.0f<p_{T}<%.0f",Str_PbPb.Data(), var_compare.Data(),Dpt_Low,Dpt_High)); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Var : %s %.2f to %.2f ", var_cut.Data(), var_cutLow, var_cutHigh)); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Prompt fraction : %.2f #pm %.2f", fMixWBkg->GetParameter(1),fMixWBkg->GetParError(1))); shiftY-=oneshift;
		// tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Prompt fraction (CS) : %.2f #pm %.2f", Pfr_cs,Pfr_csErr ) ); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("DdlErr Scale : %.2f ", ScaleErrFactorArr[m] )); shiftY-=oneshift;


	c_Ddls_PfrFitWBkg->SaveAs(Form("./plots_smear/%s/%s_pfrFitWBkg_%s%.0fto%.0f_DdlErrScale%.0fem3.png", var_compare.Data(),Str_PbPb.Data(), var_cut.Data(),var_cutLow*100,var_cutHigh*100, ScaleErrFactorArr[m]*1e3));


	// return 1;


	// TFractionFitter
	TObjArray *mc_Ddls=new TObjArray(2);
	mc_Ddls->Add(hDsDlsMCPSignal);
	mc_Ddls->Add(hDsDlsMCNPSignal);

	TH1D *h_Ddls_DataSig_temp2=(TH1D*)h_Ddls_DataSig->Clone("h_Ddls_DataSig_temp2");
	TFractionFitter *fr_fit=new TFractionFitter(h_Ddls_DataSig_temp2,mc_Ddls);
	// fr_fit->Constant(1.0,1.0);
	fr_fit->SetRangeX(1,15);

	fr_fit->Fit();
	fr_fit->Fit();
	fr_fit->Fit();

	TCanvas *c_TFracFit=new TCanvas("c_TFracFit","",600,600);
	c_TFracFit->cd();

	TH1D* result =(TH1D*) fr_fit->GetPlot();
	h_Ddls_DataSig_temp2->Draw("Ep");
	result->Draw("same")	;	

	c_TFracFit->SaveAs(Form("./plots_smear/%s/%s_pfrTFrFit_%s%.0fto%.0f_DdlErrScale%.0fem3.png", var_compare.Data(),Str_PbPb.Data(), var_cut.Data(),var_cutLow*100,var_cutHigh*100, ScaleErrFactorArr[m]*1e3));
// this doesn't right

	// manually scan

	double PfrLow=0.7;
	double PfrStep=0.01;	
	double Best_Pfr=0;
	double Best_Pfr_Prob=0;

	// return 1;

	for(int i=0; i<31;i++){

		gStyle->SetOptFit(1111);

		double Pfr=PfrLow+i*PfrStep;
		double Pfr_withCut_scan= (Pfr*Best_Scale_PromptEff)/(Pfr*Best_Scale_PromptEff + (1-Pfr)*Best_Scale_NonPromptEff);
		TH1D *h_var_MCMix=(TH1D*)h_var_PromptMC_scl->Clone("h_var_MixMC");
		h_var_MCMix->Add(h_var_PromptMC_scl,h_var_NonPromptMC_scl,Pfr_withCut_scan,1-Pfr_withCut_scan);
		h_var_MCMix->Scale(1/h_var_MCMix->Integral());
		TH1D *h_ratio_Ddls=(TH1D*)h_Ddls_DataSig->Clone("h_ratio_Ddls");
		h_ratio_Ddls->Divide(h_var_MCMix);

		TF1 *f1_line=new TF1("f1_line","[0]");
		f1_line->SetParameter(0,1);
		f1_line->FixParameter(0,1);
		f1_line->SetRange(bins_var_DrawLow,bins_var_DrawHigh);
		h_ratio_Ddls->Fit("f1_line","q");

		TCanvas *c_Ddls_Pfr=new TCanvas("c_Ddls_Pfr","",800,800);
		c_Ddls_Pfr->cd();

		TPad *pad1 = new TPad("pad1","top pad", 0.0,0.25,1.0,1.0);
		pad1->SetBottomMargin(0.0);
		pad1->Draw();
		TPad *pad2 = new TPad("pad2","bottom pad",0.0,0.0,1.0,0.25);
		pad2->SetTopMargin(0.0);
		pad2->SetBottomMargin(0.30);
		pad2->Draw();

		pad1->cd();
		h_Ddls_DataSig->SetMarkerStyle(22);
		h_Ddls_DataSig->SetMarkerColor(1);
		h_Ddls_DataSig->Draw("");
		h_var_MCMix->SetMarkerStyle(26);
		h_var_MCMix->SetMarkerColor(2);
		h_var_MCMix->Draw("same");


		tltx->SetTextSize(0.035);
		shiftY=0.035;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Ds %s ",Str_PbPb.Data(), var_compare.Data())); shiftY-=oneshift;
		// tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Ds %s %.0f<p_{T}<%.0f",Str_PbPb.Data(), var_compare.Data(),Dpt_Low,Dpt_High)); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Var : %s %.2f to %.2f ", var_cut.Data(), var_cutLow, var_cutHigh)); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Prompt fraction (CS) : %.2f ", Pfr)); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("DdlErr Scale : %.2f ", ScaleErrFactorArr[m] )); shiftY-=oneshift;

		TLegend *le_Ddls=new TLegend(0.6,0.35,0.85,0.5);
		le_Ddls->SetBorderSize(0);
		le_Ddls->AddEntry(h_Ddls_DataSig,"Data","pl");
		le_Ddls->AddEntry(h_var_MCMix,"MC","pl");
		le_Ddls->Draw("same");


		pad2->cd();

		h_ratio_Ddls->GetYaxis()->SetTitle("Data/MC");
		h_ratio_Ddls->GetYaxis()->CenterTitle();
		h_ratio_Ddls->GetYaxis()->SetLabelSize(0.05);
		h_ratio_Ddls->GetYaxis()->SetTitleSize(0.08);
		h_ratio_Ddls->GetYaxis()->SetTitleOffset(0.38);
		h_ratio_Ddls->GetXaxis()->SetLabelSize(0.12);
		h_ratio_Ddls->GetXaxis()->SetTitleSize(0.12);
		h_ratio_Ddls->SetMaximum(3);
		h_ratio_Ddls->SetMinimum(-1);
		h_ratio_Ddls->Draw();


		f1_line->Draw("same");


		c_Ddls_Pfr->SaveAs(Form("./plots_smear/%s/%s_%s%.0fto%.0f_DdlErrScale%.0fem3_pfr%.0f.png", var_compare.Data(),Str_PbPb.Data(), var_cut.Data(),var_cutLow*100,var_cutHigh*100, ScaleErrFactorArr[m]*1e3,Pfr*100));


		cout<<"Prompt Fr = "<<Pfr<<" ,GetNDF() = "<<f1_line->GetNDF()<<" GetChisquare() = "<<f1_line->GetChisquare()<<" prob = "<<TMath::Prob(f1_line->GetChisquare(),f1_line->GetNDF())<<endl;


		if(Best_Pfr_Prob < TMath::Prob(f1_line->GetChisquare(),f1_line->GetNDF()) ){
			Best_Pfr_Prob=TMath::Prob(f1_line->GetChisquare(),f1_line->GetNDF());
			Best_Pfr=Pfr;
		}


	} // end for Pfr scan


	TGraphErrors *gr_Best_PNP=new TGraphErrors();
	gr_Best_PNP->SetName("gr_Best_PNP");
	gr_Best_PNP->SetPoint(0,(Dpt_Low+Dpt_High)/2, Best_Pfr);
	gr_Best_PNP->SetPointError(0,(Dpt_Low+Dpt_High)/2-Dpt_Low, 0);

	fout->cd();
	gr_Best_PNP->Write("",TObject::kOverwrite);




	return 0;
}


int smearTree(TFile *f_mcp_new,TFile *f_mc,int nSclErrF, Float_t *ScaleErrFactorArr,  int nSmear=10, double DptLow=2,double DptHigh=40){

	TRandom3 Rdm;
	Rdm.SetSeed(0);

	TTree *t_mcp=(TTree*)f_mc->Get("t_fit");
	TH1D *hGen_pt=(TH1D*)f_mc->Get("hGen_pt");
	f_mcp_new->cd();
	hGen_pt->Write("",TObject::kOverwrite);

	Float_t Ddls_SclErrF[nSclErrF];
	Float_t Ddl_SclErrF[nSclErrF];
	Float_t DdlErr_SclErrF[nSclErrF];


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
	// t_mcp->SetBranchAddress("DtktkResmass",&DtktkResmass);
	t_mcp->SetBranchAddress("Dalpha",&Dalpha);
	t_mcp->SetBranchAddress("Dchi2cl",&Dchi2cl);
	t_mcp->SetBranchAddress("TotalWeight",&TotalWeight);
	t_mcp->SetBranchAddress("Dalpha_BS_2D",&Dalpha_BS_2D);
	t_mcp->SetBranchAddress("DlxyBS",&DlxyBS);
	t_mcp->SetBranchAddress("DlxyBSErr",&DlxyBSErr);
	t_mcp->SetBranchAddress("DlxyBSs",&DlxyBSs);

	Float_t SMrWt=(float)1.0/nSmear;

	f_mcp_new->cd();
	TTree *t_new=new TTree("t_fit","t_fit");
	t_new->Branch("Dmass",&Dmass);
	t_new->Branch("Dpt",&Dpt);
	t_new->Branch("Ddls",&Ddls);
	t_new->Branch("Ddl",&Ddl);
	t_new->Branch("DdlErr",&DdlErr);
	t_new->Branch("Dalpha",&Dalpha);
	t_new->Branch("Dchi2cl",&Dchi2cl);
	// t_new->Branch("DtktkResmass",&DtktkResmass);
	t_new->Branch("TotalWeight",&TotalWeight);
	t_new->Branch("SMrWt",&SMrWt);
	for(int i=0; i<nSclErrF; i++){
		t_new->Branch(Form("Ddls_Scl%.0fem3",ScaleErrFactorArr[i]*1e3),&Ddls_SclErrF[i]);
		t_new->Branch(Form("Ddl_Scl%.0fem3",ScaleErrFactorArr[i]*1e3),&Ddl_SclErrF[i]);
		t_new->Branch(Form("DdlErr_Scl%.0fem3",ScaleErrFactorArr[i]*1e3),&DdlErr_SclErrF[i]);
	}


	Long64_t nentries=t_mcp->GetEntries();
	double smF=0.1;
	double DdlErrTemp=0;

	for(Long64_t i=0;i<nentries; i++)
		// for(Long64_t i=0;i<500; i++)
	{
		if(i%200000==0) {cout<<setw(10)<<i<<" / "<<nentries<<endl;}
		t_mcp->GetEntry(i);

		if(Dpt<DptLow || Dpt>DptHigh) {
			continue;
		}


		for(int j=0; j<nSmear;j++){
			/*
				 for(int k=0;k<nSmrRelF;k++){
				 DdlErr_SmrRelF[k]=Rdm.Gaus(DdlErr,DdlErr*SmearRelFactorArr[k]);
				 while(DdlErr_SmrRelF[k]<0){
				 DdlErr_SmrRelF[k]=Rdm.Gaus(DdlErr,DdlErr*SmearRelFactorArr[k]);
				 }
				 Ddls_SmrRelF[k]=Ddl / DdlErr_SmrRelF[k];
			// Ddls_SmrRelF[k]=Ddl / Rdm.Gaus(DdlErr,DdlErr*SmearRelFactorArr[k]);
			// while(Ddls_SmrRelF[k]<0){
			// Ddls_SmrRelF[k]=Ddl / Rdm.Gaus(DdlErr,DdlErr*SmearRelFactorArr[k]);
			// }
			}
			for(int k=0;k<nSmrAbsF;k++){
			DdlErr_SmrAbsF[k]=Rdm.Gaus(DdlErr,SmearAbsFactorArr[k]);
			while(DdlErr_SmrAbsF[k]<0){
			DdlErr_SmrAbsF[k]=Rdm.Gaus(DdlErr,SmearAbsFactorArr[k]);
			}
			Ddls_SmrAbsF[k]=Ddl / DdlErr_SmrAbsF[k];
			// Ddls_SmrAbsF[k]=Ddl / Rdm.Gaus(DdlErr,DdlErr*SmearAbsFactorArr[k]);
			// while(Ddls_SmrAbsF[k]<0){
			// Ddls_SmrAbsF[k]=Ddl / Rdm.Gaus(DdlErr,DdlErr*SmearAbsFactorArr[k]);
			// }
			}
			*/
			for(int k=0;k<nSclErrF;k++){
				DdlErr_SclErrF[k]=DdlErr*ScaleErrFactorArr[k];
				Ddls_SclErrF[k]=Ddl / DdlErr_SclErrF[k];
				Ddl_SclErrF[k]=Ddl;
				if(ScaleErrFactorArr[k]>1){
					Ddls_SclErrF[k]=Rdm.Gaus(Ddls_SclErrF[k],sqrt(1.0*1.0-1.0/ScaleErrFactorArr[k]*1.0/ScaleErrFactorArr[k]*1.0));
					Ddl_SclErrF[k]=Rdm.Gaus(Ddl_SclErrF[k],DdlErr_SclErrF[k]*sqrt(1.0*1.0-1.0/ScaleErrFactorArr[k]*1.0/ScaleErrFactorArr[k]*1.0) );
				}
			}
			t_new->Fill();
		}

	}

	f_mcp_new->cd();
	t_new->Write("",TObject::kOverwrite);

	return 1;


}
int main(int argc, char*argv[]){


// int DdlErrScaleScan2(int isPbPb=0, TString var_compare="DdlErr", TString var_cut="Dpt" , double var_cutLow=20, double var_cutHigh=40, double bins_var_Low=0, double bins_var_High=1, double bins_var_DrawLow=0, double bins_var_DrawHigh=0.1 , int nRebin=1, TString MC_DefaultWeight="TotalWeight",int doSmear=0, int nSmear=5, int doCreateSMTree=1){
  if(argc==6){
    DdlErrScaleScan(atoi(argv[1]), argv[2], argv[3], atof(argv[4]), atof(argv[5]) );
  }else if(argc==10){
    DdlErrScaleScan(atoi(argv[1]), argv[2], argv[3], atof(argv[4]), atof(argv[5]) , atof(argv[6]) , atof(argv[7]) , atof(argv[8]) , atof(argv[9]));
  }else{
    DdlErrScaleScan();
    cout<<"wrong number of input parameters , need 6\n int isPbPb=0, TString var_compare=DdxyzErr, TString var_cut=Dtrk1Pt , double var_cutLow=0.75, double var_cutHigh=1.25"<<endl;
    return 1;
  }
  //int Fit_sideband(int isPbPb=0, TString var_compare="DdxyzErr", TString var_cut="Dtrk1Pt" , double var_cutLow=0.75, double var_cutHigh=1.25){


  return 0;
}



