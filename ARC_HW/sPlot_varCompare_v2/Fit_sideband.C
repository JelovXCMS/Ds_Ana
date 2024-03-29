#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"

// #include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
// #include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
// #include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"

#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"


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

	double textposx=0.2;
	double textposy=0.77;

  double shiftY=0.0;
  double shiftX=0.3;
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

int Fit_sideband(int isPbPb=0, TString var_compare="Ddls", TString var_cut="Dpt" , double var_cutLow=20, double var_cutHigh=40, double bins_var_Low=0, double bins_var_High=200, double bins_var_DrawLow=0, double bins_var_DrawHigh=30 , int nRebin=1, TString MC_DefaultWeight="TotalWeight"){

	// InitStyle();

	TString GenWt="";
	if(MC_DefaultWeight=="D0Weight"){
		GenWt="_D0Weight";
	}

	gSystem->Exec("mkdir -p fitout");
	gSystem->Exec("mkdir -p plots/fit");
	gSystem->Exec(Form("mkdir -p plots/%s",var_compare.Data()));

	gStyle->SetOptStat(0);
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

	if(var_compare=="Dchi2cl"){
		Dchi2cl_cut=0.02;
	}

	double *bins_var=bins_DdxyzErr_pp;
	int nbin_var=nbin_DdxyzErr_pp;
	// double bins_var_Low=0;
	// double bins_var_High=1;
	nbin_var=1000;
	nbin_var=nbin_var/nRebin;


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

		if(var_compare=="Dchi2cl"){
			Dchi2cl_cut=0.05;
		}

	}

	TString DataCuts=Form("Dpt>%f && Dpt<%f && Dalpha < %f && Dchi2cl > %f && Ddls >%f && %s > %f && %s <%f", Dpt_Low, Dpt_High, Dalpha_cut, Dchi2cl_cut, Ddls_cut, var_cut.Data(), var_cutLow,var_cut.Data(), var_cutHigh);
	// TString DataCuts=Form("Dalpha < %f && Dchi2cl > %f && Ddls >%f && %s > %f && %s <%f", Dalpha_cut, Dchi2cl_cut, Ddls_cut, var_cut.Data(), var_cutLow,var_cut.Data(), var_cutHigh);

  TFile *f_data=TFile::Open(dataName.Data());
  TFile *f_mc_Prompt=TFile::Open(mcName_Prompt.Data());
  TFile *f_mc_NonPrompt=TFile::Open(mcName_NonPrompt.Data());
  TFile *f_mc_NonPrompt_f0=TFile::Open(mcName_NonPrompt_f0.Data());

  TTree *t_Ds_MCPrompt=(TTree*)f_mc_Prompt->Get(Form("t_fit"));
  TTree *t_Ds_MCNonPrompt=(TTree*)f_mc_NonPrompt->Get(Form("t_fit"));
  TTree *t_Ds_Data=(TTree*)f_data->Get(Form("t_fit"));

	gSystem->Exec("mkdir -p FitSideBand_out");

	TFile *f_out=TFile::Open(Form("./fitout/%s_%s_%s%.0fto%.0f.root",Str_PbPb.Data(), var_compare.Data(), var_cut.Data(), var_cutLow*100, var_cutHigh*100) ,"recreate");
	
	// import mc & fit
  RooRealVar Dmass("Dmass","Dmass",DsDataFitRangeLow,DsDataFitRangeHigh);
  // RooRealVar Ddca("Ddca","Ddca",0,0.1); // temp
  // RooRealVar Ddls("Ddls","Ddls",0,200); // temp
  RooRealVar TotalWeight("TotalWeight","TotalWeight",0,1e15);
	RooRealVar Dpt("Dpt","Dpt",0);
	RooRealVar Dalpha("Dalpha","Dalpha",0);
	RooRealVar Dchi2cl ("Dchi2cl","Dchi2cl",0); 
	RooRealVar Ddls("Ddls","Ddls",0);
	RooRealVar Dtrk1Pt ("Dtrk1Pt","Dtrk1Pt",0); 
	RooRealVar Dtrk2Pt ("Dtrk2Pt","Dtrk2Pt",0); 
	RooRealVar Dtrk3Pt ("Dtrk3Pt","Dtrk3Pt",0); 

  RooDataSet RooDS_MCAll("RooDS_MCAll","RooDS_MCAll",RooArgSet(Dmass,TotalWeight,Dpt,Dalpha,Dchi2cl,Ddls,Dtrk1Pt,Dtrk2Pt,Dtrk3Pt),WeightVar(TotalWeight),Import(*t_Ds_MCPrompt));
  // RooDataSet RooDS_MC_noW("RooDS_MC_noW","RooDS_MC_noW",RooArgSet(Dmass,TotalWeight,Dpt,Dalpha,Dchi2cl,Ddls,Dtrk1Pt,Dtrk2Pt,Dtrk3Pt),Import(*t_Ds_MCPrompt));
	

	// RooDataSet RooDS_MC_cut=*(RooDataSet *)RooDS_MCAll.reduce(RooArgSet(Dmass,TotalWeight),DataCuts.Data() ); // Weight is inherited from its mother
	RooDataSet RooDS_MC_cut=*(RooDataSet *)RooDS_MCAll.reduce(RooArgSet(Dmass),DataCuts.Data() ); // Weight is inherited from its mother
	// RooDS_MC_cut.Print("v");

	// return 1;

	// RooDataSet *RooDS_MC1_noW=(RooDataSet *)RooDS_MC_noW.reduce(RooArgSet(Dmass,TotalWeight),DataCuts.Data() );
	// RooDS_MC1->setWeightVar(TotalWeight);  // decrpreacted 

  RooRealVar DsMassMean_MC("DsMassMean_MC","DsMassMean_MC",1.9690,1.965,1.971);
  RooRealVar DsWidth1_MC("DsWidth1_MC","DsWidth1_MC",0.011,0.008,0.02);
  RooRealVar DsWidth2_MC("DsWidth2_MC","DsWidth2_MC",0.0060,0.001,0.01);
  RooRealVar DsGaus1Fr_MC("DsGaus1Fr_MC","DsGaus1Fr_MC",0.65,0.3,1);

  RooGaussian Gaus1_MC("Gaus1_MC","gauss(Dmass,DsMassMean_MC,DsWidth1_MC)",Dmass,DsMassMean_MC,DsWidth1_MC);
  RooGaussian Gaus2_MC("Gaus2_MC","gauss(Dmass,DsMassMean_MC,DsWidth2_MC)",Dmass,DsMassMean_MC,DsWidth2_MC);
  RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(Gaus1_MC,Gaus2_MC),DsGaus1Fr_MC);

  SigPdf_MC.fitTo(RooDS_MC_cut,NumCPU(10));
  SigPdf_MC.fitTo(RooDS_MC_cut,NumCPU(10));

	RooPlot *massframe_MC = new RooPlot("massframe_MC","Dmass",Dmass,DsDataFitRangeLow,DsDataFitRangeHigh,nbin_DmassDraw) ;
	// RooDS_MCAll.plotOn(massframe_MC);
	// SigPdf_MC.plotOn(massframe_MC,LineColor(2));
	RooDS_MC_cut.plotOn(massframe_MC,LineColor(1));
	SigPdf_MC.plotOn(massframe_MC,LineColor(2));
	// RooDS_MC1_noW->plotOn(massframe_MC,LineColor(2));

  RooDS_MC_cut.Print("v");
	// SigPdf_MC.Print("v");
	TCanvas *c_MC= new TCanvas("c_MC","c_MC");
	c_MC->cd();
	massframe_MC->GetXaxis()->SetTitle("m_{KK#pi} (GeV/c^{2})");
	massframe_MC->Draw();


	cout<<"DataCuts : "<<DataCuts<<endl;
	cout<<"RooDS_MCAll weight? "<<RooDS_MCAll.isWeighted()<<endl;
	cout<<"RooDS_MCcut weight? "<<RooDS_MC_cut.isWeighted()<<endl;

	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Ds MC %.0f<p_{T}<%.0f",Str_PbPb.Data(),Dpt_Low,Dpt_High)); shiftY-=oneshift;
	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Var : %s %.2f to %.2f ", var_cut.Data(), var_cutLow, var_cutHigh)); shiftY-=oneshift;

	c_MC->SaveAs(Form("./plots/fit/%s_MC_%s%.0fto%.0f.png",Str_PbPb.Data(),var_cut.Data(),var_cutLow*100,var_cutHigh*100));

	// end MC fit part

	// return 1;


	// import data
  //RooDataSet RooDSAll("RooDSAll","RooDSAll",RooArgSet(Dmass,Ddls),Import(*t_DsMassData));
  //RooDataSet *RooDS=(RooDataSet*)RooDSAll.reduce(RooArgSet(Dmass));
	
	// RooDataSet *RooDS= new RooDataSet("RooDS","RooDS",RooArgSet(Dmass),Import(*t_DsMassData));

  RooDataSet RooDS_DataAll("RooDS_DataAll","RooDS_DataAll",RooArgSet(Dmass,Dpt,Dalpha,Dchi2cl,Ddls,Dtrk1Pt,Dtrk2Pt,Dtrk3Pt),Import(*t_Ds_Data));

	RooDataSet RooDS_Data_cut=*(RooDataSet *)RooDS_DataAll.reduce(RooArgSet(Dmass),DataCuts.Data() ); // Weight is inherited from its mother


  double DsMassMeanV=DsMassMean_MC.getValV();
  double DsWidth1V=DsWidth1_MC.getValV();
  double DsWidth2V=DsWidth2_MC.getValV();
  double DsGaus1FrV=DsGaus1Fr_MC.getValV();

  RooRealVar DsMassMean("DsMassMean","DsMassMean",DsMassMeanV,1.965,1.971);
  RooRealVar DsWidth1("DsWidth1","DsWidth1",DsWidth1V,0.000001,0.1);
  RooRealVar DsWidth2("DsWidth2","DsWidth2",DsWidth2V,0.000001,0.1);
  RooRealVar DsGaus1Fr("DsGaus1Fr","DsGaus1Fr",DsGaus1FrV,0.0,1);
  RooRealVar DsFloatWidth("DsFloatWidth","DsFloatWidth",0,-0.6,0.6);
  // RooRealVar DsFloatWidth("DsFloatWidth","DsFloatWidth",DsFloatWidthV,DsFloatWidthV,DsFloatWidthV);

  DsWidth1.setConstant(kTRUE);
  DsWidth2.setConstant(kTRUE);
  DsGaus1Fr.setConstant(kTRUE);
  RooFormulaVar scale_width1("scale width1","scaled width1","DsWidth1*(1+DsFloatWidth)",RooArgSet(DsWidth1,DsFloatWidth));
  RooFormulaVar scale_width2("scale width2","scaled width2","DsWidth2*(1+DsFloatWidth)",RooArgSet(DsWidth2,DsFloatWidth));
  RooGaussian Gaus1("Gaus1","gauss(Dmass,DsMassMean,scale_width1)",Dmass,DsMassMean,scale_width1);
  RooGaussian Gaus2("Gaus2","gauss(Dmass,DsMassMean,scale_width2)",Dmass,DsMassMean,scale_width2);
  RooAddPdf SigPdf("SigPdf","SigPdf",RooArgList(Gaus1,Gaus2),DsGaus1Fr);
  RooRealVar Cheb1("Cheb1","Cheb1",0,-1,1); // no input is better than with input
  RooRealVar Cheb2("Cheb2","Cheb2",0,-1,1);
  RooRealVar Cheb3("Cheb3","Cheb3",0,-1,1);
  RooChebychev *BkgPdf;
  BkgPdf=new RooChebychev("BkgPdf","BkgPdf",Dmass,RooArgList(Cheb1,Cheb2));
  // RooRealVar NumSig("NumSig","Number of Signal",1000,-1e4,4e8);
  RooRealVar NumSig("NumSig","Number of Signal",1000,0,4e8);
  RooRealVar NumBkg("NumBkg","Number of Background",10000,0,1e10);

  RooAddPdf *RooDsMixPdf= new RooAddPdf("RooDsMixPdf","RooDsMixPdf",RooArgList(SigPdf,*BkgPdf),RooArgList(NumSig,NumBkg));

	if(isPbPb==0 ){

		Cheb2.setConstant(kTRUE);
	}

  RooDsMixPdf->fitTo(RooDS_Data_cut,Extended(kTRUE),NumCPU(20));
  RooDsMixPdf->fitTo(RooDS_Data_cut,Extended(kTRUE),NumCPU(20));
  RooDsMixPdf->fitTo(RooDS_Data_cut,Extended(kTRUE),NumCPU(20));

	TCanvas *c_data=new TCanvas("c_data","c_data");
	c_data->cd();

  RooPlot* massframe=new RooPlot("massframe","Dmass",Dmass,DsDataFitRangeLow,DsDataFitRangeHigh,nbin_DmassDraw);
  // RooPlot* massframe=new RooPlot("massframe","Dmass",Dmass,DsDataFitRangeLow,DsDataFitRangeHigh,nbin_DmassDraw);
  RooDS_Data_cut.plotOn(massframe,DataError(RooAbsData::SumW2));
  RooDsMixPdf->plotOn(massframe,LineColor(2));
	massframe->GetXaxis()->SetTitle("m_{KK#pi} (GeV/c^{2})");
  massframe->Draw();

	shiftY=0;
	// tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Ds Data %.0f<p_{T}<%.0f",Str_PbPb.Data(),Dpt_Low,Dpt_High)); shiftY-=oneshift;
	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Ds Data ",Str_PbPb.Data())); shiftY-=oneshift;
	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Var : %s %.2f to %.2f ", var_cut.Data(), var_cutLow, var_cutHigh)); shiftY-=oneshift;

	c_data->SaveAs(Form("./plots/fit/%s_Data_%s%.0fto%.0f.png",Str_PbPb.Data(),var_cut.Data(),var_cutLow*100,var_cutHigh*100));

	// return 1;

  cout<<" \n\nSigPdf = "<< SigPdf.getVal()<<endl;

	RooArgSet nset(Dmass);
	cout<<"SigPdf[Dmass] = "<< SigPdf.getVal(&nset)<<endl;

	RooAbsReal* igx = SigPdf.createIntegral(Dmass) ;
	cout << "gx_Int[x] = " << igx->getVal() << endl ;

	Dmass.setRange("signal",Dmass2SigLow,Dmass2SigHigh);

	RooAbsReal* igx_sig = SigPdf.createIntegral(Dmass,NormSet(Dmass),Range("signal")) ;
	cout << "Sig_Int[x|signal]_Norm[x] = " << igx_sig->getVal() << endl ;	

	RooAbsReal* igx_bkg = BkgPdf->createIntegral(Dmass,NormSet(Dmass),Range("signal")) ;
	cout << "Bkg_Int[x|signal]_Norm[x] = " << igx_bkg->getVal() << endl ;	


	cout<<" \n------------------ \n"<<endl;

	double N_Sig_2sig=NumSig.getValV()* igx_sig->getVal();
	double N_Bkg_2sig=NumBkg.getValV()* igx_bkg->getVal();

	cout<<"N_Sig_2sig = "<<N_Sig_2sig<<endl;
	cout<<"N_Bkg_2sig = "<<N_Bkg_2sig<<endl;

// end fitting
	// return 1;


	// TString DataCuts=Form("Dpt>%f && Dpt<%f && Dalpha < %f && Dchi2cl > %f && Ddls >%f && %s > %f && %s <%f", Dpt_Low, Dpt_High, Dalpha_cut, Dchi2cl_cut, Ddls_cut, var_cut.Data(), var_cutLow,var_cut.Data(), var_cutHigh);
	// start to plot var by sideband subtraction

	

	// TH1D *h_var_sideband= new TH1D(Form("h_%s_sideband",var_compare.Data()),Form("h_%s_sideband",var_compare.Data()),nbin_var, bins_var); h_var_sideband->Sumw2();

	TH1D *h_Data2sig_SigFrac=new TH1D("h_Data2sig_SigFrac","h_Data2sig_SigFrac",1,0,1); h_Data2sig_SigFrac->Sumw2();
	h_Data2sig_SigFrac->SetBinContent(1,N_Sig_2sig/(N_Sig_2sig+N_Bkg_2sig));


	TH1D *h_var_sideband= new TH1D(Form("h_%s_sideband",var_compare.Data()),Form("h_%s_sideband",var_compare.Data()),nbin_var, bins_var_Low,bins_var_High); h_var_sideband->Sumw2();
	t_Ds_Data->Project(Form("h_%s_sideband",var_compare.Data()),var_compare.Data(),Form("( (Dmass>%f && Dmass <%f ) || ( Dmass >%f && Dmass<%f) ) && %s",DmassSideBand1Low, DmassSideBand1High, DmassSideBand2Low, DmassSideBand2High, DataCuts.Data()));

	// TH1D *h_var_Data2sig= new TH1D(Form("h_%s_Data2sig",var_compare.Data()),Form("h_%s_Data2sig",var_compare.Data()),nbin_var, bins_var); h_var_Data2sig->Sumw2();
	TH1D *h_var_Data2sig= new TH1D(Form("h_%s_Data2sig",var_compare.Data()),Form("h_%s_Data2sig",var_compare.Data()),nbin_var, bins_var_Low,bins_var_High); h_var_Data2sig->Sumw2();
	t_Ds_Data->Project(Form("h_%s_Data2sig",var_compare.Data()),var_compare.Data(),Form("(Dmass>%f && Dmass<%f) && %s",Dmass2SigLow,Dmass2SigHigh, DataCuts.Data()));

	TH1D *h_var_Data2sigMix=(TH1D*)h_var_Data2sig->Clone(Form("h_%s_Data2sigMix",var_compare.Data()));
	
	h_var_sideband->Scale(N_Bkg_2sig/h_var_sideband->Integral());
	
	h_var_Data2sig->Add(h_var_sideband,-1);
	h_var_Data2sig->Scale(1/h_var_Data2sig->Integral());

	TCanvas *c_var_Data=new TCanvas("c_varData","c_varData",800,800);
	c_var_Data->cd();
	
	// h_var_Data2sig->GetXaxis()->SetRangeUser(0,0.1);
	// h_var_Data2sig->SetMaximum();
	h_var_Data2sig->GetXaxis()->SetRangeUser(bins_var_DrawLow,bins_var_DrawHigh); // normalized before setRange
	h_var_Data2sig->SetTitle("");
	h_var_Data2sig->GetXaxis()->SetTitle(var_compare.Data());
	h_var_Data2sig->Draw();
	h_var_Data2sig->SetLineColor(1);

	h_var_sideband->Scale(1/h_var_sideband->Integral());
	// h_var_sideband->GetXaxis()->SetRangeUser(0,0.1);
	h_var_sideband->SetLineColor(2);
	h_var_sideband->Draw("same");
	// h_var_sideband->Draw("same");
	h_var_Data2sigMix->Scale(1/h_var_Data2sigMix->Integral());
	h_var_Data2sigMix->SetLineColor(4);
	h_var_Data2sigMix->Draw("same");
//	shiftY=0;
//	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form(""))




	TLegend *le_var_Data=new TLegend(0.5,0.5,0.85,0.85);
	le_var_Data->SetBorderSize(0);
	// le_var_Data->AddEntry((TObject*)0,Form("%s Ds %s %.0f<p_{T}<%.0f",Str_PbPb.Data(), var_compare.Data(),Dpt_Low,Dpt_High),"");
	le_var_Data->AddEntry((TObject*)0,Form("%s Ds %s ",Str_PbPb.Data(), var_compare.Data()),"");
	le_var_Data->AddEntry((TObject*)0,Form("Var : %s %.2f to %.2f ", var_cut.Data(), var_cutLow, var_cutHigh),"");
	le_var_Data->AddEntry((TObject*)0,Form("Sig frac = %.2f", N_Sig_2sig/(N_Sig_2sig+N_Bkg_2sig) ),"");
	le_var_Data->AddEntry(h_var_Data2sig,"Data Sig","pl");
	le_var_Data->AddEntry(h_var_sideband,"Data sideband","pl");
	le_var_Data->AddEntry(h_var_Data2sigMix,"Data Mix in 2sigma","pl");
	le_var_Data->Draw("same");

	// tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Ds %s %.0f<p_{T}<%.0f",Str_PbPb.Data(), var_compare.Data(),Dpt_Low,Dpt_High)); shiftY-=oneshift;
	// tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Var : %s %.2f to %.2f ", var_cut.Data(), var_cutLow, var_cutHigh)); shiftY-=oneshift;

	c_var_Data->SaveAs(Form("./plots/%s/%s_DataAndSideBand_%s%.0fto%.0f.png", var_compare.Data(),Str_PbPb.Data(), var_cut.Data(),var_cutLow*100,var_cutHigh*100));

	// return 1;

	// MC part

/*
	double *bins_var=bins_Ddls_pp;
	int nbin_var= nbin_Ddls_pp;

	if(var=="Dalpha"){
		bins_var=bins_Dalpha_pp;
		nbin_var=nbin_Dalpha_pp;
	}else if(var=="Dchi2cl"){
		bins_var=bins_Dchi2cl_pp;
		nbin_var=nbin_Dchi2cl_pp;
	}

	if(isPbPb){
		bins_var=bins_Ddls_PbPb;
		nbin_var= nbin_Ddls_PbPb;
		if(var=="Dalpha"){
		bins_var=bins_Dalpha_PbPb;
		nbin_var=nbin_Dalpha_PbPb;
		}else if(var=="Dchi2cl"){
		bins_var=bins_Dchi2cl_PbPb;
		nbin_var=nbin_Dchi2cl_PbPb;
		}
	}
*/

	// MC histogram 
	// TH1D *h_var_PromptMC =new TH1D(Form("h_%s_PromptMC",var_compare.Data()), Form("h_%s_PromptMC",var_compare.Data()), nbin_var, bins_var); h_var_PromptMC->Sumw2();
	TH1D *h_var_PromptMC =new TH1D(Form("h_%s_PromptMC",var_compare.Data()), Form("h_%s_PromptMC",var_compare.Data()), nbin_var, bins_var_Low,bins_var_High); h_var_PromptMC->Sumw2();
	t_Ds_MCPrompt->Project(Form("h_%s_PromptMC",var_compare.Data()),var_compare.Data(),Form("(Dmass>%f && Dmass<%f && %s)*%s",Dmass2SigLow, Dmass2SigHigh, DataCuts.Data(),MC_DefaultWeight.Data() ));

	h_var_PromptMC->Scale(1/h_var_PromptMC->Integral());


	// TH1D *h_var_NonPromptMC =new TH1D(Form("h_%s_NonPromptMC",var_compare.Data()), Form("h_%s_NonPromptMC",var_compare.Data()), nbin_var, bins_var); h_var_NonPromptMC->Sumw2();
	TH1D *h_var_NonPromptMC =new TH1D(Form("h_%s_NonPromptMC",var_compare.Data()), Form("h_%s_NonPromptMC",var_compare.Data()), nbin_var, bins_var_Low,bins_var_High); h_var_NonPromptMC->Sumw2();
	t_Ds_MCNonPrompt->Project(Form("h_%s_NonPromptMC",var_compare.Data()),var_compare.Data(),Form("(Dmass>%f && Dmass<%f && %s)*%s", Dmass2SigLow, Dmass2SigHigh, DataCuts.Data(),MC_DefaultWeight.Data()));

  h_var_NonPromptMC->Scale(1/h_var_NonPromptMC->Integral());	


	h_var_PromptMC->SetLineColor(2);
	h_var_NonPromptMC->Draw("same");
	h_var_NonPromptMC->SetLineColor(4);
	h_var_PromptMC->Draw("same");


	// return 1;
 
	double fr_withcut_prompt=0.85; // just for get quick result

	// extract fr_prompt from fit and BtoDs estimation
	bool use_fr_fromBtoD=true;

	if(use_fr_fromBtoD){
		double N_yield=NumSig.getValV();
		double LumiNevt=LumiSum;
		if(isPbPb){LumiNevt=NevtPbPb3;}


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

	// nee nonprompt efficiency here
		

//		TH1D *h_eff_nonprompt_phi=(TH1D*)f_mc_NonPrompt->Get(Form("h_RecoEff_for%s",var_compare.Data()));
//		TH1D *h_eff_nonprompt_f0=(TH1D*)f_mc_NonPrompt_f0->Get(Form("h_RecoEff_for%s",var_compare.Data()));

		TH1D *h_eff_nonprompt_phi=new TH1D("h_eff_nonprompt_phi","h_eff_nonprompt_phi",40,0,40);	h_eff_nonprompt_phi->Sumw2();
		t_Ds_MCNonPrompt->Project("h_eff_nonprompt_phi","Dpt",(TCut)Form("%s*(%s)",MC_DefaultWeight.Data(),DataCuts.Data()));

		double NonPrompt_phi_reco=h_eff_nonprompt_phi->Integral(h_eff_nonprompt_phi->FindBin(Dpt_Low+0.00001), h_eff_nonprompt_phi->FindBin(Dpt_High-0.00001));
		TH1D *hGen_pt_nonprompt_phi=(TH1D*)f_mc_NonPrompt->Get(Form("hGen_pt%s",GenWt.Data() ) );
		double NonPrompt_phi_Gen=hGen_pt_nonprompt_phi->Integral(hGen_pt_nonprompt_phi->FindBin(Dpt_Low+0.00001),hGen_pt_nonprompt_phi->FindBin(Dpt_High-0.0001) );
		

//		double NonPrompt_phi_Eff=h_eff_nonprompt_phi->GetBinContent(1);
//		double NonPrompt_f0_Eff=h_eff_nonprompt_f0->GetBinContent(1);
		double NonPrompt_phi_Eff=NonPrompt_phi_reco/NonPrompt_phi_Gen;



		cout<<"NonPrompt_phi_Eff = "<<NonPrompt_phi_Eff<<endl;
//		cout<<"NonPrompt_f0_Eff = "<<NonPrompt_f0_Eff<<endl;
		double phiFr=0.95 ;
		if(isPbPb){
			phiFr=0.896;
		}
	
		double N_NonPrompt_yield=CS_Integral*(2*LumiNevt)*( BRphi*NonPrompt_phi_Eff)/phiFr;
		// double N_NonPrompt_yield=CS_Integral*(2*LumiNevt)*( BRphi*NonPrompt_phi_Eff+ BRf0*NonPrompt_f0_Eff );



//		cout<<"LumiNevt = "<<LumiNevt<<endl;

		fr_withcut_prompt=1-(N_NonPrompt_yield/N_yield);

		cout<<"fr_withcut_prompt = "<<fr_withcut_prompt<<" , N_NonPrompt_yield = "<<N_NonPrompt_yield<<" , N_yield = "<<N_yield<<endl;


	} // end fr_prompt calculation



	// TH1D *h_var_MixMC= new TH1D("h_var_MixMC","h_var_MixMC",nbin_var, bins_var); h_var_MixMC->Sumw2();
	TH1D *h_var_MixMC= new TH1D(Form("h_%s_MixMC",var_compare.Data()),"h_var_MixMC",nbin_var, bins_var_Low,bins_var_High); h_var_MixMC->Sumw2();

	h_var_MixMC->Add(h_var_PromptMC,h_var_NonPromptMC,fr_withcut_prompt,1-fr_withcut_prompt);

	h_var_MixMC->Draw("same");

	cout<<"start drawing"<<endl;


	TH1D *h_fr_PromptMC=new TH1D("h_fr_PromptMC","",1,0,1);
	h_fr_PromptMC->SetBinContent(1,fr_withcut_prompt); 

	// save histogram before setRange

	f_out->cd();
	h_var_sideband->Write();
	h_var_Data2sig->Write();
	h_var_Data2sigMix->Write();
	h_var_PromptMC->Write();
	h_var_NonPromptMC->Write();
	h_var_MixMC->Write();

	h_Data2sig_SigFrac->Write();
	h_fr_PromptMC->Write();

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

	int rebinN=2;
	if(var_compare=="Dchi2cl"){
		rebinN=20;
	}
	h_var_Data2sig->Rebin(rebinN);
	h_var_PromptMC->Rebin(rebinN);
	h_var_NonPromptMC->Rebin(rebinN);
	h_var_MixMC->Rebin(rebinN);

	h_var_Data2sig->GetXaxis()->SetRangeUser(bins_var_DrawLow,bins_var_DrawHigh); // normalized before setRange
	h_var_Data2sig->GetXaxis()->SetTitle(var_compare.Data());
	h_var_Data2sig->SetTitle("");
	h_var_Data2sig->SetLineColor(1);
	h_var_Data2sig->SetMarkerColor(1);
	h_var_Data2sig->SetMarkerStyle(22);
	h_var_Data2sig->Draw();
	h_var_PromptMC->SetLineColor(2);
	h_var_PromptMC->SetMarkerColor(2);
	h_var_PromptMC->SetMarkerStyle(26);
	h_var_PromptMC->Draw("same");
	h_var_NonPromptMC->SetLineColor(4);
	h_var_NonPromptMC->SetMarkerColor(4);
	h_var_NonPromptMC->SetMarkerStyle(26);
	h_var_NonPromptMC->Draw("same");
	h_var_MixMC->SetLineColor(1);
	h_var_MixMC->SetMarkerColor(1);
	h_var_MixMC->SetMarkerStyle(26);
	h_var_MixMC->Draw("same");


	TLegend *le_var = new TLegend(0.5,0.28,0.85,0.65);
	le_var->SetBorderSize(0);
	le_var->AddEntry(h_var_Data2sig,"Data","pl");
	le_var->AddEntry(h_var_PromptMC,"MC Prompt D_{S}","pl");
	le_var->AddEntry(h_var_NonPromptMC,"MC NonPrompt D_{S}","pl");
	le_var->AddEntry(h_var_MixMC,"MC Mix D_{S}","pl");
	le_var->Draw("same");

	shiftY=0.03;
	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Ds %s ",Str_PbPb.Data(), var_compare.Data())); shiftY-=oneshift;
	// tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Ds %s %.0f<p_{T}<%.0f",Str_PbPb.Data(), var_compare.Data(),Dpt_Low,Dpt_High)); shiftY-=oneshift;
	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Var : %s %.2f to %.2f ", var_cut.Data(), var_cutLow, var_cutHigh)); shiftY-=oneshift;
	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Prompt fraction : %.2f ", fr_withcut_prompt)); shiftY-=oneshift;

	pad2->cd();

	TH1D *h_ratio=(TH1D*)h_var_Data2sig->Clone("h_ratio");
	h_ratio->Divide(h_var_MixMC);
	h_ratio->GetYaxis()->SetTitle("Data/MC");
	h_ratio->GetYaxis()->CenterTitle();
	h_ratio->GetYaxis()->SetLabelSize(0.05);
	h_ratio->GetYaxis()->SetTitleSize(0.08);
	h_ratio->GetYaxis()->SetTitleOffset(0.38);
	h_ratio->GetXaxis()->SetLabelSize(0.12);
	h_ratio->GetXaxis()->SetTitleSize(0.12);
	h_ratio->SetMaximum(2);
	h_ratio->SetMinimum(0);
	h_ratio->Draw();



	c_var_compare->SaveAs(Form("./plots/%s/%s_%s%.0fto%.0f.png", var_compare.Data(),Str_PbPb.Data(), var_cut.Data(),var_cutLow*100,var_cutHigh*100));


	TH1D *h_var_ratioSmooth=(TH1D*)h_var_Data2sig->Clone(Form("h_%s_ratioSmooth",var_compare.Data()));
	// int RebinN=2;
	int SmoothN=2;
	// h_var_ratioSmooth->Rebin(RebinN);
	// h_var_MixMC->Rebin(RebinN);
	h_var_ratioSmooth->Divide(h_var_MixMC);
  h_var_ratioSmooth->Smooth(SmoothN);

	h_var_ratioSmooth->Write();


	TCanvas *c_weight=new TCanvas("c_weight","c_weight");
	c_weight->cd();
	h_var_ratioSmooth->Draw();


	cout<<"Data mean = "<<h_var_Data2sig->GetMean()<<" , RMS = "<<h_var_Data2sig->GetRMS()<<" , StdDev = "<<h_var_Data2sig->GetStdDev()<<endl;
	cout<<"MC mean = "<<h_var_MixMC->GetMean()<<" , RMS = "<<h_var_MixMC->GetRMS()<<" , StdDev = "<<h_var_MixMC->GetStdDev()<<endl;

	f_out->Close();

/*
	c_Ddls_sideband->Divide(2,2);
	c_Ddls_sideband->cd(1);
	h_var_compare.Data2sig->Draw();
	c_Ddls_sideband->cd(2);
//	h_var_sideband->Draw();
	h_var_PromptMC->Draw();
	c_Ddls_sideband->cd(3);
	h_var_NonPromptMC->Draw();
	c_Ddls_sideband->cd(4);

	h_var_MixMC->SetMaximum(0.3);
	h_var_MixMC->SetMinimum(0.0);
	h_var_MixMC->Draw();
	h_var_compare.Data2sig->SetLineColor(2);
	h_var_compare.Data2sig->Draw("same");

	c_Ddls_sideband->SaveAs(Form("./plots/%s_%s.png",Str_PbPb.Data(),var_compare.Data()));

	InitStyle();
	setTDRStyle();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TCanvas *c_var_compare = new TCanvas("c_var_compare","c_var_compare",c_wtopx,c_wtopy,c_W,c_H);
	SetCanvas(c_var_compare);
	c_var_compare->cd();

  h_var_MixMC->SetMaximum(0.3);
  h_var_MixMC->SetMinimum(0.0);
	h_var_MixMC->SetTitle("");
	h_var_MixMC->GetXaxis()->SetTitle(var_compare.Data());
	h_var_MixMC->GetYaxis()->SetTitle("Normalized Distribution");
  h_var_MixMC->Draw();
  h_var_compare.Data2sig->SetLineColor(2);
	h_var_compare.Data2sig->SetMarkerColor(2);
  h_var_compare.Data2sig->Draw("same");

	TLegend *le = new TLegend(0.6,0.62,0.88,0.85);
	le->SetBorderSize(0);
	le->AddEntry((TObject*)0,Form("%s",Str_PbPb.Data()),"");
	le->AddEntry(h_var_MixMC,"MC","l");
	le->AddEntry(h_var_compare.Data2sig,"data sideband sub.","l");
	le->Draw("same");

	SavePlotDirs(c_var_compare,Form("%s_%s_var_MCDataCompare",Str_PbPb.Data(),var_compare.Data()),{"Miscellaneous","var_MCDataCompare"});



	f_out->cd();
	h_var_MixMC->Write();
	h_var_compare.Data2sig->Write();


	// cumulative for all above, considering the cut is > certain value
	TH1D *h_var_MixMC_Cum= new TH1D("h_var_MixMC_Cum","h_var_MixMC_Cum",nbin_var, bins_var); h_var_MixMC_Cum->Sumw2();
	TH1D *h_var_compare.Data2sig_Cum= new TH1D("h_var_compare.Data2sig_Cum","h_var_compare.Data2sig_Cum",nbin_var, bins_var); h_var_compare.Data2sig_Cum->Sumw2();
	double cum_MixMC=0;
	double cum_Data=0; 

	for(int i=1; i<=nbin_var; i++){
		for(int j=i; j<=nbin_var; j++){
			cum_MixMC+=h_var_MixMC->GetBinContent(j);
			cum_Data+=h_var_compare.Data2sig->GetBinContent(j);			

		} // end for j
		h_var_MixMC_Cum->SetBinContent(i,cum_MixMC);		
		h_var_compare.Data2sig_Cum->SetBinContent(i,cum_Data);		

	}// end for i 


	TH1D *h_var_Ratio_Cum= new TH1D("h_var_Ratio_Cum","h_var_Ratio_Cum",nbin_var, bins_var); h_var_Ratio_Cum->Sumw2();
	h_var_Ratio_Cum->Divide(h_var_compare.Data2sig_Cum,h_var_MixMC_Cum);

	h_var_MixMC_Cum->Write();
	h_var_compare.Data2sig_Cum->Write();
	h_var_Ratio_Cum->Write();
 


	return 1;
*/
	// cout<<"textposx = "<<textposx<<endl;
	// cout<<"textposy = "<<textposy<<endl;


	return 0;

}

int main(int argc, char*argv[]){


	if(argc==6){
		Fit_sideband(atoi(argv[1]), argv[2], argv[3], atof(argv[4]), atof(argv[5]) );
	}else if(argc==10){	
		Fit_sideband(atoi(argv[1]), argv[2], argv[3], atof(argv[4]), atof(argv[5]) , atof(argv[6]) , atof(argv[7]) , atof(argv[8]) , atof(argv[9]));
	}else{
		Fit_sideband();
		cout<<"wrong number of input parameters , need 6\n int isPbPb=0, TString var_compare=DdxyzErr, TString var_cut=Dtrk1Pt , double var_cutLow=0.75, double var_cutHigh=1.25"<<endl;
		return 1;
	}
//int Fit_sideband(int isPbPb=0, TString var_compare="DdxyzErr", TString var_cut="Dtrk1Pt" , double var_cutLow=0.75, double var_cutHigh=1.25){


	return 0;
}
