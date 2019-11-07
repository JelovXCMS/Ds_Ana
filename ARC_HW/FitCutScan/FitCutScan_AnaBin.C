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

using namespace RooFit;
using namespace std;

// double textposx=0.2;
// double textposy=0.77;

double shiftY=0.0;
double shiftX=0.35;
double oneshift=0.075;

  double DsMassMCMean=1.96834;
  double DsMass3SigWidth=0.03;



int FitCutScan_AnaBin(int isPbPb=0,int useAnaBin=1, int ibin_Dpt=7, double DptLow=6, double DptHigh=40, TString var_scan="Dalpha", double scan_val=0.05, int FixShape=0, double shapeVal=3.0 , int useDefaultCut=1, int useBkgSub=0){

// double textposx=0.2;
// double textposy=0.77;

	initParameter();
	setTDRStyle();

  double tex_upperY=0.95;

  TLatex* texCmsPre = new TLatex(0.15,tex_upperY, "#scale[1.25]{CMS} Preliminary");
  texCmsPre->SetNDC();
  texCmsPre->SetTextAlign(12);
  texCmsPre->SetTextSize(0.04);
  texCmsPre->SetTextFont(42);

  TLatex* texCmsSim = new TLatex(0.15,tex_upperY, "#scale[1.25]{CMS} Simulations");
  texCmsSim->SetNDC();
  texCmsSim->SetTextAlign(12);
  texCmsSim->SetTextSize(0.04);
  texCmsSim->SetTextFont(42);	

  TLatex* texColPbPb = new TLatex(0.95,tex_upperY, "530 #mub^{-1} (5.02 TeV PbPb)");
  texColPbPb->SetNDC();
  texColPbPb->SetTextAlign(32);
  texColPbPb->SetTextSize(0.04);
  texColPbPb->SetTextFont(42);

  // TLatex* texColpp = new TLatex(0.88,tex_upperY, "pp #sqrt{s} = 5.02 TeV");
  TLatex* texColpp = new TLatex(0.95,tex_upperY, "27.4 pb^{-1} (5.02 TeV pp )");
  texColpp->SetNDC();
  texColpp->SetTextAlign(32);
  texColpp->SetTextSize(0.04);
  texColpp->SetTextFont(42);


	double DsDataFitRangeLow=1.91;
	double DsDataFitRangeHigh=2.11;
	int nbin_DmassDraw=50;

	double DsMassMeanFix  =0; 
	double DsWidth1Fix    =0; 
	double DsWidth2Fix    =0; 
	double DsWidth3Fix    =0; 
	double DsGaus1FrFix   =0; 
	double DsGaus2FrFix   =0; 
	double DsFloatWidthFix=0;
	TString s_FixShape="";

	// int FixShape=0;//

	TString s_PbPb3="pp";
	TString s_PbPb="pp";
	if(isPbPb){
		s_PbPb3="PbPb3";
		s_PbPb="PbPb";
	}

	double Dalpha_cut=0.12;
	double Dchi2cl_cut=0.1;
	double Ddls_cut=2.5;
	if(isPbPb==0 && DptHigh<=6){
		Dchi2cl_cut=0.1;
		Ddls_cut=3.5;
	}
	if(isPbPb==0 && DptLow>=6){
		Dchi2cl_cut=0.03;
		Ddls_cut=2.4;
	}
	if(isPbPb==3){
		Dchi2cl_cut=0.25;
		Ddls_cut=4.8;
	}

	double *DalphaMax_bins=DalphaMax_bins_pp;
	double *Dchi2clMin_bins=Dchi2clMin_bins_pp;
	double *DdlsMin_bins=DdlsMin_bins_pp;
	double *bins_pt=bins_pt_pp;
//	PhiMassScan_bins

	if(isPbPb==3){
		DalphaMax_bins=DalphaMax_bins_PbPb3;
		Dchi2clMin_bins=Dchi2clMin_bins_PbPb3;
		DdlsMin_bins=DdlsMin_bins_PbPb3;
		bins_pt=bins_pt_PbPb3;
	}

	if(useAnaBin){
		Dalpha_cut=DalphaMax_bins[ibin_Dpt];			
		Dchi2cl_cut=Dchi2clMin_bins[ibin_Dpt];
		Ddls_cut=DdlsMin_bins[ibin_Dpt];
		// PhiMass_cut=DtktkResmassCutWidth;
		DptLow=bins_pt[ibin_Dpt];
		DptHigh=bins_pt[ibin_Dpt+1];
	}



	// if(isPbPb==3 && DptHigh<=10){
	// Dchi2cl_cut=0.25;
	// Ddls_cut=4.8;
	// }
	// if(isPbPb==3 && DptLow>=10){
	// Dchi2cl_cut=0.08;
	// Ddls_cut=4.7;
	// }



	double DtktkResmass_cut=DtktkResmassCutWidth;  // 0.009
	double PhiMass=DtktkResmassCutMean;

	if(!useDefaultCut) // do not change cut if produce default cut fit
	{ 

		if(var_scan=="Ddls"){
			Ddls_cut=scan_val;
		}else if(var_scan=="Dchi2cl"){
			Dchi2cl_cut=scan_val;
		}else if(var_scan=="Dalpha"){
			Dalpha_cut=scan_val;
		}else if(var_scan=="PhiMass"){
			DtktkResmass_cut=scan_val;
		}else{
			cout<<"unknow var for scan , terminate withouth running "<<var_scan<<endl;
			return 2;
	
		}
	}


	if(FixShape){
		s_FixShape="FixShape";
		char buffer[100];
		fstream file;
		TString fparName=Form("./par/%s/%s_FitPar_Dpt%.0fto%.0f_%s%.0f.txt",var_scan.Data(),s_PbPb.Data(), DptLow,DptHigh,var_scan.Data(),shapeVal*100 );
		if(useAnaBin){
			fparName=Form("./par/%s_FitPar_Dpt%.0fto%.0f_Default.txt",s_PbPb.Data(), DptLow,DptHigh);
		}
		// file.open(Form("./par/%s/%s_FitPar_Dpt%.0fto%.0f_%s%.0f.txt",var_scan.Data(),s_PbPb.Data(), DptLow,DptHigh,var_scan.Data(),shapeVal*100 ));
		file.open(fparName.Data());

		TString stemp;
		double pars[10];
		int count_par=0;
		if(!file){
			cout<<"failed in open file"<<endl;
			return 2;
		}else{
			while(file>>stemp){
				// count_par++;
				if(count_par%2==1){
					pars[count_par/2]=atof(stemp);
				}	
				cout<<"read from file :"<<stemp<<endl;
				count_par++;
			}
			DsMassMeanFix   =pars[0];
			DsWidth1Fix     =pars[1];
			DsWidth2Fix     =pars[2];
			DsWidth3Fix     =pars[3];
			DsGaus1FrFix    =pars[4];
			DsGaus2FrFix    =pars[5];
			DsFloatWidthFix =pars[6];
		}
		for(int i =0; i<10; i++){
			cout<<"par i "<<i<<" , = "<<pars[i]<<endl;	
		}
	}
	// return 1;

	TLatex *tltx=new TLatex();

	gSystem->Exec(Form("mkdir -p fitout/%s",var_scan.Data()));
	gSystem->Exec(Form("mkdir -p plots/%s",var_scan.Data()));

	gStyle->SetOptStat(0);

	TString dataName=Form("/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_cutScan/%s_fitFile.root",s_PbPb3.Data());
	TString mcName_Prompt=Form("/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_cutScan/%sMC_phiPrompt_fitFile.root",s_PbPb3.Data());
	TString mcName_NonPrompt=Form("/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_cutScan/%sMC_phiNonPrompt_fitFile.root",s_PbPb3.Data());

	TFile *f_data=TFile::Open(dataName.Data());
	TFile *f_mc_Prompt=TFile::Open(mcName_Prompt.Data());
	TFile *f_mc_NonPrompt=TFile::Open(mcName_NonPrompt.Data());
	//  TFile *f_mc_NonPrompt_f0=TFile::Open(mcName_NonPrompt_f0.Data());

	TTree *t_Ds_MCPrompt=(TTree*)f_mc_Prompt->Get(Form("t_fit"));
	TTree *t_Ds_MCNonPrompt=(TTree*)f_mc_NonPrompt->Get(Form("t_fit"));
	TTree *t_Ds_Data=(TTree*)f_data->Get(Form("t_fit"));


	TString foutName=Form("./fitout/%s/%s_%sDpt%.0fto%.0f_%s%.0f.root",var_scan.Data(),s_PbPb.Data(),s_FixShape.Data(),DptLow,DptHigh,var_scan.Data(),scan_val*100);
	if(var_scan=="PhiMass"){
		foutName=Form("./fitout/%s/%s_%sDpt%.0fto%.0f_%s%.0f.root",var_scan.Data(),s_PbPb.Data(),s_FixShape.Data(),DptLow,DptHigh,var_scan.Data(),scan_val*10000);
	}
	

	if(useDefaultCut && useAnaBin){
		foutName=Form("./fitout/%s_%sDpt%.0fto%.0f_Default.root",s_PbPb.Data(),s_FixShape.Data(),DptLow,DptHigh);	
	}

	TFile *fout=TFile::Open(foutName.Data(),"recreate");

	TString DataCuts=Form("Dpt>%f && Dpt<%f && Dalpha<%f && Dchi2cl>%f && Ddls>%f && DtktkResmass>%f && DtktkResmass<%f", DptLow,DptHigh,Dalpha_cut,Dchi2cl_cut,Ddls_cut, PhiMass-DtktkResmass_cut, PhiMass+DtktkResmass_cut);

	/// read MC

	RooRealVar Dmass("Dmass","Dmass",DsDataFitRangeLow,DsDataFitRangeHigh);
	// RooRealVar Ddca("Ddca","Ddca",0,0.1); // temp
	// RooRealVar Ddls("Ddls","Ddls",0,200); // temp
	RooRealVar TotalWeight("TotalWeight","TotalWeight",0,1e15);
	RooRealVar Dpt("Dpt","Dpt",0);
	RooRealVar Dalpha("Dalpha","Dalpha",0);
	RooRealVar Dchi2cl ("Dchi2cl","Dchi2cl",0);
	RooRealVar Ddls("Ddls","Ddls",0);
	RooRealVar DtktkResmass("DtktkResmass","DtktkResmass",0);

	RooDataSet RooDS_MCAll("RooDS_MCAll","RooDS_MCAll",RooArgSet(Dmass,TotalWeight,Dpt,Dalpha,Dchi2cl,Ddls,DtktkResmass),WeightVar(TotalWeight),Import(*t_Ds_MCPrompt));
	RooDataSet RooDS_MC_cut=*(RooDataSet *)RooDS_MCAll.reduce(RooArgSet(Dmass),DataCuts.Data() ); // Weight is inherited from its mother


	RooRealVar DsMassMean_MC("DsMassMean_MC","DsMassMean_MC",1.9690,1.965,1.971);
	RooRealVar DsWidth1_MC("DsWidth1_MC","DsWidth1_MC",0.011,0.003,0.1);
	RooRealVar DsWidth2_MC("DsWidth2_MC","DsWidth2_MC",0.0060,0.0001,0.02);
	RooRealVar DsWidth3_MC("DsWidth3_MC","DsWidth3_MC",0.0060,0.0001,0.02);
	RooRealVar DsGaus1Fr_MC("DsGaus1Fr_MC","DsGaus1Fr_MC",0.65,0.01,1);
	RooRealVar DsGaus2Fr_MC("DsGaus2Fr_MC","DsGaus2Fr_MC",0.65,0.01,1);

	RooGaussian Gaus1_MC("Gaus1_MC","gauss(Dmass,DsMassMean_MC,DsWidth1_MC)",Dmass,DsMassMean_MC,DsWidth1_MC);
	RooGaussian Gaus2_MC("Gaus2_MC","gauss(Dmass,DsMassMean_MC,DsWidth2_MC)",Dmass,DsMassMean_MC,DsWidth2_MC);
	RooGaussian Gaus3_MC("Gaus3_MC","gauss(Dmass,DsMassMean_MC,DsWidth3_MC)",Dmass,DsMassMean_MC,DsWidth3_MC);
	RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(Gaus1_MC,Gaus2_MC),DsGaus1Fr_MC);
	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(Gaus1_MC,Gaus2_MC,Gaus3_MC),RooArgList(DsGaus1Fr_MC,DsGaus2Fr_MC));

	if(FixShape){
		DsMassMean_MC.setVal(DsMassMeanFix);
		DsWidth1_MC.setVal(DsWidth1Fix);
		DsWidth2_MC.setVal(DsWidth2Fix);
		DsWidth3_MC.setVal(DsWidth3Fix);
		DsGaus1Fr_MC.setVal(DsGaus1FrFix);
		DsGaus2Fr_MC.setVal(DsGaus2FrFix);

		DsWidth1_MC.setConstant(kTRUE);
		DsWidth2_MC.setConstant(kTRUE);
		DsWidth3_MC.setConstant(kTRUE);
		DsGaus1Fr_MC.setConstant(kTRUE);
		DsGaus2Fr_MC.setConstant(kTRUE);

	}



	// SigPdf_MC.fitTo(RooDS_MC_cut,Range(1.92,2.02));
	// SigPdf_MC.fitTo(RooDS_MC_cut,Range(1.91,2.03));
	// SigPdf_MC.fitTo(RooDS_MC_cut,Range(1.91,2.03));
	// SigPdf_MC.fitTo(RooDS_MC_cut,NumCPU(10));
	SigPdf_MC.fitTo(RooDS_MC_cut,NumCPU(20));
	SigPdf_MC.fitTo(RooDS_MC_cut,NumCPU(20));
	SigPdf_MC.fitTo(RooDS_MC_cut,NumCPU(20));

	// RooPlot *massframe_MC = new RooPlot("massframe_MC","Dmass",Dmass,DsDataFitRangeLow,DsDataFitRangeHigh,nbin_DmassDraw) ;
	RooPlot *massframe_MC = Dmass.frame();
	// RooDS_MCAll.plotOn(massframe_MC);
	// SigPdf_MC.plotOn(massframe_MC,LineColor(2));
	RooDS_MC_cut.plotOn(massframe_MC,LineColor(1));
	SigPdf_MC.plotOn(massframe_MC,LineColor(2));
	// RooDS_MC1_noW->plotOn(massframe_MC,LineColor(2));

	RooDS_MC_cut.Print("v");
	TCanvas *c_MC= new TCanvas("c_MC","c_MC");
	c_MC->cd();
	massframe_MC->GetXaxis()->SetTitle("m_{KK#pi} (GeV/c^{2})");
	massframe_MC->Draw();

	cout<<"DataCuts : "<<DataCuts<<endl;
	cout<<"RooDS_MCAll weight? "<<RooDS_MCAll.isWeighted()<<endl;
	cout<<"RooDS_MCcut weight? "<<RooDS_MC_cut.isWeighted()<<endl;

	// return 1;

	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Ds MC %.0f<p_{T}<%.0f",s_PbPb.Data(),DptLow,DptHigh)); shiftY-=oneshift;
	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Cut : %.2f  ", var_scan.Data(),scan_val )); shiftY-=oneshift;

	c_MC->SaveAs(Form("./plots/%s/%s_MC%s_Dpt%.0fto%.0f_%s%.0f.png",var_scan.Data(),s_PbPb.Data(),s_FixShape.Data(),DptLow,DptHigh,var_scan.Data(),scan_val*10000));


	// Data Fit

	RooDataSet RooDS_DataAll("RooDS_DataAll","RooDS_DataAll",RooArgSet(Dmass,Dpt,Dalpha,Dchi2cl,Ddls,DtktkResmass),Import(*t_Ds_Data));

	RooDataSet RooDS_Data_cut=*(RooDataSet *)RooDS_DataAll.reduce(RooArgSet(Dmass),DataCuts.Data() ); // Weight is inherited from its mother


	double DsMassMeanV=DsMassMean_MC.getValV();
	double DsWidth1V=DsWidth1_MC.getValV();
	double DsWidth2V=DsWidth2_MC.getValV();
	double DsWidth3V=DsWidth3_MC.getValV();
	double DsGaus1FrV=DsGaus1Fr_MC.getValV();
	double DsGaus2FrV=DsGaus2Fr_MC.getValV();

	// RooRealVar DsMassMean("DsMassMean","DsMassMean",DsMassMeanV,1.965,1.971);
	RooRealVar DsMassMean("DsMassMean","DsMassMean",DsMassMeanV,1.956,1.98);
	RooRealVar DsWidth1("DsWidth1","DsWidth1",DsWidth1V,0.000001,0.1);
	RooRealVar DsWidth2("DsWidth2","DsWidth2",DsWidth2V,0.000001,0.1);
	RooRealVar DsWidth3("DsWidth3","DsWidth3",DsWidth3V,0.000001,0.1);
	RooRealVar DsGaus1Fr("DsGaus1Fr","DsGaus1Fr",DsGaus1FrV,0.0,1);
	RooRealVar DsGaus2Fr("DsGaus2Fr","DsGaus2Fr",DsGaus2FrV,0.0,1);
	RooRealVar DsFloatWidth("DsFloatWidth","DsFloatWidth",0,-0.6,0.6);
	// RooRealVar DsFloatWidth("DsFloatWidth","DsFloatWidth",DsFloatWidthV,DsFloatWidthV,DsFloatWidthV);

	DsWidth1.setConstant(kTRUE);
	DsWidth2.setConstant(kTRUE);
	DsWidth3.setConstant(kTRUE);
	DsGaus1Fr.setConstant(kTRUE);
	DsGaus2Fr.setConstant(kTRUE);


	if(FixShape){
		DsMassMean.setVal(DsMassMeanFix);
		DsWidth1.setVal(DsWidth1Fix);
		DsWidth2.setVal(DsWidth2Fix);
		DsWidth3.setVal(DsWidth3Fix);
		DsGaus1Fr.setVal(DsGaus1FrFix);
		DsGaus2Fr.setVal(DsGaus2FrFix);
		DsFloatWidth.setVal(DsFloatWidthFix);

		DsMassMean.setConstant(kTRUE);
		DsWidth1.setConstant(kTRUE);
		DsWidth2.setConstant(kTRUE);
		DsWidth3.setConstant(kTRUE);
		DsGaus1Fr.setConstant(kTRUE);
		DsGaus2Fr.setConstant(kTRUE);
		DsFloatWidth.setConstant(kTRUE);

	}





	RooFormulaVar scale_width1("scale width1","scaled width1","DsWidth1*(1+DsFloatWidth)",RooArgSet(DsWidth1,DsFloatWidth));
	RooFormulaVar scale_width2("scale width2","scaled width2","DsWidth2*(1+DsFloatWidth)",RooArgSet(DsWidth2,DsFloatWidth));
	RooFormulaVar scale_width3("scale width3","scaled width3","DsWidth3*(1+DsFloatWidth)",RooArgSet(DsWidth3,DsFloatWidth));
	RooGaussian Gaus1("Gaus1","gauss(Dmass,DsMassMean,scale_width1)",Dmass,DsMassMean,scale_width1);
	RooGaussian Gaus2("Gaus2","gauss(Dmass,DsMassMean,scale_width2)",Dmass,DsMassMean,scale_width2);
	RooGaussian Gaus3("Gaus3","gauss(Dmass,DsMassMean,scale_width3)",Dmass,DsMassMean,scale_width3);

	// RooAddPdf SigPdf("SigPdf","SigPdf",RooArgList(Gaus1,Gaus2,Gaus3),RooArgList(DsGaus1Fr,DsGaus2Fr));
	RooAddPdf SigPdf("SigPdf","SigPdf",RooArgList(Gaus1,Gaus2),DsGaus1Fr);
	RooRealVar Cheb1("Cheb1","Cheb1",0,-100,100); // no input is better than with input
	RooRealVar Cheb2("Cheb2","Cheb2",0,-100,100);
	RooRealVar Cheb3("Cheb3","Cheb3",0,-100,100);
	RooChebychev *BkgPdf;
	// BkgPdf=new RooChebychev("BkgPdf","BkgPdf",Dmass,RooArgList(Cheb1));
	BkgPdf=new RooChebychev("BkgPdf","BkgPdf",Dmass,RooArgList(Cheb1,Cheb2));
	// RooRealVar NumSig("NumSig","Number of Signal",1000,-1e4,4e8);
	double NumSigPre=10000;
	double NumBkgPre=3000000;

	RooRealVar NumSig("NumSig","Number of Signal",NumSigPre,0,4e9);
	RooRealVar NumBkg("NumBkg","Number of Background",NumBkgPre,0,1e11);

	// Cheb2.setConstant(kTRUE);
	if(isPbPb==0 || DptLow>=10)
	// if(isPbPb==0)
	{
		Cheb2.setConstant(kTRUE);
	}

	RooAddPdf *RooDsMixPdf= new RooAddPdf("RooDsMixPdf","RooDsMixPdf",RooArgList(SigPdf,*BkgPdf),RooArgList(NumSig,NumBkg));
	RooDsMixPdf->fitTo(RooDS_Data_cut,Extended(kTRUE),NumCPU(20));
	RooDsMixPdf->fitTo(RooDS_Data_cut,Extended(kTRUE),NumCPU(20));

	RooFitResult *fitresult=NULL;
	fitresult=RooDsMixPdf->fitTo(RooDS_Data_cut,Extended(kTRUE),NumCPU(20),Save());

	/// BkgSub
	RooDataSet RooDS_Data_cut2Sig=*(RooDataSet *)RooDS_Data_cut.reduce(RooArgSet(Dmass),Form("Dmass>%f&& Dmass<%f",DsMassMCMean-DsMass3SigWidth,DsMassMCMean+DsMass3SigWidth) ); // Weight is inherited from its mother
	double nCount_2Sig=RooDS_Data_cut2Sig.sumEntries();
	double ncountErr_2Sig=sqrt(nCount_2Sig);
	cout<<"cut sumEntries() = "<<RooDS_Data_cut.sumEntries()<<endl;


  Dmass.setRange("signal",DsMassMCMean-DsMass3SigWidth,DsMassMCMean+DsMass3SigWidth);
  // RooAbsReal *igx=NULL;
  RooAbsReal *igx_sig=NULL;
  RooAbsReal *igx_bkg=NULL;

  igx_sig = SigPdf.createIntegral(Dmass,NormSet(Dmass),Range("signal")) ;
  cout << "Sig_Int[x|signal]_Norm[x] = " << igx_sig->getVal() << endl ;

  igx_bkg = BkgPdf->createIntegral(Dmass,NormSet(Dmass),Range("signal")) ;
  cout << "Bkg_Int[x|signal]_Norm[x] = " << igx_bkg->getVal() << endl ;

  cout<<" \n------------------ \n"<<endl;

  double N_Sig_2sig=NumSig.getValV()* igx_sig->getVal();
  double N_Bkg_2sig=NumBkg.getValV()* igx_bkg->getVal();
  // N_BkgError_2sig=NumBkg.getValV()* igx_bkg->getPropagatedError(fitresultN, RooArgSet(Dmass));
  double N_BkgError_2sig=NumBkg.getValV()* igx_bkg->getPropagatedError(*fitresult);

  cout<<"bkg reltive error = "<<NumBkg.getError()/NumBkg.getValV()<<endl;
  cout<<"bkg reltive error  2 = "<<N_BkgError_2sig/N_Bkg_2sig<<endl;

  double N_Sig_BkgSub=nCount_2Sig-N_Bkg_2sig;
  double N_SigErr_BkgSub=N_Sig_BkgSub*sqrt(N_BkgError_2sig/N_Bkg_2sig*N_BkgError_2sig/N_Bkg_2sig + 1/nCount_2Sig);



	/// end BkgSub


	TCanvas *c_data=new TCanvas("c_data","c_data",c_wtopx,c_wtopx,c_W,800);
	c_data->cd();
//	c_data->SetCanvas
	SetCanvas(c_data);

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

	RooPlot* massframe=new RooPlot("massframe","Dmass",Dmass,DsDataFitRangeLow,DsDataFitRangeHigh,nbin_DmassDraw);
	// RooPlot* massframe=new RooPlot("massframe","Dmass",Dmass,DsDataFitRangeLow,DsDataFitRangeHigh,nbin_DmassDraw);
	RooDS_Data_cut.plotOn(massframe,DataError(RooAbsData::SumW2));
	// RooDsMixPdf->plotOn(massframe,LineColor(2));
	RooDsMixPdf->plotOn(massframe,Components(*BkgPdf),LineColor(4),LineStyle(kDashed));
	RooDsMixPdf->plotOn(massframe,LineColor(2));
	// massframe->GetXaxis()->SetTitle("m_{KK#pi} (GeV/c^{2})");
	massframe->GetXaxis()->SetTitle("");
	massframe->SetTitle("");
	massframe->Draw();


	massframe->Print("V");
	TString s_pdf="RooDsMixPdf_Norm[Dmass]";
	RooGoF goftest(massframe->getHist("h_RooDS_DataAll"),massframe->getCurve(s_pdf));
	cout<<"check 0"<<endl;
	goftest.setRange(Dmass.getMin(),Dmass.getMax());
	goftest.setRebin(5,false); // better rebin to make sure that all bins have >=5 expected events  
	int ndf=0;
	int d_ndf=fitresult->floatParsFinal().getSize();
	double pvalue_BC=0.0;
	double testStat_BC=0.0;
	goftest.BCChi2Test(pvalue_BC,testStat_BC,ndf,d_ndf);
	// file_pval << "BC: " << pvalue_BC << ", " << testStat_BC << endl;
	cout<<"BC: " << pvalue_BC << ", " << testStat_BC << endl;
	double pvalue_RooFit=0.0;
	double testStat_RooFit=0.0;
	goftest.RooFitChi2Test(pvalue_RooFit,testStat_RooFit,ndf,d_ndf);
	// file_pval << "RooFitChi2: " << pvalue_RooFit << ", " << testStat_RooFit << endl;
	cout<<"RooFitChi2: " << pvalue_RooFit << ", " << testStat_RooFit << endl;

	TH1 * htemp=RooDS_Data_cut.createHistogram("htemp",Dmass,Binning(nbin_DmassDraw));
	// RooHist *Roohtemp=massframe->getHist("h_RooDS_DataAll");
	// cout<<"max1 = "<<Roohtemp->GetMaximum();
	// cout<<"max2 = "<<massframe->GetMaximum();
	cout<<"max3 = "<<htemp->GetMaximum();
	// cout<<"min1 = "<<Roohtemp->GetMinimum();
	// cout<<"min2 = "<<massframe->GetMinimum();
	cout<<" , min3 = "<<htemp->GetMinimum();

	double hmax=htemp->GetMaximum();
	double hmin=htemp->GetMinimum();
	double hdiff=hmax-hmin;
	double maxExtra=0.5;
	double minExtra=0.15;
	double newMax=hmax+(maxExtra)*hdiff;
	double newMin=hmin-(minExtra)*hdiff; 
	if(newMax<=0){ newMin=1;}
	massframe->SetMaximum(newMax);
	massframe->SetMinimum(newMin);
	massframe->GetYaxis()->CenterTitle();
	massframe->GetYaxis()->SetTitleOffset(1.1);

	shiftY=0.07;
	shiftX=0.38;
	// tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Ds Data %.0f<p_{T}<%.0f",s_PbPb.Data(),DptLow,DptHigh)); shiftY-=oneshift;
		texCmsPre->Draw("same");
	if(isPbPb==0){
		texColpp->Draw("same");
	}else if(isPbPb==3){
		texColPbPb->Draw("same");
	}

	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("D_{S}^{#pm}  %.0f < p_{T} < %.0f GeV/c",DptLow,DptHigh)); shiftY-=oneshift;
	if(useAnaBin!=1 && useDefaultCut!=1){
		if(var_scan=="PhiMass"){
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Cut : %.4f  ", var_scan.Data(),scan_val )); shiftY-=oneshift;
		}else{
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Cut : %.2f  ", var_scan.Data(),scan_val )); shiftY-=oneshift;
		}
	}
	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Raw Yield = %.1f #pm %.1f",NumSig.getValV(),NumSig.getError() ));  shiftY-=oneshift;
	if(useBkgSub){
	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("BkgSub Yield = %.1f #pm %.1f",N_Sig_BkgSub,N_SigErr_BkgSub ));  shiftY-=oneshift;
	}
	// shiftX=0.06;
	// shiftY=-0.4;	
	// tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("BC pvale : %.2f",pvalue_BC)); shiftY-=oneshift;
	// tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("RooFitChi2 pvale : %.2f",pvalue_RooFit)); shiftY-=oneshift;
	
	if(var_scan=="PhiMass"){
	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s : %.4f", var_scan.Data(), scan_val)); shiftY-=oneshift;
	}else{
	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s : %.2f", var_scan.Data(), scan_val)); shiftY-=oneshift;
	}


  TH1D *hData=new TH1D("hData","hData",1,0,1);
  TH1D *hSig=new TH1D("hSig","hSig",1,0,1);
  TH1D *hBkg=new TH1D("hBkg","hBkg",1,0,1);

  hSig->SetLineColor(2);
  hBkg->SetLineColor(4);

  // TLegend *le_withpull=new TLegend(0.58,0.30,0.8,0.50);
  TLegend *le_withpull=new TLegend(0.58,0.40,0.8,0.65);
  le_withpull->SetBorderSize(0);
	le_withpull->SetTextSize(0.04);
  le_withpull->AddEntry(hData,"Data","pl");
  le_withpull->AddEntry(hSig,"Signal + Background","l");
  le_withpull->AddEntry(hBkg,"Background","l");
  le_withpull->Draw("same");

	pad2->cd();
	// RooPlot* massframe=new RooPlot("massframe","Dmass",Dmass,DsDataFitRangeLow,DsDataFitRangeHigh,nbin_DmassDraw);
	RooPlot *massframe_pull=new RooPlot("massframe_pull","Dmass",Dmass,DsDataFitRangeLow,DsDataFitRangeHigh,nbin_DmassDraw);
	RooHist *hpull=NULL;
	hpull=massframe->pullHist();
	massframe_pull->addPlotable(hpull,"P");

	massframe_pull->SetMaximum(2.8);
	massframe_pull->SetMinimum(-2.8);
	massframe_pull->SetTitle("");
	massframe_pull->GetXaxis()->SetTitle("m_{KK#pi} (GeV/c^2)");
	massframe_pull->GetXaxis()->CenterTitle();
	massframe_pull->GetXaxis()->SetTitleSize(0.14);
	massframe_pull->GetXaxis()->SetLabelSize(0.11); 
	massframe_pull->GetYaxis()->SetTitle("pull");
	massframe_pull->GetYaxis()->CenterTitle();
	massframe_pull->GetYaxis()->SetTitleSize(0.14);
	massframe_pull->GetYaxis()->SetTitleOffset(0.4);
	massframe_pull->GetYaxis()->SetLabelSize(0.11); 
	massframe_pull->Draw();

	TLine *Lm1=new TLine(DsDataFitRangeLow,-1,DsDataFitRangeHigh,-1);
	Lm1->SetLineColor(3);
	Lm1->SetLineStyle(7);
	Lm1->SetLineWidth(2);

	TLine *Lp1=new TLine(DsDataFitRangeLow,1,DsDataFitRangeHigh,1);
	Lp1->SetLineColor(3);
	Lp1->SetLineStyle(7);
	Lp1->SetLineWidth(2);

	Lm1->Draw("same");
	Lp1->Draw("same");


	c_data->SaveAs(Form("./plots/%s/%s_Data%s_Dpt%.0fto%.0f_%s%.0f.png",var_scan.Data(),s_PbPb.Data(),s_FixShape.Data(),DptLow,DptHigh,var_scan.Data(),scan_val*10000));
	SavePlotDirs(c_data,Form("%s_Data%s_Dpt%.0fto%.0f_%s%.0f",s_PbPb.Data(),s_FixShape.Data(),DptLow,DptHigh,var_scan.Data(),scan_val*10000), {"CutScan","Fit",Form("%s",var_scan.Data() )});

	fout->cd();
	TH1D *h_RawYield=new TH1D("h_RawYield","h_RawYield",1,0,1); h_RawYield->Sumw2();
	h_RawYield->SetBinContent(1,NumSig.getValV());
	h_RawYield->SetBinError(1,NumSig.getError());

	if(useBkgSub){	
	h_RawYield->SetBinContent(1,N_Sig_BkgSub);
	h_RawYield->SetBinError(1,N_SigErr_BkgSub);
	}


	h_RawYield->Write("",TObject::kOverwrite);



	if(FixShape==0){
		ofstream f_par;
		gSystem->Exec(Form("mkdir -p par/%s", var_scan.Data()));

		TString fparNameOut=Form("./par/%s/%s_FitPar_Dpt%.0fto%.0f_%s%.0f.txt",var_scan.Data(),s_PbPb.Data(), DptLow,DptHigh,var_scan.Data(),shapeVal*100 );
		if(useAnaBin){
			fparNameOut=Form("./par/%s_FitPar_Dpt%.0fto%.0f_Default.txt",s_PbPb.Data(), DptLow,DptHigh);
		}
		// f_par.open(Form("./par/%s/%s_FitPar_Dpt%.0fto%.0f_%s%.0f.txt",var_scan.Data(),s_PbPb.Data(),DptLow,DptHigh,var_scan.Data(),scan_val*100));
		f_par.open(fparNameOut.Data());

		f_par<<std::setprecision(20);
		f_par<<"GausMean "<<DsMassMean.getValV()<<endl;
		f_par<<"Gauswidth1 "<<DsWidth1.getValV()<<endl;
		f_par<<"Gauswidth2 "<<DsWidth2.getValV()<<endl;
		f_par<<"Gauswidth3 "<<DsWidth3.getValV()<<endl;
		f_par<<"Gausfr1 "<<DsGaus1Fr.getValV()<<endl;
		f_par<<"Gausfr2 "<<DsGaus2Fr.getValV()<<endl;
		f_par<<"Gausfloatwidth "<<DsFloatWidth.getValV()<<endl;
	}


	cout<<"is PbPb : "<<isPbPb<<endl;
	cout<<"useAnaBin : "<<useAnaBin<<endl;
	cout<<"Dpt : "<<DptLow<<" - "<<DptHigh<<endl;
	cout<<"var_scan : "<<var_scan<<endl;
	cout<<"scan_val : "<<scan_val<<endl;
	cout<<"FixShape : "<<FixShape<<endl;
	cout<<"FixShape Val : "<<shapeVal<<endl;
	cout<<"Dalpha_cut = "<<Dalpha_cut<<endl;
	cout<<"Dchi2cl_cut = "<<Dchi2cl_cut<<endl;
	cout<<"Ddls_cut = "<<Ddls_cut<<endl;


	// RooDataSet RooDS_Data_cut2Sig=*(RooDataSet *)RooDS_Data_cut.reduce(RooArgSet(Dmass),Form("Dmass>%f&& Dmass<%f",DsMassMCMean-DsMass3SigWidth,DsMassMCMean+DsMass3SigWidth) ); // Weight is inherited from its mother
	cout<<"cut sumEntries() = "<<RooDS_Data_cut.sumEntries()<<endl;
	cout<<"cut 2sigsumEntries() = "<<RooDS_Data_cut2Sig.sumEntries()<<endl;
	cout<<"all sumEntries() = "<<RooDS_DataAll.sumEntries()<<endl;


	return 0;
}


int main(int argc, char*argv[]){


	cout<<"argc = "<<argc<<endl;

	cout<<"argv[1] = "<<argv[1]<<endl;
	cout<<"argv[2] = "<<argv[2]<<endl;
	cout<<"argv[3] = "<<argv[3]<<endl;
	cout<<"argv[4] = "<<argv[4]<<endl;
	cout<<"argv[5] = "<<argv[5]<<endl;
	cout<<"argv[6] = "<<argv[6]<<endl;
	cout<<"argv[7] = "<<argv[7]<<endl;
	cout<<"argv[8] = "<<argv[8]<<endl;
	cout<<"argv[9] = "<<argv[9]<<endl;

	// int FitCutScan_AnaBin(int isPbPb=0,int useAnaBin=1, int ibin_Dpt=7, double DptLow=6, double DptHigh=40, TString var_scan="Dchi2cl", double scan_val=0.05, int FixShape=0, double shapeVal=3.0 , int useDefaultCut=0){
	if(argc==8){
		FitCutScan_AnaBin(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atof(argv[4]), atof(argv[5]), argv[6], atof(argv[7]) );	
	}else if(argc==11){
		FitCutScan_AnaBin(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atof(argv[4]), atof(argv[5]), argv[6], atof(argv[7]), atoi(argv[8]), atof(argv[9]), atoi(argv[10]) );	
	}else{
		cout<<"wrong number of inputs"<<endl;
		return 1;
	}



	return 0;
}

