#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"

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

using namespace RooStats;


using namespace RooFit;
using namespace std;
  double DsDataFitRangeLow =1.91;
  double DsDataFitRangeHigh = 2.11;

  double shiftY=0.05;
  double shiftX=0.32;
  double oneshift=0.070;

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



int Fit_sPlot(){

  TString dataName="pp_fitTest.root";
  TString mcName_Prompt="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output/MC_eff_pp_Prompt_phikkpi.root";
  TString mcName_NonPrompt="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output/MC_eff_pp_NonPrompt_phikkpi.root";

  TFile *f_data=TFile::Open(dataName.Data());
  TFile *f_mc_Prompt=TFile::Open(mcName_Prompt.Data());
  TFile *f_mc_NonPrompt=TFile::Open(mcName_NonPrompt.Data());

	double DptLow=6;
	double DptHigh=8;

  TTree *t_DsMassDsMC_Prompt=(TTree*)f_mc_Prompt->Get(Form("t_DsMass_pt%.0fto%.0f",DptLow,DptHigh));
  // TTree *t_DsMassDsMC_NonPrompt=(TTree*)f_mc_NonPrompt->Get(Form("t_DsMass_pt%.0fto%.0f",DptLow,DptHigh));
  TTree *t_DsMassData=(TTree*)f_data->Get(Form("t_forDdls"));

	
	// import mc & fit
  RooRealVar Dmass("Dmass","Dmass",DsDataFitRangeLow,DsDataFitRangeHigh);
  // RooRealVar Ddca("Ddca","Ddca",0,0.1); // temp
  RooRealVar Ddls("Ddls","Ddls",0,200); // temp
  RooRealVar D0DataWeight("D0DataWeight","D0DataWeight",0,1e15);
  RooDataSet RooDS_MC("RooDS_MC","RooDS_MC",RooArgSet(Dmass,Ddls),WeightVar(D0DataWeight),Import(*t_DsMassDsMC_Prompt));
  RooDS_MC.Print("v");

  RooRealVar DsMassMean_MC("DsMassMean_MC","DsMassMean_MC",1.9690,1.965,1.971);
  RooRealVar DsWidth1_MC("DsWidth1_MC","DsWidth1_MC",0.011,0.008,0.02);
  RooRealVar DsWidth2_MC("DsWidth2_MC","DsWidth2_MC",0.0060,0.001,0.01);
  RooRealVar DsGaus1Fr_MC("DsGaus1Fr_MC","DsGaus1Fr_MC",0.65,0.3,1);

  RooGaussian Gaus1_MC("Gaus1_MC","gauss(Dmass,DsMassMean_MC,DsWidth1_MC)",Dmass,DsMassMean_MC,DsWidth1_MC);
  RooGaussian Gaus2_MC("Gaus2_MC","gauss(Dmass,DsMassMean_MC,DsWidth2_MC)",Dmass,DsMassMean_MC,DsWidth2_MC);
  RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(Gaus1_MC,Gaus2_MC),DsGaus1Fr_MC);
  SigPdf_MC.fitTo(RooDS_MC,NumCPU(10));
  SigPdf_MC.fitTo(RooDS_MC,NumCPU(10));




	// import data
  RooDataSet RooDSAll("RooDSAll","RooDSAll",RooArgSet(Dmass,Ddls),Import(*t_DsMassData));
  RooDataSet *RooDS=(RooDataSet*)RooDSAll.reduce(RooArgSet(Dmass));
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
  BkgPdf=new RooChebychev("BkgPdf","BkgPdf",Dmass,RooArgList(Cheb1));
  RooRealVar NumSig("NumSig","Number of Signal",1000,-1e4,4e5);
  RooRealVar NumBkg("NumBkg","Number of Background",10000,0,1e8);

  RooAddPdf *RooDsMixPdf= new RooAddPdf("RooDsMixPdf","RooDsMixPdf",RooArgList(SigPdf,*BkgPdf),RooArgList(NumSig,NumBkg));
  RooDsMixPdf->fitTo(*RooDS,Extended(kTRUE),NumCPU(20));
  RooDsMixPdf->fitTo(*RooDS,Extended(kTRUE),NumCPU(20));
  RooDsMixPdf->fitTo(*RooDS,Extended(kTRUE),NumCPU(20));

	TCanvas *c_test=new TCanvas("c_test","c_test",800,600);
	c_test->cd();

  RooPlot* massframe=new RooPlot("massframe","Dmass",Dmass,DsDataFitRangeLow,DsDataFitRangeHigh,50);
  // RooPlot* massframe=new RooPlot("massframe","Dmass",Dmass,DsDataFitRangeLow,DsDataFitRangeHigh,nbin_DmassDraw);
  RooDS->plotOn(massframe,DataError(RooAbsData::SumW2));
  RooDsMixPdf->plotOn(massframe,LineColor(2));
  massframe->Draw();


	RooMsgService::instance().setSilentMode(true);


  RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",*RooDS, RooDsMixPdf, RooArgList(NumSig,NumBkg) );

   std::cout << std::endl <<  "Yield of Sig is "
               << NumSig.getVal() << ".  From sWeights it is "
               << sData->GetYieldFromSWeight("NumSig") << std::endl;
   std::cout << "Yield of Bkg is "
               << NumBkg.getVal() << ".  From sWeights it is "
               << sData->GetYieldFromSWeight("NumBkg") << std::endl
               << std::endl;

   std::cout << "import new dataset with sWeights" << std::endl;
//   ws->import(*data, Rename("dataWithSWeights"));
//	RooDataSet *data_wts=(RooDataSet)
  RooDsMixPdf->fitTo(*RooDS,Extended(kTRUE),NumCPU(20));

	// plot the distribution

	TCanvas *c_Ddls= new TCanvas("c_Ddls","c_Ddls",800,600);
	c_Ddls->cd();
   RooDataSet * dataw_s = new RooDataSet(RooDS->GetName(),RooDS->GetTitle(),RooDS,*RooDS->get(),0,"NumSig_sw") ;
	cout<<"check 1"<<endl;

   RooPlot* frame2 = Ddls.frame() ;
	cout<<"check 2"<<endl;
   dataw_s->plotOn(frame2, DataError(RooAbsData::SumW2) ) ;
	cout<<"check 3"<<endl;
   frame2->SetTitle("Ddls distribution for Sig");
   frame2->Draw() ;





	return 0;

}



