#include "MkkFit_para.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"

#include "TFitter.h"
#include "TFitResult.h"
#include <cmath>

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
#include "RooCBShape.h"
#include "RooBreitWigner.h"
#include "RooGoF.C"

using namespace RooFit;




 TF1* fitDsMassMC(TH1D* h_DsMassDsMC);
TF1* fitDsMass(TH1D *hData, TF1 *f1MC,double &FitYield,double &FitYieldErr,double &BkgSubYield,double &BkgSubYieldErr,  TString s_ExtraName="",int fixSigShape=0, TF1 *f1_default=NULL);
//	gSystem->Exec("mkdir -p plots/fit");


  double DsMassCandWidth=0.05;
  double DsMassCandFitMeanWidth=0.02;
  double shiftY=0.05;
  double shiftX=0.32;
  double oneshift=0.070;

	double DptLowGL=5;
	double DptHighGL=40;
	TString s_ppPbPb="pp";
	TString s_type="";


int FitMkkpi(int isPbPb=3, double DptLow=6, double DptHigh=40){

	DptLowGL=DptLow;
	DptHighGL=DptHighGL;

	s_ppPbPb="pp";

	if(isPbPb==3)
	{
		s_ppPbPb="PbPb";
	}



  gSystem->Exec("mkdir -p plots/fit");

//	TString s_ppPbPb="pp";

	TString s_finName=Form("./his/%s_Dpt%.0fto%.0f.root",s_ppPbPb.Data(),DptLow,DptHigh);
	TString s_foutName=Form("./his/%s_Dpt%.0fto%.0f_mkkpiFit.root",s_ppPbPb.Data(),DptLow,DptHigh);

	s_type=Form("%s_Dpt%.0fto%.0f",s_ppPbPb.Data(),DptLow,DptHigh);

	TFile *fin=TFile::Open(s_finName.Data(),"READ");
	TFile *fout=TFile::Open(s_foutName.Data(),"RECREATE");

	TH1D *h_mkkpi_default_MCP_Phi=(TH1D*)fin->Get("h_mkkpi_default_MCP_Phi"); 
	TH1D *h_mkkpi_default_MCP_f0 =(TH1D*)fin->Get("h_mkkpi_default_MCP_f0");
	TH1D *h_mkk_MCP_Phi          =(TH1D*)fin->Get("h_mkk_MCP_Phi");
	TH1D *h_mkk_MCP_f0           =(TH1D*)fin->Get("h_mkk_MCP_f0");
	TH1D *h_mkk_MCP_Phi_fineBin  =(TH1D*)fin->Get("h_mkk_MCP_Phi_fineBin");
	TH1D *h_mkk_MCP_f0_fineBin   =(TH1D*)fin->Get("h_mkk_MCP_f0_fineBin");
	TH1D *h_mkkpi_default        =(TH1D*)fin->Get("h_mkkpi_default");



	TF1 *f1MC=fitDsMassMC(h_mkkpi_default_MCP_Phi);

	cout<<"f1MC->GetParameter(1) = "<<f1MC->GetParameter(1)<<endl;

	double FitYield;
	double FitYieldErr;
	double BkgSubYield;
	double BkgSubYieldErr;

	TF1 *f1Data=fitDsMass(h_mkkpi_default,f1MC, FitYield,FitYieldErr,BkgSubYield,BkgSubYieldErr);

	// RooFitTest 

	RooRealVar Dmass("Dmass","Dmass",1.91,2.11);
	RooDataHist dh_MC("dh_MC","dh_MC",Dmass,Import(*h_mkkpi_default_MCP_Phi)); 
	RooPlot *MC_frame=Dmass.frame(Title("MC"));
	dh_MC.plotOn(MC_frame);

  RooRealVar DsMassMean_MC("DsMassMean_MC","DsMassMean_MC",1.9690,1.965,1.971);
  RooRealVar DsWidth1_MC("DsWidth1_MC","DsWidth1_MC",0.011,0.003,0.1);
  RooRealVar DsWidth2_MC("DsWidth2_MC","DsWidth2_MC",0.0060,0.0001,0.02);
  RooRealVar DsWidth3_MC("DsWidth3_MC","DsWidth3_MC",0.0060,0.0001,0.02);
   RooRealVar DsGaus1Fr_MC("DsGaus1Fr_MC","DsGaus1Fr_MC",0.65,0.01,1);
  RooRealVar DsGaus2Fr_MC("DsGaus2Fr_MC","DsGaus2Fr_MC",0.65,0.01,1);

  RooGaussian Gaus1_MC("Gaus1_MC","gauss(Dmass,DsMassMean_MC,DsWidth1_MC)",Dmass,DsMassMean_MC,DsWidth1_MC);
  RooGaussian Gaus2_MC("Gaus2_MC","gauss(Dmass,DsMassMean_MC,DsWidth2_MC)",Dmass,DsMassMean_MC,DsWidth2_MC);
  RooGaussian Gaus3_MC("Gaus3_MC","gauss(Dmass,DsMassMean_MC,DsWidth3_MC)",Dmass,DsMassMean_MC,DsWidth3_MC);
  // RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(Gaus1_MC,Gaus2_MC),DsGaus1Fr_MC);
  RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(Gaus1_MC,Gaus2_MC,Gaus3_MC),RooArgList(DsGaus1Fr_MC,DsGaus2Fr_MC));
  SigPdf_MC.fitTo(dh_MC,NumCPU(20));
  SigPdf_MC.fitTo(dh_MC,NumCPU(20));
  SigPdf_MC.plotOn(MC_frame,LineColor(2));
  TCanvas *c_MC= new TCanvas("c_MC","c_MC");
  c_MC->cd();
  MC_frame->GetXaxis()->SetTitle("m_{KK#pi} (GeV/c^{2})");
  MC_frame->Draw();

	// Data Fit

	RooDataHist RooDS_Data_cut("RooDS_Data_cut","RooDS_Data_cut",Dmass,Import(*h_mkkpi_default));

  double DsMassMeanV=DsMassMean_MC.getValV();
  double DsWidth1V=DsWidth1_MC.getValV();
  double DsWidth2V=DsWidth2_MC.getValV();
  double DsWidth3V=DsWidth3_MC.getValV();
  double DsGaus1FrV=DsGaus1Fr_MC.getValV();
  double DsGaus2FrV=DsGaus2Fr_MC.getValV();
  
  RooRealVar DsMassMean("DsMassMean","DsMassMean",DsMassMeanV,1.965,1.971);
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

 RooFormulaVar scale_width1("scale width1","scaled width1","DsWidth1*(1+DsFloatWidth)",RooArgSet(DsWidth1,DsFloatWidth));
  RooFormulaVar scale_width2("scale width2","scaled width2","DsWidth2*(1+DsFloatWidth)",RooArgSet(DsWidth2,DsFloatWidth));
  RooFormulaVar scale_width3("scale width3","scaled width3","DsWidth3*(1+DsFloatWidth)",RooArgSet(DsWidth3,DsFloatWidth));
  RooGaussian Gaus1("Gaus1","gauss(Dmass,DsMassMean,scale_width1)",Dmass,DsMassMean,scale_width1);
  RooGaussian Gaus2("Gaus2","gauss(Dmass,DsMassMean,scale_width2)",Dmass,DsMassMean,scale_width2);
  RooGaussian Gaus3("Gaus3","gauss(Dmass,DsMassMean,scale_width3)",Dmass,DsMassMean,scale_width3);

  RooAddPdf SigPdf("SigPdf","SigPdf",RooArgList(Gaus1,Gaus2,Gaus3),RooArgList(DsGaus1Fr,DsGaus2Fr));
  // RooAddPdf SigPdf("SigPdf","SigPdf",RooArgList(Gaus1,Gaus2),DsGaus1Fr);
  RooRealVar Cheb1("Cheb1","Cheb1",0,-100,100); // no input is better than with input
  RooRealVar Cheb2("Cheb2","Cheb2",0,-100,100);
  RooRealVar Cheb3("Cheb3","Cheb3",0,-100,100);
  RooChebychev *BkgPdf;
  // BkgPdf=new RooChebychev("BkgPdf","BkgPdf",Dmass,RooArgList(Cheb1));
  BkgPdf=new RooChebychev("BkgPdf","BkgPdf",Dmass,RooArgList(Cheb1,Cheb2));

	if(isPbPb==0){
		Cheb2.setConstant(kTRUE);
	}
	
  // RooRealVar NumSig("NumSig","Number of Signal",1000,-1e4,4e8);
  double NumSigPre=10000;
  double NumBkgPre=3000000;

  RooRealVar NumSig("NumSig","Number of Signal",NumSigPre,-1e3,4e9);
  RooRealVar NumBkg("NumBkg","Number of Background",NumBkgPre,-1e3,1e11);

 RooAddPdf *RooDsMixPdf= new RooAddPdf("RooDsMixPdf","RooDsMixPdf",RooArgList(SigPdf,*BkgPdf),RooArgList(NumSig,NumBkg));
  RooDsMixPdf->fitTo(RooDS_Data_cut,Extended(kTRUE),NumCPU(20));
  RooDsMixPdf->fitTo(RooDS_Data_cut,Extended(kTRUE),NumCPU(20));

  RooFitResult *fitresult=NULL;
  fitresult=RooDsMixPdf->fitTo(RooDS_Data_cut,Extended(kTRUE),NumCPU(20),Save());

  TCanvas *c_data=new TCanvas("c_data","c_data");
  c_data->cd();
  // RooPlot* massframe=new RooPlot("massframe","Dmass",Dmass,DsDataFitRangeLow,DsDataFitRangeHigh,nbin_DmassDraw);
  RooPlot* massframe=Dmass.frame();
  RooDS_Data_cut.plotOn(massframe,DataError(RooAbsData::SumW2));
  RooDsMixPdf->plotOn(massframe,LineColor(2));
  massframe->GetXaxis()->SetTitle("m_{KK#pi} (GeV/c^{2})");
  massframe->SetTitle("");
  massframe->Draw();

	gSystem->Exec("mkdir -p plots/roofit");
  // c_data->SaveAs(Form("./plots/roofit/%s_Data_Mkkpi_%s.png",s_type.Data(),s_ExtraName.Data()));
  c_data->SaveAs(Form("./plots/roofit/%s_Data_Mkkpi.png",s_type.Data()));
	// return 1;

	
	TH1D *h_mkkpi_mkkbins[nbin_mkk];
	TF1 *f1Data_mkkbins[nbin_mkk];

	TString extraName="";

	fout->cd();
	TH1D *h_NDs_mkkbin_Fit=new TH1D("h_NDs_mkkbin_Fit",";m_{KK}",nbin_mkk,bins_mkk); h_NDs_mkkbin_Fit->Sumw2();
	TH1D *h_NDs_mkkbin_BkgSub=new TH1D("h_NDs_mkkbin_BkgSub",";m_{KK}",nbin_mkk,bins_mkk); h_NDs_mkkbin_BkgSub->Sumw2();

	TH1D *h_NDs_mkkbin_Fit_roofit=new TH1D("h_NDs_mkkbin_Fit_roofit",";m_{KK}",nbin_mkk,bins_mkk); h_NDs_mkkbin_Fit_roofit->Sumw2();
	TH1D *h_NDs_mkkbin_BkgSub_roofit=new TH1D("h_NDs_mkkbin_BkgSub_roofit",";m_{KK}",nbin_mkk,bins_mkk); h_NDs_mkkbin_BkgSub_roofit->Sumw2();

	RooDataHist *RooDS_Data_cutN[nbin_mkk];
	RooArgSet nset(Dmass);
	RooAbsReal *igx=NULL;
	RooAbsReal *igx_sig=NULL;
	RooAbsReal *igx_bkg=NULL;
	double N_Sig_2sig;
	double N_Bkg_2sig;
	double N_BkgError_2sig;
	double N_All_2sig;

	double N_Sig_BkgSub;
	double N_SigErr_BkgSub;



	for(int i=0; i<nbin_mkk; i++)
	 // for(int i=0; i<1; i++)
	{
		h_mkkpi_mkkbins[i]=(TH1D*)fin->Get(Form("h_mkkpi_mkkbins_%i",i));
		extraName=Form("mkk%.0fto%.0f",bins_mkk[i]*1000,bins_mkk[i+1]*1000);
		f1Data_mkkbins[i]=fitDsMass(h_mkkpi_mkkbins[i],f1MC,FitYield,FitYieldErr,BkgSubYield,BkgSubYieldErr,extraName,1,f1Data);

		h_NDs_mkkbin_Fit->SetBinContent(i+1,FitYield);
		h_NDs_mkkbin_Fit->SetBinError(i+1,FitYieldErr);

		h_NDs_mkkbin_BkgSub->SetBinContent(i+1,BkgSubYield);
		h_NDs_mkkbin_BkgSub->SetBinError(i+1,BkgSubYieldErr);


	RooDS_Data_cutN[i]=new RooDataHist(Form("RooDS_Data_cut%i",i),Form("RooDS_Data_cut%i",i),Dmass,Import(*h_mkkpi_mkkbins[i]));
/*
  double DsMassMeanV=DsMassMean_MC.getValV();
  double DsWidth1V=DsWidth1_MC.getValV();
  double DsWidth2V=DsWidth2_MC.getValV();
  double DsWidth3V=DsWidth3_MC.getValV();
  double DsGaus1FrV=DsGaus1Fr_MC.getValV();
  double DsGaus2FrV=DsGaus2Fr_MC.getValV();
  
  RooRealVar DsMassMean("DsMassMean","DsMassMean",DsMassMeanV,1.965,1.971);
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

 RooFormulaVar scale_width1("scale width1","scaled width1","DsWidth1*(1+DsFloatWidth)",RooArgSet(DsWidth1,DsFloatWidth));
  RooFormulaVar scale_width2("scale width2","scaled width2","DsWidth2*(1+DsFloatWidth)",RooArgSet(DsWidth2,DsFloatWidth));
  RooFormulaVar scale_width3("scale width3","scaled width3","DsWidth3*(1+DsFloatWidth)",RooArgSet(DsWidth3,DsFloatWidth));
  RooGaussian Gaus1("Gaus1","gauss(Dmass,DsMassMean,scale_width1)",Dmass,DsMassMean,scale_width1);
  RooGaussian Gaus2("Gaus2","gauss(Dmass,DsMassMean,scale_width2)",Dmass,DsMassMean,scale_width2);
  RooGaussian Gaus3("Gaus3","gauss(Dmass,DsMassMean,scale_width3)",Dmass,DsMassMean,scale_width3);

  RooAddPdf SigPdf("SigPdf","SigPdf",RooArgList(Gaus1,Gaus2,Gaus3),RooArgList(DsGaus1Fr,DsGaus2Fr));
  // RooAddPdf SigPdf("SigPdf","SigPdf",RooArgList(Gaus1,Gaus2),DsGaus1Fr);
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

 RooAddPdf *RooDsMixPdf= new RooAddPdf("RooDsMixPdf","RooDsMixPdf",RooArgList(SigPdf,*BkgPdf),RooArgList(NumSig,NumBkg));
*/
	DsFloatWidth.setConstant(kTRUE);
	DsMassMean.setConstant(kTRUE);

  RooDsMixPdf->fitTo(*RooDS_Data_cutN[i],Extended(kTRUE),NumCPU(20));
   RooDsMixPdf->fitTo(*RooDS_Data_cutN[i],Extended(kTRUE),NumCPU(20));

  RooFitResult *fitresultN=NULL;
  fitresultN=RooDsMixPdf->fitTo(*RooDS_Data_cutN[i],Extended(kTRUE),NumCPU(20),Save());

  TCanvas *c_data=new TCanvas("c_data","c_data");
  c_data->cd();
  // RooPlot* massframe=new RooPlot("massframe","Dmass",Dmass,DsDataFitRangeLow,DsDataFitRangeHigh,nbin_DmassDraw);
  RooPlot* massframe=Dmass.frame();
  RooDS_Data_cutN[i]->plotOn(massframe,DataError(RooAbsData::SumW2));
  RooDsMixPdf->plotOn(massframe,LineColor(2));
  massframe->GetXaxis()->SetTitle("m_{KK#pi} (GeV/c^{2})");
  massframe->SetTitle("");
  massframe->Draw();


	// cal bkg sub
	N_All_2sig=h_mkkpi_mkkbins[i]->Integral(h_mkkpi_mkkbins[i]->FindBin(DsMassMCMean-DsMass2SigWidth),h_mkkpi_mkkbins[i]->FindBin(DsMassMCMean+DsMass2SigWidth));

  // Dmass.setRange("signal",DsMassMean-DsMass2SigWidth,DsMassMean+DsMass2SigWidth);
  Dmass.setRange("signal",h_mkkpi_mkkbins[i]->GetBinLowEdge(h_mkkpi_mkkbins[i]->FindBin(DsMassMCMean-DsMass2SigWidth)) , h_mkkpi_mkkbins[i]->GetBinLowEdge(h_mkkpi_mkkbins[i]->FindBin(DsMassMCMean+DsMass2SigWidth)) +  h_mkkpi_mkkbins[i]->GetBinWidth(1) );
		cout<<"LowEdge : "<<h_mkkpi_mkkbins[i]->GetBinLowEdge(h_mkkpi_mkkbins[i]->FindBin(DsMassMCMean-DsMass2SigWidth))<<endl;
		cout<<"HiEdge : "<<h_mkkpi_mkkbins[i]->GetBinLowEdge(h_mkkpi_mkkbins[i]->FindBin(DsMassMCMean+DsMass2SigWidth))+h_mkkpi_mkkbins[i]->GetBinWidth(1) <<endl;


  igx_sig = SigPdf.createIntegral(Dmass,NormSet(Dmass),Range("signal")) ;
  cout << "Sig_Int[x|signal]_Norm[x] = " << igx_sig->getVal() << endl ;

  igx_bkg = BkgPdf->createIntegral(Dmass,NormSet(Dmass),Range("signal")) ;
  cout << "Bkg_Int[x|signal]_Norm[x] = " << igx_bkg->getVal() << endl ;


  cout<<" \n------------------ \n"<<endl;

  N_Sig_2sig=NumSig.getValV()* igx_sig->getVal();
  N_Bkg_2sig=NumBkg.getValV()* igx_bkg->getVal();
  // N_BkgError_2sig=NumBkg.getValV()* igx_bkg->getPropagatedError(fitresultN, RooArgSet(Dmass));
  N_BkgError_2sig=NumBkg.getValV()* igx_bkg->getPropagatedError(*fitresultN);

	cout<<"bkg reltive error = "<<NumBkg.getError()/NumBkg.getValV()<<endl;
	cout<<"bkg reltive error  2 = "<<N_BkgError_2sig/N_Bkg_2sig<<endl;
	
	N_Sig_BkgSub=N_All_2sig-N_Bkg_2sig;
	// N_SigErr_BkgSub=N_Sig_BkgSub*sqrt(NumBkg.getError()/NumBkg.getValV()* NumBkg.getError()/NumBkg.getValV() + 1/N_All_2sig); 
	N_SigErr_BkgSub=N_Sig_BkgSub*sqrt(N_BkgError_2sig/N_Bkg_2sig*N_BkgError_2sig/N_Bkg_2sig + 1/N_All_2sig); 


  TLatex *tl_roofitData =new TLatex();
//  tl_roofitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s_Data_roofit",str_PbPb.Data()));  shiftY-=oneshift;
	shiftY=0;

	if(isPbPb){
	shiftY=-0.25;
	}
  tl_roofitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s D_{S} %.0f < pt < %.0f ",s_ppPbPb.Data(),DptLowGL,DptHighGL)); shiftY-=oneshift;
  tl_roofitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Fit Yield = %.1f #pm %.1f ",NumSig.getValV(),NumSig.getError() ));  shiftY-=oneshift;
  tl_roofitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("BkgSub Yield = %.1f #pm %.1f ",N_Sig_BkgSub,N_SigErr_BkgSub ));  shiftY-=oneshift;
  tl_roofitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s",extraName.Data() ));  shiftY-=oneshift;




	gSystem->Exec("mkdir -p plots/roofit");
  // c_data->SaveAs(Form("./plots/roofit/%s_Data_Mkkpi_%s.png",s_type.Data(),s_ExtraName.Data()));
  c_data->SaveAs(Form("./plots/roofit/%s_Data_Mkkpi_bin%i.png",s_type.Data(),i));

	h_NDs_mkkbin_Fit_roofit->SetBinContent(i+1,NumSig.getValV());
	h_NDs_mkkbin_Fit_roofit->SetBinError(i+1,NumSig.getError());


	h_NDs_mkkbin_BkgSub_roofit->SetBinContent(i+1,N_Sig_BkgSub);
	h_NDs_mkkbin_BkgSub_roofit->SetBinError(i+1,N_SigErr_BkgSub);


	delete tl_roofitData;
	delete c_data;
	delete fitresultN;


	}

	// cout<<"FitYield = "<<FitYield<<endl;

		h_NDs_mkkbin_Fit->Write();
		h_NDs_mkkbin_BkgSub->Write();

		h_NDs_mkkbin_Fit_roofit->Write();
		h_NDs_mkkbin_BkgSub_roofit->Write();


	TLatex *tlatex=new TLatex();
	TCanvas *c_Nds_mkkbin_Fit=new TCanvas("c_Nds_mkkbin_Fit","c_Nds_mkkbin_Fit");
	c_Nds_mkkbin_Fit->cd();
	h_NDs_mkkbin_Fit->Draw();

	shiftY=0;
  tlatex->DrawLatexNDC(textposx+0.08,textposy+shiftY,Form("%s D_{S} %.0f < pt < %.0f ",s_ppPbPb.Data(),DptLowGL,DptHighGL)); shiftY-=oneshift;

	c_Nds_mkkbin_Fit->SaveAs(Form("./plots/fit/%s_Nds_mkkbin_Fit_Dpt%.0fto%.0f.png",s_type.Data(),DptLowGL,DptHighGL));

	TCanvas *c_Nds_mkkbin_BkgSub=new TCanvas("c_Nds_mkkbin_BkgSub","c_Nds_mkkbin_BkgSub");
	c_Nds_mkkbin_BkgSub->cd();
	h_NDs_mkkbin_BkgSub->Draw();

	shiftY=0;
  tlatex->DrawLatexNDC(textposx+0.08,textposy+shiftY,Form("%s D_{S} %.0f < pt < %.0f ",s_ppPbPb.Data(),DptLowGL,DptHighGL)); shiftY-=oneshift;

	c_Nds_mkkbin_BkgSub->SaveAs(Form("./plots/fit/%s_Nds_mkkbin_BkgSub_Dpt%.0fto%.0f.png",s_type.Data(),DptLowGL,DptHighGL));



	cout<<"s_ppPbPb = "<<s_ppPbPb.Data()<<" , isPbPb = "<<isPbPb<<endl;


	return 0;

}


TF1* fitDsMassMC(TH1D* h_DsMassDsMC){

  TF1 *f1_DsMC = new TF1("f1_DsMC","[0]*( [4]*TMath::Gaus(x,[1],[2]*(1+[5]))/( sqrt(2*TMath::Pi())*[2]*(1+[5]) ) + (1-[4])*TMath::Gaus(x,[1],[3]*(1+[5]))/(sqrt(2*TMath::Pi())*[3]*(1+[5])  ) )",DsMassFitLow,DsMassFitHigh); // kTRUE nomalize Gaus
  f1_DsMC->SetParameter(0, h_DsMassDsMC->GetBinWidth(1) * h_DsMassDsMC->Integral( h_DsMassDsMC->FindBin(DsMass-DsMassCandWidth) , h_DsMassDsMC->FindBin(DsMass+DsMassCandWidth) ) ) ;
  f1_DsMC->SetParameter(1,DsMass);   // DsMean
  f1_DsMC->SetParLimits(1,DsMass-DsMassCandFitMeanWidth,DsMass+DsMassCandFitMeanWidth);
  f1_DsMC->SetParameter(2,0.01); // defaut (wider) gaus1 sigma
  f1_DsMC->SetParLimits(2,0.0001,0.1);
  f1_DsMC->SetParameter(3,0.001); // default (narrow) gaus2 sigma
  f1_DsMC->SetParLimits(3,0.00001,0.015);
  f1_DsMC->SetParameter(4,0.5);  // fraction of Gaus1 / All
  f1_DsMC->SetParLimits(4,0,1);
  f1_DsMC->SetParameter(5,0); // for data /mc discrepancy
  f1_DsMC->SetParLimits(5,0,0);
  f1_DsMC->FixParameter(5,0);

  f1_DsMC->SetLineColor(kRed);

  h_DsMassDsMC->Fit("f1_DsMC","QN0","",    DsMass-DsMassCandWidth,DsMass+DsMassCandWidth);
  h_DsMassDsMC->Fit("f1_DsMC","QN0","",    f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);
  h_DsMassDsMC->Fit("f1_DsMC","M IN0","",    f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);
  h_DsMassDsMC->Fit("f1_DsMC","L QN0","",  f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);
  h_DsMassDsMC->Fit("f1_DsMC","LN0","",    f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);
  h_DsMassDsMC->Fit("f1_DsMC","L MN0","",  f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);
  h_DsMassDsMC->Fit("f1_DsMC","L MN0","",  f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);
  h_DsMassDsMC->Fit("f1_DsMC","L M IN0","",f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);
  h_DsMassDsMC->Fit("f1_DsMC","L M IN0","",f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);
  h_DsMassDsMC->Fit("f1_DsMC","L M IN0","",f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);
  h_DsMassDsMC->Fit("f1_DsMC","L ME I0 S","",f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);

 double yieldDsMC= f1_DsMC->Integral(f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth)/ h_DsMassDsMC->GetBinWidth(1);
 double yieldDsMCErr = f1_DsMC->Integral(DsMass-DsMassCandWidth,DsMass+DsMassCandWidth)/h_DsMassDsMC->GetBinWidth(1) * f1_DsMC->GetParError(0)/f1_DsMC->GetParameter(0);

	TCanvas *c_MCfit=new TCanvas("c_MCfit","c_MCfit");
	c_MCfit->cd();
	h_DsMassDsMC->Draw();

  shiftY=0;
	TLatex *tl_binfitMC=new TLatex();
  tl_binfitMC->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s D_{S} %.0f < pt < %.0f ",s_ppPbPb.Data(),DptLowGL,DptHighGL)); shiftY-=oneshift;
 tl_binfitMC->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Mean=%.4f, frGaus1=%.3f  ",f1_DsMC->GetParameter(1),f1_DsMC->GetParameter(4) ));  shiftY-=oneshift;
  tl_binfitMC->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Sigma1=%.4f, Sigma2=%.4f",f1_DsMC->GetParameter(2), f1_DsMC->GetParameter(3)));  shiftY-=oneshift;



	f1_DsMC->Draw("same");


	c_MCfit->SaveAs(Form("./plots/fit/%s_MC_Mkkpi.png",s_type.Data()));

 return f1_DsMC;

}



TF1* fitDsMass(TH1D *h_DsMassData, TF1 *f1MC,double &FitYield,double &FitYieldErr,double &BkgSubYield,double &BkgSubYieldErr, TString s_ExtraName, int fixSigShape, TF1 *f1_default){

  // TF1* f = new TF1("fMass","[0]*([7]*([9]*TMath::Gaus(x,[1],[2]*(1+[11]))/(sqrt(2*3.1415927)*[2]*(1+[11]))+(1-[9])*TMath::Gaus(x,[1],[10]*(1+[11]))/(sqrt(2*3.1415927)*[10]*(1+[11])))+(1-[7]\
)*TMath::Gaus(x,[1],[8])/(sqrt(2*3.1415927)*[8]))+[3]+[4]*x+[5]*x*x+[6]*x*x*x", 1.7, 2.0);

	// FitYield=472;

	cout<<"check -1"<<endl;
  cout<<"f1MC->GetParameter(1) = "<<f1MC->GetParameter(1)<<endl;
	

  TF1 *f1_DsMC = new TF1("f1_DsMC","[0]*( [4]*TMath::Gaus(x,[1],[2]*(1+[5]))/( sqrt(2*TMath::Pi())*[2]*(1+[5]) ) + (1-[4])*TMath::Gaus(x,[1],[3]*(1+[5]))/(sqrt(2*TMath::Pi())*[3]*(1+[5])  ) )",DsMassFitLow,DsMassFitHigh); // kTRUE nomalize Gaus

  // f1_DsMC->SetParameter(0, h_DsMassDsMC->GetBinWidth(1) * h_DsMassDsMC->Integral( h_DsMassDsMC->FindBin(DsMass-DsMassCandWidth) , h_DsMassDsMC->FindBin(DsMass-DsMassCandWidth) ) ) ;
  f1_DsMC->SetParameter(1,f1MC->GetParameter(1));   // DsMean
  f1_DsMC->SetParLimits(1,DsMass-DsMassCandFitMeanWidth,DsMass+DsMassCandFitMeanWidth);
  f1_DsMC->SetParameter(2,f1MC->GetParameter(2)); // defaut (wider) gaus1 sigma
 f1_DsMC->SetParLimits(2,0.0001,0.1);
  f1_DsMC->SetParameter(3,f1MC->GetParameter(3)); // default (narrow) gaus2 sigma
  f1_DsMC->SetParLimits(3,0.00001,0.1);
  f1_DsMC->SetParameter(4,f1MC->GetParameter(4));  // fraction of Gaus1 / All
  f1_DsMC->SetParLimits(4,0,1);
  f1_DsMC->SetParameter(5,0); // for data /mc discrepancy
//  f1_DsMC->SetParLimits(5,0,0);
//  f1_DsMC->FixParameter(5,0);

  f1_DsMC->SetLineColor(kRed);


	cout<<"check 0"<<endl;



  TF1 *f1_DsSignal = (TF1*)f1_DsMC->Clone("f1_DsSignal");
  TF1 *f1_DsBkg= new TF1("f1_DsBkg","[0]*(1+[1]*x+[2]*(2*x*x-1))");
  TF1 *f1_DsBkg_clone=(TF1*)f1_DsBkg->Clone("f1_DsBkg_clone");  // stupid bug in root, must use clone function
  TF1 *f1_DsMix = new TF1("f1_DsMix","f1_DsSignal+f1_DsBkg");

	// f1_DsMix->SetParLimits(0,0,1000000);
  f1_DsMix->SetParameter(5,0); // width ratio difference for Data/ MC, data width = (1+5)*MC width
  f1_DsMix->SetParLimits(5,-1,1);
	f1_DsMix->SetParameter(0,0.1*h_DsMassData->Integral()*h_DsMassData->GetBinWidth(1));
  f1_DsMix->FixParameter(1,f1_DsMC->GetParameter(1));  // Ds Mass mean
  f1_DsMix->FixParameter(2,f1_DsMC->GetParameter(2));  // Ds Sigma1
  f1_DsMix->FixParameter(3,f1_DsMC->GetParameter(3));  // Ds Sigma2
  f1_DsMix->FixParameter(4,f1_DsMC->GetParameter(4));  // Ds Gaus1/all yield ratio
  f1_DsMix->SetLineColor(kRed);
//  f1_DsMix->SetParameter(6,1000);
  f1_DsMix->SetParameter(7,0);
  f1_DsMix->SetParameter(6,0.9*h_DsMassData->Integral()*h_DsMassData->GetBinWidth(1));


	if(fixSigShape==1 && f1_default!=NULL){
		cout<<"fixSigShape fit"<<endl;
		// cout<<"f1_default->GetParameter(7) = "<<f1_default->GetParameter(7)<<endl;
		// cout<<"f1_default->GetParameter(8) = "<<f1_default->GetParameter(8)<<endl;

		f1_DsMix->SetParameters(0.1*h_DsMassData->Integral()*h_DsMassData->GetBinWidth(1), f1_default->GetParameter(1), f1_default->GetParameter(2), f1_default->GetParameter(3), f1_default->GetParameter(4), f1_default->GetParameter(5), f1_default->GetParameter(6), f1_default->GetParameter(7), f1_default->GetParameter(8) );

		f1_DsMix->FixParameter(1,f1_default->GetParameter(1));  // Ds Mass mean
		f1_DsMix->FixParameter(2,f1_default->GetParameter(2));  
		f1_DsMix->FixParameter(3,f1_default->GetParameter(3));  
		f1_DsMix->FixParameter(4,f1_default->GetParameter(4));  
		f1_DsMix->FixParameter(5,f1_default->GetParameter(5));  

		f1_DsMix->FixParameter(7,f1_default->GetParameter(7));  
		f1_DsMix->FixParameter(8,f1_default->GetParameter(8));  

	}

	cout<<"check 1"<<endl;

  h_DsMassData->Fit("f1_DsMix","QN0","",    DsMassFitLow, DsMassFitHigh );
  h_DsMassData->Fit("f1_DsMix","QN0","",    DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L QN0","",  DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","LN0","",    DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L MN0","",  DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L MN0","",  DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L M IN0","",DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L M IN0","",DsMassFitLow, DsMassFitHigh);

	if(fixSigShape!=1){
  f1_DsMix->ReleaseParameter(5);
  f1_DsMix->ReleaseParameter(1);
	}

  h_DsMassData->Fit("f1_DsMix","QN0","",    DsMassFitLow, DsMassFitHigh );
  h_DsMassData->Fit("f1_DsMix","QN0","",    DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L QN0","",  DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","LN0","",    DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L MN0","",  DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L MN0","",  DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L M IN0","",DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L M I0","",DsMassFitLow, DsMassFitHigh);
// test release background

	f1_DsMix->ReleaseParameter(7);
	f1_DsMix->ReleaseParameter(8);
  h_DsMassData->Fit("f1_DsMix","L M IN0","",DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L M I0","",DsMassFitLow, DsMassFitHigh);


	cout<<"check 2"<<endl;

	f1_DsMix->SetRange(DsMassFitLow,DsMassFitHigh);
	int fitStatus=1;
  int fitIsValid=0;
  TFitResultPtr fitResut;
  double fitPrecision=1.e-8;
  while(fitStatus>0 && fitStatus!=4000){
    TFitter::SetPrecision(fitPrecision);
    fitResut=h_DsMassData->Fit("f1_DsMix","L EMI S0","",DsMassFitLow, DsMassFitHigh);
    fitStatus=fitResut->Status();
    fitIsValid=fitResut->IsValid();
    cout<<"fit precision: "<<TFitter::GetPrecision()<<"   status: "<<fitStatus<<" isValid = "<< fitIsValid <<endl;
    if(fitStatus){
      fitPrecision *= 3;
    }
    if(fitPrecision>1.e-3) {break;}
  }
	fitResut=h_DsMassData->Fit("f1_DsMix","L EMI S0","",DsMassFitLow, DsMassFitHigh);
 cout<<"f1_DsMix->GetParameter(0) = "<<f1_DsMix->GetParameter(0)<<endl;
cout<<"f1_DsSignal->GetParameter(0) = "<<f1_DsSignal->GetParameter(0)<<endl;
  f1_DsSignal->SetParameters(f1_DsMix->GetParameter(0), f1_DsMix->GetParameter(1) , f1_DsMix->GetParameter(2) , f1_DsMix->GetParameter(3) , f1_DsMix->GetParameter(4), f1_DsMix->GetParameter(5));
  f1_DsSignal->SetParError(0,f1_DsMix->GetParError(0));
  f1_DsSignal->SetParError(1,f1_DsMix->GetParError(1));
  f1_DsSignal->SetParError(2,f1_DsMix->GetParError(2));
  f1_DsSignal->SetParError(3,f1_DsMix->GetParError(3));
  f1_DsSignal->SetParError(4,f1_DsMix->GetParError(4));
  f1_DsSignal->SetParError(5,f1_DsMix->GetParError(5));


	cout<<"check 3"<<endl;

  f1_DsSignal->SetFillColor(kOrange-3);
  f1_DsSignal->SetFillStyle(3002);
  f1_DsSignal->SetLineColor(kOrange-3);
  f1_DsSignal->SetLineWidth(2);
  f1_DsSignal->SetLineStyle(2);
  f1_DsSignal->SetRange(DsMassFitLow,DsMassFitHigh); // must setrange before plot or get empty
  double yieldDsSignal= f1_DsSignal->Integral(DsMassFitLow,DsMassFitHigh)/ h_DsMassData->GetBinWidth(1);
  double yieldDsSignalErr = f1_DsSignal->Integral(DsMass-DsMassCandWidth,DsMass+DsMassCandWidth)/h_DsMassData->GetBinWidth(1) * f1_DsSignal->GetParError(0)/f1_DsSignal->GetParameter(0);

  double yieldDsSignalTry = f1_DsSignal->GetParameter(0)/h_DsMassData->GetBinWidth(1);
  double yieldDsSignalTryErr = f1_DsSignal->GetParError(0)/h_DsMassData->GetBinWidth(1);


	cout<<"check 4"<<endl;

  cout<<"parameter 0 = "<<f1_DsSignal->GetParameter(0)<<" BinWidth = "<< h_DsMassData->GetBinWidth(1) <<endl;
  cout<<"yield DsSignal = "<< yieldDsSignal <<" +- "<<yieldDsSignalErr<<endl;
  cout<<"yield DsSignal from Para 0 = "<< yieldDsSignalTry <<" +- "<<yieldDsSignalTryErr<<endl;
  f1_DsBkg->SetParameters(f1_DsMix->GetParameter(6), f1_DsMix->GetParameter(7) ,f1_DsMix->GetParameter(8));
  // f1_DsBkg->SetParameters(f1_DsMix->GetParameter(6), f1_DsMix->GetParameter(7));
  f1_DsBkg->SetParError(0,f1_DsMix->GetParError(6));
  f1_DsBkg->SetParError(1,f1_DsMix->GetParError(7));
  f1_DsBkg->SetParError(2,f1_DsMix->GetParError(8));
//  f1_DsBkg->SetParError(3,f1_DsMix->GetParError(9));

  f1_DsBkg->SetLineColor(1);
  f1_DsBkg->SetLineStyle(2);
  f1_DsBkg->SetRange(DsMassFitLow,DsMassFitHigh);


  double covmat[3][3];
  for(int icov =0; icov<3; icov++){
    for(int jcov =0; jcov<3; jcov++){
      covmat[icov][jcov]=fitResut->CovMatrix(6+icov,6+jcov);
    }
  }

  double *covmatArr=covmat[0];
	

  // double yieldBkg= f1_DsBkg->Integral(DsMassFitLow,DsMassFitHigh)/ h_DsMassData->GetBinWidth(1);
  double yieldBkg= f1_DsBkg->Integral(DsMassFitLow,DsMassFitHigh)/ h_DsMassData->GetBinWidth(1);

  cout<<"yield Bkg = " <<yieldBkg<<endl;

	// bkg sub method

	

	
	double N2Sig=h_DsMassData->Integral(h_DsMassData->FindBin(DsMassMCMean-DsMass2SigWidth),h_DsMassData->FindBin(DsMassMCMean+DsMass2SigWidth));
	double binLow2SigVal=h_DsMassData->GetBinLowEdge(h_DsMassData->FindBin(DsMassMCMean-DsMass2SigWidth));
	double binHigh2SigVal=h_DsMassData->GetBinLowEdge(h_DsMassData->FindBin(DsMassMCMean+DsMass2SigWidth))+h_DsMassData->GetBinWidth(1);
	// double Nbkg2Sig=f1_DsBkg->Integral(DsMassMCMean-DsMass2SigWidth,DsMassMCMean+DsMass2SigWidth)/ h_DsMassData->GetBinWidth(1); // wrong, the edge need to match with th1 bin
	double Nbkg2Sig=f1_DsBkg->Integral(binLow2SigVal,binHigh2SigVal)/ h_DsMassData->GetBinWidth(1);
	double Nbkg2SigError=f1_DsBkg->IntegralError(binLow2SigVal,binHigh2SigVal,f1_DsBkg->GetParameters(), covmatArr )/ h_DsMassData->GetBinWidth(1);
	double NDsbkgSub=N2Sig-Nbkg2Sig;
	double NDsbkgSubError=sqrt(N2Sig+Nbkg2SigError*Nbkg2SigError);
	cout<<"N2Sig = "<<N2Sig<<" +- "<<sqrt(N2Sig)<<endl;
	cout<<"Nbkg2Sig = "<<Nbkg2Sig<<" +- "<<Nbkg2SigError<<endl;
	cout<<"NDsbkgSub = "<<NDsbkgSub<<" +- "<<NDsbkgSubError<<endl;




// #<{(|
	TCanvas *c_DataFit=new TCanvas("c_DataFit","c_DataFit");
	c_DataFit->cd();
	h_DsMassData->Draw();
	f1_DsMix->Draw("same");
	f1_DsBkg->Draw("same");
	f1_DsSignal->Draw("same");

  shiftY=0;
  TLatex *tl_binfitData =new TLatex();
//  tl_binfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s_Data_binfit",str_PbPb.Data()));  shiftY-=oneshift;
  tl_binfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s D_{S} %.0f < pt < %.0f ",s_ppPbPb.Data(),DptLowGL,DptHighGL)); shiftY-=oneshift;
  tl_binfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Fit Yield = %.1f #pm %.1f ",yieldDsSignal,yieldDsSignalErr ));  shiftY-=oneshift;
  tl_binfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("BkgSub Yield = %.1f #pm %.1f ",NDsbkgSub,NDsbkgSubError ));  shiftY-=oneshift;
  // tl_binfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Float width = %.3f",f1_DsMix->GetParameter(5))); shiftY-=oneshift;

  if(s_ExtraName!="")
  {
  tl_binfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,s_ExtraName); shiftY-=oneshift;
  }
  // SavePlotDirs(c_binfit_Data[count_c_binfit],Form("%s_Data_binfit_pt%.0fto%.0f_%s",str_PbPb.Data(),DptLow,DptHigh,extraName.Data() ),{"SignalFit",str_PbPb,"binfit"} );

  c_DataFit->SaveAs(Form("./plots/fit/%s_Data_Mkkpi_%s.png",s_type.Data(),s_ExtraName.Data()));


	FitYield=yieldDsSignal;
	FitYieldErr=yieldDsSignalErr;
	BkgSubYield=NDsbkgSub;
	BkgSubYieldErr=NDsbkgSubError;

	delete c_DataFit;
// |)}>#

	return f1_DsMix;
}


