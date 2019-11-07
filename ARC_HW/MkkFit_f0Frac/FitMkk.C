#include "MkkFit_para.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting_simple.h"
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
#include "RooGaussModel.h"
#include "RooFFTConvPdf.h"


using namespace RooFit;


/*
double crystalball_function(double x, double alpha, double n, double sigma, double mean) {
  // evaluate the crystal ball function
  if (sigma < 0.)     return 0.;
  double z = (x - mean)/sigma; 
  // if (alpha < 0) z = -z; 
  // double abs_alpha = std::abs(alpha);

  // double abs_alpha =alpha;
  // double C = n/abs_alpha * 1./(n-1.) * std::exp(-alpha*alpha/2.);
  // double D = std::sqrt(M_PI/2.)*(1.+ROOT::Math::erf(abs_alpha/std::sqrt(2.)));
  // double N = 1./(sigma*(C+D));

  // if (z  > - abs_alpha)
  if (z  <= - alpha)
    // return std::exp(- 0.5 * z * z);
		return 0;
  else if(z > -alpha){

		// return 0;
    //double A = std::pow(n/abs_alpha,n) * std::exp(-0.5*abs_alpha*abs_alpha);
    double nDivAlpha = n/abs(alpha);
    // double AA =  std::exp(-0.5*alpha*alpha);
    // double B = nDivAlpha -abs(alpha);
    // double arg = nDivAlpha/(B-z);
    // return AA * std::pow(arg,n);
		double A=nDivAlpha*exp(-0.5*alpha*alpha);
		// double B=nDivAlpha -abs(alpha);
		double B=nDivAlpha -alpha;
		// return 
		 return A*pow((B-z),-n);

  }
	 else return 0;
}
*/
double crystalball_function(double x, double alpha, double n, double sigma, double mean) {
  // evaluate the crystal ball function
  if (sigma < 0.)     return 0.;
  double z = (x - mean)/sigma; 
  if (alpha < 0) z = -z; 
  double abs_alpha = std::abs(alpha);
  // double C = n/abs_alpha * 1./(n-1.) * std::exp(-alpha*alpha/2.);
  // double D = std::sqrt(M_PI/2.)*(1.+ROOT::Math::erf(abs_alpha/std::sqrt(2.)));
  // double N = 1./(sigma*(C+D));
  if (z  > - abs_alpha)
    return std::exp(- 0.5 * z * z);
  else {
    //double A = std::pow(n/abs_alpha,n) * std::exp(-0.5*abs_alpha*abs_alpha);
    double nDivAlpha = n/abs_alpha;
    double AA =  std::exp(-0.5*abs_alpha*abs_alpha);
    double B = nDivAlpha -abs_alpha;
    double arg = nDivAlpha/(B-z);
    return AA * std::pow(arg,n);
  }
}

double crystalball_function(const double *x, const double *p) {
  // if ((!x) || (!p)) return 0.; // just a precaution
   // [Constant] * ROOT::Math::crystalball_function(x, [Alpha], [N], [Sigma], [Mean])
  return (p[0] * crystalball_function(x[0], p[3], p[4], p[2], p[1]));
}

/*
double doublecrystalball_fun(const double *x,const double *p){
		
	return (p[0]* (crystalball_function(x[0], p[3], p[4], p[2], p[1]) +  p[8]*crystalball_function(x[0], p[6], p[7], p[5], p[1])));
}
*/

double crys_DoubleGaus_fun(double x, double alpha, double n, double sigma, double mean, double sigma2, double fr2,double dataScale) {
  // evaluate the crystal ball function
  if (sigma < 0.)     return 0.;

	sigma=sigma*(1+dataScale);
	sigma2=sigma2*(1+dataScale);
  double z = (x - mean)/sigma;
	double z1=(x - mean)/sigma2; 
  if (alpha < 0) {
		z = -z;
		z1=-z1;
	} 
  double abs_alpha = std::abs(alpha);
  // double C = n/abs_alpha * 1./(n-1.) * std::exp(-alpha*alpha/2.);
  // double D = std::sqrt(M_PI/2.)*(1.+ROOT::Math::erf(abs_alpha/std::sqrt(2.)));
  // double N = 1./(sigma*(C+D));
  if (z  > - abs_alpha)
    return std::exp(- 0.5 * z * z)+fr2*std::exp(- 0.5 * z1 * z1);
  else {
    //double A = std::pow(n/abs_alpha,n) * std::exp(-0.5*abs_alpha*abs_alpha);
    double nDivAlpha = n/abs_alpha;
    double AA =  std::exp(-0.5*abs_alpha*abs_alpha);
    double B = nDivAlpha -abs_alpha;
    double arg = nDivAlpha/(B-z);
    double arg1 = nDivAlpha/(B-z1);
    // return AA * (std::pow(arg,n)+fr2*std::pow(arg1,n) );
    return AA * std::pow(arg,n)+ fr2*std::exp(- 0.5 * z1 * z1) ;
  }
}

// Double_t bwfunf0(Double_t* x, Double_t* par)
Double_t bwfunf0(Double_t x, Double_t par0, Double_t par1, Double_t par2)
{
  Double_t arg1 = 14.0/22.0; // 2 over pi
  Double_t arg2 = par1*par1*par2*par2; //Gamma=par[1]  M=par[2]
  Double_t arg3 = ((x*x) - (par2*par2))*((x*x) - (par2*par2));
  Double_t arg4 = x*x*x*x*((par1*par1)/(par2*par2));
  return par0*arg1*arg2/(arg3 + arg4);
}



double crys_DoubleGaus_fun(const double *x, const double *p){

	return (p[0] * crys_DoubleGaus_fun(x[0], p[3], p[4], p[2], p[1],p[5],p[6],p[7]));

}



double crys_DoubleGaus_2ndCheb_fun(const double *x, const double *p){

	return (p[0] * crys_DoubleGaus_fun(x[0], p[3], p[4], p[2], p[1],p[5],p[6],p[7]) +p[8]*(1+p[9]*x[0]+p[10]*(2*x[0]*x[0]-1)+p[11]*(4*x[0]*x[0]*x[0]-3*x[0]) )  );

}


double crys_DoubleGaus_bwf0_fun(const double *x, const double *p){

	return (p[0] * crys_DoubleGaus_fun(x[0], p[3], p[4], p[2], p[1],p[5],p[6],p[7]) + bwfunf0(x[0],p[8],p[9],p[10])  );

}


/*
double crys_BW_fun(double x, double alpha, double n, double sigma, double mean, double gamma,double fr2,double dataScale) {
  // evaluate the crystal ball function


  if (sigma < 0.)     return 0.;

	sigma=sigma*1+(dataScale);
	// sigma2=sigma*2+(dataScale);

  Double_t arg1 = 14.0/22.0; // 2 over pi
  Double_t arg2 = gamma*gamma*mean*mean; //Gamma=par[1]  M=par[2]
  Double_t arg3 = ((x*x) - (mean*mean))*((x*x) - (mean*mean));
  Double_t arg4 = x*x*x*x*((gamma*gamma)/(gamma*gamma));
  double bwN=arg1*arg2/(arg3 + arg4);


  double z = (x - mean)/sigma;
	// double z1=(x - mean)/sigma2; 
  if (alpha < 0) {
		z = -z;
		// z1=-z1;
	} 
  double abs_alpha = std::abs(alpha);
  // double C = n/abs_alpha * 1./(n-1.) * std::exp(-alpha*alpha/2.);
  // double D = std::sqrt(M_PI/2.)*(1.+ROOT::Math::erf(abs_alpha/std::sqrt(2.)));
  // double N = 1./(sigma*(C+D));
  if (z  > - abs_alpha){
    // return std::exp(- 0.5 * z * z)+fr2*std::exp(- 0.5 * z1 * z1);
    return std::exp(- 0.5 * z * z)+fr2*bwN;
	 }
  else {
    //double A = std::pow(n/abs_alpha,n) * std::exp(-0.5*abs_alpha*abs_alpha);
    double nDivAlpha = n/abs_alpha;
    double AA =  std::exp(-0.5*abs_alpha*abs_alpha);
    double B = nDivAlpha -abs_alpha;
    double arg = nDivAlpha/(B-z);
    // return AA * (std::pow(arg,n)+fr2*std::pow(arg1,n) );
    return AA * std::pow(arg,n)+ fr2*bwN ;
  }
}




double crys_BW_fun(const double *x, const double *p){

  return (p[0] * crys_DoubleGaus_fun(x[0], p[3], p[4], p[2], p[1],p[5],p[6],p[7]));
      
}     
*/


 Double_t mybwfun(Double_t* x, Double_t* par)
{
  Double_t arg1 = 14.0/22.0; // 2 over pi
  Double_t arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma=par[1]  M=par[2]
  Double_t arg3 = ((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
  Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[1]*par[1])/(par[2]*par[2]));
  return par[0]*arg1*arg2/(arg3 + arg4);
}


TF1 *fitPhiMC(TH1D* h_PhiMC);
TF1 *fitPhiMass(TH1D *hData, TF1 *f1MC,double &NPhiInCut,double &NPhiInCutErr,double &NBkgInCut,double &NBkgInCutErr,  TString s_ExtraName="",int fixSigShape=0, TF1 *f1_default=NULL);

TF1 *fitPhiMass_cry2Gaus(TH1D *hData, TF1 *f1MC,double &NPhiInCut,double &NPhiInCutErr,double &NBkgInCut,double &NBkgInCutErr,  TString s_ExtraName="",int fixSigShape=0, TF1 *f1_default=NULL);


  // double PhiMassCandWidth=0.05;
  // double PhiMassCandFitMeanWidth=0.02;
  double shiftY=0.05;
  double shiftX=0.32;
  double oneshift=0.070;

  double DptLowGL=2;
  double DptHighGL=40;
   TString s_ppPbPb="pp";
  TString s_type="";


	// double PhiMass=1.01951;
	double PhiMassCandWidth=0.1;
  double PhiMassCandFitMeanWidth=0.02;
	// double PhiMassFitLow=1.01951-0.02;
	// double PhiMassFitHigh=1.01951+0.02;

	double PhiMassFitLow=0.99;
	double PhiMassFitHigh=1.04;




int FitMkk(int isPbPb=3, double DptLow=6, double DptHigh=40, int useBkgSubData=0, double DVtxPCut=0.25, double DdlsCut=4.5){


	InitStyle();
	setTDRStyle();	

  double tex_upperY=0.95;

  // TLatex* texCmsPre = new TLatex(0.15,tex_upperY, "#scale[1.25]{CMS} Preliminary");
  TLatex* texCmsPre = new TLatex(0.15,tex_upperY, "#scale[1.25]{#bf{CMS}}");
  texCmsPre->SetNDC();
  texCmsPre->SetTextAlign(12);
  texCmsPre->SetTextSize(0.045);
  texCmsPre->SetTextFont(42);

  TLatex* texCms2 = new TLatex(0.20,tex_upperY-0.08, "#scale[1.25]{CMS}");
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


  // TLatex* texColPbPb = new TLatex(0.95,tex_upperY, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  TLatex* texColPbPb = new TLatex(0.95,tex_upperY, "44 #mub^{-1} (5.02 TeV PbPb)");
  texColPbPb->SetNDC();
  texColPbPb->SetTextAlign(32);
  texColPbPb->SetTextSize(0.045);
  texColPbPb->SetTextFont(42);

  // TLatex* texColpp = new TLatex(0.95,tex_upperY, "pp #sqrt{s} = 5.02 TeV");
  TLatex* texColpp = new TLatex(0.95,tex_upperY, "38 nb^{-1} (5.02 TeV pp )");
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


	gStyle->SetOptStat(0);

   DptLowGL=DptLow;
  DptHighGL=DptHighGL;

	double phiFitLow=0.99;
	 double phiFitHigh=1.04;

  gSystem->Exec("mkdir -p plots/fitMkk");
	
	if(isPbPb==3)
	{
		s_ppPbPb="PbPb";

	}

  TString s_fMC_Prompt_Phi=Form("/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_Mkk/%s_FitTree_MC.root",s_ppPbPb.Data());
	TFile *f_MCP_Phi_tree=TFile::Open(s_fMC_Prompt_Phi.Data(),"READ");
	TTree *ntDs_MCP_Phi=(TTree*)f_MCP_Phi_tree->Get("t_fit");
	// ntDs_MCP_Phi->Print();

  // DVtxPCut=0.2;
	// DdlsCut=4.5;

	TString DataCut=Form("(Dy<1 && Dy>-1 && Dmass>%f && Dmass<%f && Dpt>%f && Dpt<%f && Dalpha<0.12 && Dchi2cl>%f && Ddls >%f && DtktkResmass>%f && DtktkResmass<%f )",DsMassMCMean-DsMass2SigWidth , DsMassMCMean+DsMass2SigWidth, DptLow, DptHigh, DVtxPCut, DdlsCut, bins_mkk[0], bins_mkk[nbin_mkk] );

	cout<<"DataCut = "<<DataCut<<endl;


	RooRealVar Dmass("Dmass","Dmass",1.7,2.3);
	RooRealVar TotalWeight("TotalWeight","TotalWeight",0);
	RooRealVar Dpt("Dpt","Dpt",0);
	RooRealVar Dy("Dy","Dy",0);
	RooRealVar Dalpha("Dalpha","Dalpha",0);
	RooRealVar Dchi2cl("Dchi2cl","Dchi2cl",0);
	RooRealVar Ddls("Ddls","Ddls",0);
	RooRealVar DtktkResmass("DtktkResmass","DtktkResmass",bins_mkk[0],bins_mkk[nbin_mkk]);
	
	RooDataSet RooDs_MCAll("RooDs_MCAll","RooDs_MCAll",RooArgSet(Dmass,Dy,TotalWeight,Dpt,Dalpha,Dchi2cl,Ddls,DtktkResmass),WeightVar(TotalWeight),Import(*ntDs_MCP_Phi) );
	RooDataSet RooDs_MC_cut=*(RooDataSet*)RooDs_MCAll.reduce(RooArgSet(DtktkResmass),DataCut.Data());

	RooPlot *massframe_MC=DtktkResmass.frame();
	// RooDs_MCAll.plotOn(massframe_MC,LineColor(1));
	RooDs_MC_cut.plotOn(massframe_MC);
	// massframe_MC->Draw();

	RooDataHist *hist=RooDs_MC_cut.binnedClone();
	RooPlot *frame1=DtktkResmass.frame(Bins(20));
	RooDs_MC_cut.plotOn(frame1);
	// hist->plotOn(frame1);

	RooHistPdf histpdf("histpdf","histpdf",DtktkResmass,*hist,0);
	histpdf.plotOn(frame1);


	RooRealVar mean("mean","mean", 0, -1000,  1000) ;
	mean.setConstant(kTRUE);
  RooRealVar sigma("sigma","sigma", 0.2, 0, 1) ;
	RooGaussModel resol("resol", "resol", DtktkResmass, mean, sigma);
  //RooHistPdf histpdf("histpdf","histpdf",x, hMC,0) ;
	
	DtktkResmass.setBins(10000, "cache");
	RooFFTConvPdf pdf("pdf", "pdf", DtktkResmass,  histpdf, resol);

	// pdf.plotOn(frame1,LineColor(2));    	
	// pdf.fitTo(*hist);
	// pdf.plotOn(frame1,LineColor(2));    	

	// frame1->Draw();



//  ntDs_MCP_Phi->Project("h_mkk_MCP_Phi_fineBin","DtktkResmass",Form("PbPbVzWeight*weight*DgenptSampleWeight*RecoD0DataWeight*(Dmass>%f && Dmass<%f && Dpt>%f && Dpt<%f && Dalpha<0.12 && Dchi2cl>%f && Ddls >%f && DtktkResmass>%f && DtktkResmass<%f && DsGen==23333)",DsMassMCMean-DsMass2SigWidth , DsMassMCMean+DsMass2SigWidth, DptLow, DptHigh, DVtxPCut, DdlsCut, bins_mkk[0], bins_mkk[nbin_mkk] ) );




  //  TString s_ppPbPb="pp";

  TString s_fFit=Form("./his/%s_Dpt%.0fto%.0f.root",s_ppPbPb.Data(),DptLow,DptHigh);
  TString s_fMC=Form("./his/%s_Dpt%.0fto%.0f_mkkpiFit.root",s_ppPbPb.Data(),DptLow,DptHigh);

  s_type=Form("%s_Dpt%.0fto%.0f",s_ppPbPb.Data(),DptLow,DptHigh);

  TFile *fFit=TFile::Open(s_fFit.Data(),"READ");
  TFile *fMC=TFile::Open(s_fMC.Data(),"READ");
	
	TFile *fout=TFile::Open(Form("./his/%s_mkkFitResult.root",s_type.Data()),"RECREATE");


   TH1D *h_mkk_MCP_Phi          =(TH1D*)fFit->Get("h_mkk_MCP_Phi");
  TH1D *h_mkk_MCP_f0           =(TH1D*)fFit->Get("h_mkk_MCP_f0");
  TH1D *h_mkk_MCP_Phi_fineBin  =(TH1D*)fFit->Get("h_mkk_MCP_Phi_fineBin");
  TH1D *h_mkk_MCP_f0_fineBin   =(TH1D*)fFit->Get("h_mkk_MCP_f0_fineBin");


	TH1D *h_NDs_mkkbin_Fit   =(TH1D*)fMC->Get("h_NDs_mkkbin_Fit");
	TH1D *h_NDs_mkkbin_Fit_roofit   =(TH1D*)fMC->Get("h_NDs_mkkbin_Fit_roofit");
	TH1D *h_NDs_mkkbin_BkgSub_roofit   =(TH1D*)fMC->Get("h_NDs_mkkbin_BkgSub_roofit");
	TH1D *h_NDs_mkkbin_BkgSub   =(TH1D*)fMC->Get("h_NDs_mkkbin_BkgSub");

	TH1D *h_NDs_mkkbin=h_NDs_mkkbin_Fit_roofit;
	TString s_BkgSub="";
	if(useBkgSubData){
			h_NDs_mkkbin=h_NDs_mkkbin_BkgSub_roofit;
			s_BkgSub="BkgSub";
	}


	// test MC/Data ratio
	// TH1D *h_mkk_MCP_Phi_overData= (TH1D*)h_mkk_MCP_Phi->Clone("h_mkk_MCP_Phi_overData");
	// h_mkk_MCP_Phi_overData->Divide(h_NDs_mkkbin_Fit);
	// h_mkk_MCP_Phi_overData->Draw();


	/////////////////
	// roofit test //
	/////////////////

	RooRealVar x("x","x",0.99,1.04);
	RooDataHist dh_MC("dh_MC","dh_MC",DtktkResmass,Import(*h_mkk_MCP_Phi_fineBin));
	RooPlot *MC_frame=DtktkResmass.frame(Title("MC"));
	dh_MC.plotOn(MC_frame);
	// dh_MC.plotOn(massframe_MC);

	RooPlot *testframe=new RooPlot("testframe","DtktkResmass",DtktkResmass,1.0,1.04,40);

	double cbMean_MC_pre=1.0196;
	double cbSigma_MC_pre=0.0066;
	double cbN_MC_pre=0.0001;
	double cbAlpha_MC_pre=-0.6;
	double sigma1_MC_pre=0.01;
	double sigma2_MC_pre=0.0024;
	double fr_cb_MC_pre=0.4;
	double fr_gaus1_MC_pre=0.4;

	double  cbSigma2_MC_pre=7.78e-3;
	double  cbSigma3_MC_pre=7.78e-3;
	double  fr_cb1_MC_pre=0.8;
	double  fr_cb2_MC_pre=0.8;


	if(isPbPb==0){
	 cbMean_MC_pre=1.0195;
	 cbSigma_MC_pre=0.002;
	 // cbN_MC_pre=1.2;
	 cbN_MC_pre=1;
	 cbAlpha_MC_pre=-1;
	 sigma1_MC_pre=0.01;
	 // sigma2_MC_pre=0.0024;
	 fr_cb_MC_pre=0.6;
	 fr_gaus1_MC_pre=0.4;
	}

	if(isPbPb==3){
	 cbMean_MC_pre=1.0195;
	 cbSigma_MC_pre=2.4e-3;
	 cbSigma2_MC_pre=7.8e-3;
	 cbSigma3_MC_pre=5.6e-3;
	 // cbN_MC_pre=1.2;
	 cbN_MC_pre=1.0;
	 cbAlpha_MC_pre=-1.7;
	 sigma1_MC_pre=7.4e-3;
	 sigma2_MC_pre=7.4e-5;
	 fr_cb_MC_pre=0.8;
	 fr_cb1_MC_pre=0.8;
	 fr_cb2_MC_pre=0.8;
	 fr_gaus1_MC_pre=0.9;
	}

// unbin fit


	RooRealVar cbMean_MC("cbMean_MC","cbMean_MC",cbMean_MC_pre,1.015,1.023);
	RooRealVar GSMean_MC("GsMean_MC","GSMean_MC",cbMean_MC_pre,1.015,1.023);
	RooRealVar cbMean2_MC("cbMean2_MC","cbMean2_MC",cbMean_MC_pre,1.015,1.023);
	RooRealVar cbSigma_MC("cbSigma_MC","cbSigma_MC",cbSigma_MC_pre,0,0.5);
	RooRealVar cbSigma2_MC("cbSigma2_MC","cbSigma2_MC",cbSigma2_MC_pre,0,0.5);
	RooRealVar cbSigma3_MC("cbSigma3_MC","cbSigma3_MC",cbSigma3_MC_pre,0,0.5);
	// RooRealVar cbNsig;
	RooRealVar cbN_MC("cbN_MC","cbN_MC",cbN_MC_pre,-1e8,1e8);
	RooRealVar cbN2_MC("cbN2_MC","cbN2_MC",cbN_MC_pre,-1e8,1e8);
	RooRealVar cbN3_MC("cbN3_MC","cbN3_MC",cbN_MC_pre,-1e8,1e8);
	RooRealVar cbAlpha_MC("cbAlpha_MC","cbAlpha_MC",cbAlpha_MC_pre,-1e6,1e6);
	RooRealVar cbAlpha2_MC("cbAlpha2_MC","cbAlpha2_MC",cbAlpha_MC_pre,-1e6,1e6);
	RooRealVar cbAlpha3_MC("cbAlpha3_MC","cbAlpha3_MC",cbAlpha_MC_pre,-1e6,1e6);

	RooCBShape cball_MC("cball_MC","crystal ball",DtktkResmass,cbMean_MC,cbSigma_MC,cbAlpha_MC,cbN_MC);	

	// RooCBShape cball2_MC("cball2_MC","crystal ball2",DtktkResmass,cbMean_MC,cbSigma2_MC,cbAlpha2_MC,cbN2_MC);	
	RooCBShape cball2_MC("cball2_MC","crystal ball2",DtktkResmass,cbMean_MC,cbSigma2_MC,cbAlpha_MC,cbN_MC);	
	RooCBShape cball3_MC("cball3_MC","crystal ball3",DtktkResmass,cbMean_MC,cbSigma3_MC,cbAlpha_MC,cbN_MC);	

	RooRealVar sigma1_MC("sigma1_MC","sigma1_MC",sigma1_MC_pre,0,0.5);
	RooRealVar sigma2_MC("sigma2_MC","sigma2_MC",sigma2_MC_pre,0,0.5);

	RooRealVar fr_cb_MC("fr_cb_MC","fr_cb_MC",fr_cb_MC_pre,0,1);
	RooRealVar fr_cb1_MC("fr_cb1_MC","fr_cb1_MC",fr_cb1_MC_pre,0,1);
	RooRealVar fr_cb2_MC("fr_cb2_MC","fr_cb2_MC",fr_cb2_MC_pre,0,1);
	RooRealVar fr_gaus1_MC("fr_gaus1_MC","fr_gaus1_MC",fr_gaus1_MC_pre,0,1);
	RooRealVar fr_gaus2_MC("fr_gaus2_MC","fr_gaus2_MC",fr_gaus1_MC_pre,0,1);

	// RooGaussian Gaus1_MC("Gaus1_MC","Gaus1_MC",DtktkResmass,GSMean_MC,sigma1_MC);
	RooGaussian Gaus1_MC("Gaus1_MC","Gaus1_MC",DtktkResmass,cbMean_MC,sigma1_MC);
	// RooGaussian Gaus1_MC("Gaus1_MC","Gaus1_MC",x,cbMean_MC,cbSigma_MC);
	RooGaussian Gaus2_MC("Gaus2_MC","Gaus2_MC",DtktkResmass,cbMean_MC,sigma2_MC);
	// RooGaussian Gaus3_MC("Gaus3_MC","Gaus3_MC",DtktkResmass,cbMean_MC,sigma2_MC);

	RooRealVar Cheb1_MC("Cheb1_MC","Cheb1_MC",0,-1e6,1e6);
	RooChebychev Cheb1_MCPdf("Cheb1_MCPdf","Cheb1_MCPdf",DtktkResmass,RooArgList(Cheb1_MC));
	RooRealVar fr_bg_MC("fr_bg_MC","fr_bg_MC",0.1,0,1);



	// RooRealVar NumCB1("NumCB1","NumCB1",1000,-1e15,1e15);
	// RooRealVar NumCB2("NumCB2","NumCB2",1000,-1e15,1e15);
	// RooRealVar NumGS1("NumGS1","NumGS1",1000,-1e15,1e15);
	// RooRealVar NumGS2("NumGS2","NumGS2",1000,-1e15,1e15);

	// RooRealVar NumCB1("NumCB1","NumCB1",1000,0,1e15);
	// RooRealVar NumCB2("NumCB2","NumCB2",1000,0,1e15);
	// RooRealVar NumGS1("NumGS1","NumGS1",1000,0,1e15);
	// RooRealVar NumGS2("NumGS2","NumGS2",1000,0,1e15);


	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,Gaus1_MC),RooArgList(NumCB1,NumGS1));
	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,Gaus1_MC,Gaus2_MC),RooArgList(NumCB1,NumGS1,NumGS2));

	RooAddPdf CB1GS1_MC("CB1GS1_MC","CB1GS1_MC",RooArgList(cball_MC,Gaus1_MC),RooArgList(fr_cb_MC));
	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(Cheb1_MCPdf , CB1GS1_MC),RooArgList(fr_bg_MC));

	RooAddPdf nGaus_MC("nGaus_MC","nGaus_MC",RooArgList(Gaus1_MC,Gaus2_MC),RooArgList(fr_gaus1_MC));
	RooAddPdf nCB_MC("nCB_MC","nCB_MC",RooArgList(cball_MC,cball2_MC),RooArgList(fr_cb_MC));
	RooAddPdf CB23_MC("CB23_MC","CB23_MC",RooArgList(cball2_MC,cball3_MC),RooArgList(fr_cb2_MC));
	RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,CB23_MC),RooArgList(fr_cb1_MC)); // this is good for pp

	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,Gaus1_MC),RooArgList(fr_cb_MC));

	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,nGaus_MC),RooArgList(fr_cb_MC,fr_gaus1_MC));

	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,cball2_MC),RooArgList(fr_cb_MC)); // good for PbPb;

	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,cball2_MC,cball3_MC),RooArgList(fr_cb_MC,fr_cb1_MC));

	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,CB23_MC),RooArgList(fr_cb1_MC)); // this is good for pp

	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,cball2_MC,Gaus1_MC),RooArgList(fr_cb_MC,fr_gaus1_MC));
	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(Gaus1_MC,nCB_MC),RooArgList(fr_gaus1_MC));

	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,Gaus1_MC,Gaus2_MC,Gaus3_MC),RooArgList(fr_cb_MC,fr_gaus1_MC,fr_gaus2_MC));

	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,Gaus1_MC,Gaus2_MC),RooArgList(fr_cb_MC,fr_gaus1_MC)); // use this for pp 

	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,Gaus1_MC),RooArgList(fr_cb_MC));

	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,cball2_MC),RooArgList(fr_cb_MC));
	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,cball2_MC,Gaus1_MC),RooArgList(fr_cb_MC,fr_gaus1_MC));

	// cbMean_MC.setConstant(kTRUE);
	// fr_cb_MC.setVal(1);
	// fr_cb_MC.setConstant(kTRUE);
	// fr_gaus1_MC.setVal(1);
	// fr_gaus1_MC.setConstant(kTRUE);

	// SigPdf_MC.fitTo(dh_MC,Minos(kTRUE));
	// SigPdf_MC.fitTo(dh_MC,Hesse(),NumCPU(10));
	// SigPdf_MC.fitTo(dh_MC,Extended());
	// SigPdf_MC.fitTo(dh_MC,Extended());
	// SigPdf_MC.fitTo(dh_MC,Range(1.01,1.03));

	massframe_MC->Draw();

	// SigPdf_MC.fitTo(RooDs_MC_cut,NumCPU(20),Extended(kTRUE),Range(1.0,1.04));
	// SigPdf_MC.fitTo(RooDs_MC_cut,NumCPU(20));

	SigPdf_MC.fitTo(RooDs_MC_cut,NumCPU(30),Hesse(kTRUE));
	SigPdf_MC.fitTo(RooDs_MC_cut,NumCPU(30),Hesse(kTRUE));
	// SigPdf_MC.fitTo(RooDs_MC_cut,NumCPU(30),Hesse(kTRUE));

	// SigPdf_MC.fitTo(RooDs_MC_cut,NumCPU(20),Extended(kTRUE));

//	cball_MC.fitTo(RooDs_MC_cut);
//	cball_MC.plotOn(MC_frame)

	SigPdf_MC.plotOn(massframe_MC,LineColor(2));
	// SigPdf_MC.plotOn(massframe_MC,Components(cball_MC),LineColor(4));
	// SigPdf_MC.plotOn(massframe_MC,Components(cball2_MC),LineColor(kYellow+2));
	// SigPdf_MC.plotOn(massframe_MC,Components(cball3_MC),LineColor(kGreen+2));
	// SigPdf_MC.plotOn(massframe_MC,Components(Gaus1_MC),LineColor(kGreen+2));
	// SigPdf_MC.plotOn(massframe_MC,Components(Gaus2_MC),LineColor(kYellow+2));

	TCanvas *c_roofit_MC=new TCanvas("c_roofit_MC","c_roofit_MC");
	c_roofit_MC->cd();

	SetCanvas(c_roofit_MC);

	// MC_frame->GetXaxis()->SetTitle("m_{KK} (GeV/c^{2})");
	massframe_MC->Draw();
/*
	TCanvas *c_test=new TCanvas("c_test","c_test");
	c_test->cd();

	SigPdf_MC.plotOn(testframe,Components(cball3_MC),LineColor(kGreen+2));
	testframe->Draw();
*/
	massframe_MC->Print("V");

/*
  RooGoF goftest_MC(massframe_MC->getHist("h_dh_data"),Data_frame->getCurve(s_pdf));
  cout<<"check 0"<<endl;
  goftest.setRange(x.getMin(),x.getMax());
  goftest.setRebin(5,true); // better rebin to make sure that all bins have >=5 expected events  
  // goftest.setRebin(10,false); // better rebin to make sure that all bins have >=5 expected events  
  int ndf=0;
  int d_ndf=fitresult->floatParsFinal().getSize();
  double pvalue_BC=0.0;
  double testStat_BC=0.0;
  goftest.BCChi2Test(pvalue_BC,testStat_BC,ndf,d_ndf);
  // file_pval << "BC: " << pvalue_BC << ", " << testStat_BC << endl;
  cout<<"BC: " << pvalue_BC << ", " << testStat_BC << endl;

  double pvalue_PC=0.0;
  double testStat_PC=0.0;
  goftest.PearsonChi2Test(pvalue_PC,testStat_PC,ndf,d_ndf);
  // file_pval << "BC: " << pvalue_BC << ", " << testStat_BC << endl;
  cout<<"PC: " << pvalue_PC << ", " << testStat_PC << endl;

  double pvalue_NC=0.0;
  double testStat_NC=0.0;
  goftest.NeymanChi2Test(pvalue_NC,testStat_NC,ndf,d_ndf);
  // file_pval << "BC: " << pvalue_BC << ", " << testStat_BC << endl;
  cout<<"NC: " << pvalue_NC << ", " << testStat_NC << endl;
 
 
  double pvalue_RooFit=0.0;
  double testStat_RooFit=0.0;
  goftest.RooFitChi2Test(pvalue_RooFit,testStat_RooFit,ndf,d_ndf);
  // file_pval << "RooFitChi2: " << pvalue_RooFit << ", " << testStat_RooFit << endl;
  cout<<"RooFitChi2: " << pvalue_RooFit << ", " << testStat_RooFit << endl;
*/













// bin fit

/*
	RooRealVar cbMean_MC("cbMean_MC","cbMean_MC",cbMean_MC_pre,1.015,1.023);
	RooRealVar cbSigma_MC("cbSigma_MC","cbSigma_MC",cbSigma_MC_pre,-0.1,0.1);
	RooRealVar cbSigma2_MC("cbSigma2_MC","cbSigma2_MC",cbSigma_MC_pre,-0.1,0.1);
	// RooRealVar cbNsig;
	RooRealVar cbN_MC("cbN_MC","cbN_MC",cbN_MC_pre,-1e8,1e8);
	RooRealVar cbN2_MC("cbN2_MC","cbN2_MC",cbN_MC_pre,-1e2,1e8);
	RooRealVar cbAlpha_MC("cbAlpha_MC","cbAlpha_MC",cbAlpha_MC_pre,-1e6,1e6);
	RooRealVar cbAlpha2_MC("cbAlpha2_MC","cbAlpha2_MC",cbAlpha_MC_pre,-1e6,-0.0000000001);

	RooCBShape cball_MC("cball_MC","crystal ball",x,cbMean_MC,cbSigma_MC,cbAlpha_MC,cbN_MC);	

	// RooCBShape cball2_MC("cball2_MC","crystal ball2",x,cbMean_MC,cbSigma2_MC,cbAlpha2_MC,cbN2_MC);	
	RooCBShape cball2_MC("cball2_MC","crystal ball2",x,cbMean_MC,cbSigma2_MC,cbAlpha_MC,cbN_MC);	

	RooRealVar sigma1_MC("sigma1_MC","sigma1_MC",sigma1_MC_pre,-0.1,0.1);
	RooRealVar sigma2_MC("sigma2_MC","sigma2_MC",sigma2_MC_pre,-0.1,0.1);

	RooRealVar fr_cb_MC("fr_cb_MC","fr_cb_MC",fr_cb_MC_pre,-2,2);
	RooRealVar fr_gaus1_MC("fr_gaus1_MC","fr_gaus1_MC",fr_gaus1_MC_pre,-2,2);
	RooRealVar fr_gaus2_MC("fr_gaus2_MC","fr_gaus2_MC",fr_gaus1_MC_pre,-2,2);

	RooGaussian Gaus1_MC("Gaus1_MC","Gaus1_MC",x,cbMean_MC,sigma1_MC);
	// RooGaussian Gaus1_MC("Gaus1_MC","Gaus1_MC",x,cbMean_MC,cbSigma_MC);
	RooGaussian Gaus2_MC("Gaus2_MC","Gaus2_MC",x,cbMean_MC,sigma2_MC);
	RooGaussian Gaus3_MC("Gaus3_MC","Gaus3_MC",x,cbMean_MC,sigma2_MC);

	RooAddPdf nGaus_MC("nGaus_MC","nGaus_MC",RooArgList(Gaus1_MC,Gaus2_MC),RooArgList(fr_gaus1_MC));

	RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,nGaus_MC),RooArgList(fr_cb_MC,fr_gaus1_MC));
	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,Gaus1_MC,Gaus2_MC,Gaus3_MC),RooArgList(fr_cb_MC,fr_gaus1_MC,fr_gaus2_MC));
	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,Gaus1_MC,Gaus2_MC),RooArgList(fr_cb_MC,fr_gaus1_MC));

	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,Gaus1_MC),RooArgList(fr_cb_MC));

	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,cball2_MC),RooArgList(fr_cb_MC));
	// RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(cball_MC,cball2_MC,Gaus1_MC),RooArgList(fr_cb_MC,fr_gaus1_MC));

	// cbMean_MC.setConstant(kTRUE);
	// fr_cb_MC.setVal(1);
	// fr_cb_MC.setConstant(kTRUE);
	// fr_gaus1_MC.setVal(1);
	// fr_gaus1_MC.setConstant(kTRUE);

	// SigPdf_MC.fitTo(dh_MC,Minos(kTRUE));
	// SigPdf_MC.fitTo(dh_MC,Hesse());
	// SigPdf_MC.fitTo(dh_MC,Extended());
	// SigPdf_MC.fitTo(dh_MC,Extended());
	// SigPdf_MC.fitTo(dh_MC,Extended());
	// SigPdf_MC.fitTo(dh_MC,Extended());
	// SigPdf_MC.fitTo(dh_MC,Extended());
	// SigPdf_MC.fitTo(dh_MC,Range(1.01,1.03));
	// SigPdf_MC.fitTo(dh_MC);


	SigPdf_MC.plotOn(MC_frame,LineColor(2));
	SigPdf_MC.plotOn(MC_frame,Components(cball_MC),LineColor(4));
	SigPdf_MC.plotOn(MC_frame,Components(Gaus1_MC),LineColor(kGreen+2));
	SigPdf_MC.plotOn(MC_frame,Components(Gaus2_MC),LineColor(kYellow+2));

	TCanvas *c_roofit_MC=new TCanvas("c_roofit_MC","c_roofit_MC");
	c_roofit_MC->cd();
	MC_frame->GetXaxis()->SetTitle("m_{KK} (GeV/c^{2})");
	MC_frame->Draw();

*/

	massframe_MC->GetXaxis()->SetTitle("m_{KK} (GeV/c^{2})");
	massframe_MC->GetXaxis()->CenterTitle();
	massframe_MC->GetYaxis()->CenterTitle();
	massframe_MC->SetTitle("");

	TLatex *tltx =new TLatex();
	shiftX=-0.00;
	shiftY=0.05;
	tltx->SetTextSize(0.042);
//  tl_binfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s_Data_binfit",str_PbPb.Data()));  shiftY-=oneshift;
  tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s MC ",s_ppPbPb.Data())); shiftY-=oneshift;
  tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%.0f < D_{S} p_{T} < %.0f GeV/c",DptLowGL,DptHighGL)); shiftY-=oneshift;
  // tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Phi Yield = %.1f #pm %.1f ",N_Sig_2sig,N_Sig_2sig*NumSig.getError()/NumSig.getValV() ));  shiftY-=oneshift;
  // tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Bkg Yield = %.1f #pm %.1f ",N_Bkg_2sig,N_Bkg_2sig*NumBkg.getError()/NumBkg.getValV() ));  shiftY-=oneshift;
  // tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Phi Ratio = %.3f #pm %.3f ",NPhiRatio,NPhiRatioErr ));  shiftY-=oneshift;
  // tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("BC pvalue = %.3f ",pvalue_BC ));  shiftY-=oneshift;
  // tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("RooFit pvalue = %.3f ",pvalue_RooFit ));  shiftY-=oneshift;
  // tltx->DrawLatexNDC(textpo+shiftX.08,textposy+shiftY,Form("BkgSub Yield = %.1f #pm %.1f ",NPhibkgSub,NPhibkgSubError ));  shiftY-=oneshift;
  // tltx->DrawLatexNDC(textpo+shiftX.28,textposy+shiftY,Form("Float width = %.3f",f1_PhiMix->GetParameter(5))); shiftY-=oneshift;

  // if(s_ExtraName!="")
  // {
  // tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,s_ExtraName); shiftY-=oneshift;
  // }
  // SavePlotDirs(c_binfit_Data[count_c_binfit],Form("%s_Data_binfit_pt%.0fto%.0f_%s",str_PbPb.Data(),DptLow,DptHigh,extraName.Data() ),{"SignalFit",str_PbPb,"binfit"} );

	
 c_roofit_MC->SaveAs(Form("./plots/fitMkk/%s_MC_Mkk_roofit.png",s_type.Data()));
 SavePlotDirs(c_roofit_MC,Form("%s_MC_Mkk_roofit%s",s_type.Data(),s_BkgSub.Data() ),{"PhiRatio"},"gcf");

	// return 1;



	// c_roofit_MC->SaveAs()


	double cbMean_MC_V=cbMean_MC.getValV();
	double cbSigma_MC_V=cbSigma_MC.getValV();
	double cbSigma2_MC_V=cbSigma2_MC.getValV();
	double cbSigma3_MC_V=cbSigma3_MC.getValV();
	double cbN_MC_V=cbN_MC.getValV();
	double cbAlpha_MC_V=cbAlpha_MC.getValV();
	double sigma1_MC_V=sigma1_MC.getValV();
	double sigma2_MC_V=sigma2_MC.getValV();
	double fr_cb_MC_V=fr_cb_MC.getValV();
	double fr_gaus1_MC_V=fr_gaus1_MC.getValV();
	double fr_cb1_MC_V=fr_cb1_MC.getValV();
	double fr_cb2_MC_V=fr_cb2_MC.getValV();

	RooRealVar cbMean("cbMean","cbMean",cbMean_MC_V,1.015,1.023);
	RooRealVar cbSigma("cbSigma","cbSigma",cbSigma_MC_V,0.0001,0.1);
	RooRealVar cbSigma2("cbSigma2","cbSigma2",cbSigma2_MC_V,0,0.5);
	RooRealVar cbSigma3("cbSigma3","cbSigma3",cbSigma3_MC_V,0,0.5);


	// RooRealVar cbNsig;
	RooRealVar cbN("cbN","cbN",cbN_MC_V,0.0001,1e8);
	RooRealVar cbAlpha("cbAlpha","cbAlpha",cbAlpha_MC_V,-1e6,-0.0000000001);


	RooRealVar sigma1("sigma1","sigma1",sigma1_MC_V,0.00001,0.1);
	RooRealVar sigma2("sigma2","sigma2",sigma2_MC_V,0.00001,0.1);

	RooRealVar FloatScale("FloatScale","FloatScale",0,-0.6,0.6);
	RooFormulaVar scale_cbSigma("scale_cbSigma","scale_cbSigma","cbSigma*(1+FloatScale)",RooArgSet(cbSigma,FloatScale));
	RooFormulaVar scale_cbSigma2("scale_cbSigma2","scale_cbSigma2","cbSigma2*(1+FloatScale)",RooArgSet(cbSigma2,FloatScale));
	RooFormulaVar scale_cbSigma3("scale_cbSigma3","scale_cbSigma3","cbSigma3*(1+FloatScale)",RooArgSet(cbSigma3,FloatScale));
	RooFormulaVar scale_Sigma1("scale_Sigma1","scale_Sigma1","sigma1*(1+FloatScale)",RooArgSet(sigma1,FloatScale));
	RooFormulaVar scale_Sigma2("scale_Sigma2","scale_Sigma2","sigma2*(1+FloatScale)",RooArgSet(sigma2,FloatScale));


	RooRealVar fr_cb("fr_cb","fr_cb",fr_cb_MC_V,0,1);
	RooRealVar fr_gaus1("fr_gaus1","fr_gaus1",fr_gaus1_MC_V,0,1);
	RooRealVar fr_cb1("fr_cb1","fr_cb1",fr_cb1_MC_V,0,1);
	RooRealVar fr_cb2("fr_cb2","fr_cb2",fr_cb2_MC_V,0,1);

	RooCBShape cball("cball","crystal ball",x,cbMean,scale_cbSigma,cbAlpha,cbN);	
	RooCBShape cball2("cball2","crystal ball2",x,cbMean,scale_cbSigma2,cbAlpha,cbN);	
	RooCBShape cball3("cball3","crystal ball3",x,cbMean,scale_cbSigma3,cbAlpha,cbN);	
	RooGaussian Gaus1("Gaus1","Gaus1",x,cbMean,scale_Sigma1);
	RooGaussian Gaus2("Gaus2","Gaus2",x,cbMean,scale_Sigma2);

	// RooAddPdf SigPdf("SigPdf","SigPdf",RooArgList(cball,Gaus1,Gaus2),RooArgList(fr_cb,fr_gaus1));
	// RooAddPdf SigPdf("SigPdf","SigPdf",RooArgList(cball,Gaus1),RooArgList(fr_cb));
	RooAddPdf CB23("CB23","CB23",RooArgList(cball2,cball3),RooArgList(fr_cb2));
	RooAddPdf SigPdf("SigPdf","SigPdf",RooArgList(cball,CB23),RooArgList(fr_cb1)); // this is good for pp
	// RooAddPdf SigPdf("SigPdf","SigPdf",RooArgList(cball,cball2,cball3),RooArgList(fr_cb1,fr_cb2)); // this is good for pp

	RooRealVar Cheb1("Cheb1","Cheb1",0,-1,1);
	RooRealVar Cheb2("Cheb2","Cheb2",0,-1e6,1e6);

	RooChebychev BkgPdf("BkgPdf","BkgPdf",x,RooArgList(Cheb1,Cheb2));
	// RooChebychev BkgPdf("BkgPdf","BkgPdf",x,RooArgList(Cheb1));


	Cheb2.setConstant(kTRUE);

	RooRealVar bwMean("bwMean","bwMean",0.99,0.97,1.01);
	RooRealVar bwWidth("bwWidth","bwWidth",0.02,0.001,0.1);

	// RooBreitWigner BkgPdf("BkgPdf","BkgPdf",x,bwMean,bwWidth);	

	RooRealVar NumSig("NumSig","Number of Signal",1000,-1e2,4e9);
	RooRealVar NumBkg("NumBkg","Number of Background",0,0,4e9);

	NumBkg.setConstant(kTRUE);

	RooAddPdf MixPdf("MixPdf","MixPdf",RooArgList(SigPdf,BkgPdf),RooArgList(NumSig,NumBkg));
	// RooAddPdf MixPdf("MixPdf","MixPdf",RooArgList(pdf,BkgPdf),RooArgList(NumSig,NumBkg));

	RooDataHist dh_data_BkgSub("dh_data","dh_data",x,Import(*h_NDs_mkkbin_BkgSub_roofit));
	// RooDataHist dh_data("dh_data","dh_data",x,Import(*h_NDs_mkkbin_Fit_roofit));
	RooDataHist dh_data_roofit("dh_data","dh_data",x,Import(*h_NDs_mkkbin_Fit_roofit));
	// RooDataHist dh_data("dh_data","dh_data",x,Import(*h_NDs_mkkbin));
	RooDataHist *dh_data_ptr=(RooDataHist*)&dh_data_roofit;
	if(useBkgSubData){
		dh_data_ptr=(RooDataHist*)&dh_data_BkgSub;
	}
	RooDataHist dh_data=(RooDataHist)*dh_data_ptr;

	RooPlot *Data_frame=x.frame(Title(""));
	dh_data.plotOn(Data_frame);

	cbMean.setConstant(kTRUE);
	cbSigma.setConstant(kTRUE);
	cbSigma2.setConstant(kTRUE);
	cbSigma3.setConstant(kTRUE);
	cbN.setConstant(kTRUE);
	cbAlpha.setConstant(kTRUE);
	sigma1.setConstant(kTRUE);	
	sigma2.setConstant(kTRUE);	
	fr_cb.setConstant(kTRUE);
	fr_cb1.setConstant(kTRUE);
	fr_cb2.setConstant(kTRUE);
	fr_gaus1.setConstant(kTRUE);

	FloatScale.setConstant(kTRUE);

	MixPdf.fitTo(dh_data);
	MixPdf.fitTo(dh_data);
	// cbMean.setConstant(kFALSE);
	MixPdf.fitTo(dh_data);
	FloatScale.setConstant(kFALSE);
	MixPdf.fitTo(dh_data);
	NumBkg.setConstant(kFALSE);
	MixPdf.fitTo(dh_data);
	MixPdf.fitTo(dh_data);
	MixPdf.fitTo(dh_data);
	RooFitResult *fitresult=MixPdf.fitTo(dh_data,Save());


	MixPdf.plotOn(Data_frame,LineColor(2));
	MixPdf.plotOn(Data_frame,Components(BkgPdf),LineColor(4));
	Data_frame->SetTitle("");
	
	// Data_frame->Draw();

cout<< " , cbMean_MC_V   = "<<    cbMean_MC_V   <<endl; 
cout<< " , cbSigma_MC_V  = "<<    cbSigma_MC_V<<endl;
cout<< " , cbSigma2_MC_V = "<<    cbSigma2_MC_V<<endl;
cout<< " , cbSigma3_MC_V = "<<    cbSigma3_MC_V<<endl;
cout<< " , cbN_MC_V      = "<<    cbN_MC_V<<endl;
cout<< " , cbAlpha_MC_V  = "<<    cbAlpha_MC_V<<endl;
cout<< " , sigma1_MC_V   = "<<    sigma1_MC_V <<endl;
cout<< " , sigma2_MC_V   = "<<    sigma2_MC_V<<endl;
cout<< " , fr_cb_MC_V    = "<<    fr_cb_MC_V<<endl;
cout<< " , fr_gaus1_MC_V = "<<    fr_gaus1_MC_V<<endl;
cout<< " , fr_cb1_MC_V   = "<<    fr_cb1_MC_V<<endl;
cout<< " , fr_cb2_MC_V   = "<<    fr_cb2_MC_V  <<endl;





	Data_frame->Print("V");
  TString s_pdf="MixPdf_Norm[x]";
  RooGoF goftest(Data_frame->getHist("h_dh_data"),Data_frame->getCurve(s_pdf));
  cout<<"check 0"<<endl;
  goftest.setRange(x.getMin(),x.getMax());
  goftest.setRebin(5,true); // better rebin to make sure that all bins have >=5 expected events  
  // goftest.setRebin(10,false); // better rebin to make sure that all bins have >=5 expected events  
  int ndf=0;
  int d_ndf=fitresult->floatParsFinal().getSize();
  double pvalue_BC=0.0;
  double testStat_BC=0.0;
  goftest.BCChi2Test(pvalue_BC,testStat_BC,ndf,d_ndf);
  // file_pval << "BC: " << pvalue_BC << ", " << testStat_BC << endl;
  cout<<"BC: " << pvalue_BC << ", " << testStat_BC << endl;

  double pvalue_PC=0.0;
  double testStat_PC=0.0;
  goftest.PearsonChi2Test(pvalue_PC,testStat_PC,ndf,d_ndf);
  // file_pval << "BC: " << pvalue_BC << ", " << testStat_BC << endl;
  cout<<"PC: " << pvalue_PC << ", " << testStat_PC << endl;

  double pvalue_NC=0.0;
  double testStat_NC=0.0;
  goftest.NeymanChi2Test(pvalue_NC,testStat_NC,ndf,d_ndf);
  // file_pval << "BC: " << pvalue_BC << ", " << testStat_BC << endl;
  cout<<"NC: " << pvalue_NC << ", " << testStat_NC << endl;
 
 
  double pvalue_RooFit=0.0;
  double testStat_RooFit=0.0;
  goftest.RooFitChi2Test(pvalue_RooFit,testStat_RooFit,ndf,d_ndf);
  // file_pval << "RooFitChi2: " << pvalue_RooFit << ", " << testStat_RooFit << endl;
  cout<<"RooFitChi2: " << pvalue_RooFit << ", " << testStat_RooFit << endl;



	TCanvas *c_roofit_Data=new TCanvas("c_roofit_Data","c_roofit_Data");
	c_roofit_Data->cd();
	SetCanvas(c_roofit_Data);
	Data_frame->GetXaxis()->SetTitle("m_{KK} (GeV/c^{2})");
	Data_frame->GetYaxis()->SetTitle("Events / 2.5 MeV/c^{2}");
	Data_frame->GetXaxis()->SetTitleSize(0.045);
	Data_frame->GetXaxis()->SetTitleOffset(1.2);
	Data_frame->GetXaxis()->SetLabelSize(0.040);
	Data_frame->GetXaxis()->CenterTitle();
	Data_frame->GetYaxis()->CenterTitle();
	Data_frame->GetYaxis()->SetTitleSize(0.045);
	Data_frame->GetYaxis()->SetTitleOffset(1.25);
	Data_frame->GetYaxis()->SetLabelSize(0.040);
	Data_frame->Draw();

	texCmsPre->Draw();
	if(isPbPb==0){
		texColpp->Draw();
	}else if(isPbPb){
		texColPbPb->Draw();
	}

	RooArgSet nset(x);
  cout<<"SigPdf[Dmass] = "<< SigPdf.getVal(&nset)<<endl;

  RooAbsReal* igx = SigPdf.createIntegral(x) ;
  cout << "gx_Int[x] = " << igx->getVal() << endl ;

  // x.setRange("signal",1.01,1.028);
  x.setRange("signal",DtktkResmassCutMean-DtktkResmassCutWidth,DtktkResmassCutMean+DtktkResmassCutWidth);
/*
  x.setRange("signal_080",DtktkResmassCutMean-PhiMassScan_bins[0],DtktkResmassCutMean+PhiMassScan_bins[0]);
  x.setRange("signal_085",DtktkResmassCutMean-PhiMassScan_bins[1],DtktkResmassCutMean+PhiMassScan_bins[1]);
  x.setRange("signal_090",DtktkResmassCutMean-PhiMassScan_bins[2],DtktkResmassCutMean+PhiMassScan_bins[2]);
  x.setRange("signal_095",DtktkResmassCutMean-PhiMassScan_bins[3],DtktkResmassCutMean+PhiMassScan_bins[3]);
  x.setRange("signal_100",DtktkResmassCutMean-PhiMassScan_bins[4],DtktkResmassCutMean+PhiMassScan_bins[4]);
  x.setRange("signal_105",DtktkResmassCutMean-PhiMassScan_bins[5],DtktkResmassCutMean+PhiMassScan_bins[5]);
  x.setRange("signal_110",DtktkResmassCutMean-PhiMassScan_bins[6],DtktkResmassCutMean+PhiMassScan_bins[6]);
*/

//

  ofstream file_phiRatio;
  file_phiRatio.open(Form("./PhiRatio_%s.out",s_ppPbPb.Data()),std::fstream::in | std::fstream::out | std::fstream::app);

	for(int i=0; i<nbin_PhiMassScan; i++){
  x.setRange(Form("signal_%.0f",PhiMassScan_bins[i]*1e4),DtktkResmassCutMean-PhiMassScan_bins[i],DtktkResmassCutMean+PhiMassScan_bins[i]);
  RooAbsReal* igx_sig = SigPdf.createIntegral(x,NormSet(x),Range(Form("signal_%.0f",PhiMassScan_bins[i]*1e4) ) ) ;
  // cout << "Sig_Int[x|signal]_Norm[x] = " << igx_sig->getVal() << endl ; 

  RooAbsReal* igx_bkg = BkgPdf.createIntegral(x,NormSet(x),Range(Form("signal_%.0f",PhiMassScan_bins[i]*1e4) ) ) ;
  // cout << "Bkg_Int[x|signal]_Norm[x] = " << igx_bkg->getVal() << endl ; 

  cout<<" \n------------------ \n"<<endl;

  double N_Sig_2sig=NumSig.getValV()* igx_sig->getVal();
  double N_Bkg_2sig=NumBkg.getValV()* igx_bkg->getVal();
	double N_SigRelErr_2sig=NumSig.getValV()*igx_sig->getPropagatedError(*fitresult)/N_Sig_2sig;
	double N_BkgRelErr_2sig=NumBkg.getValV()*igx_bkg->getPropagatedError(*fitresult)/N_Bkg_2sig;

  // cout<<"N_Sig_2sig = "<<N_Sig_2sig<<endl;
  // cout<<"N_Bkg_2sig = "<<N_Bkg_2sig<<endl;

	double NPhiRatio=N_Sig_2sig/(N_Sig_2sig+N_Bkg_2sig);
	// double NPhiRatioErr=N_Sig_2sig*NumSig.getError()/NumSig.getValV()/(N_Sig_2sig+N_Bkg_2sig);
	double NPhiRatioErr=NPhiRatio*sqrt(N_SigRelErr_2sig*N_SigRelErr_2sig+N_BkgRelErr_2sig*N_BkgRelErr_2sig);

	cout<<"i = "<<i<<" , PhiMassWindow = "<<PhiMassScan_bins[i]<<" , NPhiRatio = "<<NPhiRatio<<endl;
		if(!useBkgSubData){
			file_phiRatio<<"i = "<<i<<" , PhiMassWindow = "<<PhiMassScan_bins[i]<<" , NPhiRatio = "<<NPhiRatio<<endl;
		}


	}
//


  RooAbsReal* igx_sig = SigPdf.createIntegral(x,NormSet(x),Range("signal")) ;
  cout << "Sig_Int[x|signal]_Norm[x] = " << igx_sig->getVal() << endl ; 

  RooAbsReal* igx_bkg = BkgPdf.createIntegral(x,NormSet(x),Range("signal")) ;
  cout << "Bkg_Int[x|signal]_Norm[x] = " << igx_bkg->getVal() << endl ; 

  cout<<" \n------------------ \n"<<endl;

  double N_Sig_2sig=NumSig.getValV()* igx_sig->getVal();
  double N_Bkg_2sig=NumBkg.getValV()* igx_bkg->getVal();
	double N_SigRelErr_2sig=NumSig.getValV()*igx_sig->getPropagatedError(*fitresult)/N_Sig_2sig;
	double N_BkgRelErr_2sig=NumBkg.getValV()*igx_bkg->getPropagatedError(*fitresult)/N_Bkg_2sig;

  cout<<"N_Sig_2sig = "<<N_Sig_2sig<<endl;
  cout<<"N_Bkg_2sig = "<<N_Bkg_2sig<<endl;

	double NPhiRatio=N_Sig_2sig/(N_Sig_2sig+N_Bkg_2sig);
	// double NPhiRatioErr=N_Sig_2sig*NumSig.getError()/NumSig.getValV()/(N_Sig_2sig+N_Bkg_2sig);
	double NPhiRatioErr=NPhiRatio*sqrt(N_SigRelErr_2sig*N_SigRelErr_2sig+N_BkgRelErr_2sig*N_BkgRelErr_2sig);

	// TLatex *tltx =new TLatex();
	shiftX=-0.01;
	shiftY=0.05;
	tltx->SetTextSize(0.042);
//  tl_binfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s_Data_binfit",str_PbPb.Data()));  shiftY-=oneshift;
  tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%.0f < D_{S}^{+} p_{T} < %.0f GeV/c ",DptLowGL,DptHighGL)); shiftY-=oneshift;
  // tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("#phi Yield = %.1f #pm %.1f ",N_Sig_2sig,N_Sig_2sig*NumSig.getError()/NumSig.getValV() ));  shiftY-=oneshift;
  // tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Bkg. Yield = %.1f #pm %.1f ",N_Bkg_2sig,N_Bkg_2sig*NumBkg.getError()/NumBkg.getValV() ));  shiftY-=oneshift;
  tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("#phi fraction = %.3f #pm %.3f ",NPhiRatio,NPhiRatioErr ));  shiftY-=oneshift;

	if(useBkgSubData){
  tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s",s_BkgSub.Data() ));  shiftY-=oneshift;
	}

	TH1D *hdata=new TH1D("hdata","",1,0,1);
	TH1D *hbkg=new TH1D("hbkg","",1,0,1);

	hdata->SetLineColor(2);
	hbkg->SetLineColor(4);

	double ylow=0.5;
	double yhi=0.7;
	if(useBkgSubData){
		ylow=0.45;
		yhi=0.65;	
	}
	
	TLegend *le_Data=new TLegend(0.18,ylow,0.45,yhi);
	le_Data->SetBorderSize(0);
	le_Data->SetTextSize(0.04);
	le_Data->AddEntry(hdata,"#phi + Background","l");
	le_Data->AddEntry(hbkg,"Background","l");
	le_Data->Draw("same");


 c_roofit_Data->SaveAs(Form("./plots/fitMkk/%s_Data_Mkk_roofit%s.png",s_type.Data(),s_BkgSub.Data()));
 SavePlotDirs(c_roofit_Data,Form("%s_Data_Mkk_roofit%s",s_type.Data(),s_BkgSub.Data() ),{"PhiRatio"},"gcf");

//	SavePlotDirs()


	// RooRealVar FloatScale("FloatScale","FloatScale",0,-0.6,0.6);

	cout<<"Data FloatScale = "<<FloatScale.getValV()<<" +- "<<FloatScale.getError()<<endl;

	return 1;


  // tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("BC pvalue = %.3f ",pvalue_BC ));  shiftY-=oneshift;
  // tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Pearson pvalue = %.3f ",pvalue_PC ));  shiftY-=oneshift;
  // tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Neyman pvalue = %.3f ",pvalue_NC ));  shiftY-=oneshift;
  tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("RooFit pvalue = %.3f ",pvalue_RooFit ));  shiftY-=oneshift;
 
 // tltx->DrawLatexNDC(textpo+shiftX.08,textposy+shiftY,Form("BkgSub Yield = %.1f #pm %.1f ",NPhibkgSub,NPhibkgSubError ));  shiftY-=oneshift;
  // tltx->DrawLatexNDC(textpo+shiftX.28,textposy+shiftY,Form("Float width = %.3f",f1_PhiMix->GetParameter(5))); shiftY-=oneshift;

  // if(s_ExtraName!="")
  // {
  // tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,s_ExtraName); shiftY-=oneshift;
  // }
  // SavePlotDirs(c_binfit_Data[count_c_binfit],Form("%s_Data_binfit_pt%.0fto%.0f_%s",str_PbPb.Data(),DptLow,DptHigh,extraName.Data() ),{"SignalFit",str_PbPb,"binfit"} );

c_roofit_Data->SaveAs(Form("./plots/fitMkk/%s_Data_Mkk_roofit%s_GOF.png",s_type.Data(),s_BkgSub.Data()));






	return 1;


	/////////////////////
	// end roofit test //
	/////////////////////

	TF1 *f1_PhiMC=fitPhiMC(h_mkk_MCP_Phi_fineBin);


// for BW fit

  TF1 *f1_bw=new TF1("f1_bw",mybwfun,phiFitLow,phiFitHigh,3);
  f1_bw->SetParameters(50,0.07,0.992);
  f1_bw->SetParLimits(1,0.000001,0.15);
  f1_bw->SetParLimits(2,0.975,1.02);
  f1_bw->SetParName(0,"const");
  f1_bw->SetParName(1,"sigma");
  f1_bw->SetParName(2,"mean");
  h_mkk_MCP_Phi_fineBin->Fit("f1_bw","QNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_bw","QNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_bw","LQNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_bw","L MQNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_bw","L M INO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_bw","L M IO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_bw","L M IO","",PhiMass-0.02,PhiMass+0.02);

  f1_bw->SetRange(phiFitLow,phiFitHigh);
  f1_bw->SetRange(PhiMass-0.02,PhiMass+0.02);

	TCanvas *ctest=new TCanvas();
	ctest->cd();
	h_mkk_MCP_Phi_fineBin->Draw();
	  f1_bw->Draw("same");
  ctest->SaveAs(Form("./plots/fitMkk/%s_MC_Mkk_BW.png",s_type.Data()));
/*
	TF1 *f1_cry=new TF1("f1_cry",crystalball_function,phiFitLow,phiFitHigh,5);
   // f1_cry yield, mean, sigma, alpha, n
  f1_cry->SetParameters(1000,1.019,0.002,-1,1);
  f1_cry->SetRange(phiFitLow,phiFitHigh);
	h_mkk_MCP_Phi_fineBin->Fit("f1_cry","QNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry","QNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry","LQNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry","L MQNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry","L M INO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry","L M IO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry","L M IO","",PhiMass-0.02,PhiMass+0.02);


  TCanvas *c_cry=new TCanvas("c_cry","c_cry");
   c_cry->cd();
	h_mkk_MCP_Phi_fineBin->Draw();
  f1_cry->Draw("same");
*/

	// TF1 *f1_cry2=new TF1("f1_cry2",crys_DoubleGaus_fun,phiFitLow,phiFitHigh,8);
  // f1_cry2->SetParameters(1000,1.019,0.002,-1.5,1,0.007,0.2,0);
	TF1 *f1_cry2=new TF1("f1_cry2",crys_DoubleGaus_2ndCheb_fun,phiFitLow,phiFitHigh,12);
  // f1_cry2->SetParameters(1000,1.0195,0.0025,-1.5,1.3,0.3,0.008,0, 0,0.0,0.0);
  f1_cry2->SetParameters(1000,1.0195,0.0027,-1.6,1.2,0.3,0.008,0, 0,0.0,0.0);

	f1_cry2->SetParLimits(1,1.017,1.022);
	f1_cry2->SetParLimits(2,0.00001,0.01);
	f1_cry2->SetParLimits(3,-100,-0.1);
	// f1_cry2->SetParLimits(4,-10,1000);
	f1_cry2->SetParLimits(5,0.0,0.05);
	f1_cry2->SetParLimits(6,0.0,100);

	f1_cry2->FixParameter(8,0);
	f1_cry2->FixParameter(9,0);
	f1_cry2->FixParameter(10,0);
	f1_cry2->FixParameter(11,0);

  // f1_cry2->SetRange(1.01,1.03);
	f1_cry2->FixParameter(7,0);
	// f1_cry2->SetParLimits(5,0,0.1);
/*
	TCanvas *c_tt2=new TCanvas("c_tt2","c_tt2");
	c_tt2->cd();
	f1_cry2->Draw();
  c_tt2->SaveAs(Form("./plots/fitMkk/%s_MC_Mkk_cry2Gaus_test.png",s_type.Data()));
*/


  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","QNO","",phiFitLow,phiFitHigh); 
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","QNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","LQNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L MQNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M INO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M IO","",phiFitLow,phiFitHigh);

/*	
	f1_cry2->ReleaseParameter(8);
	f1_cry2->ReleaseParameter(9);
	f1_cry2->ReleaseParameter(10);

	f1_cry2->SetParLimits(8,0,1000000);


  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M IO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M IO","",phiFitLow,phiFitHigh);
*/

	// phiFitLow=1.005;
	// phiFitHigh=1.04;

  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M IO","",PhiMass-0.02,PhiMass+0.02);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M IO","",1.01,1.03);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M IO","",1.01,1.03);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M IO","",1.01,1.03);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M IO","",1.005,1.035);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M IO","",1.005,1.035);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M IO","",1.005,1.035);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M IO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M IO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M IO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M IO","",phiFitLow,phiFitHigh);

	TFitResultPtr fitResutTest=h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M IOS","",phiFitLow,phiFitHigh);
	cout<<"Prob = "<<fitResutTest->Prob();

  TCanvas *c_cry2=new TCanvas("c_cry2","c_cry2");
  c_cry2->cd();
  f1_cry2->SetRange(phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Draw();
	TLatex *la_cry2=new TLatex();

	double shiftY=0;
	la_cry2->DrawLatexNDC(textposx-0.05,textposy+shiftY,Form("%s MC %.0f<D_{S} p_{T}<%.0f GeV/c",s_ppPbPb.Data(),DptLowGL,DptHighGL)); shiftY-=oneshift;
	la_cry2->DrawLatexNDC(textposx-0.05,textposy+shiftY,"CrystalBall +Gaus.Fit"); shiftY-=oneshift;
	la_cry2->DrawLatexNDC(textposx-0.05,textposy+shiftY,Form("Mean : %0.4f", f1_cry2->GetParameter(1))); shiftY-=oneshift;
	la_cry2->DrawLatexNDC(textposx-0.05,textposy+shiftY,Form("CB #sigma : %0.4f", f1_cry2->GetParameter(2))); shiftY-=oneshift;
	la_cry2->DrawLatexNDC(textposx-0.05,textposy+shiftY,Form("CB #alpha : %0.4f", f1_cry2->GetParameter(3))); shiftY-=oneshift;
	la_cry2->DrawLatexNDC(textposx-0.05,textposy+shiftY,Form("CB n : %0.4f", f1_cry2->GetParameter(4))); shiftY-=oneshift;
	la_cry2->DrawLatexNDC(textposx-0.05,textposy+shiftY,Form("Gaus #sigma : %0.4f", f1_cry2->GetParameter(5))); shiftY-=oneshift;


  // f1_cry2->SetRange(1.005,1.035);
  // f1_cry2->SetRange(PhiMass-0.02,PhiMass+0.02);
  f1_cry2->Draw("same");
  // f1_cry2->Draw();

  c_cry2->SaveAs(Form("./plots/fitMkk/%s_MC_Mkk_cry2Gaus.png",s_type.Data()));


	// return 1;

//////////////////////
////// cr2+bwf0 //////
//////////////////////
/*
	TF1 *f1_cry2_bwf0=new TF1("f1_cry2_bwf0",crys_DoubleGaus_bwf0_fun,phiFitLow,phiFitHigh,11);
  // f1_cry2->SetParameters(1000,1.0195,0.0025,-1.5,1.3,0.3,0.008,0, 0,0.0,0.0);
  f1_cry2_bwf0->SetParameters(1000,1.0195,0.0027,-1.6,1.2,0.3,0.008,0, 0,0.0,0.0);

	f1_cry2_bwf0->SetParLimits(1,1.017,1.022);
	f1_cry2_bwf0->SetParLimits(2,0.00001,0.01);
	f1_cry2_bwf0->SetParLimits(3,-100,-0.1);
	// f1_cry2_bwf0->SetParLimits(4,-10,1000);
	f1_cry2_bwf0->SetParLimits(5,0.0,0.05);
	f1_cry2_bwf0->SetParLimits(6,0.0,100);

// for bwf0
	f1_cry2_bwf0->FixParameter(8,0);
	f1_cry2_bwf0->FixParameter(9,0);
	f1_cry2_bwf0->FixParameter(10,0);

  // f1_cry2_bwf0->SetRange(1.01,1.03);
	f1_cry2_bwf0->FixParameter(7,0);
	// f1_cry2_bwf0->SetParLimits(5,0,0.1);


  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2_bwf0","QNO","",phiFitLow,phiFitHigh); 
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2_bwf0","QNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2_bwf0","LQNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2_bwf0","L MQNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2_bwf0","L M INO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2_bwf0","L M IO","",phiFitLow,phiFitHigh);


	// phiFitLow=1.005;
	// phiFitHigh=1.04;

  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2_bwf0","L M IO","",PhiMass-0.02,PhiMass+0.02);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2_bwf0","L M IO","",1.01,1.03);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2_bwf0","L M IO","",1.01,1.03);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2_bwf0","L M IO","",1.01,1.03);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2_bwf0","L M IO","",1.005,1.035);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2_bwf0","L M IO","",1.005,1.035);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2_bwf0","L M IO","",1.005,1.035);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2_bwf0","L M IO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2_bwf0","L M IO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2_bwf0","L M IO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2_bwf0","L M IO","",phiFitLow,phiFitHigh);
   // |)}>#
  TCanvas *c_cry2_bwf0=new TCanvas("c_cry2_bwf0","c_cry2_bwf0");
  c_cry2_bwf0->cd();
  f1_cry2_bwf0->SetRange(phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Draw();
	TLatex *la_cry2_bwf0=new TLatex();

	shiftY=0;
	la_cry2_bwf0->DrawLatexNDC(textposx-0.05,textposy+shiftY,Form("%s MC %.0f<D_{S} p_{T}<%.0f GeV/c",s_ppPbPb.Data(),DptLowGL,DptHighGL)); shiftY-=oneshift;
	la_cry2_bwf0->DrawLatexNDC(textposx-0.05,textposy+shiftY,"CrystalBall +Gaus.Fit"); shiftY-=oneshift;
	la_cry2_bwf0->DrawLatexNDC(textposx-0.05,textposy+shiftY,Form("Mean : %0.4f", f1_cry2_bwf0->GetParameter(1))); shiftY-=oneshift;
	la_cry2_bwf0->DrawLatexNDC(textposx-0.05,textposy+shiftY,Form("CB #sigma : %0.4f", f1_cry2_bwf0->GetParameter(2))); shiftY-=oneshift;
	la_cry2_bwf0->DrawLatexNDC(textposx-0.05,textposy+shiftY,Form("CB #alpha : %0.4f", f1_cry2_bwf0->GetParameter(3))); shiftY-=oneshift;
	la_cry2_bwf0->DrawLatexNDC(textposx-0.05,textposy+shiftY,Form("CB n : %0.4f", f1_cry2_bwf0->GetParameter(4))); shiftY-=oneshift;
	la_cry2_bwf0->DrawLatexNDC(textposx-0.05,textposy+shiftY,Form("Gaus #sigma : %0.4f", f1_cry2_bwf0->GetParameter(5))); shiftY-=oneshift;


  // f1_cry2_bwf0->SetRange(1.005,1.035);
  // f1_cry2_bwf0->SetRange(PhiMass-0.02,PhiMass+0.02);
  f1_cry2_bwf0->Draw("same");
   // f1_cry2_bwf0->Draw();

  c_cry2_bwf0->SaveAs(Form("./plots/fitMkk/%s_MC_Mkk_cry2_bwf0Gaus.png",s_type.Data()));
*/

	// return 1;

	// TF1 *f1_cry2=new TF1("f1_cry2",crys_DoubleGaus_fun,phiFitLow,phiFitHigh,8);
  // f1_cry2->SetParameters(1000,1.019,0.002,-1.5,1,0.007,0.2,0);
/*	TF1 *f1_cry2=new TF1("f1_cry2",crys_DoubleGaus_2ndCheb_fun,phiFitLow,phiFitHigh,12);
  f1_cry2->SetParameters(1000,1.019,0.002,-1,1,0.3,0.004,0, 0,0.0,0.0);
	f1_cry2->FixParameter(8,0);
	f1_cry2->FixParameter(9,0);
	f1_cry2->FixParameter(10,0);
*/
  // f1_cry2->SetRange(1.01,1.03);
	 // f1_cry2->FixParameter(7,0);
	// f1_cry2->SetParLimits(5,0,0.1);
/*
	TCanvas *c_tt2=new TCanvas("c_tt2","c_tt2");
	c_tt2->cd();
	f1_cry2->Draw();
  c_tt2->SaveAs(Form("./plots/fitMkk/%s_MC_Mkk_cry2Gaus_test.png",s_type.Data()));
*/



	int testFitCry=1;

	if(testFitCry==1){

	for(int i=1; i<7; i++){
		f1_cry2->FixParameter(i,f1_cry2->GetParameter(i));
	//	f1_cry2_bwf0->FixParameter(i,f1_cry2->GetParameter(i));
	}


  h_NDs_mkkbin_Fit->Fit("f1_cry2","QNO","",phiFitLow,phiFitHigh); 
  h_NDs_mkkbin_Fit->Fit("f1_cry2","QNO","",phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Fit("f1_cry2","LQNO","",phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Fit("f1_cry2","L MQNO","",phiFitLow,phiFitHigh);
   h_NDs_mkkbin_Fit->Fit("f1_cry2","L M INO","",phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Fit("f1_cry2","L M IO","",phiFitLow,phiFitHigh);
 

  h_NDs_mkkbin_Fit->Fit("f1_cry2","L M IO","",PhiMass-0.02,PhiMass+0.02);
  h_NDs_mkkbin_Fit->Fit("f1_cry2","L M IO","",1.01,1.03);
  h_NDs_mkkbin_Fit->Fit("f1_cry2","L M IO","",1.01,1.03);
  h_NDs_mkkbin_Fit->Fit("f1_cry2","L M IO","",1.01,1.03);
  h_NDs_mkkbin_Fit->Fit("f1_cry2","L M IO","",1.005,1.035);
  h_NDs_mkkbin_Fit->Fit("f1_cry2","L M IO","",1.005,1.035);
  h_NDs_mkkbin_Fit->Fit("f1_cry2","L M IO","",1.005,1.035);
  h_NDs_mkkbin_Fit->Fit("f1_cry2","L M IO","",phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Fit("f1_cry2","L M IO","",phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Fit("f1_cry2","L M IO","",phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Fit("f1_cry2","L M IO","",phiFitLow,phiFitHigh);

	 f1_cry2->ReleaseParameter(7);
	f1_cry2->SetParLimits(7,-0.6,0.6);
	f1_cry2->ReleaseParameter(8);
	f1_cry2->SetParLimits(8,0,1000000);
	f1_cry2->ReleaseParameter(9);
	 f1_cry2->ReleaseParameter(10);

  h_NDs_mkkbin_Fit->Fit("f1_cry2","L M IO","",phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Fit("f1_cry2","L M IO","",phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Fit("f1_cry2","L M IO","",phiFitLow,phiFitHigh);

	TFitResultPtr fitResutTest2=h_NDs_mkkbin_Fit->Fit("f1_cry2","L M IOS","",phiFitLow,phiFitHigh);
	cout<<"Prob = "<<fitResutTest2->Prob();


  TCanvas *c_cry2data=new TCanvas("c_cry2data","c_cry2data");
  c_cry2data->cd();
  f1_cry2->SetRange(phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Draw();
  // f1_cry2->SetRange(1.005,1.035);
  // f1_cry2->SetRange(PhiMass-0.02,PhiMass+0.02);
  f1_cry2->Draw("same");
  // f1_cry2->Draw();

  c_cry2data->SaveAs(Form("./plots/fitMkk/%s_Data_Mkk_cry2Gaus.png",s_type.Data()));

	return 1;

//////////////
 // for bwf0 //
//////////////
/*
  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","QNO","",phiFitLow,phiFitHigh); 
  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","QNO","",phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","LQNO","",phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","L MQNO","",phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","L M INO","",phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","L M IO","",phiFitLow,phiFitHigh);
 

  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","L M IO","",PhiMass-0.02,PhiMass+0.02);
  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","L M IO","",1.01,1.03);
  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","L M IO","",1.01,1.03);
  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","L M IO","",1.01,1.03);
  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","L M IO","",1.005,1.035);
   h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","L M IO","",1.005,1.035);
  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","L M IO","",1.005,1.035);
  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","L M IO","",phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","L M IO","",phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","L M IO","",phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","L M IO","",phiFitLow,phiFitHigh);

	f1_cry2_bwf0->ReleaseParameter(7);
	f1_cry2_bwf0->SetParLimits(7,-0.6,0.6);
	f1_cry2_bwf0->ReleaseParameter(8);
	f1_cry2_bwf0->SetParLimits(8,0,1000000);
	f1_cry2_bwf0->ReleaseParameter(9);
	f1_cry2_bwf0->SetParLimits(9,0.0001,0.1);
	f1_cry2_bwf0->ReleaseParameter(10);
	f1_cry2_bwf0->SetParLimits(10,0.97,1.01);

  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","L M IO","",phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","L M IO","",phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Fit("f1_cry2_bwf0","L M IO","",phiFitLow,phiFitHigh);

  TCanvas *c_cry2_bwf0data=new TCanvas("c_cry2_bwf0data","c_cry2_bwf0data");
  c_cry2_bwf0data->cd();
  f1_cry2_bwf0->SetRange(phiFitLow,phiFitHigh);
  h_NDs_mkkbin_Fit->Draw();
  // f1_cry2_bwf0->SetRange(1.005,1.035);
   // f1_cry2_bwf0->SetRange(PhiMass-0.02,PhiMass+0.02);
  f1_cry2_bwf0->Draw("same");
  // f1_cry2_bwf0->Draw();

  c_cry2_bwf0data->SaveAs(Form("./plots/fitMkk/%s_Data_Mkk_cry2_bwf0Gaus.png",s_type.Data()));
*/



	}



/*
	 TF1 *f1_cryBW=new TF1("f1_cryBW",crys_BW_fun,phiFitLow,phiFitHigh,7);
  f1_cryBW->SetParameters(1000,1.019,0.002,-1,1,0.3,0.3,0);
  f1_cryBW->SetRange(phiFitLow,phiFitHigh);
	f1_cryBW->FixParameter(7,0);

  h_mkk_MCP_Phi_fineBin->Fit("f1_cryBW","QNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cryBW","QNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cryBW","LQNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cryBW","L MQNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cryBW","L M INO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cryBW","L M IO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cryBW","L M IO","",PhiMass-0.02,PhiMass+0.02);
  // h_mkk_MCP_Phi_fineBin->Fit("f1_cryBW","L M IO","",1.01,1.03);

 TCanvas *c_cryBW=new TCanvas("c_cryBW","c_cryBW");
  c_cryBW->cd();
  h_mkk_MCP_Phi_fineBin->Draw();
  f1_cryBW->Draw("same");
*/


/*
	TF1 *f1_2cry=new TF1("f1_2cry",doublecrystalball_fun,phiFitLow,phiFitHigh,8);
  f1_2cry->SetParameters(1000,1.019,0.002,-1,1,0.3,-1,1,0.5);
  f1_2cry->SetRange(phiFitLow,phiFitHigh);
	// f1_2cry->FixParameter(7,0);
  h_mkk_MCP_Phi_fineBin->Fit("f1_2cry","QNO","",phiFitLow,phiFitHigh);
   h_mkk_MCP_Phi_fineBin->Fit("f1_2cry","QNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_2cry","LQNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_2cry","L MQNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_2cry","L M INO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_2cry","L M IO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_2cry","L M IO","",PhiMass-0.02,PhiMass+0.02);
  // h_mkk_MCP_Phi_fineBin->Fit("f1_2cry","L M IO","",1.01,1.03);

 TCanvas *c_2cry=new TCanvas("c_2cry","c_2cry");
  c_2cry->cd();
  h_mkk_MCP_Phi_fineBin->Draw();
  f1_2cry->Draw("same");
*/

 


	// return (p[0]* (crystalball_function(x[0], p[3], p[4], p[2], p[1]) +  p[8]*crystalball_function(x[0], p[6], p[7], p[5], p[1])));


// end BW fit

	// h_NDs_mkkbin_Fit->Draw();
	// h_NDs_mkkbin_BkgSub->Draw("same");

	// h_mkk_MCP_Phi_fineBin->Draw();
	
	// return 1;
// #<{(|
	fout->cd();

  double NPhiInCut;
  double NPhiInCutErr;
  double NBkgInCut;
  double NBkgInCutErr;

//  TF1 *f1Data=fitPhiMass(h_NDs_mkkbin_Fit,f1_PhiMC, NPhiInCut, NPhiInCutErr,NBkgInCut ,NBkgInCutErr);
  // TF1 *f1Data=fitPhiMass(h_NDs_mkkbin_Fit,f1_cry2, NPhiInCut, NPhiInCutErr,NBkgInCut ,NBkgInCutErr);

	TH1D *h_DsPhiRatio_Fit=new TH1D("h_DsPhiRatio_Fit","h_DsPhiRatio_Fit",1,0,1); h_DsPhiRatio_Fit->Sumw2();
	 h_DsPhiRatio_Fit->SetBinContent(1,NPhiInCut/(NPhiInCut+NBkgInCut));
	 h_DsPhiRatio_Fit->SetBinError(1,NPhiInCutErr/(NPhiInCut+NBkgInCut));
	h_DsPhiRatio_Fit->Write();

  // TF1 *f1Data2=fitPhiMass(h_NDs_mkkbin_BkgSub,f1_PhiMC, NPhiInCut, NPhiInCutErr,NBkgInCut ,NBkgInCutErr,"BkgSub");
  // TF1 *f1Data2=fitPhiMass(h_NDs_mkkbin_BkgSub,f1_cry2, NPhiInCut, NPhiInCutErr,NBkgInCut ,NBkgInCutErr,"BkgSub");
	TH1D *h_DsPhiRatio_BkgSub=new TH1D("h_DsPhiRatio_BkgSub","h_DsPhiRatio_BkgSub",1,0,1); h_DsPhiRatio_BkgSub->Sumw2();
	h_DsPhiRatio_BkgSub->SetBinContent(1,NPhiInCut/(NPhiInCut+NBkgInCut));
	h_DsPhiRatio_BkgSub->SetBinError(1,NPhiInCutErr/(NPhiInCut+NBkgInCut));
	h_DsPhiRatio_BkgSub->Write();
// |)}>#

  // TF1 *f1Data3=fitPhiMass(h_NDs_mkkbin_Fit,f1_PhiMC, NPhiInCut, NPhiInCutErr,NBkgInCut ,NBkgInCutErr,"BW");
  // TF1 *f1Data4=fitPhiMass(h_NDs_mkkbin_BkgSub,f1_PhiMC, NPhiInCut, NPhiInCutErr,NBkgInCut ,NBkgInCutErr,"BkgSub_BW");


	TF1 *f1Data_cry2Gaus=fitPhiMass_cry2Gaus(h_NDs_mkkbin_Fit,f1_cry2, NPhiInCut, NPhiInCutErr,NBkgInCut ,NBkgInCutErr);
	TF1 *f1Data_cry2Gaus_BkgSub=fitPhiMass_cry2Gaus(h_NDs_mkkbin_BkgSub,f1_cry2, NPhiInCut, NPhiInCutErr,NBkgInCut ,NBkgInCutErr,"BkgSub");




	return 0;
}




TF1 *fitPhiMC(TH1D* h_PhiMC){

  // TF1 *f1_PhiMC = new TF1("f1_PhiMC","[0]*( [4]*TMath::Gaus(x,[1],[2]*(1+[5]))/( sqrt(2*TMath::Pi())*[2]*(1+[5]) ) + (1-[4])*TMath::Gaus(x,[1],[3]*(1+[5]))/(sqrt(2*TMath::Pi())*[3]*(1+[5])  ) )",PhiMass-PhiMassCandWidth,PhiMass+PhiMassCandWidth); // kTRUE nomalize Gaus
  TF1 *f1_PhiMC = new TF1("f1_PhiMC","[0]*( [4]*TMath::Gaus(x,[1],[2]*(1+[5]))/( sqrt(2*TMath::Pi())*[2]*(1+[5]) ) + (1-[4])*TMath::Gaus(x,[1],[3]*(1+[5]))/(sqrt(2*TMath::Pi())*[3]*(1+[5])  ) )"); // kTRUE nomalize Gaus
  f1_PhiMC->SetParameter(0, h_PhiMC->GetBinWidth(1) * h_PhiMC->Integral( h_PhiMC->FindBin(PhiMass-PhiMassCandWidth) , h_PhiMC->FindBin(PhiMass-PhiMassCandWidth) ) ) ;
   f1_PhiMC->SetParameter(1,PhiMass);   // DsMean
  f1_PhiMC->SetParLimits(1,PhiMass-PhiMassCandFitMeanWidth,PhiMass+PhiMassCandFitMeanWidth);
  f1_PhiMC->SetParameter(2,0.01); // defaut (wider) gaus1 sigma
  f1_PhiMC->SetParLimits(2,0.0001,0.1);
  f1_PhiMC->SetParameter(3,0.001); // default (narrow) gaus2 sigma
  f1_PhiMC->SetParLimits(3,0.00001,0.015);
  f1_PhiMC->SetParameter(4,0.5);  // fraction of Gaus1 / All
  f1_PhiMC->SetParLimits(4,0,1);
  f1_PhiMC->SetParameter(5,0); // for data /mc discrepancy
  f1_PhiMC->SetParLimits(5,0,0);
  f1_PhiMC->FixParameter(5,0);

  f1_PhiMC->SetLineColor(kRed);

  h_PhiMC->Fit("f1_PhiMC","QN0","",    PhiMass-PhiMassCandWidth,PhiMass+PhiMassCandWidth);
  h_PhiMC->Fit("f1_PhiMC","QN0","",    f1_PhiMC->GetParameter(1)-PhiMassCandWidth,f1_PhiMC->GetParameter(1)+PhiMassCandWidth);
  h_PhiMC->Fit("f1_PhiMC","M IN0","",    f1_PhiMC->GetParameter(1)-PhiMassCandWidth,f1_PhiMC->GetParameter(1)+PhiMassCandWidth);
  h_PhiMC->Fit("f1_PhiMC","L QN0","",  f1_PhiMC->GetParameter(1)-PhiMassCandWidth,f1_PhiMC->GetParameter(1)+PhiMassCandWidth);
  h_PhiMC->Fit("f1_PhiMC","LN0","",    f1_PhiMC->GetParameter(1)-PhiMassCandWidth,f1_PhiMC->GetParameter(1)+PhiMassCandWidth);
  h_PhiMC->Fit("f1_PhiMC","L MN0","",  f1_PhiMC->GetParameter(1)-PhiMassCandWidth,f1_PhiMC->GetParameter(1)+PhiMassCandWidth);
  h_PhiMC->Fit("f1_PhiMC","L MN0","",  f1_PhiMC->GetParameter(1)-PhiMassCandWidth,f1_PhiMC->GetParameter(1)+PhiMassCandWidth);
  h_PhiMC->Fit("f1_PhiMC","L M IN0","",f1_PhiMC->GetParameter(1)-PhiMassCandWidth,f1_PhiMC->GetParameter(1)+PhiMassCandWidth);
  h_PhiMC->Fit("f1_PhiMC","L M IN0","",f1_PhiMC->GetParameter(1)-PhiMassCandWidth,f1_PhiMC->GetParameter(1)+PhiMassCandWidth);
  h_PhiMC->Fit("f1_PhiMC","L M IN0","",f1_PhiMC->GetParameter(1)-PhiMassCandWidth,f1_PhiMC->GetParameter(1)+PhiMassCandWidth);
  h_PhiMC->Fit("f1_PhiMC","L ME I0 S","",f1_PhiMC->GetParameter(1)-PhiMassCandWidth,f1_PhiMC->GetParameter(1)+PhiMassCandWidth);

  h_PhiMC->Fit("f1_PhiMC","L ME I0 S","",PhiMass-0.02,PhiMass+0.02);

 double yieldDsMC= f1_PhiMC->Integral(f1_PhiMC->GetParameter(1)-PhiMassCandWidth,f1_PhiMC->GetParameter(1)+PhiMassCandWidth)/ h_PhiMC->GetBinWidth(1);
  double yieldDsMCErr = f1_PhiMC->Integral(PhiMass-PhiMassCandWidth,PhiMass+PhiMassCandWidth)/h_PhiMC->GetBinWidth(1) * f1_PhiMC->GetParError(0)/f1_PhiMC->GetParameter(0);

  TCanvas *c_MCfit=new TCanvas("c_MCfit","c_MCfit");
  c_MCfit->cd();
  h_PhiMC->Draw();

  shiftY=0;
  TLatex *tl_binfitMC=new TLatex();
  tl_binfitMC->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s D_{S} %.0f < pt < %.0f ",s_ppPbPb.Data(),DptLowGL,DptHighGL)); shiftY-=oneshift;
 tl_binfitMC->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Mean=%.4f, frGaus1=%.3f  ",f1_PhiMC->GetParameter(1),f1_PhiMC->GetParameter(4) ));  shiftY-=oneshift;
  tl_binfitMC->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Sigma1=%.4f, Sigma2=%.4f",f1_PhiMC->GetParameter(2), f1_PhiMC->GetParameter(3)));  shiftY-=oneshift;


	f1_PhiMC->SetRange(f1_PhiMC->GetParameter(1)-PhiMassCandWidth,f1_PhiMC->GetParameter(1)+PhiMassCandWidth);
	f1_PhiMC->SetRange(PhiMass-0.02,PhiMass+0.02);
  f1_PhiMC->Draw("same");


	cout<<"s_type = "<<s_type<<endl;

  c_MCfit->SaveAs(Form("./plots/fitMkk/%s_MC_Mkk.png",s_type.Data()));



	return f1_PhiMC;

}

TF1* fitPhiMass(TH1D *h_PhiMassData, TF1 *f1MC,double &NPhiInCut,double &NPhiInCutErr,double &NBkgInCut,double &NBkgInCutErr, TString s_ExtraName, int fixSigShape, TF1 *f1_default){

  // TF1* f = new TF1("fMass","[0]*([7]*([9]*TMath::Gaus(x,[1],[2]*(1+[11]))/(sqrt(2*3.1415927)*[2]*(1+[11]))+(1-[9])*TMath::Gaus(x,[1],[10]*(1+[11]))/(sqrt(2*3.1415927)*[10]*(1+[11])))+(1-[7]\
)*TMath::Gaus(x,[1],[8])/(sqrt(2*3.1415927)*[8]))+[3]+[4]*x+[5]*x*x+[6]*x*x*x", 1.7, 2.0);

	// FitYield=472;

	 cout<<"check -1"<<endl;
  cout<<"f1MC->GetParameter(1) = "<<f1MC->GetParameter(1)<<endl;
	

  TF1 *f1_PhiMC = new TF1("f1_PhiMC","[0]*( [4]*TMath::Gaus(x,[1],[2]*(1+[5]))/( sqrt(2*TMath::Pi())*[2]*(1+[5]) ) + (1-[4])*TMath::Gaus(x,[1],[3]*(1+[5]))/(sqrt(2*TMath::Pi())*[3]*(1+[5])  ) )"); // kTRUE nomalize Gaus

  // f1_PhiMC->SetParameter(0, h_PhiMassPhiMC->GetBinWidth(1) * h_PhiMassPhiMC->Integral( h_PhiMassPhiMC->FindBin(PhiMass-PhiMassCandWidth) , h_PhiMassPhiMC->FindBin(PhiMass-PhiMassCandWidth) ) ) ;
  f1_PhiMC->SetParameter(1,f1MC->GetParameter(1));   // PhiMean
  f1_PhiMC->SetParLimits(1,PhiMass-PhiMassCandFitMeanWidth,PhiMass+PhiMassCandFitMeanWidth);
  f1_PhiMC->SetParameter(2,f1MC->GetParameter(2)); // defaut (wider) gaus1 sigma
 f1_PhiMC->SetParLimits(2,0.0001,0.1);
  f1_PhiMC->SetParameter(3,f1MC->GetParameter(3)); // default (narrow) gaus2 sigma
  f1_PhiMC->SetParLimits(3,0.00001,0.1);
  f1_PhiMC->SetParameter(4,f1MC->GetParameter(4));  // fraction of Gaus1 / All
  f1_PhiMC->SetParLimits(4,0,1);
  f1_PhiMC->SetParameter(5,0); // for data /mc discrepancy
//  f1_PhiMC->SetParLimits(5,0,0);
//  f1_PhiMC->FixParameter(5,0);

  f1_PhiMC->SetLineColor(kRed);


	cout<<"check 0"<<endl;


	cout<<"f1MC->GetNpar() = "<<f1MC->GetNpar()<<endl;

	int Npar_Sig=f1MC->GetNpar();

  TF1 *f1_PhiSignal = (TF1*)f1_PhiMC->Clone("f1_PhiSignal");
//  TF1 *f1_PhiSignal = (TF1*)f1MC->Clone("f1_PhiSignal");
  TF1 *f1_PhiBkg= new TF1("f1_PhiBkg","[0]*(1+[1]*x+[2]*(2*x*x-1))");
  TF1 *f1_PhiBkg_clone=(TF1*)f1_PhiBkg->Clone("f1_PhiBkg_clone");  // stupid bug in root, must use clone function
  TF1 *f1_PhiMix = new TF1("f1_PhiMix","f1_PhiSignal+f1_PhiBkg");
/*
	for(int i=0 ;i <f1MC->GetNpar(); i++){
		f1_PhiMix->SetParameter(i,f1MC->GetParameter(i));
		if(i>0 && i< f1MC->GetNpar()-1 ){
		 f1_PhiMix->FixParameter(i,f1MC->GetParameter(i));
		 }
	}
*/

  f1_PhiMix->SetParameter(5,0); // width ratio difference for Data/ MC, data width = (1+5)*MC width
   f1_PhiMix->SetParLimits(5,-1,1);
	f1_PhiMix->SetParameter(0,0.9*h_PhiMassData->Integral()*h_PhiMassData->GetBinWidth(1));
  f1_PhiMix->FixParameter(1,f1_PhiMC->GetParameter(1));  // Phi Mass mean
  f1_PhiMix->FixParameter(2,f1_PhiMC->GetParameter(2));  // Phi Sigma1
  f1_PhiMix->FixParameter(3,f1_PhiMC->GetParameter(3));  // Phi Sigma2
  f1_PhiMix->FixParameter(4,f1_PhiMC->GetParameter(4));  // Phi Gaus1/all yield ratio
  f1_PhiMix->SetLineColor(kRed);
//  f1_PhiMix->SetParameter(6,1000);
  f1_PhiMix->SetParameter(7,0);
  f1_PhiMix->SetParameter(6,0.1*h_PhiMassData->Integral()*h_PhiMassData->GetBinWidth(1));
	f1_PhiMix->SetParLimits(6,0,h_PhiMassData->Integral()*h_PhiMassData->GetBinWidth(1)*10);



	if(fixSigShape==1 && f1_default!=NULL){
		cout<<"fixSigShape fit"<<endl;
		// cout<<"f1_default->GetParameter(7) = "<<f1_default->GetParameter(7)<<endl;
		// cout<<"f1_default->GetParameter(8) = "<<f1_default->GetParameter(8)<<endl;

		f1_PhiMix->SetParameters(0.1*h_PhiMassData->Integral()*h_PhiMassData->GetBinWidth(1), f1_default->GetParameter(1), f1_default->GetParameter(2), f1_default->GetParameter(3), f1_default->GetParameter(4), f1_default->GetParameter(5), f1_default->GetParameter(6), f1_default->GetParameter(7), f1_default->GetParameter(8) );

		f1_PhiMix->FixParameter(1,f1_default->GetParameter(1));  // Phi Mass mean
		f1_PhiMix->FixParameter(2,f1_default->GetParameter(2));  
		f1_PhiMix->FixParameter(3,f1_default->GetParameter(3));  
		f1_PhiMix->FixParameter(4,f1_default->GetParameter(4));  
		f1_PhiMix->FixParameter(5,f1_default->GetParameter(5));  

		f1_PhiMix->FixParameter(7,f1_default->GetParameter(7));  
		f1_PhiMix->FixParameter(8,f1_default->GetParameter(8));  

	}
 
	cout<<"check 1"<<endl;
 
   h_PhiMassData->Fit("f1_PhiMix","QN0","",    PhiMassFitLow, PhiMassFitHigh );
  h_PhiMassData->Fit("f1_PhiMix","QN0","",    PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L QN0","",  PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","LN0","",    PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L MN0","",  PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L MN0","",  PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L M IN0","",PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L M IN0","",PhiMassFitLow, PhiMassFitHigh);

	if(fixSigShape!=1){
  f1_PhiMix->ReleaseParameter(5);
  f1_PhiMix->ReleaseParameter(1);
	}

  h_PhiMassData->Fit("f1_PhiMix","QN0","",    PhiMassFitLow, PhiMassFitHigh );
  h_PhiMassData->Fit("f1_PhiMix","QN0","",    PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L QN0","",  PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","LN0","",    PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L MN0","",  PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L MN0","",  PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L M IN0","",PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L M I0","",PhiMassFitLow, PhiMassFitHigh);


	cout<<"check 2"<<endl;

	f1_PhiMix->SetRange(PhiMassFitLow,PhiMassFitHigh);
	int fitStatus=1;
  int fitIsValid=0;
  TFitResultPtr fitResut;
   double fitPrecision=1.e-8;
  while(fitStatus>0 && fitStatus!=4000){
    TFitter::SetPrecision(fitPrecision);
    fitResut=h_PhiMassData->Fit("f1_PhiMix","L EMI S0","",PhiMassFitLow, PhiMassFitHigh);
    fitStatus=fitResut->Status();
    fitIsValid=fitResut->IsValid();
    cout<<"fit precision: "<<TFitter::GetPrecision()<<"   status: "<<fitStatus<<" isValid = "<< fitIsValid <<endl;
    if(fitStatus){
      fitPrecision *= 3;
    }
     if(fitPrecision>1.e-3) {break;}
   }
	fitResut=h_PhiMassData->Fit("f1_PhiMix","L EMI S0","",PhiMassFitLow, PhiMassFitHigh);
 cout<<"f1_PhiMix->GetParameter(0) = "<<f1_PhiMix->GetParameter(0)<<endl;
cout<<"f1_PhiSignal->GetParameter(0) = "<<f1_PhiSignal->GetParameter(0)<<endl;
/* 
  for(int i=0 ;i <f1MC->GetNpar(); i++){
		f1MC->SetParameter(i, f1_PhiMix->GetParameter(i));
		f1MC->SetParError(i, f1_PhiMix->GetParError(i));
	}
*/


	f1_PhiSignal->SetParameters(f1_PhiMix->GetParameter(0), f1_PhiMix->GetParameter(1) , f1_PhiMix->GetParameter(2) , f1_PhiMix->GetParameter(3) , f1_PhiMix->GetParameter(4), f1_PhiMix->GetParameter(5));
  f1_PhiSignal->SetParError(0,f1_PhiMix->GetParError(0));
  f1_PhiSignal->SetParError(1,f1_PhiMix->GetParError(1));
  f1_PhiSignal->SetParError(2,f1_PhiMix->GetParError(2));
  f1_PhiSignal->SetParError(3,f1_PhiMix->GetParError(3));
  f1_PhiSignal->SetParError(4,f1_PhiMix->GetParError(4));
  f1_PhiSignal->SetParError(5,f1_PhiMix->GetParError(5));


	cout<<"check 3"<<endl;

  f1_PhiSignal->SetFillColor(kOrange-3);
  f1_PhiSignal->SetFillStyle(3002);
  f1_PhiSignal->SetLineColor(kOrange-3);
  f1_PhiSignal->SetLineWidth(2);
  f1_PhiSignal->SetLineStyle(2);
  f1_PhiSignal->SetRange(PhiMassFitLow,PhiMassFitHigh); // must setrange before plot or get empty
  // double yieldPhiSignal= f1_PhiSignal->Integral(PhiMassFitLow,PhiMassFitHigh)/ h_PhiMassData->GetBinWidth(1);
  double yieldPhiSignal= f1_PhiSignal->Integral(PhiMass-PhiMassCandWidth,PhiMass+PhiMassCandWidth)/ h_PhiMassData->GetBinWidth(1);
  double yieldPhiSignalErr = f1_PhiSignal->Integral(PhiMass-PhiMassCandWidth,PhiMass+PhiMassCandWidth)/h_PhiMassData->GetBinWidth(1) * f1_PhiSignal->GetParError(0)/f1_PhiSignal->GetParameter(0);

  double yieldPhiSignalTry = f1_PhiSignal->GetParameter(0)/h_PhiMassData->GetBinWidth(1);
  double yieldPhiSignalTryErr = f1_PhiSignal->GetParError(0)/h_PhiMassData->GetBinWidth(1);


	 cout<<"check 4"<<endl;
 
  cout<<"parameter 0 = "<<f1_PhiSignal->GetParameter(0)<<" BinWidth = "<< h_PhiMassData->GetBinWidth(1) <<endl;
  cout<<"yield PhiSignal = "<< yieldPhiSignal <<" +- "<<yieldPhiSignalErr<<endl;
  cout<<"yield PhiSignal from Para 0 = "<< yieldPhiSignalTry <<" +- "<<yieldPhiSignalTryErr<<endl;
/*
	for(int i=0; i<3; i++){	
	f1_PhiBkg->SetParameter(i,f1_PhiMix->GetParameter(i+Npar_Sig));
	f1_PhiBkg->SetParError(i,f1_PhiMix->GetParError(i+Npar_Sig));

	}
*/

  f1_PhiBkg->SetParameters(f1_PhiMix->GetParameter(6), f1_PhiMix->GetParameter(7) ,f1_PhiMix->GetParameter(8));
  // f1_PhiBkg->SetParameters(f1_PhiMix->GetParameter(6), f1_PhiMix->GetParameter(7));
  f1_PhiBkg->SetParError(0,f1_PhiMix->GetParError(6));
  f1_PhiBkg->SetParError(1,f1_PhiMix->GetParError(7));
  f1_PhiBkg->SetParError(2,f1_PhiMix->GetParError(8));
//  f1_PhiBkg->SetParError(3,f1_PhiMix->GetParError(9));

  f1_PhiBkg->SetLineColor(1);
  f1_PhiBkg->SetLineStyle(2);
  f1_PhiBkg->SetRange(PhiMassFitLow,PhiMassFitHigh);


  double covmat[3][3];
  for(int icov =0; icov<3; icov++){
    for(int jcov =0; jcov<3; jcov++){
      covmat[icov][jcov]=fitResut->CovMatrix(6+icov,6+jcov);
    }
  }

  double *covmatArr=covmat[0];
	

  // double yieldBkg= f1_PhiBkg->Integral(PhiMassFitLow,PhiMassFitHigh)/ h_PhiMassData->GetBinWidth(1);
  double yieldBkg= f1_PhiBkg->Integral(PhiMassFitLow,PhiMassFitHigh)/ h_PhiMassData->GetBinWidth(1);

  cout<<"yield Bkg = " <<yieldBkg<<endl;

	// bkg sub method
 
	 double PhiMassMCMean=PhiMass;
 	double PhiMass2SigWidth=PhiMassCandWidth;		

	
	double N2Sig=h_PhiMassData->Integral(h_PhiMassData->FindBin(PhiMassMCMean-PhiMass2SigWidth),h_PhiMassData->FindBin(PhiMassMCMean+PhiMass2SigWidth));
	double binLow2SigVal=h_PhiMassData->GetBinLowEdge(h_PhiMassData->FindBin(PhiMassMCMean-PhiMass2SigWidth));
	double binHigh2SigVal=h_PhiMassData->GetBinLowEdge(h_PhiMassData->FindBin(PhiMassMCMean+PhiMass2SigWidth))+h_PhiMassData->GetBinWidth(1);
	// double Nbkg2Sig=f1_PhiBkg->Integral(PhiMassMCMean-PhiMass2SigWidth,PhiMassMCMean+PhiMass2SigWidth)/ h_PhiMassData->GetBinWidth(1); // wrong, the edge need to match with th1 bin
	double Nbkg2Sig=f1_PhiBkg->Integral(binLow2SigVal,binHigh2SigVal)/ h_PhiMassData->GetBinWidth(1);
	double Nbkg2SigError=f1_PhiBkg->IntegralError(binLow2SigVal,binHigh2SigVal,f1_PhiBkg->GetParameters(), covmatArr )/ h_PhiMassData->GetBinWidth(1);
	double NPhibkgSub=N2Sig-Nbkg2Sig;
	double NPhibkgSubError=sqrt(N2Sig+Nbkg2SigError*Nbkg2SigError);
	cout<<"N2Sig = "<<N2Sig<<" +- "<<sqrt(N2Sig)<<endl;
	cout<<"Nbkg2Sig = "<<Nbkg2Sig<<" +- "<<Nbkg2SigError<<endl;
	cout<<"NPhibkgSub = "<<NPhibkgSub<<" +- "<<NPhibkgSubError<<endl;


  NPhiInCut= f1_PhiSignal->Integral(DtktkResmassCutMean-DtktkResmassCutWidth, DtktkResmassCutMean+DtktkResmassCutWidth)/ h_PhiMassData->GetBinWidth(1);
  NPhiInCutErr = f1_PhiSignal->Integral(DtktkResmassCutMean-DtktkResmassCutWidth, DtktkResmassCutMean+DtktkResmassCutWidth)/h_PhiMassData->GetBinWidth(1) * f1_PhiSignal->GetParError(0)/f1_PhiSignal->GetParameter(0);

	NBkgInCut= f1_PhiBkg->Integral(DtktkResmassCutMean-DtktkResmassCutWidth, DtktkResmassCutMean+DtktkResmassCutWidth)/ h_PhiMassData->GetBinWidth(1); 
	NBkgInCutErr= f1_PhiBkg->IntegralError(DtktkResmassCutMean-DtktkResmassCutWidth, DtktkResmassCutMean+DtktkResmassCutWidth, f1_PhiBkg->GetParameters(), covmatArr)/ h_PhiMassData->GetBinWidth(1); 

	double NTotalInCut=f1_PhiMix->Integral(DtktkResmassCutMean-DtktkResmassCutWidth, DtktkResmassCutMean+DtktkResmassCutWidth)/ h_PhiMassData->GetBinWidth(1);

	double NPhiRatio=NPhiInCut/(NPhiInCut+NBkgInCut);
	double NPhiRatioErr=NPhiInCutErr/NTotalInCut;

	cout<<"NPhiRatio = "<<NPhiRatio<<" +- "<<NPhiRatioErr<<endl;


 
 // #<{(|
	TCanvas *c_DataFit=new TCanvas("c_DataFit","c_DataFit");
	c_DataFit->cd();
	h_PhiMassData->Draw();

	f1_PhiMix->SetRange(PhiMass-PhiMassCandWidth,PhiMass+PhiMassCandWidth);
	f1_PhiBkg->SetRange(PhiMass-PhiMassCandWidth,PhiMass+PhiMassCandWidth);
	f1_PhiSignal->SetRange(PhiMass-PhiMassCandWidth,PhiMass+PhiMassCandWidth);

	f1_PhiMix->Draw("same");
	f1_PhiBkg->Draw("same");
	f1_PhiSignal->Draw("same");

  shiftY=0;
	shiftX=0.03;
  TLatex *tl_binfitData =new TLatex();
//  tl_binfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s_Data_binfit",str_PbPb.Data()));  shiftY-=oneshift;
  tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s D_{S} %.0f < pt < %.0f ",s_ppPbPb.Data(),DptLowGL,DptHighGL)); shiftY-=oneshift;
  tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Phi Yield = %.1f #pm %.1f ",NPhiInCut,NPhiInCutErr ));  shiftY-=oneshift;
  tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Bkg Yield = %.1f #pm %.1f ",NBkgInCut,NBkgInCutErr ));  shiftY-=oneshift;
  tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Phi Ratio = %.3f #pm %.3f ",NPhiRatio,NPhiRatioErr ));  shiftY-=oneshift;
  // tl_binfitData->DrawLatexNDC(textpo+shiftX.08,textposy+shiftY,Form("BkgSub Yield = %.1f #pm %.1f ",NPhibkgSub,NPhibkgSubError ));  shiftY-=oneshift;
  // tl_binfitData->DrawLatexNDC(textpo+shiftX.28,textposy+shiftY,Form("Float width = %.3f",f1_PhiMix->GetParameter(5))); shiftY-=oneshift;

  if(s_ExtraName!="")
  {
  tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,s_ExtraName); shiftY-=oneshift;
  }
  // SavePlotDirs(c_binfit_Data[count_c_binfit],Form("%s_Data_binfit_pt%.0fto%.0f_%s",str_PbPb.Data(),DptLow,DptHigh,extraName.Data() ),{"SignalFit",str_PbPb,"binfit"} );

  c_DataFit->SaveAs(Form("./plots/fitMkk/%s_Data_Mkk_%s.png",s_type.Data(),s_ExtraName.Data()));


	 cout<<"f1MC->GetNpar() = "<<f1MC->GetNpar()<<endl;

	delete c_DataFit;
// |)}>#

	return f1_PhiMix;
}



TF1* fitPhiMass_cry2Gaus(TH1D *h_PhiMassData, TF1 *f1MC,double &NPhiInCut,double &NPhiInCutErr,double &NBkgInCut,double &NBkgInCutErr, TString s_ExtraName, int fixSigShape, TF1 *f1_default){

	cout<<"check -1"<<endl;

	cout<<"f1MC->GetParameter(1) = "<<f1MC->GetParameter(1)<<endl;
/*
	TCanvas *ctt=new TCanvas("ctt","ctt");
	ctt->cd();
	f1MC->Draw();
	ctt->SaveAs(Form("./plots/fitMkk/%s_Data_Mkk_%s_test.png",s_type.Data(),s_ExtraName.Data()));
*/
//	return f1MC;
	
  TF1 *f1_PhiSignal = (TF1*)f1MC->Clone("f1_PhiSignal");

  TF1 *f1_PhiBkg= new TF1("f1_PhiBkg","[0]*(1+[1]*x+[2]*(2*x*x-1))+[3]*(4*x*x*x-3*x)");
	TF1 *f1_PhiMix= new TF1("f1_PhiMix",crys_DoubleGaus_2ndCheb_fun,PhiMassFitLow, PhiMassFitHigh,12);

	for(int i =0; i<8 ; i++){
		f1_PhiMix->SetParameter(i,f1MC->GetParameter(i));
		if(i!=0){
    f1_PhiMix->FixParameter(i,f1MC->GetParameter(i));
		}
	}
	for(int i=8; i<12; i++){
		f1_PhiMix->SetParameter(i,0);
		f1_PhiMix->FixParameter(i,0);
	}

	// f1_PhiMix->SetParameter(0,2800);
	// f1_PhiMix->FixParameter(0,2800);
	// f1_PhiMix->FixParameter(8,0);
	// f1_PhiMix->SetParLimits();


	// f1_PhiMix->SetRange(PhiMassFitLow, PhiMassFitHigh);

	// TCanvas *ctt=new TCanvas("ctt","ctt");
	// ctt->cd();
	// f1_PhiMix->Draw();
/*
	TF1 *f1_cry2=new TF1("f1_cry2",crys_DoubleGaus_2ndCheb_fun,PhiMassFitLow,PhiMassFitHigh,11);
  f1_cry2->SetParameters(1000,1.019,0.002,-1,1,0.3,0.004,0, 0,0.0,0.0);
	f1_cry2->FixParameter(8,0);
	f1_cry2->FixParameter(9,0);
	f1_cry2->FixParameter(10,0);

  // f1_cry2->SetRange(1.01,1.03);
	f1_cry2->FixParameter(7,0);

	f1_cry2->Draw();
	cout<<"check 4"<<endl;
*/
	// ctt->SaveAs(Form("./plots/fitMkk/%s_Data_Mkk_%s_test.png",s_type.Data(),s_ExtraName.Data()));

	// cout<<"PhiMassFitLow = "<<PhiMassFitLow<<" , PhiMassFitHigh = "<<PhiMassFitHigh<<endl;

	// return f1MC;

 /*
	for(int i=0 ;i <f1MC->GetNpar(); i++){
		f1_PhiMix->SetParameter(i,f1MC->GetParameter(i));
		if(i>0 && i< f1MC->GetNpar()-1 ){
		f1_PhiMix->FixParameter(i,f1MC->GetParameter(i));
		}
	}
*/

  f1_PhiMix->SetLineColor(kRed);
//  f1_PhiMix->SetParameter(6,1000);

	cout<<"check 1"<<endl;

	// PhiMassFitLow=1.01;
	// PhiMassFitHigh=1.03;

  h_PhiMassData->Fit("f1_PhiMix","QN0","",    PhiMassFitLow, PhiMassFitHigh );
  h_PhiMassData->Fit("f1_PhiMix","QN0","",    PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L QN0","",  PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","LN0","",    PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L MN0","",  PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L MN0","",  PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L M IN0","",PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L M IN0","",PhiMassFitLow, PhiMassFitHigh);

/*
	if(fixSigShape!=1){
  // f1_PhiMix->ReleaseParameter(7);
  // f1_PhiMix->ReleaseParameter(1);
  f1_PhiMix->ReleaseParameter(3);
  f1_PhiMix->ReleaseParameter(4);
  // f1_PhiMix->ReleaseParameter(6);
	}
*/


  f1_PhiMix->ReleaseParameter(7);
	f1_PhiMix->SetParLimits(7,-0.6,0.6);
  // f1_PhiMix->ReleaseParameter(8);
  // f1_PhiMix->ReleaseParameter(9);
  // f1_PhiMix->ReleaseParameter(10);

  h_PhiMassData->Fit("f1_PhiMix","QN0","",    PhiMassFitLow, PhiMassFitHigh );
  h_PhiMassData->Fit("f1_PhiMix","QN0","",    PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L QN0","",  PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","LN0","",    PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L MN0","",  PhiMassFitLow, PhiMassFitHigh);

  f1_PhiMix->ReleaseParameter(1);
	f1_PhiMix->SetParLimits(1,1.018,1.021);


  h_PhiMassData->Fit("f1_PhiMix","L MN0","",  PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L M IN0","",PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L M I0","",PhiMassFitLow, PhiMassFitHigh);


  f1_PhiMix->ReleaseParameter(8);
	// f1_PhiMix->SetParLimits(8,0,1000000);
  f1_PhiMix->ReleaseParameter(9);
	f1_PhiMix->SetParLimits(9,-100,100);

  h_PhiMassData->Fit("f1_PhiMix","L MN0","",  PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L M IN0","",PhiMassFitLow, PhiMassFitHigh);

  f1_PhiMix->ReleaseParameter(10);
	f1_PhiMix->SetParLimits(10,-100,100);

  h_PhiMassData->Fit("f1_PhiMix","L MN0","",  PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L M IN0","",PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L M IN0","",PhiMassFitLow, PhiMassFitHigh);
  h_PhiMassData->Fit("f1_PhiMix","L M IN0","",PhiMassFitLow, PhiMassFitHigh);


	cout<<"check 2"<<endl;


  TFitResultPtr fitResut;
/*
	f1_PhiMix->SetRange(PhiMassFitLow,PhiMassFitHigh);
	int fitStatus=1;
  int fitIsValid=0;
   double fitPrecision=1.e-8;
  while(fitStatus>0 && fitStatus!=4000){
    TFitter::SetPrecision(fitPrecision);
    fitResut=h_PhiMassData->Fit("f1_PhiMix","L EMI S0","",PhiMassFitLow, PhiMassFitHigh);
    fitStatus=fitResut->Status();
    fitIsValid=fitResut->IsValid();
     cout<<"fit precision: "<<TFitter::GetPrecision()<<"   status: "<<fitStatus<<" isValid = "<< fitIsValid <<endl;
    if(fitStatus){
      fitPrecision *= 3;
    }
    if(fitPrecision>1.e-3) {break;}
  }
*/
	// fitResut=h_PhiMassData->Fit("f1_PhiMix","L EMI S0","",PhiMassFitLow, PhiMassFitHigh);
	fitResut=h_PhiMassData->Fit("f1_PhiMix","L MI S0","",PhiMassFitLow, PhiMassFitHigh);

 cout<<"f1_PhiMix->GetParameter(0) = "<<f1_PhiMix->GetParameter(0)<<endl;
cout<<"f1_PhiSignal->GetParameter(0) = "<<f1_PhiSignal->GetParameter(0)<<endl;
/* 
  for(int i=0 ;i <f1MC->GetNpar(); i++){
		f1MC->SetParameter(i, f1_PhiMix->GetParameter(i));
		f1MC->SetParError(i, f1_PhiMix->GetParError(i));
	}
*/
/*

	f1_PhiSignal->SetParameters(f1_PhiMix->GetParameter(0), f1_PhiMix->GetParameter(1) , f1_PhiMix->GetParameter(2) , f1_PhiMix->GetParameter(3) , f1_PhiMix->GetParameter(4), f1_PhiMix->GetParameter(5));
  f1_PhiSignal->SetParError(0,f1_PhiMix->GetParError(0));
  f1_PhiSignal->SetParError(1,f1_PhiMix->GetParError(1));
  f1_PhiSignal->SetParError(2,f1_PhiMix->GetParError(2));
  f1_PhiSignal->SetParError(3,f1_PhiMix->GetParError(3));
  f1_PhiSignal->SetParError(4,f1_PhiMix->GetParError(4));
  f1_PhiSignal->SetParError(5,f1_PhiMix->GetParError(5));
*/

	for(int i =0; i<8 ;i++){
		f1_PhiSignal->SetParameter(i,f1_PhiMix->GetParameter(i));
		f1_PhiSignal->SetParError(i,f1_PhiMix->GetParError(i));
	}

	cout<<"check 3"<<endl;

  f1_PhiSignal->SetFillColor(kOrange-3);
  f1_PhiSignal->SetFillStyle(3002);
  f1_PhiSignal->SetLineColor(kOrange-3);
  f1_PhiSignal->SetLineWidth(2);
  f1_PhiSignal->SetLineStyle(2);
  f1_PhiSignal->SetRange(PhiMassFitLow,PhiMassFitHigh); // must setrange before plot or get empty
  // double yieldPhiSignal= f1_PhiSignal->Integral(PhiMassFitLow,PhiMassFitHigh)/ h_PhiMassData->GetBinWidth(1);
  double yieldPhiSignal= f1_PhiSignal->Integral(PhiMass-PhiMassCandWidth,PhiMass+PhiMassCandWidth)/ h_PhiMassData->GetBinWidth(1);
  double yieldPhiSignalErr = f1_PhiSignal->Integral(PhiMass-PhiMassCandWidth,PhiMass+PhiMassCandWidth)/h_PhiMassData->GetBinWidth(1) * f1_PhiSignal->GetParError(0)/f1_PhiSignal->GetParameter(0);

   double yieldPhiSignalTry = f1_PhiSignal->GetParameter(0)/h_PhiMassData->GetBinWidth(1);
  double yieldPhiSignalTryErr = f1_PhiSignal->GetParError(0)/h_PhiMassData->GetBinWidth(1);


	cout<<"check 4"<<endl;

  cout<<"parameter 0 = "<<f1_PhiSignal->GetParameter(0)<<" BinWidth = "<< h_PhiMassData->GetBinWidth(1) <<endl;
  cout<<"yield PhiSignal = "<< yieldPhiSignal <<" +- "<<yieldPhiSignalErr<<endl;
  cout<<"yield PhiSignal from Para 0 = "<< yieldPhiSignalTry <<" +- "<<yieldPhiSignalTryErr<<endl;
/*
	for(int i=0; i<3; i++){	
	f1_PhiBkg->SetParameter(i,f1_PhiMix->GetParameter(i+Npar_Sig));
	f1_PhiBkg->SetParError(i,f1_PhiMix->GetParError(i+Npar_Sig));

	}
*/

  f1_PhiBkg->SetParameters(f1_PhiMix->GetParameter(8), f1_PhiMix->GetParameter(9) ,f1_PhiMix->GetParameter(10));
  // f1_PhiBkg->SetParameters(f1_PhiMix->GetParameter(6), f1_PhiMix->GetParameter(7));
  f1_PhiBkg->SetParError(0,f1_PhiMix->GetParError(8));
  f1_PhiBkg->SetParError(1,f1_PhiMix->GetParError(9));
  f1_PhiBkg->SetParError(2,f1_PhiMix->GetParError(10));
  f1_PhiBkg->SetParError(3,f1_PhiMix->GetParError(11));
//  f1_PhiBkg->SetParError(3,f1_PhiMix->GetParError(9));

  f1_PhiBkg->SetLineColor(1);
  f1_PhiBkg->SetLineStyle(2);
  f1_PhiBkg->SetRange(PhiMassFitLow,PhiMassFitHigh);


  double covmat[4][4];
  for(int icov =0; icov<4; icov++){
    for(int jcov =0; jcov<4; jcov++){
      covmat[icov][jcov]=fitResut->CovMatrix(8+icov,8+jcov);
    }
  }

  double *covmatArr=covmat[0];
	

  // double yieldBkg= f1_PhiBkg->Integral(PhiMassFitLow,PhiMassFitHigh)/ h_PhiMassData->GetBinWidth(1);
  double yieldBkg= f1_PhiBkg->Integral(PhiMassFitLow,PhiMassFitHigh)/ h_PhiMassData->GetBinWidth(1);
  
  cout<<"yield Bkg = " <<yieldBkg<<endl;

	// bkg sub method

	 double PhiMassMCMean=PhiMass;
 	double PhiMass2SigWidth=PhiMassCandWidth;		

	
	double N2Sig=h_PhiMassData->Integral(h_PhiMassData->FindBin(PhiMassMCMean-PhiMass2SigWidth),h_PhiMassData->FindBin(PhiMassMCMean+PhiMass2SigWidth));
	double binLow2SigVal=h_PhiMassData->GetBinLowEdge(h_PhiMassData->FindBin(PhiMassMCMean-PhiMass2SigWidth));
	double binHigh2SigVal=h_PhiMassData->GetBinLowEdge(h_PhiMassData->FindBin(PhiMassMCMean+PhiMass2SigWidth))+h_PhiMassData->GetBinWidth(1);
	// double Nbkg2Sig=f1_PhiBkg->Integral(PhiMassMCMean-PhiMass2SigWidth,PhiMassMCMean+PhiMass2SigWidth)/ h_PhiMassData->GetBinWidth(1); // wrong, the edge need to match with th1 bin
	double Nbkg2Sig=f1_PhiBkg->Integral(binLow2SigVal,binHigh2SigVal)/ h_PhiMassData->GetBinWidth(1);
	double Nbkg2SigError=f1_PhiBkg->IntegralError(binLow2SigVal,binHigh2SigVal,f1_PhiBkg->GetParameters(), covmatArr )/ h_PhiMassData->GetBinWidth(1);
	double NPhibkgSub=N2Sig-Nbkg2Sig;
	double NPhibkgSubError=sqrt(N2Sig+Nbkg2SigError*Nbkg2SigError);
	cout<<"N2Sig = "<<N2Sig<<" +- "<<sqrt(N2Sig)<<endl;
	cout<<"Nbkg2Sig = "<<Nbkg2Sig<<" +- "<<Nbkg2SigError<<endl;
	cout<<"NPhibkgSub = "<<NPhibkgSub<<" +- "<<NPhibkgSubError<<endl;


  NPhiInCut= f1_PhiSignal->Integral(DtktkResmassCutMean-DtktkResmassCutWidth, DtktkResmassCutMean+DtktkResmassCutWidth)/ h_PhiMassData->GetBinWidth(1);
  NPhiInCutErr = f1_PhiSignal->Integral(DtktkResmassCutMean-DtktkResmassCutWidth, DtktkResmassCutMean+DtktkResmassCutWidth)/h_PhiMassData->GetBinWidth(1) * f1_PhiSignal->GetParError(0)/f1_PhiSignal->GetParameter(0);

	NBkgInCut= f1_PhiBkg->Integral(DtktkResmassCutMean-DtktkResmassCutWidth, DtktkResmassCutMean+DtktkResmassCutWidth)/ h_PhiMassData->GetBinWidth(1); 
	NBkgInCutErr= f1_PhiBkg->IntegralError(DtktkResmassCutMean-DtktkResmassCutWidth, DtktkResmassCutMean+DtktkResmassCutWidth, f1_PhiBkg->GetParameters(), covmatArr)/ h_PhiMassData->GetBinWidth(1); 

	double NTotalInCut=f1_PhiMix->Integral(DtktkResmassCutMean-DtktkResmassCutWidth, DtktkResmassCutMean+DtktkResmassCutWidth)/ h_PhiMassData->GetBinWidth(1);
 
	double NPhiRatio=NPhiInCut/(NPhiInCut+NBkgInCut);
	double NPhiRatioErr=NPhiInCutErr/NTotalInCut;


	cout<<"NPhiInCut = "<<NPhiInCut<<endl;
	cout<<"NBkgInCut = "<<NBkgInCut<<endl;
 
	cout<<"NPhiRatio = "<<NPhiRatio<<" +- "<<NPhiRatioErr<<endl;



// #<{(|
	TCanvas *c_DataFit=new TCanvas("c_DataFit","c_DataFit");
	c_DataFit->cd();
	h_PhiMassData->Draw();

	f1_PhiMix->SetRange(PhiMassFitLow,PhiMassFitHigh);
	f1_PhiBkg->SetRange(PhiMassFitLow,PhiMassFitHigh);
	f1_PhiSignal->SetRange(PhiMassFitLow,PhiMassFitHigh);

	f1_PhiMix->Draw("same");
	f1_PhiBkg->Draw("same");
	f1_PhiSignal->Draw("same");

  shiftY=0;
	shiftX=-0.05;
	// textposx-=0.05;
  TLatex *tl_binfitData =new TLatex();
//  tl_binfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s_Data_binfit",str_PbPb.Data()));  shiftY-=oneshift;
   tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s D_{S} %.0f < pt < %.0f ",s_ppPbPb.Data(),DptLowGL,DptHighGL)); shiftY-=oneshift;
    tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Phi Yield = %.1f #pm %.1f ",NPhiInCut,NPhiInCutErr ));  shiftY-=oneshift;
  tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Bkg Yield = %.1f #pm %.1f ",NBkgInCut,NBkgInCutErr ));  shiftY-=oneshift;
  tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Phi Ratio = %.3f #pm %.3f ",NPhiRatio,NPhiRatioErr ));  shiftY-=oneshift;
  // tl_binfitData->DrawLatexNDC(textpo+shiftX.08,textposy+shiftY,Form("BkgSub Yield = %.1f #pm %.1f ",NPhibkgSub,NPhibkgSubError ));  shiftY-=oneshift;
  // tl_binfitData->DrawLatexNDC(textpo+shiftX.28,textposy+shiftY,Form("Float width = %.3f",f1_PhiMix->GetParameter(5))); shiftY-=oneshift;

  if(s_ExtraName!="")
  {
  tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,s_ExtraName); shiftY-=oneshift;
  }
   // SavePlotDirs(c_binfit_Data[count_c_binfit],Form("%s_Data_binfit_pt%.0fto%.0f_%s",str_PbPb.Data(),DptLow,DptHigh,extraName.Data() ),{"SignalFit",str_PbPb,"binfit"} );

  c_DataFit->SaveAs(Form("./plots/fitMkk/%s_Data_Mkk_%s_cry2Gauss.png",s_type.Data(),s_ExtraName.Data()));


	cout<<"f1MC->GetNpar() = "<<f1MC->GetNpar()<<endl;

	for(int i =0; i<8 ;i++){
//		f1_PhiSignal->SetParameter(i,f1_PhiMix->GetParameter(i));
	//	f1_PhiSignal->SetParError(i,f1_PhiMix->GetParError(i));

		cout<<"par "<<i<<" , Mix =  "<<f1_PhiMix->GetParameter(i)<<" , Sig = "<<f1_PhiSignal->GetParameter(i)<<endl;

	}

 


	delete c_DataFit;
// |)}>#

	return f1_PhiMix;
}

 
