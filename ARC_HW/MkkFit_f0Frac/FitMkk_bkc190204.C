#include "MkkFit_para.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"

#include "TFitter.h"
#include "TFitResult.h"
#include <cmath>

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


double crys_DoubleGaus_fun(const double *x, const double *p){

	return (p[0] * crys_DoubleGaus_fun(x[0], p[3], p[4], p[2], p[1],p[5],p[6],p[7]));

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
	double PhiMassFitLow=1.01951-0.02;
	 double PhiMassFitHigh=1.01951-0.02;


int FitMkk(int isPbPb=0, double DptLow=2, double DptHigh=40){

   DptLowGL=DptLow;
  DptHighGL=DptHighGL;

	double phiFitLow=0.99;
	double phiFitHigh=1.04;

  gSystem->Exec("mkdir -p plots/fitMkk");

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
	TH1D *h_NDs_mkkbin_BkgSub   =(TH1D*)fMC->Get("h_NDs_mkkbin_BkgSub");

// #<{(|
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
// |)}>#
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

	TF1 *f1_cry2=new TF1("f1_cry2",crys_DoubleGaus_fun,phiFitLow,phiFitHigh,8);
  f1_cry2->SetParameters(1000,1.019,0.002,-1,1,0.3,0.004,0);
  // f1_cry2->SetRange(1.01,1.03);
	f1_cry2->FixParameter(7,0);
// #<{(|
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","QNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","QNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","LQNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L MQNO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M INO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M IO","",phiFitLow,phiFitHigh);
  h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M IO","",PhiMass-0.02,PhiMass+0.02);
  // h_mkk_MCP_Phi_fineBin->Fit("f1_cry2","L M IO","",1.01,1.03);
// |)}>#
 TCanvas *c_cry2=new TCanvas("c_cry2","c_cry2");
  c_cry2->cd();
  f1_cry2->SetRange(PhiMass-0.02,PhiMass+0.02);
  h_mkk_MCP_Phi_fineBin->Draw();
  // f1_cry2->SetRange(1.01,1.03);
  // f1_cry2->SetRange(PhiMass-0.02,PhiMass+0.02);
  f1_cry2->Draw("same");
  // f1_cry2->Draw();

  c_cry2->SaveAs(Form("./plots/fitMkk/%s_MC_Mkk_cry2Gaus.png",s_type.Data()));


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

  TF1 *f1Data=fitPhiMass(h_NDs_mkkbin_Fit,f1_PhiMC, NPhiInCut, NPhiInCutErr,NBkgInCut ,NBkgInCutErr);
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



  // TF1 *f1_PhiSignal = (TF1*)f1_PhiMC->Clone("f1_PhiSignal");
  TF1 *f1_PhiSignal = (TF1*)f1MC->Clone("f1_PhiSignal");
  TF1 *f1_PhiBkg= new TF1("f1_PhiBkg","[0]*(1+[1]*x+[2]*(2*x*x-1))");
  TF1 *f1_PhiBkg_clone=(TF1*)f1_PhiBkg->Clone("f1_PhiBkg_clone");  // stupid bug in root, must use clone function
  TF1 *f1_PhiMix = new TF1("f1_PhiMix","f1_PhiSignal+f1_PhiBkg");

	for(int i=0 ;i <f1MC->GetNpar(); i++){
		f1_PhiMix->SetParameter(i,f1_PhiMC->GetParameter(i));
		if(i>0 && )
	}
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

