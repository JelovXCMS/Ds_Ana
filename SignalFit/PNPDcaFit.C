#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"

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
#include "TFitter.h"
#include "TFitResult.h"



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

using namespace RooFit;
using namespace std;

double shiftY=0;
double oneshift=0.075;

TH1F *h_MCP_DCA;
TH1F *h_MCNP_DCA;


void normalize(TH1F* h)
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


Double_t funMix(Double_t* x_, Double_t* para)
{
	float x = x_[0];
	float APrompt = para[0];
	float ANonPrompt = para[1];
	float promptYield = 0;
	float nonPromptYield = 0;

	promptYield = h_MCP_DCA->GetBinContent(h_MCP_DCA->GetXaxis()->FindBin(x));
	nonPromptYield = h_MCNP_DCA->GetBinContent(h_MCNP_DCA->GetXaxis()->FindBin(x));

	return APrompt*promptYield+ANonPrompt*nonPromptYield;
}

/*
Double_t funNonPrompt(Double_t* x_, Double_t* para)
{
	float x = x_[0];
	float APrompt = para[0];
	float ANonPrompt = para[1];
	float nonPromptYield = 0;
	nonPromptYield = hD0DcaMCNPSignal->GetBinContent(hD0DcaMCNPSignal->GetXaxis()->FindBin(x));
	return ANonPrompt*nonPromptYield;
}
*/

TF1 *DcaFitFun(TH1F *h_Data, TH1F *h_MCP, TH1F *h_MCNP, TH1F *h_frPrompt){

	h_MCP_DCA=h_MCP;
	h_MCNP_DCA=h_MCNP;

	TCanvas *c_test= new TCanvas("c_test","c_test",800,800);
	c_test->cd();
	h_MCP->Draw("PMC PLC");
	h_MCNP->Draw("SAME PMC PLC");
	gPad->BuildLegend();

  double integralTotalYield = h_Data->Integral(1,h_Data->GetXaxis()->GetNbins(),"width");
  cout<<"integralTotalYield: "<<integralTotalYield<<endl;

	TF1* fMix = new TF1("fMix",&funMix, 0., 0.5, 2);
  fMix->SetParameters(0.5*integralTotalYield,0.5*integralTotalYield);
  fMix->SetParLimits(0,0,2*integralTotalYield);
  fMix->SetParLimits(1,0,2*integralTotalYield);

  fMix->SetLineColor(2);
  fMix->SetFillColor(kRed-9);
  fMix->SetFillStyle(1001);

	float fitRangeL=0;
	float fitRangeH=0.07; // might need to change

	int fitStatus = 1;
	TFitResultPtr fitResult;
	double fitPrecision = 1.e-6;


	while(fitStatus)
	{
		TFitter::SetPrecision(fitPrecision);
		fMix->SetParameters(0.5*integralTotalYield,0.5*integralTotalYield);
		fMix->SetParError(0,0.1*integralTotalYield);
		fMix->SetParError(1,0.1*integralTotalYield);
		fitResult = h_Data->Fit("fMix","E SNQ0", "", fitRangeL, fitRangeH);
		fitStatus = fitResult->Status();
		cout<<"fit precision: "<<TFitter::GetPrecision()<<"   status: "<<fitStatus<<endl;
		if(fitStatus)
			fitPrecision *= 10;
	}

  fMix->SetParameters(integralTotalYield,0.9);
	fMix->SetParError(0,0.1*integralTotalYield);
	fMix->SetParError(1,0.1*integralTotalYield);

	fitResult = h_Data->Fit("fMix","E S0", "", fitRangeL, fitRangeH);

	TCanvas *c_dcaFit=new TCanvas("c_dcaFit","c_dcaFit",800,800);
	c_dcaFit->cd();
	cout<<"check 1"<<endl;
	h_Data->Draw("same");
	cout<<"check 2"<<endl;
	h_Data->GetFunction("fMix")->Draw("flsame");
	cout<<"check 3"<<endl;
	h_Data->Draw("same");

	cout<<"check 4"<<endl;
  float promptDYield = fMix->GetParameter(0);
  float promptDYieldErrorDataOnly = fMix->GetParError(0);
	cout<<"check 5"<<endl;
  float bToDYield = fMix->GetParameter(1);
  float bToDYieldErrorDataOnly = fMix->GetParError(1);
  float totalYield = promptDYield+bToDYield;
  float promptFraction = promptDYield/totalYield;
	cout<<"promptFraction = "<<promptFraction<<endl;

	h_frPrompt->SetBinContent(1,promptFraction);
	return h_Data->GetFunction("fMix");


}



int PNPDcaFit(int isPbPb=0, int ibin_Dpt=4, int doBkgVariation=0, int doSigVariation=0, int doCutScan=0, int doMCstudy=0){

	cout<<"isPbPb = "<<isPbPb<<" , ibin_Dpt = "<<ibin_Dpt<<" , doBkgVariation = "<<doBkgVariation<<" , doSigVariation = "<<doSigVariation<<" , doCutScan = "<<doCutScan<<" , doMCstudy = "<<doMCstudy<<endl;

	InitStyle();
	initParameter();

	int nbin_pt=nbin_pt_pp;
	double *bins_pt=bins_pt_pp;

	double DptLow  = bins_pt[ibin_Dpt];
	double DptHigh = bins_pt[ibin_Dpt+1];

	TString mcName_Prompt="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output/MC_eff_pp_Prompt_phikkpi.root";
	TString mcName_NonPrompt="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output/MC_eff_pp_NonPrompt_phikkpi.root";
	TString dataDcaName=Form("output/FitResult_pp_pt%.0fto%.0f.root",DptLow,DptHigh);
	TString outName=Form("output/DcaFit/PNPDcaFit_pp_pt%.0fto%.0f.root",DptLow,DptHigh);
	TString str_PbPb="pp";


	if(isPbPb==3){
		nbin_pt=nbin_pt_PbPb3;
		bins_pt=bins_pt_PbPb3;

		DptLow  = bins_pt[ibin_Dpt];
		DptHigh = bins_pt[ibin_Dpt+1];

		mcName_Prompt="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output/MC_eff_PbPb3_Prompt_phikkpi.root";
		mcName_NonPrompt="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output/MC_eff_PbPb3_NonPrompt_phikkpi.root";
		dataDcaName=Form("output/FitResult_PbPb3_pt%.0fto%.0f.root",DptLow,DptHigh);
		outName=Form("output/DcaFit/PNPDca_PbPb3_pt%.0fto%.0f.root",DptLow,DptHigh);

		str_PbPb="PbPb";
	}

	TFile *f_data=TFile::Open(dataDcaName.Data(),"read");
	TFile *f_mc_Prompt=TFile::Open(mcName_Prompt.Data(),"read");
	TFile *f_mc_NonPrompt=TFile::Open(mcName_NonPrompt.Data(),"read");

	TFile *f_out=TFile::Open(outName.Data(),"RECREATE");
	f_out->cd();

	TH1F *h_dataDca=(TH1F*)f_data->Get("h_YieldDcabin");

	h_dataDca->Draw();

	TTree *t_mc_Prompt=(TTree*)f_mc_Prompt->Get(Form("t_DsMass_pt%.0fto%.0f",DptLow,DptHigh));
	TH1F *h_mc_Prompt_Dca=new TH1F("h_mc_Prompt_Dca","h_mc_Prompt_Dca",nbin_dca,bins_dca); 
	h_mc_Prompt_Dca->Sumw2();
	t_mc_Prompt->Project("h_mc_Prompt_Dca","Ddca","BasicWeight");
	// h_mc_Prompt_Dca->Scale(1/h_mc_Prompt_Dca->Integral()); // must normalize
	normalize(h_mc_Prompt_Dca);

	TTree *t_mc_NonPrompt=(TTree*)f_mc_NonPrompt->Get(Form("t_DsMass_pt%.0fto%.0f",DptLow,DptHigh));
	TH1F *h_mc_NonPrompt_Dca=new TH1F("h_mc_NonPrompt_Dca","h_mc_NonPrompt_Dca",nbin_dca,bins_dca); 
	h_mc_NonPrompt_Dca->Sumw2();
	t_mc_NonPrompt->Project("h_mc_NonPrompt_Dca","Ddca","BasicWeight");
	h_mc_NonPrompt_Dca->Scale(1/h_mc_NonPrompt_Dca->Integral()); // must normalize
	normalize(h_mc_NonPrompt_Dca);


	f_out->cd();
	TH1F *h_frPrompt=new TH1F("h_frPrompt","h_frPrompt",1,0,1); h_frPrompt->Sumw2();


	// TF1 *DcaFitFun(TH1F *h_Data, TH1F *h_MCP, TH1F *h_MCNP, TH1F *fr_Prompt){
	TF1 *f1_DcaFit=DcaFitFun(h_dataDca,h_mc_Prompt_Dca,h_mc_NonPrompt_Dca, h_frPrompt);


	f_out->cd();
	h_frPrompt->Write("",TObject::kOverwrite);


	return 0;
}
