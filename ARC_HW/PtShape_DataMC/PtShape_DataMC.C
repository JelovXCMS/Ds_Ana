
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting_simple.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"


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
#include <TGraph.h>
#include <TGraphAsymmErrors.h>


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

// Bool_t verbose=true;
Int_t StartBin=0;



int PtShape_DataMC(int startbin=1){

  setTDRStyle();
  InitStyle();

  gROOT->ForceStyle();
	gStyle->SetOptStat(0);

	TFile *fout= TFile::Open("./PtShape_DataMC.root","recreate");

	TFile *fin= TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/Results_RAA/output/RAA.root","read");
	TH1D *h_PromptDs_pp=(TH1D*)fin->Get("h_PromptDs_pp");
	TH1D *h_PromptDs_PbPb=(TH1D*)fin->Get("h_PromptDs_PbPb3");
	// this is 1/dpt , need to mutiply binwidth before normalization
	h_PromptDs_PbPb->GetXaxis()->SetRangeUser(6,40);

	MutiplyBinWidth(h_PromptDs_pp);
	MutiplyBinWidth(h_PromptDs_PbPb);

	TFile *fmc_pp=TFile::Open("/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/DsMinTree_MC_GenSampleMerge_pp_Prompt_phi.root","read");
	TFile *fmc_PbPb=TFile::Open("/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/DsMinTree_MC_GenSampleMerge_PbPb_Prompt_phi.root","read");
	
	TTree *tmc_pp=(TTree*)fmc_pp->Get("ntGen");
	TH1D *h_PromptDs_pp_mc=(TH1D*)h_PromptDs_pp->Clone("h_PromptDs_pp_mc");
	
	TTree *tmc_PbPb=(TTree*)fmc_PbPb->Get("ntGen");
	TH1D *h_PromptDs_PbPb_mc=(TH1D*)h_PromptDs_PbPb->Clone("h_PromptDs_PbPb_mc");


  TCut cutGenTrue=Form("GSignalType==1 && GcollisionId==0 && TMath::Abs(Gy)<1");
  TCut cutGenPNPrompt="GBAncestorpt<=0";
	TString Weightpp="weight*GptSampleWeight*GenFONLLWeight*GenD0DataWeight*GenDsDataWeight";	
	// TString Weightpp="weight*GptSampleWeight*GenFONLLWeight*GenD0DataWeight";	
	// TString Weightpp="weight*GptSampleWeight*GenFONLLWeight";	
	// TString Weightpp="weight*GptSampleWeight";	

	// TString WeightPbPb="weight*GptSampleWeight*PbPbVzWeight*Ncoll*GenFONLLWeight*GenD0DataWeight";
	TString WeightPbPb="weight*GptSampleWeight*PbPbVzWeight*Ncoll*GenFONLLWeight*GenD0DataWeight*GenDsDataWeight";


	tmc_pp->Project("h_PromptDs_pp_mc","Gpt",(TCut)((cutGenTrue && cutGenPNPrompt ) *Weightpp.Data() ));

	h_PromptDs_pp->Scale(1/h_PromptDs_pp->Integral());
	h_PromptDs_pp_mc->Scale(1/h_PromptDs_pp_mc->Integral());

	TH1D *h_PromptDs_pp_Ratio=(TH1D*)h_PromptDs_pp->Clone("h_PromptDs_pp_Ratio");
	h_PromptDs_pp_Ratio->Divide(h_PromptDs_pp_mc);

	cout<<" I1 = "<<h_PromptDs_pp->Integral()<<endl;
	cout<<" I2 = "<<h_PromptDs_pp_mc->Integral()<<endl;
	
  gROOT->ForceStyle();
	gStyle->SetOptStat(0);

	TCanvas *c_pp= new TCanvas("c_pp","",800,800);
	c_pp->cd();
	gPad->SetLogy();
	gPad->SetLogx();
	h_PromptDs_pp->GetXaxis()->SetTitle("D_{S}^{#pm} p_{T} (GeV/c)");
	h_PromptDs_pp->SetTitle("");
	h_PromptDs_pp->Draw();
	h_PromptDs_pp_mc->SetLineColor(2);
	h_PromptDs_pp_mc->Draw("same");

	TLatex *tlx=new TLatex();

	TLegend *le_pp=new TLegend(0.65,0.65,0.8,0.8);
	le_pp->SetBorderSize(0);
	le_pp->AddEntry(h_PromptDs_pp,"Data","l");
	le_pp->AddEntry(h_PromptDs_pp_mc,"MC","l");
	le_pp->Draw("same");

	textposx=0.65;
	textposy=0.85;
	tlx->SetTextSize(40);
	tlx->DrawLatexNDC(textposx,textposy,"pp");
	
	// c_pp->Modified();
	// c_pp->Update();

	gSystem->Exec("mkdir -p plots");

	c_pp->SaveAs("./plots/pp_ptshape_DataMC.png");
	
	// return 1;

	TCanvas *c_pp_ratio= new TCanvas("c_pp_ratio","",800,800);
	c_pp_ratio->cd();
	gPad->SetLogx();
	h_PromptDs_pp_Ratio->GetXaxis()->SetTitle("D_{S}^{#pm} p_{T} (GeV/c)");
	h_PromptDs_pp_Ratio->GetYaxis()->SetTitle("Data/MC");
	h_PromptDs_pp_Ratio->SetTitle("");
	h_PromptDs_pp_Ratio->Draw();

	tlx->DrawLatexNDC(textposx,textposy,"pp");


	TF1 *f1_pp=new TF1("f1_pp","[0]+[1]*x+[2]*x*x+[3]*x*x*x");
	// TF1 *f1_pp=new TF1("f1_pp","[0]+[1]*x+[2]*x*x");
	f1_pp->SetRange(2,40);
	h_PromptDs_pp_Ratio->Fit("f1_pp","LEMS0FI");
	h_PromptDs_pp_Ratio->Fit("f1_pp","LEMS0FI");
	h_PromptDs_pp_Ratio->Fit("f1_pp","LEMS0FI");
	// h_PromptDs_pp_Ratio->Fit("f1_pp","EMS0FI");
	// h_PromptDs_pp_Ratio->Fit("f1_pp","EMS0FI");
	// h_PromptDs_pp_Ratio->Fit("f1_pp","EMS0FI");
	f1_pp->SetLineColor(2);
	f1_pp->SetRange(2,40);
	// f1_pp->Draw("same");

	TF1 *f1_pp_log=new TF1("f1_pp_log","[0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*log(x)*log(x)*log(x)");
	// TF1 *f1_pp_log=new TF1("f1_pp_log","[0]+[1]*x+[2]*x*x");
	f1_pp_log->SetRange(2,40);
	f1_pp_log->FixParameter(3,0);
	h_PromptDs_pp_Ratio->Fit("f1_pp_log","LEMS0FI");
	h_PromptDs_pp_Ratio->Fit("f1_pp_log","LEMS0FI");
	h_PromptDs_pp_Ratio->Fit("f1_pp_log","LEMS0FI");
	// h_PromptDs_pp_Ratio->Fit("f1_pp_log","EMS0FI");
	// h_PromptDs_pp_Ratio->Fit("f1_pp_log","EMS0FI");
	// h_PromptDs_pp_Ratio->Fit("f1_pp_log","EMS0FI");
	f1_pp_log->SetLineColor(4);
	f1_pp_log->SetRange(2,40);
	f1_pp_log->Draw("same");



	c_pp_ratio->SaveAs("./plots/pp_ptshape_DataMCRatio.png");
	// return 1;

	tmc_PbPb->Project("h_PromptDs_PbPb_mc","Gpt",(TCut)((cutGenTrue && cutGenPNPrompt ) *WeightPbPb.Data() ));

	h_PromptDs_PbPb->Scale(1/h_PromptDs_PbPb->Integral());
	h_PromptDs_PbPb_mc->Scale(1/h_PromptDs_PbPb_mc->Integral());

	TH1D *h_PromptDs_PbPb_Ratio=(TH1D*)h_PromptDs_PbPb->Clone("h_PromptDs_PbPb_Ratio");
	h_PromptDs_PbPb_Ratio->Divide(h_PromptDs_PbPb_mc);

	cout<<" I1 = "<<h_PromptDs_PbPb->Integral()<<endl;
	cout<<" I2 = "<<h_PromptDs_PbPb_mc->Integral()<<endl;


	TCanvas *c_PbPb= new TCanvas("c_PbPb","",800,800);
	c_PbPb->cd();
	gPad->SetLogy();
	gPad->SetLogx();
	h_PromptDs_PbPb->GetXaxis()->SetTitle("D_{S}^{#pm} p_{T} (GeV/c)");
	h_PromptDs_PbPb->SetTitle("");
	h_PromptDs_PbPb->Draw();
	h_PromptDs_PbPb_mc->SetLineColor(2);
	h_PromptDs_PbPb_mc->Draw("same");


	TLegend *le_PbPb=new TLegend(0.65,0.65,0.8,0.8);
	le_PbPb->SetBorderSize(0);
	le_PbPb->AddEntry(h_PromptDs_PbPb,"Data","l");
	le_PbPb->AddEntry(h_PromptDs_PbPb_mc,"MC","l");
	le_PbPb->Draw("same");

	textposx=0.65;
	textposy=0.85;
	tlx->SetTextSize(40);
	tlx->DrawLatexNDC(textposx,textposy,"PbPb");
	
	// c_PbPb->Modified();
	// c_PbPb->Update();

	gSystem->Exec("mkdir -p plots");

	c_PbPb->SaveAs("./plots/PbPb_ptshape_DataMC.png");
	



	TCanvas *c_PbPb_ratio= new TCanvas("c_PbPb_ratio","",800,800);
	c_PbPb_ratio->cd();
	gPad->SetLogx();
	h_PromptDs_PbPb_Ratio->GetXaxis()->SetTitle("D_{S}^{#pm} p_{T} (GeV/c)");
	h_PromptDs_PbPb_Ratio->GetYaxis()->SetTitle("Data/MC");
	h_PromptDs_PbPb_Ratio->SetTitle("");
	h_PromptDs_PbPb_Ratio->Draw();
	tlx->DrawLatexNDC(textposx,textposy,"PbPb");
	// c_PbPb_ratio->SaveAs("./plots/PbPb_ptshape_DataMCRatio.png");


	// TF1 *f1_PbPb=new TF1("f1_PbPb","[0]+[1]*x+[2]*x*x+[3]*x*x*x");
	TF1 *f1_PbPb=new TF1("f1_PbPb","[0]+[1]*x");
	f1_PbPb->SetRange(6,40);
	// f1_PbPb->FixParameter(3,0);
	// f1_PbPb->FixParameter(2,0);
	h_PromptDs_PbPb_Ratio->Fit("f1_PbPb","EMS0F");
	h_PromptDs_PbPb_Ratio->Fit("f1_PbPb","EMS0F");
	h_PromptDs_PbPb_Ratio->Fit("f1_PbPb","EMS0F");
	f1_PbPb->SetLineColor(2);
	f1_PbPb->SetRange(6,40);
	// f1_PbPb->Draw("same");


	TF1 *f1_PbPb_log=new TF1("f1_PbPb_log","[0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*log(x)*log(x)*log(x)");
	f1_PbPb_log->SetRange(6,40);
	f1_PbPb_log->FixParameter(3,0);
	// f1_PbPb_log->FixParameter(2,0);
	h_PromptDs_PbPb_Ratio->Fit("f1_PbPb_log","EMS0F");
	h_PromptDs_PbPb_Ratio->Fit("f1_PbPb_log","EMS0F");
	h_PromptDs_PbPb_Ratio->Fit("f1_PbPb_log","EMS0F");
	f1_PbPb_log->SetLineColor(4);
	f1_PbPb_log->SetRange(6,40);
	f1_PbPb_log->Draw("same");




	c_PbPb_ratio->SaveAs("./plots/PbPb_ptshape_DataMCRatio.png");


/*
	// log x plot and fit
 	double bins_pt_pp_log[nbin_pt_pp+1];
	for(int i=0; i<nbin_pt_pp+1; i++){
		bins_pt_pp_log[i]=log(bins_pt_pp[i]);
		cout<<"bin "<<i<<" = "<<bins_pt_pp_log[i]<<endl;
	}
	TH1D *h_PromptDs_pp_Ratio_logx=new TH1D("h_PromptDs_pp_logx",";D_{S}^{#pm} p_{T} (GeV)", nbin_pt_pp,bins_pt_pp_log); 
	for(int i=0; i<nbin_pt_pp; i++){
		h_PromptDs_pp_Ratio_logx->SetBinContent(i+1,h_PromptDs_pp_Ratio->GetBinContent(i));
		h_PromptDs_pp_Ratio_logx->SetBinError(i+1,h_PromptDs_pp_Ratio->GetBinError(i));
	}

	TCanvas *c_pp_ratio_logx=new TCanvas("c_pp_ratio_logx","",800,800);
	c_pp_ratio_logx->cd();
	c_pp_ratio_logx->SetLogy();
	h_PromptDs_pp_Ratio_logx->Draw();
*/


	fout->cd();
	h_PromptDs_pp->Write();
	h_PromptDs_PbPb->Write();
	h_PromptDs_pp_mc->Write();
	h_PromptDs_PbPb_mc->Write();
	h_PromptDs_pp_Ratio->Write();
	h_PromptDs_PbPb_Ratio->Write();
	// f1_pp->Write();
	// f1_PbPb->Write();

	ofstream file_fit;
	file_fit.open("./PtShape_DataMC_fitPar.txt");	
	file_fit<<"pp"<<endl;
	file_fit<<"par0 : "<<f1_pp->GetParameter(0)<<endl;
	file_fit<<"par1 : "<<f1_pp->GetParameter(1)<<endl;
	file_fit<<"par2 : "<<f1_pp->GetParameter(2)<<endl;
	file_fit<<"par3 : "<<f1_pp->GetParameter(3)<<endl;
	file_fit<<"\nPbPb"<<endl;
	file_fit<<"par0 : "<<f1_PbPb->GetParameter(0)<<endl;
	file_fit<<"par1 : "<<f1_PbPb->GetParameter(1)<<endl;

	file_fit<<"\n\nlog(x) polynomials"<<endl;
	file_fit<<"pp"<<endl;
	file_fit<<"par0 : "<<f1_pp->GetParameter(0)<<endl;
	file_fit<<"par1 : "<<f1_pp->GetParameter(1)<<endl;
	file_fit<<"par2 : "<<f1_pp->GetParameter(2)<<endl;
	file_fit<<"par3 : "<<f1_pp->GetParameter(3)<<endl;
	file_fit<<"\nPbPb"<<endl;
	file_fit<<"par0 : "<<f1_PbPb->GetParameter(0)<<endl;
	file_fit<<"par1 : "<<f1_PbPb->GetParameter(1)<<endl;
	file_fit<<"par2 : "<<f1_PbPb->GetParameter(2)<<endl;



	return 0;

}

