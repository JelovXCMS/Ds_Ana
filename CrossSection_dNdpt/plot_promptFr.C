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

using namespace std;

  double shiftY=0.05;
  double shiftX=0.32;
  double oneshift=0.070;




int plot_promptFr(){

	InitStyle();
	setTDRStyle();

	TFile *fpp=TFile::Open("./output/PromptDsCrossSectionPP_FixShape.root");
	TFile *fPbPb3=TFile::Open("./output/PromptDsdNdptPbPb3_FixShape.root");	

	TH1D *h_CSfr_PromptDs_pp=(TH1D*)fpp->Get("h_CSfr_PromptDs_pp");
	TH1D *h_CSfr_PromptDs_PbPb3=(TH1D*)fPbPb3->Get("h_CSfr_PromptDs_PbPb3");

	h_CSfr_PromptDs_PbPb3->SetBinContent(2,0);
	h_CSfr_PromptDs_PbPb3->SetBinError(2,0);

	TCanvas *c_fr=new TCanvas("c_fr","c_fr",c_wtopx,c_wtopy,c_W,c_H);
	c_fr->cd();

	gPad->SetLogx();
	gStyle->SetOptStat(0);

	c_fr->SetTopMargin(c_Tmg);
	c_fr->SetBottomMargin(c_Bmg);
	c_fr->SetLeftMargin(c_Lmg);
	c_fr->SetRightMargin(c_Rmg);

	h_CSfr_PromptDs_pp->SetMaximum(1);
	h_CSfr_PromptDs_pp->SetMinimum(0);
	h_CSfr_PromptDs_pp->SetTitle("");
	h_CSfr_PromptDs_pp->GetXaxis()->SetTitle("D_{S} p_{T} (GeV)");
	h_CSfr_PromptDs_pp->GetXaxis()->SetTitleOffset(0.8);
	h_CSfr_PromptDs_pp->GetYaxis()->SetTitle("Prompt D_{S} fraction");
	h_CSfr_PromptDs_pp->GetYaxis()->SetTitleOffset(1.1);
	h_CSfr_PromptDs_pp->Draw();
	h_CSfr_PromptDs_pp->SetLineColor(2);
	h_CSfr_PromptDs_pp->SetMarkerColor(2);	

	h_CSfr_PromptDs_PbPb3->SetLineColor(4);
	h_CSfr_PromptDs_PbPb3->SetMarkerColor(4);
	h_CSfr_PromptDs_PbPb3->Draw("same");

	TLegend *le=new TLegend(0.65,0.20,0.85,0.40);
	le->SetBorderSize(0);
	le->AddEntry(h_CSfr_PromptDs_pp,"pp","l");
	le->AddEntry(h_CSfr_PromptDs_PbPb3,"PbPb","l");

	le->Draw("same");

	lumiTextSize=0.45;
	writeExtraText=1;


	CMS_lumi(c_fr,3,0);


	SavePlotDirs(c_fr,"PromptFr",{"CSdNdpt_Result"});


	return 0;

}
