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
#include <RooFitResult.h>
#include <RooWorkspace.h>
#include <RooMCStudy.h>
#include "RooHist.h"
#include "RooConstVar.h"



#include "TFitter.h"
#include "TFitResult.h"

using namespace RooFit;
using namespace std;


int plotBfr(){

	InitStyle();
	setTDRStyle();

	TFile *f5TeV=TFile::Open("bfraction5TeV.root");

	TH1D *h_B0_5TeV_frac=(TH1D*)f5TeV->Get("h_B0_5TeV_frac");
	TH1D *h_Bp_5TeV_frac=(TH1D*)f5TeV->Get("h_Bp_5TeV_frac");
	TH1D *h_Bs_5TeV_frac=(TH1D*)f5TeV->Get("h_Bs_5TeV_frac");
	
	TFile *f7TeV=TFile::Open("bfraction7TeV.root");

	TH1D *h_B0_7TeV_frac=(TH1D*)f7TeV->Get("h_B0_7TeV_frac");
	TH1D *h_Bp_7TeV_frac=(TH1D*)f7TeV->Get("h_Bp_7TeV_frac");
	TH1D *h_Bs_7TeV_frac=(TH1D*)f7TeV->Get("h_Bs_7TeV_frac");


	TCanvas *c_bfr=new TCanvas("c_bfr","c_bfr",c_wtopx,c_wtopy,c_W,c_H);	
	c_bfr->cd();

	gStyle->SetOptStat(0);
	gStyle->SetMarkerSize(0.4);
	gStyle->SetTextFont(42);

	gStyle->SetLegendFont(42);

  c_bfr->cd();
  c_bfr->SetFillColor(0);
  c_bfr->SetBorderMode(0);
  c_bfr->SetFrameFillStyle(0);
  c_bfr->SetFrameBorderMode(0);
  c_bfr->SetLeftMargin( c_Lmg );
  c_bfr->SetRightMargin( c_Rmg );
  c_bfr->SetTopMargin( c_Tmg );
  c_bfr->SetBottomMargin( c_Bmg );
  c_bfr->SetTickx(0);
  c_bfr->SetTicky(0);

	h_B0_5TeV_frac->SetMaximum(0.9);
	h_B0_5TeV_frac->SetMinimum(0);
	h_B0_5TeV_frac->SetTitle("");
	h_B0_5TeV_frac->GetXaxis()->SetTitle("B Rapidity");
	h_B0_5TeV_frac->GetXaxis()->CenterTitle();
	h_B0_5TeV_frac->GetYaxis()->SetTitle("fraction");
	h_B0_5TeV_frac->GetYaxis()->CenterTitle();

	
	double MarkerSize=0.6;

	h_B0_5TeV_frac->Draw();
	h_B0_5TeV_frac->SetMarkerStyle(21);
	h_B0_5TeV_frac->SetMarkerSize(MarkerSize);
	h_Bp_5TeV_frac->SetMarkerColor(kBlack);
	h_Bp_5TeV_frac->SetLineColor(kBlack);
	h_B0_7TeV_frac->SetMarkerStyle(26);
	h_B0_7TeV_frac->SetMarkerSize(MarkerSize);
	h_Bp_7TeV_frac->SetMarkerColor(kBlack-2);
	h_Bp_7TeV_frac->SetLineColor(kBlack-2);
	h_B0_7TeV_frac->Draw("SAME");

	h_Bp_5TeV_frac->Draw("SAME");
	h_Bp_5TeV_frac->SetMarkerStyle(21);
	h_Bp_5TeV_frac->SetMarkerSize(MarkerSize);
	h_Bp_5TeV_frac->SetMarkerColor(kRed);
	h_Bp_5TeV_frac->SetLineColor(kRed);
	h_Bp_7TeV_frac->SetMarkerStyle(26);
	h_Bp_7TeV_frac->SetMarkerSize(MarkerSize);
	h_Bp_7TeV_frac->SetMarkerColor(kRed-2);
	h_Bp_7TeV_frac->SetLineColor(kRed-2);
	h_Bp_7TeV_frac->Draw("SAME");

	h_Bs_5TeV_frac->Draw("SAME");
	h_Bs_5TeV_frac->SetMarkerStyle(21);
	h_Bs_5TeV_frac->SetMarkerSize(MarkerSize);
	h_Bs_5TeV_frac->SetMarkerColor(kBlue);
	h_Bs_5TeV_frac->SetLineColor(kBlue);
	h_Bs_7TeV_frac->SetMarkerStyle(26);
	h_Bs_7TeV_frac->SetMarkerSize(MarkerSize);
	h_Bs_7TeV_frac->SetMarkerColor(kBlue-2);
	h_Bs_7TeV_frac->SetLineColor(kBlue-2);
	h_Bs_7TeV_frac->Draw("SAME");


	TLatex *tl=new TLatex();
	tl->DrawLatexNDC(0.3,0.85,"B mesons production fraction (PYTHIA)");

	TLegend *le_5TeV=new TLegend(0.2,0.65,0.4,0.80);
	le_5TeV->SetBorderSize(0);
	le_5TeV->AddEntry(h_B0_5TeV_frac,"B^{0} 5TeV","p");
	le_5TeV->AddEntry(h_Bp_5TeV_frac,"B^{+} 5TeV","p");
	le_5TeV->AddEntry(h_Bs_5TeV_frac,"B_{S} 5TeV","p");
	le_5TeV->Draw("same");

	TLegend *le_7TeV=new TLegend(0.6,0.65,0.8,0.80);
	le_7TeV->SetBorderSize(0);
	le_7TeV->AddEntry(h_B0_7TeV_frac,"B^{0} 7TeV","p");
	le_7TeV->AddEntry(h_Bp_7TeV_frac,"B^{+} 7TeV","p");
	le_7TeV->AddEntry(h_Bs_7TeV_frac,"B_{S} 7TeV","p");
	le_7TeV->Draw("same");


	SavePlotDirs(c_bfr,"bfractions",{"Miscellaneous","bfraction"});


return 0;
}
