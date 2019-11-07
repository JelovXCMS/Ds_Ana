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



int SameSignTest(){

	InitStyle();

	TFile *fin=TFile::Open("/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/DsMinTree_pp_Data_MBall_kpkptest.root");
	TTree *ntDs=(TTree*)fin->Get("ntDs");

	TH1D *h_DsMass=new TH1D("h_DsMass","h_DsMass",40,1.91,2.11);	h_DsMass->Sumw2();

	ntDs->Project("h_DsMass","Dmass","Dpt>3 && Dpt<20 && Dalpha<0.12 && Dchi2cl>0.03 && Ddls>2.0 && DtktkResmass>1.01051 && DtktkResmass<1.02851");

	gStyle->SetOptStat(0);

	TCanvas *c_DsMass= new TCanvas("c_DsMass","c_DsMass",800,800);
	c_DsMass->cd();
	h_DsMass->SetTitle("");
	h_DsMass->GetXaxis()->SetTitle("m_{KK#pi} (GeV)");
	h_DsMass->GetXaxis()->CenterTitle();
	h_DsMass->Draw();
	

	TLatex *tex_kpkp=new TLatex(0.55,0.72,"KK same sign bkg");
	tex_kpkp->SetNDC();
	tex_kpkp->Draw();

	SavePlotDirs(c_DsMass,"SameSignTest",{"Miscellaneous","SameSignTest"});


	return 0;

}


