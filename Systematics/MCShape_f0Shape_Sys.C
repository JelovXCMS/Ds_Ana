#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"


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
#include <TLine.h>


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

// using namespace RooFit;
using namespace std;

// Int_t StartBin=0;

TCanvas *cav[500];
Int_t c_count=0;

double shiftY=0;
double oneshift=0.075;


int MCShape_f0Shape_Sys(Int_t isPbPb=3, Int_t startbin=0, Int_t start_var=0){
  if(isPbPb==3){startbin=2;}
  // StartBin=startbin;

  InitStyle();
  initParameter();

	gStyle->SetOptStat(0);

  TString str_PbPb="pp";
	TString str_PbPbtext="pp";
  int nbin_pt=nbin_pt_pp;
  double *bins_pt=bins_pt_pp;

  TString outfile="output/MCShapeSys_pp.root";
  // TString infile="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/CrossSection_dNdpt/output/PromptDsCrossSectionPP.root";
  TString infile="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/CrossSection_dNdpt/output/PromptDsCrossSectionPP_FixShape.root";

  double *DalphaMaxScan_bins=DalphaMaxScan_bins_pp;
  double *DdlsMinScan_bins=DdlsMinScan_bins_pp;
  double *Dchi2clMinScan_bins=Dchi2clMinScan_bins_pp;

  double *bins_DalphaMaxScan=bins_DalphaMaxScan_pp;
  double *bins_DdlsMinScan=bins_DdlsMinScan_pp;
  double *bins_Dchi2clMinScan=bins_Dchi2clMinScan_pp;

	double *DdlsMin_bins=DdlsMin_bins_pp;
	double *DalphaMax_bins=DalphaMax_bins_pp;
	double *Dchi2clMin_bins=Dchi2clMin_bins_pp;

  if(isPbPb==3){
    str_PbPb="PbPb3";
		str_PbPbtext="PbPb";
    nbin_pt=nbin_pt_PbPb3;
    bins_pt=bins_pt_PbPb3;
    outfile="output/MCShapeSys_PbPb3.root";
    // infile="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/CrossSection_dNdpt/output/PromptDsdNdptPbPb3.root";
    infile="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/CrossSection_dNdpt/output/PromptDsdNdptPbPb3_FixShape.root";

    DalphaMaxScan_bins=DalphaMaxScan_bins_PbPb3;
    DdlsMinScan_bins=DdlsMinScan_bins_PbPb3;
    Dchi2clMinScan_bins=Dchi2clMinScan_bins_PbPb3;

    bins_DalphaMaxScan=bins_DalphaMaxScan_PbPb3;
    bins_DdlsMinScan=bins_DdlsMinScan_PbPb3;
    bins_Dchi2clMinScan=bins_Dchi2clMinScan_PbPb3;

		DdlsMin_bins=DdlsMin_bins_PbPb3;
		DalphaMax_bins=DalphaMax_bins_PbPb3;
		Dchi2clMin_bins=Dchi2clMin_bins_PbPb3;

  }


	cout<<"start marco"<<endl;
	TFile *fin=TFile::Open(infile.Data(),"READ");
	TFile *fout=TFile::Open(outfile.Data(),"RECREATE");
	fout->cd();

	// saving results
	TH1D *h_PromptDs_MCShape_SysRel=new TH1D("h_PromptDs_MCShape_SysRel","h_PromptDs_MCShape_SysRel",nbin_pt,bins_pt);
	h_PromptDs_MCShape_SysRel->Sumw2();


	TH1D *h_PromptDs_f0Eff_SysRel=new TH1D("h_PromptDs_f0Eff_SysRel","h_PromptDs_f0Eff_SysRel",nbin_pt,bins_pt);
	h_PromptDs_f0Eff_SysRel->Sumw2();
	
	TH1D *h_PromptDs_PhiRatio_SysRel=new TH1D("h_PromptDs_PhiRatio_SysRel","h_PromptDs_PhiRatio_SysRel",nbin_pt,bins_pt);
	h_PromptDs_PhiRatio_SysRel->Sumw2();


	TLatex *tlt=new TLatex();
	cout<<"reading file done "<<endl;

	TLine *tl=new TLine(0.02,1,0.22,1);
	tl->SetLineColor(2);
	tl->SetLineStyle(2);

	TH1D *h_PromptDs=(TH1D*)fin->Get(Form("h_PromptDs_%s",str_PbPb.Data())); // central value

	TH1D *h_PromptDs_MCS_Pythia=(TH1D*)fin->Get(Form("h_PromptDs_MCS_Pythia%s",str_PbPb.Data()));
	TH1D *h_PromptDs_MCS_FONLL=(TH1D*)fin->Get(Form("h_PromptDs_MCS_FONLL%s",str_PbPb.Data()));

	TH1D *h_PromptDs_f0Effup=(TH1D*)fin->Get(Form("h_PromptDs_f0Effup%s",str_PbPb.Data()));
	TH1D *h_PromptDs_f0Effdown=(TH1D*)fin->Get(Form("h_PromptDs_f0Effdown%s",str_PbPb.Data()));


	double CSdNdpt=0;
	double CSdNdpt_Pythia=0;
	double CSdNdpt_FONLL=0;
	double CSdNdpt_Max=0;
	double CSdNdptRel=0;

	double CSdNdpt_f0Effup=0;
	double CSdNdpt_f0Effdown=0;
	double CSdNdptRel_f0Eff=0;

	double PhiRatioErr_pp=0.014/0.954;
	double PhiRatioErr_PbPb=0.009/0.896;

	double PhiRatioErr=PhiRatioErr_pp;
	if(isPbPb){
		PhiRatioErr=PhiRatioErr_PbPb;
	}

	for(int i=startbin; i<nbin_pt; i++){
		CSdNdpt=h_PromptDs->GetBinContent(i+1);
		CSdNdpt_Pythia=h_PromptDs_MCS_Pythia->GetBinContent(i+1);
		CSdNdpt_FONLL=h_PromptDs_MCS_FONLL->GetBinContent(i+1);
		// CSdNdpt_Max=abs(CSdNdpt - CSdNdpt_Pythia) > abs(CSdNdpt - CSdNdpt_FONLL) ? 	abs(CSdNdpt - CSdNdpt_Pythia) : abs(CSdNdpt - CSdNdpt_FONLL);
		CSdNdpt_Max= abs(CSdNdpt - CSdNdpt_FONLL);
		CSdNdptRel=CSdNdpt_Max/CSdNdpt;
		h_PromptDs_MCShape_SysRel->SetBinContent(i+1,CSdNdptRel);

		cout<<setprecision(1)<<std::fixed;

		cout<<" \n \n --- pt : "<<bins_pt[i]<<" - "<<bins_pt[i+1]<<endl;
		cout<<"CSdNdptRel_MCptShape = "<<CSdNdptRel*100<<"% , CSdNdpt = "<<CSdNdpt<<" , CSdNdpt_Pythia = "<<CSdNdpt_Pythia<<" CSdNdpt_FONLL = "<<CSdNdpt_FONLL<<endl;

		CSdNdpt_f0Effup=h_PromptDs_f0Effup->GetBinContent(i+1);
		CSdNdpt_f0Effdown=h_PromptDs_f0Effdown->GetBinContent(i+1);

		CSdNdptRel_f0Eff= abs(CSdNdpt_f0Effup - CSdNdpt )> abs(CSdNdpt_f0Effdown - CSdNdpt ) ? abs(CSdNdpt_f0Effup - CSdNdpt )/CSdNdpt : abs(CSdNdpt_f0Effdown - CSdNdpt )/CSdNdpt; 

		h_PromptDs_f0Eff_SysRel->SetBinContent(i+1, CSdNdptRel_f0Eff);
		cout<<"CSdNdptRel_f0Eff : up = "<< abs(CSdNdpt_f0Effup - CSdNdpt )/CSdNdpt*100  <<" % , down = "<<abs(CSdNdpt_f0Effdown - CSdNdpt )/CSdNdpt*100<<"%"<<endl;

	
		h_PromptDs_PhiRatio_SysRel->SetBinContent(i+1,PhiRatioErr);

	}

	fout->cd();
	h_PromptDs_MCShape_SysRel->Write("",TObject::kOverwrite);
	h_PromptDs_f0Eff_SysRel->Write("",TObject::kOverwrite);
	h_PromptDs_PhiRatio_SysRel->Write("",TObject::kOverwrite);

	return 0;
}
