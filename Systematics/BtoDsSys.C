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


int BtoDsSys(Int_t isPbPb=3, Int_t startbin=0, Int_t start_var=0){
  if(isPbPb==3){startbin=2;}
  // StartBin=startbin;

	// use Raa 1.5 factor or not


  InitStyle();
  initParameter();

	gStyle->SetOptStat(0);

  TString str_PbPb="pp";
	TString str_PbPbtext="pp";
  int nbin_pt=nbin_pt_pp;
  double *bins_pt=bins_pt_pp;

  TString outfile="output/BtoDsSys_pp.root";
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
    outfile="output/BtoDsSys_PbPb3.root";
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
	TH1D *h_PromptDs_BtoDs_SysRel=new TH1D("h_PromptDs_BtoDs_SysRel","h_PromptDs_BtoDs_SysRel",nbin_pt,bins_pt);
	h_PromptDs_BtoDs_SysRel->Sumw2();

	TH1D *h_PromptDs_BtoDs_SysRel_ForRaa=new TH1D("h_PromptDs_BtoDs_SysRel_ForRaa","h_PromptDs_BtoDs_SysRel_ForRaa",nbin_pt,bins_pt);
	h_PromptDs_BtoDs_SysRel_ForRaa->Sumw2();

	TH1D *h_PromptDs_BtoDs_SysRel_ForRaaCancel=new TH1D("h_PromptDs_BtoDs_SysRel_ForRaaCancel","h_PromptDs_BtoDs_SysRel_ForRaaCancel",nbin_pt,bins_pt);
	h_PromptDs_BtoDs_SysRel_ForRaaCancel->Sumw2();




	TLatex *tlt=new TLatex();
	cout<<"reading file done "<<endl;

	TLine *tl=new TLine(0.02,1,0.22,1);
	tl->SetLineColor(2);
	tl->SetLineStyle(2);

	TH1D *h_PromptDs=(TH1D*)fin->Get(Form("h_PromptDs_%s",str_PbPb.Data())); // central value

	TH1D *h_PromptDs_BtoDs_FONLL=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_FONLL%s",str_PbPb.Data()));
	TH1D *h_PromptDs_BtoDs_FONLL_SysRel=new TH1D("h_PromptDs_BtoDs_FONLL_SysRel","h_PromptDs_BtoDs_FONLL_SysRel",nbin_pt,bins_pt); h_PromptDs_BtoDs_FONLL_SysRel->Sumw2();

	TH1D *h_PromptDs_BtoDs_Raaup=NULL;
	TH1D *h_PromptDs_BtoDs_Raadown=NULL;

	if(isPbPb){
		h_PromptDs_BtoDs_Raaup=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_Raaup%s",str_PbPb.Data()));
		h_PromptDs_BtoDs_Raadown=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_Raadown%s",str_PbPb.Data()));
	}

	TH1D *h_PromptDs_BtoDs_FONLLShape=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_FONLLShape%s",str_PbPb.Data()));

	TH1D *h_PromptDs_BtoDs_Errup=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_Errup%s",str_PbPb.Data()));
	TH1D *h_PromptDs_BtoDs_Errdown=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_Errdown%s",str_PbPb.Data()));

	TH1D *h_PromptDs_BtoDs_ScaleBRup=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_ScaleBRup%s",str_PbPb.Data()));
	TH1D *h_PromptDs_BtoDs_ScaleBRdown=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_ScaleBRdown%s",str_PbPb.Data()));


	TH1D *h_PromptDs_BtoDs_ScaleBR_B0toDmax=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_ScaleBR_B0toDmax%s",str_PbPb.Data()));
	TH1D *h_PromptDs_BtoDs_ScaleBR_B0toDmin=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_ScaleBR_B0toDmin%s",str_PbPb.Data()));

	TH1D *h_PromptDs_BtoDs_ScaleBR_BptoDmax=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_ScaleBR_BptoDmax%s",str_PbPb.Data()));
	TH1D *h_PromptDs_BtoDs_ScaleBR_BptoDmin=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_ScaleBR_BptoDmin%s",str_PbPb.Data()));

	TH1D *h_PromptDs_BtoDs_ScaleBR_BstoDmax=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_ScaleBR_BstoDmax%s",str_PbPb.Data()));
	TH1D *h_PromptDs_BtoDs_ScaleBR_BstoDmin=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_ScaleBR_BstoDmin%s",str_PbPb.Data()));

	TH1D *h_PromptDs_BtoDs_ScaleBR_B0toDsmax=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_ScaleBR_B0toDsmax%s",str_PbPb.Data()));
	TH1D *h_PromptDs_BtoDs_ScaleBR_B0toDsmin=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_ScaleBR_B0toDsmin%s",str_PbPb.Data()));

	TH1D *h_PromptDs_BtoDs_ScaleBR_BptoDsmax=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_ScaleBR_BptoDsmax%s",str_PbPb.Data()));
	TH1D *h_PromptDs_BtoDs_ScaleBR_BptoDsmin=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_ScaleBR_BptoDsmin%s",str_PbPb.Data()));

	TH1D *h_PromptDs_BtoDs_ScaleBR_BstoDsmax=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_ScaleBR_BstoDsmax%s",str_PbPb.Data()));
	TH1D *h_PromptDs_BtoDs_ScaleBR_BstoDsmin=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_ScaleBR_BstoDsmin%s",str_PbPb.Data()));




	TH1D *h_PromptDs_BtoDs_ScaleFrZ=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_ScaleFrZ%s",str_PbPb.Data()));
	TH1D *h_PromptDs_BtoDs_ScaleFrp=(TH1D*)fin->Get(Form("h_PromptDs_BtoDs_ScaleFrp%s",str_PbPb.Data()));


	double CSdNdpt=0;
	double CSdNdptRel=0;
	double CSdNdptRel_ForRaa=0;
	double CSdNdptRel_ForRaaCancel=0;

	double CSdNdpt_BtoDs_FONLL=0;
	double CSdNdpt_BtoDs_FONLL_SysRel=0;

	double CSdNdpt_BtoDs_Raaup=0;
	double CSdNdpt_BtoDs_Raadown=0;
	double CSdNdpt_BtoDs_Raa_SysRel=0;

	double CSdNdpt_BtoDs_FONLLShape=0;
	double CSdNdpt_BtoDs_FONLLShape_SysRel=0;

  double CSdNdpt_BtoDs_Errup=0;
	double CSdNdpt_BtoDs_Errdown=0;
	double CSdNdpt_BtoDs_Err_SysRel=0;

  double CSdNdpt_BtoDs_ScaleBRup=0;
	double CSdNdpt_BtoDs_ScaleBRdown=0;
	double CSdNdpt_BtoDs_ScaleBR_SysRel=0;

  double CSdNdpt_BtoDs_ScaleBR_B0toDmax=0;
	double CSdNdpt_BtoDs_ScaleBR_B0toDmin=0;
	double CSdNdpt_BtoDs_ScaleBR_B0toD_SysRel=0;

  double CSdNdpt_BtoDs_ScaleBR_BptoDmax=0;
	double CSdNdpt_BtoDs_ScaleBR_BptoDmin=0;
	double CSdNdpt_BtoDs_ScaleBR_BptoD_SysRel=0;

  double CSdNdpt_BtoDs_ScaleBR_BstoDmax=0;
	double CSdNdpt_BtoDs_ScaleBR_BstoDmin=0;
	double CSdNdpt_BtoDs_ScaleBR_BstoD_SysRel=0;

  double CSdNdpt_BtoDs_ScaleBR_B0toDsmax=0;
	double CSdNdpt_BtoDs_ScaleBR_B0toDsmin=0;
	double CSdNdpt_BtoDs_ScaleBR_B0toDs_SysRel=0;

  double CSdNdpt_BtoDs_ScaleBR_BptoDsmax=0;
	double CSdNdpt_BtoDs_ScaleBR_BptoDsmin=0;
	double CSdNdpt_BtoDs_ScaleBR_BptoDs_SysRel=0;

  double CSdNdpt_BtoDs_ScaleBR_BstoDsmax=0;
	double CSdNdpt_BtoDs_ScaleBR_BstoDsmin=0;
	double CSdNdpt_BtoDs_ScaleBR_BstoDs_SysRel=0;


	double CSdNdpt_BtoDs_ScaleBRSum_SysRel=0;


  double CSdNdpt_BtoDs_ScaleFrZ=0;
	double CSdNdpt_BtoDs_ScaleFrp=0;
	double CSdNdpt_BtoDs_ScaleFr_SysRel=0;

	
	TH1D *h_PromptDs_BtoDs_SysRel_FONLL=new TH1D("h_PromptDs_BtoDs_SysRel_FONLL","h_PromptDs_BtoDs_SysRel_FONLL",nbin_pt,bins_pt);  h_PromptDs_BtoDs_SysRel_FONLL->Sumw2();

	TH1D *h_PromptDs_BtoDs_SysRel_Raa=new TH1D("h_PromptDs_BtoDs_SysRel_Raa","h_PromptDs_BtoDs_SysRel_Raa",nbin_pt,bins_pt);  h_PromptDs_BtoDs_SysRel_Raa->Sumw2();

	TH1D *h_PromptDs_BtoDs_SysRel_Err=new TH1D("h_PromptDs_BtoDs_SysRel_Err","h_PromptDs_BtoDs_SysRel_Err",nbin_pt,bins_pt);  h_PromptDs_BtoDs_SysRel_Err->Sumw2();

	TH1D *h_PromptDs_BtoDs_SysRel_FONLLShape=new TH1D("h_PromptDs_BtoDs_SysRel_FONLLShape","h_PromptDs_BtoDs_SysRel_FONLLShape",nbin_pt,bins_pt);  h_PromptDs_BtoDs_SysRel_FONLLShape->Sumw2();

	TH1D *h_PromptDs_BtoDs_SysRel_ScaleBR=new TH1D("h_PromptDs_BtoDs_SysRel_ScaleBR","h_PromptDs_BtoDs_SysRel_ScaleBR",nbin_pt,bins_pt);  h_PromptDs_BtoDs_SysRel_ScaleBR->Sumw2();

	TH1D *h_PromptDs_BtoDs_SysRel_ScaleFr=new TH1D("h_PromptDs_BtoDs_SysRel_ScaleFr","h_PromptDs_BtoDs_SysRel_ScaleFr",nbin_pt,bins_pt);  h_PromptDs_BtoDs_SysRel_ScaleFr->Sumw2();



	for(int i=startbin; i<nbin_pt; i++){
		// cout<<"\n\n # pt : "<<bins_pt[i]<<" - "<<bins_pt[i+1]<<endl;
		CSdNdpt=h_PromptDs->GetBinContent(i+1);

		CSdNdpt_BtoDs_FONLL=h_PromptDs_BtoDs_FONLL->GetBinContent(i+1);
		CSdNdpt_BtoDs_FONLL_SysRel=abs(CSdNdpt_BtoDs_FONLL-CSdNdpt)/CSdNdpt;
		h_PromptDs_BtoDs_FONLL_SysRel->SetBinContent(i+1,CSdNdpt_BtoDs_FONLL_SysRel);

		if(isPbPb && use_BtoD_Raa){
		CSdNdpt_BtoDs_Raaup=h_PromptDs_BtoDs_Raaup->GetBinContent(i+1);
		CSdNdpt_BtoDs_Raadown=h_PromptDs_BtoDs_Raadown->GetBinContent(i+1);
		CSdNdpt_BtoDs_Raa_SysRel= abs(CSdNdpt_BtoDs_Raaup-CSdNdpt) > abs(CSdNdpt_BtoDs_Raadown - CSdNdpt )? abs(CSdNdpt_BtoDs_Raaup-CSdNdpt)/CSdNdpt : abs(CSdNdpt_BtoDs_Raadown - CSdNdpt )/CSdNdpt	;
		}else{
		CSdNdpt_BtoDs_Raa_SysRel=0;
		}


		CSdNdpt_BtoDs_FONLLShape=h_PromptDs_BtoDs_FONLLShape->GetBinContent(i+1);
		CSdNdpt_BtoDs_FONLLShape_SysRel= abs(CSdNdpt_BtoDs_FONLLShape-CSdNdpt)/CSdNdpt ;

		CSdNdpt_BtoDs_Errup=h_PromptDs_BtoDs_Errup->GetBinContent(i+1);
		CSdNdpt_BtoDs_Errdown=h_PromptDs_BtoDs_Errdown->GetBinContent(i+1);
		CSdNdpt_BtoDs_Err_SysRel= abs(CSdNdpt_BtoDs_Errup-CSdNdpt) > abs(CSdNdpt_BtoDs_Errdown - CSdNdpt )? abs(CSdNdpt_BtoDs_Errup-CSdNdpt)/CSdNdpt : abs(CSdNdpt_BtoDs_Errdown - CSdNdpt )/CSdNdpt	;

		CSdNdpt_BtoDs_ScaleBRup=h_PromptDs_BtoDs_ScaleBRup->GetBinContent(i+1);
		CSdNdpt_BtoDs_ScaleBRdown=h_PromptDs_BtoDs_ScaleBRdown->GetBinContent(i+1);
		CSdNdpt_BtoDs_ScaleBR_SysRel= abs(CSdNdpt_BtoDs_ScaleBRup-CSdNdpt) > abs(CSdNdpt_BtoDs_ScaleBRdown - CSdNdpt )? abs(CSdNdpt_BtoDs_ScaleBRup-CSdNdpt)/CSdNdpt : abs(CSdNdpt_BtoDs_ScaleBRdown - CSdNdpt )/CSdNdpt	;

		// scale BR seperately

		CSdNdpt_BtoDs_ScaleBR_B0toDmax=h_PromptDs_BtoDs_ScaleBR_B0toDmax->GetBinContent(i+1);
		CSdNdpt_BtoDs_ScaleBR_B0toDmin=h_PromptDs_BtoDs_ScaleBR_B0toDmin->GetBinContent(i+1);
		CSdNdpt_BtoDs_ScaleBR_B0toD_SysRel= abs(CSdNdpt_BtoDs_ScaleBR_B0toDmax-CSdNdpt) > abs(CSdNdpt_BtoDs_ScaleBR_B0toDmin - CSdNdpt )? abs(CSdNdpt_BtoDs_ScaleBR_B0toDmax-CSdNdpt)/CSdNdpt : abs(CSdNdpt_BtoDs_ScaleBR_B0toDmin - CSdNdpt )/CSdNdpt	;

		CSdNdpt_BtoDs_ScaleBR_BptoDmax=h_PromptDs_BtoDs_ScaleBR_BptoDmax->GetBinContent(i+1);
		CSdNdpt_BtoDs_ScaleBR_BptoDmin=h_PromptDs_BtoDs_ScaleBR_BptoDmin->GetBinContent(i+1);
		CSdNdpt_BtoDs_ScaleBR_BptoD_SysRel= abs(CSdNdpt_BtoDs_ScaleBR_BptoDmax-CSdNdpt) > abs(CSdNdpt_BtoDs_ScaleBR_BptoDmin - CSdNdpt )? abs(CSdNdpt_BtoDs_ScaleBR_BptoDmax-CSdNdpt)/CSdNdpt : abs(CSdNdpt_BtoDs_ScaleBR_BptoDmin - CSdNdpt )/CSdNdpt	;

		CSdNdpt_BtoDs_ScaleBR_BstoDmax=h_PromptDs_BtoDs_ScaleBR_BstoDmax->GetBinContent(i+1);
		CSdNdpt_BtoDs_ScaleBR_BstoDmin=h_PromptDs_BtoDs_ScaleBR_BstoDmin->GetBinContent(i+1);
		CSdNdpt_BtoDs_ScaleBR_BstoD_SysRel= abs(CSdNdpt_BtoDs_ScaleBR_BstoDmax-CSdNdpt) > abs(CSdNdpt_BtoDs_ScaleBR_BstoDmin - CSdNdpt )? abs(CSdNdpt_BtoDs_ScaleBR_BstoDmax-CSdNdpt)/CSdNdpt : abs(CSdNdpt_BtoDs_ScaleBR_BstoDmin - CSdNdpt )/CSdNdpt	;

		CSdNdpt_BtoDs_ScaleBR_B0toDsmax=h_PromptDs_BtoDs_ScaleBR_B0toDsmax->GetBinContent(i+1);
		CSdNdpt_BtoDs_ScaleBR_B0toDsmin=h_PromptDs_BtoDs_ScaleBR_B0toDsmin->GetBinContent(i+1);
		CSdNdpt_BtoDs_ScaleBR_B0toDs_SysRel= abs(CSdNdpt_BtoDs_ScaleBR_B0toDsmax-CSdNdpt) > abs(CSdNdpt_BtoDs_ScaleBR_B0toDsmin - CSdNdpt )? abs(CSdNdpt_BtoDs_ScaleBR_B0toDsmax-CSdNdpt)/CSdNdpt : abs(CSdNdpt_BtoDs_ScaleBR_B0toDsmin - CSdNdpt )/CSdNdpt	;

		CSdNdpt_BtoDs_ScaleBR_BptoDsmax=h_PromptDs_BtoDs_ScaleBR_BptoDsmax->GetBinContent(i+1);
		CSdNdpt_BtoDs_ScaleBR_BptoDsmin=h_PromptDs_BtoDs_ScaleBR_BptoDsmin->GetBinContent(i+1);
		CSdNdpt_BtoDs_ScaleBR_BptoDs_SysRel= abs(CSdNdpt_BtoDs_ScaleBR_BptoDsmax-CSdNdpt) > abs(CSdNdpt_BtoDs_ScaleBR_BptoDsmin - CSdNdpt )? abs(CSdNdpt_BtoDs_ScaleBR_BptoDsmax-CSdNdpt)/CSdNdpt : abs(CSdNdpt_BtoDs_ScaleBR_BptoDsmin - CSdNdpt )/CSdNdpt	;

		CSdNdpt_BtoDs_ScaleBR_BstoDsmax=h_PromptDs_BtoDs_ScaleBR_BstoDsmax->GetBinContent(i+1);
		CSdNdpt_BtoDs_ScaleBR_BstoDsmin=h_PromptDs_BtoDs_ScaleBR_BstoDsmin->GetBinContent(i+1);
		CSdNdpt_BtoDs_ScaleBR_BstoDs_SysRel= abs(CSdNdpt_BtoDs_ScaleBR_BstoDsmax-CSdNdpt) > abs(CSdNdpt_BtoDs_ScaleBR_BstoDsmin - CSdNdpt )? abs(CSdNdpt_BtoDs_ScaleBR_BstoDsmax-CSdNdpt)/CSdNdpt : abs(CSdNdpt_BtoDs_ScaleBR_BstoDsmin - CSdNdpt )/CSdNdpt	;



	// end BR 





		CSdNdpt_BtoDs_ScaleFrZ=h_PromptDs_BtoDs_ScaleFrZ->GetBinContent(i+1);
		CSdNdpt_BtoDs_ScaleFrp=h_PromptDs_BtoDs_ScaleFrp->GetBinContent(i+1);
		CSdNdpt_BtoDs_ScaleFr_SysRel= abs(CSdNdpt_BtoDs_ScaleFrZ-CSdNdpt) > abs(CSdNdpt_BtoDs_ScaleFrp - CSdNdpt )? abs(CSdNdpt_BtoDs_ScaleFrZ-CSdNdpt)/CSdNdpt : abs(CSdNdpt_BtoDs_ScaleFrp - CSdNdpt )/CSdNdpt	;

		// CSdNdpt_Max=abs(CSdNdpt - CSdNdpt_Pythia) > abs(CSdNdpt - CSdNdpt_FONLL) ? 	abs(CSdNdpt - CSdNdpt_Pythia) : abs(CSdNdpt - CSdNdpt_FONLL);
		// CSdNdptRel=sqrt( CSdNdpt_BtoDs_Raa_SysRel*CSdNdpt_BtoDs_Raa_SysRel + CSdNdpt_BtoDs_FONLLShape_SysRel*CSdNdpt_BtoDs_FONLLShape_SysRel + CSdNdpt_BtoDs_Err_SysRel*CSdNdpt_BtoDs_Err_SysRel + CSdNdpt_BtoDs_ScaleBR_SysRel*CSdNdpt_BtoDs_ScaleBR_SysRel+ CSdNdpt_BtoDs_ScaleFr_SysRel*CSdNdpt_BtoDs_ScaleFr_SysRel  );
		CSdNdptRel=sqrt( CSdNdpt_BtoDs_Raa_SysRel*CSdNdpt_BtoDs_Raa_SysRel + CSdNdpt_BtoDs_FONLLShape_SysRel*CSdNdpt_BtoDs_FONLLShape_SysRel + CSdNdpt_BtoDs_Err_SysRel*CSdNdpt_BtoDs_Err_SysRel + CSdNdpt_BtoDs_ScaleFr_SysRel*CSdNdpt_BtoDs_ScaleFr_SysRel  +pow(CSdNdpt_BtoDs_ScaleBR_B0toD_SysRel,2) +pow(CSdNdpt_BtoDs_ScaleBR_BptoD_SysRel,2) +pow(CSdNdpt_BtoDs_ScaleBR_BstoD_SysRel,2) +pow(CSdNdpt_BtoDs_ScaleBR_B0toDs_SysRel,2) +pow(CSdNdpt_BtoDs_ScaleBR_BptoDs_SysRel,2) +pow(CSdNdpt_BtoDs_ScaleBR_BstoDs_SysRel,2));


		CSdNdpt_BtoDs_ScaleBRSum_SysRel=sqrt(pow(CSdNdpt_BtoDs_ScaleBR_B0toD_SysRel,2) +pow(CSdNdpt_BtoDs_ScaleBR_BptoD_SysRel,2) +pow(CSdNdpt_BtoDs_ScaleBR_BstoD_SysRel,2) +pow(CSdNdpt_BtoDs_ScaleBR_B0toDs_SysRel,2) +pow(CSdNdpt_BtoDs_ScaleBR_BptoDs_SysRel,2) +pow(CSdNdpt_BtoDs_ScaleBR_BstoDs_SysRel,2));


		h_PromptDs_BtoDs_SysRel->SetBinContent(i+1,CSdNdptRel);
	
		// CSdNdptRel_ForRaa=sqrt(CSdNdpt_BtoDs_Raa_SysRel*CSdNdpt_BtoDs_Raa_SysRel + CSdNdpt_BtoDs_FONLLShape_SysRel*CSdNdpt_BtoDs_FONLLShape_SysRel + CSdNdpt_BtoDs_Err_SysRel*CSdNdpt_BtoDs_Err_SysRel);
		CSdNdptRel_ForRaa=sqrt( CSdNdpt_BtoDs_Raa_SysRel*CSdNdpt_BtoDs_Raa_SysRel + CSdNdpt_BtoDs_FONLLShape_SysRel*CSdNdpt_BtoDs_FONLLShape_SysRel + CSdNdpt_BtoDs_Err_SysRel*CSdNdpt_BtoDs_Err_SysRel);
		CSdNdptRel_ForRaaCancel=sqrt( CSdNdpt_BtoDs_ScaleFr_SysRel*CSdNdpt_BtoDs_ScaleFr_SysRel  +pow(CSdNdpt_BtoDs_ScaleBR_B0toD_SysRel,2) +pow(CSdNdpt_BtoDs_ScaleBR_BptoD_SysRel,2) +pow(CSdNdpt_BtoDs_ScaleBR_BstoD_SysRel,2) +pow(CSdNdpt_BtoDs_ScaleBR_B0toDs_SysRel,2) +pow(CSdNdpt_BtoDs_ScaleBR_BptoDs_SysRel,2) +pow(CSdNdpt_BtoDs_ScaleBR_BstoDs_SysRel,2) );
		h_PromptDs_BtoDs_SysRel_ForRaa->SetBinContent(i+1,CSdNdptRel_ForRaa);
		h_PromptDs_BtoDs_SysRel_ForRaaCancel->SetBinContent(i+1,CSdNdptRel_ForRaaCancel);


		h_PromptDs_BtoDs_SysRel_FONLL->SetBinContent(i+1,CSdNdpt_BtoDs_FONLL_SysRel);
		h_PromptDs_BtoDs_SysRel_Raa->SetBinContent(i+1,CSdNdpt_BtoDs_Raa_SysRel);
		h_PromptDs_BtoDs_SysRel_FONLLShape->SetBinContent(i+1,CSdNdpt_BtoDs_FONLLShape_SysRel);
		h_PromptDs_BtoDs_SysRel_Err->SetBinContent(i+1,CSdNdpt_BtoDs_Err_SysRel);
		h_PromptDs_BtoDs_SysRel_ScaleBR->SetBinContent(i+1,CSdNdpt_BtoDs_ScaleBR_SysRel);
		h_PromptDs_BtoDs_SysRel_ScaleFr->SetBinContent(i+1,CSdNdpt_BtoDs_ScaleFr_SysRel);

		cout<<setprecision(1)<<std::fixed;

		cout<<"\n\n -------  pt : "<<bins_pt[i]<<" - "<<bins_pt[i+1]<<endl;
		cout<<setw(40)<<"CSdNdptRel All = "<<CSdNdptRel*100<<"% , FONLL rel = "<<CSdNdpt_BtoDs_FONLL_SysRel*100<<"% , CSdNdpt = "<<CSdNdpt<<endl;
		cout<<setw(40)<<"CSdNdpt_BtoDs_ScaleFr_SysRel = "<<CSdNdpt_BtoDs_ScaleFr_SysRel*100<<"%" <<endl;
		cout<<setw(40)<<"CSdNdpt_BtoDs_ScaleBRSum_SysRel = "<<CSdNdpt_BtoDs_ScaleBRSum_SysRel*100<<"%" <<endl;
		cout<<setw(40)<<"CSdNdpt_BtoDs_Err_SysRel = "<<CSdNdpt_BtoDs_Err_SysRel*100<<"%" <<endl;
		cout<<setw(40)<<"CSdNdpt_BtoDs_FONLLShape_SysRel = "<<CSdNdpt_BtoDs_FONLLShape_SysRel*100<<"%" <<endl;

		int print_extra=0;
		if(print_extra){
		cout<<"----------\n extra  "<<endl;
		cout<<setw(40)<<"CSdNdptRel_ForRaa = "<<CSdNdptRel_ForRaa*100<<"% "<<endl;
		cout<<setw(40)<<"CSdNdpt_BtoDs_Raa_SysRel = "<<CSdNdpt_BtoDs_Raa_SysRel*100 <<"%" <<endl;
		cout<<setw(40)<<"CSdNdpt_BtoDs_Err_SysRel = "<<CSdNdpt_BtoDs_Err_SysRel*100<<"%" <<endl;
		// cousetw(40)<<t<<"CSdNdpt_BtoDs_ScaleBR_SysRel = "<<CSdNdpt_BtoDs_ScaleBR_SysRel*100<<"%" <<endl;

		cout<<setw(40)<<"CSdNdpt_BtoDs_ScaleBR_B0toD_SysRel = "<<CSdNdpt_BtoDs_ScaleBR_B0toD_SysRel*100 <<"%" <<endl;
		cout<<setw(40)<<"CSdNdpt_BtoDs_ScaleBR_BptoD_SysRel = "<<CSdNdpt_BtoDs_ScaleBR_BptoD_SysRel*100 <<"%" <<endl;
		cout<<setw(40)<<"CSdNdpt_BtoDs_ScaleBR_BstoD_SysRel = "<<CSdNdpt_BtoDs_ScaleBR_BstoD_SysRel*100<<"%" <<endl;
		cout<<setw(40)<<"CSdNdpt_BtoDs_ScaleBR_B0toDs_SysRel = "<<CSdNdpt_BtoDs_ScaleBR_B0toDs_SysRel*100 <<"%" <<endl;
		cout<<setw(40)<<"CSdNdpt_BtoDs_ScaleBR_BptoDs_SysRel = "<<CSdNdpt_BtoDs_ScaleBR_BptoDs_SysRel*100 <<"%" <<endl;
		cout<<setw(40)<<"CSdNdpt_BtoDs_ScaleBR_BstoDs_SysRel = "<<CSdNdpt_BtoDs_ScaleBR_BstoDs_SysRel*100 <<"%" <<endl;
		}


	}



	// for(int i=startbin; i<nbin_pt; i++){

	// }




	fout->cd();


	h_PromptDs_BtoDs_SysRel->Write("",TObject::kOverwrite);
	h_PromptDs_BtoDs_SysRel_ForRaa->Write("",TObject::kOverwrite);
	h_PromptDs_BtoDs_SysRel_ForRaaCancel->Write("",TObject::kOverwrite);


	h_PromptDs_BtoDs_SysRel_FONLL->Write("",TObject::kOverwrite);
	h_PromptDs_BtoDs_SysRel_Raa->Write("",TObject::kOverwrite);
	h_PromptDs_BtoDs_SysRel_FONLLShape->Write("",TObject::kOverwrite);
	h_PromptDs_BtoDs_SysRel_Err->Write("",TObject::kOverwrite);
	h_PromptDs_BtoDs_SysRel_ScaleBR->Write("",TObject::kOverwrite);
	h_PromptDs_BtoDs_SysRel_ScaleFr->Write("",TObject::kOverwrite);

	return 0;
}
