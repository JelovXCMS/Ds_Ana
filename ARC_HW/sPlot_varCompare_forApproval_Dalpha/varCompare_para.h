#ifndef _VARCOMPARE_PARAMETERS_H_
#define _VARCOMPARE_PARAMETERS_H_
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"


#include <TLatex.h>
#include <TCut.h>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <iostream>
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



using namespace std;


  Float_t Dalpha_cut_pp=0.2;
  Float_t Dchi2cl_cut_pp=0.05;
  Float_t Ddls_cut_pp=2.5;

  Float_t Ddls_maxcut_pp=14;
  Float_t Dalpha_maxcut_pp=0.2;
  Float_t Dchi2cl_mincut_pp=0.02;

  Float_t Dalpha_cut_PbPb3=0.2;
  Float_t Dchi2cl_cut_PbPb3=0.2;
  Float_t Ddls_cut_PbPb3=3.5;

  Float_t Ddls_maxcut_PbPb3=14;
  Float_t Ddls_mincut_PbPb3=3.5;
  Float_t Dalpha_maxcut_PbPb3=0.2;
  Float_t Dchi2cl_mincut_PbPb3=0.05;


	double Dpt_Low_pp=2;
	double Dpt_Hight_pp=40;
	double Dpt_Low_PbPb=6;
	double Dpt_Hight_PbPb=40;

  double bins_Ddls_pp[]={1.5,2,2.5,3,3.5,4,5,6,8,10,14};
  const int nbin_Ddls_pp= sizeof(bins_Ddls_pp)/sizeof(bins_Ddls_pp[1])-1;
	
  double bins_Dalpha_pp[]={0,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20};
  const int nbin_Dalpha_pp= sizeof(bins_Dalpha_pp)/sizeof(bins_Dalpha_pp[1])-1;

  double bins_Dchi2cl_pp[]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
  const int nbin_Dchi2cl_pp= sizeof(bins_Dchi2cl_pp)/sizeof(bins_Dchi2cl_pp[1])-1;


  double bins_Ddls_PbPb[]={3.5,4,4.5,5,6,8,10,14};
  const int nbin_Ddls_PbPb= sizeof(bins_Ddls_PbPb)/sizeof(bins_Ddls_PbPb[1])-1;
	
  double bins_Dalpha_PbPb[]={0,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20};
  const int nbin_Dalpha_PbPb= sizeof(bins_Dalpha_PbPb)/sizeof(bins_Dalpha_PbPb[1])-1;

  double bins_Dchi2cl_PbPb[]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
  const int nbin_Dchi2cl_PbPb= sizeof(bins_Dchi2cl_PbPb)/sizeof(bins_Dchi2cl_PbPb[1])-1;
	
	double bins_DdxyzErr_pp[]={0,0.01,0.015,0.0175,0.02,0.0225,0.025,0.0275,0.03,0.0325,0.035,0.0375,0.04,0.045,0.05,0.06,0.08,0.1,1};
	const int nbin_DdxyzErr_pp=sizeof(bins_DdxyzErr_pp)/sizeof(bins_DdxyzErr_pp[1])-1;



#endif

