#ifndef _ANA_BFEEDDOWN_PARAMETERS_H_
#define _ANA_BFEEDDOWN_PARAMETERS_H_

#include <TLatex.h>
#include <TCut.h>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <iostream>

// TMVA cut need to be read by initParameter function , not default cut

using namespace std;

bool verbose=true;

// for plot
  int c_W = 800;
  int c_H = 600;

	int c_wtopx=50;
	int c_wtopy=50;

  float c_Tmg = 0.08;
  float c_Bmg = 0.12; 
  float c_Lmg = 0.15;
  float c_Rmg = 0.04;


const float PI = 3.14159265359;

	double NevtPbPb3=292362708; // Trigger & Event selection , without PvZ cut, others use 2.94, part of my file failed in Dfiner but Jian refuse to resubmit
	double TAA0to100=5.607*1e-9;  // 5.607/ mb , 5.58 in Hao's Anafrom https://twiki.cern.ch/twiki/pub/CMS/HiCentrality2016/AN-15-080_temp_20161206.pdf
	double LumiSum=0.0382;  //(pb)
  double BRphi=0.0227;
  double BRf0=0.0115;

  double bins_pt_Gen[]={0,1.8,3.8,5.7,9.5,19,100};
  const int nbin_pt_Gen=sizeof(bins_pt_Gen)/sizeof(bins_pt_Gen[0]) -1;

  double bins_pt_GenD0[]={0,2,4,10,20,40,60,100};
  const int nbin_pt_GenD0=sizeof(bins_pt_GenD0)/sizeof(bins_pt_GenD0[0]) -1;

  // double bins_pt_Norm[]={0,2,3,4,5,6,8,10,20,40,100};
  double bins_pt_Norm[]={0,2,3,4,5,6,8,10,20,40,60,100}; // for consistency with GenD0
  const int nbin_pt_Norm=sizeof(bins_pt_Norm)/sizeof(bins_pt_Norm[0])-1;


// for dca bin
	// const int nbin_dca=13;
	const int nbin_dca=13;
	float bins_dca[nbin_dca+1]; // initialize in function
	float firstBinDcaWidth=0.001;
	// float binDcaWidthRatio=1.27;
	float binDcaWidthRatio=1.27;

//const int nPtBins = 9;
//Double_t ptBins[nPtBins+1] = {1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0, 25.0, 40.0};


//const int nPtBins = 11;
//Double_t ptBins[nPtBins+1] = {1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0, 40.0, 60.0,100.0};

//const int nPtBins = 12;
//Double_t ptBins[nPtBins+1] = { 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0 ,40.0, 60.0,80.0,100.0};


// original bin from hao

// general parameter part
  double DtktkResmassCutMean=1.01951;
  double DtktkResmassCutWidth=0.009;

	const int nbin_PhiMassScan=7;

	double PhiMassScan_Min=0.008;
	double PhiMassScan_Max=0.011;

	double PhiMassScan_bins[nbin_PhiMassScan]={0.008,0.0085,0.0090,0.0095,0.0100,0.0105,0.0110};
	double bins_PhiMassScan[nbin_PhiMassScan+1]={0.00775,0.00825,0.00875,0.00925,0.00975,0.01025,0.01075,0.01125};


  double bins_pt_pp[]={2,3,4,5,6,8,10,20,40}; // for test
  const int nbin_pt_pp=sizeof(bins_pt_pp)/sizeof(bins_pt_pp[0])-1;

	double DataFitRangeLow=1.91;
	double DataFitRangeHigh=2.11;

	double DsMass=1.96828;  //fitted 1.97
  double DsMassRange=0.03;

  double PhiMass=1.01951; // not used
  double PhiMassRange=0.009;

  double TrketaAcc=1.5;

// pp part
  double TrkptAcc_pp=0.7;
//  double TrkptAcc=TrkptAcc_pp;

  double Dchi2clMin_bins_pp[nbin_pt_pp]={0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05};
  double DdlsMin_bins_pp[nbin_pt_pp]={4,4,4,4,4,4,4,4};
  double DalphaMax_bins_pp[nbin_pt_pp]={0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12};

	double Dtrk1PtMin_bins_pp[nbin_pt_pp]={0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7};
  double Dtrk2PtMin_bins_pp[nbin_pt_pp]={0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7};
  double Dtrk3PtMin_bins_pp[nbin_pt_pp]={0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7};

  double Dchi2clMinLoose_pp=0.02;  
  double DalphaMaxLoose_pp=0.12; // fix Dalpha cut in TMVA 
  double DdlsMinLoose_pp=1.5;
	double DtrkPtMinLoose_pp=0.7;

	TString mycut_PP=Form("Dchi2cl>%f && Dalpha<%f && Ddls>%f&& Dtrk1Pt>%f && Dtrk2Pt>%f && Dtrk3Pt>%f && TMath::Abs(DtktkResmass-1.01951) <%f",Dchi2clMinLoose_pp, DalphaMaxLoose_pp, DdlsMinLoose_pp, DtrkPtMinLoose_pp, DtrkPtMinLoose_pp,DtrkPtMinLoose_pp, PhiMassRange);
	TString mycuts_PP = Form("(%s)&&DsGen==23333 && DgencollisionId==0 ",mycut_PP.Data());
	TString mycutb_PP = Form("(%s)&&abs(Dmass-1.97)>0.04&&abs(Dmass-1.97)<0.6",mycut_PP.Data());


	const int nbin_Dchi2clMinScan=5;
	const int nbin_DdlsMinScan=5;
	const int nbin_DalphaMaxScan=5;

	double Dchi2clMinScan_Min_pp=0.02;
	double Dchi2clMinScan_Max_pp=0.42;

	double DdlsMinScan_Min_pp=1.5;
	double DdlsMinScan_Max_pp=5.5;

	double DalphaMaxScan_Min_pp=0.08;
	double DalphaMaxScan_Max_pp=0.2;

	double Dchi2clMinScan_bins_pp[nbin_Dchi2clMinScan];
	double DdlsMinScan_bins_pp[nbin_DdlsMinScan];
	double DalphaMaxScan_bins_pp[nbin_DalphaMaxScan];


	double bins_Dchi2clMinScan_pp[nbin_Dchi2clMinScan+1];
	double bins_DdlsMinScan_pp[nbin_DdlsMinScan+1];
	double bins_DalphaMaxScan_pp[nbin_DalphaMaxScan+1];

/* 
 	TString ppMCFile[nbin_pt_pp]={                                                                                                        "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_phikkpi_pt4.root"
 , "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_phikkpi_pt4.root"
 , "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_phikkpi_pt4.root"
 , "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_phikkpi_pt4.root"};
*/ 

	TString ppDataFile="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_Data/DsMinTree_pp_Data_MBAll.root";

	TString MCFile_Prompt_phikkpi_pp_Merge="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_MergeWeight/root_output/DsMinTree_MC_GenSampleMerge_pp_Prompt_phi.root";
	TString MCFile_Prompt_f0980kkpi_pp_Merge="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_MergeWeight/root_output/DsMinTree_MC_GenSampleMerge_pp_Prompt_f0980.root";
	TString MCFile_NonPrompt_phikkpi_pp_Merge="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_MergeWeight/root_output/DsMinTree_MC_GenSampleMerge_pp_NonPrompt_phi.root";
	TString MCFile_NonPrompt_f0980kkpi_pp_Merge="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_MergeWeight/root_output/DsMinTree_MC_GenSampleMerge_pp_NonPrompt_f0980.root";

/*
	TString MCFile_Prompt_phikkpi_pp[nbin_pt_pp]={
  "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_Prompt_phikkpi_pt1p8.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_Prompt_phikkpi_pt1p8.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_Prompt_phikkpi_pt3p8.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_Prompt_phikkpi_pt5p7.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_Prompt_phikkpi_pt5p7.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_Prompt_phikkpi_pt9p5.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_Prompt_phikkpi_pt19.root" 
	};

	TString MCFile_NonPrompt_phikkpi_pp[nbin_pt_pp]={
	"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_NonPrompt_phikkpi_pt1p8.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_NonPrompt_phikkpi_pt1p8.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_NonPrompt_phikkpi_pt3p8.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_NonPrompt_phikkpi_pt5p7.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_NonPrompt_phikkpi_pt5p7.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_NonPrompt_phikkpi_pt9p5.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_NonPrompt_phikkpi_pt19.root" 
	};

	TString MCFile_Prompt_f0980kkpi_pp[nbin_pt_pp]={
	"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_Prompt_f0980kkpi_pt1p8.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_Prompt_f0980kkpi_pt1p8.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_Prompt_f0980kkpi_pt3p8.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_Prompt_f0980kkpi_pt5p7.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_Prompt_f0980kkpi_pt5p7.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_Prompt_f0980kkpi_pt9p5.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_Prompt_f0980kkpi_pt19.root" 
	};

	TString MCFile_NonPrompt_f0980kkpi_pp[nbin_pt_pp]={
	"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_NonPrompt_f0980kkpi_pt1p8.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_NonPrompt_f0980kkpi_pt1p8.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_NonPrompt_f0980kkpi_pt3p8.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_NonPrompt_f0980kkpi_pt5p7.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_NonPrompt_f0980kkpi_pt5p7.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_NonPrompt_f0980kkpi_pt9p5.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_NonPrompt_f0980kkpi_pt19.root" 
	};
*/


// PbPb part 1 // for cent 30-100


  double bins_pt_PbPb1[]={3,4,6,8,10,20,40}; // for test
  const int nbin_pt_PbPb1=sizeof(bins_pt_PbPb1)/sizeof(bins_pt_PbPb1[0])-1;

	TString PbPb1MCFile[nbin_pt_PbPb1]={
  "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/PbPb_MC/DsMinTree_PbPb_Prompt_phikkpi_pt1p8.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/PbPb_MC/DsMinTree_PbPb_MC_Ds_phikkpi_pt4.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/PbPb_MC/DsMinTree_PbPb_MC_Ds_phikkpi_pt4.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/PbPb_MC/DsMinTree_PbPb_MC_Ds_phikkpi_pt4.root"
};

  TString PbPb1DataFile="";


	TString MCFile_Prompt_phikkpi_PbPb_Merge="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_MergeWeight/root_output/DsMinTree_MC_GenSampleMerge_PbPb_Prompt_phi.root";
	TString MCFile_Prompt_f0980kkpi_PbPb_Merge="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_MergeWeight/root_output/DsMinTree_MC_GenSampleMerge_PbPb_Prompt_f0980.root";
	TString MCFile_NonPrompt_phikkpi_PbPb_Merge="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_MergeWeight/root_output/DsMinTree_MC_GenSampleMerge_PbPb_NonPrompt_phi.root";
	TString MCFile_NonPrompt_f0980kkpi_PbPb_Merge="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_MergeWeight/root_output/DsMinTree_MC_GenSampleMerge_PbPb_NonPrompt_f0980.root";



  double TrkptAcc_PbPb1=0.7;
  double Dchi2clMin_bins_PbPb1[nbin_pt_PbPb1]={0.05,0.05,0.05,0.05,0.05,0.05};
  double DalphaMax_bins_PbPb1[nbin_pt_PbPb1]={0.12,0.12,0.12,0.12,0.12,0.12};
  double DdlsMin_bins_PbPb1[nbin_pt_PbPb1]={3.5,3.5,3.5,3.5,3.5,3.5};

  double Dtrk1PtMin_bins_PbPb1[nbin_pt_PbPb1]={1,1,1,1,1,1};
  double Dtrk2PtMin_bins_PbPb1[nbin_pt_PbPb1]={1,1,1,1,1,1};
  double Dtrk3PtMin_bins_PbPb1[nbin_pt_PbPb1]={1,1,1,1,1,1};

  double Dchi2clMinLoose_PbPb1=0.05;
  double DalphaMaxLoose_PbPb1=0.2;
  double DdlsMinLoose_PbPb1=2.5;
	double DtrkPtMinLoose_PbPb1=1;

	int hiBinLow_PbPb1=60;
	int hiBinHigh_PbPb1=201;

	TCut cuthiBin_PbPb1="hiBin>60";


	double Dchi2clMinScan_Min_PbPb1=0.05;
	double Dchi2clMinScan_Max_PbPb1=0.45;

	double DdlsMinScan_Min_PbPb1=2.5;
	double DdlsMinScan_Max_PbPb1=6.5;

	double DalphaMaxScan_Min_PbPb1=0.08;
	double DalphaMaxScan_Max_PbPb1=0.2;

	double Dchi2clMinScan_bins_PbPb1[nbin_Dchi2clMinScan];
	double DdlsMinScan_bins_PbPb1[nbin_DdlsMinScan];
	double DalphaMaxScan_bins_PbPb1[nbin_DalphaMaxScan];

 



	// PbPb part2 // for cent 0-30


  double bins_pt_PbPb2[]={3,4,6,8,10}; // for test
  const int nbin_pt_PbPb2=sizeof(bins_pt_PbPb2)/sizeof(bins_pt_PbPb2[0])-1;

	TString PbPb2MCFile[nbin_pt_PbPb2]={
  "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/PbPb_MC/DsMinTree_PbPb_Prompt_phikkpi_pt1p8.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/PbPb_MC/DsMinTree_PbPb_MC_Ds_phikkpi_pt4.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/PbPb_MC/DsMinTree_PbPb_MC_Ds_phikkpi_pt4.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/PbPb_MC/DsMinTree_PbPb_MC_Ds_phikkpi_pt4.root"
};

  double TrkptAcc_PbPb2=0.7;
  double Dchi2clMin_bins_PbPb2[nbin_pt_PbPb2]={0.05,0.05,0.05,0.05};
  double DalphaMax_bins_PbPb2[nbin_pt_PbPb2]={0.12,0.12,0.12,0.12};
  double DdlsMin_bins_PbPb2[nbin_pt_PbPb2]={3.5,3.5,3.5,3.5};

  double Dtrk1PtMin_bins_PbPb2[nbin_pt_PbPb2]={1,1,1,1};
  double Dtrk2PtMin_bins_PbPb2[nbin_pt_PbPb2]={1,1,1,1};
  double Dtrk3PtMin_bins_PbPb2[nbin_pt_PbPb2]={1,1,1,1};

  double Dchi2clMinLoose_PbPb2=0.05;
  double DalphaMaxLoose_PbPb2=0.2; 
  double DdlsMinLoose_PbPb2=2.5;
  double DtrkPtMinLoose_PbPb2=1;

	TCut cuthiBin_PbPb2="hiBin<60";

	int hiBinLow_PbPb2=-1;
	int hiBinHigh_PbPb2=60;


	double Dchi2clMinScan_Min_PbPb2=0.05;
	double Dchi2clMinScan_Max_PbPb2=0.45;

	double DdlsMinScan_Min_PbPb2=2.5;
	double DdlsMinScan_Max_PbPb2=6.5;

	double DalphaMaxScan_Min_PbPb2=0.08;
	double DalphaMaxScan_Max_PbPb2=0.2;

	double Dchi2clMinScan_bins_PbPb2[nbin_Dchi2clMinScan];
	double DdlsMinScan_bins_PbPb2[nbin_DdlsMinScan];
	double DalphaMaxScan_bins_PbPb2[nbin_DalphaMaxScan];


	// PbPb part3 // for cent 0-100

  double bins_pt_PbPb3[]={4,5,6,8,10,20,40}; // for test
  const int nbin_pt_PbPb3=sizeof(bins_pt_PbPb3)/sizeof(bins_pt_PbPb3[0])-1;


  // TString PbPb3DataFile="DsMinTree_PbPb_Data_HIMB3_GJ_part2.root"; // only for test
  TString PbPb3DataFile_GJMB3part2="/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/GJ/DsMinTree_PbPb_Data_HIMB3_GJ_part2.root"; // only for test
  TString PbPb3DataFile_test="/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/GJ/Ds_PbPb_Data_HIMB234_GJ_forTMVAtest.root"; // only for test
/*
	TString PbPb3MCFile[nbin_pt_PbPb3]={
  "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/PbPb_MC/DsMinTree_PbPb_Prompt_phikkpi_pt1p8.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/PbPb_MC/DsMinTree_PbPb_MC_Ds_phikkpi_pt4.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/PbPb_MC/DsMinTree_PbPb_MC_Ds_phikkpi_pt4.root"
 ,"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/PbPb_MC/DsMinTree_PbPb_MC_Ds_phikkpi_pt4.root"
};
*/
  double TrkptAcc_PbPb3=0.7;
  double Dchi2clMin_bins_PbPb3[nbin_pt_PbPb3]={0.05,0.05,0.05,0.05,0.05,0.05};
  double DalphaMax_bins_PbPb3[nbin_pt_PbPb3]={0.12,0.12,0.12,0.12,0.12,0.12};
  double DdlsMin_bins_PbPb3[nbin_pt_PbPb3]={3.5,3.5,3.5,3.5,3.5,3.5};

  double Dtrk1PtMin_bins_PbPb3[nbin_pt_PbPb3]={1,1,1,1,1,1};
  double Dtrk2PtMin_bins_PbPb3[nbin_pt_PbPb3]={1,1,1,1,1,1};
  double Dtrk3PtMin_bins_PbPb3[nbin_pt_PbPb3]={1,1,1,1,1,1};

  double Dchi2clMinLoose_PbPb3=0.05;
  double DalphaMaxLoose_PbPb3=0.12; // just changed, not apply and generate new MC_eff & fitfile yet
  double DdlsMinLoose_PbPb3=2.5;
	double DtrkPtMinLoose_PbPb3=1;


	TString mycut_PbPb3=Form("Dchi2cl>%f && Dalpha<%f && Ddls>%f&& Dtrk1Pt>%f && Dtrk2Pt>%f && Dtrk3Pt>%f && TMath::Abs(DtktkResmass-1.01951) <%f",Dchi2clMinLoose_PbPb3, DalphaMaxLoose_PbPb3, DdlsMinLoose_PbPb3, DtrkPtMinLoose_PbPb3, DtrkPtMinLoose_PbPb3,DtrkPtMinLoose_PbPb3, PhiMassRange);
	TString mycuts_PbPb3 = Form("(%s)&&DsGen==23333 && DgencollisionId==0 ",mycut_PbPb3.Data());
	TString mycutb_PbPb3 = Form("(%s)&&abs(Dmass-1.97)>0.04&&abs(Dmass-1.97)<0.6",mycut_PbPb3.Data());


TCut cuthiBin_PbPb3="hiBin<200";

	int hiBinLow_PbPb3=0;
	int hiBinHigh_PbPb3=200;


	double Dchi2clMinScan_Min_PbPb3=0.05;
	double Dchi2clMinScan_Max_PbPb3=0.45;

	double DdlsMinScan_Min_PbPb3=2.5;
	double DdlsMinScan_Max_PbPb3=6.5;

	double DalphaMaxScan_Min_PbPb3=0.08;
	double DalphaMaxScan_Max_PbPb3=0.2;

	double Dchi2clMinScan_bins_PbPb3[nbin_Dchi2clMinScan];
	double DdlsMinScan_bins_PbPb3[nbin_DdlsMinScan];
	double DalphaMaxScan_bins_PbPb3[nbin_DalphaMaxScan];

	
	double bins_Dchi2clMinScan_PbPb3[nbin_Dchi2clMinScan+1];
	double bins_DdlsMinScan_PbPb3[nbin_DdlsMinScan+1];
	double bins_DalphaMaxScan_PbPb3[nbin_DalphaMaxScan+1];


/*
  const int nPtBins = 9;
  Double_t ptBins[nPtBins+1] = {2.,4.,6.,8.,10.,12.,20.,40.,60.,100};



Double_t ffls3dcut[nPtBins] = {6.00,  5.86,  5.86,  4.86,  4.54,  4.42,  4.06,  3.50,  3.00};
Double_t vprobcut[nPtBins] =  {0.250, 0.224, 0.224, 0.170, 0.125, 0.091, 0.069, 0.054, 0.050};
const int nPhiBins = 5;


const int nCentBins = 6;
Int_t centBins[nCentBins+1] = {0, 10, 30, 50, 70, 90, 100};
Double_t EPm_resolution_v2_etagap[nCentBins] = {0.685732, 0.859684, 0.805492, 0.566930, 0.211378, 0.0307577};
Double_t EPp_resolution_v2_etagap[nCentBins] = {0.685895, 0.859866, 0.805762, 0.567147, 0.210694, 0.0329058};

Double_t minhisto=1.7;
Double_t maxhisto=2.0;
Double_t nbinsmasshisto=60;
Double_t binwidthmass=(maxhisto-minhisto)/nbinsmasshisto;
const int nMassBins = 60;
Double_t fstMassBin = minhisto;
Double_t widMassBin = binwidthmass;
Double_t massBins[nMassBins+1];

const int nDcaBins = 12;
Double_t fstDcaBin = 0.001;
Double_t widDcaBin = 1.27;
//Double_t dcaBins[nDcaBins+1] = {0, 0.001, 0.002, 0.003, 0.004, 0.006, 0.008, 0.012, 0.016, 0.024, 0.032, 0.045, 0.070};
Double_t dcaBins[nDcaBins+1] = {0, 0.001, 0.00227, 0.0038829, 0.00593128, 0.008, 0.0118366, 0.0160324, 0.0213612, 0.0281287, 0.0367235, 0.0476388, 0.07};
const int nD0Bins = 20;
Double_t fstD0Bin = 3.5;
Double_t widD0Bin = 5;
Double_t d0Bins[nD0Bins+1];

TString tfname[3][2] = {{"v2_inpl","v2_outpl"},{"v3_inpl","v3_outpl"},{"inclusive",""}};

const int nDcaxyBins = 20;
Double_t DcaxyBins[nDcaxyBins+1] = {-0.0734,-0.0562,-0.0428,-0.0320,-0.0236,-0.0170,-0.0118,-0.0078,-0.0046,-0.002,0.0,0.002,0.0046,0.0078,0.0118,0.0170,0.0236,0.0320,0.0428,0.0562,0.0734};

const int nFonllBins = 400;
Double_t fstFonllBins = 0;
Double_t lstFonllBins = 100;
Double_t widFonllBins = (lstFonllBins-fstFonllBins)/nFonllBins;

  const int nBinZ = 50;
  Double_t binsZ[nBinZ+1];
*/
/*
void initBins(bool resetDca=false)
{
  if(resetDca)
    {
      dcaBins[0] = 0;
      for(int i=1;i<nDcaBins+1;i++) dcaBins[i] = dcaBins[i-1]+fstDcaBin*pow(widDcaBin,i-1);
    }
  for(int i=0;i<nMassBins+1;i++) massBins[i] = fstMassBin+i*widMassBin;
  for(int i=0;i<nD0Bins+1;i++) d0Bins[i] = fstD0Bin+i*widD0Bin;

  float firstBinZ = -75;
  float binWidthZ = 3;
  for(int i=0; i<=nBinZ; i++){
    binsZ[i] = firstBinZ+binWidthZ*i;
	}

}
*/
void DrawCmsTlatex(TString collision, Float_t tsize=0.04)
{
  TLatex* texCms = new TLatex(0.18,0.93, "#scale[1.25]{CMS} #bf{#it{Preliminary}}");
  texCms->SetNDC();
  texCms->SetTextAlign(12);
  texCms->SetTextSize(tsize);
  texCms->SetTextFont(42);
  texCms->Draw();

  TLatex* texCol = new TLatex(0.96,0.93, Form("%s #sqrt{s_{NN}} = 5.02 TeV",collision.Data()));
  texCol->SetNDC();
  texCol->SetTextAlign(32);
  texCol->SetTextSize(tsize);
  texCol->SetTextFont(42);
  texCol->Draw();
}

void initParameter(){



	for(int i=0; i<nbin_Dchi2clMinScan; i++){
		Dchi2clMinScan_bins_pp[i]=Dchi2clMinScan_Min_pp+(double)i*(Dchi2clMinScan_Max_pp-Dchi2clMinScan_Min_pp)/(double)(nbin_Dchi2clMinScan-1);
		Dchi2clMinScan_bins_PbPb1[i]=Dchi2clMinScan_Min_PbPb1+(double)i*(Dchi2clMinScan_Max_PbPb1-Dchi2clMinScan_Min_PbPb1)/(double)(nbin_Dchi2clMinScan-1);
		Dchi2clMinScan_bins_PbPb2[i]=Dchi2clMinScan_Min_PbPb2+(double)i*(Dchi2clMinScan_Max_PbPb2-Dchi2clMinScan_Min_PbPb2)/(double)(nbin_Dchi2clMinScan-1);
		Dchi2clMinScan_bins_PbPb3[i]=Dchi2clMinScan_Min_PbPb3+(double)i*(Dchi2clMinScan_Max_PbPb3-Dchi2clMinScan_Min_PbPb3)/(double)(nbin_Dchi2clMinScan-1);

		bins_Dchi2clMinScan_pp[i]=Dchi2clMinScan_bins_pp[i]-0.5*(Dchi2clMinScan_Max_pp-Dchi2clMinScan_Min_pp)/(double)(nbin_Dchi2clMinScan-1);
		bins_Dchi2clMinScan_PbPb3[i]=Dchi2clMinScan_bins_PbPb3[i]-0.5*(Dchi2clMinScan_Max_PbPb3-Dchi2clMinScan_Min_PbPb3)/(double)(nbin_Dchi2clMinScan-1);		
	}
		bins_Dchi2clMinScan_pp[nbin_Dchi2clMinScan]=Dchi2clMinScan_bins_pp[nbin_Dchi2clMinScan-1]+0.5*(Dchi2clMinScan_Max_pp-Dchi2clMinScan_Min_pp)/(double)(nbin_Dchi2clMinScan-1);
		bins_Dchi2clMinScan_PbPb3[nbin_Dchi2clMinScan]=Dchi2clMinScan_bins_PbPb3[nbin_Dchi2clMinScan-1]+0.5*(Dchi2clMinScan_Max_PbPb3-Dchi2clMinScan_Min_PbPb3)/(double)(nbin_Dchi2clMinScan-1);		



	for(int i=0; i<nbin_DdlsMinScan; i++){
		DdlsMinScan_bins_pp[i]=DdlsMinScan_Min_pp+(double)i*(DdlsMinScan_Max_pp-DdlsMinScan_Min_pp)/(double)(nbin_DdlsMinScan-1);
		DdlsMinScan_bins_PbPb1[i]=DdlsMinScan_Min_PbPb1+(double)i*(DdlsMinScan_Max_PbPb1-DdlsMinScan_Min_PbPb1)/(double)(nbin_DdlsMinScan-1);
		DdlsMinScan_bins_PbPb2[i]=DdlsMinScan_Min_PbPb2+(double)i*(DdlsMinScan_Max_PbPb2-DdlsMinScan_Min_PbPb2)/(double)(nbin_DdlsMinScan-1);
		DdlsMinScan_bins_PbPb3[i]=DdlsMinScan_Min_PbPb3+(double)i*(DdlsMinScan_Max_PbPb3-DdlsMinScan_Min_PbPb3)/(double)(nbin_DdlsMinScan-1);

		bins_DdlsMinScan_pp[i]=DdlsMinScan_bins_pp[i]-0.5*(DdlsMinScan_Max_pp-DdlsMinScan_Min_pp)/(double)(nbin_DdlsMinScan-1);
		bins_DdlsMinScan_PbPb3[i]=DdlsMinScan_bins_PbPb3[i]-0.5*(DdlsMinScan_Max_PbPb3-DdlsMinScan_Min_PbPb3)/(double)(nbin_DdlsMinScan-1);		
	}
		bins_DdlsMinScan_pp[nbin_DdlsMinScan]=DdlsMinScan_bins_pp[nbin_DdlsMinScan-1]+0.5*(DdlsMinScan_Max_pp-DdlsMinScan_Min_pp)/(double)(nbin_DdlsMinScan-1);
		bins_DdlsMinScan_PbPb3[nbin_DdlsMinScan]=DdlsMinScan_bins_PbPb3[nbin_DdlsMinScan-1]+0.5*(DdlsMinScan_Max_PbPb3-DdlsMinScan_Min_PbPb3)/(double)(nbin_DdlsMinScan-1);		



	for(int i=0; i<nbin_DalphaMaxScan; i++){
		DalphaMaxScan_bins_pp[i]=DalphaMaxScan_Min_pp+(double)i*(DalphaMaxScan_Max_pp-DalphaMaxScan_Min_pp)/(double)(nbin_DalphaMaxScan-1);
		DalphaMaxScan_bins_PbPb1[i]=DalphaMaxScan_Min_PbPb1+(double)i*(DalphaMaxScan_Max_PbPb1-DalphaMaxScan_Min_PbPb1)/(double)(nbin_DalphaMaxScan-1);
		DalphaMaxScan_bins_PbPb2[i]=DalphaMaxScan_Min_PbPb2+(double)i*(DalphaMaxScan_Max_PbPb2-DalphaMaxScan_Min_PbPb2)/(double)(nbin_DalphaMaxScan-1);
		DalphaMaxScan_bins_PbPb3[i]=DalphaMaxScan_Min_PbPb3+(double)i*(DalphaMaxScan_Max_PbPb3-DalphaMaxScan_Min_PbPb3)/(double)(nbin_DalphaMaxScan-1);


		bins_DalphaMaxScan_pp[i]=DalphaMaxScan_bins_pp[i]-0.5*(DalphaMaxScan_Max_pp-DalphaMaxScan_Min_pp)/(double)(nbin_DalphaMaxScan-1);
		bins_DalphaMaxScan_PbPb3[i]=DalphaMaxScan_bins_PbPb3[i]-0.5*(DalphaMaxScan_Max_PbPb3-DalphaMaxScan_Min_PbPb3)/(double)(nbin_DalphaMaxScan-1);		
	}
		bins_DalphaMaxScan_pp[nbin_DalphaMaxScan]=DalphaMaxScan_bins_pp[nbin_DalphaMaxScan-1]+0.5*(DalphaMaxScan_Max_pp-DalphaMaxScan_Min_pp)/(double)(nbin_DalphaMaxScan-1);
		bins_DalphaMaxScan_PbPb3[nbin_DalphaMaxScan]=DalphaMaxScan_bins_PbPb3[nbin_DalphaMaxScan-1]+0.5*(DalphaMaxScan_Max_PbPb3-DalphaMaxScan_Min_PbPb3)/(double)(nbin_DalphaMaxScan-1);		



  char buffer[100];
  fstream file;
  file.open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/TMVA/readxml/TMVA_result_isPbPb0.txt",ios::in);
  int ibin=0;
  int jbin=0;
  if(!file){
    cout<<"failed in open file"<<endl;
		return;
  }else{
    while(!file.eof()){
    file.getline(buffer,sizeof(buffer));
    if(strcmp(buffer,"Dchi2cl")==0){
      file.getline(buffer,sizeof(buffer));
      Dchi2clMin_bins_pp[ibin]=atof(buffer);
			if(verbose){
			cout<<"ibin : "<< ibin<<" , Dchi2clMin_bins_pp = "<<atof(buffer)<<endl;
			}
      ibin++;
    }
    if(strcmp(buffer,"Ddls")==0){
      file.getline(buffer,sizeof(buffer));
      DdlsMin_bins_pp[jbin]=atof(buffer);
			if(verbose){
			cout<<"jbin : "<< jbin<<" , DdlsMin_bins_pp = "<<atof(buffer)<<endl;
			}
      jbin++;
      }
    }
    file.close();
		cout<<"reading TMVA pp cut done"<<endl;
  }

  fstream file_PbPb3;
  file_PbPb3.open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/TMVA/readxml/TMVA_result_isPbPb3.txt",ios::in);
  ibin=0;
  jbin=0;
  if(!file_PbPb3){
    cout<<"failed in open file"<<endl;
		return;
  }else{
    while(!file_PbPb3.eof()){
    file_PbPb3.getline(buffer,sizeof(buffer));
    if(strcmp(buffer,"Dchi2cl")==0){
      file_PbPb3.getline(buffer,sizeof(buffer));
      Dchi2clMin_bins_PbPb3[ibin]=atof(buffer);
			if(verbose){
			cout<<"ibin : "<< ibin<<" , Dchi2clMin_bins_PbPb3 = "<<atof(buffer)<<endl;
			}

      ibin++;
    }
    if(strcmp(buffer,"Ddls")==0){
      file_PbPb3.getline(buffer,sizeof(buffer));
      DdlsMin_bins_PbPb3[jbin]=atof(buffer);
			if(verbose){
			cout<<"jbin : "<< jbin<<" , DdlsMin_bins_PbPb3 = "<<atof(buffer)<<endl;
			}
      jbin++;
      }
    }
    file_PbPb3.close();
		cout<<"reading TMVA PbPb cut done"<<endl;
  }

	// initialize dca bin
	bins_dca[0]=0;
	for(int i=1; i<=nbin_dca; i++){
		bins_dca[i]=bins_dca[i-1]+firstBinDcaWidth*pow(binDcaWidthRatio,i-1);
		cout<<"bin i : "<<i<<" , dca = "<< bins_dca[i]<<endl;
	}
	cout<<"bin dca init done"<<endl;


	return;

}




#endif
