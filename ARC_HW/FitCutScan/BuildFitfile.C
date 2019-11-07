#include <iostream>
#include <fstream>
#include <iomanip>

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
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/DsMinTreeLoad.h"

#include "varCompare_para.h"

using namespace std;


int BuildFitfile(Int_t isPbPb=3, Int_t isReal=1,int PNPrompt=0,int DsChannel=0, TString s_fin="",TString foutTag="", double Dpt_Low=0, double Dpt_High=0){

	int nbin_pt=nbin_pt_pp;
	double *bins_pt=bins_pt_pp;


	Bool_t REAL=true;

	if(Dpt_High==0 || Dpt_High<=Dpt_Low ){
		cout<<"using Default Dpt setting"<<endl;
		if(isPbPb==0){	Dpt_Low=Dpt_Low_pp ; Dpt_High=Dpt_Hight_pp;}
		else if(isPbPb>=0){ Dpt_Low=Dpt_Low_PbPb ; Dpt_High=Dpt_Hight_PbPb;
		}
	}
  TString Str_PNPrompt="Prompt";
  if(PNPrompt==1) {Str_PNPrompt="NonPrompt";}

  Int_t DsGenTrue=23333;
  Int_t GSignalTypeTrue=1;
  TString Str_DsChannel="phi";
	if(isReal==0){
		foutTag="MC_phi";
		REAL=false;
  if(DsChannel==1){
    Str_DsChannel="f0" ;
    DsGenTrue=24433;
    GSignalTypeTrue=2;
		foutTag="MC_f0";
  }else if(DsChannel==2){
    Str_DsChannel="kstar" ;
    DsGenTrue=25544;
    GSignalTypeTrue=3;
		foutTag="MC_kstar";
  }
		foutTag+=Str_PNPrompt;		
	}

  TCut cutGenTrue=Form("GSignalType==%i && GcollisionId==0 && TMath::Abs(Gy)<1",GSignalTypeTrue);
  TCut cutGenPNPrompt="GBAncestorpt<=0";
	TCut cutGpt=Form("Gpt>%f && Gpt<%f",Dpt_Low,Dpt_High);
  if(PNPrompt==1){cutGenPNPrompt="GBAncestorpt>0";}

  TCut cuthiBin="";
  TString GenWeight_Pythia="weight*GptSampleWeight";
  TString GenWeight_FONLL="weight*GptSampleWeight*GenFONLLWeight";
  TString GenWeight_D0Data="weight*GptSampleWeight*GenFONLLWeight*GenD0DataWeight";
  TString GenWeight_DsData="weight*GptSampleWeight*GenFONLLWeight*GenD0DataWeight*GenDsDataWeight";
  // TString RecoWeight="weight*DgenptSampleWeight";
  if(isPbPb){
    GenWeight_Pythia="weight*GptSampleWeight*PbPbVzWeight*Ncoll";
    GenWeight_FONLL="weight*GptSampleWeight*PbPbVzWeight*Ncoll*GenFONLLRaaWeight";
    GenWeight_D0Data="weight*GptSampleWeight*PbPbVzWeight*Ncoll*GenFONLLWeight*GenD0DataWeight";
    GenWeight_DsData="weight*GptSampleWeight*PbPbVzWeight*Ncoll*GenFONLLWeight*GenD0DataWeight*GenDsDataWeight";
  }




	Float_t DmassLow=1.91;
	Float_t DmassHigh=2.11;
	Float_t PhiMassLow=1.0045;
	Float_t PhiMassHigh=1.0345;


	Float_t Dalpha_cut=Dalpha_cut_pp;
	Float_t Dchi2cl_cut=Dchi2cl_cut_pp;
	Float_t Ddls_cut=Ddls_cut_pp;

	Float_t Ddls_maxcut=Ddls_maxcut_pp;
	Float_t Dalpha_maxcut=Dalpha_maxcut_pp;
	Float_t Dchi2cl_mincut=Dchi2cl_mincut_pp;

	Float_t TrkPtCut=0.7;

	

	initParameter();
  // double *bins_pt=bins_pt_pp;
  // int nbin_pt=nbin_pt_pp;
  double *Dchi2clMinScan_bins=Dchi2clMinScan_bins_pp;
  double *DalphaMaxScan_bins=DalphaMaxScan_bins_pp;
  double *DdlsMinScan_bins=DdlsMinScan_bins_pp;

  TString Str_isPbPb="pp";
	TString fin_name="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_Data/DsMinTree_pp_Data_MBAll.root";
	// if(!REAL) {fin_name="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_MergeWeight/root_output/DsMinTree_MC_GenSampleMerge_pp_Prompt_phi.root";}
 
	if(isPbPb==3) {
    Str_isPbPb="PbPb3";
    bins_pt=bins_pt_PbPb3;
    nbin_pt=nbin_pt_PbPb3;
    // TrkptAcc=TrkptAcc_PbPb3;
    // cuthiBin=cuthiBin_PbPb3;
    Dchi2clMinScan_bins=Dchi2clMinScan_bins_PbPb3;
    DalphaMaxScan_bins=DalphaMaxScan_bins_PbPb3;
    DdlsMinScan_bins=DdlsMinScan_bins_PbPb3;
		TrkPtCut=1;
		fin_name="/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/Ds_PbPb_Data_HIMBAll.root";
		// fin_name="/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/Ds_PbPb_Data_HIMB5_TrkOnly.root";
    // infileData=PbPb3DataFile;
  }

	if(s_fin!=""){
		fin_name=s_fin;
	}

	gSystem->Exec("mkdir -p /scratch/halstead/p/peng43/Ds_phikkpi/FitFile_cutScan");

	TString fout_name=Form("/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_cutScan/%s%s_fitFile.root",Str_isPbPb.Data(),foutTag.Data());


	Float_t Dmass;
	Float_t Dpt;
	Float_t Ddls;
	Float_t Ddl;
	Float_t DdlErr;
	Float_t Dalpha;
	Float_t Dchi2cl;
	// Float_t Ddca;
	Float_t Dalpha_BS_2D;
	Float_t DlxyBS;
	Float_t DlxyBSErr;
	Float_t DlxyBSs;

	Float_t Dtrk1Pt;
	Float_t Dtrk2Pt;
	Float_t Dtrk3Pt;
	Float_t DdxyzErr;
	Float_t DtktkResmass;

  Float_t weight;
  Float_t DgenptSampleWeight;
  Float_t PbPbVzWeight;
  Float_t RecoFONLLWeight;
  Float_t RecoFONLLRaaWeight;
  Float_t RecoD0DataWeight;
  Float_t RecoDsDataWeight;


	Float_t TotalWeight;

  TString MCFile;
  if(isPbPb==0 && PNPrompt==0 && DsChannel==0){MCFile=MCFile_Prompt_phikkpi_pp_Merge;  }
  if(isPbPb==0 && PNPrompt==1 && DsChannel==0){MCFile=MCFile_NonPrompt_phikkpi_pp_Merge;  }
  if(isPbPb==0 && PNPrompt==0 && DsChannel==1){MCFile=MCFile_Prompt_f0980kkpi_pp_Merge;  }
  if(isPbPb==0 && PNPrompt==1 && DsChannel==1){MCFile=MCFile_NonPrompt_f0980kkpi_pp_Merge;  }

  if(isPbPb && PNPrompt==0 && DsChannel==0){MCFile=MCFile_Prompt_phikkpi_PbPb_Merge;  }
  if(isPbPb && PNPrompt==1 && DsChannel==0){MCFile=MCFile_NonPrompt_phikkpi_PbPb_Merge;  }
  if(isPbPb && PNPrompt==0 && DsChannel==1){MCFile=MCFile_Prompt_f0980kkpi_PbPb_Merge;  }
  if(isPbPb && PNPrompt==1 && DsChannel==1){MCFile=MCFile_NonPrompt_f0980kkpi_PbPb_Merge;  }

	if(isReal==0){
		fin_name=MCFile;
	}



	TFile *fin=TFile::Open(fin_name.Data(),"read");
	TTree *ntDsData=(TTree*)fin->Get("ntDs");
	// TTree *ntGen;


	DsMinTreeLoad DsData;
	DsData.SetBranch(ntDsData,REAL,isPbPb);

	if(!REAL){		
  ntDsData->SetBranchAddress("weight",&weight);
  ntDsData->SetBranchAddress("DgenptSampleWeight",&DgenptSampleWeight);
  ntDsData->SetBranchAddress("RecoFONLLWeight",&RecoFONLLWeight);
  ntDsData->SetBranchAddress("RecoD0DataWeight",&RecoD0DataWeight);
  ntDsData->SetBranchAddress("RecoDsDataWeight",&RecoDsDataWeight);
  if(isPbPb){
  ntDsData->SetBranchAddress("PbPbVzWeight",&PbPbVzWeight);
  ntDsData->SetBranchAddress("RecoFONLLRaaWeight",&RecoFONLLRaaWeight);
  }
	}



	TFile *fout=TFile::Open(fout_name.Data(),"recreate");


	TTree *t_fit= new TTree("t_fit","t_fit");
	t_fit->Branch("Dmass",&Dmass);
	t_fit->Branch("Dpt",&Dpt);
	t_fit->Branch("Dalpha",&Dalpha);
	t_fit->Branch("Dchi2cl",&Dchi2cl);
	t_fit->Branch("Ddls",&Ddls);
	t_fit->Branch("Ddl",&Ddl);
	t_fit->Branch("DdlErr",&DdlErr);
	t_fit->Branch("DtktkResmass",&DtktkResmass);
	// t_fit->Branch("Ddl",&Ddl);
	// t_fit->Branch("DdlErr",&DdlErr);
	// t_fit->Branch("Dtrk1Pt",&Dtrk1Pt);
	// t_fit->Branch("Dtrk2Pt",&Dtrk2Pt);
	// t_fit->Branch("Dtrk3Pt",&Dtrk3Pt);
	t_fit->Branch("DdxyzErr",&DdxyzErr);

	if(!isPbPb && !REAL){
		t_fit->Branch("Dalpha_BS_2D",&Dalpha_BS_2D);
		t_fit->Branch("DlxyBS",&DlxyBS);
		t_fit->Branch("DlxyBSErr",&DlxyBSErr);
		t_fit->Branch("DlxyBSs",&DlxyBSs);
	}

	if(!REAL){
		t_fit->Branch("TotalWeight",&TotalWeight);
	}

/*
	TTree *t_forDdls=new TTree("t_forDdls","t_forDdls");
	t_forDdls->Branch("Dmass",&Dmass);
	t_forDdls->Branch("Ddls",&Ddls);

	TTree *t_forDalpha=new TTree("t_forDalpha","t_forDalpha");
	t_forDalpha->Branch("Dmass",&Dmass);
	t_forDalpha->Branch("Dalpha",&Dalpha);

	TTree *t_forDchi2cl=new TTree("t_forDchi2cl","t_forDchi2cl");
	t_forDchi2cl->Branch("Dmass",&Dmass);
	t_forDchi2cl->Branch("Dchi2cl",&Dchi2cl);
*/

//	cout<<"bins_pt[ibin_Dpt] = "<<bins_pt[ibin_Dpt]<<endl;
//	return 1;



	Long64_t nentries= ntDsData->GetEntries();

	for(Long64_t i=0; i<nentries; i++){
		ntDsData->GetEntry(i);
    if(i%200000==0) {cout<<setw(10)<<i<<" / "<<nentries<<endl;}
    Dmass=DsData.Dmass;
    Ddls=DsData.Ddls;
    Ddl=DsData.Ddl;
    DdlErr=DsData.DdlErr;
    DlxyBS=DsData.DlxyBS;
    DlxyBSErr=DsData.DlxyBSErr;
    DlxyBSs=DsData.DlxyBSs;
		Dalpha=DsData.Dalpha;
		Dalpha_BS_2D=DsData.Dalpha_BS_2D;
		Dchi2cl=DsData.Dchi2cl;
		Dpt=DsData.Dpt;
	  Dtrk1Pt=DsData.Dtrk1Pt;
	  Dtrk2Pt=DsData.Dtrk2Pt;
	  Dtrk3Pt=DsData.Dtrk3Pt;
	  DdxyzErr=DsData.DdxyzErr;
		DtktkResmass=DsData.DtktkResmass;
		
		// TotalWeight=weight*DgenptSampleWeight*RecoFONLLWeight*RecoD0DataWeight;
		TotalWeight=weight*DgenptSampleWeight*RecoFONLLWeight*RecoD0DataWeight*RecoDsDataWeight;
		if(isPbPb){TotalWeight=PbPbVzWeight*DsData.Ncoll * weight*DgenptSampleWeight*RecoFONLLWeight*RecoD0DataWeight*RecoDsDataWeight;}

			// cout<<"DsData.Dpt = "<<DsData.Dpt<<" DsData.Dmass = "<<DsData.Dmass<<" DsData.DtktkResmass = "<<DsData.DtktkResmass<<" Dalpha = "<<Dalpha<<" Dchi2cl = "<<Dchi2cl<<" , Ddls = "<<Ddls<<endl;


		// if(DsData.Dpt>Dpt_Low && DsData.Dpt < Dpt_High && DsData.Dmass>DmassLow && DsData.Dmass<DmassHigh && DsData.DtktkResmass>PhiMassLow && DsData.DtktkResmass < PhiMassHigh && Dalpha<Dalpha_cut && Dchi2cl>Dchi2cl_cut && Ddls > Ddls_cut && Dtrk1Pt>TrkPtCut && Dtrk2Pt>TrkPtCut && Dtrk3Pt>TrkPtCut){
		if(DsData.Dpt>Dpt_Low && DsData.Dpt < Dpt_High && DsData.Dmass>DmassLow && DsData.Dmass<DmassHigh && DsData.DtktkResmass>PhiMassLow && DsData.DtktkResmass < PhiMassHigh && Dtrk1Pt>TrkPtCut && Dtrk2Pt>TrkPtCut && Dtrk3Pt>TrkPtCut){
		
			// cout<<"pass cut"<<endl;

			if(isReal){
				if(Dalpha<Dalpha_cut && Dchi2cl>Dchi2cl_cut && Ddls > Ddls_cut){
				t_fit->Fill();	
				}
			}else if(DsData.DsGen==DsGenTrue && DsData.DgencollisionId==0 && (  (PNPrompt==0 && DsData.DgenBAncestorpt<=0) || (PNPrompt==1 && DsData.DgenBAncestorpt>0 )  ) ){
				t_fit->Fill();
			}


 
		}// end if general cut



	} // end for loop ientries

	fout->cd();

//	t_forDdls->Write();
//	t_forDalpha->Write();
//	t_forDchi2cl->Write();

	t_fit->Write("",TObject::kOverwrite);


// ntGen part
	// ntGen part
	if(!REAL){
		TTree *ntGen=(TTree*)fin->Get("ntGen");
		TH1D *hGen_pt=new TH1D("hGen_pt","; Gpt",nbin_pt,bins_pt); hGen_pt->Sumw2();
		ntGen->Project("hGen_pt","Gpt",(TCut)((cutGenTrue && cutGenPNPrompt && cuthiBin && cutGpt)*GenWeight_DsData.Data()) );	
		fout->cd();
		hGen_pt->Write();
		// copy the Gen Tree
		// TTree *ntGen_new=(TTree*)ntGen->CloneTree();
		// ntGen_new->Write("ntGen",TObject::kOverwrite);

	// for quick check MC

	TH1D *hDpt= new TH1D("hDpt","; Dpt",nbin_pt,bins_pt); hDpt->Sumw2();
	t_fit->Project("hDpt","Dpt","TotalWeight");
	hDpt->Write();

	TH1D *hEffDpt=(TH1D*)hDpt->Clone("hEffDpt");
	hEffDpt->Divide(hGen_pt);
	hEffDpt->Write();


	}





	cout<<"fin  : "<<fin_name<<endl;
	cout<<"fout : "<<fout_name<<endl;
	cout<<"cut : "<<endl;
	cout<<"Dpt : "<<Dpt_Low<<" to "<<Dpt_High<<endl;
	cout<<"DMass : "<<DmassLow<<" to "<<DmassHigh<<endl;
	cout<<"PhiMass : "<<PhiMassLow<<" to "<<PhiMassHigh<<endl;
	cout<<"Dalpha : "<<Dalpha_cut<<endl;
	cout<<"Dchi2cl : "<<Dchi2cl_cut<<endl;
	cout<<"Ddls : "<<Ddls_cut<<endl;






	fout->Close();

	return 0;

}
