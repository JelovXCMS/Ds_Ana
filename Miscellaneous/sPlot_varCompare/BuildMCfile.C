/*
	with three kinds of weight to :
	1.build MC fit file (tree with three kinds of weight)
	2.MC efficiency for different weights
	3.store var distribtuion with assigend binning (?) --> just from tree , build later


*/

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


int BuildMCfile(Int_t isPbPb=0,int PNPrompt=1,int DsChannel=1 , double Dpt_Low=0, double Dpt_High=0){


  Bool_t REAL=false;

  if(Dpt_High==0 || Dpt_High<=Dpt_Low ){
    cout<<"using Default Dpt setting"<<endl;
    if(isPbPb==0){  Dpt_Low=Dpt_Low_pp ; Dpt_High=Dpt_Hight_pp;}
    else if(isPbPb>=0){   Dpt_Low=Dpt_Low_PbPb ; Dpt_High=Dpt_Hight_PbPb;}
  }
  initParameter();

  TString Str_PbPb="pp";
  if(isPbPb){Str_PbPb="PbPb";}
  TString Str_PNPrompt="Prompt";
  if(PNPrompt==1) {Str_PNPrompt="NonPrompt";}

  Int_t DsGenTrue=23333;
  Int_t GSignalTypeTrue=1;

  TString Str_DsChannel="phikkpi";
  if(DsChannel==1){
    Str_DsChannel="f0kkpi" ;
    DsGenTrue=24433;
    GSignalTypeTrue=2;
  }else if(DsChannel==2){
    Str_DsChannel="kstarkkpi" ;
    DsGenTrue=25544;
    GSignalTypeTrue=3;
  }



  Float_t Dalpha_cut=Dalpha_cut_pp;
  Float_t Dchi2cl_cut=Dchi2cl_cut_pp;
  Float_t Ddls_cut=Ddls_cut_pp;

  Float_t Ddls_maxcut=Ddls_maxcut_pp;
  Float_t Dalpha_maxcut=Dalpha_maxcut_pp;
  Float_t Dchi2cl_mincut=Dchi2cl_mincut_pp;

  Float_t weight;
  Float_t DgenptSampleWeight;
  Float_t PbPbVzWeight;
  Float_t RecoFONLLWeight;
  Float_t RecoFONLLRaaWeight;
  Float_t RecoD0DataWeight;


  TString MCFile;
  if(isPbPb==0 && PNPrompt==0 && DsChannel==0){MCFile=MCFile_Prompt_phikkpi_pp_Merge;  }
  if(isPbPb==0 && PNPrompt==1 && DsChannel==0){MCFile=MCFile_NonPrompt_phikkpi_pp_Merge;  }
  if(isPbPb==0 && PNPrompt==0 && DsChannel==1){MCFile=MCFile_Prompt_f0980kkpi_pp_Merge;  }
  if(isPbPb==0 && PNPrompt==1 && DsChannel==1){MCFile=MCFile_NonPrompt_f0980kkpi_pp_Merge;  }

  if(isPbPb && PNPrompt==0 && DsChannel==0){MCFile=MCFile_Prompt_phikkpi_PbPb_Merge;  }
  if(isPbPb && PNPrompt==1 && DsChannel==0){MCFile=MCFile_NonPrompt_phikkpi_PbPb_Merge;  }
  if(isPbPb && PNPrompt==0 && DsChannel==1){MCFile=MCFile_Prompt_f0980kkpi_PbPb_Merge;  }
  if(isPbPb && PNPrompt==1 && DsChannel==1){MCFile=MCFile_NonPrompt_f0980kkpi_PbPb_Merge;  }

  TString inputFile=MCFile;

  TFile *fin=TFile::Open(inputFile.Data());

  TTree *ntDs=(TTree*)fin->Get("ntDs");
  TTree *ntGen=(TTree*)fin->Get("ntGen");

  TCut cutGenTrue=Form("GSignalType==%i && GcollisionId==0 && TMath::Abs(Gy)<1",GSignalTypeTrue);
  TCut cutGenPNPrompt="GBAncestorpt<=0";
  if(PNPrompt==1){cutGenPNPrompt="GBAncestorpt>0";}

  TCut cuthiBin="";
  TString GenWeight_Pythia="weight*GptSampleWeight";
  TString GenWeight_FONLL="weight*GptSampleWeight*GenFONLLWeight";
  TString GenWeight_D0Data="weight*GptSampleWeight*GenFONLLWeight*GenD0DataWeight";
  // TString RecoWeight="weight*DgenptSampleWeight";
  if(isPbPb){
    GenWeight_Pythia="weight*GptSampleWeight*PbPbVzWeight*Ncoll";
    GenWeight_FONLL="weight*GptSampleWeight*PbPbVzWeight*Ncoll*GenFONLLRaaWeight";
    GenWeight_D0Data="weight*GptSampleWeight*PbPbVzWeight*Ncoll*GenFONLLWeight*GenD0DataWeight";
    // RecoWeight="weight*DgenptSampleWeight*PbPbVzWeight*Ncoll";
  }

  TString outfileName=Form("./output/MC_%s_%s_%s.root",Str_PbPb.Data(), Str_PNPrompt.Data(), Str_DsChannel.Data());
	// for fout
  TFile *fout=TFile::Open(outfileName.Data(),"recreate");

	TH1D *h_GenAll = new TH1D("h_GenAll","h_GenAll",1,Dpt_Low,Dpt_High); h_GenAll->Sumw2();
	TH1D *h_GenAll_Pythia = new TH1D("h_GenAll_Pythia","h_GenAll_Pythia",1,Dpt_Low,Dpt_High); h_GenAll_Pythia->Sumw2();
	TH1D *h_GenAll_FONLL = new TH1D("h_GenAll_FONLL","h_GenAll_FONLL",1,Dpt_Low,Dpt_High); h_GenAll_FONLL->Sumw2();

	TH1D *h_Reco_forDdls = new TH1D("h_Reco_forDdls","h_Reco_forDdls",1,Dpt_Low,Dpt_High); h_Reco_forDdls->Sumw2();
	TH1D *h_Reco_forDdls_Pythia = new TH1D("h_Reco_forDdls_Pythia","h_Reco_forDdls_Pythia",1,Dpt_Low,Dpt_High); h_Reco_forDdls_Pythia->Sumw2();
	TH1D *h_Reco_forDdls_FONLL = new TH1D("h_Reco_forDdls_FONLL","h_Reco_forDdls_FONLL",1,Dpt_Low,Dpt_High); h_Reco_forDdls_FONLL->Sumw2();

	TH1D *h_Reco_forDalpha = new TH1D("h_Reco_forDalpha","h_Reco_forDalpha",1,Dpt_Low,Dpt_High); h_Reco_forDalpha->Sumw2();
	TH1D *h_Reco_forDalpha_Pythia = new TH1D("h_Reco_forDalpha_Pythia","h_Reco_forDalpha_Pythia",1,Dpt_Low,Dpt_High); h_Reco_forDalpha_Pythia->Sumw2();
	TH1D *h_Reco_forDalpha_FONLL = new TH1D("h_Reco_forDalpha_FONLL","h_Reco_forDalpha_FONLL",1,Dpt_Low,Dpt_High); h_Reco_forDalpha_FONLL->Sumw2();

	TH1D *h_Reco_forDchi2cl = new TH1D("h_Reco_forDchi2cl","h_Reco_forDchi2cl",1,Dpt_Low,Dpt_High); h_Reco_forDchi2cl->Sumw2();
	TH1D *h_Reco_forDchi2cl_Pythia = new TH1D("h_Reco_forDchi2cl_Pythia","h_Reco_forDchi2cl_Pythia",1,Dpt_Low,Dpt_High); h_Reco_forDchi2cl_Pythia->Sumw2();
	TH1D *h_Reco_forDchi2cl_FONLL = new TH1D("h_Reco_forDchi2cl_FONLL","h_Reco_forDchi2cl_FONLL",1,Dpt_Low,Dpt_High); h_Reco_forDchi2cl_FONLL->Sumw2();

// eff histo

 	TH1D *h_RecoEff_forDdls = new TH1D("h_RecoEff_forDdls","h_RecoEff_forDdls",1,Dpt_Low,Dpt_High); h_RecoEff_forDdls->Sumw2();
	TH1D *h_RecoEff_forDdls_Pythia = new TH1D("h_RecoEff_forDdls_Pythia","h_RecoEff_forDdls_Pythia",1,Dpt_Low,Dpt_High); h_RecoEff_forDdls_Pythia->Sumw2();
	TH1D *h_RecoEff_forDdls_FONLL = new TH1D("h_RecoEff_forDdls_FONLL","h_RecoEff_forDdls_FONLL",1,Dpt_Low,Dpt_High); h_RecoEff_forDdls_FONLL->Sumw2();

	TH1D *h_RecoEff_forDalpha = new TH1D("h_RecoEff_forDalpha","h_RecoEff_forDalpha",1,Dpt_Low,Dpt_High); h_RecoEff_forDalpha->Sumw2();
	TH1D *h_RecoEff_forDalpha_Pythia = new TH1D("h_RecoEff_forDalpha_Pythia","h_RecoEff_forDalpha_Pythia",1,Dpt_Low,Dpt_High); h_RecoEff_forDalpha_Pythia->Sumw2();
	TH1D *h_RecoEff_forDalpha_FONLL = new TH1D("h_RecoEff_forDalpha_FONLL","h_RecoEff_forDalpha_FONLL",1,Dpt_Low,Dpt_High); h_RecoEff_forDalpha_FONLL->Sumw2();

	TH1D *h_RecoEff_forDchi2cl = new TH1D("h_RecoEff_forDchi2cl","h_RecoEff_forDchi2cl",1,Dpt_Low,Dpt_High); h_RecoEff_forDchi2cl->Sumw2();
	TH1D *h_RecoEff_forDchi2cl_Pythia = new TH1D("h_RecoEff_forDchi2cl_Pythia","h_RecoEff_forDchi2cl_Pythia",1,Dpt_Low,Dpt_High); h_RecoEff_forDchi2cl_Pythia->Sumw2();
	TH1D *h_RecoEff_forDchi2cl_FONLL = new TH1D("h_RecoEff_forDchi2cl_FONLL","h_RecoEff_forDchi2cl_FONLL",1,Dpt_Low,Dpt_High); h_RecoEff_forDchi2cl_FONLL->Sumw2();



  ntGen->Project("h_GenAll","Gpt",(TCut)((cutGenTrue && cutGenPNPrompt && cuthiBin)*GenWeight_D0Data.Data()) );
  ntGen->Project("h_GenAll_Pythia","Gpt",(TCut)((cutGenTrue && cutGenPNPrompt && cuthiBin)*GenWeight_Pythia.Data()) );
  ntGen->Project("h_GenAll_FONLL","Gpt",(TCut)((cutGenTrue && cutGenPNPrompt && cuthiBin)*GenWeight_FONLL.Data()) );



	// for tree

  DsMinTreeLoad DsMCReco;
  DsMCReco.SetBranch(ntDs, REAL, isPbPb); // false for MC, to read more branch info

  ntDs->SetBranchAddress("weight",&weight);
  ntDs->SetBranchAddress("DgenptSampleWeight",&DgenptSampleWeight);
  ntDs->SetBranchAddress("RecoFONLLWeight",&RecoFONLLWeight);
  ntDs->SetBranchAddress("RecoD0DataWeight",&RecoD0DataWeight);
  if(isPbPb){
  ntDs->SetBranchAddress("PbPbVzWeight",&PbPbVzWeight);
  ntDs->SetBranchAddress("RecoFONLLRaaWeight",&RecoFONLLRaaWeight);
  }

  cout<<"check4"<<endl;

  Long64_t nentries_Reco=ntDs->GetEntries();

	Float_t Dpt;
	Float_t Dmass;
  Float_t Ddls;
  Float_t Dalpha;
  Float_t Dchi2cl;
  Float_t Ddca;

  Float_t D0DataWeight;
  Float_t FONLLWeight;
  Float_t PythiaWeight;

	fout->cd();

  TTree *t_forDdls=new TTree("t_forDdls","t_forDdls");
  t_forDdls->Branch("Dmass",&Dmass);
  t_forDdls->Branch("Ddls",&Ddls);
	t_forDdls->Branch("D0DataWeight",&D0DataWeight);
	t_forDdls->Branch("FONLLWeight",&FONLLWeight);
	t_forDdls->Branch("PythiaWeight",&PythiaWeight);
	
  TTree *t_forDalpha=new TTree("t_forDalpha","t_forDalpha");
  t_forDalpha->Branch("Dmass",&Dmass);
  t_forDalpha->Branch("Dalpha",&Dalpha);
	t_forDalpha->Branch("D0DataWeight",&D0DataWeight);
	t_forDalpha->Branch("FONLLWeight",&FONLLWeight);
	t_forDalpha->Branch("PythiaWeight",&PythiaWeight);

  TTree *t_forDchi2cl=new TTree("t_forDchi2cl","t_forDchi2cl");
  t_forDchi2cl->Branch("Dmass",&Dmass);
  t_forDchi2cl->Branch("Dchi2cl",&Dchi2cl);
	t_forDchi2cl->Branch("D0DataWeight",&D0DataWeight);
	t_forDchi2cl->Branch("FONLLWeight",&FONLLWeight);
	t_forDchi2cl->Branch("PythiaWeight",&PythiaWeight);


  cout<<"nentries_Reco = "<<nentries_Reco<<endl;

  for(Long64_t ientry=0 ; ientry<nentries_Reco; ientry++){
    if(ientry%50000==0){
      cout<<setw(9)<<ientry<<" / "<<nentries_Reco<<endl;
    }

    ntDs->GetEntry(ientry);
 		Dpt=DsMCReco.Dpt;
    Dmass=DsMCReco.Dmass;
    Ddca=DsMCReco.Ddca;
    Ddls=DsMCReco.Ddls;
    Dalpha=DsMCReco.Dalpha;
    Dchi2cl=DsMCReco.Dchi2cl;

    // cout<<"Dmass = "<<Dmass<<endl;
    D0DataWeight=weight*DgenptSampleWeight*RecoFONLLWeight*RecoD0DataWeight;
    FONLLWeight=weight*DgenptSampleWeight*RecoFONLLWeight;
    PythiaWeight=weight*DgenptSampleWeight;
    // cout<<"weight= "<<weight<<endl;
    // cout<<"DgenptSampleWeight = "<<DgenptSampleWeight<<endl;
    // D0DataWeight=1;

    // cout<<"D0DataWeight = "<<D0DataWeight<<endl;

    if(isPbPb){
      D0DataWeight=weight*DgenptSampleWeight*PbPbVzWeight*DsMCReco.Ncoll*RecoFONLLWeight*RecoD0DataWeight;
      FONLLWeight=weight*DgenptSampleWeight*PbPbVzWeight*DsMCReco.Ncoll*RecoFONLLRaaWeight;
      PythiaWeight=weight*DgenptSampleWeight*PbPbVzWeight*DsMCReco.Ncoll;
    }

//		cout<<"Dpt = "<<DsMCReco.Dpt<<" , "<<endl;

    if(DsMCReco.Dpt>Dpt_Low && DsMCReco.Dpt < Dpt_High && DsMCReco.Dmass>1.90 && DsMCReco.Dmass<2.12 && DsMCReco.DtktkResmass>1.0105 && DsMCReco.DtktkResmass < 1.0285 && DsMCReco.DsGen==DsGenTrue && DsMCReco.DgencollisionId==0 && (( PNPrompt==0 && DsMCReco.DgenBAncestorpt<=0 ) || ( PNPrompt==1 && DsMCReco.DgenBAncestorpt>0 ))  ){
    // general cuts 1
			
      if(Dalpha < Dalpha_maxcut && Dchi2cl > Dchi2cl_mincut && Ddls <Ddls_maxcut){
      // general cuts 2

			// cout<<"inside g cut2 "<<endl;

      if(DsMCReco.Dalpha < Dalpha_cut && DsMCReco.Dchi2cl > Dchi2cl_cut){
        t_forDdls->Fill();
				h_Reco_forDdls->Fill(Dpt,D0DataWeight);
				h_Reco_forDdls_Pythia->Fill(Dpt,PythiaWeight);
				h_Reco_forDdls_FONLL->Fill(Dpt,FONLLWeight);
      }

      if(DsMCReco.Dchi2cl > Dchi2cl_cut && Ddls >Ddls_cut){
        t_forDalpha->Fill();
				h_Reco_forDalpha->Fill(Dpt,D0DataWeight);
				h_Reco_forDalpha_Pythia->Fill(Dpt,PythiaWeight);
				h_Reco_forDalpha_FONLL->Fill(Dpt,FONLLWeight);

      }

      if(DsMCReco.Dalpha < Dalpha_cut && Ddls >Ddls_cut){
        t_forDchi2cl->Fill();
				h_Reco_forDchi2cl->Fill(Dpt,D0DataWeight);
				h_Reco_forDchi2cl_Pythia->Fill(Dpt,PythiaWeight);
				h_Reco_forDchi2cl_FONLL->Fill(Dpt,FONLLWeight);

      }

      }// end if general cut2
    }// end if general cut


	}// end for loop


	// calculate eff


	h_RecoEff_forDdls->Divide(h_Reco_forDdls, h_GenAll, 1,1,"B");
	h_RecoEff_forDdls_Pythia->Divide(h_Reco_forDdls_Pythia, h_GenAll_Pythia, 1,1,"B");
	h_RecoEff_forDdls_FONLL->Divide(h_Reco_forDdls_FONLL, h_GenAll_FONLL, 1,1,"B");

	h_RecoEff_forDalpha->Divide(h_Reco_forDalpha, h_GenAll, 1,1,"B");
	h_RecoEff_forDalpha_Pythia->Divide(h_Reco_forDalpha_Pythia, h_GenAll_Pythia, 1,1,"B");
	h_RecoEff_forDalpha_FONLL->Divide(h_Reco_forDalpha_FONLL, h_GenAll_FONLL, 1,1,"B");

	h_RecoEff_forDchi2cl->Divide(h_Reco_forDchi2cl, h_GenAll, 1,1,"B");
	h_RecoEff_forDchi2cl_Pythia->Divide(h_Reco_forDchi2cl_Pythia, h_GenAll_Pythia, 1,1,"B");
	h_RecoEff_forDchi2cl_FONLL->Divide(h_Reco_forDchi2cl_FONLL, h_GenAll_FONLL, 1,1,"B");

  fout->cd();


	h_GenAll->Write();
	h_GenAll_Pythia->Write();
	h_GenAll_FONLL->Write();

	h_Reco_forDdls->Write();
	h_Reco_forDdls_Pythia->Write();
	h_Reco_forDdls_FONLL->Write();

	h_Reco_forDalpha->Write();
	h_Reco_forDalpha_Pythia->Write();
	h_Reco_forDalpha_FONLL->Write();

	h_Reco_forDchi2cl->Write();
	h_Reco_forDchi2cl_Pythia->Write();
	h_Reco_forDchi2cl_FONLL->Write();

	h_RecoEff_forDdls->Write();
	h_RecoEff_forDdls_Pythia->Write();
	h_RecoEff_forDdls_FONLL->Write();

	h_RecoEff_forDalpha->Write();
	h_RecoEff_forDalpha_Pythia->Write();
	h_RecoEff_forDalpha_FONLL->Write();

	h_RecoEff_forDchi2cl->Write();
	h_RecoEff_forDchi2cl_Pythia->Write();
	h_RecoEff_forDchi2cl_FONLL->Write();




  t_forDdls->Write();
  t_forDalpha->Write();
  t_forDchi2cl->Write();

  return 0;



}
