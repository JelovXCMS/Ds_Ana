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

#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/Miscellaneous/sPlot_varCompare/varCompare_para.h"

using namespace std;






int BuildFitfile(TString infile="", TString outfile="", Int_t isPbPb=0, double Dpt_Low=0, double Dpt_High=0){

	Bool_t REAL=true;

	if(Dpt_High==0 || Dpt_High<=Dpt_Low ){
		cout<<"using Default Dpt setting"<<endl;
		if(isPbPb==0){	Dpt_Low=Dpt_Low_pp ; Dpt_High=Dpt_Hight_pp;}
		else if(isPbPb>=0){   Dpt_Low=Dpt_Low_PbPb ; Dpt_High=Dpt_Hight_PbPb;}
	}


	Float_t Dalpha_cut=Dalpha_cut_pp;
	Float_t Dchi2cl_cut=Dchi2cl_cut_pp;
	Float_t Ddls_cut=Ddls_cut_pp;

	Float_t Ddls_maxcut=Ddls_maxcut_pp;
	Float_t Ddls_mincut=1.5;
	Float_t Dalpha_maxcut=Dalpha_maxcut_pp;
	Float_t Dchi2cl_mincut=Dchi2cl_mincut_pp;
	

	initParameter();
  double *bins_pt=bins_pt_pp;
  int nbin_pt=nbin_pt_pp;
  double *Dchi2clMinScan_bins=Dchi2clMinScan_bins_pp;
  double *DalphaMaxScan_bins=DalphaMaxScan_bins_pp;
  double *DdlsMinScan_bins=DdlsMinScan_bins_pp;
  TString Str_isPbPb="pp";
	
	// TString fin_name="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_Data/DsMinTree_pp_Data_MBAll.root";
//	if(!REAL) {fin_name="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_MergeWeight/root_output/DsMinTree_MC_GenSampleMerge_pp_Prompt_phi.root";}

	if(isPbPb==3) {
    Str_isPbPb="PbPb3";
    bins_pt=bins_pt_PbPb3;
    nbin_pt=nbin_pt_PbPb3;
    // TrkptAcc=TrkptAcc_PbPb3;
    // cuthiBin=cuthiBin_PbPb3;
    Dchi2clMinScan_bins=Dchi2clMinScan_bins_PbPb3;
    DalphaMaxScan_bins=DalphaMaxScan_bins_PbPb3;
    DdlsMinScan_bins=DdlsMinScan_bins_PbPb3;

		Dalpha_cut=Dalpha_cut_PbPb3;
		Dchi2cl_cut=Dchi2cl_cut_PbPb3;
		Ddls_cut=Ddls_cut_PbPb3;

		Ddls_maxcut=Ddls_maxcut_PbPb3;
		Ddls_mincut=Ddls_mincut_PbPb3;
		Dalpha_maxcut=Dalpha_maxcut_PbPb3;
		Dchi2cl_mincut=Dchi2cl_mincut_PbPb3;
	

    // infileData=PbPb3DataFile;

  }

	Float_t Dmass;
	Float_t Ddls;
	Float_t Dalpha;
	Float_t Dchi2cl;
	Float_t Ddca;

	TFile *fin=TFile::Open(infile.Data(),"read");
	TTree *ntDsData=(TTree*)fin->Get("ntDs");

	DsMinTreeLoad DsData;
	DsData.SetBranch(ntDsData,REAL,isPbPb);


	TFile *fout=TFile::Open(outfile.Data(),"recreate");

	TTree *t_forDdls=new TTree("t_forDdls","t_forDdls");
	t_forDdls->Branch("Dmass",&Dmass);
	t_forDdls->Branch("Ddls",&Ddls);

	TTree *t_forDalpha=new TTree("t_forDalpha","t_forDalpha");
	t_forDalpha->Branch("Dmass",&Dmass);
	t_forDalpha->Branch("Dalpha",&Dalpha);

	TTree *t_forDchi2cl=new TTree("t_forDchi2cl","t_forDchi2cl");
	t_forDchi2cl->Branch("Dmass",&Dmass);
	t_forDchi2cl->Branch("Dchi2cl",&Dchi2cl);


//	cout<<"bins_pt[ibin_Dpt] = "<<bins_pt[ibin_Dpt]<<endl;
//	return 1;

	Long64_t nentries= ntDsData->GetEntries();

	for(Long64_t i=0; i<nentries; i++){
		ntDsData->GetEntry(i);
    if(i%200000==0) {cout<<setw(10)<<i<<" / "<<nentries<<endl;}
    Dmass=DsData.Dmass;
    Ddls=DsData.Ddls;
		Dalpha=DsData.Dalpha;
		Dchi2cl=DsData.Dchi2cl;
		if(DsData.Dpt>Dpt_Low && DsData.Dpt<Dpt_High && DsData.Dmass>1.90 && DsData.Dmass<2.12 && DsData.DtktkResmass>1.0105 && DsData.DtktkResmass < 1.0285){
		// general cuts 1
			if(Dalpha < Dalpha_maxcut && Dchi2cl > Dchi2cl_mincut && Ddls <Ddls_maxcut && Ddls > Ddls_mincut){
			// general cuts 2 

			if(DsData.Dalpha < Dalpha_cut && DsData.Dchi2cl > Dchi2cl_cut){
				t_forDdls->Fill();
			}

			if(DsData.Dchi2cl > Dchi2cl_cut && Ddls >Ddls_cut){
				t_forDalpha->Fill();
			}

			if(DsData.Dalpha < Dalpha_cut && Ddls >Ddls_cut){
				t_forDchi2cl->Fill();
			}

			}// end if general cut2 
		}// end if general cut



	}// end for loop ientries

	fout->cd();

	t_forDdls->Write();
	t_forDalpha->Write();
	t_forDchi2cl->Write();


	fout->Close();

	return 0;

}


int main(int argc, char *argv[]){

    if(argc==4)
    {
			BuildFitfile(argv[1], argv[2],atoi(argv[3]));
		}else if(argc==6){
			BuildFitfile(argv[1], argv[2],atoi(argv[3]),atof(argv[4]),atof(argv[5]) );
		}else{
	
			cout<<"input format error ,terminate"<<endl;
			return 1;
		}


	return 0;

}

