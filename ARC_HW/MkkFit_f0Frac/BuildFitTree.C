
#include "MkkFit_para.h"

int BuildFitTree(int isPbPb=0,int isMC=1){

	TString s_fMC_Prompt_Phi=MCFile_Prompt_phikkpi_pp_Merge;
	TString s_ppPbPb="pp";
	TString s_MCData="Data";
	if(isMC){
		s_MCData="MC";
	}
	if(isPbPb){
		s_fMC_Prompt_Phi=MCFile_Prompt_phikkpi_PbPb_Merge;
		s_ppPbPb="PbPb";
	}
	TFile *f_MCP_Phi=TFile::Open(s_fMC_Prompt_Phi.Data());
	TTree *ntDs_MCP_Phi=(TTree*)f_MCP_Phi->Get("ntDs");

	gSystem->Exec("mkdir -p /scratch/halstead/p/peng43/Ds_phikkpi/FitFile_Mkk/");

	TFile *fout=TFile::Open(Form("/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_Mkk/%s_FitTree_%s.root",s_ppPbPb.Data(),s_MCData.Data()),"recreate");
	fout->cd();
	// TTree *ntGen=(TTree*)f_MCP_Phi->Get("ntGen");
	// TTree *ntGen_new=(TTree*)ntGen->CloneTree();
	// ntGen_new->Write();


	// create TTree
	Float_t Dmass;
	Float_t Dpt;
	Float_t Dy;
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

	Int_t DsGen;
	Int_t DgencollisionId;
	Float_t DgenBAncestorpt;

	Float_t Ncoll;

	Float_t weight;
	Float_t DgenptSampleWeight;
	Float_t PbPbVzWeight;
	Float_t RecoFONLLWeight;
	Float_t RecoFONLLRaaWeight;
	Float_t RecoD0DataWeight;
	Float_t RecoDsDataWeight;

	Float_t TotalWeight;

	ntDs_MCP_Phi->SetBranchAddress("Dmass",&Dmass);
	ntDs_MCP_Phi->SetBranchAddress("Dpt",&Dpt);
	ntDs_MCP_Phi->SetBranchAddress("Dy",&Dy);
	ntDs_MCP_Phi->SetBranchAddress("Ddls",&Ddls);
	// ntDs_MCP_Phi->SetBranchAddress("Ddl",&Ddl);
	ntDs_MCP_Phi->SetBranchAddress("DdlErr",&DdlErr);
	ntDs_MCP_Phi->SetBranchAddress("Dalpha",&Dalpha);
	ntDs_MCP_Phi->SetBranchAddress("Dchi2cl",&Dchi2cl);
	ntDs_MCP_Phi->SetBranchAddress("Dtrk1Pt",&Dtrk1Pt);
	ntDs_MCP_Phi->SetBranchAddress("Dtrk2Pt",&Dtrk2Pt);
	ntDs_MCP_Phi->SetBranchAddress("Dtrk3Pt",&Dtrk3Pt);
	ntDs_MCP_Phi->SetBranchAddress("DtktkResmass",&DtktkResmass);


	// ntDs_MCP_Phi->SetBranchAddress("Dalpha_BS_2D",&Dalpha_BS_2D);
	// ntDs_MCP_Phi->SetBranchAddress("DlxyBS",&DlxyBS);
	// ntDs_MCP_Phi->SetBranchAddress("DlxyBSs",&DlxyBSs);
	// ntDs_MCP_Phi->SetBranchAddress("DlxyBSErr",&DlxyBSErr);

	if(isPbPb){
	ntDs_MCP_Phi->SetBranchAddress("Ncoll",&Ncoll);
	}

	if(isMC){
	ntDs_MCP_Phi->SetBranchAddress("DsGen",&DsGen);
	ntDs_MCP_Phi->SetBranchAddress("DgencollisionId",&DgencollisionId);
	ntDs_MCP_Phi->SetBranchAddress("DgenBAncestorpt",&DgenBAncestorpt);

	ntDs_MCP_Phi->SetBranchAddress("weight",&weight);
	ntDs_MCP_Phi->SetBranchAddress("DgenptSampleWeight",&DgenptSampleWeight);
	ntDs_MCP_Phi->SetBranchAddress("RecoFONLLWeight",&RecoFONLLWeight);
	ntDs_MCP_Phi->SetBranchAddress("RecoD0DataWeight",&RecoD0DataWeight);
	ntDs_MCP_Phi->SetBranchAddress("RecoDsDataWeight",&RecoDsDataWeight);

	if(isPbPb){
		ntDs_MCP_Phi->SetBranchAddress("PbPbVzWeight",&PbPbVzWeight);
		ntDs_MCP_Phi->SetBranchAddress("RecoFONLLRaaWeight",&RecoFONLLRaaWeight);
	}
	}

	fout->cd();
	TTree *t_fit=new TTree("t_fit","t_fit");
	t_fit->Branch("Dmass",&Dmass);
	t_fit->Branch("Dpt",&Dpt);
	t_fit->Branch("Dy",&Dy);
	t_fit->Branch("Dalpha",&Dalpha);
	t_fit->Branch("Dchi2cl",&Dchi2cl);
	t_fit->Branch("Ddls",&Ddls);
	// t_fit->Branch("Ddl",&Ddl);
	t_fit->Branch("DdlErr",&DdlErr);
	if(isPbPb){
		t_fit->Branch("Ncoll",&Ncoll);
	}
	t_fit->Branch("DtktkResmass",&DtktkResmass);

	// t_fit->Branch("Dalpha_BS_2D",&Dalpha_BS_2D);
	// t_fit->Branch("DlxyBS",&DlxyBS);
	// t_fit->Branch("DlxyBSErr",&DlxyBSErr);
	// t_fit->Branch("DlxyBSs",&DlxyBSs);

	t_fit->Branch("TotalWeight",&TotalWeight);

	Long64_t nEntries=ntDs_MCP_Phi->GetEntries();

	for(int i =0; i<nEntries; i++){
			ntDs_MCP_Phi->GetEntry(i);

			if(isMC){			
		  TotalWeight=weight*DgenptSampleWeight*RecoFONLLWeight*RecoD0DataWeight*RecoDsDataWeight;
			if(isPbPb){
		  TotalWeight=Ncoll*weight*DgenptSampleWeight*RecoFONLLRaaWeight*RecoD0DataWeight*RecoDsDataWeight*PbPbVzWeight;
			}
			}

			t_fit->Fill();


	}

	fout->cd(); 
	t_fit->Write();
	fout->Close();


	return 1;

}

