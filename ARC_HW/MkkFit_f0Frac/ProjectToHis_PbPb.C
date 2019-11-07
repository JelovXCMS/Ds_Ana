
#include "MkkFit_para.h"

int ProjectToHis_PbPb(int isPbPb=3, double DptLow=6, double DptHigh=40, double DVtxPCut=0.2, double DdlsCut=4.5, TString finName="/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/GJ/Ds_PbPb_Data_HIMB2_GJ.root" , TString foutName="PbPb_MB2_GJ"){

	TString s_ppPbPb="pp";

	// int bins_phi=

  gSystem->Exec("mkdir -p his");
  gSystem->Exec("mkdir -p plots");
  gSystem->Exec("mkdir -p plots/Ori");

	TString s_finName="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_Data/DsMinTree_pp_Data_MBAll.root";

	if(isPbPb){
	//	s_finName="/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/GJ/Ds_PbPb_Data_HIMB2_GJ.root";
		s_finName=finName;
		s_ppPbPb="PbPb";
	}

	// TString s_foutName=Form("./his/%s_Dpt%.0fto%.0f.root",s_ppPbPb.Data(),DptLow,DptHigh);
	TString s_foutName=Form("./his/%s_Dpt%.0fto%.0f.root",foutName.Data(),DptLow,DptHigh);

	TString s_fMC_Prompt_Phi=MCFile_Prompt_phikkpi_pp_Merge;
	TString s_fMC_Prompt_f0=MCFile_Prompt_f0980kkpi_pp_Merge;
	
	if(isPbPb){
		s_fMC_Prompt_Phi=MCFile_Prompt_phikkpi_PbPb_Merge;
		s_fMC_Prompt_f0=MCFile_Prompt_f0980kkpi_PbPb_Merge;
	}
	

	// TString s_ppPbPb="pp";

	TFile *fin=TFile::Open(s_finName.Data(),"read");
	TTree *ntDs=(TTree*)fin->Get("ntDs");

	TFile *f_MCP_Phi=TFile::Open(s_fMC_Prompt_Phi.Data());
	TTree *ntDs_MCP_Phi=(TTree*)f_MCP_Phi->Get("ntDs");
	TFile *f_MCP_f0=TFile::Open(s_fMC_Prompt_f0.Data());
	TTree *ntDs_MCP_f0=(TTree*)f_MCP_f0->Get("ntDs");

	TFile *fout=TFile::Open(s_foutName.Data(),"recreate");
	fout->cd();


	if(foutName=="PbPb_MB2_GJ")
	{
/*

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

	Float_t DsGen;
	Float_t DgencollisionId;
	Float_t DgenBAncestorpt;

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
  ntDs_MCP_Phi->SetBranchAddress("DdlErr",&DdlErr);
  ntDs_MCP_Phi->SetBranchAddress("Dalpha",&Dalpha);
  ntDs_MCP_Phi->SetBranchAddress("Dchi2cl",&Dchi2cl);
  ntDs_MCP_Phi->SetBranchAddress("Dtrk1Pt",&Dtrk1Pt);
  ntDs_MCP_Phi->SetBranchAddress("Dtrk2Pt",&Dtrk2Pt);
  ntDs_MCP_Phi->SetBranchAddress("Dtrk3Pt",&Dtrk3Pt);

  ntDs_MCP_Phi->SetBranchAddress("DsGen",&DsGen);
  ntDs_MCP_Phi->SetBranchAddress("DgencollisionId",&DgencollisionId);
  ntDs_MCP_Phi->SetBranchAddress("DgenBAncestorpt",&DgenBAncestorpt);

  ntDs_MCP_Phi->SetBranchAddress("weight",&weight);
  ntDs_MCP_Phi->SetBranchAddress("DgenptSampleWeight",&DgenptSampleWeight);
  ntDs_MCP_Phi->SetBranchAddress("RecoFONLLWeight",&RecoFONLLWeight);
  ntDs_MCP_Phi->SetBranchAddress("RecoD0DataWeight",&RecoD0DataWeight);
  ntDs_MCP_Phi->SetBranchAddress("RecoDsDataWeight",&RecoDsDataWeight);
  ntDs_MCP_Phi->SetBranchAddress("PbPbVzWeight",&PbPbVzWeight);
  ntDs_MCP_Phi->SetBranchAddress("RecoFONLLRaaWeight",&RecoFONLLRaaWeight);

	fout->cd();
	TTree *t_fit=new TTree("t_fit","t_fit");
  t_fit->Branch("Dmass",&Dmass);
  t_fit->Branch("Dpt",&Dpt);
  t_fit->Branch("Dalpha",&Dalpha);
  t_fit->Branch("Dchi2cl",&Dchi2cl);
  t_fit->Branch("Ddls",&Ddls);
  t_fit->Branch("Ddl",&Ddl);
  t_fit->Branch("DdlErr",&DdlErr);
  t_fit->Branch("DtktkResmass",&DtktkResmass);
	
    t_fit->Branch("Dalpha_BS_2D",&Dalpha_BS_2D);
    t_fit->Branch("DlxyBS",&DlxyBS);
    t_fit->Branch("DlxyBSErr",&DlxyBSErr);
    t_fit->Branch("DlxyBSs",&DlxyBSs);

    t_fit->Branch("TotalWeight",&TotalWeight);

*/

// create TH1D


	TH1D *h_mkkpi_default_MCP_Phi=new TH1D("h_mkkpi_default_MCP_Phi","; m_{KK#pi}",nbin_DsMass*2, DsMassFitLow,DsMassFitHigh);
	h_mkkpi_default_MCP_Phi->Sumw2();

	ntDs_MCP_Phi->Project(Form("h_mkkpi_default_MCP_Phi"),"Dmass",Form("PbPbVzWeight*weight*DgenptSampleWeight*RecoD0DataWeight*(Dmass>%f && Dmass<%f && Dpt>%f && Dpt<%f && Dalpha<0.12 && Dchi2cl>%f && Ddls >%f && DtktkResmass>%f && DtktkResmass<%f)",DsMassFitLow,DsMassFitHigh, DptLow, DptHigh, DVtxPCut, DdlsCut, DtktkResmassCutMean-DtktkResmassCutWidth, DtktkResmassCutMean+DtktkResmassCutWidth ) );

	h_mkkpi_default_MCP_Phi->Write();

	TCanvas *ctest=new TCanvas("ctest","ctest");
  ctest->cd();
	h_mkkpi_default_MCP_Phi->Draw();
	ctest->SaveAs(Form("plots/Ori/OriHis_MCP_Phi_%s_pt%.0fto%.0f_mkkDefault.png",s_ppPbPb.Data(),DptLow,DptHigh));

	TH1D *h_mkkpi_default_MCP_f0=new TH1D("h_mkkpi_default_MCP_f0","; m_{KK#pi}",nbin_DsMass*2, DsMassFitLow,DsMassFitHigh);
	h_mkkpi_default_MCP_f0->Sumw2();

	ntDs_MCP_f0->Project(Form("h_mkkpi_default_MCP_f0"),"Dmass",Form("PbPbVzWeight*weight*DgenptSampleWeight*RecoD0DataWeight*(Dmass>%f && Dmass<%f && Dpt>%f && Dpt<%f && Dalpha<0.12 && Dchi2cl>%f && Ddls >%f && DtktkResmass>%f && DtktkResmass<%f)",DsMassFitLow,DsMassFitHigh, DptLow, DptHigh, DVtxPCut, DdlsCut, DtktkResmassCutMean-DtktkResmassCutWidth, DtktkResmassCutMean+DtktkResmassCutWidth ) );

	h_mkkpi_default_MCP_f0->Write();

  ctest->cd();
	h_mkkpi_default_MCP_f0->Draw();
	ctest->SaveAs(Form("plots/Ori/OriHis_MCP_f0_%s_pt%.0fto%.0f_mkkDefault.png",s_ppPbPb.Data(),DptLow,DptHigh));

	// mkk his for template fit

	TH1D *h_mkk_MCP_Phi=new TH1D("h_mkk_MCP_Phi",";m_{KK}",nbin_mkk,bins_mkk);
	h_mkk_MCP_Phi->Sumw2();
	ntDs_MCP_Phi->Project("h_mkk_MCP_Phi","DtktkResmass",Form("PbPbVzWeight*weight*DgenptSampleWeight*RecoD0DataWeight*(Dmass>%f && Dmass<%f && Dpt>%f && Dpt<%f && Dalpha<0.12 && Dchi2cl>%f && Ddls >%f && DtktkResmass>%f && DtktkResmass<%f && DsGen==23333)",DsMassMCMean-DsMass2SigWidth , DsMassMCMean+DsMass2SigWidth, DptLow, DptHigh, DVtxPCut, DdlsCut, bins_mkk[0], bins_mkk[nbin_mkk] ) );
	h_mkk_MCP_Phi->Write();
	
	ctest->cd();
	h_mkk_MCP_Phi->Draw();
	ctest->SaveAs(Form("plots/Ori/OriHis_Mkk_MCP_Phi_%s_pt%.0fto%.0f.png",s_ppPbPb.Data(),DptLow,DptHigh));

	TH1D *h_mkk_MCP_f0=new TH1D("h_mkk_MCP_f0",";m_{KK}",nbin_mkk,bins_mkk);
	h_mkk_MCP_f0->Sumw2();
	ntDs_MCP_f0->Project("h_mkk_MCP_f0","DtktkResmass",Form("PbPbVzWeight*weight*DgenptSampleWeight*RecoD0DataWeight*(Dmass>%f && Dmass<%f && Dpt>%f && Dpt<%f && Dalpha<0.12 && Dchi2cl>%f && Ddls >%f && DtktkResmass>%f && DtktkResmass<%f && DsGen==24433)",DsMassMCMean-DsMass2SigWidth , DsMassMCMean+DsMass2SigWidth, DptLow, DptHigh, DVtxPCut, DdlsCut, bins_mkk[0], bins_mkk[nbin_mkk] ) );
	h_mkk_MCP_f0->Write();
	
	ctest->cd();
	h_mkk_MCP_f0->Draw();
	ctest->SaveAs(Form("plots/Ori/OriHis_Mkk_MCP_f0_%s_pt%.0fto%.0f.png",s_ppPbPb.Data(),DptLow,DptHigh));


	TH1D *h_mkk_MCP_Phi_fineBin=new TH1D("h_mkk_MCP_Phi_fineBin",";m_{KK}",nbin_mkk*2,bin_mkkLow,bin_mkkHigh);
	h_mkk_MCP_Phi_fineBin->Sumw2();
	ntDs_MCP_Phi->Project("h_mkk_MCP_Phi_fineBin","DtktkResmass",Form("PbPbVzWeight*weight*DgenptSampleWeight*RecoD0DataWeight*(Dmass>%f && Dmass<%f && Dpt>%f && Dpt<%f && Dalpha<0.12 && Dchi2cl>%f && Ddls >%f && DtktkResmass>%f && DtktkResmass<%f && DsGen==23333)",DsMassMCMean-DsMass2SigWidth , DsMassMCMean+DsMass2SigWidth, DptLow, DptHigh, DVtxPCut, DdlsCut, bins_mkk[0], bins_mkk[nbin_mkk] ) );
	h_mkk_MCP_Phi_fineBin->Write();
	
	ctest->cd();
	h_mkk_MCP_Phi_fineBin->Draw();
	ctest->SaveAs(Form("plots/Ori/OriHis_Mkk_MCP_Phi_fineBin_%s_pt%.0fto%.0f.png",s_ppPbPb.Data(),DptLow,DptHigh));

	TH1D *h_mkk_MCP_f0_fineBin=new TH1D("h_mkk_MCP_f0_fineBin",";m_{KK}",nbin_mkk*2,bin_mkkLow,bin_mkkHigh);
	h_mkk_MCP_f0_fineBin->Sumw2();
	ntDs_MCP_f0->Project("h_mkk_MCP_f0_fineBin","DtktkResmass",Form("PbPbVzWeight*weight*DgenptSampleWeight*RecoD0DataWeight*(Dmass>%f && Dmass<%f && Dpt>%f && Dpt<%f && Dalpha<0.12 && Dchi2cl>%f && Ddls >%f && DtktkResmass>%f && DtktkResmass<%f && DsGen==24433)",DsMassMCMean-DsMass2SigWidth , DsMassMCMean+DsMass2SigWidth, DptLow, DptHigh, DVtxPCut, DdlsCut, bins_mkk[0], bins_mkk[nbin_mkk] ) );
	h_mkk_MCP_f0_fineBin->Write();
	
	ctest->cd();
	h_mkk_MCP_f0_fineBin->Draw();
	ctest->SaveAs(Form("plots/Ori/OriHis_Mkk_MCP_f0_fineBin_%s_pt%.0fto%.0f.png",s_ppPbPb.Data(),DptLow,DptHigh));

	}// only run MC part once , in MB2_GJ


	TH1D *h_mkkpi_mkkbins[nbin_mkk];
	TH1D *h_mkkpi_default=new TH1D("h_mkkpi_default","; m_{KK#pi}",nbin_DsMass, DsMassFitLow,DsMassFitHigh);
	h_mkkpi_default->Sumw2();

	ntDs->Project(Form("h_mkkpi_default"),"Dmass",Form("Dmass>%f && Dmass<%f && Dpt>%f && Dpt<%f && Dalpha<0.12 && Dchi2cl>%f && Ddls >%f && DtktkResmass>%f && DtktkResmass<%f",DsMassFitLow,DsMassFitHigh, DptLow, DptHigh, DVtxPCut, DdlsCut, DtktkResmassCutMean-DtktkResmassCutWidth, DtktkResmassCutMean+DtktkResmassCutWidth ) );

	h_mkkpi_default->Write();

	TCanvas *ctest=new TCanvas("ctest","ctest");
  ctest->cd();
	h_mkkpi_default->Draw();
	ctest->SaveAs(Form("plots/Ori/OriHis_%s_pt%.0fto%.0f_mkkDefault.png",s_ppPbPb.Data(),DptLow,DptHigh));

	delete ctest;	

	for(int i=0; i<nbin_mkk ; i++){

		cout<<"project bin : "<<i<<endl;
		h_mkkpi_mkkbins[i]=new TH1D(Form("h_mkkpi_mkkbins_%i",i),"; m_{KK#pi}", nbin_DsMass , DsMassFitLow, DsMassFitHigh );
		h_mkkpi_mkkbins[i]->Sumw2();
		ntDs->Project(Form("h_mkkpi_mkkbins_%i",i),"Dmass",Form("Dmass>%f && Dmass<%f && Dpt>%f && Dpt<%f && Dalpha<0.12 && Dchi2cl>%f && Ddls>%f && DtktkResmass>%f && DtktkResmass<%f",DsMassFitLow,DsMassFitHigh, DptLow, DptHigh, DVtxPCut,DdlsCut , bins_mkk[i], bins_mkk[i+1] ) );
		TCanvas *ctest=new TCanvas("ctest","ctest");
		ctest->cd();
		h_mkkpi_mkkbins[i]->Draw();

		// ctest->SaveAs(Form("plots/Ori/OriHis_%i.pdf",i));
		ctest->SaveAs(Form("plots/Ori/OriHis_%s_pt%.0fto%.0f_mkk%.0fto%.0f.png",s_ppPbPb.Data(),DptLow,DptHigh,bins_mkk[i]*1000,bins_mkk[i+1]*1000));

		h_mkkpi_mkkbins[i]->Write();

		delete ctest;

	}

	fout->Close();



	return 0;
}



int main(int argc, char*argv[]){


	ProjectToHis_PbPb(atoi(argv[1]), atof(argv[2]) , atof(argv[3]) , atof(argv[4]) , atof(argv[5]) , argv[6] , argv[7]  );


	return 0;
}
