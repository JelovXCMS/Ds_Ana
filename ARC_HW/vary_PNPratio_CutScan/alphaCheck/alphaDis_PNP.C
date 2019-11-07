#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"


int alphaDis_PNP(int isPbPb=0, double ptLow=6,double ptHigh=10){

  TString s_fMC_Prompt_Phi=MCFile_Prompt_phikkpi_pp_Merge;
  TString s_fMC_Prompt_f0=MCFile_Prompt_f0980kkpi_pp_Merge;
  TString s_fMC_NonPrompt_Phi=MCFile_NonPrompt_phikkpi_pp_Merge;
  TString s_fMC_NonPrompt_f0=MCFile_NonPrompt_f0980kkpi_pp_Merge;


	int nbin_alpha=20;
	double binLow_alpha=0;
	double binHigh_alpha=0.2;

	int nbin_TrkPt=60;
	double binLow_TrkPt=0;
	double binHigh_TrkPt=6;


	TFile *fMCP_Phi=TFile::Open(s_fMC_Prompt_Phi.Data(),"read");
	TTree *nt_MCP_Phi=(TTree*)fMCP_Phi->Get("ntDs");
	TH1D *h_Trk2Pt_MCP_Phi=new TH1D("h_Trk2Pt_MCP_Phi",";Trk2Pt",nbin_TrkPt,binLow_TrkPt,binHigh_TrkPt);
	TH1D *h_Trk3Pt_MCP_Phi=new TH1D("h_Trk3Pt_MCP_Phi",";Trk3Pt",nbin_TrkPt,binLow_TrkPt,binHigh_TrkPt);
	nt_MCP_Phi->Project("h_Trk2Pt_MCP_Phi","Dtrk2Pt",Form("weight*DgenptSampleWeight*RecoD0DataWeight*(DsGen==23333 && Dpt<5)"));
	nt_MCP_Phi->Project("h_Trk3Pt_MCP_Phi","Dtrk3Pt",Form("weight*DgenptSampleWeight*RecoD0DataWeight*(DsGen==23333 && Dpt<5)"));
	// h_Trk2Pt_MCP_Phi->Draw();

	h_Trk2Pt_MCP_Phi->Scale(1/h_Trk2Pt_MCP_Phi->Integral());
	h_Trk3Pt_MCP_Phi->Scale(1/h_Trk3Pt_MCP_Phi->Integral());


	TFile *fMCP_f0=TFile::Open(s_fMC_Prompt_f0.Data(),"read");
	TTree *nt_MCP_f0=(TTree*)fMCP_f0->Get("ntDs");
	TH1D *h_Trk2Pt_MCP_f0=new TH1D("h_Trk2Pt_MCP_f0",";Trk2Pt",nbin_TrkPt,binLow_TrkPt,binHigh_TrkPt);
	TH1D *h_Trk3Pt_MCP_f0=new TH1D("h_Trk3Pt_MCP_f0",";Trk3Pt",nbin_TrkPt,binLow_TrkPt,binHigh_TrkPt);
	nt_MCP_f0->Project("h_Trk2Pt_MCP_f0","Dtrk2Pt",Form("weight*DgenptSampleWeight*RecoD0DataWeight*(DsGen==24433 && Dpt<5 )"));
	nt_MCP_f0->Project("h_Trk3Pt_MCP_f0","Dtrk3Pt",Form("weight*DgenptSampleWeight*RecoD0DataWeight*(DsGen==24433 && Dpt<5 )"));
	
	h_Trk2Pt_MCP_f0->Scale(1/h_Trk2Pt_MCP_f0->Integral());
	h_Trk3Pt_MCP_f0->Scale(1/h_Trk3Pt_MCP_f0->Integral());

	h_Trk2Pt_MCP_f0->SetLineColor(2);
	h_Trk3Pt_MCP_f0->SetLineColor(2);

	// h_Trk2Pt_MCP_f0->Draw();

	TCanvas *c_Trk2Pt=new TCanvas("c_Trk2Pt","c_Trk2Pt",600,600);
	c_Trk2Pt->cd();
	h_Trk2Pt_MCP_Phi->Draw();
	h_Trk2Pt_MCP_f0->Draw("same");

	TCanvas *c_Trk3Pt=new TCanvas("c_Trk3Pt","c_Trk3Pt",600,600);
	c_Trk3Pt->cd();
	h_Trk3Pt_MCP_Phi->Draw();
	h_Trk3Pt_MCP_f0->Draw("same");




	return 0;


}
