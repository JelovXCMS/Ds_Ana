#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"

int plotVar_PNP(){

	gStyle->SetOptStat(0);
	gSystem->Exec("mkdir -p plots_VarPNP");

	TFile *f_Prompt=TFile::Open(MCFile_Prompt_phikkpi_pp_Merge.Data());
	TTree *nt_Prompt=(TTree*)f_Prompt->Get("ntDs");



	TH1D *h_Prompt_Ddls=new TH1D("h_Prompt_Ddls",";Decay Length Significance",100,0,20);

	nt_Prompt->Project("h_Prompt_Ddls","Ddls","weight*DgenptSampleWeight*RecoD0DataWeight*(Ddls<20 && DsGen==23333 && DgencollisionId==0&& Dpt>6 && Dpt<8)");

	TFile *f_NonPrompt=TFile::Open(MCFile_NonPrompt_phikkpi_pp_Merge.Data());
	TTree *nt_NonPrompt=(TTree*)f_NonPrompt->Get("ntDs");


	TH1D *h_NonPrompt_Ddls=new TH1D("h_NonPrompt_Ddls",";Decay Length Significance",100,0,20);

	nt_NonPrompt->Project("h_NonPrompt_Ddls","Ddls","weight*DgenptSampleWeight*RecoD0DataWeight*(Ddls<20 && DsGen==23333 && DgencollisionId==0&& Dpt>6 && Dpt<8)");


	TCanvas *c_PNP_Ddls=new TCanvas("c_PNP_Ddls","c_PNP_Ddls") ;
	c_PNP_Ddls->cd();

	// h_Prompt_Ddls->Scale(1/h_Prompt_Ddls->Integral());
	h_Prompt_Ddls->SetLineColor(2);
	h_Prompt_Ddls->SetMarkerColor(2);
	h_Prompt_Ddls->SetMarkerStyle(22);
	h_Prompt_Ddls->Draw();


	// h_NonPrompt_Ddls->Scale(1/h_NonPrompt_Ddls->Integral());
	h_NonPrompt_Ddls->SetLineColor(4);
	h_NonPrompt_Ddls->SetMarkerColor(4);
	h_NonPrompt_Ddls->SetMarkerStyle(23);
	h_NonPrompt_Ddls->Draw("same");

	TLegend *le=new TLegend(0.65,0.65,0.85,0.85);
	le->SetBorderSize(0);
	le->AddEntry((TObject*)0,"pp 6<p_{T}<8 GeV/c","");
	le->AddEntry(h_Prompt_Ddls,"Prompt","l");
	le->AddEntry(h_NonPrompt_Ddls,"NonPrompt","l");
	le->Draw("same");

	c_PNP_Ddls->SaveAs("plots_VarPNP/Ddls_pp_pt6to8.png");

	return 0;

}
