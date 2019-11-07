#include "../include/uti.h"
#include "../include/DsMinTreeLoad.h"
#include "../include/parameters.h"

#include <TString.h>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

void MC_dca(int isPbPb=0){

	initParameter();

	TFile *f_Prompt_phi = TFile::Open("./output/MC_eff_PbPb3_Prompt_phikkpi.root","read");
	TTree *t_Dmass_pt6to8_Prompt_phi=(TTree*)f_Prompt_phi->Get("t_DsMass_pt6to8");

	TH1F *h_dca_pt6to8_Prompt_phi=new TH1F("h_dca_pt6to8_Prompt_phi","h_dca_pt6to8_Prompt_phi",nbin_dca,bins_dca); 
	h_dca_pt6to8_Prompt_phi->Sumw2();
	t_Dmass_pt6to8_Prompt_phi->Project("h_dca_pt6to8_Prompt_phi","Ddca","BasicWeight");
	h_dca_pt6to8_Prompt_phi->Scale(1/h_dca_pt6to8_Prompt_phi->Integral());

	TFile *f_NonPrompt_phi = TFile::Open("./output/MC_eff_PbPb3_NonPrompt_phikkpi.root","read");
	TTree *t_Dmass_pt6to8_NonPrompt_phi=(TTree*)f_NonPrompt_phi->Get("t_DsMass_pt6to8");

	TH1F *h_dca_pt6to8_NonPrompt_phi=new TH1F("h_dca_pt6to8_NonPrompt_phi","h_dca_pt6to8_NonPrompt_phi",nbin_dca,bins_dca); 
	h_dca_pt6to8_NonPrompt_phi->Sumw2();
	t_Dmass_pt6to8_NonPrompt_phi->Project("h_dca_pt6to8_NonPrompt_phi","Ddca","BasicWeight");
	h_dca_pt6to8_NonPrompt_phi->Scale(1/h_dca_pt6to8_NonPrompt_phi->Integral());

	TCanvas *c_dca_PNP=new TCanvas("c_dca_PNP","c_dca_PNP",800,800);
	c_dca_PNP->cd();
	h_dca_pt6to8_Prompt_phi->SetLineColor(2);
	h_dca_pt6to8_Prompt_phi->Draw("PLC PMC");
	h_dca_pt6to8_NonPrompt_phi->Draw("SAME PLC PMC");
	gPad->BuildLegend();
	c_dca_PNP->SaveAs("plots/pdf/dna_PNP_example.pdf");



}
