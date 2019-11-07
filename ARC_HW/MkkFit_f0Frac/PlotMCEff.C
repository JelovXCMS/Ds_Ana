#include "MkkFit_para.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"

int PlotMCEff(int isPbPb=3){

	TString str_PbPb="pp";
	TString s_ppPbPb="pp";
	if(isPbPb==3){
		str_PbPb="PbPb3";
		s_ppPbPb="PbPb";
}


	TString fname_P_phi=Form("../../MC_Eff/output%s/MC_eff_%s_Prompt_phikkpi.root", s_CutSet.Data(),str_PbPb.Data());
	TFile *f_P_phi=TFile::Open(fname_P_phi.Data());
  TH1D *h_RecoNormEff_P_phi=(TH1D*)f_P_phi->Get("h_RecoNormEff");
	h_RecoNormEff_P_phi->Draw();
  TH1D *h_RecoNormEffBr_P_phi=(TH1D*)h_RecoNormEff_P_phi->Clone("h_RecoNormEffBr_P_phi");
  h_RecoNormEffBr_P_phi->Scale(BRphi);

	// h_RecoNormEffBr_P_phi->Draw();	


	TString fname_P_f0=Form("../../MC_Eff/output/MC_eff_%s_Prompt_f0kkpi.root", str_PbPb.Data());
	TFile *f_P_f0=TFile::Open(fname_P_f0.Data());
  TH1D *h_RecoNormEff_P_f0=(TH1D*)f_P_f0->Get("h_RecoNormEff");
	h_RecoNormEff_P_f0->Draw();
  TH1D *h_RecoNormEffBr_P_f0=(TH1D*)h_RecoNormEff_P_f0->Clone("h_RecoNormEffBr_P_f0");
  h_RecoNormEffBr_P_f0->Scale(BRf0);

	// h_RecoNormEffBr_P_f0->Draw("same");	

	TH1D *h_RecoNormEffBr_P_all=(TH1D*)h_RecoNormEffBr_P_phi->Clone("h_RecoNormEffBr_P_phiRatio");
	h_RecoNormEffBr_P_all->Add(h_RecoNormEffBr_P_f0);

  TH1D *h_RecoNormEffBr_P_phiRatio=(TH1D*)h_RecoNormEffBr_P_phi->Clone("h_RecoNormEffBr_P_phiRatio");
	h_RecoNormEffBr_P_phiRatio->Divide(h_RecoNormEffBr_P_all);


	TCanvas *ctest=new TCanvas();
	ctest->cd();
	gStyle->SetOptStat(0);
	if(isPbPb==3){
		h_RecoNormEffBr_P_phiRatio->GetXaxis()->SetRangeUser(6,40);
	}
	h_RecoNormEffBr_P_phiRatio->SetTitle("");
	h_RecoNormEffBr_P_phiRatio->GetXaxis()->SetTitle("D_{S} p_{T} (GeV/c)");
	h_RecoNormEffBr_P_phiRatio->GetYaxis()->SetTitle("ratio");
	h_RecoNormEffBr_P_phiRatio->Draw();

	TLatex *tlatex = new TLatex();	
	tlatex->DrawLatexNDC(textposx,textposy,Form("%s Ds #phi channel ratio",s_ppPbPb.Data()));


  gSystem->Exec("mkdir -p plots/MCphiRatio");

	ctest->SaveAs(Form("plots/MCphiRatio/%s_MCphiRatio_P_prompt.png",s_ppPbPb.Data()));

	return 0;
}
