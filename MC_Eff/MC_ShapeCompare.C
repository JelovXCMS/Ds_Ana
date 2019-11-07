#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"


void MC_ShapeCompare(){

	InitStyle();
	initParameter();

	TFile *f_Prompt_phi= TFile::Open("./output/MC_eff_pp_Prompt_phikkpi.root");
	TH1D *h_Prompt_phi_ptbin[nbin_pt_pp];

	TFile *f_Prompt_f0= TFile::Open("./output/MC_eff_pp_Prompt_f0kkpi.root");
	TH1D *h_Prompt_f0_ptbin[nbin_pt_pp];

	TFile *f_NonPrompt_phi= TFile::Open("./output/MC_eff_pp_NonPrompt_phikkpi.root");
	TH1D *h_NonPrompt_phi_ptbin[nbin_pt_pp];

	for(int i=0;i<nbin_pt_pp; i++){
		h_Prompt_phi_ptbin[i]=(TH1D*)f_Prompt_phi->Get(Form("h_DsMass_pt%.0fto%.0f",bins_pt_pp[i],bins_pt_pp[i+1]));
		h_Prompt_f0_ptbin[i]=(TH1D*)f_Prompt_f0->Get(Form("h_DsMass_pt%.0fto%.0f",bins_pt_pp[i],bins_pt_pp[i+1]));
		h_NonPrompt_phi_ptbin[i]=(TH1D*)f_NonPrompt_phi->Get(Form("h_DsMass_pt%.0fto%.0f",bins_pt_pp[i],bins_pt_pp[i+1]));
		Normalize({h_Prompt_phi_ptbin[i],h_NonPrompt_phi_ptbin[i], h_Prompt_f0_ptbin[i]});
		Draw({h_Prompt_phi_ptbin[i],h_NonPrompt_phi_ptbin[i], h_Prompt_f0_ptbin[i]});
		DrawCompare(h_Prompt_phi_ptbin[i],h_NonPrompt_phi_ptbin[i]);
		DrawCompare(h_Prompt_phi_ptbin[i],h_Prompt_f0_ptbin[i]);
		TCanvas *c_test=Draw({h_Prompt_phi_ptbin[i],h_NonPrompt_phi_ptbin[i], h_Prompt_f0_ptbin[i]});
		SavePlotDirs(c_test,"MCShape",{"MC_study","MC_shape","test"}); // works, need to change save name
	}

	// later using fit to compare fraction of Gauss, width ratio


}
