int test(){

	TFile *f_old=TFile::Open("./MC_eff_PbPb3_Prompt_phikkpi.root","read");
	TFile *f_new=TFile::Open("./old_DdlErrScale/MC_eff_PbPb3_Prompt_phikkpi.root","read");
	
	TH1D *h_eff_old=(TH1D*)f_old->Get("h_RecoNormEff");
	TH1D *h_eff_new=(TH1D*)f_new->Get("h_RecoNormEff");

	TCanvas *c_eff=new TCanvas("c_eff","c_eff",600,600);
	c_eff->cd();
	h_eff_old->Draw();
	h_eff_new->SetLineColor(2);
	h_eff_new->Draw("same");


	return 0;
}
