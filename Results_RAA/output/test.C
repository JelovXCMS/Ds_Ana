int test(){

	TFile *f_old=new TFile("./old_noScale/RAA.root");
	TFile *f_new=new TFile("./old_DdlErrScale/RAA.root");

	TH1D *h_CS_pp_old=(TH1D*)f_old->Get("h_PromptDs_PbPb3");
	TH1D *h_CS_pp_new=(TH1D*)f_new->Get("h_PromptDs_PbPb3");

	TCanvas *c_CS_pp=new TCanvas("c_CS_pp","",600,600);
	gPad->SetLogy();
	gStyle->SetOptStat(0);
	c_CS_pp->cd();
	h_CS_pp_old->SetTitle("");
	h_CS_pp_old->GetXaxis()->SetTitle("p_{T}");
	h_CS_pp_old->GetYaxis()->SetTitle("dN/dp_{T}");
	h_CS_pp_old->SetLineColor(1);
	h_CS_pp_old->SetMarkerColor(1);
	h_CS_pp_old->Draw();
	h_CS_pp_new->SetLineColor(2);
	h_CS_pp_new->SetMarkerColor(2);
	h_CS_pp_new->Draw("same");

	TLegend *le=new TLegend(0.5,0.6,0.8,0.8);
	le->SetBorderSize(0);
	le->AddEntry(h_CS_pp_old,"no scale","p");
	le->AddEntry(h_CS_pp_new,"Best DdlErr scale","p");
	le->Draw("same");



	return 1;

}
