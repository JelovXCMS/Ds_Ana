// only fro showing plot, not used by other code


void DsDptShapeRatio(){

	TFile *f_D0=TFile::Open("DptSampleWeight_pp_NonPrompt.root");
	TH1D *h_D0Gpt_weighted=(TH1D*)f_D0->Get("h_Gpt_fromHaoWeight");
	h_D0Gpt_weighted->Scale(1/h_D0Gpt_weighted->Integral());
	TH1D *h_D0Gpt_weightedFineBin=(TH1D*)f_D0->Get("h_Gpt_fromHaoWeightFineBin");
	h_D0Gpt_weightedFineBin->Scale(1/h_D0Gpt_weightedFineBin->Integral());

	TFile *f_Ds=TFile::Open("DsptSampleWeight_pp_NonPrompt_phi.root");
	TH1D *h_DsGpt_weighted=(TH1D*)f_Ds->Get("h_Gpt_weighted2");
	h_DsGpt_weighted->Scale(1/h_DsGpt_weighted->Integral());
	TH1D *h_DsGpt_weightedFineBin=(TH1D*)f_Ds->Get("h_Gpt_weighted2FineBin");
	h_DsGpt_weightedFineBin->Scale(1/h_DsGpt_weightedFineBin->Integral());

	TH1D *h_DsD0Ratio=(TH1D*)h_DsGpt_weighted->Clone("h_DsD0Ratio");
	h_DsD0Ratio->SetTitle("h_DsD0Ratio");
	h_DsD0Ratio->Divide(h_D0Gpt_weighted);


	gStyle->SetOptStat(0);
	TCanvas *c0=new TCanvas("c0","c0",800,800);
	c0->cd();
	gPad->SetLogx();
	h_DsD0Ratio->GetXaxis()->SetTitle("p_{T} GeV");
	h_DsD0Ratio->GetYaxis()->SetTitle("ratio");
	h_DsD0Ratio->GetXaxis()->SetRangeUser(0,40);
	h_DsD0Ratio->Draw();
	c0->SaveAs("DsD0Ratio.pdf");


	TCanvas *c1=new TCanvas("c1","c1",800,800);
	c1->cd();
	gPad->SetLogx();
	gPad->SetLogy();
	h_D0Gpt_weightedFineBin->SetTitle("Gpt");
	h_D0Gpt_weightedFineBin->GetXaxis()->SetTitle("p_{T} GeV");
	h_D0Gpt_weightedFineBin->SetLineColor(2);
	h_D0Gpt_weightedFineBin->Draw();
	h_DsGpt_weightedFineBin->SetLineColor(4);
	h_DsGpt_weightedFineBin->Draw("SAME");
	TLegend *le=new TLegend(0.15,0.15,0.40,0.40);
	le->SetBorderSize(0);
	le->AddEntry(h_D0Gpt_weightedFineBin,"D0 pt","l");
	le->AddEntry(h_DsGpt_weightedFineBin,"Ds pt","l");
	le->Draw("SAME");
	c1->SaveAs("h_DsD0Gpt.pdf");

}
