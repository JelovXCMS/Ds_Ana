int compare_MCData_relPterr(){

	gStyle->SetOptStat(0);

	TFile *f_data= new TFile("f_trkrelpterr.root");
	TFile *f_MC= new TFile("MC_trkrelpterr.root");

	TH1D *h_data_trkrelptr=(TH1D*)f_data->Get("h_trksrelpterr");
	TH1D *h_MC_trkrelptr=(TH1D*)f_MC->Get("h_trksrelpterr");

	h_data_trkrelptr->Scale(1/h_data_trkrelptr->Integral());
	h_MC_trkrelptr->Scale(1/h_MC_trkrelptr->Integral());

	TCanvas *c_1=new TCanvas("c_1","c_1",800,600);
	c_1->cd();
	h_data_trkrelptr->SetTitle("");
	h_data_trkrelptr->GetXaxis()->SetTitle("rel. p_{T} err.");
	h_data_trkrelptr->Draw();
	h_MC_trkrelptr->SetLineColor(2);
	h_MC_trkrelptr->Draw("same");	
	
	TLegend *le=new TLegend(0.7,0.7,0.85,0.85);
	le->SetBorderSize(0);
	le->AddEntry(h_data_trkrelptr,"Data","l");
	le->AddEntry(h_MC_trkrelptr,"MC","l");
	le->Draw("same");	




	return 0;
}
