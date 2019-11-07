
// test for PbPb , bin 8-10 first

void compare_wide(){


	TFile *f_normal=TFile::Open("../output_sig99eff50test/FitResult_FixShape_PbPb3_pt10to20.root");
	// TFile *f_normal=TFile::Open("../output_sig99eff50test/FitResult_PbPb3_pt8to10.root");
	// TFile *f_wide=TFile::Open("../output_wide_sig99eff50test/FitResult_PbPb3_pt8to10.root");
	TFile *f_wide=TFile::Open("../output_wide_sig99eff50test/FitResult_FixShape_PbPb3_pt10to20.root");

	const int nbin=5;
	double bins[nbin+1]={0,0.1,0.2,0.3,0.4,0.5};

	TH1F *h_normal=new TH1F("h_normal","h_normal",nbin,bins); h_normal->Sumw2();
	TH1F *h_wide=new TH1F("h_wide","h_wide",nbin,bins); h_wide->Sumw2();

	for(int i=0; i<nbin; i++){
		TH1F *h_norm_temp=(TH1F*)f_normal->Get(Form("h_RawRooFitYield_Dchi2clMinScan%i", i*10+5)) ; 
		h_normal->SetBinContent(i+1,h_norm_temp->GetBinContent(1));
		h_normal->SetBinError(i+1,h_norm_temp->GetBinError(1));

		TH1F *h_wide_temp=(TH1F*)f_wide->Get(Form("h_RawRooFitYield_Dchi2clMinScan%i", i*10+5)) ; 
		h_wide->SetBinContent(i+1,h_wide_temp->GetBinContent(1));
		h_wide->SetBinError(i+1,h_wide_temp->GetBinError(1));


		delete h_norm_temp;
		delete h_wide_temp;

	}	

	gStyle->SetOptStat(0);

	TCanvas *c_yield=new TCanvas("c_yield","c_yield",800,800);
	c_yield->cd();

	h_normal->SetTitle("");
  h_normal->GetXaxis()->SetTitle("vertex probability");	
  h_normal->GetYaxis()->SetTitle("raw Yield");		
	h_normal->Draw();
	h_normal->SetMinimum(0);
	h_wide->SetLineColor(2);
	h_wide->Draw("same");

	TLegend *le=new TLegend(0.2,0.2,0.45,0.45);
	le->SetBorderSize(0);
	le->AddEntry((TObject*)0,"PbPb 10-20 GeV","");
	le->AddEntry(h_normal,"Default","l");
	le->AddEntry(h_wide,"Wide Range","l");
	le->Draw("same");


	c_yield->SaveAs("yield_compare.png");

	TH1F *h_ratio=new TH1F("h_ratio","h_ratio",nbin,bins); h_ratio->Sumw2();

	h_ratio->Divide(h_wide,h_normal,1,1,"B");

	TCanvas *c_ratio= new TCanvas("c_ratio","c_ratio",800,800);
	c_ratio->cd();
	h_ratio->SetMinimum(0.8);
	h_ratio->SetMaximum(1.2);
	h_ratio->SetTitle("");
	h_ratio->GetXaxis()->SetTitle("vertex probability");
	h_ratio->GetYaxis()->SetTitle("Wide/Default ratio");
	h_ratio->Draw();

	c_ratio->SaveAs("ratio_compare.png");







}
