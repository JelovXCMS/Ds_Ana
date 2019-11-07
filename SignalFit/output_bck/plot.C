void plot(){

	TFile *fin=TFile::Open("FitResult_pp.root");
	TH1D *h1=(TH1D*)fin->Get("h_DsMassData_pt8to10");
	TF1 *f1=(TF1*)fin->Get("myfunc");
	
	TCanvas *c_myfunc = new TCanvas("c_myfunc","c_myfunc",800,800);
	c_myfunc->cd();
	h1->Draw();
	// f1->Draw("same");
	// cout<<"f1->GetParameter(0) = "<<f1->GetParameter(0)<<endl;
	TF1 *f3=(TF1*)h1->GetFunction("f1_DsMix");
	TCanvas *c_test2 =new TCanvas("c_test2","c_test2",800,800);
	c_test2->cd();
	h1->Draw();
	f3->Draw("same");

	TF1 *f2=(TF1*)fin->Get("f1_DsMix");
	TCanvas *c_f1DsMix = new TCanvas("c_f1DsMix","c_f1DsMix",800,800);
	c_f1DsMix->cd();
	h1->Draw();
	f2->SetRange(1.91,2.11);
	f2->Draw("same"); // this not work...
	cout<<"f2->GetParameter(0) = "<<f2->GetParameter(0)<<endl;
	return;
}
