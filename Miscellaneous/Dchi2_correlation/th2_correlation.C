int th2_correlation(TString Svar="Dtrk1PtErr/Dtrk1Pt", int nbin=100, double binlow =0, double binhi=0.1){

	TFile *fin= new TFile("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/offpp_MC/DsNtuple_pp_Prompt_Phi_pt6.root");

	TTree *ntD=(TTree*)fin->Get("ntDPhikkpi");
	TTree *ntHi=(TTree*)fin->Get("ntHi");

	ntD->AddFriend(ntHi);

	TH2D *h2=new TH2D("h2","h2",100,0,1,nbin,binlow,binhi); h2->Sumw2();

	ntD->Draw(Form("%s:Dchi2cl>>h2",Svar.Data()),"weight");
	cout<<Svar<<" corr =  " <<h2->GetCorrelationFactor()<<endl;

	TCanvas *c_test=new TCanvas("c_test","c_test",600,600);
	c_test->cd();
	h2->Draw();
	c_test->SaveAs(Form("%s.png",Svar.Data()));



	return 0;

}
