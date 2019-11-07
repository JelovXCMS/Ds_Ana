


int test_th2(){



	TFile *fout=new TFile("test.root","recreate");
	fout->cd();
	TH2D *h2= new TH2D("h2","h2",100,0,100,100,0,100);
	h2->Sumw2();

	for(int i=0; i<80; i++){

		h2->Fill(i,100-i);

	}

	double corr=h2->GetCorrelationFactor();
	cout<<"corr = "<<corr<<endl;

	h2->Write();



	return 0;

}
