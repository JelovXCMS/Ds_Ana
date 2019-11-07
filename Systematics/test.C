

int test(){

	TRandom *r3 = new TRandom3();
	r3->GetSeed();
	
	TH1D *h1=new TH1D("h1","h1",10,0,10);

	
	for(int i=0; i<10;i++){
		h1->SetBinContent(i+1,5+i*5+5*r3->Gaus(0,1));
		// h1->SetBinError(i+1, h1->GetBinContent(i+1)*r3->Gaus(0,0.2));
		h1->SetBinError(i+1, sqrt(h1->GetBinContent(i+1)));

	  cout<<"i : "<<h1->GetBinContent(i+1)<<" +- "<<h1->GetBinError(i+1)<<endl;
	
	}

	//h1->Draw();

	TF1 *f1=new TF1("f1","[0]+x*[1]");
	h1->Fit("f1","EMIS0",0,10);
	h1->Fit("f1","EMIS0",0,10);
	TFitResultPtr r=h1->Fit("f1","EMIS0",0,10);;
	
	h1->Draw();
	f1->SetRange(0,10);
	f1->Draw("same");

	cout<<"f1 at 5 = "<<f1->Eval(5)<<endl;

	double x[1]={5};
	double err[1];
	r->GetConfidenceIntervals(1,1,1,x,err,0.683,false);
	cout<<"err = "<<err[0]<<endl;
	


	return 1;
}
