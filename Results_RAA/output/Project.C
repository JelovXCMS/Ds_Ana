int Project(){

	double ErrFac=sqrt(5);

	TFile *fin=TFile::Open("RAA.root","read");
	TFile *fout=TFile::Open("RAA_Project.root","recreate");
	fout->cd();	

	TH1D *h_RAA=(TH1D*)fin->Get("h_RAA");
	TH1D *h_RAA_Sta=(TH1D*)fin->Get("h_RAA_Sta");

	TH1D *h_RAA_new=(TH1D*)h_RAA->Clone("h_RAA_new");
	TH1D *h_RAA_Sta_new=(TH1D*)h_RAA_Sta->Clone("h_RAA_Sta_new");

	TH1D *h_DsoverD0_pp_Sta=(TH1D*)fin->Get("h_DsoverD0_pp_Sta");
	TH1D *h_DsoverD0_pp_Sys=(TH1D*)fin->Get("h_DsoverD0_pp_Sys");

	TH1D *h_DsoverD0_pp_Sta_new=(TH1D*)h_DsoverD0_pp_Sta->Clone("h_DsoverD0_pp_Sta_new");
	TH1D *h_DsoverD0_pp_Sys_new=(TH1D*)h_DsoverD0_pp_Sys->Clone("h_DsoverD0_pp_Sys_new");

	
	TH1D *h_DsoverD0_PbPb_Sta=(TH1D*)fin->Get("h_DsoverD0_PbPb_Sta");
	TH1D *h_DsoverD0_PbPb_Sys=(TH1D*)fin->Get("h_DsoverD0_PbPb_Sys");

	TH1D *h_DsoverD0_PbPb_Sta_new=(TH1D*)h_DsoverD0_PbPb_Sta->Clone("h_DsoverD0_PbPb_Sta_new");
	TH1D *h_DsoverD0_PbPb_Sys_new=(TH1D*)h_DsoverD0_PbPb_Sys->Clone("h_DsoverD0_PbPb_Sys_new");

	h_DsoverD0_PbPb_Sys->Draw();


	for(int i =0; i<h_RAA->GetNbinsX(); i++){

		double RAA= h_RAA->GetBinContent(i+1);
		double RAA_Sys= h_RAA->GetBinError(i+1)/ h_RAA->GetBinContent(i+1);
		double RAA_Sta= h_RAA_Sta->GetBinError(i+1)/ h_RAA->GetBinContent(i+1);
		double RAA_Sta_new=RAA_Sta/ErrFac;
		double RAA_Sys_new=sqrt( (RAA_Sys*RAA_Sys-0.18*0.18)/ErrFac+0.18*0.18);

		h_RAA_new->SetBinContent(i+1, RAA);
		h_RAA_new->SetBinError(i+1, RAA*RAA_Sys_new);
		h_RAA_Sta_new->SetBinContent(i+1, RAA);
		h_RAA_Sta_new->SetBinError(i+1, RAA*RAA_Sta_new);
		
		cout<<"h_RAA : "<<h_RAA->GetBinContent(i+1)<<" , Sys = "<<RAA_Sys<<" , RAA_Sta = "<<RAA_Sta<<" new Sys = "<<RAA_Sys_new<<" RAA_Sta_new = "<<RAA_Sta_new<<endl;

	}

	cout<<" pp "<<endl;

	for(int i =0; i<h_DsoverD0_pp_Sys->GetNbinsX(); i++){
		double DsD0= h_DsoverD0_pp_Sys->GetBinContent(i+1);
		double DsD0_Sys= h_DsoverD0_pp_Sys->GetBinError(i+1)/h_DsoverD0_pp_Sys->GetBinContent(i+1);
		double DsD0_Sta= h_DsoverD0_pp_Sta->GetBinError(i+1)/h_DsoverD0_pp_Sta->GetBinContent(i+1);
		double DsD0_Sys_new = sqrt( (DsD0_Sys*DsD0_Sys-0.04*0.04)/ErrFac + 0.04*0.04 ) ;
		double DsD0_Sta_new = DsD0_Sta/ErrFac;

		h_DsoverD0_pp_Sta_new->SetBinError(i+1,DsD0_Sta_new);
		h_DsoverD0_pp_Sys_new->SetBinError(i+1,DsD0_Sys_new);


		cout<<"DsD0 = "<<DsD0<<" , DsD0_Sys = "<<DsD0_Sys<<" , DsD0_Sta = "<<DsD0_Sta<<" , DsD0_Sys_new = "<<DsD0_Sys_new<<" ,DsD0_Sta_new = "<<DsD0_Sta_new<<endl;


	}

	cout<<" PbPb "<<endl;

	for(int i =0; i<h_DsoverD0_PbPb_Sys->GetNbinsX(); i++){
		double DsD0= h_DsoverD0_PbPb_Sys->GetBinContent(i+1);
		double DsD0_Sys= h_DsoverD0_PbPb_Sys->GetBinError(i+1)/h_DsoverD0_PbPb_Sys->GetBinContent(i+1);
		double DsD0_Sta= h_DsoverD0_PbPb_Sta->GetBinError(i+1)/h_DsoverD0_PbPb_Sta->GetBinContent(i+1);
		double DsD0_Sys_new = sqrt( (DsD0_Sys*DsD0_Sys-0.05*0.05)/ErrFac + 0.05*0.05 ) ;
		double DsD0_Sta_new = DsD0_Sta/ErrFac;

		cout<<"DsD0 = "<<DsD0<<" , DsD0_Sys = "<<DsD0_Sys<<" , DsD0_Sta = "<<DsD0_Sta<<" , DsD0_Sys_new = "<<DsD0_Sys_new<<" ,DsD0_Sta_new = "<<DsD0_Sta_new<<endl;

		h_DsoverD0_PbPb_Sta_new->SetBinError(i+1,DsD0_Sta_new);
		h_DsoverD0_PbPb_Sys_new->SetBinError(i+1,DsD0_Sys_new);


	}

	h_RAA->Write();
	h_RAA_Sta->Write();
	h_RAA_new->Write();
	h_RAA_Sta_new->Write();

	h_DsoverD0_pp_Sta->Write();
	h_DsoverD0_pp_Sys->Write();
	h_DsoverD0_pp_Sta_new->Write();
	h_DsoverD0_pp_Sys_new->Write();

	h_DsoverD0_PbPb_Sta->Write();
	h_DsoverD0_PbPb_Sys->Write();
	h_DsoverD0_PbPb_Sta_new->Write();
	h_DsoverD0_PbPb_Sys_new->Write();

	return 0;
}
