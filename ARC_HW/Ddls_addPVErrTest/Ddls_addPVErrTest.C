#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"

int Ddls_addPVErrTest(){

	gSystem->Exec("mkdir -p plots");
	gStyle->SetOptStat(0);
	TFile *fin=TFile::Open("/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_offMC_test190202/DsMinTree_pp_offMC_Prompt_Phi_pt6.root","READ");
	TTree *ntDs=(TTree*)fin->Get("ntDs");

	TH1D *hDdls=(TH1D*)fin->Get("h_Ddls");
	TH1D *hDdls_Gaus=(TH1D*)fin->Get("h_Ddls_Gaus");

	hDdls->GetXaxis()->SetRangeUser(2,100);
	// hDdls->Draw();
	hDdls_Gaus->GetXaxis()->SetRangeUser(2,100);
	// hDdls_Gaus->Draw("same");

	TH1D *hDdls_Cum=(TH1D*)hDdls->Clone("hDdls_Cum");
	TH1D *hDdls_Gaus_Cum=(TH1D*)hDdls_Gaus->Clone("hDdls_Gaus_Cum");

	int nbin_h=hDdls_Cum->GetNbinsX();
	cout<<"nbin_h= "<<nbin_h<<endl;
	
	double Cum_Ddls=0;
	double Cum_Ddls_Gaus=0;

	for(int i =nbin_h; i>=0; i--){
		Cum_Ddls+=hDdls->GetBinContent(i+1);
		hDdls_Cum->SetBinContent(i+1,Cum_Ddls);
		hDdls_Cum->SetBinError(i+1,sqrt(Cum_Ddls));

		Cum_Ddls_Gaus+=hDdls_Gaus->GetBinContent(i+1);
		hDdls_Gaus_Cum->SetBinContent(i+1,Cum_Ddls_Gaus);
		hDdls_Gaus_Cum->SetBinError(i+1,sqrt(Cum_Ddls_Gaus));
	}



	hDdls_Cum->Draw();
	hDdls_Gaus_Cum->SetLineColor(2);
	hDdls_Gaus_Cum->Draw("Same");




	TH1D *h_Gaus_EffRatio=(TH1D*)hDdls_Gaus_Cum->Clone("h_Gaus_EffRatio");
	h_Gaus_EffRatio->Divide(hDdls_Cum);


	TLatex *tla=new TLatex();

	TCanvas *c_Ddls_compare = new TCanvas("c_Ddls_compare","c_Ddls_compare");
	c_Ddls_compare->cd();
	hDdls->SetTitle("");
	hDdls->GetXaxis()->SetTitle("Decay Length Significance");
	hDdls->GetXaxis()->SetRangeUser(2,10);
	hDdls->Draw();
	hDdls_Gaus->SetLineColor(2);
	hDdls_Gaus->Draw("same");

	TLegend *tle =new TLegend(0.55,0.55,0.85,0.85);
	tle->SetBorderSize(0);
	tle->AddEntry((TObject*)0,"pp : 8<p_{T}<10 GeV/c","");
	tle->AddEntry(hDdls,"Default","pl");
	tle->AddEntry(hDdls_Gaus,"Gaus. smear with PVErr.","pl");
	tle->Draw("same");

	c_Ddls_compare->SaveAs("plots/Ddls_GausCompare.png");


	TCanvas *c_Ddls_GausEffratio = new TCanvas("c_Ddls_GausEffratio","c_Ddls_GausEffratio");
	c_Ddls_GausEffratio->cd();

	h_Gaus_EffRatio->GetXaxis()->SetRangeUser(2,10);
	h_Gaus_EffRatio->GetXaxis()->SetTitle("Decay Length Significance");
	h_Gaus_EffRatio->GetYaxis()->SetTitle("Eff. Ratio : smeared/default");
	h_Gaus_EffRatio->SetTitle("");
	h_Gaus_EffRatio->SetMinimum(0.7);
	h_Gaus_EffRatio->SetMaximum(1.3);
	// h_Gaus_EffRatio->GetXaxis()->Rebin(2);
	 h_Gaus_EffRatio->Draw();

	tla->DrawLatexNDC(textposx,textposy,"pp : 8<pt<10 GeV/c");

	c_Ddls_GausEffratio->SaveAs("plots/Ddls_GausEffRatio.png");



	return 1;


///////////////////////

		
	int nbin=92;
	double binLow=1.6;
	double binHigh=100;

	double firstBinWidth=0.2;
	double BinWRatio=1.1;
	double cbin=binLow;
	double cbinw=firstBinWidth;
	int nbinc=0;

/*
	while(cbin<binHigh){
		nbinc++;
		cbin=cbin+cbinw;
		cbinw=cbinw*BinWRatio;	
		cout<<cbin<<" , ";
		// cout<<"cbin = "<<cbin<<endl;
	}
*/

	// double bins[45]={2.05 , 2.104 , 2.16232 , 2.22531 , 2.29333 , 2.3668 , 2.44614 , 2.53183 , 2.62438 , 2.72433 , 2.83227 , 2.94886 , 3.07476 , 3.21075 , 3.35761 , 3.51621 , 3.68751 , 3.87251 , 4.07231 , 4.2881 , 4.52115 , 4.77284 , 5.04466 , 5.33824 , 5.6553 , 5.99772 , 6.36754 , 6.76694 , 7.1983 , 7.66416 , 8.16729 , 8.71068 , 9.29753 , 9.93133 , 10.6158 , 11.3551 , 12.1535 , 13.0158 , 13.9471 , 14.9528 , 16.0391 , 17.2122 , 18.4792 , 19.8475 , 21.3253};
	// int nbin1=44;

	// double bins[]={2.08 , 2.1656 , 2.25719 , 2.3552 , 2.46006 , 2.57226 , 2.69232 , 2.82078 , 2.95824 , 3.10532 , 3.26269 , 3.43108 , 3.61125 , 3.80404 , 4.01032 , 4.23104 , 4.46722 , 4.71992 , 4.99032 , 5.27964 , 5.58921 , 5.92046 , 6.27489 , 6.65413 , 7.05992 , 7.49412 , 7.95871 , 8.45582 , 8.98772 , 9.55686 , 10.1658 , 10.8175 , 11.5147 , 12.2607 , 13.059 , 13.9131 , 14.827 , 15.8049 , 16.8512 , 17.9708 , 19.1688 , 20.4506};
	// int nbin1=41;

	// int nbin1=25;
	// double bins[]={2 ,2.2 , 2.42 , 2.662 , 2.9282 , 3.22102 , 3.54312 , 3.89743 , 4.28718 , 4.7159 , 5.18748 , 5.70623 , 6.27686 , 6.90454 , 7.595 , 8.3545 , 9.18995 , 10.1089 , 11.1198 , 12.2318 , 13.455 , 14.8005 , 16.2805 , 17.9086 , 19.6995 , 21.6694};


	cout<<"nbinc = "<<nbinc<<endl;

	double DptLow=8;
	double DptHigh=10;

	double DvtxP=0.05;
	

	TH1D *h_Ddls=new TH1D("h_Ddls","; Decay Length Significance",nbin,binLow,binHigh);
	TH1D *h_DdlsPPVE=new TH1D("h_DdlsPPVE","; Decay Length Significance",nbin,binLow,binHigh);
	TH1D *h_DdlsMPVE=new TH1D("h_DdlsMPVE","; Decay Length Significance",nbin,binLow,binHigh);

	// TH1D *h_Ddls=new TH1D("h_Ddls","; Decay Length Significance",nbin1,bins);
	// TH1D *h_DdlsPPVE=new TH1D("h_DdlsPPVE","; Decay Length Significance",nbin1,bins);
	// TH1D *h_DdlsMPVE=new TH1D("h_DdlsMPVE","; Decay Length Significance",nbin1,bins);


	ntDs->Project("h_Ddls","Ddls",(TCut)Form("Dpt>%f && Dpt<%f && Dchi2cl>%f && Dalpha<0.12 && DsGen==23333 && DgencollisionId==0 && DgenBAncestorpt<=0", DptLow,DptHigh, DvtxP));
	ntDs->Project("h_DdlsPPVE","DdlsPPVE",(TCut)Form("Dpt>%f && Dpt<%f && Dchi2cl>%f && Dalpha<0.12 && DsGen==23333 && DgencollisionId==0 && DgenBAncestorpt<=0", DptLow,DptHigh, DvtxP));
	ntDs->Project("h_DdlsMPVE","DdlsMPVE",(TCut)Form("Dpt>%f && Dpt<%f && Dchi2cl>%f && Dalpha<0.12 && DsGen==23333 && DgencollisionId==0 && DgenBAncestorpt<=0", DptLow,DptHigh, DvtxP));
	
	h_Ddls->Draw();
	h_DdlsPPVE->SetLineColor(1);
	h_DdlsPPVE->Draw("same");
	h_DdlsMPVE->SetLineColor(2);
	h_DdlsMPVE->Draw("same");


	TH1D *h_Ddls_GausS=new TH1D("h_Ddls_GausS",";Decay Length Significance",nbin,binLow,binHigh);
	
	for(int i =0; i<nbin; i++){
		h_Ddls_GausS->Fill(h_Ddls->GetBinCenter(i+1), h_Ddls->GetBinContent(i+1) + gRandom->Gaus(0,0.5));
	}



	
	TH1D *h_DdlsPPVE_Ratio=(TH1D*)h_DdlsPPVE->Clone("h_DdlsPPVE_Ratio");
	h_DdlsPPVE_Ratio->Divide(h_Ddls);

	TH1D *h_Ddls_Cum=new TH1D("h_Ddls_Cum","; Decay Length Significance",nbin,binLow,binHigh);
	double Ddls_Cum=0;

	TH1D *h_DdlsPPVE_Cum=new TH1D("h_DdlsPPVE_Cum","; Decay Length Significance",nbin,binLow,binHigh);
	double DdlsPPVE_Cum=0;

	for(int i=nbin; i>=0; i--){
		Ddls_Cum+=h_Ddls->GetBinContent(i+1);
		h_Ddls_Cum->SetBinContent(i+1,Ddls_Cum);
		h_Ddls_Cum->SetBinError(i+1,sqrt(Ddls_Cum));

		DdlsPPVE_Cum+=h_DdlsPPVE->GetBinContent(i+1);
		h_DdlsPPVE_Cum->SetBinContent(i+1,DdlsPPVE_Cum);
		h_DdlsPPVE_Cum->SetBinError(i+1,sqrt(DdlsPPVE_Cum));

	}

	h_Ddls_Cum->SetLineColor(1);
	h_Ddls_Cum->Draw();
	h_DdlsPPVE_Cum->SetLineColor(2);
	h_DdlsPPVE_Cum->Draw("same");


  TH1D *h_DdlsPPVE_CumRatio=(TH1D*)h_DdlsPPVE_Cum->Clone("h_DdlsPPVE_CumRatio");
  h_DdlsPPVE_CumRatio->Divide(h_Ddls_Cum);

	h_DdlsPPVE_CumRatio->Draw();

	// h_DdlsPPVE_Ratio->Draw();

	h_Ddls->Draw();
	// h_Ddls_GausS->Draw();


	return 0;
}
