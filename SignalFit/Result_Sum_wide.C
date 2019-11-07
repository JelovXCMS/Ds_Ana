// plot ratio of deafult fit & wide fit
// plot Ds/Dp double ratio
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"


int Result_Sum_wide(){

	gSystem->Exec(Form("mkdir -p output_wide%s/plots",s_CutSet.Data()));
	gStyle->SetOptStat(0);

	TLatex *tlatex=new TLatex();
	double shiftY=0;
	double shiftX=0;

	double DptLow=0;
	double DptHigh=0;


  TString inName_pp_default=Form("./output%s/FitResult_FixShapeTrkPtScan_pp",s_CutSet.Data());
  TString inName_pp_wide=Form("./output_wide%s/FitResult_pp",s_CutSet.Data());
//  TString str_PbPb="pp";
	TFile *f_pp_default[nbin_pt_pp];
	TFile *f_pp_wide[nbin_pt_pp];

	TH1F *h_RawYield_pp_default_ptbin[nbin_pt_pp];
	TH1F *h_RawYield_pp_default=new TH1F("h_RawYield_pp_default","h_RawYield_pp_default",nbin_pt_pp,bins_pt_pp);

	TH1F *h_RawYield_pp_wide_ptbin[nbin_pt_pp];
	TH1F *h_RawYield_pp_wide=new TH1F("h_RawYield_pp_wide","h_RawYield_pp_wide",nbin_pt_pp,bins_pt_pp);

	TH1F *h_RawYield_pp_wideDp_ptbin[nbin_pt_pp];
	TH1F *h_RawYield_pp_wideDp=new TH1F("h_RawYield_pp_wideDp","h_RawYield_pp_wideDp",nbin_pt_pp,bins_pt_pp);

	TH1F *h_RawYield_pp_DsOverDp=new TH1F("h_RawYield_pp_DsOverDp","h_RawYield_pp_DsOverDp",nbin_pt_pp,bins_pt_pp);
	// TH1F *h_RawYield_pp_WideOverDefault=new TH1F("h_RawYield_pp_WideOverDefault","h_RawYield_pp_WideOverDefault",nbin_pt_pp,bins_pt_pp);

	for(int i =0; i<nbin_pt_pp; i++){

		DptLow=bins_pt_pp[i];
		DptHigh=bins_pt_pp[i+1];	

		f_pp_default[i]=TFile::Open(Form("%s_pt%.0fto%.0f.root",inName_pp_default.Data(),DptLow,DptHigh));
		// cout<<"fin name = "<<Form("%s_pt%.0fto%.0f.root",inName_pp_default.Data(),DptLow,DptHigh)<<endl;

		h_RawYield_pp_default_ptbin[i]=(TH1F*)f_pp_default[i]->Get("h_RawRooFitYield");	
		h_RawYield_pp_default->SetBinContent(i+1,h_RawYield_pp_default_ptbin[i]->GetBinContent(1));
		h_RawYield_pp_default->SetBinError(i+1,h_RawYield_pp_default_ptbin[i]->GetBinError(1));

		// cout<<"bin i : "<<i<<" , Yield = "<<h_RawYield_pp_default_ptbin[i]->GetBinContent(1)<<endl;

		f_pp_wide[i]=TFile::Open(Form("%s_pt%.0fto%.0f.root",inName_pp_wide.Data(),DptLow,DptHigh));
		// cout<<"fin name = "<<Form("%s_pt%.0fto%.0f.root",inName_pp_wide.Data(),DptLow,DptHigh)<<endl;

		h_RawYield_pp_wide_ptbin[i]=(TH1F*)f_pp_wide[i]->Get("h_RawRooFitYield");	
		h_RawYield_pp_wide->SetBinContent(i+1,h_RawYield_pp_wide_ptbin[i]->GetBinContent(1));
		h_RawYield_pp_wide->SetBinError(i+1,h_RawYield_pp_wide_ptbin[i]->GetBinError(1));

		// cout<<"bin i : "<<i<<" , Yield = "<<h_RawYield_pp_wide_ptbin[i]->GetBinContent(1)<<endl;

		h_RawYield_pp_wideDp_ptbin[i]=(TH1F*)f_pp_wide[i]->Get("h_RawRooFitYieldDp");	
		h_RawYield_pp_wideDp->SetBinContent(i+1,h_RawYield_pp_wideDp_ptbin[i]->GetBinContent(1));
		h_RawYield_pp_wideDp->SetBinError(i+1,h_RawYield_pp_wideDp_ptbin[i]->GetBinError(1));

		cout<<"bin i : "<<i<<" ,Dp Yield = "<<h_RawYield_pp_wideDp_ptbin[i]->GetBinContent(1)<<" ,relErr = "<<h_RawYield_pp_wideDp_ptbin[i]->GetBinError(1) / h_RawYield_pp_wideDp_ptbin[i]->GetBinContent(1) <<endl;
		cout<<"bin i : "<<i<<" ,Ds Yield = "<<h_RawYield_pp_wide_ptbin[i]->GetBinContent(1)<<" ,relErr = "<<h_RawYield_pp_wide_ptbin[i]->GetBinError(1) / h_RawYield_pp_wide_ptbin[i]->GetBinContent(1) <<endl;

		double DsOverDp=h_RawYield_pp_wide->GetBinContent(i+1)/h_RawYield_pp_wideDp->GetBinContent(i+1);
		double DsOverDpErr=DsOverDp*sqrt(pow(h_RawYield_pp_wide->GetBinError(i+1)/h_RawYield_pp_wide->GetBinContent(i+1) ,2) + pow(h_RawYield_pp_wideDp->GetBinError(i+1)/h_RawYield_pp_wideDp->GetBinContent(i+1) ,2)  );

		h_RawYield_pp_DsOverDp->SetBinContent(i+1,DsOverDp);
		h_RawYield_pp_DsOverDp->SetBinError(i+1,DsOverDpErr);
		

//		h_RawYield_pp_default_ptbin[i]->Draw();

	//	return 1;
	}	

	TCanvas *c_pp_DsOverDp = new TCanvas("c_pp_DsOverDp","c_pp_DsOverDp");
	c_pp_DsOverDp->cd();
	
  h_RawYield_pp_DsOverDp->GetXaxis()->SetRangeUser(2,40);
	h_RawYield_pp_DsOverDp->SetTitle("");
	h_RawYield_pp_DsOverDp->GetXaxis()->SetTitle("M_{KK#pi} p_{T} (GeV/c)");
	h_RawYield_pp_DsOverDp->GetYaxis()->SetTitle("Yield Ratio D_{S}/D^{+}");
	h_RawYield_pp_DsOverDp->SetMarkerStyle(21);
	h_RawYield_pp_DsOverDp->Draw();
	
	tlatex->DrawLatexNDC(textposx,textposy,"pp D_{S}/D^{+} Yield Ratio");
	c_pp_DsOverDp->SaveAs(Form("./output_wide%s/plots/pp_DsOverDp.png",s_CutSet.Data()));


	TCanvas *c_pp_WideOverDefault = new TCanvas("c_pp_WideOverDefault","c_pp_WideOverDefault");
	c_pp_WideOverDefault->cd();	

	TH1F *h_RawYield_pp_WideOverDefault=(TH1F*)h_RawYield_pp_wide->Clone("h_RawYield_pp_WideOverDefault");
	h_RawYield_pp_WideOverDefault->Divide(h_RawYield_pp_default);
	h_RawYield_pp_WideOverDefault->SetTitle("");
	h_RawYield_pp_WideOverDefault->GetXaxis()->SetTitle("M_{KK#pi} p_{T} (GeV/c)");
	h_RawYield_pp_WideOverDefault->GetYaxis()->SetTitle("Yield Ratio Wide/Default");

	h_RawYield_pp_WideOverDefault->GetXaxis()->SetRangeUser(2,40);
	h_RawYield_pp_WideOverDefault->SetMarkerStyle(21);
	h_RawYield_pp_WideOverDefault->Draw();
	h_RawYield_pp_WideOverDefault->SetMaximum(1.4);
	h_RawYield_pp_WideOverDefault->SetMinimum(0.6);

	tlatex->DrawLatexNDC(textposx,textposy,"pp Wide/Default Yield Ratio");

	c_pp_WideOverDefault->SaveAs(Form("./output_wide%s/plots/pp_WideOverDefault.png",s_CutSet.Data()));

	// h_RawYield_pp_wideDp->Draw();


  TString inName_PbPb3_default=Form("./output%s/FitResult_FixShapeTrkPtScan_PbPb3",s_CutSet.Data());
  TString inName_PbPb3_wide=Form("./output_wide%s/FitResult_PbPb3",s_CutSet.Data());
//  TString str_PbPb="PbPb3";
	TFile *f_PbPb3_default[nbin_pt_PbPb3];
	TFile *f_PbPb3_wide[nbin_pt_PbPb3];

	TH1F *h_RawYield_PbPb3_default_ptbin[nbin_pt_PbPb3];
	TH1F *h_RawYield_PbPb3_default=new TH1F("h_RawYield_PbPb3_default","h_RawYield_PbPb3_default",nbin_pt_PbPb3,bins_pt_PbPb3);

	TH1F *h_RawYield_PbPb3_wide_ptbin[nbin_pt_PbPb3];
	TH1F *h_RawYield_PbPb3_wide=new TH1F("h_RawYield_PbPb3_wide","h_RawYield_PbPb3_wide",nbin_pt_PbPb3,bins_pt_PbPb3);

	TH1F *h_RawYield_PbPb3_wideDp_ptbin[nbin_pt_PbPb3];
	TH1F *h_RawYield_PbPb3_wideDp=new TH1F("h_RawYield_PbPb3_wideDp","h_RawYield_PbPb3_wideDp",nbin_pt_PbPb3,bins_pt_PbPb3);

	TH1F *h_RawYield_PbPb3_DsOverDp=new TH1F("h_RawYield_PbPb3_DsOverDp","h_RawYield_PbPb3_DsOverDp",nbin_pt_PbPb3,bins_pt_PbPb3);
	// TH1F *h_RawYield_PbPb3_WideOverDefault=new TH1F("h_RawYield_PbPb3_WideOverDefault","h_RawYield_PbPb3_WideOverDefault",nbin_pt_PbPb3,bins_pt_PbPb3);

	for(int i =2; i<nbin_pt_PbPb3; i++){

		DptLow=bins_pt_PbPb3[i];
		DptHigh=bins_pt_PbPb3[i+1];	

		f_PbPb3_default[i]=TFile::Open(Form("%s_pt%.0fto%.0f.root",inName_PbPb3_default.Data(),DptLow,DptHigh));
		// cout<<"fin name = "<<Form("%s_pt%.0fto%.0f.root",inName_PbPb3_default.Data(),DptLow,DptHigh)<<endl;

		h_RawYield_PbPb3_default_ptbin[i]=(TH1F*)f_PbPb3_default[i]->Get("h_RawRooFitYield");	
		h_RawYield_PbPb3_default->SetBinContent(i+1,h_RawYield_PbPb3_default_ptbin[i]->GetBinContent(1));
		h_RawYield_PbPb3_default->SetBinError(i+1,h_RawYield_PbPb3_default_ptbin[i]->GetBinError(1));

		// cout<<"bin i : "<<i<<" , Yield = "<<h_RawYield_PbPb3_default_ptbin[i]->GetBinContent(1)<<endl;

		f_PbPb3_wide[i]=TFile::Open(Form("%s_pt%.0fto%.0f.root",inName_PbPb3_wide.Data(),DptLow,DptHigh));
		// cout<<"fin name = "<<Form("%s_pt%.0fto%.0f.root",inName_PbPb3_wide.Data(),DptLow,DptHigh)<<endl;

		h_RawYield_PbPb3_wide_ptbin[i]=(TH1F*)f_PbPb3_wide[i]->Get("h_RawRooFitYield");	
		h_RawYield_PbPb3_wide->SetBinContent(i+1,h_RawYield_PbPb3_wide_ptbin[i]->GetBinContent(1));
		h_RawYield_PbPb3_wide->SetBinError(i+1,h_RawYield_PbPb3_wide_ptbin[i]->GetBinError(1));

		// cout<<"bin i : "<<i<<" , Yield = "<<h_RawYield_PbPb3_wide_ptbin[i]->GetBinContent(1)<<endl;

		h_RawYield_PbPb3_wideDp_ptbin[i]=(TH1F*)f_PbPb3_wide[i]->Get("h_RawRooFitYieldDp");	
		h_RawYield_PbPb3_wideDp->SetBinContent(i+1,h_RawYield_PbPb3_wideDp_ptbin[i]->GetBinContent(1));
		h_RawYield_PbPb3_wideDp->SetBinError(i+1,h_RawYield_PbPb3_wideDp_ptbin[i]->GetBinError(1));

		cout<<"bin i : "<<i<<" ,Dp Yield = "<<h_RawYield_PbPb3_wideDp_ptbin[i]->GetBinContent(1)<<" ,relErr = "<<h_RawYield_PbPb3_wideDp_ptbin[i]->GetBinError(1) / h_RawYield_PbPb3_wideDp_ptbin[i]->GetBinContent(1) <<endl;
		cout<<"bin i : "<<i<<" ,Ds Yield = "<<h_RawYield_PbPb3_wide_ptbin[i]->GetBinContent(1)<<" ,relErr = "<<h_RawYield_PbPb3_wide_ptbin[i]->GetBinError(1) / h_RawYield_PbPb3_wide_ptbin[i]->GetBinContent(1) <<endl;

		double DsOverDp=h_RawYield_PbPb3_wide->GetBinContent(i+1)/h_RawYield_PbPb3_wideDp->GetBinContent(i+1);
		double DsOverDpErr=DsOverDp*sqrt(pow(h_RawYield_PbPb3_wide->GetBinError(i+1)/h_RawYield_PbPb3_wide->GetBinContent(i+1) ,2) + pow(h_RawYield_PbPb3_wideDp->GetBinError(i+1)/h_RawYield_PbPb3_wideDp->GetBinContent(i+1) ,2)  );

		h_RawYield_PbPb3_DsOverDp->SetBinContent(i+1,DsOverDp);
		h_RawYield_PbPb3_DsOverDp->SetBinError(i+1,DsOverDpErr);
	

		cout<<"DsOverDp = "<<DsOverDp<<" +- "<<DsOverDpErr<<endl;	

//		h_RawYield_PbPb3_default_ptbin[i]->Draw();

	//	return 1;
	}	

	TCanvas *c_PbPb3_DsOverDp = new TCanvas("c_PbPb3_DsOverDp","c_PbPb3_DsOverDp");
	c_PbPb3_DsOverDp->cd();
	
  h_RawYield_PbPb3_DsOverDp->GetXaxis()->SetRangeUser(6,40);
	h_RawYield_PbPb3_DsOverDp->SetTitle("");
	h_RawYield_PbPb3_DsOverDp->GetXaxis()->SetTitle("M_{KK#pi} p_{T} (GeV/c)");
	h_RawYield_PbPb3_DsOverDp->GetYaxis()->SetTitle("Yield Ratio D_{S}/D^{+}");
	h_RawYield_PbPb3_DsOverDp->SetMarkerStyle(21);
	h_RawYield_PbPb3_DsOverDp->Draw();
	
	tlatex->DrawLatexNDC(textposx,textposy,"PbPb D_{S}/D^{+} Yield Ratio");
	c_PbPb3_DsOverDp->SaveAs(Form("./output_wide%s/plots/PbPb3_DsOverDp.png",s_CutSet.Data()));

	TCanvas *c_PbPb3_WideOverDefault = new TCanvas("c_PbPb3_WideOverDefault","c_PbPb3_WideOverDefault");
	c_PbPb3_WideOverDefault->cd();	

	TH1F *h_RawYield_PbPb3_WideOverDefault=(TH1F*)h_RawYield_PbPb3_wide->Clone("h_RawYield_PbPb3_WideOverDefault");
	h_RawYield_PbPb3_WideOverDefault->Divide(h_RawYield_PbPb3_default);
	h_RawYield_PbPb3_WideOverDefault->SetTitle("");
	h_RawYield_PbPb3_WideOverDefault->GetXaxis()->SetTitle("M_{KK#pi} p_{T} (GeV/c)");
	h_RawYield_PbPb3_WideOverDefault->GetYaxis()->SetTitle("Yield Ratio Wide/Default");

	h_RawYield_PbPb3_WideOverDefault->GetXaxis()->SetRangeUser(6,40);
	h_RawYield_PbPb3_WideOverDefault->SetMarkerStyle(21);
	h_RawYield_PbPb3_WideOverDefault->Draw();
	h_RawYield_PbPb3_WideOverDefault->SetMaximum(1.4);
	h_RawYield_PbPb3_WideOverDefault->SetMinimum(0.6);

	tlatex->DrawLatexNDC(textposx,textposy,"PbPb Wide/Default Yield Ratio");

	c_PbPb3_WideOverDefault->SaveAs(Form("./output_wide%s/plots/PbPb3_WideOverDefault.png",s_CutSet.Data()));


	// DsOverDp pp/PbPb double ratio

	TH1D *h_RawYield_PbPbpp_DsOverDp=new TH1D("h_RawYield_PbPbpp_DsOverDp",";M_{KK#pi} p_{T} (GeV/c);D_{S}/D^{+} PbPb/pp",nbin_pt_PbPb3,bins_pt_PbPb3);

	for(int i=2; i<nbin_pt_PbPb3; i++){

		double ratio=h_RawYield_PbPb3_DsOverDp->GetBinContent(i+1)/ h_RawYield_pp_DsOverDp->GetBinContent(i+1+2);
		cout<<"ratio = "<<ratio<<endl;
		double ppRelErr=h_RawYield_pp_DsOverDp->GetBinError(i+1+2)/ h_RawYield_pp_DsOverDp->GetBinContent(i+1+2);
		double PbPbRelErr=h_RawYield_PbPb3_DsOverDp->GetBinError(i+1)/ h_RawYield_PbPb3_DsOverDp->GetBinContent(i+1);
		double ratioErr = ratio*sqrt(pow(ppRelErr ,2) +pow(PbPbRelErr,2));
	
		h_RawYield_PbPbpp_DsOverDp->SetBinContent(i+1,ratio);
		h_RawYield_PbPbpp_DsOverDp->SetBinError(i+1,ratioErr);

	}


	TCanvas *c_PbPbpp_DsOverDp=new TCanvas("PbPbpp_DsOverDp","PbPbpp_DsOverDp");
	c_PbPbpp_DsOverDp->cd();

	h_RawYield_PbPbpp_DsOverDp->GetXaxis()->SetRangeUser(6,40);
	h_RawYield_PbPbpp_DsOverDp->SetMarkerStyle(21);
	h_RawYield_PbPbpp_DsOverDp->Draw();

	tlatex->DrawLatexNDC(textposx,textposy,"RawYield D_{S}/D^{+} PbPb/pp Double Ratio");

	c_PbPbpp_DsOverDp->SaveAs(Form("./output_wide%s/plots/PbPbpp_DsOverDp.png",s_CutSet.Data()));




	return 0;

}



