#include "../include/uti.h"
#include "../include/DsMinTreeLoad.h"
#include "../include/parameters.h"


#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting_simple.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"


#include <TString.h>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>
#include <TLatex.h>


void Plot_Eff(Int_t isPbPb=0, Int_t DsChannel=0, int prelim=0){

	InitStyle();
  double tex_upperY=0.93;

	int nbin_pt=nbin_pt_pp;
	if(isPbPb>0){
		nbin_pt=nbin_pt_PbPb3;
	}

	TString s_pre="";
	if(prelim){
		s_pre="_pre";
	}

  TLatex* texCmsPre = new TLatex(0.15,tex_upperY, "#scale[1.25]{#bf{CMS}} #it{Preliminary}");
  texCmsPre->SetNDC();
  texCmsPre->SetTextAlign(12);
  texCmsPre->SetTextSize(0.035);
  texCmsPre->SetTextFont(42);

  TLatex* texCmsSim = new TLatex(0.15,tex_upperY, "#scale[1.25]{#bf{CMS}} #it{Simulation}");
	if(prelim){
		texCmsSim = new TLatex(0.15,tex_upperY, "#scale[1.25]{#bf{CMS}} #it{Simulation Preliminary} ");
	}
  texCmsSim->SetNDC();
  texCmsSim->SetTextAlign(12);
  texCmsSim->SetTextSize(0.035);
  texCmsSim->SetTextFont(42);


  TLatex* texColPbPb = new TLatex(0.90,tex_upperY, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  texColPbPb->SetNDC();
  texColPbPb->SetTextAlign(32);
  texColPbPb->SetTextSize(0.035);
  texColPbPb->SetTextFont(42);

  TLatex* texColpp = new TLatex(0.90,tex_upperY, "pp #sqrt{s} = 5.02 TeV");
  texColpp->SetNDC();
  texColpp->SetTextAlign(32);
  texColpp->SetTextSize(0.035);
  texColpp->SetTextFont(42);

	TString str_PbPb="pp";
	if(isPbPb==3){
		str_PbPb="PbPb3";
	}
	TString str_DsChannel="phikkpi";
	// TString str_DsChannel_plot="D_{S}^{#pm} #rightarrow #phi #pi^{#pm}";
	TString str_DsChannel_plot="D_{S}^{+} #rightarrow #phi#pi^{+} #rightarrow K^{+}K^{-}#pi^{+}";
	if(DsChannel==1){
		str_DsChannel="f0kkpi";
		str_DsChannel_plot="D_{S}^{#pm} #rightarrow f0 #pi^{#pm}";
	}

	TString fname_NP=Form("./output%s/MC_eff_%s_NonPrompt_%s.root",s_CutSet.Data(), str_PbPb.Data(),str_DsChannel.Data());
	TString fname_P=Form("./output%s/MC_eff_%s_Prompt_%s.root",s_CutSet.Data(), str_PbPb.Data(),str_DsChannel.Data());

	TFile *f_P=TFile::Open(fname_P,"READ");
	TFile *f_NP=TFile::Open(fname_NP,"READ");

	TH1D *h_GenACCEff_P=(TH1D*)f_P->Get("h_GenAccEff");
	TH1D *h_GenACCEff_NP=(TH1D*)f_NP->Get("h_GenAccEff");

	TH1D *h_RecoNormEff_P=(TH1D*)f_P->Get("h_RecoNormEff");
	TH1D *h_RecoNormEff_NP=(TH1D*)f_NP->Get("h_RecoNormEff");

	TH1D *h_RecoNormEff_FONLL_P=(TH1D*)f_P->Get("h_RecoNormEff_FONLL");
	TH1D *h_RecoNormEff_FONLL_NP=(TH1D*)f_NP->Get("h_RecoNormEff_FONLL");


	TH1D *h_EffRatio_P=(TH1D*)h_RecoNormEff_FONLL_P->Clone();
	h_EffRatio_P->Divide(h_RecoNormEff_FONLL_P,h_RecoNormEff_P,1,1,"B");

	TH1D *h_EffRatio_NP=(TH1D*)h_RecoNormEff_FONLL_NP->Clone();
	h_EffRatio_NP->Divide(h_RecoNormEff_FONLL_NP,h_RecoNormEff_NP,1,1,"B");


	// taking the uncorrelated error for eff ratio
	for(int i=0; i<nbin_pt; i++){
		double eff_P=h_RecoNormEff_P->GetBinContent(i+1);
		double effErr_P=h_RecoNormEff_P->GetBinError(i+1);
		double eff_FONLL_P=h_RecoNormEff_FONLL_P->GetBinContent(i+1);
		double effErr_FONLL_P=h_RecoNormEff_FONLL_P->GetBinError(i+1);

		double SumErr2= pow(effErr_P/eff_P,2) + pow(effErr_FONLL_P/eff_FONLL_P,2) - 2.0*effErr_P*effErr_P/eff_P/eff_FONLL_P;
		if(SumErr2<0){
			SumErr2= pow(effErr_P/eff_P,2) + pow(effErr_FONLL_P/eff_FONLL_P,2) - 2.0*effErr_FONLL_P*effErr_FONLL_P/eff_P/eff_FONLL_P;
			if(SumErr2<0){cout<<"SumErr2 < 0"<<endl;}
		}
		if(SumErr2>=0){
			h_EffRatio_P->SetBinError(i+1, eff_FONLL_P/eff_P*sqrt(SumErr2));
		}

		double eff_NP=h_RecoNormEff_NP->GetBinContent(i+1);
		double effErr_NP=h_RecoNormEff_NP->GetBinError(i+1);
		double eff_FONLL_NP=h_RecoNormEff_FONLL_NP->GetBinContent(i+1);
		double effErr_FONLL_NP=h_RecoNormEff_FONLL_NP->GetBinError(i+1);

		SumErr2= pow(effErr_NP/eff_NP,2) + pow(effErr_FONLL_NP/eff_FONLL_NP,2) - 2.0*effErr_NP*effErr_NP/eff_NP/eff_FONLL_NP;
		if(SumErr2<0){
			SumErr2= pow(effErr_NP/eff_NP,2) + pow(effErr_FONLL_NP/eff_FONLL_NP,2) - 2.0*effErr_FONLL_NP*effErr_FONLL_NP/eff_NP/eff_FONLL_NP;
			if(SumErr2<0){cout<<"SumErr2 < 0"<<endl;}
		}
		if(SumErr2>=0){
			h_EffRatio_NP->SetBinError(i+1, eff_FONLL_NP/eff_NP*sqrt(SumErr2));
		}




	}



	gStyle->SetOptStat(0);


	TCanvas *c_Acc= new TCanvas("c_Acc","c_Acc",800,800);
	c_Acc->cd();
	gPad->SetLogx();
	h_GenACCEff_P->SetMaximum(1.02);
	h_GenACCEff_P->SetMinimum(0);
	h_GenACCEff_P->SetTitle("");
	h_GenACCEff_P->GetXaxis()->SetTitle("p_{T} (GeV)");
	h_GenACCEff_P->GetXaxis()->CenterTitle();
	h_GenACCEff_P->GetXaxis()->SetTitleOffset(1.30);
	h_GenACCEff_P->GetYaxis()->SetTitle("acceptance");
	h_GenACCEff_P->GetYaxis()->CenterTitle();
	h_GenACCEff_P->SetLineColor(2);
	h_GenACCEff_P->SetMarkerColor(2);
	h_GenACCEff_P->Draw("SAME");
	h_GenACCEff_NP->SetLineColor(4);
	h_GenACCEff_NP->SetMarkerColor(4);
	h_GenACCEff_NP->Draw("SAME");

	TLegend *le_Acc=new TLegend(0.60,0.2,0.85,0.5,NULL,"brNDC");
	le_Acc->SetBorderSize(0);
	// le_Acc->AddEntry((TObject*)0,"D_{D}^{#pm} #rightarrow #phi #pi","");
	le_Acc->AddEntry((TObject*)0,str_DsChannel_plot.Data(),"");
	le_Acc->AddEntry(h_GenACCEff_P,"Prompt D_{S}^{#pm}","lp");
	le_Acc->AddEntry(h_GenACCEff_NP,"NonPrompt D_{S}^{#pm}","lp");

	le_Acc->Draw("SAME");

	texCmsSim->Draw("SAME");
	if(isPbPb==0){
	texColpp->Draw("SAME");
	}else{
	texColPbPb->Draw("SAME");
	}

	SavePlotDirs(c_Acc,Form("Acc_%s_%s%s",str_PbPb.Data(),str_DsChannel.Data(),s_pre.Data()),{"MC"});


	TCanvas *c_Eff= new TCanvas("c_Eff","c_Eff",800,800);
	c_Eff->cd();
	gPad->SetLogx();
	// SetCanvas(c_Eff);

	if(isPbPb==3){
	h_RecoNormEff_P->GetXaxis()->SetRangeUser(6,40);
	h_RecoNormEff_NP->GetXaxis()->SetRangeUser(6,40);
	}
	h_RecoNormEff_P->SetTitle("");
	// h_RecoNormEff_P->SetMaximum(h_RecoNormEff_NP->GetMaximum()+0.01);
	h_RecoNormEff_P->SetMaximum(1.02);
	h_RecoNormEff_P->SetMinimum(0);
	h_RecoNormEff_P->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h_RecoNormEff_P->GetXaxis()->CenterTitle();
	h_RecoNormEff_P->GetYaxis()->SetTitle("#alpha #times #epsilon_{Reco} #times #epsilon_{sel}");
	h_RecoNormEff_P->GetYaxis()->CenterTitle();
	h_RecoNormEff_P->SetLineColor(2);
	h_RecoNormEff_P->SetMarkerColor(2);
	h_RecoNormEff_P->Draw("SAME");
	h_RecoNormEff_NP->SetLineColor(4);
	h_RecoNormEff_NP->SetMarkerColor(4);
	h_RecoNormEff_NP->Draw("SAME");

	TLegend *le_Eff=new TLegend(0.20,0.55,0.50,0.85,NULL,"brNDC");
	le_Eff->SetBorderSize(0);
	// le_Eff->AddEntry((TObject*)0,"D_{D}^{#pm} #rightarrow #phi #pi","");
	le_Eff->AddEntry((TObject*)0,str_DsChannel_plot.Data(),"");
	le_Eff->AddEntry(h_RecoNormEff_P,"Prompt D_{S}^{+}","lp");
	le_Eff->AddEntry(h_RecoNormEff_NP,"NonPrompt D_{S}^{+}","lp");

	le_Eff->Draw("SAME");

	texCmsSim->Draw("SAME");
	if(isPbPb==0){
	texColpp->Draw("SAME");
	}else{
	texColPbPb->Draw("SAME");
	}

	SavePlotDirs(c_Eff,Form("Eff_%s_%s%s",str_PbPb.Data(),str_DsChannel.Data(),s_pre.Data()),{"MC"},"gfc");


	double BR=BRphi;
	if(DsChannel==1){
		BR=BRf0;
	}
	TH1D *h_RecoNormEffBr_P=(TH1D*)h_RecoNormEff_P->Clone("h_RecoNormEffBr_P");
	h_RecoNormEffBr_P->Scale(BR);
	TH1D *h_RecoNormEffBr_NP=(TH1D*)h_RecoNormEff_NP->Clone("h_RecoNormEffBr_NP");
	h_RecoNormEffBr_NP->Scale(BR);

	TCanvas *c_EffBr= new TCanvas("c_EffBr","c_EffBr",800,800);
	c_EffBr->cd();
	gPad->SetLogx();

	if(isPbPb==3){
	h_RecoNormEffBr_P->GetXaxis()->SetRangeUser(6,40);
	h_RecoNormEffBr_NP->GetXaxis()->SetRangeUser(6,40);
	}
	h_RecoNormEffBr_P->SetTitle("");
	h_RecoNormEffBr_P->SetMaximum(h_RecoNormEffBr_NP->GetMaximum()*1.1);
	h_RecoNormEffBr_P->GetXaxis()->SetTitle("p_{T} (GeV)");
	h_RecoNormEffBr_P->GetXaxis()->CenterTitle();
	h_RecoNormEffBr_P->GetYaxis()->SetTitle("BR #times acceptance #times #epsilon_{Reco} #times #epsilon_{sel}");
	h_RecoNormEffBr_P->GetYaxis()->CenterTitle();
	h_RecoNormEffBr_P->SetLineColor(2);
	h_RecoNormEffBr_P->SetMarkerColor(2);
	h_RecoNormEffBr_P->Draw("SAME");
	h_RecoNormEffBr_NP->SetLineColor(4);
	h_RecoNormEffBr_NP->SetMarkerColor(4);
	h_RecoNormEffBr_NP->Draw("SAME");

	TLegend *le_EffBr=new TLegend(0.60,0.15,0.85,0.4,NULL,"brNDC");
	le_EffBr->SetBorderSize(0);
	// le_EffBr->AddEntry((TObject*)0,"D_{D}^{#pm} #rightarrow #phi #pi","");
	le_EffBr->AddEntry((TObject*)0,str_DsChannel_plot.Data(),"");
	le_EffBr->AddEntry(h_RecoNormEffBr_P,"Prompt D_{S}^{#pm}","lp");
	le_EffBr->AddEntry(h_RecoNormEffBr_NP,"NonPrompt D_{S}^{#pm}","lp");

	le_EffBr->Draw("SAME");

	texCmsSim->Draw("SAME");
	if(isPbPb==0){
	texColpp->Draw("SAME");
	}else{
	texColPbPb->Draw("SAME");
	}

	SavePlotDirs(c_EffBr,Form("EffBr_%s_%s%s",str_PbPb.Data(),str_DsChannel.Data(),s_pre.Data()),{"MC"});




	TCanvas *c_EffRatio= new TCanvas("c_EffRatio","c_EffRatio",800,800);
	c_EffRatio->cd();
	gPad->SetLogx();
	h_EffRatio_P->SetTitle("");
	h_EffRatio_P->SetMaximum(1.25);
	h_EffRatio_P->SetMinimum(0.6);
	if(isPbPb==3){
	h_EffRatio_P->GetXaxis()->SetRangeUser(6,40);
	h_EffRatio_NP->GetXaxis()->SetRangeUser(6,40);
	}
	h_EffRatio_P->GetXaxis()->SetTitle("p_{T} (GeV)");
	h_EffRatio_P->GetXaxis()->CenterTitle();
	h_EffRatio_P->GetYaxis()->SetTitle("Efficiency Ratio of Alternative to Default");
	h_EffRatio_P->GetYaxis()->CenterTitle();
	h_EffRatio_P->SetLineColor(2);
	h_EffRatio_P->SetMarkerColor(2);
	h_EffRatio_P->Draw("SAME");
	h_EffRatio_NP->SetLineColor(4);
	h_EffRatio_NP->SetMarkerColor(4);
	h_EffRatio_NP->Draw("SAME");

	TLegend *le_EffRatio=new TLegend(0.60,0.15,0.85,0.4,NULL,"brNDC");
	le_EffRatio->SetBorderSize(0);
	// le_EffRatio->AddEntry((TObject*)0,"D_{D}^{#pm} #rightarrow #phi #pi","");
	le_EffRatio->AddEntry((TObject*)0,str_DsChannel_plot.Data(),"");
	le_EffRatio->AddEntry(h_EffRatio_P,"Prompt D_{S}^{#pm}","lp");
	le_EffRatio->AddEntry(h_EffRatio_NP,"NonPrompt D_{S}^{#pm}","lp");

	le_EffRatio->Draw("SAME");

	texCmsSim->Draw("SAME");
	if(isPbPb==0){
	texColpp->Draw("SAME");
	}else{
	texColPbPb->Draw("SAME");
	}

	SavePlotDirs(c_EffRatio,Form("EffRatio_%s_%s%s",str_PbPb.Data(),str_DsChannel.Data(),s_pre.Data()),{"MC"});






}
