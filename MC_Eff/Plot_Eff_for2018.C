#include "../include/uti.h"
#include "../include/DsMinTreeLoad.h"
#include "../include/parameters.h"


#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"


#include <TString.h>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>
#include <TLatex.h>


void Plot_Eff_for2018(Int_t isPbPb=3, Int_t DsChannel=0){

	InitStyle();
  double tex_upperY=0.93;

  TLatex* texCmsPre = new TLatex(0.15,tex_upperY, "#scale[1.25]{CMS} Preliminary");
  texCmsPre->SetNDC();
  texCmsPre->SetTextAlign(12);
  texCmsPre->SetTextSize(0.035);
  texCmsPre->SetTextFont(42);

  TLatex* texCmsSim = new TLatex(0.15,tex_upperY, "#scale[1.25]{CMS} Simulations");
  texCmsSim->SetNDC();
  texCmsSim->SetTextAlign(12);
  texCmsSim->SetTextSize(0.035);
  texCmsSim->SetTextFont(42);


  TLatex* texColPbPb = new TLatex(0.88,tex_upperY, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  texColPbPb->SetNDC();
  texColPbPb->SetTextAlign(32);
  texColPbPb->SetTextSize(0.035);
  texColPbPb->SetTextFont(42);

  TLatex* texColpp = new TLatex(0.88,tex_upperY, "pp #sqrt{s} = 5.02 TeV");
  texColpp->SetNDC();
  texColpp->SetTextAlign(32);
  texColpp->SetTextSize(0.035);
  texColpp->SetTextFont(42);

	TString str_PbPb="pp";
	if(isPbPb==3){
		str_PbPb="PbPb3";
	}
	TString str_DsChannel="phikkpi";
	TString str_DsChannel_plot="D_{S}^{#pm} #rightarrow #phi #pi^{#pm}";
	if(DsChannel==1){
		str_DsChannel="f0kkpi";
		str_DsChannel_plot="D_{S}^{#pm} #rightarrow f0 #pi^{#pm}";
	}

	TString fname_NP=Form("./output/MC_eff_%s_NonPrompt_%s.root", str_PbPb.Data(),str_DsChannel.Data());
	TString fname_P=Form("./output/MC_eff_%s_Prompt_%s.root", str_PbPb.Data(),str_DsChannel.Data());

	TFile *f_P=TFile::Open(fname_P,"READ");
	TFile *f_NP=TFile::Open(fname_NP,"READ");

	TH1D *h_GenACCEff_P=(TH1D*)f_P->Get("h_GenAccEff");
	TH1D *h_GenACCEff_NP=(TH1D*)f_NP->Get("h_GenAccEff");

	TH1D *h_RecoNormEff_P=(TH1D*)f_P->Get("h_RecoNormEff");
	TH1D *h_RecoNormEff_NP=(TH1D*)f_NP->Get("h_RecoNormEff");

	TH1D *h_RecoNormEff_FONLL_P=(TH1D*)f_P->Get("h_RecoNormEff_FONLL");
	TH1D *h_RecoNormEff_FONLL_NP=(TH1D*)f_NP->Get("h_RecoNormEff_FONLL");


	TH1D *h_EffRatio_P=(TH1D*)h_RecoNormEff_FONLL_P->Clone();
	h_EffRatio_P->Divide(h_RecoNormEff_P);

	TH1D *h_EffRatio_NP=(TH1D*)h_RecoNormEff_FONLL_NP->Clone();
	h_EffRatio_NP->Divide(h_RecoNormEff_FONLL_NP,h_RecoNormEff_NP,1,1,"B");


	gStyle->SetOptStat(0);


	TCanvas *c_Acc= new TCanvas("c_Acc","c_Acc",800,800);
	c_Acc->cd();
	gPad->SetLogx();
	h_GenACCEff_P->SetMaximum(1.05);
	h_GenACCEff_P->SetTitle("");
	h_GenACCEff_P->GetXaxis()->SetTitle("p_{T} (GeV)");
	h_GenACCEff_P->GetXaxis()->CenterTitle();
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

	// SavePlotDirs(c_Acc,Form("Acc_%s_%s",str_PbPb.Data(),str_DsChannel.Data()),{"MC"});


	TCanvas *c_Eff= new TCanvas("c_Eff","c_Eff",800,800);
	c_Eff->cd();
	gPad->SetLogx();
	h_RecoNormEff_P->GetXaxis()->SetRangeUser(5,40);
	h_RecoNormEff_P->SetTitle("");
	h_RecoNormEff_P->SetMaximum(h_RecoNormEff_NP->GetMaximum()+0.01);
	h_RecoNormEff_P->GetXaxis()->SetTitle("p_{T} (GeV)");
	h_RecoNormEff_P->GetXaxis()->CenterTitle();
	h_RecoNormEff_P->GetYaxis()->SetTitle("acceptance #times #epsilon_{Reco} #times #epsilon_{sel}");
	h_RecoNormEff_P->GetYaxis()->CenterTitle();
	h_RecoNormEff_P->SetLineColor(4);
	h_RecoNormEff_P->SetMarkerColor(4);
	h_RecoNormEff_P->Draw("SAME");
	// h_RecoNormEff_NP->SetLineColor(4);
	// h_RecoNormEff_NP->SetMarkerColor(4);
	// h_RecoNormEff_NP->Draw("SAME");

	TLegend *le_Eff=new TLegend(0.60,0.15,0.85,0.4,NULL,"brNDC");
	le_Eff->SetBorderSize(0);
	// le_Eff->AddEntry((TObject*)0,"D_{D}^{#pm} #rightarrow #phi #pi","");
	le_Eff->AddEntry((TObject*)0,str_DsChannel_plot.Data(),"");
	le_Eff->AddEntry(h_RecoNormEff_P,"Prompt D_{S}^{#pm}","");
	// le_Eff->AddEntry(h_RecoNormEff_NP,"NonPrompt D_{S}^{#pm}","lp");

	le_Eff->Draw("SAME");

	// texCmsSim->Draw("SAME");
	// if(isPbPb==0){
	// texColpp->Draw("SAME");
	// }else{
	texColPbPb->Draw("SAME");
	// }

	for(int i=1; i<=h_RecoNormEff_P->GetNbinsX() ;i++){
		cout<<"i = "<<i<<" , Eff = "<<h_RecoNormEff_P->GetBinContent(i)<<" ,bincenter : "<<h_RecoNormEff_P->GetBinCenter(i)<<endl;
	}

	// SavePlotDirs(c_Eff,Form("Eff_%s_%s",str_PbPb.Data(),str_DsChannel.Data()),{"MC"});


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

	// SavePlotDirs(c_EffRatio,Form("EffRatio_%s_%s",str_PbPb.Data(),str_DsChannel.Data()),{"MC"});






}
