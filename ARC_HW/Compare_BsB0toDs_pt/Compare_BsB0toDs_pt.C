#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TSystem.h>

int Compare_BsB0toDs_pt(){


	TFile *fin=TFile::Open("/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/DsMinTree_MC_GenSampleMerge_pp_NonPrompt_phi.root");
	TTree *ntGen=(TTree*)fin->Get("ntGen");


	gStyle->SetOptStat(0);

	int nbin_pt=38;
	double bin_pt_low=2;
	double bin_pt_high=40;

	TH1D *h_BstoDs=new TH1D("h_BstoDs","h_BstoDs",nbin_pt,bin_pt_low,bin_pt_high); h_BstoDs->Sumw2();
	ntGen->Project("h_BstoDs","Gpt","weight*GptSampleWeight*(abs(Gy)<1 && GSignalType==1 && GcollisionId==0 && abs(GpdgId)==431 && abs(GBAncestorpdgId)==531 && Gpt>2 && Gpt<40 )  ");
	
	h_BstoDs->Scale(1/h_BstoDs->Integral());

	TH1D *h_B0toDs=new TH1D("h_B0toDs","h_B0toDs",nbin_pt,bin_pt_low,bin_pt_high); h_B0toDs->Sumw2();
	ntGen->Project("h_B0toDs","Gpt","weight*GptSampleWeight*(abs(Gy)<1 && GSignalType==1 && GcollisionId==0 && abs(GpdgId)==431 && abs(GBAncestorpdgId)==511 && Gpt>2 && Gpt<40 )  ");
	
	h_B0toDs->Scale(1/h_B0toDs->Integral());

	TH1D *h_BptoDs=new TH1D("h_BptoDs","h_BptoDs",nbin_pt,bin_pt_low,bin_pt_high); h_BptoDs->Sumw2();
	ntGen->Project("h_BptoDs","Gpt","weight*GptSampleWeight*(abs(Gy)<1 && GSignalType==1 && GcollisionId==0 && abs(GpdgId)==431 && abs(GBAncestorpdgId)==521 && Gpt>2 && Gpt<40 )  ");
	
	h_BptoDs->Scale(1/h_BptoDs->Integral());




	TH1D *h_B0BsRatio=new TH1D("h_B0BsRatio","h_B0BsRatio",nbin_pt,bin_pt_low,bin_pt_high); 
	h_B0BsRatio->Divide(h_B0toDs,h_BstoDs);

  TH1D *h_BpBsRatio=new TH1D("h_BpBsRatio","h_BpBsRatio",nbin_pt,bin_pt_low,bin_pt_high); 
	h_BpBsRatio->Divide(h_BptoDs,h_BstoDs);




	TCanvas *c_Dspt=new TCanvas("c_Dspt");
	c_Dspt->cd();

	h_BstoDs->SetTitle("");
	h_BstoDs->SetLineColor(1);
	h_BstoDs->Draw();
	h_B0toDs->SetLineColor(2);
	h_B0toDs->Draw("same");
	h_BptoDs->SetLineColor(4);
	h_BptoDs->Draw("same");


	c_Dspt->BuildLegend();

	TCanvas *c_ratio=new TCanvas("c_ratio");
	c_ratio->cd();
	// h_B0BsRatio->Smooth(2);
	h_B0BsRatio->SetLineColor(2);
	h_B0BsRatio->SetTitle("");
	h_B0BsRatio->Draw();
	h_BpBsRatio->SetLineColor(4);
	h_BpBsRatio->Draw("same");
	c_ratio->BuildLegend();

	// TLegend

	double c_Rmg=0.04;
	double c_Lmg=0.15;

	TCanvas *c_DsptWRatio=new TCanvas("c_DsptWRatio","c_DsptWRatio",800,800);
	c_DsptWRatio->cd();

  TPad *pad1 = new TPad("pad1","top pad",0.0,0.25,1.0,1.0);
  pad1->SetTopMargin(0.08);
  pad1->SetBottomMargin(0.0);
  pad1->SetRightMargin(c_Rmg);
  pad1->SetLeftMargin(c_Lmg);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","bottom pad",0.0,0.00,1.0,0.25);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.30);
  pad2->SetRightMargin(c_Rmg);
  pad2->SetLeftMargin(c_Lmg);
  pad2->Draw();

  pad1->cd();

	gPad->SetLogy();
	gPad->SetLogx();
	h_BstoDs->SetMarkerStyle(26);
	h_B0toDs->SetMarkerStyle(26);
	h_BptoDs->SetMarkerStyle(26);
	h_B0toDs->SetMarkerColor(2);
	h_BptoDs->SetMarkerColor(4);



	h_BstoDs->Draw();
	h_BstoDs->GetYaxis()->SetTitle("Normalized D_{S} #frac{d#sigma}{dpt}");
	h_BstoDs->GetYaxis()->CenterTitle();
	h_B0toDs->SetLineColor(2);
	h_B0toDs->Draw("same");
	h_BptoDs->SetLineColor(4);
	h_BptoDs->Draw("same");

	TLegend *le=new TLegend(0.2,0.1,0.5,0.4);
	le->SetTextSize(0.04);
	le->SetBorderSize(0);
	le->AddEntry(h_BstoDs,"B_{S} to D_{S}","lp");
	le->AddEntry(h_B0toDs,"B^{0} to D_{S}","lp");
	le->AddEntry(h_BptoDs,"B^{+} to D_{S}","lp");

	le->Draw("same");

	pad2->cd();

	gPad->SetLogy(0);
	gPad->SetLogx();

	h_B0BsRatio->SetMarkerStyle(26);
	h_BpBsRatio->SetMarkerStyle(26);
	h_B0BsRatio->SetMarkerColor(2);
	h_BpBsRatio->SetMarkerColor(4);




  h_B0BsRatio->GetXaxis()->SetTitleSize(0.12);
  h_B0BsRatio->GetXaxis()->SetTitleOffset(1.2);
  h_B0BsRatio->GetXaxis()->SetLabelSize(0.12);
  h_B0BsRatio->GetYaxis()->SetTitleSize(0.1);
  h_B0BsRatio->GetYaxis()->SetLabelSize(0.08);
  h_B0BsRatio->GetYaxis()->SetTitleOffset(0.58);

	h_B0BsRatio->GetXaxis()->SetTitle("D_{S} p_{T} (GeV/c) ");
	h_B0BsRatio->GetXaxis()->CenterTitle();
	h_B0BsRatio->GetYaxis()->CenterTitle();
	h_B0BsRatio->SetMaximum(1.7);
	h_B0BsRatio->SetMinimum(0.5);
	h_B0BsRatio->GetYaxis()->SetTitle("#frac{B^{0}(B^{+}) to D_{S} }{B_{S} to D_{S}} ");
	h_B0BsRatio->SetLineColor(2);
	h_B0BsRatio->Draw();
	h_BpBsRatio->SetLineColor(4);
	h_BpBsRatio->Draw("same");


	gSystem->Exec("mkdir -p plots");
	c_DsptWRatio->SaveAs("./plots/DsptFromBmesons.png");
	c_DsptWRatio->SaveAs("./plots/DsptFromBmesons.pdf");

	return 0;


}
