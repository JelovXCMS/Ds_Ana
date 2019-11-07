#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"

int compareDs(){

	TFile *fD0=new TFile("fPPMB_Chi2_DoubleRatio.root");
	TFile *fDs=new TFile("CutScanSys_FixShape_pp.root");

	TH1F *hD0_chi2cl=(TH1F*)fD0->Get("hDoubleRatio");
	TH1D *hDs_chi2cl=(TH1D*)fDs->Get("h_PromptDs_Dchi2clMinScan_ratio_pt10to20");
	TH1D *hPhi_chi2cl=(TH1D*)fDs->Get("h_PromptDs_Reschi2clScan_ratio_pt10to20");


	TH1D *hDs_StatError=(TH1D*)fDs->Get("h_PromptDs_Dchi2clMinScan_Stat_pt10to20");
	TH1D *hPhi_StatError=(TH1D*)fDs->Get("h_PromptDs_Reschi2clScan_Stat_pt10to20");

	hDs_chi2cl->SetBinError(1,hDs_StatError->GetBinError(1));
	hPhi_chi2cl->SetBinError(1,hPhi_StatError->GetBinError(1));

	TCanvas *c2=new TCanvas("c2","c2",800,600);
	c2->cd();

	hDs_chi2cl->SetMinimum(0.85);
	hDs_chi2cl->SetMaximum(1.4);
	hDs_chi2cl->SetTitle("");
	hDs_chi2cl->GetXaxis()->SetTitle("vertex probability");
	hDs_chi2cl->GetYaxis()->SetTitle("#sigma_{varied}/#sigma_{loose}");
	hDs_chi2cl->GetXaxis()->SetTitleOffset(1.0);
	hDs_chi2cl->GetYaxis()->SetTitleOffset(1.0);
	hDs_chi2cl->SetLineColor(2);
	hDs_chi2cl->SetMarkerColor(2);
	hDs_chi2cl->SetMarkerStyle(26);
	hDs_chi2cl->Draw("same");
	hPhi_chi2cl->SetLineColor(4);
	hPhi_chi2cl->SetMarkerColor(4);
	hPhi_chi2cl->SetMarkerStyle(27);
	hPhi_chi2cl->Draw("same");
	hD0_chi2cl->SetMarkerStyle(28);
	hD0_chi2cl->Draw("same");

	TLegend *le=new TLegend(0.2,0.65,0.5,0.85);
	le->SetBorderSize(0);
	le->AddEntry(hDs_chi2cl,"D_{S}","pl");
	le->AddEntry(hPhi_chi2cl,"#phi","pl");
	le->AddEntry(hD0_chi2cl,"D^{0}","pl");
	le->Draw("same");

	TLatex *tl=new TLatex();
	tl->DrawLatexNDC(0.5,0.75,"pp p_{T}:10-20 GeV/c");

	gPad->Update();

	// TGraphErrors *gr_D0=new TGraphErrors();

	







	return 0;

}
