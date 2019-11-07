#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting_simple.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/CMS_lumi.C"
#include "varCompare_para.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TFitter.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TGraphErrors.h"

// double textposx=0.2;
// double textposy=0.77;

double shiftY=0.06;
double shiftX=0.3;
double oneshift=0.075;

int reWeight=0;
TString s_reWeight="DdxyzErrWeight";



int Plot_PromptFr3(int isPbPb=0,TString var_scan="Ddls", int useReweight=0, TString ReweightTree="DdxyzErrWeight"){

  initParameter();
  setTDRStyle();
  InitStyle();
	gStyle->SetOptStat(0);
	
	TString s_weight="";
	if(useReweight){
		s_weight="_"+ReweightTree;
	} 
	TString s_FixShape="FixShape";
	TString s_PFrScan="_PFrScan";
	s_PFrScan="";

	int nbin_pt=nbin_pt_pp;
	double *bins_pt=bins_pt_pp;
	TString s_ppPbPb="pp";
	int startBin=0;
	double Xlow=0;
	double Xhigh=1;

	TFile *fD0=new TFile("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BtoDs_Results/output/BtoDs.root","read");

	TH1D *h_CS_PromptD0_pp=(TH1D*)fD0->Get("hPromptDCrossSectionPP_AnaBin");
	TH1D *h_CS_BtoD0_pp=(TH1D*)fD0->Get("hBtoDCrossSectionPP_AnaBin");

	TH1D *h_CS_PromptD0_PbPb=(TH1D*)fD0->Get("hPromptDdNdPtPbPb_AnaBin");
	TH1D *h_CS_BtoD0_PbPb=(TH1D*)fD0->Get("hBtoDdNdPtPbPb_AnaBin");

	TH1D *h_CS_PromptD0=h_CS_PromptD0_pp;
	TH1D *h_CS_BtoD0=h_CS_BtoD0_pp;	

	if(isPbPb){
		nbin_pt=nbin_pt_PbPb3;
		bins_pt=bins_pt_PbPb3;
		s_ppPbPb="PbPb";
		startBin=2;
		Xlow=bins_pt[startBin];
		Xhigh=bins_pt[nbin_pt];

		h_CS_PromptD0=h_CS_PromptD0_PbPb;
		h_CS_BtoD0=h_CS_BtoD0_PbPb;

	}


	TH1D *h_CS_PromptD0_fr=(TH1D*)h_CS_PromptD0->Clone("h_CS_PromptD0_fr");
	TH1D *h_CS_Sum=(TH1D*)h_CS_PromptD0->Clone("h_CS_Sum");
	h_CS_Sum->Add(h_CS_BtoD0);
	h_CS_PromptD0_fr->Divide(h_CS_Sum);
	h_CS_PromptD0_fr->SetLineColor(4);
	h_CS_PromptD0_fr->SetMarkerColor(4);
	h_CS_PromptD0_fr->SetMarkerStyle(22);

	TString s_fin="";
	TFile *fin[nbin_pt];
	TGraphAsymmErrors *gr_PFrDefault[nbin_pt];
	TGraphAsymmErrors *gr_PFrDefault_Sys[nbin_pt];
	TGraphErrors *gr_ScaleBestPfr[nbin_pt];
	TGraphAsymmErrors *gr_ScaleBest[nbin_pt];
	TGraphErrors *gr_Fit_PNP[nbin_pt];
	// TGraphErrors *gr_PFrBest[nbin_pt];

	TCanvas *c_fr=new TCanvas("c_fr","",800,800);
	c_fr->cd();
	gPad->SetLogx();
	TH1D *htemp=new TH1D("htemp",";p_{T} (GeV); Prompt D_{s} fraction",nbin_pt,bins_pt);
	htemp->GetXaxis()->SetRangeUser(Xlow,Xhigh);
	htemp->SetMaximum(1.3);
	htemp->SetMinimum(0.95);
	htemp->GetXaxis()->SetTitleOffset(1);
	htemp->GetYaxis()->SetTitleOffset(1);
	htemp->Draw();

	// for(int i=0; i<1; i++)
	for(int i=startBin; i<nbin_pt ; i++)
	{
		// s_fin=Form("./RatioOut_smear/%s%s/%s_%s_Dpt%.0fto%.0f%s%s.root",var_scan.Data(),s_weight.Data(),s_ppPbPb.Data(),var_scan.Data(),bins_pt[i],bins_pt[i+1],s_FixShape.Data(),s_PFrScan.Data());
		s_fin=Form("../sPlot_varCompare_v2/DdlErrScaleOut/%s_Dpt%.0fto%.0f_DdlErrScale.root",s_ppPbPb.Data(),bins_pt[i],bins_pt[i+1]);

		fin[i]=TFile::Open(s_fin.Data(),"read");
/*
		gr_PFrDefault[i]=(TGraphAsymmErrors*)fin[i]->Get("gr_PFrDefault")	;
		gr_PFrDefault[i]->Draw("p");

		gr_PFrDefault_Sys[i]=(TGraphAsymmErrors*)fin[i]->Get("gr_PFrDefault_Sys")	;
		gr_PFrDefault_Sys[i]->SetLineColor(1);
		gr_PFrDefault_Sys[i]->SetFillColor(kGray+2);
		gr_PFrDefault_Sys[i]->SetFillStyle(3002);
		gr_PFrDefault_Sys[i]->Draw("pF2");

		gr_Fit_PNP[i]=(TGraphErrors*)fin[i]->Get("gr_Fit_PNP");
		gr_Fit_PNP[i]->SetLineColor(2);
		gr_Fit_PNP[i]->SetMarkerColor(2);
		gr_Fit_PNP[i]->Draw("p");

*/
		// gr_PFrBest[i]=(TGraphErrors*)fin[i]->Get("gr_PFrBest");
		// gr_PFrBest[i]->SetLineColor(2);
		// gr_PFrBest[i]->SetMarkerColor(2);
		// gr_PFrBest[i]->Draw("p");
/*
		gr_ScaleBestPfr[i]=(TGraphErrors*)fin[i]->Get("gr_ScaleBestPfr");
		gr_ScaleBestPfr[i]->SetLineColor(2);
		gr_ScaleBestPfr[i]->SetMarkerColor(2);
		gr_ScaleBestPfr[i]->Draw("p");
*/	

	
	}

/*

	TLegend *le=new TLegend(0.25,0.25,0.5,0.5);
	le->SetBorderSize(0);
	le->AddEntry((TObject*)0,s_ppPbPb.Data(),"");
	le->AddEntry(gr_PFrDefault[nbin_pt-1],"Default Fr.","p");
	le->AddEntry(gr_Fit_PNP[nbin_pt-1],Form("Best Fr.%s",var_scan.Data()), "p");
	if(useReweight){
		le->AddEntry((TObject*)0,ReweightTree.Data(), "");
	}
	le->Draw("same");

	gSystem->Exec(Form("mkdir -p Ratioplots_BestScale/%s%s",var_scan.Data(),s_weight.Data()));
	c_fr->SaveAs(Form("./Ratioplots_BestScale/%s%s/%s_%s_%s%s_BestScalePfr.png",var_scan.Data(),s_weight.Data(),s_ppPbPb.Data(),var_scan.Data(),s_FixShape.Data(),s_PFrScan.Data()));

	if(isPbPb){
		h_CS_PromptD0_fr->GetXaxis()->SetRangeUser(6,40);
	}

	h_CS_PromptD0_fr->Draw("same");
	le->AddEntry(h_CS_PromptD0_fr,"D^{0} fraction","p");
	le->Draw("same");
	c_fr->SaveAs(Form("./Ratioplots_BestScale/%s%s/%s_%s_%s%s_BestScalePfr_withD0.png",var_scan.Data(),s_weight.Data(),s_ppPbPb.Data(),var_scan.Data(),s_FixShape.Data(),s_PFrScan.Data()));

*/
	TCanvas *c_scl=new TCanvas("c_scl","",800,800);
	c_scl->cd();
	gPad->SetLogx();
	TH1D *htemp1=new TH1D("htemp1",";p_{T} (GeV); Best DdlErr scale",nbin_pt,bins_pt);
	htemp1->GetXaxis()->SetRangeUser(Xlow,Xhigh);
	htemp1->SetMaximum(1.5);
	htemp1->SetMinimum(0.7);
	htemp1->GetXaxis()->SetTitleOffset(1);
	htemp1->GetYaxis()->SetTitleOffset(1);
	htemp1->Draw();

	TLatex *tlx=new TLatex();
	tlx->SetTextSize(25);
	tlx->DrawLatexNDC(0.25,0.77,s_ppPbPb.Data());


	for(int i=startBin; i<nbin_pt ; i++)
	{

		gr_ScaleBest[i]=(TGraphAsymmErrors*)fin[i]->Get("gr_Best_Scale");
		gr_ScaleBest[i]->Draw("p");
		double x0,y0;
		gr_ScaleBest[i]->GetPoint(0,x0,y0);
		cout<<"scale = "<<y0<<endl;
	}

	gSystem->Exec(Form("mkdir -p Ratioplots_BestScale/%s%s",var_scan.Data(),s_weight.Data()));
	c_scl->SaveAs(Form("./Ratioplots_BestScale/%s%s/%s_%s_%s%s_BestScale.png",var_scan.Data(),s_weight.Data(),s_ppPbPb.Data(),var_scan.Data(),s_FixShape.Data(),s_PFrScan.Data()));





	return 0;
}
