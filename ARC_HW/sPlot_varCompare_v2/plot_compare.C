#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"

// #include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
// #include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
// #include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"

#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"
#include <TString.h>
#include <TCanvas.h>
#include "TBranch.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLatex.h"
#include <TH2F.h>
#include <TCut.h>
#include <TStyle.h>

#include <math.h>

#include <string>
#include <TChain.h>
#include <TFile.h>
#include "TH1.h"
#include "TROOT.h"
#include "TPad.h"
#include "TRandom.h"
#include "THStack.h"
#include "TH2.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TAxis.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <TCut.h>
#include <RooFit.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooHistPdf.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooGlobalFunc.h>
#include <RooGaussian.h>
#include <RooPolynomial.h>
#include <RooGenericPdf.h>
#include <RooFormulaVar.h>
#include <RooChebychev.h>
#include <RooExponential.h>
#include <RooFitResult.h>
#include <RooWorkspace.h>
#include <RooMCStudy.h>
#include "RooHist.h"
#include "RooConstVar.h"
#include "RooMsgService.h"

#include "RooStats/SPlot.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"

#include "varCompare_para.h"


	// double bins_varcut[]={0.75,2,3, 4, 5, 6, 8, 12};
	// double bins_varcut[]={1,2,3, 4, 5, 6, 8, 12};
	// double bins_varcut[]={2,3,4,5,6,8,10,20,40};
	double bins_varcut[]={2,3,4,5,6,8,10,20,40};
	// double bins_varcut[]={2,4,6,10,40};
	// double bins_varcut[]={6,10,40};
	int nbin_varcut=sizeof(bins_varcut)/sizeof(bins_varcut[0])-1;

// int Fit_sideband(int isPbPb=0, TString var_compare="DdxyzErr", TString var_cut="Dtrk1Pt" , double var_cutLow=0.75, double var_cutHigh=1.25){

int plot_compare(int isPbPb=0, TString var_compare="DdxyzErr", TString var_cut="Dtrk1Pt"){

	gStyle->SetOptStat(0);

	TString Str_ppPbPb="pp";
	if(isPbPb==3){
		Str_ppPbPb="PbPb";
	}

	// TFile *fout

	TH1D *h_Data_var=new TH1D("h_Data_var",Form("Data;%s;%s",var_cut.Data(),var_compare.Data()), nbin_varcut,bins_varcut); 
	h_Data_var->Sumw2();
	TH1D *h_PromptMC_var=new TH1D("h_PromptMC_var",Form("Prompt MC;%s;%s",var_cut.Data(),var_compare.Data()), nbin_varcut,bins_varcut); h_PromptMC_var->Sumw2();
	TH1D *h_NonPromptMC_var=new TH1D("h_NonPromptMC_var",Form("NonPrompt MC;%s;%s",var_cut.Data(),var_compare.Data()), nbin_varcut,bins_varcut); h_NonPromptMC_var->Sumw2();
	TH1D *h_MixMC_var=new TH1D("h_MixMC_var",Form("Mix MC;%s;%s",var_cut.Data(),var_compare.Data()), nbin_varcut,bins_varcut); h_MixMC_var->Sumw2();


	TFile *fin[nbin_varcut];
	TH1D  *h_Data[nbin_varcut];
	TH1D  *h_PromptMC[nbin_varcut];
	TH1D  *h_NonPromptMC[nbin_varcut];
	TH1D  *h_MixMC[nbin_varcut];


	for(int i=0; i<nbin_varcut; i++){

		fin[i]=TFile::Open(Form("./fitout/%s_%s_%s%.0fto%.0f.root",Str_ppPbPb.Data(),var_compare.Data(),var_cut.Data(), bins_varcut[i]*100, bins_varcut[i+1]*100 ),"READ");
		h_Data[i]=(TH1D*)fin[i]->Get(Form("h_%s_Data2sig",var_compare.Data()));
		cout<<"i = "<<i<<" mean = "<<h_Data[i]->GetMean()<<" , Std = "<<h_Data[i]->GetStdDev()<<endl;
		h_Data_var->SetBinContent(i+1, h_Data[i]->GetMean());
		// h_Data_var->SetBinError(i+1, h_Data[i]->GetStdDev());
		h_Data_var->SetBinError(i+1, h_Data[i]->GetMeanError());
	
		h_PromptMC[i]=(TH1D*)fin[i]->Get(Form("h_%s_PromptMC",var_compare.Data()));
		h_PromptMC_var->SetBinContent(i+1, h_PromptMC[i]->GetMean());
		// h_PromptMC_var->SetBinError(i+1, h_PromptMC[i]->GetStdDev());
		h_PromptMC_var->SetBinError(i+1, h_PromptMC[i]->GetMeanError());
	
		h_NonPromptMC[i]=(TH1D*)fin[i]->Get(Form("h_%s_NonPromptMC",var_compare.Data()));
		h_NonPromptMC_var->SetBinContent(i+1, h_NonPromptMC[i]->GetMean());
		// h_NonPromptMC_var->SetBinError(i+1, h_NonPromptMC[i]->GetStdDev());
		h_NonPromptMC_var->SetBinError(i+1, h_NonPromptMC[i]->GetMeanError());

		h_MixMC[i]=(TH1D*)fin[i]->Get(Form("h_%s_MixMC",var_compare.Data()));
		h_MixMC_var->SetBinContent(i+1, h_MixMC[i]->GetMean());
		// h_MixMC_var->SetBinError(i+1, h_MixMC[i]->GetStdDev());
		h_MixMC_var->SetBinError(i+1, h_MixMC[i]->GetMeanError());
	
	
	}

	TCanvas *c_var= new TCanvas("c_var","c_var",600,600);
	c_var->cd();

	if(var_cut=="Dpt"){
		gPad->SetLogx();
	}

	h_Data_var->SetTitle("");
	h_Data_var->GetYaxis()->SetTitle(Form("%s",var_compare.Data()));
	h_Data_var->GetYaxis()->SetTitleOffset(1.4);

	h_Data_var->SetLineColor(1);	
	h_Data_var->SetMarkerColor(1);	
	h_Data_var->SetMarkerStyle(29);	

	h_PromptMC_var->SetLineColor(2);	
	h_PromptMC_var->SetMarkerColor(2);	
	h_PromptMC_var->SetMarkerStyle(26);	

	h_NonPromptMC_var->SetLineColor(4);	
	h_NonPromptMC_var->SetMarkerColor(4);	
	h_NonPromptMC_var->SetMarkerStyle(26);	

	h_MixMC_var->SetLineColor(1);	
	h_MixMC_var->SetMarkerColor(1);	
	h_MixMC_var->SetMarkerStyle(26);	

	h_MixMC_var->SetLineColor(2);	
	h_MixMC_var->SetMarkerColor(2);	
	h_MixMC_var->SetMarkerStyle(26);	

	if(h_MixMC_var->GetMaximum()>h_Data_var->GetMaximum()){
		h_Data_var->SetMaximum(h_MixMC_var->GetMaximum()*1.1);
	}
	h_Data_var->Draw();
	// h_PromptMC_var->Draw("same");
	// h_NonPromptMC_var->Draw("same");
	h_MixMC_var->Draw("same");

	double lex1=0.6;
	double lex2=0.85;
	double ley1=0.2;
	double ley2=0.45;

	if(isPbPb){
		ley1=0.2;
		ley2=0.45;	

	}

	TLegend *le=new TLegend(lex1,ley1,lex2,ley2);
	le->SetBorderSize(0);
	le->AddEntry(h_Data_var,Form("%s Data",Str_ppPbPb.Data()),"pl");
	// le->AddEntry(h_PromptMC_var,"Prompt MC","pl");
	// le->AddEntry(h_NonPromptMC_var,"NonPrompt MC","pl");
	le->AddEntry(h_MixMC_var,"Mix MC","pl");
	le->Draw("same");


	c_var->SaveAs(Form("./plots/%s/%s_%s_compare.png",var_compare.Data(),Str_ppPbPb.Data(),var_cut.Data()));


	h_MixMC_var->SetLineColor(1);	
	h_MixMC_var->SetMarkerColor(1);	
	h_MixMC_var->SetMarkerStyle(26);	

	h_PromptMC_var->Draw("same");
	h_NonPromptMC_var->Draw("same");
	le->AddEntry(h_PromptMC_var,"Prompt MC","pl");
	le->AddEntry(h_NonPromptMC_var,"NonPrompt MC","pl");
	le->Draw("same");
	c_var->SaveAs(Form("./plots/%s/%s_%s_compareWithPNP.png",var_compare.Data(),Str_ppPbPb.Data(),var_cut.Data()));

	return 0;

}


int main(int argc, char *argv[]){

	if(argc==4){
		plot_compare(atoi(argv[1]), argv[2], argv[3] );
	}else{
		cout<<"incorrect number of argc, need 4 "<<endl;
		return 1;
	}
	

// int plot_compare(int isPbPb=0, TString var_compare="DdxyzErr", TString var_cut="Dtrk1Pt"){

	

	return 0;
}

