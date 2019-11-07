#include "../varCompare_para.h"


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

using namespace RooFit;


int reweightTest(int isPbPb=0, TString var_compare="DdlErr", TString var_cut="Dpt" , double var_cutLow=6, double var_cutHigh=10, double bins_var_Low=0, double bins_var_High=1, double bins_var_DrawLow=0, double bins_var_DrawHigh=0.1 ){

	TString Str_PbPb="pp";
	TString Str_PbPb3="pp";
	if(isPbPb){
		Str_PbPb="PbPb";
		Str_PbPb3="PbPb3";
	}	

  TString S_dataName=Form("../rootF/%s_fitFile.root",Str_PbPb3.Data());
  TString S_MCP=Form("../rootF/%sMC_phiPrompt_fitFile.root",Str_PbPb3.Data());
  TString S_MCNP=Form("../rootF/%sMC_phiNonPrompt_fitFile.root",Str_PbPb3.Data());
	TString S_fit=Form("../fitout/%s_%s_%s%.0fto%.0f.root",Str_PbPb.Data(),var_compare.Data(),var_cut.Data(),var_cutLow*100,var_cutHigh*100);


	TFile *f_data=TFile::Open(S_dataName);
	TFile *f_MCP=TFile::Open(S_MCP);
	TFile *f_MCNP=TFile::Open(S_MCNP);
	TFile *f_fit=TFile::Open(S_fit);

	TH1D *h_fr_PromptMC=(TH1D*)f_data->Get("h_fr_PromptMC");

	TH1D *h_var_MixMC=(TH1D*)f_fit->Get(Form("h_%s_MixMC",var_compare.Data()));


	TH1D *h_var_Data2sig=(TH1D*)f_fit->Get(Form("h_%s_Data2sig",var_compare.Data()));

	int rebinN=4;
	h_var_MixMC->Rebin(rebinN);
	h_var_Data2sig->Rebin(rebinN);

	TCanvas *c_tf1=new TCanvas("c_tf1","c_tf1");
	c_tf1->cd();

	TH1D *h_var_DataMCRatio=(TH1D*)h_var_Data2sig->Clone(Form("h_%s_DataMCRatio",var_compare.Data()));
	h_var_DataMCRatio->Divide(h_var_MixMC);

	h_var_DataMCRatio->SetMaximum(5);
	h_var_DataMCRatio->SetMinimum(-1);
	// h_var_DataMCRatio->GetXaxis()->SetRangeUser(0,0.5);

	// h_var_DataMCRatio->Rebin(4);
	h_var_DataMCRatio->GetXaxis()->SetRangeUser(0,0.1);
	h_var_DataMCRatio->Smooth(3,"R");
	h_var_DataMCRatio->Draw();



	// TF1 *f1_poly2=new TF1("f1_poly2","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x");
	// TF1 *f1_poly2=new TF1("f1_poly2","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x");
	TF1 *f1_poly2=new TF1("f1_poly2","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
	// TF1 *f1_poly2=new TF1("f1_poly2","[0]*(1+[1]*x+[2]*(2*x*x-1)+[3]*(4*x*x*x-3*x)+[4]*(8*x*x*x*x-8*x*x+1)  ) ");
	// f1_poly2->SetRange(0.005,0.2);
	f1_poly2->SetRange(0.01,0.05);
	h_var_DataMCRatio->Fit("f1_poly2","EMIS0");
	h_var_DataMCRatio->Fit("f1_poly2","EMIS0");
	h_var_DataMCRatio->Fit("f1_poly2","EMIS0");
	h_var_DataMCRatio->Fit("f1_poly2","EMIS0");
	h_var_DataMCRatio->Fit("f1_poly2","MS0");
	h_var_DataMCRatio->Fit("f1_poly2","MS0");

	f1_poly2->SetLineColor(2);

	f1_poly2->SetRange(0.0,0.1);
	f1_poly2->Draw("same");


	cout<<"test"<<endl;


	TCanvas *c_roofit=new TCanvas("c_roofit","c_roofit");
	c_roofit->cd();

	RooRealVar x("x","x",0.011,0.05);
	RooDataHist dh("dh","dh",x, Import(*h_var_DataMCRatio));
	RooPlot *frame=x.frame(Title("title"));
	dh.plotOn(frame);

	RooRealVar Cheb1("Cheb1","Cheb1",0,-1e6,1e6);
	RooRealVar Cheb2("Cheb2","Cheb2",0,-1e6,1e6);
	RooRealVar Cheb3("Cheb3","Cheb3",0,-1e6,1e6);
	RooRealVar Cheb4("Cheb4","Cheb4",0,-1e6,1e6);
	RooRealVar Cheb5("Cheb5","Cheb5",0,-1e6,1e6);
	
	RooChebychev *ChebPdf=new RooChebychev("ChebPdf","ChebPdf",x,RooArgList(Cheb1,Cheb2,Cheb3,Cheb4,Cheb5));
	// Cheb3.setConstant(kTRUE);
	// Cheb4.setConstant(kTRUE);
	Cheb5.setConstant(kTRUE);

	ChebPdf->fitTo(dh,SumW2Error(kTRUE));
	ChebPdf->fitTo(dh,SumW2Error(kTRUE));
	ChebPdf->fitTo(dh,SumW2Error(kTRUE));
	ChebPdf->fitTo(dh,SumW2Error(kTRUE));
	ChebPdf->fitTo(dh,SumW2Error(kTRUE));
	// ChebPdf->fitTo(dh);
	// ChebPdf->fitTo(dh);
	ChebPdf->plotOn(frame,LineColor(2));

	frame->SetMaximum(5);
	frame->SetMinimum(-1);
	frame->Draw();


	return 0;

}

int main(int argc, char*argv[]){

  if(argc==6){
    reweightTest(atoi(argv[1]), argv[2], argv[3], atof(argv[4]), atof(argv[5]) );
  }else if(argc==10){
    reweightTest(atoi(argv[1]), argv[2], argv[3], atof(argv[4]), atof(argv[5]) , atof(argv[6]) , atof(argv[7]) , atof(argv[8]) , atof(argv[9]));
  }else{
    reweightTest();
    cout<<"wrong number of input parameters , need 6\n int isPbPb=0, TString var_compare=DdxyzErr, TString var_cut=Dtrk1Pt , double var_cutLow=0.75, double var_cutHigh=1.25"<<endl;
    return 1;
  }
//int Fit_sideband(int isPbPb=0, TString var_compare="DdxyzErr", TString var_cut="Dtrk1Pt" , double var_cutLow=0.75, double var_cutHigh=1.25){


  return 0;
}

