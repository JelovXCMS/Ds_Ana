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

	// double bins_varcut[]={2,4,6,10,40};
	// double bins_varcut[]={6,10,40};
	// int nbin_varcut=sizeof(bins_varcut)/sizeof(bins_varcut[0])-1;

int makeReweight_DdxyzErr(int isPbPb=3, TString var_compare="DdxyzErr", TString var_cut="Dpt" , double var_cutLow=6, double var_cutHigh=10, double bins_var_Low=0, double bins_var_High=1, double bins_var_DrawLow=0, double bins_var_DrawHigh=0.1 ){

	double bins_varcut_pp[]={2,4,6,10,40};
	double bins_varcut_PbPb[]={6,10,40};
	double *bins_varcut=bins_varcut_pp;

	// int nbin_varcut=sizeof(bins_varcut)/sizeof(bins_varcut[0])-1;
	int nbin_varcut=4;
	TString Str_PbPb="pp";
	TString Str_PbPb3="pp";
	if(isPbPb){
		Str_PbPb="PbPb";
		Str_PbPb3="PbPb3";
	  bins_varcut=bins_varcut_PbPb;
		nbin_varcut=2;
	}	

  TString S_dataName=Form("../rootF/%s_fitFile.root",Str_PbPb3.Data());

	TString S_fDir="/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_cutScan/";
	// TString S_fDir="../rootF/";

  TString S_MCP=Form("%s%sMC_phiPrompt_fitFile.root",S_fDir.Data(),Str_PbPb3.Data());
  TString S_MCNP=Form("%s%sMC_phiNonPrompt_fitFile.root",S_fDir.Data(),Str_PbPb3.Data());
  TString S_MCP_new=Form("%s%sMC_phiPrompt_fitFile_%sWeight.root",S_fDir.Data(),Str_PbPb3.Data(),var_compare.Data());
  TString S_MCNP_new=Form("%s%sMC_phiNonPrompt_fitFile_%sWeight.root",S_fDir.Data(),Str_PbPb3.Data(),var_compare.Data());



	TFile *f_data=TFile::Open(S_dataName);
	TFile *f_MCP=TFile::Open(S_MCP);
	TFile *f_MCNP=TFile::Open(S_MCNP);

	TFile *f_MCP_new=TFile::Open(S_MCP_new,"recreate");
	TFile *f_MCNP_new=TFile::Open(S_MCNP_new,"recreate");
//	read ratio histograms

	TFile *f_fit[nbin_varcut];
	TH1D *h_ratio[nbin_varcut];

	// TString S_fit=Form("../fitout/%s_%s_%s%.0fto%.0f.root",Str_PbPb.Data(),var_compare.Data(),var_cut.Data(),var_cutLow*100,var_cutHigh*100);

	TCanvas *c_test[nbin_varcut];

	for(int i=0; i<nbin_varcut; i++){
		f_fit[i]=TFile::Open(Form("../fitout/%s_%s_%s%.0fto%.0f.root",Str_PbPb.Data(),var_compare.Data(),var_cut.Data(),bins_varcut[i]*100,bins_varcut[i+1]*100));
		h_ratio[i]=(TH1D*)f_fit[i]->Get(Form("h_%s_ratioSmooth",var_compare.Data()));
		
		c_test[i]=new TCanvas(Form("c_test%i",i), Form("c_test%i",i));
		c_test[i]->cd();
		h_ratio[i]->SetMaximum(5);
		h_ratio[i]->SetMinimum(-1);
		h_ratio[i]->GetXaxis()->SetRangeUser(0,0.05);
		h_ratio[i]->Draw();

	}
	double WeightVarMin=0.002;
	double WeightVarMax=0.1;
	double DptLowest=2;
	double DptHighest=40;

	// setbranch
	Float_t Dpt;
	Float_t WeightVar;
//	double TotalWeight;
	Float_t VarWeight;	


	TH1D *hGen_pt_MCP=(TH1D*)f_MCP->Get("hGen_pt");

	TTree *t_fit_MCP=(TTree*)f_MCP->Get("t_fit");
	t_fit_MCP->SetBranchAddress("Dpt", &Dpt);
	t_fit_MCP->SetBranchAddress(Form("%s",var_compare.Data()), &WeightVar);

	// TFile *f_MCP_new=TFile::Open(Form("../rootF/%sMC_phiPrompt_fitFile_%sWeight.root",Str_PbPb3.Data(),var_compare.Data()),"recreate");
	f_MCP_new->cd();
	TTree *t_fit_MCP_new=(TTree*)t_fit_MCP->CloneTree(0);

	// t_fit_new->Branch("DdlErrWeight",&DdlErrWeight,"DdlErrWeight/F");	
	t_fit_MCP_new->Branch(Form("%sWeight",var_compare.Data()),&VarWeight);	


	// loop events to add weight

	// double NorFactor[nbin_varcut];
	double NorFactor[4];
	NorFactor[0]	=1;
	NorFactor[1]	=1;
	NorFactor[2]	=1;
	NorFactor[3]	=1;
	if(isPbPb==0){
	NorFactor[0]=1.78111/1.89362;
	NorFactor[1]=2.13961/2.21607;
	NorFactor[2]=1.31111/1.29016;
	NorFactor[3]=3.52482/3.31269;
	}else if(isPbPb==3){
	NorFactor[0]=0.490175/0.588831;
	NorFactor[1]=0.149462/0.148137;
	NorFactor[0]=0.571728/0.680722;
	NorFactor[1]=0.161897/0.15997;
	}


	long nentries=t_fit_MCP->GetEntries();
	// nentries=10000;
	for(int i=0; i<nentries;i++){
		t_fit_MCP->GetEntry(i);
		VarWeight=1;
		if(Dpt>DptLowest && Dpt<DptHighest && WeightVar >WeightVarMin && WeightVar< WeightVarMax){
		for(int j=0; j<nbin_varcut; j++){
			if(Dpt>bins_varcut[j] && Dpt<=bins_varcut[j+1]){
				VarWeight=h_ratio[j]->GetBinContent(h_ratio[j]->FindBin(WeightVar))* NorFactor[j];
				if(VarWeight<0.001){	VarWeight=0.001;}
				if(VarWeight>5){	VarWeight=5;}
				// cout<<"DdlErrWeight = "<<DdlErrWeight<<endl;	
				break;
			}
		} // for j<nbin_varcut
		} // if Dpt>2 && <40

		t_fit_MCP_new->Fill();


	}

	f_MCP_new->cd();
	t_fit_MCP_new->Write();
	hGen_pt_MCP->Write();

// check the total integral after weight

	TH1D *h_before_MCP[nbin_varcut];
	TH1D *h_after_MCP[nbin_varcut];
	TCanvas *c_weight_MCP[nbin_varcut];

	int nbin_test=200;
	double binLow_test=0;
	double binHigh_test=200;
	TString s_test="Ddls";



	for(int i=0; i<nbin_varcut;i++){
		h_before_MCP[i]=new TH1D(Form("h_before_MCP%i",i),"",nbin_test,binLow_test,binHigh_test); h_before_MCP[i]->Sumw2();
		h_after_MCP[i]=new TH1D(Form("h_after_MCP%i",i),"",nbin_test,binLow_test,binHigh_test); h_after_MCP[i]->Sumw2();
	
		t_fit_MCP_new->Project(Form("h_before_MCP%i",i),s_test.Data(),Form("TotalWeight*(Dpt>%f && Dpt<%f)",bins_varcut[i],bins_varcut[i+1]));
		t_fit_MCP_new->Project(Form("h_after_MCP%i",i),s_test.Data(),Form("%sWeight*TotalWeight*(Dpt>%f && Dpt<%f)",var_compare.Data(),bins_varcut[i],bins_varcut[i+1]));


		c_weight_MCP[i]=new TCanvas(Form("c_weight_MCP%i",i),Form("c_weight_MCP%i",i));
		c_weight_MCP[i]->cd();
		h_before_MCP[i]->Draw();
		h_after_MCP[i]->SetLineColor(2);
		h_after_MCP[i]->Draw("same");
	
		cout<<"h_before_MCP Integral = "<<h_before_MCP[i]->Integral()<<" , h_after_MCP Integral = "<<h_after_MCP[i]->Integral()<<endl;

		cout<<"h_before_MCP Mean = "<<h_before_MCP[i]->GetMean()<<" , h_after_MCP Mean = "<<h_after_MCP[i]->GetMean()<<endl;

	}



	// NP part

	TH1D *hGen_pt_MCNP=(TH1D*)f_MCNP->Get("hGen_pt");
	TTree *t_fit_MCNP=(TTree*)f_MCNP->Get("t_fit");
	t_fit_MCNP->SetBranchAddress("Dpt", &Dpt);
	t_fit_MCNP->SetBranchAddress(Form("%s",var_compare.Data()), &WeightVar);

	// TFile *f_MCNP_new=TFile::Open(Form("../rootF/%sMC_phiNonPrompt_fitFile_%sWeight.root",Str_PbPb3.Data(), var_compare.Data()),"recreate");
	f_MCNP_new->cd();
	TTree *t_fit_MCNP_new=(TTree*)t_fit_MCNP->CloneTree(0);

	// t_fit_new->Branch("DdlErrWeight",&DdlErrWeight,"DdlErrWeight/F");	
	t_fit_MCNP_new->Branch(Form("%sWeight",var_compare.Data()),&VarWeight);	


	// loop events to add weight

	// double NorFactor[nbin_varcut];
	// double NorFactor[4];
	// NorFactor[0]=1.78111/1.82891;
	// NorFactor[1]=2.13961/2.17648;
	// NorFactor[2]=1.31111/1.2707;
	// NorFactor[3]=3.52482/3.45746;
	if(isPbPb==0){
	NorFactor[0]=1.963222/2.08734;
	NorFactor[1]=2.96439/3.13352;
	NorFactor[2]=2.82868/2.84464;
	NorFactor[3]=1.08609/1.04122;
	}else if(isPbPb==3){
	NorFactor[0]=0.00524482/0.00622198;
	NorFactor[1]=0.00228868/0.00221477;
	NorFactor[0]=0.00715918/0.00845727;
	NorFactor[1]=0.0025658/0.00246642;
	}
	


	nentries=t_fit_MCNP->GetEntries();
	// nentries=10000;
	for(int i=0; i<nentries;i++){
		t_fit_MCNP->GetEntry(i);
		VarWeight=1;
		if(Dpt>DptLowest && Dpt<DptHighest && WeightVar >WeightVarMin && WeightVar< WeightVarMax){
		for(int j=0; j<nbin_varcut; j++){
			if(Dpt>bins_varcut[j] && Dpt<=bins_varcut[j+1]){
				VarWeight=h_ratio[j]->GetBinContent(h_ratio[j]->FindBin(WeightVar))* NorFactor[j];
				if(VarWeight<=0.001){VarWeight=0.001;}
				if(VarWeight>5){VarWeight=5;}
				// cout<<"DdlErrWeight = "<<DdlErrWeight<<endl;	
				break;
			}
		} // for j<nbin_varcut
		} // if Dpt>2 && <40

		t_fit_MCNP_new->Fill();


	}

	f_MCNP_new->cd();
	t_fit_MCNP_new->Write();
	hGen_pt_MCNP->Write();

// check the total integral after weight

	TH1D *h_before_MCNP[nbin_varcut];
	TH1D *h_after_MCNP[nbin_varcut];
	TCanvas *c_weight_MCNP[nbin_varcut];

/*
	int nbin_test=200;
	double binLow_test=0;
	double binHigh_test=200;
	TString s_test="Ddls";
*/


	for(int i=0; i<nbin_varcut;i++){
		h_before_MCNP[i]=new TH1D(Form("h_before_MCNP%i",i),"",nbin_test,binLow_test,binHigh_test); h_before_MCNP[i]->Sumw2();
		h_after_MCNP[i]=new TH1D(Form("h_after_MCNP%i",i),"",nbin_test,binLow_test,binHigh_test); h_after_MCNP[i]->Sumw2();
	
		t_fit_MCNP_new->Project(Form("h_before_MCNP%i",i),s_test.Data(),Form("TotalWeight*(Dpt>%f && Dpt<%f)",bins_varcut[i],bins_varcut[i+1]));
		t_fit_MCNP_new->Project(Form("h_after_MCNP%i",i),s_test.Data(),Form("%sWeight*TotalWeight*(Dpt>%f && Dpt<%f)",var_compare.Data(),bins_varcut[i],bins_varcut[i+1]));


		c_weight_MCNP[i]=new TCanvas(Form("c_weight_MCNP%i",i),Form("c_weight_MCNP%i",i));
		c_weight_MCNP[i]->cd();
		h_before_MCNP[i]->Draw();
		h_after_MCNP[i]->SetLineColor(2);
		h_after_MCNP[i]->Draw("same");
	
		cout<<"h_before_MCNP Integral = "<<h_before_MCNP[i]->Integral()<<" , h_after_MCNP Integral = "<<h_after_MCNP[i]->Integral()<<endl;

		cout<<"h_before_MCNP Mean = "<<h_before_MCNP[i]->GetMean()<<" , h_after_MCNP Mean = "<<h_after_MCNP[i]->GetMean()<<endl;

	}









	// f_MCP_new->Close();


	// TFile *f_fit=TFile::Open(S_fit);

	// TH1D *h_fr_PromptMC=(TH1D*)f_data->Get("h_fr_PromptMC");

	// TH1D *h_var_MixMC=(TH1D*)f_fit->Get(Form("h_%s_MixMC",var_compare.Data()));


	// TH1D *h_var_Data2sig=(TH1D*)f_fit->Get(Form("h_%s_Data2sig",var_compare.Data()));



/*
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
*/

	return 0;

}

int main(int argc, char*argv[]){

  if(argc==6){
    makeReweight_DdxyzErr(atoi(argv[1]), argv[2], argv[3], atof(argv[4]), atof(argv[5]) );
  }else if(argc==10){
    makeReweight_DdxyzErr(atoi(argv[1]), argv[2], argv[3], atof(argv[4]), atof(argv[5]) , atof(argv[6]) , atof(argv[7]) , atof(argv[8]) , atof(argv[9]));
  }else{
     makeReweight_DdxyzErr();
    cout<<"wrong number of input parameters , need 6\n int isPbPb=0, TString var_compare=DdxyzErr, TString var_cut=Dtrk1Pt , double var_cutLow=0.75, double var_cutHigh=1.25"<<endl;
    return 1;
  }
//int Fit_sideband(int isPbPb=0, TString var_compare="DdxyzErr", TString var_cut="Dtrk1Pt" , double var_cutLow=0.75, double var_cutHigh=1.25){


  return 0;
}

