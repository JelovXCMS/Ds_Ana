#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"

// #include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
// #include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
// #include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"

#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"

#include <TGraphErrors.h>
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
#include <TFitResult.h>
#include <TFitResultPtr.h>

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
#include "RooConstVar.h"

#include "varCompare_para.h"

#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooAbsPdf.h"
#include <iomanip>


using namespace RooStats;


using namespace RooFit;
using namespace std;
double DsDataFitRangeLow =1.91;
double DsDataFitRangeHigh = 2.11;
Int_t  nbin_DmassDraw=50;

double textposx=0.2;
double textposy=0.77;

double shiftY=0.0;
double shiftX=0.3;
double oneshift=0.075;

TCanvas *c_binfit_MC[500];
TCanvas *c_binfit_Data[500];
TCanvas *c_roofit_MC[500] ;
TCanvas *c_roofit_Data[500];
TCanvas *c_roofit_Data_pull[500];
TCanvas *c_roofit_Data_withpull[500];

TCanvas *c_roofit_MCstudy[500];
TCanvas *c_roofit_GoF[500];

TCanvas *c_roofit_Data_Dca[500][30];

int count_c_binfit=0;
int count_c_rootfit=0;

double FloatWidthErr_Norm=0;
double FloatWidthVal_Norm=0;
double FloatWidth_DsMass_Norm=1.965;


double DmassSideBand1Low=1.91;
double DmassSideBand1High=1.93;
double DmassSideBand2Low=2.01;
double DmassSideBand2High=2.03;
double Dmass2SigLow=1.949;
double Dmass2SigHigh=1.989;

Float_t SmearRelFactorArr[]={0.01,0.05,0.1,0.2,0.3};
Float_t SmearAbsFactorArr[]={0.0001,0.001,0.005,0.02,0.03,0.05,0.08,0.1,0.15};
Float_t ScaleErrFactorArr[]={0.7,0.8,0.85,0.9,0.95,0.98,1,1.02,1.05,1.1,1.15,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9};
// Float_t SmearFactorArr[]={0.01,0.05};
const int nSmrRelF=sizeof(SmearRelFactorArr)/sizeof(SmearRelFactorArr[0]);
const int nSmrAbsF=sizeof(SmearAbsFactorArr)/sizeof(SmearAbsFactorArr[0]);
const int nSclErrF=sizeof(ScaleErrFactorArr)/sizeof(ScaleErrFactorArr[0]);
Float_t Ddls_SmrRelF[nSmrRelF];
Float_t DdlErr_SmrRelF[nSmrRelF];
Float_t Ddls_SmrAbsF[nSmrAbsF];
Float_t DdlErr_SmrAbsF[nSmrAbsF];
Float_t Ddls_SclErrF[nSclErrF];
Float_t Ddl_SclErrF[nSclErrF];
Float_t DdlErr_SclErrF[nSclErrF];



int CompAndTune_DataMC_Ddls(int isPbPb=3, TString var_compare="Ddls", TString var_cut="Dpt" , double var_cutLow=6, double var_cutHigh=8, double bins_var_Low=0, double bins_var_High=20, double bins_var_DrawLow=0, double bins_var_DrawHigh=20 , int nRebin=1, TString MC_DefaultWeight="TotalWeight*SMrWt",int doSmear=1, int nSmear=10){

	// InitStyle();

	TString var_compareMC="Ddls_BestScale";

	TString GenWt="";
	if(MC_DefaultWeight=="D0Weight"){
		GenWt="_D0Weight";
	}

	gSystem->Exec("mkdir -p fitout");
	gSystem->Exec("mkdir -p plots/fit");
	gSystem->Exec(Form("mkdir -p plots/%s",var_compare.Data()));

	gStyle->SetOptStat(0);
	TString Str_PbPb="pp";
	if(isPbPb){Str_PbPb="PbPb";}

	TLatex *tltx=new TLatex();

	// TString dataName="../sPlot_varCompare_v2/rootF/pp_fitFile.root";
	// TString mcName_Prompt="../sPlot_varCompare_v2/rootF/ppMC_phiPrompt_fitFile.root";
	// TString mcName_NonPrompt="../sPlot_varCompare_v2/rootF/ppMC_phiNonPrompt_fitFile.root";

	TString dataName="./rootF/pp_fitFile.root";
	// TString mcName_Prompt="./rootF/ppMC_phiPrompt_fitFile.root";
	// TString mcName_NonPrompt="./rootF/ppMC_phiNonPrompt_fitFile.root";
	TString mcName_Prompt="./rootF/ppMC_phiPrompt_fitFile_DdlsScl.root";
	TString mcName_NonPrompt="./rootF/ppMC_phiNonPrompt_fitFile_DdlsScl.root";
	double Dpt_Low=Dpt_Low_pp;
	double Dpt_High=Dpt_Hight_pp;

	double Dalpha_cut=Dalpha_cut_pp;
	double Dchi2cl_cut=Dchi2cl_cut_pp;
	double Ddls_cut=Ddls_cut_pp;

	double *bins_var=bins_DdxyzErr_pp;
	int nbin_var=nbin_DdxyzErr_pp;
	// double bins_var_Low=0;
	// double bins_var_High=1;
	nbin_var=1000;
	nbin_var=nbin_var/nRebin;

	TString s_PbPb="pp";

	if(isPbPb){
		Dpt_Low=Dpt_Low_PbPb;
		Dpt_High=Dpt_Hight_PbPb;

		//dataName="../sPlot_varCompare_v2/rootF/PbPb3_fitFile.root";
		// mcName_Prompt="../sPlot_varCompare_v2/rootF/PbPb3MC_phiPrompt_fitFile.root";
		// mcName_NonPrompt="../sPlot_varCompare_v2/rootF/PbPb3MC_phiNonPrompt_fitFile.root";
	  dataName="./rootF/PbPb3_fitFile.root";
	  mcName_Prompt="./rootF/PbPb3MC_phiPrompt_fitFile_DdlsScl.root";
	  mcName_NonPrompt="./rootF/PbPb3MC_phiNonPrompt_fitFile_DdlsScl.root";

		Dalpha_cut=Dalpha_cut_PbPb3;
		Dchi2cl_cut=Dchi2cl_cut_PbPb3;
		Ddls_cut=Ddls_cut_PbPb3;
		s_PbPb="PbPb";

		if(var_compare=="Dchi2cl"){
			Dchi2cl_cut=0.05;
		}

	}

	if(var_compare=="Dchi2cl"){
		Dchi2cl_cut=0.02;
	}

	if(var_cut=="Dpt"){
		Dpt_Low=var_cutLow;
		Dpt_High=var_cutHigh;
	}

	// set Dalpha_cut=0.12 for Ddls
	Dalpha_cut=0.12;


// tighter cut for PbPb low pt
 if(!isPbPb && Dpt_High<=4){
    Dchi2cl_cut=0.15;
    Ddls_cut=3.5;
  }else if(!isPbPb && Dpt_High<=6 && Dpt_Low>4){
    Dchi2cl_cut=0.1;
    Ddls_cut=2.9;
	}else if(isPbPb && Dpt_High<=10){
		Dchi2cl_cut=0.3;
		Ddls_cut=5.0;
	}



	TString sPlotDataName=Form("./fitout/%s_Dpt%.0fto%.0f.root",Str_PbPb.Data(),Dpt_Low*100,Dpt_High*100);
// pp_Dpt600to800.root
	cout<<"sPlotDataName = "<<sPlotDataName<<endl;




	TString DataCuts=Form("Dpt>%f && Dpt<%f && Dalpha < %f && Dchi2cl > %f && Ddls >%f && %s > %f && %s <%f", Dpt_Low, Dpt_High, Dalpha_cut, Dchi2cl_cut, Ddls_cut, var_cut.Data(), var_cutLow,var_cut.Data(), var_cutHigh);
	TString DataCutsMC=Form("Dpt>%f && Dpt<%f && Dalpha < %f && Dchi2cl > %f && Ddls_BestScale >%f && %s > %f && %s <%f", Dpt_Low, Dpt_High, Dalpha_cut, Dchi2cl_cut, Ddls_cut, var_cut.Data(), var_cutLow,var_cut.Data(), var_cutHigh);
	// TString DataCuts=Form("Dalpha < %f && Dchi2cl > %f && Ddls >%f && %s > %f && %s <%f", Dalpha_cut, Dchi2cl_cut, Ddls_cut, var_cut.Data(), var_cutLow,var_cut.Data(), var_cutHigh);

	TFile *f_data=TFile::Open(dataName.Data());
	TFile *f_mc_Prompt=TFile::Open(mcName_Prompt.Data());
	TFile *f_mc_NonPrompt=TFile::Open(mcName_NonPrompt.Data());

	TTree *t_Ds_MCPrompt=(TTree*)f_mc_Prompt->Get(Form("t_fit"));
	TTree *t_Ds_MCNonPrompt=(TTree*)f_mc_NonPrompt->Get(Form("t_fit"));
	TTree *t_Ds_Data=(TTree*)f_data->Get(Form("t_fit"));

	TFile *f_sPlot=TFile::Open(sPlotDataName);
	if(!f_sPlot){
		cout<<"fit sPlot file : "<<sPlotDataName<<" not exist, terminate"<<endl;
		return 2;
	}
	TTree *t_sPlot=(TTree*)f_sPlot->Get("tree_sPlot");
	TCanvas *c_Ddls=new TCanvas("c_Ddls");
	c_Ddls->cd();
	// t_sPlot->Draw("Dalpha","NumSig_sw");
	// t_sPlot->Draw("Dmass","NumSig_sw");
	// t_sPlot->Draw("Dmass","NumBkg_sw");

	gSystem->Exec("mkdir -p Compare_DataMC_out_Ddls");
	TFile *f_out=TFile::Open(Form("./Compare_DataMC_out_Ddls/%s_%s%.0fto%.0f.root",Str_PbPb.Data(), var_cut.Data(), var_cutLow*100, var_cutHigh*100) ,"recreate");

	TH1D *h_DataRawYield_test=(TH1D*)f_sPlot->Get("h_DataRawYield");
	TH1D *h_Data2sig_SigFrac_test=(TH1D*)f_sPlot->Get("h_Data2sig_SigFrac");

	h_DataRawYield_test->Write();
	h_Data2sig_SigFrac_test->Write();

	double Yield = h_DataRawYield_test->GetBinContent(1);
	double SigFrac = h_Data2sig_SigFrac_test->GetBinContent(1);

	cout<<"Yield = "<<Yield<<" , SigFrac = "<<SigFrac<<endl;


	gSystem->Exec("mkdir -p FitSideBand_out");


	// compare sideband with sPlot
	
	int nbin=200;
	double binLow=Ddls_cut;
	double binHigh=30;

	double binstep=0.1;
	double cur=binLow;
	int nbincnt=0;
	vector<double> binsV;
	while(cur<=binHigh+0.2){
		cout<<"cur  = "<<cur<<endl;
		binsV.emplace_back(cur);
		cur+=binstep;
		nbincnt++;
		if(cur>7){ binstep=0.2;}
		if(cur>10){ binstep=0.4;}
		if(cur>12){ binstep=0.6;}
		if(cur>15){ binstep=1;}
	}

cout<<"nbincnt = "<<nbincnt<<endl;	
	if(binsV.size()%2==0){
		binsV.emplace_back(cur+binstep);
	}


		nbin=binsV.size()-1;
	double binsArr[binsV.size()];
	for(unsigned int i=0; i<binsV.size(); i++){
			binsArr[i]=binsV[i];
	}


	TH1D *h_Ddls_sPlot=new TH1D("h_Ddls_sPlot","h_Ddls_sPlot",nbin,binLow,binHigh); h_Ddls_sPlot->Sumw2();
	// TH1D *h_Ddls_sPlot=new TH1D("h_Ddls_sPlot","h_Ddls_sPlot",nbin,binsArr); h_Ddls_sPlot->Sumw2();
	t_sPlot->Draw("Ddls>>h_Ddls_sPlot",Form("NumSig_sw*(%s)",DataCuts.Data()));


	// return 1;

	// TString DataCuts=Form("Dpt>%f && Dpt<%f && Dalpha < %f && Dchi2cl > %f && Ddls >%f && %s > %f && %s <%f", Dpt_Low, Dpt_High, Dalpha_cut, Dchi2cl_cut, Ddls_cut, var_cut.Data(), var_cutLow,var_cut.Data(), var_cutHigh);

	h_Ddls_sPlot->Scale(1/h_Ddls_sPlot->Integral());

	f_out->cd();
	h_Ddls_sPlot->Write();

	c_Ddls->cd();

	h_Ddls_sPlot->Draw();
	// c_Ddls->BuildLegend();

	// build Mix MC
	TFile *f_NonPromptDs=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BtoDs_Results/output/BtoDs.root");
	// TH1D *hBtoDsCrossSectionPP_AnaBin=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin");
	TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight");

	// TH1D *hBtoDsdNdPtPbPb_AnaBin=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin");
	TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight");

	TH1D *hBtoDs_AnaBin_pythiaWeight=hBtoDsCrossSectionPP_AnaBin_pythiaWeight;  // this is differential cross section
	if(isPbPb){ hBtoDs_AnaBin_pythiaWeight=hBtoDsdNdPtPbPb_AnaBin_pythiaWeight; }

	MutiplyBinWidth(hBtoDs_AnaBin_pythiaWeight); 
	double CS_Integral=hBtoDs_AnaBin_pythiaWeight->Integral(hBtoDs_AnaBin_pythiaWeight->GetXaxis()->FindBin(Dpt_Low+0.00001), hBtoDs_AnaBin_pythiaWeight->GetXaxis()->FindBin(Dpt_High-0.00001) );	

	cout<<"CS_Integral = "<<CS_Integral<<endl;


			// prompt fraction cal
			double fr_withcut_prompt=0.85; // just for get quick result

			// extract fr_prompt from fit and BtoDs estimation
			bool use_fr_fromBtoD=true;
			if(isPbPb){use_fr_fromBtoD=false;}

			if(use_fr_fromBtoD){

				double N_yield=Yield;
				double LumiNevt=LumiSum;
				if(isPbPb){LumiNevt=NevtPbPb3;}

				// nonprompt eff
				TH1D *h_eff_nonprompt_phi=new TH1D("h_eff_nonprompt_phi","h_eff_nonprompt_phi",40,0,40);	h_eff_nonprompt_phi->Sumw2();
				t_Ds_MCNonPrompt->Project("h_eff_nonprompt_phi","Dpt",(TCut)Form("%s*(%s)",MC_DefaultWeight.Data(),DataCutsMC.Data()));
				double NonPrompt_phi_reco=h_eff_nonprompt_phi->Integral(h_eff_nonprompt_phi->FindBin(Dpt_Low+0.00001), h_eff_nonprompt_phi->FindBin(Dpt_High-0.00001));
				TH1D *hGen_pt_nonprompt_phi=(TH1D*)f_mc_NonPrompt->Get(Form("hGen_pt%s",GenWt.Data() ) );
				double NonPrompt_phi_Gen=hGen_pt_nonprompt_phi->Integral(hGen_pt_nonprompt_phi->FindBin(Dpt_Low+0.00001),hGen_pt_nonprompt_phi->FindBin(Dpt_High-0.0001) );

				double NonPrompt_phi_Eff=NonPrompt_phi_reco/NonPrompt_phi_Gen;
			
				cout<<"NonPrompt_phi_Eff = "<<NonPrompt_phi_Eff<<endl;
				//		cout<<"NonPrompt_f0_Eff = "<<NonPrompt_f0_Eff<<endl;

				// prompt eff , only for comparison
				TH1D *h_eff_prompt_phi=new TH1D("h_eff_prompt_phi","h_eff_prompt_phi",40,0,40);	h_eff_prompt_phi->Sumw2();
				t_Ds_MCPrompt->Project("h_eff_prompt_phi","Dpt",(TCut)Form("%s*(%s)",MC_DefaultWeight.Data(),DataCutsMC.Data()));
				double Prompt_phi_reco=h_eff_prompt_phi->Integral(h_eff_prompt_phi->FindBin(Dpt_Low+0.00001), h_eff_prompt_phi->FindBin(Dpt_High-0.00001));
				TH1D *hGen_pt_prompt_phi=(TH1D*)f_mc_Prompt->Get(Form("hGen_pt%s",GenWt.Data() ) );
				double Prompt_phi_Gen=hGen_pt_prompt_phi->Integral(hGen_pt_prompt_phi->FindBin(Dpt_Low+0.00001),hGen_pt_prompt_phi->FindBin(Dpt_High-0.0001) );

				double Prompt_phi_Eff=Prompt_phi_reco/Prompt_phi_Gen;
			
				cout<<"Prompt_phi_Eff = "<<Prompt_phi_Eff<<endl;


				 double phiFr=0.96 ;
				if(isPbPb){
					phiFr=0.9;
				}

				double N_NonPrompt_yield=CS_Integral*(2*LumiNevt)*( BRphi*NonPrompt_phi_Eff)/phiFr;
				// double N_NonPrompt_yield=CS_Integral*(2*LumiNevt)*( BRphi*NonPrompt_phi_Eff+ BRf0*NonPrompt_f0_Eff );


				fr_withcut_prompt=1-(N_NonPrompt_yield/N_yield);

				cout<<"fr_withcut_prompt = "<<fr_withcut_prompt<<" , N_NonPrompt_yield = "<<N_NonPrompt_yield<<" , N_yield = "<<N_yield<<endl;

				delete h_eff_nonprompt_phi;


			} // end fr_prompt calculation



	// Draw MC plots

  TH1D *h_var_PromptMC =new TH1D(Form("h_%s_PromptMC",var_compare.Data()), Form("h_%s_PromptMC",var_compare.Data()),nbin,binLow,binHigh); h_var_PromptMC->Sumw2();
  // TH1D *h_var_PromptMC =new TH1D(Form("h_%s_PromptMC",var_compare.Data()), Form("h_%s_PromptMC",var_compare.Data()),nbin,binsArr); h_var_PromptMC->Sumw2();
 
 t_Ds_MCPrompt->Project(Form("h_%s_PromptMC",var_compare.Data()),var_compareMC.Data(),Form("(Dmass>%f && Dmass<%f && %s)*%s",Dmass2SigLow, Dmass2SigHigh, DataCutsMC.Data(),MC_DefaultWeight.Data() ));

  h_var_PromptMC->Scale(1/h_var_PromptMC->Integral());

  TH1D *h_var_NonPromptMC =new TH1D(Form("h_%s_NonPromptMC",var_compare.Data()), Form("h_%s_NonPromptMC",var_compare.Data()),nbin,binLow,binHigh); h_var_NonPromptMC->Sumw2();
  // TH1D *h_var_NonPromptMC =new TH1D(Form("h_%s_NonPromptMC",var_compare.Data()), Form("h_%s_NonPromptMC",var_compare.Data()),nbin,binsArr); h_var_NonPromptMC->Sumw2();
  t_Ds_MCNonPrompt->Project(Form("h_%s_NonPromptMC",var_compare.Data()),var_compareMC.Data(),Form("(Dmass>%f && Dmass<%f && %s)*%s", Dmass2SigLow, Dmass2SigHigh, DataCutsMC.Data(),MC_DefaultWeight.Data()));

  h_var_NonPromptMC->Scale(1/h_var_NonPromptMC->Integral());

	TCanvas *c_MCPNP=new TCanvas();
	c_MCPNP->cd();
	h_var_PromptMC->Draw();
	h_var_NonPromptMC->SetLineColor(2);
	h_var_NonPromptMC->Draw("same");



	TH1D *h_var_MixMC= new TH1D(Form("h_%s_MixMC",var_compare.Data()),"h_var_MixMC",nbin,binLow,binHigh); h_var_MixMC->Sumw2();
	// TH1D *h_var_MixMC= new TH1D(Form("h_%s_MixMC",var_compare.Data()),"h_var_MixMC",nbin,binsArr); h_var_MixMC->Sumw2();

	h_var_MixMC->Add(h_var_PromptMC,h_var_NonPromptMC,fr_withcut_prompt,1-fr_withcut_prompt);

	TH1D *h_var_MixMC_fit=(TH1D*)h_var_MixMC->Clone("h_var_MixMC_fit");


	h_var_MixMC->Draw("same");

	gSystem->Exec("mkdir -p plots_DataMC_Ddls");
	c_MCPNP->SaveAs(Form("./plots_DataMC_Ddls/MCMixPNP_%s_%s_%s%.0fto%.0f.png",var_compare.Data(), Str_PbPb.Data(), var_cut.Data(), var_cutLow*100, var_cutHigh*100));

	cout<<"start drawing"<<endl;


	TH1D *h_fr_PromptMC=new TH1D("h_fr_PromptMC","",1,0,1);
	h_fr_PromptMC->SetBinContent(1,fr_withcut_prompt); 

	// save histogram before setRange

	f_out->cd();
	// h_var_sideband->Write();
	// h_var_Data2sig->Write();
	// h_var_Data2sigMix->Write();
	h_var_PromptMC->Write();
	h_var_NonPromptMC->Write();
	h_var_MixMC->Write();

	// h_Data2sig_SigFrac->Write();
	h_fr_PromptMC->Write();

	TCanvas *c_var_compare= new TCanvas("c_var_compare","c_var_compare",800,800);
	c_var_compare->cd();

	TPad *pad1 = new TPad("pad1","top pad", 0.0,0.25,1.0,1.0);
	pad1->SetBottomMargin(0.0);
	pad1->Draw();
	TPad *pad2 = new TPad("pad2","bottom pad",0.0,0.0,1.0,0.25);
	pad2->SetTopMargin(0.0);
	pad2->SetBottomMargin(0.30);
	pad2->Draw();

	pad1->cd();

	int rebinN=4;
	if(var_compare=="Dchi2cl"){
		rebinN=20;
	}

	
	h_Ddls_sPlot->Rebin(rebinN);
	h_Ddls_sPlot->Smooth(3);
	h_var_PromptMC->Rebin(rebinN);
	h_var_NonPromptMC->Rebin(rebinN);
	h_var_MixMC->Rebin(rebinN);

	TH1D *h_Ddls_sPlot_divBinWid=(TH1D*)h_Ddls_sPlot->Clone("h_Ddls_sPlot_divBinWid");
	TH1D *h_var_PromptMC_divBinWid=(TH1D*)h_var_PromptMC->Clone("h_var_PromptMC_divBinWid");
	TH1D *h_var_NonPromptMC_divBinWid=(TH1D*)h_var_NonPromptMC->Clone("h_var_NonPromptMC_divBinWid");
	TH1D *h_var_MixMC_divBinWid=(TH1D*)h_var_MixMC->Clone("h_var_MixMC_divBinWid");

	divideBinWidth(h_Ddls_sPlot_divBinWid);
	divideBinWidth(h_var_PromptMC_divBinWid);
	divideBinWidth(h_var_NonPromptMC_divBinWid);
	divideBinWidth(h_var_MixMC_divBinWid);

	double h_max=h_Ddls_sPlot_divBinWid->GetMaximum();
	double h_min=h_Ddls_sPlot_divBinWid->GetMinimum();

	h_Ddls_sPlot_divBinWid->SetMaximum(h_max +  0.4*(h_max-h_min ) );
	h_Ddls_sPlot_divBinWid->SetMinimum(h_min -  0.1*(h_max-h_min ) );

	// h_Ddls_sPlot->GetXaxis()->SetRangeUser(bins_var_DrawLow,bins_var_DrawHigh); // normalized before setRange

	h_Ddls_sPlot_divBinWid->GetXaxis()->SetRangeUser(binLow,binHigh); // normalized before setRange
	h_Ddls_sPlot_divBinWid->GetXaxis()->SetTitle(var_compare.Data());
	h_Ddls_sPlot_divBinWid->SetTitle("");
	h_Ddls_sPlot_divBinWid->SetLineColor(1);
	h_Ddls_sPlot_divBinWid->SetMarkerColor(1);
	h_Ddls_sPlot_divBinWid->SetMarkerStyle(22);
	h_Ddls_sPlot_divBinWid->Draw();
	h_var_PromptMC_divBinWid->SetLineColor(2);
	h_var_PromptMC_divBinWid->SetMarkerColor(2);
	h_var_PromptMC_divBinWid->SetMarkerStyle(26);
	h_var_PromptMC_divBinWid->Draw("same");
	h_var_NonPromptMC_divBinWid->SetLineColor(4);
	h_var_NonPromptMC_divBinWid->SetMarkerColor(4);
	h_var_NonPromptMC_divBinWid->SetMarkerStyle(26);
	h_var_NonPromptMC_divBinWid->Draw("same");
	h_var_MixMC_divBinWid->SetLineColor(kMagenta+3);
	h_var_MixMC_divBinWid->SetMarkerColor(kMagenta+3);
	h_var_MixMC_divBinWid->SetMarkerStyle(26);
	h_var_MixMC_divBinWid->Draw("same");

	gPad->SetLogx();

	TLegend *le_var = new TLegend(0.55,0.45,0.88,0.73);
	le_var->SetBorderSize(0);
	le_var->SetTextSize(0.04);
	le_var->AddEntry(h_Ddls_sPlot_divBinWid,"Data","pl");
	le_var->AddEntry(h_var_PromptMC_divBinWid,"MC Prompt D_{S}","pl");
	le_var->AddEntry(h_var_NonPromptMC_divBinWid,"MC NonPrompt D_{S}","pl");
	le_var->AddEntry(h_var_MixMC_divBinWid,"MC Mix D_{S}","pl");
	le_var->Draw("same");

	shiftY=0.07;
	shiftX=0.35;
	tltx->SetTextSize(0.04);
	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Ds %s ",Str_PbPb.Data(), var_compare.Data())); shiftY-=oneshift;
	// tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s Ds %s %.0f<p_{T}<%.0f",Str_PbPb.Data(), var_compare.Data(),Dpt_Low,Dpt_High)); shiftY-=oneshift;
	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Var : %s %.2f to %.2f ", var_cut.Data(), var_cutLow, var_cutHigh)); shiftY-=oneshift;
	// tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Prompt fraction : %.2f ", fr_withcut_prompt)); shiftY-=oneshift;

	pad2->cd();


	TH1D *h_ratio=(TH1D*)h_Ddls_sPlot->Clone("h_ratio");
	h_ratio->Divide(h_var_MixMC);
	h_ratio->SetTitle("");
	h_ratio->GetYaxis()->SetTitle("Data/MC");
	h_ratio->GetXaxis()->SetTitle(var_compare.Data());
	h_ratio->GetYaxis()->CenterTitle();
	h_ratio->GetYaxis()->SetLabelSize(0.05);
	h_ratio->GetYaxis()->SetTitleSize(0.08);
	h_ratio->GetYaxis()->SetTitleOffset(0.38);
	h_ratio->GetXaxis()->SetLabelSize(0.12);
	h_ratio->GetXaxis()->SetTitleSize(0.12);
	h_ratio->SetMaximum(2.1);
	h_ratio->SetMinimum(0.01);
	h_ratio->Draw();

	gPad->SetLogx();


	gSystem->Exec("mkdir -p plots_sPlotWeight_Ddls");
	// c_var_compare->SaveAs(Form("./plots_DataMC/MCDataCompare_%s_%s_%s%.0fto%.0f.png",var_compare.Data(), Str_PbPb.Data(), var_cut.Data(), var_cutLow*100, var_cutHigh*100));
	c_var_compare->SaveAs(Form("./plots_sPlotWeight_Ddls/MCDataCompare_%s_%s_%s%.0fto%.0f.png",var_compare.Data(), Str_PbPb.Data(), var_cut.Data(), var_cutLow*100, var_cutHigh*100));
	// c_var_compare->SaveAs(Form("./plots/%s/%s_%s%.0fto%.0f.png", var_compare.Data(),Str_PbPb.Data(), var_cut.Data(),var_cutLow*100,var_cutHigh*100));


	TH1D *h_var_ratioSmooth=(TH1D*)h_ratio->Clone(Form("h_%s_ratioSmooth",var_compare.Data()));
	// int RebinN=2;
	int SmoothN=2;
	// h_var_ratioSmooth->Rebin(RebinN);
	// h_var_MixMC->Rebin(RebinN);
	// h_var_ratioSmooth->Divide(h_var_MixMC);
	h_var_ratioSmooth->Smooth(SmoothN);

	h_var_ratioSmooth->Write();


	TCanvas *c_weight=new TCanvas("c_weight","c_weight");
	c_weight->cd();
	// h_var_ratioSmooth->Draw();


	cout<<"Data mean = "<<h_Ddls_sPlot->GetMean()<<" , RMS = "<<h_Ddls_sPlot->GetRMS()<<" , StdDev = "<<h_Ddls_sPlot->GetStdDev()<<endl;
	cout<<"MC mean = "<<h_var_MixMC->GetMean()<<" , RMS = "<<h_var_MixMC->GetRMS()<<" , StdDev = "<<h_var_MixMC->GetStdDev()<<endl;

// Using TGraph for weight so we can use interpolation.

	TGraphErrors *gr_weight=new TGraphErrors(h_var_ratioSmooth);
	gr_weight->SetMarkerColor(2);
	// gr_weight->Draw("ap");
	// h_var_ratioSmooth->Draw("same");
	gr_weight->Draw("ap");
	
	c_weight->SaveAs(Form("./plots_DataMC_Ddls/MCWeight_%s_%s_%s%.0fto%.0f.png",var_compare.Data(), Str_PbPb.Data(), var_cut.Data(), var_cutLow*100, var_cutHigh*100));

	f_out->cd();
	gr_weight->Write("gr_weight");

// check Prompt and Nonprompt total integral after apply weight, need to scale the weight to 1

  // TH1D *h_var_PromptMC =new TH1D(Form("h_%s_PromptMC",var_compare.Data()), Form("h_%s_PromptMC",var_compare.Data()),nbin,binLow,binHigh); h_var_PromptMC->Sumw2();
//  t_Ds_MCPrompt->Project(Form("h_%s_PromptMC",var_compare.Data()),var_compare.Data(),Form("(Dmass>%f && Dmass<%f && %s)*%s",Dmass2SigLow, Dmass2SigHigh, DataCuts.Data(),MC_DefaultWeight.Data() ));
	// TString DataCuts=Form("Dpt>%f && Dpt<%f && Dalpha < %f && Dchi2cl > %f && Ddls >%f && %s > %f && %s <%f", Dpt_Low, Dpt_High, Dalpha_cut, Dchi2cl_cut, Ddls_cut, var_cut.Data(), var_cutLow,var_cut.Data(), var_cutHigh);


/*
	// bool doFit=true;
	// if(doFit){

	TH1D *h_Ddls_sPlot_fit=new TH1D("h_Ddls_sPlot_fit","h_Ddls_sPlot_fit",nbin,binLow,binHigh); h_Ddls_sPlot_fit->Sumw2();
	t_sPlot->Draw("Ddls>>h_Ddls_sPlot_fit","NumSig_sw");

	h_Ddls_sPlot_fit->Rebin(2);
	h_var_MixMC_fit->Rebin(2);

	h_Ddls_sPlot_fit->Scale(1/h_Ddls_sPlot_fit->Integral());
	h_var_MixMC_fit->Scale(1/h_Ddls_sPlot_fit->Integral());

	TCanvas *c_f1Fit=new TCanvas("c_f1Fit");
	c_f1Fit->cd();

	TF1 *f1_Ddls =new TF1("f1_Ddls","[0]*pow(x,[1])*exp([2]*x)+ [3]*pow(x,[4])*exp([5]*x)"); 
	f1_Ddls->SetRange(0.0,0.2) ;

	if(Dpt_High<=4){
		f1_Ddls->FixParameter(3,0);
		f1_Ddls->FixParameter(4,0);
		f1_Ddls->FixParameter(5,0);
	}	

	h_Ddls_sPlot_fit->Fit("f1_Ddls","QN0","",0.0,0.2);
	h_Ddls_sPlot_fit->Fit("f1_Ddls","QN0","",0.0,0.2);
	h_Ddls_sPlot_fit->Fit("f1_Ddls","QN0","",0.0,0.2);
	h_Ddls_sPlot_fit->Fit("f1_Ddls","0","",0.0,0.2);
	h_Ddls_sPlot_fit->Fit("f1_Ddls","0","",0.0,0.2);
	h_Ddls_sPlot_fit->Fit("f1_Ddls","0","",0.0,0.2);
	h_Ddls_sPlot_fit->Fit("f1_Ddls","0","",0.0,0.2);
	// h_Ddls_sPlot_fit->Fit("f1_Ddls","0ILM","",0.0,0.2);
	// h_Ddls_sPlot_fit->Fit("f1_Ddls","0ILM","",0.0,0.2);
	// h_Ddls_sPlot_fit->Fit("f1_Ddls","0ILM","",0.0,0.2);

	double f1_data_p0=f1_Ddls->GetParameter(0);
	double f1_data_p1=f1_Ddls->GetParameter(1);
	double f1_data_p2=f1_Ddls->GetParameter(2);
	double f1_data_p3=f1_Ddls->GetParameter(3);
	double f1_data_p4=f1_Ddls->GetParameter(4);
	double f1_data_p5=f1_Ddls->GetParameter(5);


	h_Ddls_sPlot_fit->Draw();
	f1_Ddls->Draw("same");


	TCanvas *c_f1Fit_MC=new TCanvas("c_f1Fit_MC");
	c_f1Fit_MC->cd();
	
	// TF1 *f1_Ddls_MC =new TF1("f1_Ddls_MC","[0]*pow(x,[1])*exp([2]*x)");
	TF1 *f1_Ddls_MC =new TF1("f1_Ddls_MC","[0]*pow(x,[1])*exp([2]*x) + [3]*pow(x,[4])*exp([5]*x)");
	f1_Ddls_MC->SetRange(0.0,0.2) ;
	h_var_MixMC_fit->Fit("f1_Ddls_MC","QN0","",0.0,0.2);
	h_var_MixMC_fit->Fit("f1_Ddls_MC","QN0","",0.0,0.2);
	h_var_MixMC_fit->Fit("f1_Ddls_MC","QN0","",0.0,0.2);
	h_var_MixMC_fit->Fit("f1_Ddls_MC","0","",0.0,0.2);
	h_var_MixMC_fit->Fit("f1_Ddls_MC","0","",0.0,0.2);
	h_var_MixMC_fit->Fit("f1_Ddls_MC","0","",0.0,0.2);
  // h_var_MixMC_fit->Fit("f1_Ddls_MC","0ILM","",0.0,0.2);
	// h_var_MixMC_fit->Fit("f1_Ddls_MC","0ILM","",0.0,0.2);
	// h_var_MixMC_fit->Fit("f1_Ddls_MC","0ILM","",0.0,0.2);
	


	
	h_var_MixMC_fit->Draw();
	f1_Ddls_MC->Draw("same");

	double f1_MC_p0=f1_Ddls_MC->GetParameter(0);
	double f1_MC_p1=f1_Ddls_MC->GetParameter(1);
	double f1_MC_p2=f1_Ddls_MC->GetParameter(2);
	double f1_MC_p3=f1_Ddls_MC->GetParameter(3);
	double f1_MC_p4=f1_Ddls_MC->GetParameter(4);
	double f1_MC_p5=f1_Ddls_MC->GetParameter(5);



	TCanvas *c_fitRatio=new TCanvas("c_fitRatio");
	c_fitRatio->cd();

	f1_Ddls->Write();
	f1_Ddls_MC->Write();


		
	// TF1 *f1_DdlsRatio =new TF1("f1_DdlsRatio","([0]*pow(x,[1])*exp([2]*x)  )/ ( [3]*pow(x,[4])*exp([5]*x) )"); 
	TF1 *f1_DdlsRatio =new TF1("f1_DdlsRatio","( [0]*pow(x,[1])*exp([2]*x) + [3]*pow(x,[4])*exp([5]*x)  )/ ( [6]*pow(x,[7])*exp([8]*x)+ [9]*pow(x,[10])*exp([11]*x) )"); 
	f1_DdlsRatio->SetRange(0.0,0.2);
	// f1_DdlsRatio->SetParameters(f1_data_p0,f1_data_p1,f1_data_p2,f1_data_p3,f1_data_p4,f1_data_p5,f1_MC_p0,f1_MC_p1,f1_MC_p2,f1_MC_p3,f1_MC_p4,f1_MC_p5);
	double pars[12]={f1_data_p0,f1_data_p1,f1_data_p2,f1_data_p3,f1_data_p4,f1_data_p5,f1_MC_p0,f1_MC_p1,f1_MC_p2,f1_MC_p3,f1_MC_p4,f1_MC_p5};
	f1_DdlsRatio->SetParameters(pars);
	f1_DdlsRatio->Draw();

	
	f1_DdlsRatio->Write();
	f1_DdlsRatio->SetTitle("");
	cout<<"f1_data_p0 = "<<f1_data_p0<<endl;

	TCanvas *c_FitWRatio=new TCanvas("c_FitWRatio","c_FitWRatio",800,800);
	c_FitWRatio->cd();
	
  TPad *pad1a = new TPad("pad1a","top pad",0.0,0.25,1.0,1.0);
  pad1a->SetTopMargin(0.08);
  pad1a->SetBottomMargin(0.0);
  // pad1a->SetRightMargin(c_Rmg);
  // pad1a->SetLeftMargin(c_Lmg);
  pad1a->Draw();
  TPad *pad2a = new TPad("pad2a","bottom pad",0.0,0.00,1.0,0.25);
  pad2a->SetTopMargin(0.0);
  pad2a->SetBottomMargin(0.30);
  // pad2a->SetRightMargin(c_Rmg);
  // pad2a->SetLeftMargin(c_Lmg);
  pad2a->Draw();

	pad1a->cd();


	h_Ddls_sPlot_fit->SetTitle("");
	h_Ddls_sPlot_fit->SetLineColor(1);
	f1_Ddls->SetLineColor(1);
	h_Ddls_sPlot_fit->Draw();
	f1_Ddls->Draw("same");

	h_var_MixMC_fit->SetLineColor(2);
	f1_Ddls_MC->SetLineColor(2);
	h_var_MixMC_fit->Draw("same");
	f1_Ddls_MC->Draw("same");

	TLegend *le_fit=new TLegend(0.6,0.6,0.85,0.85);
	le_fit->SetBorderSize(0);
	le_fit->AddEntry((TObject*)0,Form("%s",Str_PbPb.Data()),"" );
	le_fit->AddEntry((TObject*)0,Form("%.0f < D_{S} p_{T} < %.0f",Dpt_Low,Dpt_High ),"" );
	le_fit->AddEntry(h_Ddls_sPlot_fit,"Data","l");
	le_fit->AddEntry(h_var_MixMC_fit,"MC","l");
	le_fit->Draw("same");

	pad2a->cd();
	f1_DdlsRatio->Draw();
	f1_DdlsRatio->SetMaximum(2.5);
	f1_DdlsRatio->SetMinimum(0);
	f1_DdlsRatio->GetXaxis()->SetTitle("Decay Length Significance");
	f1_DdlsRatio->GetYaxis()->SetTitle("Data/MC");
	f1_DdlsRatio->GetYaxis()->SetTitleSize(0.10);
	f1_DdlsRatio->GetYaxis()->CenterTitle();
	f1_DdlsRatio->GetXaxis()->CenterTitle();
	f1_DdlsRatio->GetYaxis()->SetTitleOffset(0.4);
	f1_DdlsRatio->GetYaxis()->SetLabelSize(0.08);
	f1_DdlsRatio->GetXaxis()->SetTitleSize(0.12);
	f1_DdlsRatio->GetXaxis()->SetTitleOffset(1.0);
	f1_DdlsRatio->GetXaxis()->SetLabelSize(0.12);

	// c_FitWRatio->SaveAs();	
	// c_FitWRatio->SaveAs(Form("./plots_DataMC/MCFitWeight_%s_%s_%s%.0fto%.0f.png",var_compare.Data(), Str_PbPb.Data(), var_cut.Data(), var_cutLow*100, var_cutHigh*100));

	c_FitWRatio->SaveAs(Form("./plots_sPlotWeight/MCDataFitCompare_%s_%s_%s%.0fto%.0f.png",var_compare.Data(), Str_PbPb.Data(), var_cut.Data(), var_cutLow*100, var_cutHigh*100));
	c_FitWRatio->SaveAs(Form("./plots_sPlotWeight/MCDataFitCompare_%s_%s_%s%.0fto%.0f.pdf",var_compare.Data(), Str_PbPb.Data(), var_cut.Data(), var_cutLow*100, var_cutHigh*100));

	TGraph *gr_FitWeight=new TGraph();
	double begin=0.0001;
	double end=0.1999;
	double steps=0.0001;
	double current=begin;
	int istep=0;
	while(current<end){
		double FitWeight=f1_DdlsRatio->Eval(current);
		if(FitWeight<0){FitWeight=0;}
		// cout<<"x = "<<current<<" , weight = "<<FitWeight<<endl;
		gr_FitWeight->SetPoint(istep,current,FitWeight);
		istep++;
		current+=steps;
	}	


	gr_FitWeight->Write("gr_Fitweight");



	bool doRoofit=false;
	if(doRoofit){
	RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

	RooRealVar x("x","Dalpha",0.0,0.2);
	RooDataHist dh_Data("dh_Data","dh_Data",x,Import(*h_Ddls_sPlot_fit));
	RooPlot *Dalpha_frame=x.frame();


	RooRealVar p0("p0","p0",0,100);
	RooRealVar p1("p1","p1",-100,100);
	RooRealVar p2("p2","p2",-1000,1000);

	RooGenericPdf gPdf("gPdf","p0*pow(x,p1)*exp(x*p2)",RooArgSet(x,p0,p1,p2));


	// p0.setVal(3.8);
	// p1.setVal(0.82);
	// p2.setVal(-0.61);

	gPdf.fitTo(dh_Data,NumCPU(30));
	gPdf.fitTo(dh_Data,NumCPU(30));
	gPdf.fitTo(dh_Data,NumCPU(30));
	gPdf.fitTo(dh_Data,NumCPU(30));
	gPdf.fitTo(dh_Data,NumCPU(30));
	gPdf.fitTo(dh_Data,NumCPU(30));
	
	dh_Data.plotOn(Dalpha_frame);
	gPdf.plotOn(Dalpha_frame,LineColor(2));


	
	TCanvas *c_roofit=new TCanvas("c_roofit");
	c_roofit->cd();
	Dalpha_frame->Draw();	
	} // end if doRoofit
*/
	// } // end if dofit






	Float_t Dpt;
	Float_t Dalpha;
	Float_t Dchi2cl;
	Float_t Ddls;
	Float_t Dmass;
	Float_t TotalWeight;

	t_Ds_MCPrompt->SetBranchAddress("Dpt",&Dpt);
	t_Ds_MCPrompt->SetBranchAddress("Dalpha",&Dalpha);
	t_Ds_MCPrompt->SetBranchAddress("Dchi2cl",&Dchi2cl);
	// t_Ds_MCPrompt->SetBranchAddress("Ddls",&Ddls);
	t_Ds_MCPrompt->SetBranchAddress("Ddls_BestScale",&Ddls);
	t_Ds_MCPrompt->SetBranchAddress("Dmass",&Dmass);
	t_Ds_MCPrompt->SetBranchAddress("TotalWeight",&TotalWeight);


	nbin=200;
	TH1D *h_MCP_Ddls_Ori=new TH1D("h_MCP_Ddls_Ori","h_MCP_Ddls_Ori",nbin,binLow,binHigh); h_MCP_Ddls_Ori->Sumw2();
	TH1D *h_MCP_Ddls_Wet=new TH1D("h_MCP_Ddls_Wet","h_MCP_Ddls_Wet",nbin,binLow,binHigh); h_MCP_Ddls_Wet->Sumw2();
	TH1D *h_MCP_Ddls_FitWet=new TH1D("h_MCP_Ddls_FitWet","h_MCP_Ddls_FitWet",nbin,binLow,binHigh); h_MCP_Ddls_FitWet->Sumw2();

	long int nEntries=t_Ds_MCPrompt->GetEntries();
	for(int i =0; i<nEntries; i++){
		t_Ds_MCPrompt->GetEntry(i);
		// if(Dpt> Dpt_Low && Dpt < Dpt_High && Dalpha< Dalpha_cut && Dchi2cl > Dchi2cl_cut && Ddls> Ddls_cut && Dpt > var_cutLow && Dpt < var_cutHigh)
		// remove Ddls cut
		if(Dpt> Dpt_Low && Dpt < Dpt_High && Dalpha< Dalpha_cut && Dchi2cl > Dchi2cl_cut && Dpt > var_cutLow && Dpt < var_cutHigh){
			h_MCP_Ddls_Ori->Fill(Ddls,TotalWeight);
			double alpha_weight= gr_weight->Eval(Ddls);
			if(alpha_weight<0) {alpha_weight=0;}
			h_MCP_Ddls_Wet->Fill(Ddls,TotalWeight*alpha_weight);
			// h_MCP_Ddls_FitWet->Fill(Ddls,TotalWeight*gr_FitWeight->Eval(Ddls));
			
		}// end if cut

	} // end for  i<nEntries


		cout<<"Ori integral = "<<h_MCP_Ddls_Ori->Integral()<<endl;
		cout<<"Wet integral = "<<h_MCP_Ddls_Wet->Integral()<<endl;
		// cout<<"FitWet integral = "<<h_MCP_Ddls_FitWet->Integral()<<endl;
	TH1D *h_MCP_WtScale=new TH1D("h_MCP_WtScale","h_MCP_WtScale",1,0,1);
	h_MCP_WtScale->SetBinContent(1,h_MCP_Ddls_Ori->Integral()/h_MCP_Ddls_Wet->Integral());
	h_MCP_WtScale->Write();

	// TH1D *h_MCP_FitWtScale=new TH1D("h_MCP_FitWtScale","h_MCP_FitWtScale",1,0,1);
	// h_MCP_FitWtScale->SetBinContent(1,h_MCP_Ddls_Ori->Integral()/h_MCP_Ddls_FitWet->Integral());
	// h_MCP_FitWtScale->Write();

	TCanvas *c_MCPWet=new TCanvas("c_MCPWet");
	c_MCPWet->cd();


		h_MCP_Ddls_Ori->Scale(1/h_MCP_Ddls_Ori->Integral());
		h_MCP_Ddls_Ori->Rebin(rebinN);
		h_MCP_Ddls_Ori->Draw();
		h_MCP_Ddls_Wet->Draw("same");
		// h_MCP_Ddls_FitWet->SetLineColor(kMagenta+3);
		// h_MCP_Ddls_FitWet->Draw("same");
		
		c_MCPWet->BuildLegend();



	t_Ds_MCNonPrompt->SetBranchAddress("Dpt",&Dpt);
	t_Ds_MCNonPrompt->SetBranchAddress("Dalpha",&Dalpha);
	t_Ds_MCNonPrompt->SetBranchAddress("Dchi2cl",&Dchi2cl);
	// t_Ds_MCNonPrompt->SetBranchAddress("Ddls",&Ddls);
	t_Ds_MCNonPrompt->SetBranchAddress("Ddls_BestScale",&Ddls);
	t_Ds_MCNonPrompt->SetBranchAddress("Dmass",&Dmass);
	t_Ds_MCNonPrompt->SetBranchAddress("TotalWeight",&TotalWeight);

	TH1D *h_MCNP_Ddls_Ori=new TH1D("h_MCNP_Ddls_Ori","h_MCNP_Ddls_Ori",nbin,binLow,binHigh); h_MCNP_Ddls_Ori->Sumw2();
	TH1D *h_MCNP_Ddls_Wet=new TH1D("h_MCNP_Ddls_Wet","h_MCNP_Ddls_Wet",nbin,binLow,binHigh); h_MCNP_Ddls_Wet->Sumw2();
	TH1D *h_MCNP_Ddls_FitWet=new TH1D("h_MCNP_Ddls_FitWet","h_MCNP_Ddls_FitWet",nbin,binLow,binHigh); h_MCNP_Ddls_FitWet->Sumw2();


	nEntries=t_Ds_MCNonPrompt->GetEntries();
	for(int i =0; i<nEntries; i++){
		t_Ds_MCNonPrompt->GetEntry(i);
		// if(Dpt> Dpt_Low && Dpt < Dpt_High && Dalpha< Dalpha_cut && Dchi2cl > Dchi2cl_cut && Ddls> Ddls_cut && Dpt > var_cutLow && Dpt < var_cutHigh)
		// remove Ddls cut
		if(Dpt> Dpt_Low && Dpt < Dpt_High && Dalpha< Dalpha_cut && Dchi2cl > Dchi2cl_cut && Dpt > var_cutLow && Dpt < var_cutHigh){
			h_MCNP_Ddls_Ori->Fill(Ddls,TotalWeight);
			double alpha_weight= gr_weight->Eval(Ddls);
			if(alpha_weight<0) {alpha_weight=0;}
			h_MCNP_Ddls_Wet->Fill(Ddls,TotalWeight*alpha_weight);
			// h_MCNP_Ddls_FitWet->Fill(Ddls,TotalWeight*gr_FitWeight->Eval(Ddls));
			
		}// end if cut

	} // end for  i<nEntries


		cout<<"Ori integral = "<<h_MCNP_Ddls_Ori->Integral()<<endl;
		cout<<"Wet integral = "<<h_MCNP_Ddls_Wet->Integral()<<endl;
		// cout<<"FitWet integral = "<<h_MCNP_Ddls_FitWet->Integral()<<endl;
	TH1D *h_MCNP_WtScale=new TH1D("h_MCNP_WtScale","h_MCNP_WtScale",1,0,1);
	h_MCNP_WtScale->SetBinContent(1,h_MCNP_Ddls_Ori->Integral()/h_MCNP_Ddls_Wet->Integral());
	h_MCNP_WtScale->Write();

	// TH1D *h_MCNP_FitWtScale=new TH1D("h_MCNP_FitWtScale","h_MCNP_FitWtScale",1,0,1);
	// h_MCNP_FitWtScale->SetBinContent(1,h_MCNP_Ddls_Ori->Integral()/h_MCNP_Ddls_FitWet->Integral());
	// h_MCNP_FitWtScale->Write();

	TCanvas *c_MCNPWet=new TCanvas("c_MCNPWet");
	c_MCNPWet->cd();
		
		h_MCNP_Ddls_Ori->Scale(1/h_MCNP_Ddls_Ori->Integral());
		h_MCNP_Ddls_Ori->Rebin(rebinN);
		h_MCNP_Ddls_Ori->Draw();
		h_MCNP_Ddls_Wet->Draw("same");
		c_MCNPWet->BuildLegend();
		

		TH1D *h_MCMix_Ddls_Wet=new TH1D("h_MCMix_Ddls_Wet","h_MCMix_Ddls_Wet",nbin,binLow,binHigh); h_MCMix_Ddls_Wet->Sumw2();
		h_MCNP_Ddls_Wet->Scale(1/h_MCNP_Ddls_Wet->Integral());
		h_MCP_Ddls_Wet->Scale(1/h_MCP_Ddls_Wet->Integral());

		h_MCMix_Ddls_Wet->Add(h_MCP_Ddls_Wet,h_MCNP_Ddls_Wet, fr_withcut_prompt,1-fr_withcut_prompt);


		TCanvas *c_WetedMC=new TCanvas("c_WetedMC");
		c_WetedMC->cd(); 

		h_MCMix_Ddls_Wet->Rebin(rebinN)	;
		h_MCMix_Ddls_Wet->SetLineColor(kMagenta+3);
		h_MCMix_Ddls_Wet->Draw();
		h_max=h_MCMix_Ddls_Wet->GetMaximum();
		h_min=h_MCMix_Ddls_Wet->GetMinimum();
		h_MCMix_Ddls_Wet->SetMaximum(h_max+0.3*(h_max-h_min));
		h_MCMix_Ddls_Wet->SetMinimum(h_min-0.1*(h_max-h_min));

		h_MCP_Ddls_Wet->SetLineColor(2);
		h_MCP_Ddls_Wet->Rebin(rebinN);
		h_MCP_Ddls_Wet->Draw("same");
		h_MCNP_Ddls_Wet->SetLineColor(4);
		h_MCNP_Ddls_Wet->Rebin(rebinN);
		h_MCNP_Ddls_Wet->Draw("same");
		h_Ddls_sPlot->Draw("same");
		c_WetedMC->BuildLegend();


	c_WetedMC->SaveAs(Form("./plots_DataMC_Ddls/WetedMC_%s_%s_%s%.0fto%.0f.png",var_compare.Data(), Str_PbPb.Data(), var_cut.Data(), var_cutLow*100, var_cutHigh*100));



	// fitting Data & MC for weight


return 1;

		return 0;

}

int main(int argc, char*argv[]){


	if(argc==6){
		CompAndTune_DataMC_Ddls(atoi(argv[1]), argv[2], argv[3], atof(argv[4]), atof(argv[5]) );
	}else if(argc==10){	
		CompAndTune_DataMC_Ddls(atoi(argv[1]), argv[2], argv[3], atof(argv[4]), atof(argv[5]) , atof(argv[6]) , atof(argv[7]) , atof(argv[8]) , atof(argv[9]));
	}else{
		CompAndTune_DataMC_Ddls();
		cout<<"wrong number of input parameters , need 6\n int isPbPb=0, TString var_compare=DdxyzErr, TString var_cut=Dtrk1Pt , double var_cutLow=0.75, double var_cutHigh=1.25"<<endl;
		return 1;
	}
	//int Fit_sideband(int isPbPb=0, TString var_compare="DdxyzErr", TString var_cut="Dtrk1Pt" , double var_cutLow=0.75, double var_cutHigh=1.25){


	return 0;
}

