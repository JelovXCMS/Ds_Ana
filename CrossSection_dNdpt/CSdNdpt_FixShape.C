// error propergation https://en.wikipedia.org/wiki/Propagation_of_uncertainty


#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"


#include <TString.h>
#include <TCanvas.h>
#include "TBranch.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1F.h"
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
#include <RooFitResult.h>

#include "TFitter.h"
#include "TFitResult.h"

using namespace RooFit;
using namespace std;



TH1D *h_fr_PromptDs; // global for later use
// Bool_t verbose=true;

Int_t StartBin=0;

Int_t c_count=0;

TCanvas *c_PromptDs[500];

TCanvas *c_test[500];
Int_t count_test=0;


double shiftY=0;
double oneshift=0.075;

int isPbPb_Glo=-1;
// double PhiRatio_PbPb_=0.896;
//double PhiRatio_pp_=0.954;
int useDataPhiRatio_Glo=-1;

TH1D *h_CSfr_PromptDs_fun(TH1D *h_PromptDsCrossSectionPP_temp, TH1D* hBtoDsCrossSectionPP_temp, Int_t isPbPb=0,TString extraname=""){

  Int_t nbin_pt=nbin_pt_pp;
	double *bins_pt=bins_pt_pp;
  TString str_PbPb="pp";

  if(isPbPb==3){
    nbin_pt=nbin_pt_PbPb3;
    bins_pt=bins_pt_PbPb3;
    str_PbPb="PbPb3";
  }


	TH1D *h_CSfr_PromptDs_temp=new TH1D(Form("h_CSfr_PromptDs_%s%s",extraname.Data(),str_PbPb.Data()),Form("h_CSfr_PromptDs_%s%s",extraname.Data(),str_PbPb.Data()),nbin_pt,bins_pt); h_CSfr_PromptDs_temp->Sumw2();
  for(int ibin_pt=StartBin; ibin_pt<nbin_pt; ibin_pt++){
		double CS_PromptDs=h_PromptDsCrossSectionPP_temp->GetBinContent(ibin_pt+1);
		double CSErr_PromptDs=h_PromptDsCrossSectionPP_temp->GetBinError(ibin_pt+1);

		double CS_NonPromptDs=hBtoDsCrossSectionPP_temp->GetBinContent(ibin_pt+1);
		double CSErr_NonPromptDs=hBtoDsCrossSectionPP_temp->GetBinError(ibin_pt+1);

		h_CSfr_PromptDs_temp->SetBinContent(ibin_pt+1, CS_PromptDs/(CS_PromptDs+CS_NonPromptDs) ); 
		h_CSfr_PromptDs_temp->SetBinError(ibin_pt+1, ErrorPro_AoverAplusB(CSErr_PromptDs,CSErr_NonPromptDs,CS_PromptDs,CS_NonPromptDs));

	}

	return h_CSfr_PromptDs_temp;

}

TH1D *h_DsoverD0_fun(TH1D *h_PromptDsCrossSectionPP_temp, TH1D *hPromptDCrossSectionPP_temp, Int_t isPbPb=0,TString extraname=""){

  Int_t nbin_pt=nbin_pt_pp;
	double *bins_pt=bins_pt_pp;
  TString str_PbPb="pp";

  if(isPbPb==3){
    nbin_pt=nbin_pt_PbPb3;
    bins_pt=bins_pt_PbPb3;
    str_PbPb="PbPb3";
  }


	TH1D *h_DsoverD0PP_temp=new TH1D(Form("h_DsoverD0_%s%s",extraname.Data(),str_PbPb.Data()),Form("h_DsoverD0_%s%s",extraname.Data(),str_PbPb.Data()),nbin_pt,bins_pt); h_DsoverD0PP_temp->Sumw2();

  for(int ibin_pt=StartBin; ibin_pt<nbin_pt; ibin_pt++){
		double CS_PromptDs=h_PromptDsCrossSectionPP_temp->GetBinContent(ibin_pt+1);
		double CSErr_PromptDs=h_PromptDsCrossSectionPP_temp->GetBinError(ibin_pt+1);

		double CS_PromptD=hPromptDCrossSectionPP_temp->GetBinContent(ibin_pt+1);
		double CSErr_PromptD=hPromptDCrossSectionPP_temp->GetBinError(ibin_pt+1);

		double CSErr_DsoverD0=ErrorPro_AoverB(CSErr_PromptDs,CSErr_PromptD,CS_PromptDs,CS_PromptD);

		h_DsoverD0PP_temp->SetBinContent(ibin_pt+1, CS_PromptDs/CS_PromptD);
		h_DsoverD0PP_temp->SetBinError(ibin_pt+1, CSErr_DsoverD0);
	}
	return h_DsoverD0PP_temp;
}


TH1D *h_PromptDsCal_fun(TH1D *h_RawFitYield_temp1, TH1D* hBtoDsCrossSectionPP_temp1, Int_t isPbPb=0 ,TString str_h_eff="h_RecoNormEff", TString extraname="",TString DirName="",Bool_t fixPNPRatio=0, Double_t f0EffRatio=1, TH1D *h_fr_PromptDs_temp=h_fr_PromptDs, Int_t MC_weight=0, Bool_t doVarScan=false, TString Var="Dalpha", Int_t BR_option=0, Double_t LumiSumVar=0, Double_t BRphiVar=0, Double_t Brf0Var=0 )
{

	// statistic error quote fit only now

	vector<TString> savedirs;
	savedirs.push_back("CSdNdpt_Result");
	if(DirName!=""){
	savedirs.push_back(DirName);
	}

	double LumiNevt=LumiSum;
	Int_t nbin_pt=nbin_pt_pp;
	double *bins_pt=bins_pt_pp;
	TString str_eff_Prompt_phi=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_pp_Prompt_phikkpi.root",s_CutSet.Data());
	TString str_eff_NonPrompt_phi=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_pp_NonPrompt_phikkpi.root",s_CutSet.Data());
	TString str_eff_Prompt_f0=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_pp_Prompt_f0kkpi.root",s_CutSet.Data());
	TString str_eff_NonPrompt_f0=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_pp_NonPrompt_f0kkpi.root",s_CutSet.Data());
	TString str_PbPb="pp";

	if(isPbPb==3){
		LumiNevt=NevtPbPb3;	
		nbin_pt=nbin_pt_PbPb3;
		bins_pt=bins_pt_PbPb3;
	  str_eff_Prompt_phi=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_PbPb3_Prompt_phikkpi.root",s_CutSet.Data());
	  str_eff_NonPrompt_phi=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_PbPb3_NonPrompt_phikkpi.root",s_CutSet.Data());
	  str_eff_Prompt_f0=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_PbPb3_Prompt_f0kkpi.root",s_CutSet.Data());
	  str_eff_NonPrompt_f0=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_PbPb3_NonPrompt_f0kkpi.root",s_CutSet.Data());
		str_PbPb="PbPb3";
	}

	TH1D *h_PromptDsCrossSectionPP_temp=new TH1D(Form("h_PromptDs_%s%s",extraname.Data(),str_PbPb.Data()),Form("h_PromptDs_%s_%s",extraname.Data(),str_PbPb.Data()),nbin_pt,bins_pt); // need to modify this
	h_PromptDsCrossSectionPP_temp->Sumw2();

	TH1D *h_RawFitYield_temp=(TH1D*)h_RawFitYield_temp1->Clone("h_RawFitYield_temp"); // do not modify original Th1
	divideBinWidth(h_RawFitYield_temp);  // the calculataion is based on 1/dpt

	TH1D *hBtoDsCrossSectionPP_temp=(TH1D*)hBtoDsCrossSectionPP_temp1->Clone("hBtoDsCrossSectionPP_temp"); // do not modify original Th1 // this should be in anabin and already divided by bin width 

	// check binning 
	if(h_RawFitYield_temp->GetNbinsX()!= hBtoDsCrossSectionPP_temp->GetNbinsX() ){
		cout<<"  ---- Warning !!! the binning of hBtoDsCrossSectionPP_temp != h_RawFitYield_temp !!!! "<<endl;
	}

	// Efficiency
  TFile *f_Prompt_phikkpi= TFile::Open(str_eff_Prompt_phi.Data());
  TFile *f_NonPrompt_phikkpi= TFile::Open(str_eff_NonPrompt_phi.Data());
  TFile *f_Prompt_f0kkpi= TFile::Open(str_eff_Prompt_f0.Data());
  TFile *f_NonPrompt_f0kkpi= TFile::Open(str_eff_NonPrompt_f0.Data());

	// TString str_h_eff="h_RecoNormEff"; // if MC option ... change eff input here
  TH1D *h_Prompt_phikkpi=(TH1D*)f_Prompt_phikkpi->Get(str_h_eff);
  TH1D *h_NonPrompt_phikkpi=(TH1D*)f_NonPrompt_phikkpi->Get(str_h_eff);
  TH1D *h_Prompt_f0kkpi=(TH1D*)f_Prompt_f0kkpi->Get(str_h_eff);
  TH1D *h_NonPrompt_f0kkpi=(TH1D*)f_NonPrompt_f0kkpi->Get(str_h_eff);


	for(int ibin_pt=StartBin; ibin_pt<nbin_pt; ibin_pt++){

		double DptLow=bins_pt[ibin_pt];
		double DptHigh=bins_pt[ibin_pt+1];

		double N_yield=h_RawFitYield_temp->GetBinContent(ibin_pt+1);	
		double NErr_yield=h_RawFitYield_temp->GetBinError(ibin_pt+1);
		double NErrRel_yield=NErr_yield/N_yield;

    double Eff_Prompt_phi=h_Prompt_phikkpi->GetBinContent(ibin_pt+1);
    double Eff_NonPrompt_phi=h_NonPrompt_phikkpi->GetBinContent(ibin_pt+1);
    double Eff_Prompt_f0=h_Prompt_f0kkpi->GetBinContent(ibin_pt+1)*f0EffRatio;
    double Eff_NonPrompt_f0=h_NonPrompt_f0kkpi->GetBinContent(ibin_pt+1)*f0EffRatio;

    double EffErr_Prompt_phi=h_Prompt_phikkpi->GetBinError(ibin_pt+1);
    double EffErr_NonPrompt_phi=h_NonPrompt_phikkpi->GetBinError(ibin_pt+1);
    double EffErr_Prompt_f0=h_Prompt_f0kkpi->GetBinError(ibin_pt+1)*f0EffRatio;
    double EffErr_NonPrompt_f0=h_NonPrompt_f0kkpi->GetBinError(ibin_pt+1)*f0EffRatio;

		double CS_BtoDs=hBtoDsCrossSectionPP_temp->GetBinContent(ibin_pt+1);
		double CSErr_BtoDs=hBtoDsCrossSectionPP_temp->GetBinError(ibin_pt+1);

		double CS_PromptDs= (N_yield/(2*LumiNevt) - (CS_BtoDs* (BRphi*Eff_NonPrompt_phi + BRf0*Eff_NonPrompt_f0) ) ) / ( BRphi*Eff_Prompt_phi + BRf0*Eff_Prompt_f0 );

		double PhiRatio=1;
		if(isPbPb_Glo==0){
			PhiRatio=PhiRatio_pp;
		}else if(isPbPb_Glo>=3){
			PhiRatio=PhiRatio_PbPb;
		}else {
			cout<<"unknow isPbPb_Glo = "<<isPbPb_Glo<<endl;
			// return 99;
		}

		if(useDataPhiRatio_Glo){
			cout<<"useDataPhiRatio ,Ratio = "<<PhiRatio<<endl;
			CS_PromptDs = (N_yield*PhiRatio/(2*LumiNevt) - (CS_BtoDs* (BRphi*Eff_NonPrompt_phi ) ) ) / ( BRphi*Eff_Prompt_phi );
		}


		cout<<"\nfr pre = "<<CS_PromptDs/(CS_PromptDs+CS_BtoDs)<<" CS = "<<CS_PromptDs<<endl;

		// fr fix method for systematics
		double fr_Prompt=0;
		if(fixPNPRatio==1){
			cout<<"fixing PNP ratio "<<endl;
			fr_Prompt=h_fr_PromptDs_temp->GetBinContent(ibin_pt+1);
			CS_PromptDs=(N_yield/(2*LumiNevt)) / ( ( BRphi*Eff_Prompt_phi + BRf0*Eff_Prompt_f0 ) + (1.0-fr_Prompt)/(fr_Prompt)*(BRphi*Eff_NonPrompt_phi + BRf0*Eff_NonPrompt_f0) );
			cout<<"fr_Prompt = "<<fr_Prompt<<" CS = "<<CS_PromptDs<<endl;

			if(useDataPhiRatio_Glo){
  		  CS_PromptDs=(N_yield*PhiRatio/(2*LumiNevt)) / ( ( BRphi*Eff_Prompt_phi ) + (1.0-fr_Prompt)/(fr_Prompt)*(BRphi*Eff_NonPrompt_phi ) );		

			}

		}

		// double CSErr_PromptDs_fitonly=NErr_yield/2/LumiNevt/( BRphi*Eff_Prompt_phi + BRf0*Eff_Prompt_f0 );
		double CSErr_PromptDs_fitonly=CS_PromptDs*NErrRel_yield;


		double Err_BRPromptSum_test=sqrt(pow(BRphi*EffErr_Prompt_phi,2) + pow(BRf0*EffErr_Prompt_f0,2)); // assume uncorelated
		double Err_BRNonPromptSum_test=sqrt(pow(BRphi*EffErr_NonPrompt_phi,2) + pow(BRf0*EffErr_NonPrompt_f0,2)); // assume uncorelated
		double Err_CSBtoDs_and_eff_test=CS_BtoDs* (BRphi*Eff_NonPrompt_phi + BRf0*Eff_NonPrompt_f0) * sqrt( pow(CSErr_BtoDs/CS_BtoDs ,2) + pow(Err_BRNonPromptSum_test/(BRphi*Eff_NonPrompt_phi + BRf0*Eff_NonPrompt_f0) ,2)   );



		double Err_BRPromptSum=ErrorPro_aAPlusbB(EffErr_Prompt_phi,EffErr_Prompt_f0,BRphi,BRf0); // assume uncorelated
		double Err_BRNonPromptSum=ErrorPro_aAPlusbB(EffErr_NonPrompt_phi,EffErr_NonPrompt_f0,BRphi,BRf0); // assume uncorelated
		double Err_CSBtoDs_and_eff=ErrorPro_AtimesB(CSErr_BtoDs,Err_BRNonPromptSum,CS_BtoDs,(BRphi*Eff_NonPrompt_phi + BRf0*Eff_NonPrompt_f0));

		// double Err_NumeratorAll= ErrorPro_aAPlusbB_Corr(NErr_yield/2/LumiNevt, Err_CSBtoDs_and_eff); // treat totalDs and BtoDs CS (from fit result are totally correlated)
		double Err_NumeratorAll= ErrorPro_aAPlusbB(NErr_yield/2/LumiNevt, Err_CSBtoDs_and_eff); // treat totalDs and BtoDs CS (from fit result are totally correlated)
		double CSErr_PromptDs=ErrorPro_AoverB(Err_NumeratorAll,Err_BRPromptSum, N_yield/2/LumiNevt-(CS_BtoDs* (BRphi*Eff_NonPrompt_phi + BRf0*Eff_NonPrompt_f0) ), BRphi*Eff_Prompt_phi + BRf0*Eff_Prompt_f0 );


		if(CS_PromptDs<0){
			cout<<" warning !!CS_PromptDs <0 , set to 0"<<endl;
			CS_PromptDs=0;
		}

		h_PromptDsCrossSectionPP_temp->SetBinContent(ibin_pt+1, CS_PromptDs);
		h_PromptDsCrossSectionPP_temp->SetBinError(ibin_pt+1, CSErr_PromptDs_fitonly);
		// h_PromptDsCrossSectionPP_temp->SetBinError(ibin_pt+1, CSErr_PromptDs);

		// majority of  errors coming from 1. rawfit , 2 rawfit - NBtoDs , since fraction of NonpromptDs is ~50% is current cut, could cut on DCA to reduce Nonprompt ?

		if(verbose==true){
		cout<<"\n\n ibin_pt = "<<ibin_pt<<" pt : "<<DptLow<<" - "<<DptHigh<<" GeV "<<endl;
		cout<<"N_yield = "<< N_yield<<" +- "<<NErr_yield<<endl;
		cout<<"Eff_Prompt_phi = " <<Eff_Prompt_phi<<endl;
		cout<<"Eff_NonPrompt_phi = " <<Eff_NonPrompt_phi<<endl;
		cout<<"Eff_Prompt_f0 = " <<Eff_Prompt_f0<<endl;
		cout<<"Eff_NonPrompt_f0 = " <<Eff_NonPrompt_f0<<endl;
		cout<<"CS_BtoDs = "<<CS_BtoDs<<", CSErr_BtoDs = "<<CSErr_BtoDs<< endl;
		cout<<"N_NonPromptDs = "<<CS_BtoDs*2*LumiNevt<<endl;
		cout<<"N_NonPromptDs Raw= "<<CS_BtoDs*2*LumiNevt*(BRphi*Eff_NonPrompt_phi + BRf0*Eff_NonPrompt_f0)<<endl;

		cout<<" --- BRPromptEff = "<<( BRphi*Eff_Prompt_phi + BRf0*Eff_Prompt_f0 )<<" , Err_BRPromptSum_test = "<<Err_BRPromptSum_test <<" , Err_BRPromptSum = "<<Err_BRPromptSum<<endl;
		cout<<"BRNonPromptEff = "<< (BRphi*Eff_NonPrompt_phi + BRf0*Eff_NonPrompt_f0) <<" , Err_BRNonPromptSum_test = "<<Err_BRNonPromptSum_test <<" , Err_BRNonPromptSum = "<<Err_BRNonPromptSum<<endl;
		cout<<" , CSBtoDs_and_eff = "<< (CS_BtoDs* (BRphi*Eff_NonPrompt_phi + BRf0*Eff_NonPrompt_f0) )<<" , Err_CSBtoDs_and_eff_test = "<<Err_CSBtoDs_and_eff_test<<" , Err_CSBtoDs_and_eff = "<<Err_CSBtoDs_and_eff<<endl;
		cout<<" NumeratorAll (N_PromptDs) = "<<(N_yield/(2*LumiNevt) - (CS_BtoDs* (BRphi*Eff_NonPrompt_phi + BRf0*Eff_NonPrompt_f0) ) ) <<" ,Err_NumeratorAll = "<<Err_NumeratorAll<<endl;
		cout<<"CS_PromptDs = "<<CS_PromptDs<<" CSErr_PromptDs  = "<<CSErr_PromptDs<<" , CSErr_PromptDs_fitonly = "<<CSErr_PromptDs_fitonly<<endl;
		}

	}

	// TCanvas *c_PromptDsCrossSectionPP= new TCanvas("c_PromptDsCrossSectionPP","c_PromptDsCrossSectionPP",800,800);
	c_PromptDs[c_count]= new TCanvas(Form("c_PromptDs_%i",c_count),Form("c_PromptDs_%i",c_count),800,800);
	c_PromptDs[c_count]->cd();
	h_PromptDsCrossSectionPP_temp->Draw();
	c_PromptDs[c_count]->SetLogy();
	SavePlotDirs(c_PromptDs[c_count],Form("PromptDs_%s_%s",str_PbPb.Data(),extraname.Data()),savedirs);


	delete h_RawFitYield_temp;
	delete hBtoDsCrossSectionPP_temp;

	c_count++;

	return h_PromptDsCrossSectionPP_temp;

}


void CSdNdpt_FixShape(Int_t isPbPb=0, Int_t startbin=0, Bool_t doCutScan=0, int useDataPhiRatio=1){


	useDataPhiRatio_Glo=useDataPhiRatio;
	isPbPb_Glo=isPbPb;

	bool doDauPtScan=0;

// in parameters.h
//	double LumiSum=0.0382;  //(pb)
//	double BRphi=0.0227;
//	double BRf0=0.0115;
	// if(isPbPb==3){startbin=4;}
	if(isPbPb==3){startbin=2;}
	if(isPbPb==0){startbin=0;}
	StartBin=startbin;

	InitStyle();
	initParameter();

	TFile *f_NonPromptDs=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BtoDs_Results/output/BtoDs.root");
	TH1D *hPromptDCrossSectionPP_AnaBin=(TH1D*)f_NonPromptDs->Get("hPromptDCrossSectionPP_AnaBin");
	TH1D *hBtoDCrossSectionPP_AnaBin=(TH1D*)f_NonPromptDs->Get("hBtoDCrossSectionPP_AnaBin");
	TH1D *hBtoDsCrossSectionPP_AnaBin=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin");
	TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight");
	TH1D *hTotalErrBtoDsCrossSectionPP_AnaBin_pythiaWeight=(TH1D*)f_NonPromptDs->Get("hTotalErrBtoDsCrossSectionPP_AnaBin_pythiaWeight");
                                                                                                           
	TH1D *hPromptDdNdPtPbPb_AnaBin=(TH1D*)f_NonPromptDs->Get("hPromptDdNdPtPbPb_AnaBin");                      
	TH1D *hBtoDdNdPtPbPb_AnaBin=(TH1D*)f_NonPromptDs->Get("hBtoDdNdPtPbPb_AnaBin");
	TH1D *hBtoDsdNdPtPbPb_AnaBin=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin");
	TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight");
	TH1D *hTotalErrBtoDsdNdPtPbPb_AnaBin_pythiaWeight=(TH1D*)f_NonPromptDs->Get("hTotalErrBtoDsdNdPtPbPb_AnaBin_pythiaWeight");

	TH1D *hBtoDs_AnaBin_pythiaWeight=hBtoDsCrossSectionPP_AnaBin_pythiaWeight;
	TH1D *hPromptD_AnaBin=hPromptDCrossSectionPP_AnaBin;

	TString str_PbPb="pp";
	int nbin_pt=nbin_pt_pp;
	double *bins_pt=bins_pt_pp;

	TString outfile="output/PromptDsCrossSectionPP_FixShape.root";
	TString rawyieldfile=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/SignalFit/output%s/RawFitYield_FixShapeTrkPtScan_pp.root",s_CutSet.Data());
	// TString rawyieldfile=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/SignalFit/output%s/RawFitYield_FixShape_pp.root",s_CutSet.Data());

	double *DalphaMaxScan_bins=DalphaMaxScan_bins_pp;
	double *DdlsMinScan_bins=DdlsMinScan_bins_pp;
	double *Dchi2clMinScan_bins=Dchi2clMinScan_bins_pp;

	double *bins_DalphaMaxScan=bins_DalphaMaxScan_pp;
	double *bins_DdlsMinScan=bins_DdlsMinScan_pp;
	double *bins_Dchi2clMinScan=bins_Dchi2clMinScan_pp;


	if(isPbPb==3){
		str_PbPb="PbPb3";
		nbin_pt=nbin_pt_PbPb3;
		bins_pt=bins_pt_PbPb3;
		outfile="output/PromptDsdNdptPbPb3_FixShape.root";
		// rawyieldfile=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/SignalFit/output%s/RawFitYield_FixShape_PbPb3.root",s_CutSet.Data());
		rawyieldfile=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/SignalFit/output%s/RawFitYield_FixShapeTrkPtScan_PbPb3.root",s_CutSet.Data());
		hBtoDs_AnaBin_pythiaWeight=hBtoDsdNdPtPbPb_AnaBin_pythiaWeight;
		// hBtoDs_AnaBin_pythiaWeight=hBtoDsdNdPtPbPb_AnaBin;
		hPromptD_AnaBin=hPromptDdNdPtPbPb_AnaBin;

	  DalphaMaxScan_bins=DalphaMaxScan_bins_PbPb3;
	  DdlsMinScan_bins=DdlsMinScan_bins_PbPb3;
	  Dchi2clMinScan_bins=Dchi2clMinScan_bins_PbPb3;

	  bins_DalphaMaxScan=bins_DalphaMaxScan_PbPb3;
	  bins_DdlsMinScan=bins_DdlsMinScan_PbPb3;
	  bins_Dchi2clMinScan=bins_Dchi2clMinScan_PbPb3;


	}


	TFile *f_out=new TFile(outfile.Data(),"RECREATE");

	TFile *f_fitYield=TFile::Open(rawyieldfile.Data(),"READ");
	TH1D *h_RawRooFitYield=(TH1D*)f_fitYield->Get("h_RawRooFitYield");

	f_out->cd();


// TH1D *h_PromptDsCal_fun(TH1D *h_RawFitYield_temp1, TH1D* hBtoDsCrossSectionPP_temp1, Int_t isPbPb=0 ,TString str_h_eff="h_RecoNormEff", TString extraname="",TString DirName="",Bool_t fixPNPRatio=0, TH1D *h_fr_PromptDs_temp=h_fr_PromptDs, Int_t MC_weight=0, Bool_t doVarScan=false, TString Var="Dalpha", Int_t BR_option=0, Double_t LumiSumVar=0, Double_t BRphiVar=0, Double_t Brf0Var=0 )

	TH1D *h_PromptDs=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_AnaBin_pythiaWeight, isPbPb);
	TH1D *h_DsoverD0=h_DsoverD0_fun(h_PromptDs,hPromptD_AnaBin, isPbPb);
	TH1D *h_CSfr_PromptDs=h_CSfr_PromptDs_fun(h_PromptDs,hBtoDs_AnaBin_pythiaWeight, isPbPb);

	h_fr_PromptDs=h_CSfr_PromptDs; // set the global fr for later use

	// test fix ratio 
	// TH1D *h_PromptDs_test=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_AnaBin_pythiaWeight, isPbPb,"h_RecoNormEff","","",1);
	// h_PromptDs_test->SetName("h_PromptDs_test");

	TCanvas *c_test_hPromptDdNdPtPbPb_AnaBin=new TCanvas("c_test_hPromptDdNdPtPbPb_AnaBin","c_test_hPromptDdNdPtPbPb_AnaBin",800,800);
	c_test_hPromptDdNdPtPbPb_AnaBin->cd();
	hPromptD_AnaBin->Draw();


	f_out->cd();

	h_PromptDs->Write("",TObject::kOverwrite);
  h_DsoverD0->Write("",TObject::kOverwrite);
	h_CSfr_PromptDs->Write("",TObject::kOverwrite);

	// h_PromptDs_test->Write("",TObject::kOverwrite);

	gStyle->SetOptStat(0);

	TCanvas *c_DsoverD0= new TCanvas("c_DsoverD0","c_DsoverD0",800,800);
	c_DsoverD0->cd();
	h_DsoverD0->Draw();
	SavePlotDirs(c_DsoverD0,Form("DsoverD0_%s",str_PbPb.Data()), {"CSdNdpt_Result",str_PbPb});
	// c_DsoverD0PP->SaveAs(c_DsoverD0PP,"DsoverD0PP.pdf",{"Result"});	

	TCanvas *c_CSfr_PromptDs= new TCanvas("c_CSfr_PromptDs","c_CSfr_PromptDs",800,800);
	c_CSfr_PromptDs->cd();
	h_CSfr_PromptDs->Draw();
	SavePlotDirs(c_CSfr_PromptDs,Form("CSfr_PromptDs_%s",str_PbPb.Data()), {"CSdNdpt_Result",str_PbPb});


////-- f0 shape efficiency syst --////


// TH1D *h_PromptDsCal_fun(TH1D *h_RawFitYield_temp1, TH1D* hBtoDsCrossSectionPP_temp1, Int_t isPbPb=0 ,TString str_h_eff="h_RecoNormEff", TString extraname="",TString DirName="",Bool_t fixPNPRatio=0, Double_t f0EffRatio=1, TH1D *h_fr_PromptDs_temp=h_fr_PromptDs, Int_t MC_weight=0, Bool_t doVarScan=false, TString Var="Dalpha", Int_t BR_option=0, Double_t LumiSumVar=0, Double_t BRphiVar=0, Double_t Brf0Var=0 )
	f_out->cd();
	double f0Effup=2.22467;
	double f0Effdown=0.317062;

	TH1D *h_PromptDs_f0Effup=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_AnaBin_pythiaWeight, isPbPb, "h_RecoNormEff","f0Effup","f0Eff",1,f0Effup);
	TH1D *h_PromptDs_f0Effdown=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_AnaBin_pythiaWeight, isPbPb, "h_RecoNormEff","f0Effdown","f0Eff",1,f0Effdown);

	f_out->cd();
	h_PromptDs_f0Effup->Write("",TObject::kOverwrite);
	h_PromptDs_f0Effdown->Write("",TObject::kOverwrite);


////-- f0 shape efficiency syst --////


////-- MC pt shape --////

	f_out->cd();

	TH1D *h_PromptDs_MCSPythia=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_AnaBin_pythiaWeight, isPbPb, "h_RecoNormEff_Pythia","MCS_Pythia","MCS"); // MC shape by default Pythia
	TH1D *h_PromptDs_MCSFONLL=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_AnaBin_pythiaWeight, isPbPb, "h_RecoNormEff_FONLL","MCS_FONLL","MCS"); // MC shape by default Pythia

	f_out->cd();
	h_PromptDs_MCSPythia->Write("",TObject::kOverwrite);
	h_PromptDs_MCSFONLL->Write("",TObject::kOverwrite);

////-- end MC pt shape --////


//// -- BtoDs syst --////

	// FONLL 
	TH1D *hBtoDsCrossSectionPP_FONLL=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_FONLL");
	TH1D *hBtoDsdNdPtPbPb_FONLL=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_FONLL");
	TH1D *hBtoDs_FONLL=hBtoDsCrossSectionPP_FONLL;
	if(isPbPb){
	hBtoDs_FONLL=hBtoDsdNdPtPbPb_FONLL;
	}

  TH1D *h_PromptDs_BtoDsFONLL=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_FONLL, isPbPb, "h_RecoNormEff","BtoDs_FONLL","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_BtoDsFONLL->Write("",TObject::kOverwrite);

	// Raa up & down
	TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Raaup=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Raaup");
	TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Raadown=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Raadown");
	
	if(isPbPb){
	TH1D *h_PromptDs_Raaup=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Raaup,isPbPb, "h_RecoNormEff", "BtoDs_Raaup", "BtoDsSys");
	TH1D *h_PromptDs_Raadown=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Raadown,isPbPb, "h_RecoNormEff", "BtoDs_Raadown", "BtoDsSys");
	f_out->cd();
	h_PromptDs_Raaup->Write("",TObject::kOverwrite);
	h_PromptDs_Raadown->Write("",TObject::kOverwrite);
	}

	// FONLL shape with Data CS

	TH1D *hBtoDsCrossSectionPP_FONLLShape=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_FONLLShape");
	TH1D *hBtoDsdNdPtPbPb_FONLLShape=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_FONLLShape");
	TH1D *hBtoDs_FONLLShape=hBtoDsCrossSectionPP_FONLLShape;
	if(isPbPb){
	hBtoDs_FONLLShape=hBtoDsdNdPtPbPb_FONLLShape;
	}

  TH1D *h_PromptDs_BtoDsFONLLShape=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_FONLLShape, isPbPb, "h_RecoNormEff","BtoDs_FONLLShape","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_BtoDsFONLLShape->Write("",TObject::kOverwrite);

	// BtoD total err up & down

	TH1D *hBtoDsCrossSectionPP_Errup=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_Errup");
	TH1D *hBtoDsdNdPtPbPb_Errup=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Errup");
	TH1D *hBtoDs_Errup=hBtoDsCrossSectionPP_Errup;
	if(isPbPb){
	hBtoDs_Errup=hBtoDsdNdPtPbPb_Errup;
	}

  TH1D *h_PromptDs_Errup=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_Errup, isPbPb, "h_RecoNormEff","BtoDs_Errup","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_Errup->Write("",TObject::kOverwrite);


	TH1D *hBtoDsCrossSectionPP_Errdown=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_Errdown");
	TH1D *hBtoDsdNdPtPbPb_Errdown=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Errdown");
	TH1D *hBtoDs_Errdown=hBtoDsCrossSectionPP_Errdown;
	if(isPbPb){
	hBtoDs_Errdown=hBtoDsdNdPtPbPb_Errdown;
	}

  TH1D *h_PromptDs_Errdown=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_Errdown, isPbPb, "h_RecoNormEff","BtoDs_Errdown","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_Errdown->Write("",TObject::kOverwrite);


	// Scale Br up & down

	TH1D *hBtoDsCrossSectionPP_ScaleBRup=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBRup");
	TH1D *hBtoDsdNdPtPbPb_ScaleBRup=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBRup");
	TH1D *hBtoDs_ScaleBRup=hBtoDsCrossSectionPP_ScaleBRup;
	if(isPbPb){
	hBtoDs_ScaleBRup=hBtoDsdNdPtPbPb_ScaleBRup;
	}

  TH1D *h_PromptDs_ScaleBRup=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_ScaleBRup, isPbPb, "h_RecoNormEff","BtoDs_ScaleBRup","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_ScaleBRup->Write("",TObject::kOverwrite);

	TH1D *hBtoDsCrossSectionPP_ScaleBRdown=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBRdown");
	TH1D *hBtoDsdNdPtPbPb_ScaleBRdown=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBRdown");
	TH1D *hBtoDs_ScaleBRdown=hBtoDsCrossSectionPP_ScaleBRdown;
	if(isPbPb){
	hBtoDs_ScaleBRdown=hBtoDsdNdPtPbPb_ScaleBRdown;
	}

  TH1D *h_PromptDs_ScaleBRdown=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_ScaleBRdown, isPbPb, "h_RecoNormEff","BtoDs_ScaleBRdown","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_ScaleBRdown->Write("",TObject::kOverwrite);


	// Scale Br up & down seperately 


	TH1D *hBtoDsCrossSectionPP_ScaleBR_B0toDmax=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDmax");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_B0toDmax=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDmax");
	TH1D *hBtoDs_ScaleBR_B0toDmax=hBtoDsCrossSectionPP_ScaleBR_B0toDmax;
	if(isPbPb){
	hBtoDs_ScaleBR_B0toDmax=hBtoDsdNdPtPbPb_ScaleBR_B0toDmax;
	}

  TH1D *h_PromptDs_ScaleBR_B0toDmax=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_ScaleBR_B0toDmax, isPbPb, "h_RecoNormEff","BtoDs_ScaleBR_B0toDmax","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_ScaleBR_B0toDmax->Write("",TObject::kOverwrite);

	TH1D *hBtoDsCrossSectionPP_ScaleBR_B0toDmin=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDmin");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_B0toDmin=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDmin");
	TH1D *hBtoDs_ScaleBR_B0toDmin=hBtoDsCrossSectionPP_ScaleBR_B0toDmin;
	if(isPbPb){
	hBtoDs_ScaleBR_B0toDmin=hBtoDsdNdPtPbPb_ScaleBR_B0toDmin;
	}

  TH1D *h_PromptDs_ScaleBR_B0toDmin=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_ScaleBR_B0toDmin, isPbPb, "h_RecoNormEff","BtoDs_ScaleBR_B0toDmin","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_ScaleBR_B0toDmin->Write("",TObject::kOverwrite);

///

	TH1D *hBtoDsCrossSectionPP_ScaleBR_BptoDmax=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDmax");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_BptoDmax=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDmax");
	TH1D *hBtoDs_ScaleBR_BptoDmax=hBtoDsCrossSectionPP_ScaleBR_BptoDmax;
	if(isPbPb){
	hBtoDs_ScaleBR_BptoDmax=hBtoDsdNdPtPbPb_ScaleBR_BptoDmax;
	}

  TH1D *h_PromptDs_ScaleBR_BptoDmax=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_ScaleBR_BptoDmax, isPbPb, "h_RecoNormEff","BtoDs_ScaleBR_BptoDmax","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_ScaleBR_BptoDmax->Write("",TObject::kOverwrite);

	TH1D *hBtoDsCrossSectionPP_ScaleBR_BptoDmin=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDmin");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_BptoDmin=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDmin");
	TH1D *hBtoDs_ScaleBR_BptoDmin=hBtoDsCrossSectionPP_ScaleBR_BptoDmin;
	if(isPbPb){
	hBtoDs_ScaleBR_BptoDmin=hBtoDsdNdPtPbPb_ScaleBR_BptoDmin;
	}

  TH1D *h_PromptDs_ScaleBR_BptoDmin=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_ScaleBR_BptoDmin, isPbPb, "h_RecoNormEff","BtoDs_ScaleBR_BptoDmin","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_ScaleBR_BptoDmin->Write("",TObject::kOverwrite);

//


	TH1D *hBtoDsCrossSectionPP_ScaleBR_BstoDmax=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDmax");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_BstoDmax=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDmax");
	TH1D *hBtoDs_ScaleBR_BstoDmax=hBtoDsCrossSectionPP_ScaleBR_BstoDmax;
	if(isPbPb){
	hBtoDs_ScaleBR_BstoDmax=hBtoDsdNdPtPbPb_ScaleBR_BstoDmax;
	}

  TH1D *h_PromptDs_ScaleBR_BstoDmax=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_ScaleBR_BstoDmax, isPbPb, "h_RecoNormEff","BtoDs_ScaleBR_BstoDmax","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_ScaleBR_BstoDmax->Write("",TObject::kOverwrite);

	TH1D *hBtoDsCrossSectionPP_ScaleBR_BstoDmin=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDmin");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_BstoDmin=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDmin");
	TH1D *hBtoDs_ScaleBR_BstoDmin=hBtoDsCrossSectionPP_ScaleBR_BstoDmin;
	if(isPbPb){
	hBtoDs_ScaleBR_BstoDmin=hBtoDsdNdPtPbPb_ScaleBR_BstoDmin;
	}

  TH1D *h_PromptDs_ScaleBR_BstoDmin=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_ScaleBR_BstoDmin, isPbPb, "h_RecoNormEff","BtoDs_ScaleBR_BstoDmin","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_ScaleBR_BstoDmin->Write("",TObject::kOverwrite);

// Ds


	TH1D *hBtoDsCrossSectionPP_ScaleBR_B0toDsmax=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDsmax");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_B0toDsmax=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDsmax");
	TH1D *hBtoDs_ScaleBR_B0toDsmax=hBtoDsCrossSectionPP_ScaleBR_B0toDsmax;
	if(isPbPb){
	hBtoDs_ScaleBR_B0toDsmax=hBtoDsdNdPtPbPb_ScaleBR_B0toDsmax;
	}

  TH1D *h_PromptDs_ScaleBR_B0toDsmax=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_ScaleBR_B0toDsmax, isPbPb, "h_RecoNormEff","BtoDs_ScaleBR_B0toDsmax","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_ScaleBR_B0toDsmax->Write("",TObject::kOverwrite);

	TH1D *hBtoDsCrossSectionPP_ScaleBR_B0toDsmin=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDsmin");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_B0toDsmin=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDsmin");
	TH1D *hBtoDs_ScaleBR_B0toDsmin=hBtoDsCrossSectionPP_ScaleBR_B0toDsmin;
	if(isPbPb){
	hBtoDs_ScaleBR_B0toDsmin=hBtoDsdNdPtPbPb_ScaleBR_B0toDsmin;
	}

  TH1D *h_PromptDs_ScaleBR_B0toDsmin=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_ScaleBR_B0toDsmin, isPbPb, "h_RecoNormEff","BtoDs_ScaleBR_B0toDsmin","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_ScaleBR_B0toDsmin->Write("",TObject::kOverwrite);

///

	TH1D *hBtoDsCrossSectionPP_ScaleBR_BptoDsmax=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDsmax");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_BptoDsmax=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDsmax");
	TH1D *hBtoDs_ScaleBR_BptoDsmax=hBtoDsCrossSectionPP_ScaleBR_BptoDsmax;
	if(isPbPb){
	hBtoDs_ScaleBR_BptoDsmax=hBtoDsdNdPtPbPb_ScaleBR_BptoDsmax;
	}

  TH1D *h_PromptDs_ScaleBR_BptoDsmax=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_ScaleBR_BptoDsmax, isPbPb, "h_RecoNormEff","BtoDs_ScaleBR_BptoDsmax","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_ScaleBR_BptoDsmax->Write("",TObject::kOverwrite);

	TH1D *hBtoDsCrossSectionPP_ScaleBR_BptoDsmin=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDsmin");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_BptoDsmin=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDsmin");
	TH1D *hBtoDs_ScaleBR_BptoDsmin=hBtoDsCrossSectionPP_ScaleBR_BptoDsmin;
	if(isPbPb){
	hBtoDs_ScaleBR_BptoDsmin=hBtoDsdNdPtPbPb_ScaleBR_BptoDsmin;
	}

  TH1D *h_PromptDs_ScaleBR_BptoDsmin=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_ScaleBR_BptoDsmin, isPbPb, "h_RecoNormEff","BtoDs_ScaleBR_BptoDsmin","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_ScaleBR_BptoDsmin->Write("",TObject::kOverwrite);

//


	TH1D *hBtoDsCrossSectionPP_ScaleBR_BstoDsmax=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDsmax");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_BstoDsmax=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDsmax");
	TH1D *hBtoDs_ScaleBR_BstoDsmax=hBtoDsCrossSectionPP_ScaleBR_BstoDsmax;
	if(isPbPb){
	hBtoDs_ScaleBR_BstoDsmax=hBtoDsdNdPtPbPb_ScaleBR_BstoDsmax;
	}

  TH1D *h_PromptDs_ScaleBR_BstoDsmax=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_ScaleBR_BstoDsmax, isPbPb, "h_RecoNormEff","BtoDs_ScaleBR_BstoDsmax","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_ScaleBR_BstoDsmax->Write("",TObject::kOverwrite);

	TH1D *hBtoDsCrossSectionPP_ScaleBR_BstoDsmin=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDsmin");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_BstoDsmin=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDsmin");
	TH1D *hBtoDs_ScaleBR_BstoDsmin=hBtoDsCrossSectionPP_ScaleBR_BstoDsmin;
	if(isPbPb){
	hBtoDs_ScaleBR_BstoDsmin=hBtoDsdNdPtPbPb_ScaleBR_BstoDsmin;
	}

  TH1D *h_PromptDs_ScaleBR_BstoDsmin=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_ScaleBR_BstoDsmin, isPbPb, "h_RecoNormEff","BtoDs_ScaleBR_BstoDsmin","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_ScaleBR_BstoDsmin->Write("",TObject::kOverwrite);




	// Scale Fr_Z & Fr_p

	TH1D *hBtoDsCrossSectionPP_ScaleFrZ=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleFrZ");
	TH1D *hBtoDsdNdPtPbPb_ScaleFrZ=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleFrZ");
	TH1D *hBtoDs_ScaleFrZ=hBtoDsCrossSectionPP_ScaleFrZ;
	if(isPbPb){
	hBtoDs_ScaleFrZ=hBtoDsdNdPtPbPb_ScaleFrZ;
	}

  TH1D *h_PromptDs_ScaleFrZ=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_ScaleFrZ, isPbPb, "h_RecoNormEff","BtoDs_ScaleFrZ","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_ScaleFrZ->Write("",TObject::kOverwrite);



	TH1D *hBtoDsCrossSectionPP_ScaleFrp=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleFrp");
	TH1D *hBtoDsdNdPtPbPb_ScaleFrp=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleFrp");
	TH1D *hBtoDs_ScaleFrp=hBtoDsCrossSectionPP_ScaleFrp;
	if(isPbPb){
	hBtoDs_ScaleFrp=hBtoDsdNdPtPbPb_ScaleFrp;
	}

  TH1D *h_PromptDs_ScaleFrp=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_ScaleFrp, isPbPb, "h_RecoNormEff","BtoDs_ScaleFrp","BtoDsSys"); // MC shape by default Pythia
	f_out->cd();
	h_PromptDs_ScaleFrp->Write("",TObject::kOverwrite);


//// -- end BtoDs syst --////




////--  cut scan result --////

	// defaut , fix PNP ratio, method2, using BtoD yield

// TH1D *h_PromptDsCal_fun(TH1D *h_RawFitYield_temp1, TH1D* hBtoDsCrossSectionPP_temp1, Int_t isPbPb=0 ,TString str_h_eff="h_RecoNormEff", TString extraname="",TString DirName="",Bool_t fixPNPRatio=0, TH1D *h_fr_PromptDs_temp=h_fr_PromptDs, Int_t MC_weight=0, Bool_t doVarScan=false, TString Var="Dalpha", Int_t BR_option=0, Double_t LumiSumVar=0, Double_t BRphiVar=0, Double_t Brf0Var=0 )
	// TH1D *h_PromptDs=h_PromptDsCal_fun(h_RawRooFitYield,hBtoDs_AnaBin_pythiaWeight, isPbPb);
	// TH1D *h_RawRooFitYield=(TH1D*)f_fitYield->Get("h_RawRooFitYield");


	if(doCutScan){

	//-- DalphaMaxScan --//
	cout<<"\n\n ----- start DalphaMaxScan ----- \n"<<endl;
	f_out->cd();
	TH1D *h_PromptDs_DalphaMaxScan[nbin_DalphaMaxScan];
	TH1D *h_PromptDs_DalphaMaxScan2[nbin_DalphaMaxScan];
	TH1D *h_RawRooFitYield_DalphaMaxScan[nbin_DalphaMaxScan];
	for(int ibin=0; ibin<nbin_DalphaMaxScan; ibin++){
		cout<<"\n -- DalphaMaxScan bin "<<ibin<<" --"<<endl;
		h_RawRooFitYield_DalphaMaxScan[ibin]=(TH1D*)f_fitYield->Get(Form("h_RawRooFitYield_DalphaMaxScan_%i",ibin));
		h_PromptDs_DalphaMaxScan[ibin]=h_PromptDsCal_fun(h_RawRooFitYield_DalphaMaxScan[ibin],hBtoDs_AnaBin_pythiaWeight, isPbPb,Form("h_RecoNormEff_Dalpha_%i",ibin),Form("DalphaMaxScan_%i",ibin),"DalphaMaxScan",1);
		h_PromptDs_DalphaMaxScan[ibin]->SetName(Form("h_PromptDs_DalphaMaxScan_%i",ibin));	
		f_out->cd();
		h_PromptDs_DalphaMaxScan[ibin]->Write("",TObject::kOverwrite);

		h_PromptDs_DalphaMaxScan2[ibin]=h_PromptDsCal_fun(h_RawRooFitYield_DalphaMaxScan[ibin],hBtoDs_AnaBin_pythiaWeight, isPbPb,Form("h_RecoNormEff_Dalpha_%i",ibin),Form("DalphaMaxScan2_%i",ibin),"DalphaMaxScan2");
		h_PromptDs_DalphaMaxScan2[ibin]->SetName(Form("h_PromptDs_DalphaMaxScan2_%i",ibin));	
		f_out->cd();
		h_PromptDs_DalphaMaxScan2[ibin]->Write("",TObject::kOverwrite);

	} // end for ibin<nbin_DalphaMaxScan
		cout<<"end for ibin<nbin_DalphaMaxScan"<<endl;
	TH1D *h_PromptDs_DalphaMaxScan_pt[nbin_pt];
	TH1D *h_PromptDs_DalphaMaxScan2_pt[nbin_pt];
	for(int ibin_pt=StartBin; ibin_pt<nbin_pt; ibin_pt++){
		h_PromptDs_DalphaMaxScan_pt[ibin_pt]=new TH1D(Form("h_PromptDs_DalphaMaxScan_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), Form("h_PromptDs_DalphaMaxScan_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), nbin_DalphaMaxScan,bins_DalphaMaxScan );	
		h_PromptDs_DalphaMaxScan2_pt[ibin_pt]=new TH1D(Form("h_PromptDs_DalphaMaxScan2_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), Form("h_PromptDs_DalphaMaxScan2_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), nbin_DalphaMaxScan,bins_DalphaMaxScan );	
		for(int ibin=0; ibin<nbin_DalphaMaxScan; ibin++){
			cout<<"DalphaMaxScan_bins = "<<DalphaMaxScan_bins[ibin]<<" , bins_DalphaMaxScan = "<<bins_DalphaMaxScan[ibin]<<endl;
			h_PromptDs_DalphaMaxScan_pt[ibin_pt]->SetBinContent(ibin+1,h_PromptDs_DalphaMaxScan[ibin]->GetBinContent(ibin_pt+1));
			h_PromptDs_DalphaMaxScan_pt[ibin_pt]->SetBinError(ibin+1,h_PromptDs_DalphaMaxScan[ibin]->GetBinError(ibin_pt+1));

			h_PromptDs_DalphaMaxScan2_pt[ibin_pt]->SetBinContent(ibin+1,h_PromptDs_DalphaMaxScan2[ibin]->GetBinContent(ibin_pt+1));
			h_PromptDs_DalphaMaxScan2_pt[ibin_pt]->SetBinError(ibin+1,h_PromptDs_DalphaMaxScan2[ibin]->GetBinError(ibin_pt+1));
		}
		f_out->cd();
		h_PromptDs_DalphaMaxScan_pt[ibin_pt]->Write("",TObject::kOverwrite);
		h_PromptDs_DalphaMaxScan2_pt[ibin_pt]->Write("",TObject::kOverwrite);

		c_test[count_test]=new TCanvas(Form("c_test_%i",count_test),Form("c_test_%i",count_test),800,800);
		c_test[count_test]->cd();
		h_PromptDs_DalphaMaxScan_pt[ibin_pt]->Draw();

 	  shiftY=0;
  	TLatex *tl_DalphaMaxScan =new TLatex();
  	tl_DalphaMaxScan->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%.0f < pt < %.0f ",bins_pt[ibin_pt],bins_pt[ibin_pt+1])); shiftY-=oneshift;
  	tl_DalphaMaxScan->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s DalphaMaxScan ",str_PbPb.Data())); shiftY-=oneshift;

		SavePlotDirs(c_test[count_test],Form("PromptDs_%s_pt%.0fto%.0f_DalphaMaxScan",str_PbPb.Data(),bins_pt[ibin_pt],bins_pt[ibin_pt+1]),{"CSdNdpt_Result","DalphaMaxScan"});
		count_test++;

	} // end ibin_pt;


	//-- Dchi2clMinScan --//
	cout<<"\n\n ----- start Dchi2clMinScan ----- \n"<<endl;
	f_out->cd();
	TH1D *h_PromptDs_Dchi2clMinScan[nbin_Dchi2clMinScan];
	TH1D *h_PromptDs_Dchi2clMinScan2[nbin_Dchi2clMinScan];
	TH1D *h_RawRooFitYield_Dchi2clMinScan[nbin_Dchi2clMinScan];
	for(int ibin=0; ibin<nbin_Dchi2clMinScan; ibin++){
		cout<<"\n -- Dchi2clMinScan bin "<<ibin<<" --"<<endl;
		h_RawRooFitYield_Dchi2clMinScan[ibin]=(TH1D*)f_fitYield->Get(Form("h_RawRooFitYield_Dchi2clMinScan_%i",ibin));
		h_PromptDs_Dchi2clMinScan[ibin]=h_PromptDsCal_fun(h_RawRooFitYield_Dchi2clMinScan[ibin],hBtoDs_AnaBin_pythiaWeight, isPbPb,Form("h_RecoNormEff_Dchi2cl_%i",ibin),Form("Dchi2clMinScan_%i",ibin),"Dchi2clMinScan",1);
		h_PromptDs_Dchi2clMinScan[ibin]->SetName(Form("h_PromptDs_Dchi2clMinScan_%i",ibin));	
		f_out->cd();
		h_PromptDs_Dchi2clMinScan[ibin]->Write("",TObject::kOverwrite);

		h_PromptDs_Dchi2clMinScan2[ibin]=h_PromptDsCal_fun(h_RawRooFitYield_Dchi2clMinScan[ibin],hBtoDs_AnaBin_pythiaWeight, isPbPb,Form("h_RecoNormEff_Dchi2cl_%i",ibin),Form("Dchi2clMinScan2_%i",ibin),"Dchi2clMinScan2");
		h_PromptDs_Dchi2clMinScan2[ibin]->SetName(Form("h_PromptDs_Dchi2clMinScan2_%i",ibin));	
		f_out->cd();
		h_PromptDs_Dchi2clMinScan2[ibin]->Write("",TObject::kOverwrite);

	} // end for ibin<nbin_Dchi2clMinScan
		cout<<"end for ibin<nbin_Dchi2clMinScan"<<endl;
	TH1D *h_PromptDs_Dchi2clMinScan_pt[nbin_pt];
	TH1D *h_PromptDs_Dchi2clMinScan2_pt[nbin_pt];
	for(int ibin_pt=StartBin; ibin_pt<nbin_pt; ibin_pt++){
		h_PromptDs_Dchi2clMinScan_pt[ibin_pt]=new TH1D(Form("h_PromptDs_Dchi2clMinScan_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), Form("h_PromptDs_Dchi2clMinScan_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), nbin_Dchi2clMinScan,bins_Dchi2clMinScan );	
		h_PromptDs_Dchi2clMinScan2_pt[ibin_pt]=new TH1D(Form("h_PromptDs_Dchi2clMinScan2_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), Form("h_PromptDs_Dchi2clMinScan2_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), nbin_Dchi2clMinScan,bins_Dchi2clMinScan );	
		for(int ibin=0; ibin<nbin_Dchi2clMinScan; ibin++){
			cout<<"Dchi2clMinScan_bins = "<<Dchi2clMinScan_bins[ibin]<<" , bins_Dchi2clMinScan = "<<bins_Dchi2clMinScan[ibin]<<endl;
			h_PromptDs_Dchi2clMinScan_pt[ibin_pt]->SetBinContent(ibin+1,h_PromptDs_Dchi2clMinScan[ibin]->GetBinContent(ibin_pt+1));
			h_PromptDs_Dchi2clMinScan_pt[ibin_pt]->SetBinError(ibin+1,h_PromptDs_Dchi2clMinScan[ibin]->GetBinError(ibin_pt+1));

			h_PromptDs_Dchi2clMinScan2_pt[ibin_pt]->SetBinContent(ibin+1,h_PromptDs_Dchi2clMinScan2[ibin]->GetBinContent(ibin_pt+1));
			h_PromptDs_Dchi2clMinScan2_pt[ibin_pt]->SetBinError(ibin+1,h_PromptDs_Dchi2clMinScan2[ibin]->GetBinError(ibin_pt+1));
		}
		f_out->cd();
		h_PromptDs_Dchi2clMinScan_pt[ibin_pt]->Write("",TObject::kOverwrite);
		h_PromptDs_Dchi2clMinScan2_pt[ibin_pt]->Write("",TObject::kOverwrite);

		c_test[count_test]=new TCanvas(Form("c_test_%i",count_test),Form("c_test_%i",count_test),800,800);
		c_test[count_test]->cd();
		h_PromptDs_Dchi2clMinScan_pt[ibin_pt]->Draw();

 	  shiftY=0;
  	TLatex *tl_Dchi2clMinScan =new TLatex();
  	tl_Dchi2clMinScan->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%.0f < pt < %.0f ",bins_pt[ibin_pt],bins_pt[ibin_pt+1])); shiftY-=oneshift;
  	tl_Dchi2clMinScan->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s Dchi2clMinScan ",str_PbPb.Data())); shiftY-=oneshift;

		SavePlotDirs(c_test[count_test],Form("PromptDs_%s_pt%.0fto%.0f_Dchi2clMinScan",str_PbPb.Data(),bins_pt[ibin_pt],bins_pt[ibin_pt+1]),{"CSdNdpt_Result","Dchi2clMinScan"});
		count_test++;

	} // end ibin_pt;


	//-- DdlsMinScan --//
	cout<<"\n\n ----- start DdlsMinScan ----- \n"<<endl;
	f_out->cd();
	TH1D *h_PromptDs_DdlsMinScan[nbin_DdlsMinScan];
	TH1D *h_PromptDs_DdlsMinScan2[nbin_DdlsMinScan];
	TH1D *h_RawRooFitYield_DdlsMinScan[nbin_DdlsMinScan];
	for(int ibin=0; ibin<nbin_DdlsMinScan; ibin++){
		cout<<"\n -- DdlsMinScan bin "<<ibin<<" --"<<endl;
		h_RawRooFitYield_DdlsMinScan[ibin]=(TH1D*)f_fitYield->Get(Form("h_RawRooFitYield_DdlsMinScan_%i",ibin));
		h_PromptDs_DdlsMinScan[ibin]=h_PromptDsCal_fun(h_RawRooFitYield_DdlsMinScan[ibin],hBtoDs_AnaBin_pythiaWeight, isPbPb,Form("h_RecoNormEff_Ddls_%i",ibin),Form("DdlsMinScan_%i",ibin),"DdlsMinScan",1);
		h_PromptDs_DdlsMinScan[ibin]->SetName(Form("h_PromptDs_DdlsMinScan_%i",ibin));	
		f_out->cd();
		h_PromptDs_DdlsMinScan[ibin]->Write("",TObject::kOverwrite);

		h_PromptDs_DdlsMinScan2[ibin]=h_PromptDsCal_fun(h_RawRooFitYield_DdlsMinScan[ibin],hBtoDs_AnaBin_pythiaWeight, isPbPb,Form("h_RecoNormEff_Ddls_%i",ibin),Form("DdlsMinScan2_%i",ibin),"DdlsMinScan2");
		h_PromptDs_DdlsMinScan2[ibin]->SetName(Form("h_PromptDs_DdlsMinScan2_%i",ibin));	
		f_out->cd();
		h_PromptDs_DdlsMinScan2[ibin]->Write("",TObject::kOverwrite);

	} // end for ibin<nbin_DdlsMinScan
		cout<<"end for ibin<nbin_DdlsMinScan"<<endl;
	TH1D *h_PromptDs_DdlsMinScan_pt[nbin_pt];
	TH1D *h_PromptDs_DdlsMinScan2_pt[nbin_pt];
	for(int ibin_pt=StartBin; ibin_pt<nbin_pt; ibin_pt++){
		h_PromptDs_DdlsMinScan_pt[ibin_pt]=new TH1D(Form("h_PromptDs_DdlsMinScan_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), Form("h_PromptDs_DdlsMinScan_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), nbin_DdlsMinScan,bins_DdlsMinScan );	
		h_PromptDs_DdlsMinScan2_pt[ibin_pt]=new TH1D(Form("h_PromptDs_DdlsMinScan2_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), Form("h_PromptDs_DdlsMinScan2_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), nbin_DdlsMinScan,bins_DdlsMinScan );	
		for(int ibin=0; ibin<nbin_DdlsMinScan; ibin++){
			cout<<"DdlsMinScan_bins = "<<DdlsMinScan_bins[ibin]<<" , bins_DdlsMinScan = "<<bins_DdlsMinScan[ibin]<<endl;
			h_PromptDs_DdlsMinScan_pt[ibin_pt]->SetBinContent(ibin+1,h_PromptDs_DdlsMinScan[ibin]->GetBinContent(ibin_pt+1));
			h_PromptDs_DdlsMinScan_pt[ibin_pt]->SetBinError(ibin+1,h_PromptDs_DdlsMinScan[ibin]->GetBinError(ibin_pt+1));

			h_PromptDs_DdlsMinScan2_pt[ibin_pt]->SetBinContent(ibin+1,h_PromptDs_DdlsMinScan2[ibin]->GetBinContent(ibin_pt+1));
			h_PromptDs_DdlsMinScan2_pt[ibin_pt]->SetBinError(ibin+1,h_PromptDs_DdlsMinScan2[ibin]->GetBinError(ibin_pt+1));
		}
		f_out->cd();
		h_PromptDs_DdlsMinScan_pt[ibin_pt]->Write("",TObject::kOverwrite);
		h_PromptDs_DdlsMinScan2_pt[ibin_pt]->Write("",TObject::kOverwrite);

		c_test[count_test]=new TCanvas(Form("c_test_%i",count_test),Form("c_test_%i",count_test),800,800);
		c_test[count_test]->cd();
		h_PromptDs_DdlsMinScan_pt[ibin_pt]->Draw();

 	  shiftY=0;
  	TLatex *tl_DdlsMinScan =new TLatex();
  	tl_DdlsMinScan->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%.0f < pt < %.0f ",bins_pt[ibin_pt],bins_pt[ibin_pt+1])); shiftY-=oneshift;
  	tl_DdlsMinScan->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s DdlsMinScan ",str_PbPb.Data())); shiftY-=oneshift;

		SavePlotDirs(c_test[count_test],Form("PromptDs_%s_pt%.0fto%.0f_DdlsMinScan",str_PbPb.Data(),bins_pt[ibin_pt],bins_pt[ibin_pt+1]),{"CSdNdpt_Result","DdlsMinScan"});
		count_test++;

	} // end ibin_pt;

	// -- end DdlsMinScan


	//-- PhiMassScan --//
	cout<<"\n\n ----- start PhiMassScan ----- \n"<<endl;
	f_out->cd();
	TH1D *h_PromptDs_PhiMassScan[nbin_PhiMassScan];
	TH1D *h_PromptDs_PhiMassScan2[nbin_PhiMassScan];
	TH1D *h_RawRooFitYield_PhiMassScan[nbin_PhiMassScan];
	for(int ibin=0; ibin<nbin_PhiMassScan; ibin++){
		cout<<"\n -- PhiMassScan bin "<<ibin<<" --"<<endl;
		h_RawRooFitYield_PhiMassScan[ibin]=(TH1D*)f_fitYield->Get(Form("h_RawRooFitYield_PhiMassScan_%i",ibin));
		h_PromptDs_PhiMassScan[ibin]=h_PromptDsCal_fun(h_RawRooFitYield_PhiMassScan[ibin],hBtoDs_AnaBin_pythiaWeight, isPbPb,Form("h_RecoNormEff_PhiMass_%i",ibin),Form("PhiMassScan_%i",ibin),"PhiMassScan",1);
		h_PromptDs_PhiMassScan[ibin]->SetName(Form("h_PromptDs_PhiMassScan_%i",ibin));	
		f_out->cd();
		h_PromptDs_PhiMassScan[ibin]->Write("",TObject::kOverwrite);

		h_PromptDs_PhiMassScan2[ibin]=h_PromptDsCal_fun(h_RawRooFitYield_PhiMassScan[ibin],hBtoDs_AnaBin_pythiaWeight, isPbPb,Form("h_RecoNormEff_PhiMass_%i",ibin),Form("PhiMassScan2_%i",ibin),"PhiMassScan2");
		h_PromptDs_PhiMassScan2[ibin]->SetName(Form("h_PromptDs_PhiMassScan2_%i",ibin));	
		f_out->cd();
		h_PromptDs_PhiMassScan2[ibin]->Write("",TObject::kOverwrite);

	} // end for ibin<nbin_PhiMassScan
		cout<<"end for ibin<nbin_PhiMassScan"<<endl;
	TH1D *h_PromptDs_PhiMassScan_pt[nbin_pt];
	TH1D *h_PromptDs_PhiMassScan2_pt[nbin_pt];
	for(int ibin_pt=StartBin; ibin_pt<nbin_pt; ibin_pt++){
		h_PromptDs_PhiMassScan_pt[ibin_pt]=new TH1D(Form("h_PromptDs_PhiMassScan_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), Form("h_PromptDs_PhiMassScan_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), nbin_PhiMassScan,bins_PhiMassScan );	
		h_PromptDs_PhiMassScan2_pt[ibin_pt]=new TH1D(Form("h_PromptDs_PhiMassScan2_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), Form("h_PromptDs_PhiMassScan2_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), nbin_PhiMassScan,bins_PhiMassScan );	
		for(int ibin=0; ibin<nbin_PhiMassScan; ibin++){
			cout<<"PhiMassScan_bins = "<<PhiMassScan_bins[ibin]<<" , bins_PhiMassScan = "<<bins_PhiMassScan[ibin]<<endl;
			h_PromptDs_PhiMassScan_pt[ibin_pt]->SetBinContent(ibin+1,h_PromptDs_PhiMassScan[ibin]->GetBinContent(ibin_pt+1));
			h_PromptDs_PhiMassScan_pt[ibin_pt]->SetBinError(ibin+1,h_PromptDs_PhiMassScan[ibin]->GetBinError(ibin_pt+1));

			h_PromptDs_PhiMassScan2_pt[ibin_pt]->SetBinContent(ibin+1,h_PromptDs_PhiMassScan2[ibin]->GetBinContent(ibin_pt+1));
			h_PromptDs_PhiMassScan2_pt[ibin_pt]->SetBinError(ibin+1,h_PromptDs_PhiMassScan2[ibin]->GetBinError(ibin_pt+1));
		}
		f_out->cd();
		h_PromptDs_PhiMassScan_pt[ibin_pt]->Write("",TObject::kOverwrite);
		h_PromptDs_PhiMassScan2_pt[ibin_pt]->Write("",TObject::kOverwrite);

		c_test[count_test]=new TCanvas(Form("c_test_%i",count_test),Form("c_test_%i",count_test),800,800);
		c_test[count_test]->cd();
		h_PromptDs_PhiMassScan_pt[ibin_pt]->Draw();

 	  shiftY=0;
  	TLatex *tl_PhiMassScan =new TLatex();
  	tl_PhiMassScan->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%.0f < pt < %.0f ",bins_pt[ibin_pt],bins_pt[ibin_pt+1])); shiftY-=oneshift;
  	tl_PhiMassScan->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s PhiMassScan ",str_PbPb.Data())); shiftY-=oneshift;

		SavePlotDirs(c_test[count_test],Form("PromptDs_%s_pt%.0fto%.0f_PhiMassScan",str_PbPb.Data(),bins_pt[ibin_pt],bins_pt[ibin_pt+1]),{"CSdNdpt_Result","PhiMassScan"});
		count_test++;

	} // end ibin_pt;

	// -- end PhiMassScan

	if(doDauPtScan){	
	
	//-- KaonPtScan --//
	cout<<"\n\n ----- start KaonPtScan ----- \n"<<endl;
	f_out->cd();
	TH1D *h_PromptDs_KaonPtScan[nbin_DauPtScan];
	TH1D *h_PromptDs_KaonPtScan2[nbin_DauPtScan];
	TH1D *h_RawRooFitYield_KaonPtScan[nbin_DauPtScan];
	for(int ibin=0; ibin<nbin_DauPtScan; ibin++){
		cout<<"\n -- KaonPtScan bin "<<ibin<<" --"<<endl;
		h_RawRooFitYield_KaonPtScan[ibin]=(TH1D*)f_fitYield->Get(Form("h_RawRooFitYield_KaonPtScan_%i",ibin));
		h_PromptDs_KaonPtScan[ibin]=h_PromptDsCal_fun(h_RawRooFitYield_KaonPtScan[ibin],hBtoDs_AnaBin_pythiaWeight, isPbPb,Form("h_RecoNormEff_KaonPt_%i",ibin),Form("KaonPtScan_%i",ibin),"KaonPtScan",1);
		h_PromptDs_KaonPtScan[ibin]->SetName(Form("h_PromptDs_KaonPtScan_%i",ibin));	
		f_out->cd();
		h_PromptDs_KaonPtScan[ibin]->Write("",TObject::kOverwrite);

		h_PromptDs_KaonPtScan2[ibin]=h_PromptDsCal_fun(h_RawRooFitYield_KaonPtScan[ibin],hBtoDs_AnaBin_pythiaWeight, isPbPb,Form("h_RecoNormEff_KaonPt_%i",ibin),Form("KaonPtScan2_%i",ibin),"KaonPtScan2");
		h_PromptDs_KaonPtScan2[ibin]->SetName(Form("h_PromptDs_KaonPtScan2_%i",ibin));	
		f_out->cd();
		h_PromptDs_KaonPtScan2[ibin]->Write("",TObject::kOverwrite);

	} // end for ibin<nbin_DauPtScan
		cout<<"end for ibin<nbin_DauPtScan"<<endl;
	TH1D *h_PromptDs_KaonPtScan_pt[nbin_pt];
	TH1D *h_PromptDs_KaonPtScan2_pt[nbin_pt];
	for(int ibin_pt=StartBin; ibin_pt<nbin_pt; ibin_pt++){
		h_PromptDs_KaonPtScan_pt[ibin_pt]=new TH1D(Form("h_PromptDs_KaonPtScan_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), Form("h_PromptDs_KaonPtScan_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), nbin_DauPtScan,bins_DauPtScan );	
		h_PromptDs_KaonPtScan2_pt[ibin_pt]=new TH1D(Form("h_PromptDs_KaonPtScan2_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), Form("h_PromptDs_KaonPtScan2_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), nbin_DauPtScan,bins_DauPtScan );	
		for(int ibin=0; ibin<nbin_DauPtScan; ibin++){
//			// cout<<"KaonPtScan_bins = "<<KaonPtScan_bins[ibin]<<" , bins_KaonPtScan = "<<bins_KaonPtScan[ibin]<<endl;
			h_PromptDs_KaonPtScan_pt[ibin_pt]->SetBinContent(ibin+1,h_PromptDs_KaonPtScan[ibin]->GetBinContent(ibin_pt+1));
			h_PromptDs_KaonPtScan_pt[ibin_pt]->SetBinError(ibin+1,h_PromptDs_KaonPtScan[ibin]->GetBinError(ibin_pt+1));

			h_PromptDs_KaonPtScan2_pt[ibin_pt]->SetBinContent(ibin+1,h_PromptDs_KaonPtScan2[ibin]->GetBinContent(ibin_pt+1));
			h_PromptDs_KaonPtScan2_pt[ibin_pt]->SetBinError(ibin+1,h_PromptDs_KaonPtScan2[ibin]->GetBinError(ibin_pt+1));
		}
		f_out->cd();
		h_PromptDs_KaonPtScan_pt[ibin_pt]->Write("",TObject::kOverwrite);
		h_PromptDs_KaonPtScan2_pt[ibin_pt]->Write("",TObject::kOverwrite);

		c_test[count_test]=new TCanvas(Form("c_test_%i",count_test),Form("c_test_%i",count_test),800,800);
		c_test[count_test]->cd();
		h_PromptDs_KaonPtScan_pt[ibin_pt]->Draw();

 	  shiftY=0;
  	TLatex *tl_KaonPtScan =new TLatex();
  	tl_KaonPtScan->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%.0f < pt < %.0f ",bins_pt[ibin_pt],bins_pt[ibin_pt+1])); shiftY-=oneshift;
  	tl_KaonPtScan->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s KaonPtScan ",str_PbPb.Data())); shiftY-=oneshift;

		SavePlotDirs(c_test[count_test],Form("PromptDs_%s_pt%.0fto%.0f_KaonPtScan",str_PbPb.Data(),bins_pt[ibin_pt],bins_pt[ibin_pt+1]),{"CSdNdpt_Result","KaonPtScan"});
		count_test++;

	} // end ibin_pt;

	// -- end KaonPtScan

	//-- PionPtScan --//
	cout<<"\n\n ----- start PionPtScan ----- \n"<<endl;
	f_out->cd();
	TH1D *h_PromptDs_PionPtScan[nbin_DauPtScan];
	TH1D *h_PromptDs_PionPtScan2[nbin_DauPtScan];
	TH1D *h_RawRooFitYield_PionPtScan[nbin_DauPtScan];
	for(int ibin=0; ibin<nbin_DauPtScan; ibin++){
		cout<<"\n -- PionPtScan bin "<<ibin<<" --"<<endl;
		h_RawRooFitYield_PionPtScan[ibin]=(TH1D*)f_fitYield->Get(Form("h_RawRooFitYield_PionPtScan_%i",ibin));
		h_PromptDs_PionPtScan[ibin]=h_PromptDsCal_fun(h_RawRooFitYield_PionPtScan[ibin],hBtoDs_AnaBin_pythiaWeight, isPbPb,Form("h_RecoNormEff_PionPt_%i",ibin),Form("PionPtScan_%i",ibin),"PionPtScan",1);
		h_PromptDs_PionPtScan[ibin]->SetName(Form("h_PromptDs_PionPtScan_%i",ibin));	
		f_out->cd();
		h_PromptDs_PionPtScan[ibin]->Write("",TObject::kOverwrite);

		h_PromptDs_PionPtScan2[ibin]=h_PromptDsCal_fun(h_RawRooFitYield_PionPtScan[ibin],hBtoDs_AnaBin_pythiaWeight, isPbPb,Form("h_RecoNormEff_PionPt_%i",ibin),Form("PionPtScan2_%i",ibin),"PionPtScan2");
		h_PromptDs_PionPtScan2[ibin]->SetName(Form("h_PromptDs_PionPtScan2_%i",ibin));	
		f_out->cd();
		h_PromptDs_PionPtScan2[ibin]->Write("",TObject::kOverwrite);

	} // end for ibin<nbin_DauPtScan
		cout<<"end for ibin<nbin_DauPtScan"<<endl;
	TH1D *h_PromptDs_PionPtScan_pt[nbin_pt];
	TH1D *h_PromptDs_PionPtScan2_pt[nbin_pt];
	for(int ibin_pt=StartBin; ibin_pt<nbin_pt; ibin_pt++){
		h_PromptDs_PionPtScan_pt[ibin_pt]=new TH1D(Form("h_PromptDs_PionPtScan_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), Form("h_PromptDs_PionPtScan_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), nbin_DauPtScan,bins_DauPtScan );	
		h_PromptDs_PionPtScan2_pt[ibin_pt]=new TH1D(Form("h_PromptDs_PionPtScan2_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), Form("h_PromptDs_PionPtScan2_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), nbin_DauPtScan,bins_DauPtScan );	
		for(int ibin=0; ibin<nbin_DauPtScan; ibin++){
//			// cout<<"PionPtScan_bins = "<<PionPtScan_bins[ibin]<<" , bins_PionPtScan = "<<bins_PionPtScan[ibin]<<endl;
			h_PromptDs_PionPtScan_pt[ibin_pt]->SetBinContent(ibin+1,h_PromptDs_PionPtScan[ibin]->GetBinContent(ibin_pt+1));
			h_PromptDs_PionPtScan_pt[ibin_pt]->SetBinError(ibin+1,h_PromptDs_PionPtScan[ibin]->GetBinError(ibin_pt+1));

			h_PromptDs_PionPtScan2_pt[ibin_pt]->SetBinContent(ibin+1,h_PromptDs_PionPtScan2[ibin]->GetBinContent(ibin_pt+1));
			h_PromptDs_PionPtScan2_pt[ibin_pt]->SetBinError(ibin+1,h_PromptDs_PionPtScan2[ibin]->GetBinError(ibin_pt+1));
		}
		f_out->cd();
		h_PromptDs_PionPtScan_pt[ibin_pt]->Write("",TObject::kOverwrite);
		h_PromptDs_PionPtScan2_pt[ibin_pt]->Write("",TObject::kOverwrite);

		c_test[count_test]=new TCanvas(Form("c_test_%i",count_test),Form("c_test_%i",count_test),800,800);
		c_test[count_test]->cd();
		h_PromptDs_PionPtScan_pt[ibin_pt]->Draw();

 	  shiftY=0;
  	TLatex *tl_PionPtScan =new TLatex();
  	tl_PionPtScan->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%.0f < pt < %.0f ",bins_pt[ibin_pt],bins_pt[ibin_pt+1])); shiftY-=oneshift;
  	tl_PionPtScan->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s PionPtScan ",str_PbPb.Data())); shiftY-=oneshift;

		SavePlotDirs(c_test[count_test],Form("PromptDs_%s_pt%.0fto%.0f_PionPtScan",str_PbPb.Data(),bins_pt[ibin_pt],bins_pt[ibin_pt+1]),{"CSdNdpt_Result","PionPtScan"});
		count_test++;

	} // end ibin_pt;

	// -- end PionPtScan

	//-- AllDauPtScan --//
	cout<<"\n\n ----- start AllDauPtScan ----- \n"<<endl;
	f_out->cd();
	TH1D *h_PromptDs_AllDauPtScan[nbin_DauPtScan];
	TH1D *h_PromptDs_AllDauPtScan2[nbin_DauPtScan];
	TH1D *h_RawRooFitYield_AllDauPtScan[nbin_DauPtScan];
	for(int ibin=0; ibin<nbin_DauPtScan; ibin++){
		cout<<"\n -- AllDauPtScan bin "<<ibin<<" --"<<endl;
		h_RawRooFitYield_AllDauPtScan[ibin]=(TH1D*)f_fitYield->Get(Form("h_RawRooFitYield_AllDauPtScan_%i",ibin));
		h_PromptDs_AllDauPtScan[ibin]=h_PromptDsCal_fun(h_RawRooFitYield_AllDauPtScan[ibin],hBtoDs_AnaBin_pythiaWeight, isPbPb,Form("h_RecoNormEff_AllDauPt_%i",ibin),Form("AllDauPtScan_%i",ibin),"AllDauPtScan",1);
		h_PromptDs_AllDauPtScan[ibin]->SetName(Form("h_PromptDs_AllDauPtScan_%i",ibin));	
		f_out->cd();
		h_PromptDs_AllDauPtScan[ibin]->Write("",TObject::kOverwrite);

		h_PromptDs_AllDauPtScan2[ibin]=h_PromptDsCal_fun(h_RawRooFitYield_AllDauPtScan[ibin],hBtoDs_AnaBin_pythiaWeight, isPbPb,Form("h_RecoNormEff_AllDauPt_%i",ibin),Form("AllDauPtScan2_%i",ibin),"AllDauPtScan2");
		h_PromptDs_AllDauPtScan2[ibin]->SetName(Form("h_PromptDs_AllDauPtScan2_%i",ibin));	
		f_out->cd();
		h_PromptDs_AllDauPtScan2[ibin]->Write("",TObject::kOverwrite);

	} // end for ibin<nbin_DauPtScan
		cout<<"end for ibin<nbin_DauPtScan"<<endl;
	TH1D *h_PromptDs_AllDauPtScan_pt[nbin_pt];
	TH1D *h_PromptDs_AllDauPtScan2_pt[nbin_pt];
	for(int ibin_pt=StartBin; ibin_pt<nbin_pt; ibin_pt++){
		h_PromptDs_AllDauPtScan_pt[ibin_pt]=new TH1D(Form("h_PromptDs_AllDauPtScan_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), Form("h_PromptDs_AllDauPtScan_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), nbin_DauPtScan,bins_DauPtScan );	
		h_PromptDs_AllDauPtScan2_pt[ibin_pt]=new TH1D(Form("h_PromptDs_AllDauPtScan2_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), Form("h_PromptDs_AllDauPtScan2_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), nbin_DauPtScan,bins_DauPtScan );	
		for(int ibin=0; ibin<nbin_DauPtScan; ibin++){
//			// cout<<"AllDauPtScan_bins = "<<AllDauPtScan_bins[ibin]<<" , bins_AllDauPtScan = "<<bins_AllDauPtScan[ibin]<<endl;
			h_PromptDs_AllDauPtScan_pt[ibin_pt]->SetBinContent(ibin+1,h_PromptDs_AllDauPtScan[ibin]->GetBinContent(ibin_pt+1));
			h_PromptDs_AllDauPtScan_pt[ibin_pt]->SetBinError(ibin+1,h_PromptDs_AllDauPtScan[ibin]->GetBinError(ibin_pt+1));

			h_PromptDs_AllDauPtScan2_pt[ibin_pt]->SetBinContent(ibin+1,h_PromptDs_AllDauPtScan2[ibin]->GetBinContent(ibin_pt+1));
			h_PromptDs_AllDauPtScan2_pt[ibin_pt]->SetBinError(ibin+1,h_PromptDs_AllDauPtScan2[ibin]->GetBinError(ibin_pt+1));
		}
		f_out->cd();
		h_PromptDs_AllDauPtScan_pt[ibin_pt]->Write("",TObject::kOverwrite);
		h_PromptDs_AllDauPtScan2_pt[ibin_pt]->Write("",TObject::kOverwrite);

		c_test[count_test]=new TCanvas(Form("c_test_%i",count_test),Form("c_test_%i",count_test),800,800);
		c_test[count_test]->cd();
		h_PromptDs_AllDauPtScan_pt[ibin_pt]->Draw();

 	  shiftY=0;
  	TLatex *tl_AllDauPtScan =new TLatex();
  	tl_AllDauPtScan->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%.0f < pt < %.0f ",bins_pt[ibin_pt],bins_pt[ibin_pt+1])); shiftY-=oneshift;
  	tl_AllDauPtScan->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s AllDauPtScan ",str_PbPb.Data())); shiftY-=oneshift;

		SavePlotDirs(c_test[count_test],Form("PromptDs_%s_pt%.0fto%.0f_AllDauPtScan",str_PbPb.Data(),bins_pt[ibin_pt],bins_pt[ibin_pt+1]),{"CSdNdpt_Result","AllDauPtScan"});
		count_test++;

	} // end ibin_pt;

	// -- end AllDauPtScan

	} // end if doDauPtScan




	//-- Reschi2clScan --//
/*
	cout<<"\n\n ----- start Reschi2clScan ----- \n"<<endl;
	f_out->cd();
	TH1D *h_PromptDs_Reschi2clScan[nbin_Reschi2clScan];
	TH1D *h_PromptDs_Reschi2clScan2[nbin_Reschi2clScan];
	TH1D *h_RawRooFitYield_Reschi2clScan[nbin_Reschi2clScan];
	for(int ibin=0; ibin<nbin_Reschi2clScan; ibin++){
		cout<<"\n -- Reschi2clScan bin "<<ibin<<" --"<<endl;
		h_RawRooFitYield_Reschi2clScan[ibin]=(TH1D*)f_fitYield->Get(Form("h_RawRooFitYield_Reschi2clScan_%i",ibin));
		h_PromptDs_Reschi2clScan[ibin]=h_PromptDsCal_fun(h_RawRooFitYield_Reschi2clScan[ibin],hBtoDs_AnaBin_pythiaWeight, isPbPb,Form("h_RecoNormEff_Reschi2cl_%i",ibin),Form("Reschi2clScan_%i",ibin),"Reschi2clScan",1);
		h_PromptDs_Reschi2clScan[ibin]->SetName(Form("h_PromptDs_Reschi2clScan_%i",ibin));	
		f_out->cd();
		h_PromptDs_Reschi2clScan[ibin]->Write("",TObject::kOverwrite);

		h_PromptDs_Reschi2clScan2[ibin]=h_PromptDsCal_fun(h_RawRooFitYield_Reschi2clScan[ibin],hBtoDs_AnaBin_pythiaWeight, isPbPb,Form("h_RecoNormEff_Reschi2cl_%i",ibin),Form("Reschi2clScan2_%i",ibin),"Reschi2clScan2");
		h_PromptDs_Reschi2clScan2[ibin]->SetName(Form("h_PromptDs_Reschi2clScan2_%i",ibin));	
		f_out->cd();
		h_PromptDs_Reschi2clScan2[ibin]->Write("",TObject::kOverwrite);

	} // end for ibin<nbin_Reschi2clScan
		cout<<"end for ibin<nbin_Reschi2clScan"<<endl;
	TH1D *h_PromptDs_Reschi2clScan_pt[nbin_pt];
	TH1D *h_PromptDs_Reschi2clScan2_pt[nbin_pt];
	for(int ibin_pt=StartBin; ibin_pt<nbin_pt; ibin_pt++){
		h_PromptDs_Reschi2clScan_pt[ibin_pt]=new TH1D(Form("h_PromptDs_Reschi2clScan_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), Form("h_PromptDs_Reschi2clScan_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), nbin_Reschi2clScan,bins_Reschi2clScan );	
		h_PromptDs_Reschi2clScan2_pt[ibin_pt]=new TH1D(Form("h_PromptDs_Reschi2clScan2_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), Form("h_PromptDs_Reschi2clScan2_pt%.0fto%.0f", bins_pt[ibin_pt],bins_pt[ibin_pt+1]), nbin_Reschi2clScan,bins_Reschi2clScan );	
		for(int ibin=0; ibin<nbin_Reschi2clScan; ibin++){
			cout<<"Reschi2clScan_bins = "<<Reschi2clScan_bins[ibin]<<" , bins_Reschi2clScan = "<<bins_Reschi2clScan[ibin]<<endl;
			h_PromptDs_Reschi2clScan_pt[ibin_pt]->SetBinContent(ibin+1,h_PromptDs_Reschi2clScan[ibin]->GetBinContent(ibin_pt+1));
			h_PromptDs_Reschi2clScan_pt[ibin_pt]->SetBinError(ibin+1,h_PromptDs_Reschi2clScan[ibin]->GetBinError(ibin_pt+1));

			h_PromptDs_Reschi2clScan2_pt[ibin_pt]->SetBinContent(ibin+1,h_PromptDs_Reschi2clScan2[ibin]->GetBinContent(ibin_pt+1));
			h_PromptDs_Reschi2clScan2_pt[ibin_pt]->SetBinError(ibin+1,h_PromptDs_Reschi2clScan2[ibin]->GetBinError(ibin_pt+1));
		}
		f_out->cd();
		h_PromptDs_Reschi2clScan_pt[ibin_pt]->Write("",TObject::kOverwrite);
		h_PromptDs_Reschi2clScan2_pt[ibin_pt]->Write("",TObject::kOverwrite);

		c_test[count_test]=new TCanvas(Form("c_test_%i",count_test),Form("c_test_%i",count_test),800,800);
		c_test[count_test]->cd();
		h_PromptDs_Reschi2clScan_pt[ibin_pt]->Draw();

 	  shiftY=0;
  	TLatex *tl_Reschi2clScan =new TLatex();
  	tl_Reschi2clScan->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%.0f < pt < %.0f ",bins_pt[ibin_pt],bins_pt[ibin_pt+1])); shiftY-=oneshift;
  	tl_Reschi2clScan->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s Reschi2clScan ",str_PbPb.Data())); shiftY-=oneshift;

		SavePlotDirs(c_test[count_test],Form("PromptDs_%s_pt%.0fto%.0f_Reschi2clScan",str_PbPb.Data(),bins_pt[ibin_pt],bins_pt[ibin_pt+1]),{"CSdNdpt_Result","Reschi2clScan"});
		count_test++;

	} // end ibin_pt;

	// -- end Reschi2clScan
*/



	} // end doCutScan




	return;


}
