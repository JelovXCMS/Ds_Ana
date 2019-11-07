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

using namespace std;

  double shiftY=0;
	double oneshift=0.075;


// calculate effificncy in different Dpt bin with different Dpt Gen sample 
// auto combine together in this file
// must use set bin content & set bin error


// later rewrite plan : loop entry, pass to selection function,
// selection function desing: dsGen , basic selection, selection scan
// add prompt/ nonprompt/ phi, f0

Int_t findBin(double *bins_pt,int nbin_pt, double inVal){

		Int_t DptBin=-1;
	for(int ibin=0; ibin<nbin_pt; ibin++){
		if(inVal>=bins_pt[ibin] && inVal<bins_pt[ibin+1]){
			DptBin=ibin;
			break;
		}
	}
	return DptBin;

}


Int_t Selection(DsMinTreeLoad &DsClass, Int_t ibin_Dpt=0, Int_t isPbPb=0, Bool_t REAL=true,Int_t PNPrompt=0 ,Int_t DsGenTrue=23333, Bool_t doVarScan=false, TString Var="Dalpha", Int_t ibin_Var=0, Bool_t looseCut=false){


	// cout<<"DsMC.Dtrk1Pt = "<<DsClass.Dtrk1Pt<<endl;

	int hiBinLow=-999;
	int hiBinHigh=999;

	double DptLow=0;
	double DptHigh=100;
	double Dchi2clMin=0.05;
	double DalphaMax=0.15;
	double DdlsMin=2;
	double Dtrk1PtMin=0.07;
	double Dtrk2PtMin=0.07;
	double Dtrk3PtMin=0.07;
	double PhimassCutWidth=DtktkResmassCutWidth;


	int nbin_pt=0;
	double *Dchi2clMinScan_bins=Dchi2clMinScan_bins_pp;
	double *DalphaMaxScan_bins=DalphaMaxScan_bins_pp;
	double *DdlsMinScan_bins=DdlsMinScan_bins_pp;

	if(isPbPb==0){
		nbin_pt=nbin_pt_pp;
		if(ibin_Dpt>nbin_pt){
			cout<<"ibin_Dpt out of range , return 0"<<endl;
			return 0;
		}

		if(ibin_Dpt==-1){
			ibin_Dpt=findBin(bins_pt_pp , nbin_pt , DsClass.Dpt);
			// cout<<"DsClass.Dpt = "<<DsClass.Dpt<<" , ibin_Dpt = "<<ibin_Dpt<<endl;
		}else{
			DptLow=bins_pt_pp[ibin_Dpt];
			DptHigh=bins_pt_pp[ibin_Dpt+1];
		}

		hiBinLow = -999;
		hiBinHigh= 999;

		if(looseCut){
			Dchi2clMin=Dchi2clMinLoose_pp;
			DalphaMax=DalphaMaxLoose_pp;
			DdlsMin=DdlsMinLoose_pp;
			Dtrk1PtMin=DtrkPtMinLoose_pp;
			Dtrk2PtMin=DtrkPtMinLoose_pp;
			Dtrk3PtMin=DtrkPtMinLoose_pp;
		}else{

			Dchi2clMin=Dchi2clMin_bins_pp[ibin_Dpt];
			DalphaMax=DalphaMax_bins_pp[ibin_Dpt];
			DdlsMin=DdlsMin_bins_pp[ibin_Dpt];
			Dtrk1PtMin=Dtrk1PtMin_bins_pp[ibin_Dpt];
			Dtrk2PtMin=Dtrk2PtMin_bins_pp[ibin_Dpt];
			Dtrk3PtMin=Dtrk3PtMin_bins_pp[ibin_Dpt];

			Dchi2clMinScan_bins=Dchi2clMinScan_bins_pp;
			DalphaMaxScan_bins=DalphaMaxScan_bins_pp;
			DdlsMinScan_bins=DdlsMinScan_bins_pp;

		}


	}else if (isPbPb==3){
		nbin_pt=nbin_pt_PbPb3;
		if(ibin_Dpt>nbin_pt){
			cout<<"ibin_Dpt out of range , return 0"<<endl;
			return 0;
		}

		if(ibin_Dpt==-1){
      ibin_Dpt=findBin(bins_pt_PbPb3 , nbin_pt , DsClass.Dpt);
		}else{
			DptLow=bins_pt_PbPb3[ibin_Dpt];
			DptHigh=bins_pt_PbPb3[ibin_Dpt+1];
		}

		hiBinLow = hiBinLow_PbPb3;
		hiBinHigh= hiBinHigh_PbPb3;

		if(looseCut){
			Dchi2clMin=Dchi2clMinLoose_PbPb3;
			DalphaMax=DalphaMaxLoose_PbPb3;
			DdlsMin=DdlsMinLoose_PbPb3;
			Dtrk1PtMin=DtrkPtMinLoose_PbPb3;
			Dtrk2PtMin=DtrkPtMinLoose_PbPb3;
			Dtrk3PtMin=DtrkPtMinLoose_PbPb3;
		}else{

			Dchi2clMin=Dchi2clMin_bins_PbPb3[ibin_Dpt];
			DalphaMax=DalphaMax_bins_PbPb3[ibin_Dpt];
			DdlsMin=DdlsMin_bins_PbPb3[ibin_Dpt];
			Dtrk1PtMin=Dtrk1PtMin_bins_PbPb3[ibin_Dpt];
			Dtrk2PtMin=Dtrk2PtMin_bins_PbPb3[ibin_Dpt];
			Dtrk3PtMin=Dtrk3PtMin_bins_PbPb3[ibin_Dpt];

			Dchi2clMinScan_bins=Dchi2clMinScan_bins_PbPb3;
			DalphaMaxScan_bins=DalphaMaxScan_bins_PbPb3;
			DdlsMinScan_bins=DdlsMinScan_bins_PbPb3;
		}
	}else{
		cout<<"isPbPb != 0 or 3 , return 0"<<endl;
		return 0;
	}


	if(doVarScan==true){
		if(Var=="Dalpha"){
			DalphaMax=DalphaMaxScan_bins[ibin_Var];
		}
		if(Var=="Dchi2cl"){
			Dchi2clMin=Dchi2clMinScan_bins[ibin_Var];
		}
		if(Var=="Ddls"){
			DdlsMin=DdlsMinScan_bins[ibin_Var];
		}
		if(Var=="PhiMass"){
			PhimassCutWidth=PhiMassScan_bins[ibin_Var];
		}
	}

	// start selection

	if( TMath::Abs(DsClass.Dmass- DsMass)< DsMassRange
			&& DsClass.Dpt >= DptLow && DsClass.Dpt < DptHigh
			&& DsClass.hiBin >= hiBinLow && DsClass.hiBin <= hiBinHigh 
			&& TMath::Abs(DsClass.DtktkResmass - DtktkResmassCutMean) < PhimassCutWidth
			&&	DsClass.Dalpha    <   DalphaMax 
			&& DsClass.Dchi2cl   >   Dchi2clMin 
			&& DsClass.Ddls      >   DdlsMin
			&& DsClass.Dtrk1Pt   >   Dtrk1PtMin
			&& DsClass.Dtrk2Pt   >   Dtrk2PtMin
			&& DsClass.Dtrk3Pt   >   Dtrk3PtMin 
			&& DsClass.DsGen==DsGenTrue && DsClass.DgencollisionId==0 
			&& ( ( PNPrompt==0 && DsClass.DgenBAncestorpt<=0 ) || ( PNPrompt==1 && DsClass.DgenBAncestorpt>0 ) ) )
			{
			return 1; // pass cut
			}else {
			return 0;
			}

			//test could read or not
			// cout<<"DsMC.Dtrk1Pt = "<<DsClass.Dtrk1Pt<<endl;
			// cout<<"DalphaMax_bins = "<<DalphaMax_bins[ibin_Dpt]<<endl;

	return 0;

}

/*
Int_t SelectionLoose(DsMinTreeLoad &DsClass, Int_t ibin_Dpt=0, Int_t isPbPb=0, Bool_t REAL=true,Int_t PNPrompt=0 ,Int_t DsGenTrue=23333, Bool_t doVarScan=false, TString Var="Dalpha", Int_t ibin_Var=0){


	// cout<<"DsMC.Dtrk1Pt = "<<DsClass.Dtrk1Pt<<endl;

	int hiBinLow=-999;
	int hiBinHigh=999;

	double DptLow=6;
	double DptHigh=8;
	double Dchi2clMin=0.05;
	double DalphaMax=0.15;
	double DdlsMin=2;
	double Dtrk1PtMin=0.07;
	double Dtrk2PtMin=0.07;
	double Dtrk3PtMin=0.07;

	int nbin_pt=0;

	if(isPbPb==0){
		nbin_pt=nbin_pt_pp;
		if(ibin_Dpt>nbin_pt){
			cout<<"ibin_Dpt out of range , return 0"<<endl;
			return 0;
		}

		DptLow=bins_pt_pp[ibin_Dpt];
		DptHigh=bins_pt_pp[ibin_Dpt+1];

		hiBinLow = -999;
		hiBinHigh= 999;


		Dchi2clMin=Dchi2clMinLoose_pp;
		DalphaMax=DalphaMaxLoose_pp;
		DdlsMin=DdlsMinLoose_pp;
		Dtrk1PtMin=DtrkPtMinLoose_pp;
		Dtrk2PtMin=DtrkPtMinLoose_pp;
		Dtrk3PtMin=DtrkPtMinLoose_pp;
		// if(Var=="Dchi2cl"){Dchi2clMin=Dchi2clScan_bins[ibin_Var];} change the cut for variable scan
	}else if (isPbPb==1){
		cout<<"not done PbPb code yet, return 0"<<endl;
		return 0;

	}

	// start selection

	if( TMath::Abs(DsClass.Dmass- DsMass)< DsMassRange
			&& DsClass.Dpt > DptLow && DsClass.Dpt < DptHigh
			&& DsClass.hiBin > hiBinLow && DsClass.hiBin < hiBinHigh 
			&& TMath::Abs(DsClass.DtktkResmass - DtktkResmassCutMean) < DtktkResmassCutWidth
			&&	DsClass.Dalpha    <   DalphaMax 
			&& DsClass.Dchi2cl   >   Dchi2clMin 
			&& DsClass.Ddls      >   DdlsMin
			&& DsClass.Dtrk1Pt   >   Dtrk1PtMin
			&& DsClass.Dtrk2Pt   >   Dtrk2PtMin
			&& DsClass.Dtrk3Pt   >   Dtrk3PtMin 
			&& DsClass.DsGen==DsGenTrue && DsClass.DgencollisionId==0 )
	{
		// cout<<"DgenBAncestorpt = "<<DsClass.DgenBAncestorpt<<endl;
		if(PNPrompt==0 && DsClass.DgenBAncestorpt<=0) { //cout<<"hello Prompt"<<endl;
			return 1;}
		else if(PNPrompt==1 && DsClass.DgenBAncestorpt>=0) { //cout<<"hello NonPrompt"<<endl; 
			return 1;}
		else {return 0;}
	}else{
		return 0;
	}

	//test could read or not
	// cout<<"DsMC.Dtrk1Pt = "<<DsClass.Dtrk1Pt<<endl;
	// cout<<"DalphaMax_bins = "<<DalphaMax_bins[ibin_Dpt]<<endl;

	return 0;

}
*/

int MC_EffCalLoop_noNcoll(int isPbPb=0, int PNPrompt = 0, int DsChannel=0  ){

	InitStyle();
	initParameter();
	Bool_t REAL=false;

	cout<<"running , isPbPb = "<<isPbPb<<" , PNPrompt = "<<PNPrompt<<" , DsChannel = "<<DsChannel<<endl;
	// DsChannel : 0 Ds phikkpi, 1 Ds_f0kkpi
	// isPbPb :  0=pp 1=PbPb 30-100,2=0-30, 3=0-100
	// parameters later changed to use included parameter.h
	// isPbPb=0; // 1=30-100, 2=0-30

	double *bins_pt=bins_pt_pp;
	int nbin_pt=nbin_pt_pp;
	double TrkptAcc=TrkptAcc_pp;

	double *DalphaMaxScan_bins=DalphaMaxScan_bins_pp;
	double *Dchi2clMinScan_bins=Dchi2clMinScan_bins_pp;
	double *DdlsMinScan_bins=DdlsMinScan_bins_pp;

	TCut cuthiBin="";

	TString GenWeight_Pythia="weight*GptSampleWeight";
	TString GenWeight_FONLL="weight*GptSampleWeight*GenFONLLWeight";
	TString GenWeight_D0Data="weight*GptSampleWeight*GenFONLLWeight*GenD0DataWeight";
	// TString RecoWeight="weight*DgenptSampleWeight";
	if(isPbPb){
		GenWeight_Pythia="weight*GptSampleWeight*PbPbVzWeight";
		GenWeight_FONLL="weight*GptSampleWeight*PbPbVzWeight*GenFONLLRaaWeight";
		GenWeight_D0Data="weight*GptSampleWeight*PbPbVzWeight*GenFONLLWeight*GenD0DataWeight";
		// RecoWeight="weight*DgenptSampleWeight*PbPbVzWeight*Ncoll";
	}

	TString Str_isPbPb="pp";
	TString Str_PbPb="pp";
	if(isPbPb==1) {
		Str_isPbPb="PbPb1";
		bins_pt=bins_pt_PbPb1;
		nbin_pt=nbin_pt_PbPb1;
		TrkptAcc=TrkptAcc_PbPb1;
		cuthiBin=cuthiBin_PbPb1;
		DalphaMaxScan_bins=DalphaMaxScan_bins_PbPb1;
		Dchi2clMinScan_bins=Dchi2clMinScan_bins_PbPb1;
    DdlsMinScan_bins=DdlsMinScan_bins_PbPb1;
	}
	else if(isPbPb==2) {
		Str_isPbPb="PbPb2";
		bins_pt=bins_pt_PbPb2;
		nbin_pt=nbin_pt_PbPb2;
		TrkptAcc=TrkptAcc_PbPb2;
		cuthiBin=cuthiBin_PbPb2;
		DalphaMaxScan_bins=DalphaMaxScan_bins_PbPb2;
		Dchi2clMinScan_bins=Dchi2clMinScan_bins_PbPb2;
    DdlsMinScan_bins=DdlsMinScan_bins_PbPb2;
	}
	else if(isPbPb==3) {
		Str_isPbPb="PbPb3";
		Str_PbPb="PbPb 0-100%";
		bins_pt=bins_pt_PbPb3;
		nbin_pt=nbin_pt_PbPb3;
		TrkptAcc=TrkptAcc_PbPb3;
		cuthiBin=cuthiBin_PbPb3;
    DalphaMaxScan_bins=DalphaMaxScan_bins_PbPb3;
		Dchi2clMinScan_bins=Dchi2clMinScan_bins_PbPb3;
    DdlsMinScan_bins=DdlsMinScan_bins_PbPb3;
	}

	TString Str_PNPrompt="Prompt";
	if(PNPrompt==1) {Str_PNPrompt="NonPrompt";}

	Int_t DsGenTrue=23333;	
	Int_t GSignalTypeTrue=1;

	TString Str_DsChannel="phikkpi";
	if(DsChannel==1){
		Str_DsChannel="f0kkpi" ;
		DsGenTrue=24433;
		GSignalTypeTrue=2;
	}else if(DsChannel==2){
		Str_DsChannel="kstarkkpi" ;
		DsGenTrue=25544;
		GSignalTypeTrue=3;
	}


	TCut cutGenTrue=Form("GSignalType==%i && GcollisionId==0 && TMath::Abs(Gy)<1",GSignalTypeTrue);
	TCut cutGenPNPrompt="GBAncestorpt<=0";
	if(PNPrompt==1){cutGenPNPrompt="GBAncestorpt>0";}
	TCut cutGenMass=Form("TMath::Abs(Gmass-%f)<%f",DsMass,DsMassRange); // not stored
	TCut cutGenAcc=Form("Gtk2pt >%f && GRestk1pt>%f && GRestk2pt>%f && TMath::Abs(Gtk2eta)<%f && TMath::Abs(GRestk1eta)<%f && TMath::Abs(GRestk2eta)<%f  ",TrkptAcc,TrkptAcc,TrkptAcc,TrketaAcc,TrketaAcc,TrketaAcc);  // must careful using tr2, Restrk1 & Restk2

	/*
		 TString *MCFile;
		 if(isPbPb==0 && PNPrompt==0 && DsChannel==0){MCFile=MCFile_Prompt_phikkpi_pp;  }
		 if(isPbPb==0 && PNPrompt==1 && DsChannel==0){MCFile=MCFile_NonPrompt_phikkpi_pp;  }
		 if(isPbPb==0 && PNPrompt==0 && DsChannel==1){MCFile=MCFile_Prompt_f0980kkpi_pp;  }
		 if(isPbPb==0 && PNPrompt==1 && DsChannel==1){MCFile=MCFile_NonPrompt_f0980kkpi_pp;  }
		 */

	TString MCFile;
	if(isPbPb==0 && PNPrompt==0 && DsChannel==0){MCFile=MCFile_Prompt_phikkpi_pp_Merge;  }
	if(isPbPb==0 && PNPrompt==1 && DsChannel==0){MCFile=MCFile_NonPrompt_phikkpi_pp_Merge;  }
	if(isPbPb==0 && PNPrompt==0 && DsChannel==1){MCFile=MCFile_Prompt_f0980kkpi_pp_Merge;  }
	if(isPbPb==0 && PNPrompt==1 && DsChannel==1){MCFile=MCFile_NonPrompt_f0980kkpi_pp_Merge;  }

	if(isPbPb && PNPrompt==0 && DsChannel==0){MCFile=MCFile_Prompt_phikkpi_PbPb_Merge;  }
	if(isPbPb && PNPrompt==1 && DsChannel==0){MCFile=MCFile_NonPrompt_phikkpi_PbPb_Merge;  }
	if(isPbPb && PNPrompt==0 && DsChannel==1){MCFile=MCFile_Prompt_f0980kkpi_PbPb_Merge;  }
	if(isPbPb && PNPrompt==1 && DsChannel==1){MCFile=MCFile_NonPrompt_f0980kkpi_PbPb_Merge;  }


  gSystem->MakeDirectory(Form("output%s",s_CutSet.Data()));
	// int ibin_Dpt=0;
	TString outfileName=Form("./output%s/MC_eff_%s_%s_%s.root",s_CutSet.Data(),Str_isPbPb.Data(), Str_PNPrompt.Data(), Str_DsChannel.Data());

	TFile *fout= new TFile(outfileName.Data(),"RECREATE");
	fout->cd();

	TH1D *h_GenAll = new TH1D("h_GenAll","h_GenAll",nbin_pt,bins_pt); h_GenAll->Sumw2();
	TH1D *h_GenAcc = new TH1D("h_GenAcc","h_GenAcc",nbin_pt,bins_pt); h_GenAcc->Sumw2();
	TH1D *h_GenAccEff = new TH1D("h_GenAccEff","h_GenAccEff",nbin_pt,bins_pt); h_GenAccEff->Sumw2();

	TH1D *h_GenAll_Pythia = new TH1D("h_GenAll_Pythia","h_GenAll_Pythia",nbin_pt,bins_pt); h_GenAll_Pythia->Sumw2();
	TH1D *h_GenAcc_Pythia = new TH1D("h_GenAcc_Pythia","h_GenAcc_Pythia",nbin_pt,bins_pt); h_GenAcc_Pythia->Sumw2();
	TH1D *h_GenAccEff_Pythia = new TH1D("h_GenAccEff_Pythia","h_GenAccEff_Pythia",nbin_pt,bins_pt); h_GenAccEff_Pythia->Sumw2();

	TH1D *h_GenAll_FONLL = new TH1D("h_GenAll_FONLL","h_GenAll_FONLL",nbin_pt,bins_pt); h_GenAll_FONLL->Sumw2();
	TH1D *h_GenAcc_FONLL = new TH1D("h_GenAcc_FONLL","h_GenAcc_FONLL",nbin_pt,bins_pt); h_GenAcc_FONLL->Sumw2();
	TH1D *h_GenAccEff_FONLL = new TH1D("h_GenAccEff_FONLL","h_GenAccEff_FONLL",nbin_pt,bins_pt); h_GenAccEff_FONLL->Sumw2();

	TH1D *h_RecoLooseEff = new TH1D("h_RecoLooseEff","h_RecoLooseEff",nbin_pt,bins_pt); h_RecoLooseEff->Sumw2(); // loose cut
	TH1D *h_RecoNormEff = new TH1D("h_RecoNormEff","h_RecoNormEff",nbin_pt,bins_pt); h_RecoNormEff->Sumw2();

	TH1D *h_RecoNormEff_Pythia = new TH1D("h_RecoNormEff_Pythia","h_RecoNormEff_Pythia",nbin_pt,bins_pt); h_RecoNormEff_Pythia->Sumw2();
	TH1D *h_RecoNormEff_FONLL = new TH1D("h_RecoNormEff_FONLL","h_RecoNormEff_FONLL",nbin_pt,bins_pt); h_RecoNormEff_FONLL->Sumw2();


	TH1D *h_RecoLoose = new TH1D("h_RecoLoose","h_RecoLoose",nbin_pt,bins_pt); h_RecoLoose->Sumw2(); // loose cut
	TH1D *h_RecoNorm = new TH1D("h_RecoNorm","h_RecoNorm",nbin_pt,bins_pt); h_RecoNorm->Sumw2();

	TH1D *h_RecoNorm_Pythia = new TH1D("h_RecoNorm_Pythia","h_RecoNorm_Pythia",nbin_pt,bins_pt); h_RecoNorm_Pythia->Sumw2();
	TH1D *h_RecoNorm_FONLL = new TH1D("h_RecoNorm_FONLL","h_RecoNorm_FONLL",nbin_pt,bins_pt); h_RecoNorm_FONLL->Sumw2();



	TH1D *h_DsMass[nbin_pt];
	TH1D *h_DsMassLoose[nbin_pt];

	TTree *t_DsMass[nbin_pt];
	TTree *t_DsMassLoose[nbin_pt];

	TH1D *h_RecoNorm_DalphaMaxScan[nbin_DalphaMaxScan];
	TH1D *h_RecoNorm_Dchi2clMinScan[nbin_Dchi2clMinScan];
	TH1D *h_RecoNorm_DdlsMinScan[nbin_DdlsMinScan];

	TH1D *h_RecoNorm_PhiMassScan[nbin_PhiMassScan];

	TH1D *h_RecoNormEff_DalphaMaxScan[nbin_DalphaMaxScan];
	TH1D *h_RecoNormEff_Dchi2clMinScan[nbin_Dchi2clMinScan];
	TH1D *h_RecoNormEff_DdlsMinScan[nbin_DdlsMinScan];
	TH1D *h_RecoNormEff_PhiMassScan[nbin_PhiMassScan];

	TH1D *h_DsMass_DalphaMaxScan[nbin_pt][nbin_DalphaMaxScan];
	TH1D *h_DsMass_Dchi2clMinScan[nbin_pt][nbin_Dchi2clMinScan];
	TH1D *h_DsMass_DdlsMinScan[nbin_pt][nbin_DdlsMinScan];
	TH1D *h_DsMass_PhiMassScan[nbin_pt][nbin_PhiMassScan];


	TTree *t_DsMass_DalphaMaxScan[nbin_pt][nbin_DalphaMaxScan];
	TTree *t_DsMass_Dchi2clMinScan[nbin_pt][nbin_Dchi2clMinScan];
	TTree *t_DsMass_DdlsMinScan[nbin_pt][nbin_DdlsMinScan];
	TTree *t_DsMass_PhiMassScan[nbin_pt][nbin_PhiMassScan];

	Float_t Dmass;
	Float_t Ddca;
	Float_t D0DataWeight;
	Float_t FONLLWeight;
	Float_t PythiaWeight;


	int nbin_DsMass=100;
	double binLow_DsMass=1.91;
	double binHi_DsMass=2.11;

	cout<<"check1"<<endl;

	for(int ibin_pt=0 ; ibin_pt< nbin_pt; ibin_pt++){

		double DptLow=bins_pt[ibin_pt];
		double DptHigh=bins_pt[ibin_pt+1];

		h_DsMass[ibin_pt]= new TH1D(Form("h_DsMass_pt%.0fto%.0f",DptLow,DptHigh),Form("h_DsMass_pt%.0fto%.0f",DptLow,DptHigh),nbin_DsMass,binLow_DsMass,binHi_DsMass); h_DsMass[ibin_pt]->Sumw2();

		h_DsMassLoose[ibin_pt]= new TH1D(Form("h_DsMassLoose_pt%.0fto%.0f",DptLow,DptHigh),Form("h_DsMassLoose_pt%.0fto%.0f",DptLow,DptHigh),nbin_DsMass,binLow_DsMass,binHi_DsMass); h_DsMassLoose[ibin_pt]->Sumw2();

    t_DsMass[ibin_pt]=new TTree(Form("t_DsMass_pt%.0fto%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1]),Form("t_DsMass_pt%.0fto%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1]));
    t_DsMass[ibin_pt]->Branch("Dmass",&Dmass);
    t_DsMass[ibin_pt]->Branch("Ddca",&Ddca);
    t_DsMass[ibin_pt]->Branch("D0DataWeight",&D0DataWeight);

    t_DsMassLoose[ibin_pt]=new TTree(Form("t_DsMassLoose_pt%.0fto%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1]),Form("t_DsMassLoose_pt%.0fto%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1]));
    t_DsMassLoose[ibin_pt]->Branch("Dmass",&Dmass);
    t_DsMassLoose[ibin_pt]->Branch("D0DataWeight",&D0DataWeight);

// systematics
		for(int ibin_Dalpha = 0; ibin_Dalpha<nbin_DalphaMaxScan; ibin_Dalpha++){
      h_DsMass_DalphaMaxScan[ibin_pt][ibin_Dalpha]=new TH1D(Form("h_DsMass_pt%.0fto%.0f_Dalpha_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_Dalpha),Form("h_DsMass_pt%.0fto%.0f_Dalpha_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_Dalpha) ,nbin_DsMass,binLow_DsMass,binHi_DsMass );
      h_DsMass_DalphaMaxScan[ibin_pt][ibin_Dalpha]->Sumw2();

      t_DsMass_DalphaMaxScan[ibin_pt][ibin_Dalpha]=new TTree(Form("t_DsMass_pt%.0fto%.0f_Dalpha_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_Dalpha),Form("t_DsMass_pt%.0fto%.0f_Dalpha_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_Dalpha));
      t_DsMass_DalphaMaxScan[ibin_pt][ibin_Dalpha]->Branch("Dmass",&Dmass);
      t_DsMass_DalphaMaxScan[ibin_pt][ibin_Dalpha]->Branch("D0DataWeight",&D0DataWeight);
    }

    for(int ibin_Ddls = 0; ibin_Ddls<nbin_DdlsMinScan; ibin_Ddls++){
      h_DsMass_DdlsMinScan[ibin_pt][ibin_Ddls]=new TH1D(Form("h_DsMass_pt%.0fto%.0f_Ddls_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_Ddls),Form("h_DsMass_pt%.0fto%.0f_Ddls_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_Ddls) ,nbin_DsMass,binLow_DsMass,binHi_DsMass );
      h_DsMass_DdlsMinScan[ibin_pt][ibin_Ddls]->Sumw2();

      t_DsMass_DdlsMinScan[ibin_pt][ibin_Ddls]=new TTree(Form("t_DsMass_pt%.0fto%.0f_Ddls_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_Ddls),Form("t_DsMass_pt%.0fto%.0f_Ddls_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_Ddls));
      t_DsMass_DdlsMinScan[ibin_pt][ibin_Ddls]->Branch("Dmass",&Dmass);
      t_DsMass_DdlsMinScan[ibin_pt][ibin_Ddls]->Branch("D0DataWeight",&D0DataWeight);
    }

    for(int ibin_Dchi2cl = 0; ibin_Dchi2cl<nbin_Dchi2clMinScan; ibin_Dchi2cl++){
      h_DsMass_Dchi2clMinScan[ibin_pt][ibin_Dchi2cl]=new TH1D(Form("h_DsMass_pt%.0fto%.0f_Dchi2cl_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_Dchi2cl),Form("h_DsMass_pt%.0fto%.0f_Dchi2cl_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_Dchi2cl) ,nbin_DsMass,binLow_DsMass,binHi_DsMass );
      h_DsMass_Dchi2clMinScan[ibin_pt][ibin_Dchi2cl]->Sumw2();

      t_DsMass_Dchi2clMinScan[ibin_pt][ibin_Dchi2cl]=new TTree(Form("t_DsMass_pt%.0fto%.0f_Dchi2cl_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_Dchi2cl),Form("t_DsMass_pt%.0fto%.0f_Dchi2cl_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_Dchi2cl));
      t_DsMass_Dchi2clMinScan[ibin_pt][ibin_Dchi2cl]->Branch("Dmass",&Dmass);
      t_DsMass_Dchi2clMinScan[ibin_pt][ibin_Dchi2cl]->Branch("D0DataWeight",&D0DataWeight);
    }


    for(int ibin_PhiMass = 0; ibin_PhiMass<nbin_PhiMassScan; ibin_PhiMass++){
      h_DsMass_PhiMassScan[ibin_pt][ibin_PhiMass]=new TH1D(Form("h_DsMass_pt%.0fto%.0f_PhiMass_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_PhiMass),Form("h_DsMass_pt%.0fto%.0f_PhiMass_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_PhiMass) ,nbin_DsMass,binLow_DsMass,binHi_DsMass );
      h_DsMass_PhiMassScan[ibin_pt][ibin_PhiMass]->Sumw2();

      t_DsMass_PhiMassScan[ibin_pt][ibin_PhiMass]=new TTree(Form("t_DsMass_pt%.0fto%.0f_PhiMass_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_PhiMass),Form("t_DsMass_pt%.0fto%.0f_PhiMass_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_PhiMass));
      t_DsMass_PhiMassScan[ibin_pt][ibin_PhiMass]->Branch("Dmass",&Dmass);
      t_DsMass_PhiMassScan[ibin_pt][ibin_PhiMass]->Branch("D0DataWeight",&D0DataWeight);
    } // end for ibin scan Phi


	}	// end for ibin_Dpt;		 


	cout<<"check2"<<endl;

	// for systematic eff

		for(int ibin_Dalpha = 0; ibin_Dalpha<nbin_DalphaMaxScan; ibin_Dalpha++){
			h_RecoNorm_DalphaMaxScan[ibin_Dalpha]=new TH1D(Form("h_RecoNorm_Dalpha_%i",ibin_Dalpha), Form("h_RecoNorm_Dalpha_%i",ibin_Dalpha), nbin_pt, bins_pt); h_RecoNorm_DalphaMaxScan[ibin_Dalpha]->Sumw2();
			h_RecoNormEff_DalphaMaxScan[ibin_Dalpha]=new TH1D(Form("h_RecoNormEff_Dalpha_%i",ibin_Dalpha), Form("h_RecoNormEff_Dalpha_%i",ibin_Dalpha), nbin_pt, bins_pt); h_RecoNormEff_DalphaMaxScan[ibin_Dalpha]->Sumw2();
	}
    for(int ibin_Ddls = 0; ibin_Ddls<nbin_DdlsMinScan; ibin_Ddls++){
			h_RecoNorm_DdlsMinScan[ibin_Ddls]=new TH1D(Form("h_RecoNorm_Ddls_%i",ibin_Ddls), Form("h_RecoNorm_Ddls_%i",ibin_Ddls), nbin_pt, bins_pt); h_RecoNorm_DdlsMinScan[ibin_Ddls]->Sumw2();
			h_RecoNormEff_DdlsMinScan[ibin_Ddls]=new TH1D(Form("h_RecoNormEff_Ddls_%i",ibin_Ddls), Form("h_RecoNormEff_Ddls_%i",ibin_Ddls), nbin_pt, bins_pt); h_RecoNormEff_DdlsMinScan[ibin_Ddls]->Sumw2();
	}
    for(int ibin_Dchi2cl = 0; ibin_Dchi2cl<nbin_Dchi2clMinScan; ibin_Dchi2cl++){
			h_RecoNorm_Dchi2clMinScan[ibin_Dchi2cl]=new TH1D(Form("h_RecoNorm_Dchi2cl_%i",ibin_Dchi2cl), Form("h_RecoNorm_Dchi2cl_%i",ibin_Dchi2cl), nbin_pt, bins_pt); h_RecoNorm_Dchi2clMinScan[ibin_Dchi2cl]->Sumw2();
			h_RecoNormEff_Dchi2clMinScan[ibin_Dchi2cl]=new TH1D(Form("h_RecoNormEff_Dchi2cl_%i",ibin_Dchi2cl), Form("h_RecoNormEff_Dchi2cl_%i",ibin_Dchi2cl), nbin_pt, bins_pt); h_RecoNormEff_Dchi2clMinScan[ibin_Dchi2cl]->Sumw2();
	}

    for(int ibin_PhiMass = 0; ibin_PhiMass<nbin_PhiMassScan; ibin_PhiMass++){
			h_RecoNorm_PhiMassScan[ibin_PhiMass]=new TH1D(Form("h_RecoNorm_PhiMass_%i",ibin_PhiMass), Form("h_RecoNorm_PhiMass_%i",ibin_PhiMass), nbin_pt, bins_pt); h_RecoNorm_PhiMassScan[ibin_PhiMass]->Sumw2();
			h_RecoNormEff_PhiMassScan[ibin_PhiMass]=new TH1D(Form("h_RecoNormEff_PhiMass_%i",ibin_PhiMass), Form("h_RecoNormEff_PhiMass_%i",ibin_PhiMass), nbin_pt, bins_pt); h_RecoNormEff_PhiMassScan[ibin_PhiMass]->Sumw2();
	}


  cout<<"check3"<<endl;


	TString inputFile=MCFile;           

	TFile *fin=TFile::Open(inputFile.Data());

	fout->cd();
	TTree *ntDs=(TTree*)fin->Get("ntDs");
	TTree *ntGen=(TTree*)fin->Get("ntGen");

	// TCut cutDMass=Form("TMath::Abs(Dmass-%f)<%f",DsMass,DsMassRange);
	// Dy already cut to 1 in DsMinTree

	// TCut cutGenPt=Form("Gpt> %f && Gpt< %f",DptLow,DptHigh);

/*
	TCut cutGenAll=(TCut)(cutGenTrue && cutGenPNPrompt && cuthiBin)*GenWeight.Data()	;
	cout<<"cutGenAll = "<<cutGenAll<<endl;

	TCut cutGenAccAll=(TCut)(cutGenAcc && cutGenTrue && cutGenPNPrompt && cuthiBin)*GenWeight.Data()	;
	cout<<"cutGenAccAll = "<<cutGenAccAll<<endl;
*/


	ntGen->Project("h_GenAll","Gpt",(TCut)((cutGenTrue && cutGenPNPrompt && cuthiBin)*GenWeight_D0Data.Data()) ); 
	ntGen->Project("h_GenAcc","Gpt",(TCut)((cutGenTrue && cutGenAcc && cutGenPNPrompt && cuthiBin)*GenWeight_D0Data.Data()) );

	ntGen->Project("h_GenAll_Pythia","Gpt",(TCut)((cutGenTrue && cutGenPNPrompt && cuthiBin)*GenWeight_Pythia.Data()) ); 
	ntGen->Project("h_GenAcc_Pythia","Gpt",(TCut)((cutGenTrue && cutGenAcc && cutGenPNPrompt && cuthiBin)*GenWeight_Pythia.Data()) );

	ntGen->Project("h_GenAll_FONLL","Gpt",(TCut)((cutGenTrue && cutGenPNPrompt && cuthiBin)*GenWeight_FONLL.Data()) ); 
	ntGen->Project("h_GenAcc_FONLL","Gpt",(TCut)((cutGenTrue && cutGenAcc && cutGenPNPrompt && cuthiBin)*GenWeight_FONLL.Data()) );


//	ntGen->Project("h_GenAll","Gpt", "((GSignalType==1 && GcollisionId==0 && TMath::Abs(Gy)<1)&&(GBAncestorpt<=0))*weight*GptSampleWeight" ); // not work 
//	ntGen->Project("h_GenAcc","Gpt",cutGenAccAll ); // work
	//	ntGen->Draw("Gpt>>h_GenAll","((GSignalType==1 && GcollisionId==0 && TMath::Abs(Gy)<1)&&(GBAncestorpt<=0))*weight*GptSampleWeight"); // work

	DsMinTreeLoad DsMCReco;
	DsMCReco.SetBranch(ntDs, REAL, isPbPb); // false for MC, to read more branch info

	Float_t weight;
	Float_t DgenptSampleWeight;
	Float_t PbPbVzWeight;
	Float_t RecoFONLLWeight;
	Float_t RecoFONLLRaaWeight;
	Float_t RecoD0DataWeight; 


	ntDs->SetBranchAddress("weight",&weight);
	ntDs->SetBranchAddress("DgenptSampleWeight",&DgenptSampleWeight);
	ntDs->SetBranchAddress("RecoFONLLWeight",&RecoFONLLWeight);
	ntDs->SetBranchAddress("RecoD0DataWeight",&RecoD0DataWeight);
	if(isPbPb){
		ntDs->SetBranchAddress("PbPbVzWeight",&PbPbVzWeight);
	ntDs->SetBranchAddress("RecoFONLLRaaWeight",&RecoFONLLRaaWeight);
	}	

  cout<<"check4"<<endl;

	Long64_t nentries_Reco=ntDs->GetEntries();

	cout<<"nentries_Reco = "<<nentries_Reco<<endl;

	for(Long64_t ientry=0 ; ientry<nentries_Reco; ientry++){
		if(ientry%500000==0){
			cout<<setw(9)<<ientry<<" / "<<nentries_Reco<<endl;
		}
	
		ntDs->GetEntry(ientry);
		// Int_t Selection(DsMinTreeLoad &DsClass, Int_t ibin_Dpt=0, Int_t isPbPb=0, Bool_t REAL=true,Int_t PNPrompt=0 Int_t DsGenTrue=23333, Bool_t doVarScan=false, TString Var="Dalpha", Int_t ibin_Var=0, Bool_t looseCut=false)
		Dmass=DsMCReco.Dmass;
		Ddca=DsMCReco.Ddca;
		// cout<<"Dmass = "<<Dmass<<endl;
		D0DataWeight=weight*DgenptSampleWeight*RecoFONLLWeight*RecoD0DataWeight;
		FONLLWeight=weight*DgenptSampleWeight*RecoFONLLWeight;
		PythiaWeight=weight*DgenptSampleWeight;
		// cout<<"weight= "<<weight<<endl;
		// cout<<"DgenptSampleWeight = "<<DgenptSampleWeight<<endl;
		// D0DataWeight=1;

		// cout<<"D0DataWeight = "<<D0DataWeight<<endl;

		if(isPbPb){
			D0DataWeight=weight*DgenptSampleWeight*PbPbVzWeight*RecoFONLLWeight*RecoD0DataWeight;
			FONLLWeight=weight*DgenptSampleWeight*PbPbVzWeight*RecoFONLLRaaWeight;
			PythiaWeight=weight*DgenptSampleWeight*PbPbVzWeight;
		}

		bool pass_norm=false;
		bool pass_phi=false;

		// nominal eff
		if( Selection(DsMCReco , -1 , isPbPb, REAL, PNPrompt, DsGenTrue, false,"",0,false) == 1 ) {
			// cout<<"pass cut , Dpt ="<<DsMCReco.Dpt<<endl; 
			h_RecoNorm->Fill(DsMCReco.Dpt, D0DataWeight);
			h_RecoNorm_Pythia->Fill(DsMCReco.Dpt, PythiaWeight);
			h_RecoNorm_FONLL->Fill(DsMCReco.Dpt, FONLLWeight);
			pass_norm=true;

		}
		if( Selection(DsMCReco , -1 , isPbPb, REAL, PNPrompt, DsGenTrue, false,"",0,true) == 1 ) {
			h_RecoLoose->Fill(DsMCReco.Dpt, D0DataWeight);
		}

		// systemtic eff
		for(int i=0; i<nbin_DalphaMaxScan; i++){
			if( Selection(DsMCReco , -1 , isPbPb, REAL, PNPrompt, DsGenTrue, true,"Dalpha",i,false) == 1 ) {
				h_RecoNorm_DalphaMaxScan[i]->Fill(DsMCReco.Dpt, D0DataWeight);
			}
		}
		for(int i=0; i<nbin_Dchi2clMinScan; i++){
			if( Selection(DsMCReco , -1 , isPbPb, REAL, PNPrompt, DsGenTrue, true,"Dchi2cl",i,false) == 1 ) {
				h_RecoNorm_Dchi2clMinScan[i]->Fill(DsMCReco.Dpt, D0DataWeight);
			}
		}
		for(int i=0; i<nbin_DdlsMinScan; i++){
			if( Selection(DsMCReco , -1 , isPbPb, REAL, PNPrompt, DsGenTrue, true,"Ddls",i,false) == 1 ) {
				h_RecoNorm_DdlsMinScan[i]->Fill(DsMCReco.Dpt, D0DataWeight);
			}
		}
		for(int i=0; i<nbin_PhiMassScan; i++){
			if( Selection(DsMCReco , -1 , isPbPb, REAL, PNPrompt, DsGenTrue, true,"PhiMass",i,false) == 1 ) {
				h_RecoNorm_PhiMassScan[i]->Fill(DsMCReco.Dpt, D0DataWeight);
				if(i==2){pass_phi=true;}		

			}
		}

		// if(pass_norm != pass_phi){
			// cout<<"norm = "<< pass_norm<<" ,phi = "<<pass_phi<<" , PhiMass = "<<DsMCReco.DtktkResmass<<endl;
		// }
		

		// DsMass ditribution for ptbins and systematic Var bins
	  for(int ibin_pt=0 ; ibin_pt< nbin_pt; ibin_pt++){
			double DptLow=bins_pt[ibin_pt];
			double DptHigh=bins_pt[ibin_pt+1];
			if( Selection(DsMCReco , ibin_pt , isPbPb, REAL, PNPrompt, DsGenTrue, false,"",0,false) == 1 ) {
				// cout<<"pass ibin_pt"<<endl; // this works
				t_DsMass[ibin_pt]->Fill();
				h_DsMass[ibin_pt]->Fill(DsMCReco.Dmass, D0DataWeight);
			}

			if( Selection(DsMCReco , ibin_pt , isPbPb, REAL, PNPrompt, DsGenTrue, false,"",0,true) == 1 ) {
				t_DsMassLoose[ibin_pt]->Fill();
				h_DsMassLoose[ibin_pt]->Fill(DsMCReco.Dmass, D0DataWeight);  			
			}
			// systematic var scan 
			
			for(int i=0; i<nbin_DalphaMaxScan; i++){
				if( Selection(DsMCReco , ibin_pt , isPbPb, REAL, PNPrompt, DsGenTrue, true,"Dalpha",i,false) == 1 ) {
					t_DsMass_DalphaMaxScan[ibin_pt][i]->Fill();	
					h_DsMass_DalphaMaxScan[ibin_pt][i]->Fill(DsMCReco.Dmass, D0DataWeight);
				}
			}

			for(int i=0; i<nbin_Dchi2clMinScan; i++){
				if( Selection(DsMCReco , ibin_pt , isPbPb, REAL, PNPrompt, DsGenTrue, true,"Dchi2cl",i,false) == 1 ) {
					t_DsMass_Dchi2clMinScan[ibin_pt][i]->Fill();	
					h_DsMass_Dchi2clMinScan[ibin_pt][i]->Fill(DsMCReco.Dmass, D0DataWeight);
				}
			}

			for(int i=0; i<nbin_DdlsMinScan; i++){
				if( Selection(DsMCReco , ibin_pt , isPbPb, REAL, PNPrompt, DsGenTrue, true,"Ddls",i,false) == 1 ) {
					t_DsMass_DdlsMinScan[ibin_pt][i]->Fill();	
					h_DsMass_DdlsMinScan[ibin_pt][i]->Fill(DsMCReco.Dmass, D0DataWeight);
				}
			}

			for(int i=0; i<nbin_PhiMassScan; i++){
				if( Selection(DsMCReco , ibin_pt , isPbPb, REAL, PNPrompt, DsGenTrue, true,"PhiMass",i,false) == 1 ) {
					t_DsMass_PhiMassScan[ibin_pt][i]->Fill();	
					h_DsMass_PhiMassScan[ibin_pt][i]->Fill(DsMCReco.Dmass, D0DataWeight);
				}
			}



		}// end for ibin_pt

	} // end loop ientry of DsMCReco

  cout<<"check5"<<endl;
	

	// Eff Calculation & write to output

	fout->cd();

	h_GenAll->Write("",TObject::kOverwrite);
	h_GenAcc->Write("",TObject::kOverwrite);
	h_RecoNorm->Write("",TObject::kOverwrite);
	h_RecoLoose->Write("",TObject::kOverwrite);

	h_GenAll_FONLL->Write("",TObject::kOverwrite);
	h_GenAcc_FONLL->Write("",TObject::kOverwrite);
	h_RecoNorm_FONLL->Write("",TObject::kOverwrite);

	h_GenAll_Pythia->Write("",TObject::kOverwrite);
	h_GenAcc_Pythia->Write("",TObject::kOverwrite);
	h_RecoNorm_Pythia->Write("",TObject::kOverwrite);


	h_GenAccEff->Divide(h_GenAcc, h_GenAll, 1,1,"B");     h_GenAccEff->SetTitle("h_GenAccEff"); h_GenAccEff->Write("",TObject::kOverwrite);
	h_RecoNormEff->Divide(h_RecoNorm, h_GenAll, 1,1,"B"); h_RecoNormEff->SetTitle("h_RecoNormEff"); h_RecoNormEff->Write("",TObject::kOverwrite);


	h_GenAccEff_FONLL->Divide(h_GenAcc_FONLL, h_GenAll_FONLL, 1,1,"B");     h_GenAccEff_FONLL->SetTitle("h_GenAccEff_FONLL"); h_GenAccEff_FONLL->Write("",TObject::kOverwrite);
	h_RecoNormEff_FONLL->Divide(h_RecoNorm_FONLL, h_GenAll_FONLL, 1,1,"B"); h_RecoNormEff_FONLL->SetTitle("h_RecoNormEff_FONLL"); h_RecoNormEff_FONLL->Write("",TObject::kOverwrite);

	h_GenAccEff_Pythia->Divide(h_GenAcc_Pythia, h_GenAll_Pythia, 1,1,"B");     h_GenAccEff_Pythia->SetTitle("h_GenAccEff_Pythia"); h_GenAccEff_Pythia->Write("",TObject::kOverwrite);
	h_RecoNormEff_Pythia->Divide(h_RecoNorm_Pythia, h_GenAll_Pythia, 1,1,"B"); h_RecoNormEff_Pythia->SetTitle("h_RecoNormEff_Pythia"); h_RecoNormEff_Pythia->Write("",TObject::kOverwrite);



	h_RecoLooseEff->Divide(h_RecoLoose, h_GenAll, 1,1,"B"); h_RecoLooseEff->SetTitle("h_RecoLooseEff"); h_RecoLooseEff->Write("",TObject::kOverwrite);

	gStyle->SetOptStat(0);

	TCanvas *c_effNorm=new TCanvas("c_effNorm","c_effNorm",800,800);
	c_effNorm->cd();
	h_RecoNormEff->SetTitle("");
	h_RecoNormEff->GetXaxis()->SetTitle("D_{S} p_{T} GeV");
	h_RecoNormEff->GetYaxis()->SetTitle("Efficiency");
	h_RecoNormEff->Draw();
	
	shiftY=-0.3;
	TLatex *tl_eff= new TLatex();
	tl_eff->DrawLatexNDC(textposx+0.30,textposy+shiftY,Form("%s",Str_PbPb.Data()));  shiftY-=oneshift;
	tl_eff->DrawLatexNDC(textposx+0.30,textposy+shiftY,Form("%s Ds to %s",Str_PNPrompt.Data(), Str_DsChannel.Data() ));  shiftY-=oneshift;

	c_effNorm->SaveAs(Form("./plots/%s/EffNorm_%s_%s_%s.pdf",Str_isPbPb.Data(),Str_isPbPb.Data(),Str_PNPrompt.Data(),Str_DsChannel.Data()));

	TCanvas *c_effNormCompare = new TCanvas("c_effNormCompare","c_effNormCompare",800,800);
	c_effNormCompare->cd();
  h_RecoNormEff->Draw();
	h_RecoNormEff_FONLL->SetLineColor(2);
	h_RecoNormEff_FONLL->SetMarkerColor(2);
	h_RecoNormEff_FONLL->Draw("SAME");
	h_RecoNormEff_Pythia->SetLineColor(4);
	h_RecoNormEff_Pythia->SetMarkerColor(4);
	h_RecoNormEff_Pythia->Draw("SAME");
	gPad->BuildLegend();
	c_effNormCompare->SaveAs(Form("./plots/%s/EffNormCompare_%s_%s_%s.pdf",Str_isPbPb.Data(),Str_isPbPb.Data(),Str_PNPrompt.Data(),Str_DsChannel.Data()));


	TCanvas *c_AccEff=new TCanvas("c_AccEff","c_AccEff",800,800);
	c_AccEff->cd();
	h_GenAccEff->Draw();
	c_AccEff->SaveAs(Form("./plots/%s/AccEff_%s_%s_%s.pdf",Str_isPbPb.Data(),Str_isPbPb.Data(),Str_PNPrompt.Data(),Str_DsChannel.Data()));

  cout<<"check6"<<endl;

	for(int ibin_pt=0; ibin_pt<nbin_pt; ibin_pt++){
		h_DsMass[ibin_pt]->Write("",TObject::kOverwrite);
		h_DsMassLoose[ibin_pt]->Write("",TObject::kOverwrite);

		t_DsMass[ibin_pt]->Write("",TObject::kOverwrite);
		t_DsMassLoose[ibin_pt]->Write("",TObject::kOverwrite);
	}

  cout<<"check7"<<endl;
	

	TCanvas *c_eff_DalphaMaxScan=new TCanvas("c_eff_DalphaMaxScan","c_eff_DalphaMaxScan",1200,800);
	c_eff_DalphaMaxScan->Divide(3,2);
	for(int i=0; i<nbin_DalphaMaxScan; i++){
		h_RecoNormEff_DalphaMaxScan[i]->Divide(h_RecoNorm_DalphaMaxScan[i],h_GenAll, 1,1 ,"B"); h_RecoNormEff_DalphaMaxScan[i]->SetTitle(Form("h_RecoNormEff_Dalpha%.0f",DalphaMaxScan_bins[i]*100)); 
		h_RecoNormEff_DalphaMaxScan[i]->Write("",TObject::kOverwrite);
		h_RecoNorm_DalphaMaxScan[i]->Write("",TObject::kOverwrite);

		c_eff_DalphaMaxScan->cd(i+1);
		h_RecoNormEff_DalphaMaxScan[i]->Draw();
		for(int ibin_pt=0; ibin_pt<nbin_pt; ibin_pt++){
			h_DsMass_DalphaMaxScan[ibin_pt][i]->Write("",TObject::kOverwrite);	
			t_DsMass_DalphaMaxScan[ibin_pt][i]->Write("",TObject::kOverwrite);	
		}
	}
	c_eff_DalphaMaxScan->SaveAs(Form("./plots/%s/CutScan/Eff_DalphaMaxScan_%s_%s_%s.pdf",Str_isPbPb.Data(),Str_isPbPb.Data(),Str_PNPrompt.Data(),Str_DsChannel.Data()));

  cout<<"check8"<<endl;

	TCanvas *c_eff_Dchi2clMinScan=new TCanvas("c_eff_Dchi2clMinScan","c_eff_Dchi2clMinScan",1200,800);
	c_eff_Dchi2clMinScan->Divide(3,2);
	for(int i=0; i<nbin_Dchi2clMinScan; i++){
		h_RecoNormEff_Dchi2clMinScan[i]->Divide(h_RecoNorm_Dchi2clMinScan[i],h_GenAll, 1,1 ,"B"); h_RecoNormEff_Dchi2clMinScan[i]->SetTitle(Form("h_RecoNormEff_Dchi2cl%.0f",Dchi2clMinScan_bins[i]*10)); 
		h_RecoNormEff_Dchi2clMinScan[i]->Write("",TObject::kOverwrite);
		h_RecoNorm_Dchi2clMinScan[i]->Write("",TObject::kOverwrite);
		
		c_eff_Dchi2clMinScan->cd(i+1);
		h_RecoNormEff_Dchi2clMinScan[i]->Draw();
		for(int ibin_pt=0; ibin_pt<nbin_pt; ibin_pt++){
			h_DsMass_Dchi2clMinScan[ibin_pt][i]->Write("",TObject::kOverwrite);	
			t_DsMass_Dchi2clMinScan[ibin_pt][i]->Write("",TObject::kOverwrite);	
		}
	}
	c_eff_Dchi2clMinScan->SaveAs(Form("./plots/%s/CutScan/Eff_Dchi2clMinScan_%s_%s_%s.pdf",Str_isPbPb.Data(),Str_isPbPb.Data(),Str_PNPrompt.Data(),Str_DsChannel.Data()));

  cout<<"check9"<<endl;

	TCanvas *c_eff_DdlsMinScan=new TCanvas("c_eff_DdlsMinScan","c_eff_DdlsMinScan",1200,800);
	c_eff_DdlsMinScan->Divide(3,2);
	for(int i=0; i<nbin_DdlsMinScan; i++){
		h_RecoNormEff_DdlsMinScan[i]->Divide(h_RecoNorm_DdlsMinScan[i],h_GenAll, 1,1 ,"B"); h_RecoNormEff_DdlsMinScan[i]->SetTitle(Form("h_RecoNormEff_Ddls%.0f",DdlsMinScan_bins[i]*100)); 
		h_RecoNormEff_DdlsMinScan[i]->Write("",TObject::kOverwrite);
		h_RecoNorm_DdlsMinScan[i]->Write("",TObject::kOverwrite);
		c_eff_DdlsMinScan->cd(i+1);
		h_RecoNormEff_DdlsMinScan[i]->Draw();
		for(int ibin_pt=0; ibin_pt<nbin_pt; ibin_pt++){
			h_DsMass_DdlsMinScan[ibin_pt][i]->Write("",TObject::kOverwrite);	
			t_DsMass_DdlsMinScan[ibin_pt][i]->Write("",TObject::kOverwrite);	
		}
	}
	c_eff_DdlsMinScan->SaveAs(Form("./plots/%s/CutScan/Eff_DdlsMinScan_%s_%s_%s.pdf",Str_isPbPb.Data(),Str_isPbPb.Data(),Str_PNPrompt.Data(),Str_DsChannel.Data()));


	TCanvas *c_eff_PhiMassScan=new TCanvas("c_eff_PhiMassScan","c_eff_PhiMassScan",1200,800);
	c_eff_PhiMassScan->Divide(3,2);
	for(int i=0; i<nbin_PhiMassScan; i++){
		h_RecoNormEff_PhiMassScan[i]->Divide(h_RecoNorm_PhiMassScan[i],h_GenAll, 1,1 ,"B"); h_RecoNormEff_PhiMassScan[i]->SetTitle(Form("h_RecoNormEff_PhiMass_%i",i)); 
		h_RecoNormEff_PhiMassScan[i]->Write("",TObject::kOverwrite);
		h_RecoNorm_PhiMassScan[i]->Write("",TObject::kOverwrite);
		c_eff_PhiMassScan->cd(i+1);
		h_RecoNormEff_PhiMassScan[i]->Draw();
		for(int ibin_pt=0; ibin_pt<nbin_pt; ibin_pt++){
			h_DsMass_PhiMassScan[ibin_pt][i]->Write("",TObject::kOverwrite);	
			t_DsMass_PhiMassScan[ibin_pt][i]->Write("",TObject::kOverwrite);	
		}
	}
	c_eff_PhiMassScan->SaveAs(Form("./plots/%s/CutScan/Eff_PhiMassScan_%s_%s_%s.pdf",Str_isPbPb.Data(),Str_isPbPb.Data(),Str_PNPrompt.Data(),Str_DsChannel.Data()));


//
	// cout<<"DtktkResmassCutWidth = "<<DtktkResmassCutWidth<<endl;

	cout<<"running , isPbPb = "<<isPbPb<<" , PNPrompt = "<<PNPrompt<<" , DsChannel = "<<DsChannel<<endl;
	return 1;

}
