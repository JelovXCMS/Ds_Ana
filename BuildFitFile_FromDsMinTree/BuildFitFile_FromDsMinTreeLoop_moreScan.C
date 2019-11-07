/*
/// ideal fitter procedure, input divided by bins, read default cuts file/,
add systematic variable , variable cut array, doSystematic,

Generate TH1 & RooDataSet for latter use

fiiting , start from read MC shape (Ds, D+) with Cut

fitting method , 1. TH1 fit,(always do this first) 
2. RooFit unbin

fitting systematic options: on/off, method #

for 2,Roofit RooFit unbin : apply TH1 Fit first to get initial fit parameters , ex background & #Sig
Exclude Signal region to fit background first.
then Exclude D+ , to fit , in the end, add D+ , fit everything together
(optional , add Null hypo and get Likelihood for Significance calculation) 

output , RawYield
save fitting pdf: MC & Data 
save fitting result as txt, latter read


adding dca into ttree



// roofit save and read : https://root.cern.ch/root/html/tutorials/roofit/rf402_datahandling.C.html


// use new version of Root to avoid error , local setting CMSenvNew (6.10)

// make cuts & variation array a parameter file , ex, parameterPP  parameterPbPb30-100
copy & replace the parameter file before run, written in .sh 

*/

#include <iostream>
#include <fstream>
#include <iomanip>

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

//#include <sstream>

#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/DsMinTreeLoad.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
using namespace RooFit;


//#include "/home/peng43/work/Project/GspDiBjet2017/helpers/plotting.h"
// #include "/home/peng43/work/Project/GspDiBjet2017/helpers/InitStyle.h"

//// https://root.cern.ch/root/html/tutorials/roofit/rf402_datahandling.C.html

using namespace std;


Bool_t istest = false;
// change Bool isPbPb to int isPbPb, 1 = 30-100, 2=0-30, 3=0-100;
// REAL has no meaning here


/*
	 Int_t Selection(DsMinTreeLoad &DsClass, Int_t ibin_Dpt=0, Int_t isPbPb=0, Bool_t REAL=true, Bool_t doVarScan=false, TString Var="Dalpha", Int_t ibin_Var=0){


	 int hiBinLow=-999;
	 int hiBinHigh=999;

	 double DptLow=6;
	 double DptHigh=8;
	 double DalphaMax=0.15;
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

	 Dchi2clMin=Dchi2clMin_bins_pp[ibin_Dpt];
	 DalphaMax=DalphaMax_bins_pp[ibin_Dpt];
	 DdlsMin=DdlsMin_bins_pp[ibin_Dpt];
	 Dtrk1PtMin=Dtrk1PtMin_bins_pp[ibin_Dpt];
	 Dtrk2PtMin=Dtrk2PtMin_bins_pp[ibin_Dpt];
	 Dtrk3PtMin=Dtrk3PtMin_bins_pp[ibin_Dpt];
	 hiBinLow = -999;
	 hiBinHigh= 999;
// if(Var=="Dchi2cl"){Dchi2clMin=Dchi2clScan_bins[ibin_Var];}
}else if (isPbPb==1){
cout<<"not done yet, return 0"<<endl;
return 0;
}

if(DsClass.Dalpha< DalphaMax  &&  )	{

return 1;

}



//	}

//test could read or not
cout<<"DsMC.Dtrk1Pt = "<<DsClass.Dtrk1Pt<<endl;
cout<<"DalphaMax_bins = "<<DalphaMax_bins[ibin_Dpt]<<endl;




return 1;

}
*/

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


Int_t Selection(DsMinTreeLoad &DsClass, Int_t ibin_Dpt=0, Int_t isPbPb=0, Bool_t REAL=true,Int_t PNPrompt=0 ,Int_t DsGenTrue=23333, Bool_t doVarScan=false, TString Var="Dalpha", Int_t ibin_Var=0, Bool_t looseCut=false, Double_t tktkResmassCutWidth=0 ){

	double tktkResmasschi2cl=Reschi2clCut;	
	double PhimassCutWidth=DtktkResmassCutWidth;
	if(tktkResmassCutWidth!=0){
		PhimassCutWidth=tktkResmassCutWidth;
	}

	// cout<<"DsMC.Dtrk1Pt = "<<DsClass.Dtrk1Pt<<endl;

	int hiBinLow=-999;
	int hiBinHigh=999;

//	 the default is not used, load value in paramete.h later
	double DptLow=0;
	double DptHigh=100;
	double Dchi2clMin=0.05;
	double DalphaMax=0.15;
	double DdlsMin=2;
	double Dtrk1PtMin=0.7;
	double Dtrk2PtMin=0.7;
	double Dtrk3PtMin=0.7;


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

	}else if(isPbPb==3){
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



	}else {
		cout<<"not done PbPb==1,2 code yet, return 0"<<endl;
		return 0;
	} // end if else isPbPb=0,1,2,3

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
		if(Var=="Reschi2cl"){
			tktkResmasschi2cl=Reschi2clScan_bins[ibin_Var];
		}

		if(Var=="KaonPt"){
			Dtrk1PtMin=DauPtScan_bins[ibin_Var];
			Dtrk2PtMin=DauPtScan_bins[ibin_Var];
		}
		if(Var=="PionPt"){
			Dtrk3PtMin=DauPtScan_bins[ibin_Var];
		}
		if(Var=="AllDauPt"){
			Dtrk1PtMin=DauPtScan_bins[ibin_Var];
			Dtrk2PtMin=DauPtScan_bins[ibin_Var];
			Dtrk3PtMin=DauPtScan_bins[ibin_Var];
		}


	} // end do VarScan==true

	// start selection

	if(    DsClass.Dpt > DptLow && DsClass.Dpt < DptHigh
			&& DsClass.hiBin >= hiBinLow && DsClass.hiBin <= hiBinHigh 
			&& TMath::Abs(DsClass.DtktkResmass - DtktkResmassCutMean) < PhimassCutWidth
			&& DsClass.DtktkRes_chi2cl >  tktkResmasschi2cl
			&& DsClass.Dalpha    <   DalphaMax 
			&& DsClass.Dchi2cl   >   Dchi2clMin 
			&& DsClass.Ddls      >   DdlsMin
			&& DsClass.Dtrk1Pt   >   Dtrk1PtMin
			&& DsClass.Dtrk2Pt   >   Dtrk2PtMin
			&& DsClass.Dtrk3Pt   >   Dtrk3PtMin 
			//			&& DsClass.DsGen==DsGenTrue && DsClass.DgencollisionId==0 
			&& DsClass.Dmass > DataFitRangeLow && DsClass.Dmass<DataFitRangeHigh
		)
	{ return 1;
		// cout<<"DgenBAncestorpt = "<<DsClass.DgenBAncestorpt<<endl;
		// if(PNPrompt==0 && DsClass.DgenBAncestorpt<=0) { //cout<<"hello Prompt"<<endl;
		// return 1;}
		// else if(PNPrompt==1 && DsClass.DgenBAncestorpt>=0) { //cout<<"hello NonPrompt"<<endl; 
		// return 1;}
		// else {return 0;}
	}else{
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

double DptLow=0;
double DptHigh=100;
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

if(    DsClass.Dpt > DptLow && DsClass.Dpt < DptHigh
&& DsClass.hiBin > hiBinLow && DsClass.hiBin < hiBinHigh 
&& TMath::Abs(DsClass.DtktkResmass - DtktkResmassCutMean) < DtktkResmassCutWidth
&& DsClass.Dalpha    <   DalphaMax 
&& DsClass.Dchi2cl   >   Dchi2clMin 
&& DsClass.Ddls      >   DdlsMin
&& DsClass.Dtrk1Pt   >   Dtrk1PtMin
&& DsClass.Dtrk2Pt   >   Dtrk2PtMin
&& DsClass.Dtrk3Pt   >   Dtrk3PtMin 
//			&& DsClass.DsGen==DsGenTrue && DsClass.DgencollisionId==0 
//			TMath::Abs(DsClass.Dmass- DsMass)< DsMassRange
)
{ return 1;
// cout<<"DgenBAncestorpt = "<<DsClass.DgenBAncestorpt<<endl;
// if(PNPrompt==0 && DsClass.DgenBAncestorpt<=0) { //cout<<"hello Prompt"<<endl;
// return 1;}
// else if(PNPrompt==1 && DsClass.DgenBAncestorpt>=0) { //cout<<"hello NonPrompt"<<endl; 
// return 1;}
// else {return 0;}
}else{
return 0;
}

//test could read or not
// cout<<"DsMC.Dtrk1Pt = "<<DsClass.Dtrk1Pt<<endl;
// cout<<"DalphaMax_bins = "<<DalphaMax_bins[ibin_Dpt]<<endl;

return 0;

}
*/


int BuildFitFile_FromDsMinTreeLoop_moreScan(TString infile="", TString outfile="", Bool_t REAL=true, Int_t isPbPb=3, Int_t Ndiv=0, Int_t idiv=0, Int_t ibin_Dpt=2, Int_t FitMethod=1 , Int_t doSysFit=0, Int_t doSysCut=0, TString SysCutVar="", Int_t SysCutVarN=2, Float_t *SysCutVarBin=NULL  )
{

	initParameter();


	// TString infileData=ppDataFile;

	double *bins_pt=bins_pt_pp;
	int nbin_pt=nbin_pt_pp;
	// double TrkptAcc=TrkptAcc_pp;
	TCut cuthiBin="";

	double *Dchi2clMinScan_bins=Dchi2clMinScan_bins_pp;
	double *DalphaMaxScan_bins=DalphaMaxScan_bins_pp;
	double *DdlsMinScan_bins=DdlsMinScan_bins_pp;


	TString Str_isPbPb="pp";
	if(isPbPb==1) {
		Str_isPbPb="PbPb1";
		bins_pt=bins_pt_PbPb1;
		nbin_pt=nbin_pt_PbPb1;
		// TrkptAcc=TrkptAcc_PbPb1;
		cuthiBin=cuthiBin_PbPb1;
		Dchi2clMinScan_bins=Dchi2clMinScan_bins_PbPb1;
		DalphaMaxScan_bins=DalphaMaxScan_bins_PbPb1;
		DdlsMinScan_bins=DdlsMinScan_bins_PbPb1;

	}
	else if(isPbPb==2) {
		Str_isPbPb="PbPb2";
		bins_pt=bins_pt_PbPb2;
		nbin_pt=nbin_pt_PbPb2;
		// TrkptAcc=TrkptAcc_PbPb2;
		cuthiBin=cuthiBin_PbPb2;
		Dchi2clMinScan_bins=Dchi2clMinScan_bins_PbPb2;
		DalphaMaxScan_bins=DalphaMaxScan_bins_PbPb2;
		DdlsMinScan_bins=DdlsMinScan_bins_PbPb2;

	}
	else if(isPbPb==3) {
		Str_isPbPb="PbPb3";
		bins_pt=bins_pt_PbPb3;
		nbin_pt=nbin_pt_PbPb3;
		// TrkptAcc=TrkptAcc_PbPb3;
		cuthiBin=cuthiBin_PbPb3;
		Dchi2clMinScan_bins=Dchi2clMinScan_bins_PbPb3;
		DalphaMaxScan_bins=DalphaMaxScan_bins_PbPb3;
		DdlsMinScan_bins=DdlsMinScan_bins_PbPb3;

	  // infileData=PbPb3DataFile;

	}

	TString outfileRoot=Form("./output/FitFile_%s.root",Str_isPbPb.Data()); // test all in one, rooDataset might exceed limit in PbPb.. 
	// if(Ndiv>0 && idiv>=0){
		// outfileRoot=Form("./output/FitFile_%s_idiv%i.root",Str_isPbPb.Data(),idiv);
	// }



	TString infileData=infile;
	outfileRoot=Form("%s%s.root",outfile.Data(),s_CutSet.Data());
  if(Ndiv>0 && idiv>=0){
	  outfileRoot=Form("%s_%idiv%i.root",outfile.Data(),Ndiv,idiv);
	}
		cout<<"outfileRoot = "<<outfileRoot<<endl;

	TFile* outf = TFile::Open(Form("%s", outfileRoot.Data()),"recreate");
	/*
		 TString Str_PNPrompt="Prompt";
		 if(PNPrompt==1) {Str_PNPrompt="NonPrompt";}

		 Int_t DsGenTrue=23333;
		 Int_t GSignalTypeTrue=1;

		 TString Str_DsChannel="phikkpi";
		 if(DsChannel==1){
		 Str_DsChannel="f0kkpi" ;
		 DsGenTrue=24433;
		 GSignalTypeTrue=2;
		 }
		 */

	/*
		 if(ibin_Dpt > nbin_pt){
		 cout<<"ibin_Dpt = "<<ibin_Dpt<<" out of range of nbin_pt = "<<nbin_pt<<" , exit"<<endl;
		 return 0;		
		 }

	// use default pp setting first
	double DptLow=bins_pt[ibin_Dpt];	
	double DptHigh=bins_pt[ibin_Dpt+1];

	// load cuts & filenames
	double Dchi2clMin=Dchi2clMin_bins_pp[ibin_Dpt];
	double DalphaMax=DalphaMax_bins_pp[ibin_Dpt];
	double DdlsMin=DdlsMin_bins_pp[ibin_Dpt];

	int hiBinLow = -999;
	int hiBinHigh= 999;

	double TrkptAcc=TrkptAcc_pp;

	// add trk pt and other cuts later

	TString infileDsMC=ppMCFile[ibin_Dpt];
	*/

	// TString outfileRoot=Form("./output/FitFile_%s.root",Str_isPbPb.Data()); // test all in one, rooDataset might exceed limit in PbPb.. 

	// TString outfileDat=Form("./output/SignalFit_pp_pt%.0fto%.0f.root",DptLow,DptHigh);
	// TString outfileRoot=Form("./output/FitFile_pp_pt%.0fto%.0f.root",DptLow,DptHigh);

	// else cout<<"--- Processing - MC";
	if(isPbPb) cout<<" - PbPb centpart "<<isPbPb;
	else cout<<" - pp";
	cout<<endl;

	//	if(FfromMnt) ifname= Form("root://xrootd.rcac.purdue.edu//store",infile.Data());
	if (!TFile::Open(infileData))   { cout << " fail to open file" << endl; return 0;}
	// if (!TFile::Open(infileDsMC))   { cout << " fail to open file" << endl; return 0;}


	cout<<"hello"<<endl;
	TFile* fdata = TFile::Open(infileData);
	TTree* ntDsData =(TTree*)fdata->Get("ntDs");

	// TFile* fDsMC = TFile::Open(infileDsMC);
	// TTree* ntDsMC =(TTree*)fDsMC->Get("ntDs");

	DsMinTreeLoad DsData;
	DsData.SetBranch(ntDsData,REAL,isPbPb);	

	// DsMinTreeLoad DsMC;
	// DsMC.SetBranch(ntDsMC, false); // Bool_t REAL, MC set false to read additional varaiable 	

	// tree variables
	/*
	*/

	// TFile* outf = TFile::Open(Form("%s", outfileRoot.Data()),"recreate");

	int nbin_DsMass=100;  // could set into parameter
	double binLow_DsMass=1.91;
	double binHi_DsMass=2.11;

	// build phi mass fit later
	// int nbin_phiMass=100;
	// double binLow_phiMass=0.97;
	// double binHi_phiMass=1.13;

	// TH1F *h_phiMass= new TH1F("h_phiMass","h_phiMass",nbin_phiMass,binLow_phiMass,binHi_phiMass); h_phiMass->Sumw2();
	// TH1F *h_phiMass_unfit= new TH1F("h_phiMass_unfit","h_phiMass_unfit",nbin_phiMass,binLow_phiMass,binHi_phiMass);

	// TH1F *h_DsMassData=new TH1F(Form("h_DsMassData_pt%.0fto%.0f",DptLow,DptHigh),Form("h_DsMassData_pt%.0fto%.0f",DptLow,DptHigh) ,nbin_DsMass,binLow_DsMass,binHi_DsMass );     h_DsMassData->Sumw2();
	// TH1F *h_DsMassDsMC=new TH1F(Form("h_DsMassDsMC_pt%.0fto%.0f",DptLow,DptHigh),Form("h_DsMassDsMC_pt%.0fto%.0f",DptLow,DptHigh) ,nbin_DsMass,binLow_DsMass,binHi_DsMass );     h_DsMassDsMC->Sumw2();

	// TH1F *h_DsRawYield=new TH1F("h_DsRawYield","h_DsRawYield" , nbin_pt, bins_pt); h_DsRawYield->Sumw2();

	outf->cd();

	TH1F *h_tktkMassData[nbin_pt];
	TH1F *h_DsMassData[nbin_pt];
	TH1F *h_DsMassData_DalphaMaxScan[nbin_pt][nbin_DalphaMaxScan];
	TH1F *h_DsMassData_Dchi2clMinScan[nbin_pt][nbin_Dchi2clMinScan];
	TH1F *h_DsMassData_DdlsMinScan[nbin_pt][nbin_DdlsMinScan];
	TH1F *h_DsMassData_PhiMassScan[nbin_pt][nbin_PhiMassScan];
	TH1F *h_DsMassData_Reschi2clScan[nbin_pt][nbin_Reschi2clScan];

	outf->cd();
	TTree *t_DsMassData[nbin_pt];
	TTree *t_DsMassData_DalphaMaxScan[nbin_pt][nbin_DalphaMaxScan];
	TTree *t_DsMassData_Dchi2clMinScan[nbin_pt][nbin_Dchi2clMinScan];
	TTree *t_DsMassData_DdlsMinScan[nbin_pt][nbin_DdlsMinScan];
	TTree *t_DsMassData_PhiMassScan[nbin_pt][nbin_PhiMassScan];
	TTree *t_DsMassData_Reschi2clScan[nbin_pt][nbin_Reschi2clScan];

	
	// more scan

	TH1F  *h_DsMassData_KaonPtScan[nbin_pt][nbin_DauPtScan];
	TTree *t_DsMassData_KaonPtScan[nbin_pt][nbin_DauPtScan];

	TH1F  *h_DsMassData_PionPtScan[nbin_pt][nbin_DauPtScan];
	TTree *t_DsMassData_PionPtScan[nbin_pt][nbin_DauPtScan];

	TH1F  *h_DsMassData_AllDauPtScan[nbin_pt][nbin_DauPtScan];
	TTree *t_DsMassData_AllDauPtScan[nbin_pt][nbin_DauPtScan];


	TH1F *h_DsMassDataLoose[nbin_pt];

	Float_t Dmass;
	Float_t Ddca;

	for(int ibin_pt = 0; ibin_pt<nbin_pt; ibin_pt++){

	outf->cd();
		h_DsMassDataLoose[ibin_pt]=new TH1F(Form("h_DsMassDataLoose_pt%.0fto%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1]),Form("h_DsMassDataLoose_pt%.0fto%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1]) ,nbin_DsMass,binLow_DsMass,binHi_DsMass );   
		h_DsMassDataLoose[ibin_pt]->Sumw2();


		h_DsMassData[ibin_pt]=new TH1F(Form("h_DsMassData_pt%.0fto%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1]),Form("h_DsMassData_pt%.0fto%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1]) ,nbin_DsMass,binLow_DsMass,binHi_DsMass );   
		h_DsMassData[ibin_pt]->Sumw2();

		
		h_tktkMassData[ibin_pt]=new TH1F(Form("h_tktkMassData_pt%.0fto%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1]),Form("h_tktkMassData_pt%.0fto%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1]) ,150,0.98,1.28 );   
		h_tktkMassData[ibin_pt]->Sumw2();


		t_DsMassData[ibin_pt]=new TTree(Form("t_DsMassData_pt%.0fto%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1]),Form("t_DsMassData_pt%.0fto%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1]));  
		t_DsMassData[ibin_pt]->Branch("Dmass",&Dmass);
		t_DsMassData[ibin_pt]->Branch("Ddca",&Ddca);

		// systematic 
		for(int ibin_Dalpha = 0; ibin_Dalpha<nbin_DalphaMaxScan; ibin_Dalpha++){
			h_DsMassData_DalphaMaxScan[ibin_pt][ibin_Dalpha]=new TH1F(Form("h_DsMassData_pt%.0fto%.0f_Dalpha%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1],DalphaMaxScan_bins[ibin_Dalpha]*100),Form("h_DsMassData_pt%.0fto%.0f_Dalpha%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1],DalphaMaxScan_bins[ibin_Dalpha]*100) ,nbin_DsMass,binLow_DsMass,binHi_DsMass );   
			h_DsMassData_DalphaMaxScan[ibin_pt][ibin_Dalpha]->Sumw2();

			t_DsMassData_DalphaMaxScan[ibin_pt][ibin_Dalpha]=new TTree(Form("t_DsMassData_pt%.0fto%.0f_Dalpha%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1],DalphaMaxScan_bins[ibin_Dalpha]*100),Form("t_DsMassData_pt%.0fto%.0f_Dalpha%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1],DalphaMaxScan_bins[ibin_Dalpha]*100));  
			t_DsMassData_DalphaMaxScan[ibin_pt][ibin_Dalpha]->Branch("Dmass",&Dmass);
			t_DsMassData_DalphaMaxScan[ibin_pt][ibin_Dalpha]->Branch("Ddca",&Ddca);
		}

		for(int ibin_Ddls = 0; ibin_Ddls<nbin_DdlsMinScan; ibin_Ddls++){
			h_DsMassData_DdlsMinScan[ibin_pt][ibin_Ddls]=new TH1F(Form("h_DsMassData_pt%.0fto%.0f_Ddls%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1],DdlsMinScan_bins[ibin_Ddls]*10),Form("h_DsMassData_pt%.0fto%.0f_Ddls%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1],DdlsMinScan_bins[ibin_Ddls]*10) ,nbin_DsMass,binLow_DsMass,binHi_DsMass );   
			h_DsMassData_DdlsMinScan[ibin_pt][ibin_Ddls]->Sumw2();

			t_DsMassData_DdlsMinScan[ibin_pt][ibin_Ddls]=new TTree(Form("t_DsMassData_pt%.0fto%.0f_Ddls%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1],DdlsMinScan_bins[ibin_Ddls]*10),Form("t_DsMassData_pt%.0fto%.0f_Ddls%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1],DdlsMinScan_bins[ibin_Ddls]*10));  
			t_DsMassData_DdlsMinScan[ibin_pt][ibin_Ddls]->Branch("Dmass",&Dmass);
			t_DsMassData_DdlsMinScan[ibin_pt][ibin_Ddls]->Branch("Ddca",&Ddca);
		}

		for(int ibin_Dchi2cl = 0; ibin_Dchi2cl<nbin_Dchi2clMinScan; ibin_Dchi2cl++){
			h_DsMassData_Dchi2clMinScan[ibin_pt][ibin_Dchi2cl]=new TH1F(Form("h_DsMassData_pt%.0fto%.0f_Dchi2cl%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1],Dchi2clMinScan_bins[ibin_Dchi2cl]*100),Form("h_DsMassData_pt%.0fto%.0f_Dchi2cl%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1],Dchi2clMinScan_bins[ibin_Dchi2cl]*100) ,nbin_DsMass,binLow_DsMass,binHi_DsMass );   
			h_DsMassData_Dchi2clMinScan[ibin_pt][ibin_Dchi2cl]->Sumw2();

			t_DsMassData_Dchi2clMinScan[ibin_pt][ibin_Dchi2cl]=new TTree(Form("t_DsMassData_pt%.0fto%.0f_Dchi2cl%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1],Dchi2clMinScan_bins[ibin_Dchi2cl]*100),Form("t_DsMassData_pt%.0fto%.0f_Dchi2cl%.0f",bins_pt[ibin_pt],bins_pt[ibin_pt+1],Dchi2clMinScan_bins[ibin_Dchi2cl]*100));  
			t_DsMassData_Dchi2clMinScan[ibin_pt][ibin_Dchi2cl]->Branch("Dmass",&Dmass);
			t_DsMassData_Dchi2clMinScan[ibin_pt][ibin_Dchi2cl]->Branch("Ddca",&Ddca);
		}

		for(int ibin_PhiMass = 0; ibin_PhiMass<nbin_PhiMassScan; ibin_PhiMass++){
			h_DsMassData_PhiMassScan[ibin_pt][ibin_PhiMass]=new TH1F(Form("h_DsMassData_pt%.0fto%.0f_PhiMass_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_PhiMass),Form("h_DsMassData_pt%.0fto%.0f_PhiMass_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_PhiMass) ,nbin_DsMass,binLow_DsMass,binHi_DsMass );   
			h_DsMassData_PhiMassScan[ibin_pt][ibin_PhiMass]->Sumw2();

			t_DsMassData_PhiMassScan[ibin_pt][ibin_PhiMass]=new TTree(Form("t_DsMassData_pt%.0fto%.0f_PhiMass_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_PhiMass),Form("t_DsMassData_pt%.0fto%.0f_PhiMass_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_PhiMass));  
			t_DsMassData_PhiMassScan[ibin_pt][ibin_PhiMass]->Branch("Dmass",&Dmass);
			t_DsMassData_PhiMassScan[ibin_pt][ibin_PhiMass]->Branch("Ddca",&Ddca);
		}

		for(int ibin_Reschi2cl = 0; ibin_Reschi2cl<nbin_Reschi2clScan; ibin_Reschi2cl++){
			h_DsMassData_Reschi2clScan[ibin_pt][ibin_Reschi2cl]=new TH1F(Form("h_DsMassData_pt%.0fto%.0f_Reschi2cl_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_Reschi2cl),Form("h_DsMassData_pt%.0fto%.0f_Reschi2cl_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_Reschi2cl) ,nbin_DsMass,binLow_DsMass,binHi_DsMass );   
			h_DsMassData_Reschi2clScan[ibin_pt][ibin_Reschi2cl]->Sumw2();

			t_DsMassData_Reschi2clScan[ibin_pt][ibin_Reschi2cl]=new TTree(Form("t_DsMassData_pt%.0fto%.0f_Reschi2cl_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_Reschi2cl),Form("t_DsMassData_pt%.0fto%.0f_Reschi2cl_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_Reschi2cl));  
			t_DsMassData_Reschi2clScan[ibin_pt][ibin_Reschi2cl]->Branch("Dmass",&Dmass);
			t_DsMassData_Reschi2clScan[ibin_pt][ibin_Reschi2cl]->Branch("Ddca",&Ddca);
		}

	// more scan
	
		for(int ibin_KaonPt = 0; ibin_KaonPt<nbin_DauPtScan; ibin_KaonPt++){
			h_DsMassData_KaonPtScan[ibin_pt][ibin_KaonPt]=new TH1F(Form("h_DsMassData_pt%.0fto%.0f_KaonPt_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_KaonPt),Form("h_DsMassData_pt%.0fto%.0f_KaonPt_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_KaonPt) ,nbin_DsMass,binLow_DsMass,binHi_DsMass );   
			h_DsMassData_KaonPtScan[ibin_pt][ibin_KaonPt]->Sumw2();

			t_DsMassData_KaonPtScan[ibin_pt][ibin_KaonPt]=new TTree(Form("t_DsMassData_pt%.0fto%.0f_KaonPt_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_KaonPt),Form("t_DsMassData_pt%.0fto%.0f_KaonPt_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_KaonPt));  
			t_DsMassData_KaonPtScan[ibin_pt][ibin_KaonPt]->Branch("Dmass",&Dmass);
			t_DsMassData_KaonPtScan[ibin_pt][ibin_KaonPt]->Branch("Ddca",&Ddca);
		}
	
		for(int ibin_PionPt = 0; ibin_PionPt<nbin_DauPtScan; ibin_PionPt++){
			h_DsMassData_PionPtScan[ibin_pt][ibin_PionPt]=new TH1F(Form("h_DsMassData_pt%.0fto%.0f_PionPt_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_PionPt),Form("h_DsMassData_pt%.0fto%.0f_PionPt_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_PionPt) ,nbin_DsMass,binLow_DsMass,binHi_DsMass );   
			h_DsMassData_PionPtScan[ibin_pt][ibin_PionPt]->Sumw2();

			t_DsMassData_PionPtScan[ibin_pt][ibin_PionPt]=new TTree(Form("t_DsMassData_pt%.0fto%.0f_PionPt_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_PionPt),Form("t_DsMassData_pt%.0fto%.0f_PionPt_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_PionPt));  
			t_DsMassData_PionPtScan[ibin_pt][ibin_PionPt]->Branch("Dmass",&Dmass);
			t_DsMassData_PionPtScan[ibin_pt][ibin_PionPt]->Branch("Ddca",&Ddca);
		}
	
		for(int ibin_AllDauPt = 0; ibin_AllDauPt<nbin_DauPtScan; ibin_AllDauPt++){
			h_DsMassData_AllDauPtScan[ibin_pt][ibin_AllDauPt]=new TH1F(Form("h_DsMassData_pt%.0fto%.0f_AllDauPt_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_AllDauPt),Form("h_DsMassData_pt%.0fto%.0f_AllDauPt_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_AllDauPt) ,nbin_DsMass,binLow_DsMass,binHi_DsMass );   
			h_DsMassData_AllDauPtScan[ibin_pt][ibin_AllDauPt]->Sumw2();

			t_DsMassData_AllDauPtScan[ibin_pt][ibin_AllDauPt]=new TTree(Form("t_DsMassData_pt%.0fto%.0f_AllDauPt_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_AllDauPt),Form("t_DsMassData_pt%.0fto%.0f_AllDauPt_%i",bins_pt[ibin_pt],bins_pt[ibin_pt+1],ibin_AllDauPt));  
			t_DsMassData_AllDauPtScan[ibin_pt][ibin_AllDauPt]->Branch("Dmass",&Dmass);
			t_DsMassData_AllDauPtScan[ibin_pt][ibin_AllDauPt]->Branch("Ddca",&Ddca);
		}





	}  // end for ibin_pt declare of TH1F & TTree


	// TH1F *h_DsMassDsMC=new TH1F("h_DsMassDsMC","h_DsMassDsMC" ,nbin_DsMass,binLow_DsMass,binHi_DsMass );     h_DsMassDsMC->Sumw2();

	//	float DtktkResmassCutMean=1.01951; //move to parameter.h
	//	float DtktkResmassCutWidth=0.09;

	// float DtrkPtCut=0.7;
	// float DdlsCut=2;
	// float DalphaCut=0.2;

	// later cut would be in array with ptbin

	// if(isPbPb){
	// DtrkPtCut=1;
	// }


	//	TCut cut_phiMass=Form("DtktkResmass>%f && DtktkResmass<%f ", DtktkResmassCutMean-DtktkResmassCutWidth, DtktkResmassCutMean+DtktkResmassCutWidth);


	/*
		 RooRealVar RooDmass("RooDmass","m_{KK#pi} (GeV)",binLow_DsMass,binHi_DsMass);
	//	RooDataSet *data[nbin_pt];
	RooDataSet *RooDsData =new RooDataSet(Form("RooDsData_pt%.0fto%.0f",DptLow,DptHigh),"RooDsData",RooArgSet(RooDmass));
	RooDataSet *RooDsMC =new RooDataSet(Form("RooDsMC_pt%.0fto%.0f",DptLow,DptHigh),"RooDsMC",RooArgSet(RooDmass));
	*/
	/*
		 for(int ibin_pt=0; ibin_pt<nbin_pt; ibin_pt++){
		 data[ibin_pt]=new RooDataSet(Form("data_pt%i",ibin_pt),"data",RooArgSet(RooDmass));
		 }
		 */

	//	TCut cut_Dpt="Dpt>8 && Dpt<15";
	/*
	// loop for DsMC   // these part might be rewritten to a function
	Long64_t nentriesDsMC = ntDsMC->GetEntries();
	for(int i=0; i<nentriesDsMC; i++){
	ntDsMC->GetEntry(i);

	Selection(DsMC);

	if(i%1000==0) cout<<setw(7)<<i<<" / "<<nentriesDsMC<<endl;

	// later try to make it as function of selection, if(select(i)) then save, and the select could be as function of variable scan

	//		cout<<"DsGen = "<<DsMC.DsGen<<endl;
	if(DsMC.DsGen!=23333)	continue; // or desired code

	if(DsMC.Dtrk1Pt <DtrkPtCut || DsMC.Dtrk2Pt<DtrkPtCut ||  DsMC.Dtrk3Pt<DtrkPtCut ) continue; //change later
	if(DsMC.Ddls<DdlsMin) continue;
	if(DsMC.Dchi2cl<Dchi2clMin) continue;
	if(DsMC.Dalpha>DalphaMax) continue;
	// if(DsMC.DtktkRes_chi2cl<0.05) continue;
	if(DsMC.DtktkResmass< (DtktkResmassCutMean - DtktkResmassCutWidth)  || DsMC.DtktkResmass> (DtktkResmassCutMean + DtktkResmassCutWidth)) continue;

	if(DsMC.Dpt> DptLow && DsMC.Dpt < DptHigh){
	if(DsMC.Dmass>binLow_DsMass && DsMC.Dmass<binHi_DsMass){
	h_DsMassDsMC->Fill(DsMC.Dmass);
	RooDmass=DsMC.Dmass;
	RooDsMC->add(RooArgSet(RooDmass));
	}
	}
	}// end for DsMC loop 
	*/

	// loop for Data
	Long64_t nentries= ntDsData->GetEntries();
	// nentries=1000000; // for test
	Long64_t startEntry=0;
	Long64_t endEntry=nentries;
	Long64_t Nmin=1000000;

	outf->cd();
	if(Ndiv>0 && idiv >=0){ // further divide files to run parallel , also need to change the output name
		Long64_t NtoProcess=(nentries-Nmin )/(Ndiv-1);
		startEntry=idiv*NtoProcess + 1;
		if(idiv==0){startEntry-=1;}
		endEntry=(idiv+1)*NtoProcess;
		if(idiv==Ndiv-1){endEntry=nentries;}
	}	

	cout<<"Ndiv = "<<Ndiv<<" , idiv = "<<idiv<<endl;
	cout<<"startEntry = "<<startEntry<<" ,endEntry = "<<endEntry<<" ,nentries = "<<nentries<<endl; 

	for(Long64_t i=startEntry; i<endEntry; i++){
		ntDsData->GetEntry(i);
		if(i%200000==0) {cout<<setw(10)<<i<<" / "<<nentries<<endl;}
		Dmass=DsData.Dmass;
		Ddca=DsData.Ddca;

		for(int ibin_pt=0; ibin_pt<nbin_pt; ibin_pt++){

			// Int_t Selection(DsMinTreeLoad &DsClass, Int_t ibin_Dpt=0, Int_t isPbPb=0, Bool_t REAL=true,Int_t PNPrompt=0 ,Int_t DsGenTrue=23333, Bool_t doVarScan=false, TString Var="Dalpha", Int_t ibin_Var=0, Bool_t looseCut=false, Double_t tktkResmassCutWidth=0)
			if( Selection(DsData , ibin_pt , isPbPb, true, 0, 23333, false) == 1 ) {
				h_DsMassData[ibin_pt]->Fill(DsData.Dmass);
				t_DsMassData[ibin_pt]->Fill();	
			}
			if( Selection(DsData , ibin_pt , isPbPb, true, 0, 23333, false, "",0,true) == 1 ) {
				h_DsMassDataLoose[ibin_pt]->Fill(DsData.Dmass);	
			}

			if( Selection(DsData , ibin_pt , isPbPb, true, 0, 23333, false, "",0,true,0.1) == 1 ) {
				h_tktkMassData[ibin_pt]->Fill(DsData.DtktkResmass);	
			}

			// Int_t Selection(DsMinTreeLoad &DsClass, Int_t ibin_Dpt=0, Int_t isPbPb=0, Bool_t REAL=true,Int_t PNPrompt=0 ,Int_t DsGenTrue=23333, Bool_t doVarScan=false, TString Var="Dalpha", Int_t ibin_Var=0){

			for(int i=0; i<nbin_DalphaMaxScan; i++){
				if( Selection(DsData , ibin_pt , isPbPb, true, 0, 23333, true, "Dalpha", i) == 1 ){
					h_DsMassData_DalphaMaxScan[ibin_pt][i]->Fill(Dmass);
					t_DsMassData_DalphaMaxScan[ibin_pt][i]->Fill();
				}
			} // end for i<nbin_DalphaMaxScan

			for(int i=0; i<nbin_Dchi2clMinScan; i++){
				if( Selection(DsData , ibin_pt , isPbPb, true, 0, 23333, true, "Dchi2cl", i) == 1 ){
					h_DsMassData_Dchi2clMinScan[ibin_pt][i]->Fill(Dmass);
					t_DsMassData_Dchi2clMinScan[ibin_pt][i]->Fill();
				}
			} // end for i<nbin_Dchi2clMinScan

			for(int i=0; i<nbin_DdlsMinScan; i++){
				if( Selection(DsData , ibin_pt , isPbPb, true, 0, 23333, true, "Ddls", i) == 1 ){
					h_DsMassData_DdlsMinScan[ibin_pt][i]->Fill(Dmass);
					t_DsMassData_DdlsMinScan[ibin_pt][i]->Fill();
				}
			} // end for i<nbin_DdlsMinScan

			for(int i=0; i<nbin_PhiMassScan; i++){
				if( Selection(DsData , ibin_pt , isPbPb, true, 0, 23333, true, "PhiMass", i) == 1 ){
					h_DsMassData_PhiMassScan[ibin_pt][i]->Fill(Dmass);
					t_DsMassData_PhiMassScan[ibin_pt][i]->Fill();
				}
			} // end for i<nbin_PhiMassScan

			for(int i=0; i<nbin_Reschi2clScan; i++){
				if( Selection(DsData , ibin_pt , isPbPb, true, 0, 23333, true, "Reschi2cl", i) == 1 ){
					h_DsMassData_Reschi2clScan[ibin_pt][i]->Fill(Dmass);
					t_DsMassData_Reschi2clScan[ibin_pt][i]->Fill();
				}
			} // end for i<nbin_Reschi2clScan

			for(int i=0; i<nbin_DauPtScan; i++){
				if( Selection(DsData , ibin_pt , isPbPb, true, 0, 23333, true, "KaonPt", i) == 1 ){
					h_DsMassData_KaonPtScan[ibin_pt][i]->Fill(Dmass);
					t_DsMassData_KaonPtScan[ibin_pt][i]->Fill();
				}
			} // end for i<nbin_DauPtScan

			for(int i=0; i<nbin_DauPtScan; i++){
				if( Selection(DsData , ibin_pt , isPbPb, true, 0, 23333, true, "PionPt", i) == 1 ){
					h_DsMassData_PionPtScan[ibin_pt][i]->Fill(Dmass);
					t_DsMassData_PionPtScan[ibin_pt][i]->Fill();
				}
			} // end for i<nbin_DauPtScan

			for(int i=0; i<nbin_DauPtScan; i++){
				if( Selection(DsData , ibin_pt , isPbPb, true, 0, 23333, true, "AllDauPt", i) == 1 ){
					h_DsMassData_AllDauPtScan[ibin_pt][i]->Fill(Dmass);
					t_DsMassData_AllDauPtScan[ibin_pt][i]->Fill();
				}
			} // end for i<nbin_DauPtScan





		}// end for ibin_pt

		} // end DsData loop events

		outf->cd();
		//	outf->Write();
		for(int ibin_pt = 0; ibin_pt<nbin_pt; ibin_pt++){
			h_DsMassData[ibin_pt]->Write("",TObject::kOverwrite);
			h_DsMassDataLoose[ibin_pt]->Write("",TObject::kOverwrite);
			t_DsMassData[ibin_pt]->Write("",TObject::kOverwrite);

			h_tktkMassData[ibin_pt]->Write("",TObject::kOverwrite);

			for(int i=0; i<nbin_DalphaMaxScan; i++){	
				h_DsMassData_DalphaMaxScan[ibin_pt][i]->Write("",TObject::kOverwrite);
				t_DsMassData_DalphaMaxScan[ibin_pt][i]->Write("",TObject::kOverwrite);
			}

			for(int i=0; i<nbin_DdlsMinScan; i++){	
				h_DsMassData_DdlsMinScan[ibin_pt][i]->Write("",TObject::kOverwrite);
				t_DsMassData_DdlsMinScan[ibin_pt][i]->Write("",TObject::kOverwrite);
			}

			for(int i=0; i<nbin_Dchi2clMinScan; i++){	
				h_DsMassData_Dchi2clMinScan[ibin_pt][i]->Write("",TObject::kOverwrite);
				t_DsMassData_Dchi2clMinScan[ibin_pt][i]->Write("",TObject::kOverwrite);
			}

			for(int i=0; i<nbin_PhiMassScan; i++){	
				h_DsMassData_PhiMassScan[ibin_pt][i]->Write("",TObject::kOverwrite);
				t_DsMassData_PhiMassScan[ibin_pt][i]->Write("",TObject::kOverwrite);
			}

			for(int i=0; i<nbin_Reschi2clScan; i++){	
				h_DsMassData_Reschi2clScan[ibin_pt][i]->Write("",TObject::kOverwrite);
				t_DsMassData_Reschi2clScan[ibin_pt][i]->Write("",TObject::kOverwrite);
			}

			for(int i=0; i<nbin_DauPtScan; i++){	
				h_DsMassData_KaonPtScan[ibin_pt][i]->Write("",TObject::kOverwrite);
				t_DsMassData_KaonPtScan[ibin_pt][i]->Write("",TObject::kOverwrite);
			}

			for(int i=0; i<nbin_DauPtScan; i++){	
				h_DsMassData_PionPtScan[ibin_pt][i]->Write("",TObject::kOverwrite);
				t_DsMassData_PionPtScan[ibin_pt][i]->Write("",TObject::kOverwrite);
			}

			for(int i=0; i<nbin_DauPtScan; i++){	
				h_DsMassData_AllDauPtScan[ibin_pt][i]->Write("",TObject::kOverwrite);
				t_DsMassData_AllDauPtScan[ibin_pt][i]->Write("",TObject::kOverwrite);
			}



		} // end for ibin_pt

		//	 RooDsData->Write("",TObject::kOverwrite);
		cout<<"--- Writing finished"<<endl;
		outf->Close();

		cout<<"--- In/Output files"<<endl;
		cout<<infileData<<endl;
		cout<<outfile<<endl;
		cout<<endl;

		if(istest){
			cout<<"is test mode"<<endl;
		}
		return 0;
	}

	int main(int argc, char *argv[])
	{
		int bf_return=-1;
		if(argc==3)
		{
			bf_return=BuildFitFile_FromDsMinTreeLoop_moreScan(argv[1], argv[2]);
		}
		else if(argc==5)
		{
			cout<<"real = "<<atoi(argv[3])<<" ,isPbPb = "<<atoi(argv[4])<<endl;
			bf_return=BuildFitFile_FromDsMinTreeLoop_moreScan(argv[1], argv[2],atoi(argv[3]), atoi(argv[4]));
		}		
		else if(argc==7)
		{
			cout<<"real = "<<atoi(argv[3])<<" ,isPbPb = "<<atoi(argv[4])<<endl;
			bf_return=BuildFitFile_FromDsMinTreeLoop_moreScan(argv[1], argv[2],atoi(argv[3]), atoi(argv[4]),  atoi(argv[5]) , atoi(argv[6]));
		}
		else
		{
			std::cout << "Usage: BuildFitFile_FromDsMinTreeLoop.exe <input_collection> <output_file> <isReal> <isPbPb>" << std::endl;
			return 1;
		}

		if(bf_return==0){
			std::cout<< "finish BuildFitFile"<<endl;
		}


		return 0;
	}

