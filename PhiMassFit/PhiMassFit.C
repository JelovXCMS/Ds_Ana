#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
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

  double shiftY=-0.3;
	double oneshift=0.075;


void PhibinnedFit(TH1D *h_tktkMass,TH1D *h_PhiWidthData , Int_t ibin, Int_t isPbPb=0){
// single gaussian fit

	int nbin_pt=nbin_pt_pp;
	double *bins_pt=bins_pt_pp;
	TString str_PbPb="pp";


	if(isPbPb==3){
		nbin_pt=nbin_pt_PbPb3;
		bins_pt=bins_pt_PbPb3;
		str_PbPb="PbPb3";
	}

	double DptLow=bins_pt[ibin];
	double DptHigh=bins_pt[ibin+1];

	double PhiMassFitRangeLow=1.005;
	double PhiMassFitRangeHigh=1.035;

	TF1 *f1_phi= new TF1("f1_phi","[0]*( TMath::Gaus(x,[1],[2])) / ( sqrt(2*TMath::Pi()*[2])) + [3]*(1+[4]*x+[5]*(2*x*x-1))");
	f1_phi->SetParameter(0 , h_tktkMass->GetBinWidth(1)*h_tktkMass->Integral( h_tktkMass->FindBin(PhiMass-0.005), h_tktkMass->FindBin(PhiMass+0.005))*1/3);
	f1_phi->SetParameter(1 , PhiMass );
	f1_phi->SetParLimits(1 , PhiMass-0.002, PhiMass+0.002);
	f1_phi->SetParameter(2 , 0.0025);
	f1_phi->SetParLimits(2 , 0.0005 , 0.005);
	f1_phi->SetParameter(3 , h_tktkMass->GetBinWidth(1)*h_tktkMass->Integral( h_tktkMass->FindBin(PhiMassFitRangeLow), h_tktkMass->FindBin(PhiMassFitRangeHigh)) ) ;
	f1_phi->SetParameter(4 , 0.2);
	f1_phi->SetParameter(5 , 0 );
	
	f1_phi->SetLineColor(kRed);
	f1_phi->SetRange(PhiMassFitRangeLow, PhiMassFitRangeHigh);

	h_tktkMass->Fit("f1_phi", "QN0" ,"", PhiMassFitRangeLow, PhiMassFitRangeHigh );
	h_tktkMass->Fit("f1_phi", "QN0" ,"", PhiMassFitRangeLow, PhiMassFitRangeHigh );
	h_tktkMass->Fit("f1_phi", "L QN0" ,"", PhiMassFitRangeLow, PhiMassFitRangeHigh );
	h_tktkMass->Fit("f1_phi", "L N0" ,"", PhiMassFitRangeLow, PhiMassFitRangeHigh );
	h_tktkMass->Fit("f1_phi", "L MN0" ,"", PhiMassFitRangeLow, PhiMassFitRangeHigh );
	h_tktkMass->Fit("f1_phi", "L M IN0" ,"", PhiMassFitRangeLow, PhiMassFitRangeHigh );
	h_tktkMass->Fit("f1_phi", "L M IN0" ,"", PhiMassFitRangeLow, PhiMassFitRangeHigh );
	// h_tktkMass->Fit("f1_phi", "L EMI 0" ,"", PhiMassFitRangeLow, PhiMassFitRangeHigh );
	// h_tktkMass->Fit("f1_phi", "L EMI S0" ,"", PhiMassFitRangeLow, PhiMassFitRangeHigh );

	h_PhiWidthData->SetBinContent(ibin+1, f1_phi->GetParameter(2) );
	h_PhiWidthData->SetBinError(ibin+1, f1_phi->GetParError(2) );

	cout<<"bin i = "<<ibin<<" , Dpt : "<<DptLow<<" - "<<DptHigh<<" , PhiWidth = "<< f1_phi->GetParameter(2)<<endl;

	TCanvas *c_PhibinfitData= new TCanvas("c_PhibinfitData" , "c_PhibinfitData" ,800,800);
	c_PhibinfitData->cd();
	h_tktkMass->GetXaxis()->SetRangeUser(PhiMassFitRangeLow,PhiMassFitRangeHigh);
	h_tktkMass->SetTitle("");
	h_tktkMass->GetXaxis()->SetTitle("m_{KK} (GeV)");
	h_tktkMass->Draw();
	f1_phi->Draw("SAME");

	shiftY=-0.3;
	TLatex *tl_PhibinfitData=new TLatex();
	tl_PhibinfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s_Data_binfit",str_PbPb.Data()));  shiftY-=oneshift;
  tl_PhibinfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%.0f < pt < %.0f ",DptLow,DptHigh)); shiftY-=oneshift;
  tl_PhibinfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("#phi mean : %.5f #pm %.5f ",f1_phi->GetParameter(1) , f1_phi->GetParError(1))); shiftY-=oneshift;
  tl_PhibinfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("#phi width : %.5f #pm %.5f ",f1_phi->GetParameter(2) , f1_phi->GetParError(2))); shiftY-=oneshift;
	SavePlotDirs(c_PhibinfitData,Form("%s_Data_binfit_pt%.0fto%.0f",str_PbPb.Data(),DptLow,DptHigh),{"PhiMassFit",str_PbPb});

	delete tl_PhibinfitData;
	delete c_PhibinfitData;
	delete f1_phi;

	return;

}

int PhiMassFit(int isPbPb=3){

	  InitStyle();
		initParameter();

	gStyle->SetOptStat(0);


  int nbin_pt=nbin_pt_pp;
  double *bins_pt=bins_pt_pp;

  double DptLow;
  double DptHigh;

  TString dataName="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree/output/FitFile_pp.root";
  TString mcName="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output/MC_eff_pp_Prompt_phikkpi.root";
  TString outName=Form("output/PhiFitResult_pp.root");
  TString str_PbPb="pp";

  if(isPbPb==3){
    nbin_pt=nbin_pt_PbPb3;
    bins_pt=bins_pt_PbPb3;

    // dataName="/scratch/halstead/p/peng43/Ds_phikkpi/BuildFitFile_FromDsMinTree/DsFitFile_PbPb3_Data.root";
    dataName="/scratch/halstead/p/peng43/Ds_phikkpi/BuildFitFile_FromDsMinTree/DsFitFile_PbPb3_Data.root";
    mcName="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output/MC_eff_PbPb3_Prompt_phikkpi.root";
    outName=Form("output/PhiFitResult_PbPb3.root");
    str_PbPb="PbPb";
  }
  // may need to consider shape contribution from nonprompt later, no significant difference for PNP, phi f0 channel

  TFile *f_data=TFile::Open(dataName.Data());
  TFile *f_mc=TFile::Open(mcName.Data());

  TFile *f_out=TFile::Open(outName.Data(),"RECREATE");
  f_out->cd();

	TH1D *h_PhiWidthData=new TH1D("h_PhiWidthData","h_PhiWidthData",nbin_pt,bins_pt);
	TH1D *h_tktkMassData[nbin_pt];

	for(int i=0; i<nbin_pt; i++){
		DptLow=bins_pt[i];
		DptHigh=bins_pt[i+1];
		h_tktkMassData[i]=(TH1D*)f_data->Get(Form("h_tktkMassData_pt%.0fto%.0f",bins_pt[i],bins_pt[i+1]));	
		PhibinnedFit(h_tktkMassData[i], h_PhiWidthData, i , isPbPb );
	}

	cout<<"check y"<<endl;

	f_out->cd();
	h_PhiWidthData->Write("",TObject::kOverwrite);


	return 0;
}


