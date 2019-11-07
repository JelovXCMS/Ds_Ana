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
#include <TLine.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include <math.h>
#include "TFitter.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

int testFit(){

	TFile *fin=TFile::Open("CutScanSys_FixShape_PbPb3.root","read");
	TH1D *h1=(TH1D*)fin->Get("h_PromptDs_DdlsMinScan_ratio_pt10to20");
//	h1->Draw();
	TGraphErrors *Gr=new TGraphErrors();	


	for(int i =0; i<h1->GetNbinsX();i++){
		cout<<" i "<<h1->GetBinContent(i+1)<<" bincenter = "<<h1->GetBinCenter(i+1)<<endl;
		if(i!=0){
			Gr->SetPoint(i, h1->GetBinCenter(i+1), h1->GetBinContent(i+1));
			Gr->SetPointError(i, 0. ,h1->GetBinError(i+1));
		}
	}

	// h1->GetXaxis()->SetRangeUser(2.5,6.5);

	double cut_Ddls=4.656;

	TF1 *f1=new TF1("f1","[0]+(x-[2])*[1]");
	f1->SetLineColor(2);
	// f1->SetRange(2.5,6.5);
	f1->SetParameter(0,1);
	// f1->FixParameter(0,1);
	f1->SetParameter(1,0);
	f1->SetParameter(2,cut_Ddls);
	f1->FixParameter(2,cut_Ddls);
	f1->SetParLimits(1,-1,0);


	// h1->Fit("f1","QN0");
	// h1->Fit("f1","EMS0");
	// h1->Fit("f1","EMS0");
	// h1->Fit("f1","EMS0","",2.5,6.5);
	// h1->Fit("f1","LEMS0","",2.5,6.5);

/*
	h1->Fit("f1","LEMS0","",2,8);
	h1->Fit("f1","LEMS0I","",2,8);
	h1->Fit("f1","EMS0I","",2,8);

	f1->SetRange(0,8);

	h1->Draw();
	f1->Draw("same");

*/
	// Gr->Fit("f1","LEMS0","",2,8);
	// Gr->Fit("f1","LEMS0I","",2,8);
	// Gr->Fit("f1","LEMS0","");
	// Gr->Fit("f1","LEMS0I","");
	// Gr->Fit("f1","LEMS0I","");
	// Gr->Fit("f1","EMS0I","",2,8);
	// Gr->Fit("f1","EMS0","",0,8);

	TFitResultPtr fitR;
	for(int i =0; i<3 ; i++){
		fitR=Gr->Fit("f1","LEMS0F","",2,8); 
	}
	double x0[1]={0};
	double x0err[1];
	fitR->GetConfidenceIntervals(1,1,1,x0,x0err,0.683,false);	

	cout<<"x0 = "<<x0[0]<<" +- "<<x0err[0]*100<<endl;

	TCanvas *cg=new TCanvas("cr","cr");
	cg->cd();
	Gr->SetMinimum(0.5);
	Gr->SetMaximum(1.4);
	Gr->SetMarkerSize(0.5);
	Gr->SetMarkerStyle(21);
	// Gr->Draw();
	Gr->Draw("AP");
	f1->SetRange(0,8);
	f1->Draw("same");


	return 0;
}
