// var the f0 shape to check the efficiency change

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

#include "TFitter.h"
#include "TFitResult.h"

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
#include <RooWorkspace.h>
#include <RooMCStudy.h>
#include <RooBreitWigner.h>
#include <RooVoigtian.h>



double KAON_Mass=0.493677;

double f0_MassMin=2*KAON_Mass;

double MaxBinVal=0.99;
// using namespace RooFit;
using namespace std;
using namespace RooFit;


Double_t mybwfun(Double_t* x, Double_t* par)
{
	Double_t arg1 = 14.0/22.0; // 2 over pi
	Double_t arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma=par[1]  M=par[2]
	Double_t arg3 = ((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
	Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[1]*par[1])/(par[2]*par[2]));
	return par[0]*arg1*arg2/(arg3 + arg4);
}



Double_t rightExp(Double_t *x, Double_t *par){

	if ( x[0] >  -MaxBinVal){
		return par[0]*TMath::Exp(-par[1]*(x[0]-MaxBinVal));
	}else if(x[0]>=f0_MassMin){
		return par[2]*(x[0]-f0_MassMin);
	}else {
		return 0;
	}

}



int f0ShapeVar(){

	// double f0_MassMin=2*KAON_Mass;

	InitStyle();

	int nbin=144;
	double binLow=0.97;
	double binHigh=1.23;

	double binLow_test=0.99;

	cout<<"check 0"<<endl;

	TFile *fin=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_MergeWeight/root_output/DsMinTree_MC_GenSampleMerge_pp_Prompt_f0980.root");
	TTree *tin=(TTree*)fin->Get("ntDs");
	TH1D *h_f0mass=new TH1D("h_f0mass","h_f0mass",nbin,binLow,binHigh);

	tin->Project("h_f0mass","DtktkResmass","(DsGen==24433 && DgencollisionId==0 && Dgenpt>3 && Dgenpt<10 )*(weight*DgenptSampleWeight)");

	double Fr_ori=h_f0mass->Integral( h_f0mass->GetXaxis()->FindBin(DtktkResmassCutMean-DtktkResmassCutWidth), h_f0mass->GetXaxis()->FindBin(DtktkResmassCutMean + DtktkResmassCutWidth) ) / h_f0mass->Integral() ;


	TFile *f_980w10=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/pp_MC_DsMinTree/DsMinTree_pp_MC_Ds_f0m980w10kkpi_pt4.root");
	TTree *t_980w10=(TTree*)f_980w10->Get("ntDs");
	TH1D *h_f0mass_980w10=new TH1D("h_f0mass_980w10","h_f0mass_980w10",nbin,binLow,binHigh);

	t_980w10->Project("h_f0mass_980w10","DtktkResmass","(DsGen==24433 && DgencollisionId==0 && Dgenpt>3 && Dgenpt<10 )");

	cout<<"check 1"<<endl;

	TFile *f_990w100=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/pp_MC_DsMinTree/DsMinTree_pp_MC_Ds_f0m990w100kkpi_pt4.root");
	TTree *t_990w100=(TTree*)f_990w100->Get("ntDs");
	TH1D *h_f0mass_990w100=new TH1D("h_f0mass_990w100","h_f0mass_990w100",nbin,binLow,binHigh);

	t_990w100->Project("h_f0mass_990w100","DtktkResmass","(DsGen==24433 && DgencollisionId==0 && Dgenpt>3 && Dgenpt<10 )");



	int MaxBin=h_f0mass->GetMaximumBin();
	MaxBinVal=h_f0mass->GetBinLowEdge(MaxBin);
	// MaxBinVal=h_f0mass->GetBinCenter(MaxBin);

	double fitbinHigh=1.10;

	//Roofit test
	RooRealVar x("x","x",binLow,binHigh); 
	RooDataHist dh("dh","dh",x,Import(*h_f0mass));

	RooPlot *frame = x.frame();
	dh.plotOn(frame);


	// end roofit test

	TF1 *f1_bw=new TF1("f1_bw",mybwfun,MaxBinVal,fitbinHigh,3);
	f1_bw->SetParameters(50,0.07,0.992);
	f1_bw->SetParLimits(1,0.000001,0.15);
	f1_bw->SetParLimits(2,0.975,1.02);
	f1_bw->SetParName(0,"const");
	f1_bw->SetParName(1,"sigma");
	f1_bw->SetParName(2,"mean");
	h_f0mass->Fit("f1_bw","QNO","",MaxBinVal,fitbinHigh);
	h_f0mass->Fit("f1_bw","QNO","",MaxBinVal,fitbinHigh);
	h_f0mass->Fit("f1_bw","LQNO","",MaxBinVal,fitbinHigh);
	h_f0mass->Fit("f1_bw","L MQNO","",MaxBinVal,fitbinHigh);
	h_f0mass->Fit("f1_bw","L M INO","",MaxBinVal,fitbinHigh);
	h_f0mass->Fit("f1_bw","L M IO","",MaxBinVal,fitbinHigh);

	f1_bw->SetRange(MaxBinVal,fitbinHigh);


	TF1 *f1_pol2=new TF1("f1_pol2","[0]+[1]*x+[2]*x*x+[3]*x*x*x");
	/*
		 h_f0mass->Fit("f1_pol2","QNO","",binLow,binHigh);
		 h_f0mass->Fit("f1_pol2","QNO","",binLow,binHigh);
		 h_f0mass->Fit("f1_pol2","LQNO","",binLow,binHigh);
		 h_f0mass->Fit("f1_pol2","L MQNO","",binLow,binHigh);
		 h_f0mass->Fit("f1_pol2","L M INO","",binLow,binHigh);
		 h_f0mass->Fit("f1_pol2","L M IO","",binLow,binHigh);
		 */
	// f1_pol2->Draw("SAME");

	TF1 *f1_exp=new TF1("f1_exp","[0]+[1]*TMath::Exp(-[2]*(x-[3]))");
	f1_exp->SetParameters(0,1,1,MaxBinVal);
	f1_exp->FixParameter(0,0);
	f1_exp->FixParameter(3,MaxBinVal);

	/*
		 h_f0mass->Fit("f1_exp","QNO","",MaxBinVal,fitbinHigh);
		 h_f0mass->Fit("f1_exp","QNO","",MaxBinVal,fitbinHigh);
		 h_f0mass->Fit("f1_exp","LQNO","",MaxBinVal,fitbinHigh);
		 h_f0mass->Fit("f1_exp","L MQNO","",MaxBinVal,fitbinHigh);
		 h_f0mass->Fit("f1_exp","L M INO","",MaxBinVal,fitbinHigh);
		 h_f0mass->Fit("f1_exp","L M IO","",MaxBinVal,fitbinHigh);
		 */
	f1_exp->SetRange(MaxBinVal,fitbinHigh);


	cout<<"\n\n ---- end f1_exp ---- \n\n "<<endl;


	int MaxBin_980w10=h_f0mass_980w10->GetMaximumBin();
	double MaxBinVal_980w10=h_f0mass_980w10->GetBinLowEdge(MaxBin_980w10);
	// double MaxBinVal_980w10=h_f0mass_980w10->GetBinCenter(MaxBin_980w10);

	TF1 *f1_exp2=new TF1("f1_exp2","[0]+[1]*TMath::Exp(-[2]*(x-[3]))");
	f1_exp2->SetParameters(0,1,1,MaxBinVal_980w10);
	f1_exp2->FixParameter(0,0);
	f1_exp2->FixParameter(3,MaxBinVal_980w10);
	/*
		 h_f0mass_980w10->Fit("f1_exp2","QNO","",MaxBinVal_980w10,fitbinHigh);
		 h_f0mass_980w10->Fit("f1_exp2","QNO","",MaxBinVal_980w10,fitbinHigh);
		 h_f0mass_980w10->Fit("f1_exp2","LQNO","",MaxBinVal_980w10,fitbinHigh);
		 h_f0mass_980w10->Fit("f1_exp2","L MQNO","",MaxBinVal_980w10,fitbinHigh);
		 h_f0mass_980w10->Fit("f1_exp2","L M INO","",MaxBinVal_980w10,fitbinHigh);
		 h_f0mass_980w10->Fit("f1_exp2","L M IO","",MaxBinVal_980w10,fitbinHigh);
		 */
	f1_exp2->SetRange(MaxBinVal_980w10,fitbinHigh);


	TF1 *f1_bw2=new TF1("f1_bw2",mybwfun,MaxBinVal_980w10,fitbinHigh,3);
	f1_bw2->SetParameters(50,0.07,0.990);
	f1_bw2->SetParLimits(1,0.000001,0.15);
	f1_bw2->SetParLimits(2,0.975,1.02);
	h_f0mass_980w10->Fit("f1_bw2","QNO","",MaxBinVal_980w10,fitbinHigh);
	h_f0mass_980w10->Fit("f1_bw2","QNO","",MaxBinVal_980w10,fitbinHigh);
	h_f0mass_980w10->Fit("f1_bw2","LQNO","",MaxBinVal_980w10,fitbinHigh);
	h_f0mass_980w10->Fit("f1_bw2","L MQNO","",MaxBinVal_980w10,fitbinHigh);
	h_f0mass_980w10->Fit("f1_bw2","L M INO","",MaxBinVal_980w10,fitbinHigh);
	h_f0mass_980w10->Fit("f1_bw2","L M IO","",MaxBinVal_980w10,fitbinHigh);

	f1_bw2->SetRange(MaxBinVal_980w10,fitbinHigh);


	int MaxBin_990w100=h_f0mass_990w100->GetMaximumBin();
	double MaxBinVal_990w100=h_f0mass_990w100->GetBinLowEdge(MaxBin_990w100);
	// double MaxBinVal_990w100=h_f0mass_990w100->GetBinCenter(MaxBin_990w100);

	TF1 *f1_exp3=new TF1("f1_exp3","[0]+[1]*TMath::Exp(-[2]*(x-[3]))");
	f1_exp3->SetParameters(0,1,1,MaxBinVal_990w100);
	f1_exp3->FixParameter(0,0);
	f1_exp3->FixParameter(3,MaxBinVal_990w100);
	/*
		 h_f0mass_990w100->Fit("f1_exp3","QNO","",MaxBinVal_990w100,fitbinHigh);
		 h_f0mass_990w100->Fit("f1_exp3","QNO","",MaxBinVal_990w100,fitbinHigh);
		 h_f0mass_990w100->Fit("f1_exp3","LQNO","",MaxBinVal_990w100,fitbinHigh);
		 h_f0mass_990w100->Fit("f1_exp3","L MQNO","",MaxBinVal_990w100,fitbinHigh);
		 h_f0mass_990w100->Fit("f1_exp3","L M INO","",MaxBinVal_990w100,fitbinHigh);
		 h_f0mass_990w100->Fit("f1_exp3","L M IO","",MaxBinVal_990w100,fitbinHigh);
		 */
	f1_exp3->SetRange(MaxBinVal_990w100,fitbinHigh);


	TF1 *f1_bw3=new TF1("f1_bw3",mybwfun,MaxBinVal_990w100,fitbinHigh,3);
	f1_bw3->SetParameters(50,0.07,0.990);
	f1_bw3->SetParLimits(1,0.000001,0.15);
	f1_bw3->SetParLimits(2,0.975,1.02);
	h_f0mass_990w100->Fit("f1_bw3","QNO","",MaxBinVal_990w100,fitbinHigh);
	h_f0mass_990w100->Fit("f1_bw3","QNO","",MaxBinVal_990w100,fitbinHigh);
	h_f0mass_990w100->Fit("f1_bw3","LQNO","",MaxBinVal_990w100,fitbinHigh);
	h_f0mass_990w100->Fit("f1_bw3","L MQNO","",MaxBinVal_990w100,fitbinHigh);
	h_f0mass_990w100->Fit("f1_bw3","L M INO","",MaxBinVal_990w100,fitbinHigh);
	h_f0mass_990w100->Fit("f1_bw3","L M IO","",MaxBinVal_990w100,fitbinHigh);

	f1_bw3->SetRange(MaxBinVal_990w100,fitbinHigh);


	gStyle->SetOptStat(0);

	TCanvas *c_test=new TCanvas("c_test","c_test",800,800);
	c_test->cd();

	h_f0mass_980w10->SetTitle("");
	h_f0mass_980w10->GetXaxis()->SetRangeUser(0.97,1.15);
	h_f0mass_980w10->GetXaxis()->SetTitle("m_{KK}");
  h_f0mass_980w10->SetMarkerColor(4);	
	h_f0mass_980w10->SetLineColor(4);
	h_f0mass_980w10->Draw("SAME");
	f1_bw2->SetLineColor(4);
	f1_bw2->Draw("SAME");

//	h_f0mass_990w100->SetLineColor(6);
//	h_f0mass_990w100->Draw("SAME");
	// f1_bw3->Draw("SAME");

	h_f0mass->SetMarkerColor(1);
	h_f0mass->SetLineColor(1);
	h_f0mass->Draw("SAME");
	f1_bw->SetLineColor(1);
	f1_bw->Draw("SAME");

	TLegend *le_test=new TLegend(0.5,0.6,0.85,0.85);
	le_test->SetBorderSize(0);
	le_test->AddEntry((TObject*)0,"PYTHIA f0 mass and fit","");
	le_test->AddEntry(h_f0mass,"Default","l");
	le_test->AddEntry(h_f0mass_980w10,"Mean 980 width 10","l");
	le_test->Draw("SAME");


	SavePlotDirs(c_test,"f0massFit",{"Systematics","f0Shape"});

	// f1_exp->Draw("SAME");
	// f1_exp2->Draw("SAME");
	// f1_exp3->Draw("SAME");



	// return 1;



	/*
		 TF1 *f1_fun=new TF1("f1_fun",rightExp,binLow,binHigh,3);
		 f1_exp->SetParameters(1,1,1);

		 h_f0mass->Fit("f1_fun","QNO","",binLow_test,binHigh);
		 h_f0mass->Fit("f1_fun","QNO","",binLow_test,binHigh);
		 h_f0mass->Fit("f1_fun","LQNO","",binLow_test,binHigh);
		 h_f0mass->Fit("f1_fun","L MQNO","",binLow_test,binHigh);
		 h_f0mass->Fit("f1_fun","L M INO","",binLow_test,binHigh);
		 h_f0mass->Fit("f1_fun","L M IO","",binLow_test,binHigh);

		 f1_fun->SetRange(binLow,binHigh);

		 f1_fun->Draw("SAME");
		 */

	// bw function method

	TCanvas *c_default_bw=new TCanvas("c_default_bw","c_default_bw",800,800);
	c_default_bw->cd();

	cout<<"start bw method"<<endl;

	int nbin_bw=800;
	double binLow_bw=0.97;
	double binHigh_bw=1.37;
	// binHigh_bw=1.1;

	f1_bw->SetRange(MaxBinVal,binHigh_bw); 

	TH1D *h_default_bw=new TH1D("h_default_bw","h_default_bw",nbin_bw,binLow_bw,binHigh_bw); h_default_bw->Sumw2();
	double Max_bw=f1_bw->Eval(MaxBinVal);
	int MaxBin_bw=h_default_bw->GetXaxis()->FindBin(MaxBinVal);

	cout<<"Max_bw = "<<Max_bw<<""<<endl;

	TF1 *f1_pol1=new TF1("f1_pol1","[0]*(x-[1])");
	f1_pol1->SetParameter(0, Max_bw/(MaxBinVal-f0_MassMin) );
	f1_pol1->SetParameter(1, f0_MassMin);

	f1_pol1->SetRange(f0_MassMin,MaxBinVal);
	//	f1_pol1->Draw("SAME");

	TCanvas *c_test_bw= new TCanvas("c_test_bw","c_test_bw",800,800);
	c_test_bw->cd();
	h_f0mass->Draw("SAME");
	f1_pol1->Draw("SAME");
	f1_bw->Draw("SAME");

	double binVal=0;
	double binHeight=0;

	for(int i=0; i<nbin_bw; i++){
		binVal=h_default_bw->GetBinCenter(i+1);

		if(binVal> f0_MassMin && binVal<MaxBinVal){
			binHeight=f1_pol1->Eval(binVal);	
			//	cout<<"i = "<<i<<" , binVal = "<<binVal<<" ,  binHight = "<<binHight<<endl;
		}else if(binVal >=MaxBinVal && binVal> f0_MassMin ){
			binHeight=f1_bw->Eval(binVal);
		}else{
			binHeight=0;
		}
		cout<<"i = "<<i<<" , binVal = "<<binVal<<" ,  binHeight = "<<binHeight<<endl;
		h_default_bw->SetBinContent(i+1,binHeight);

	}

	c_default_bw->cd();
	h_default_bw->Scale(1/h_default_bw->Integral());
	h_default_bw->Draw();

	double Fr_default_bw=0;
	Fr_default_bw=h_default_bw->Integral(h_default_bw->GetXaxis()->FindBin(DtktkResmassCutMean-DtktkResmassCutWidth), h_default_bw->GetXaxis()->FindBin(DtktkResmassCutMean + DtktkResmassCutWidth));
	cout<<"Fr_default_bw = "<<Fr_default_bw<<endl;
	cout<<"Fr_ori = "<<Fr_ori<<endl; // difference due to we do not save the tail part of MC, to compare , need to include all shape


	// Loop for variation

	double f0_MassMeanLow=0.97;
	double f0_MassMeanHigh=1.01;
	double f0_WidthLow=0.010;
	double f0_WidthHigh=0.100;

	int nbin_f0_MassScan=41;
	int nbin_f0_WidthScan=91;

	double MassCurrent=0;
	double WidthCurrent=0;

	double FrCurrent=0;
	double FrMax=0;
	double FrMin=1;
	int massbin_Max=0;
	int massbin_Min=0;
	int widthbin_Max=0;
	int widthbin_Min=0;

	TH1D *h_scan[nbin_f0_MassScan][nbin_f0_WidthScan];

	for(int ibinWidth=0; ibinWidth<nbin_f0_WidthScan; ibinWidth++){
		WidthCurrent=f0_WidthLow+(double)ibinWidth*(f0_WidthHigh-f0_WidthLow)/(nbin_f0_WidthScan -1);
		cout<<"WidthCurrent = "<<WidthCurrent<<endl;

		for(int ibinMass=0; ibinMass<nbin_f0_MassScan; ibinMass++){
			h_scan[ibinMass][ibinWidth]= new TH1D(Form("h_scan_Mass%i_Width%i", ibinMass, ibinWidth),  Form("h_scan_Mass%i_Width%i", ibinMass, ibinWidth), nbin_bw, binLow_bw, binHigh_bw);

			MassCurrent=f0_MassMeanLow+(double)ibinMass*(f0_MassMeanHigh-f0_MassMeanLow)/(nbin_f0_MassScan-1);
			cout<<"MassCurrent = "<<MassCurrent<<endl;

				// set f1
				f1_bw->SetParameter(1,WidthCurrent);
				f1_bw->SetParameter(2,MassCurrent);

				f1_pol1->SetParameter(0, Max_bw/(MassCurrent-f0_MassMin) );


			for(int ibinH=0; ibinH<nbin_bw; ibinH++){

				binVal=h_scan[ibinMass][ibinWidth]->GetBinCenter(ibinH+1);
				if(binVal> f0_MassMin && binVal<MassCurrent){
					binHeight=f1_pol1->Eval(binVal);	
					//	cout<<"i = "<<i<<" , binVal = "<<binVal<<" ,  binHight = "<<binHight<<endl;
				}else if(binVal >=MassCurrent && binVal > f0_MassMin){
					binHeight=f1_bw->Eval(binVal);
				}else{
					binHeight=0;
				}
				//		cout<<"i = "<<ibinH<<" , binVal = "<<binVal<<" ,  binHeight = "<<binHeight<<endl;
				h_scan[ibinMass][ibinWidth]->SetBinContent(ibinH+1,binHeight);


			}// end for ibinH

			// normalize
			h_scan[ibinMass][ibinWidth]->Scale(1/h_scan[ibinMass][ibinWidth]->Integral());

			cout<<"Fr : Mass = "<<MassCurrent<<" , Width = "<<WidthCurrent<<" Fr = "<<h_scan[ibinMass][ibinWidth]->Integral(h_scan[ibinMass][ibinWidth]->GetXaxis()->FindBin(DtktkResmassCutMean-DtktkResmassCutWidth), h_scan[ibinMass][ibinWidth]->GetXaxis()->FindBin(DtktkResmassCutMean+DtktkResmassCutWidth) )/ h_scan[ibinMass][ibinWidth]->Integral() <<endl;
			FrCurrent=h_scan[ibinMass][ibinWidth]->Integral(h_scan[ibinMass][ibinWidth]->GetXaxis()->FindBin(DtktkResmassCutMean-DtktkResmassCutWidth), h_scan[ibinMass][ibinWidth]->GetXaxis()->FindBin(DtktkResmassCutMean+DtktkResmassCutWidth) )/ h_scan[ibinMass][ibinWidth]->Integral();

			if(FrCurrent>FrMax){
				FrMax=FrCurrent;
				massbin_Max=ibinMass;
				widthbin_Max=ibinWidth;
			}

			if(FrCurrent<FrMin){
				FrMin=FrCurrent;
				massbin_Min=ibinMass;
				widthbin_Min=ibinWidth;
			}


		} // end for int ibinMass<nbin_f0_MassScan
	} // end for int ibinWidth<nbin_f0_WidthScan

	cout<<"FrMax = "<<FrMax<<" , mass_Max = "<<f0_MassMeanLow+(double)massbin_Max*(f0_MassMeanHigh-f0_MassMeanLow)/(nbin_f0_MassScan-1)<<" , widthbin_Max = "<<f0_WidthLow+(double)widthbin_Max*(f0_WidthHigh-f0_WidthLow)/(nbin_f0_WidthScan -1)<<endl;
	cout<<"FrMax = "<<FrMin<<" , mass_Min = "<<f0_MassMeanLow+(double)massbin_Min*(f0_MassMeanHigh-f0_MassMeanLow)/(nbin_f0_MassScan-1)<<" , widthbin_Min = "<<f0_WidthLow+(double)widthbin_Min*(f0_WidthHigh-f0_WidthLow)/(nbin_f0_WidthScan -1)<<endl;
	// cout<<"FrMin = "<<FrMin<<" , massbin_Min = "<<massbin_Min<<" , widthbin_Min = "<<widthbin_Min<<endl;

	double MaxVal=h_default_bw->GetMaximum();
	if(h_scan[massbin_Max][widthbin_Max]->GetMaximum() > MaxVal ){MaxVal= h_scan[massbin_Max][widthbin_Max]->GetMaximum() ;}
	if(h_scan[massbin_Min][widthbin_Min]->GetMaximum() > MaxVal ){MaxVal= h_scan[massbin_Min][widthbin_Min]->GetMaximum() ;}

	h_default_bw->SetMaximum(MaxVal);

	TCanvas *c_FrResult=new TCanvas("c_FrResult","c_FrResult",800,800);
	c_FrResult->cd();
	gStyle->SetOptStat(0);
	h_default_bw->GetXaxis()->SetRangeUser(0.95,1.2);
	h_default_bw->SetTitle("");
	h_default_bw->GetXaxis()->SetTitle("M_{KK} GeV");
	h_default_bw->SetLineColor(1);
	h_default_bw->Draw("SAME");
	
	h_scan[massbin_Max][widthbin_Max]->SetLineColor(2);
	h_scan[massbin_Max][widthbin_Max]->Draw("SAME");
	h_scan[massbin_Min][widthbin_Min]->SetLineColor(4);
	h_scan[massbin_Min][widthbin_Min]->Draw("SAME");

	TLegend *le_FrResult=new TLegend(0.5,0.6,0.82,0.82, NULL,"brNDC");
	le_FrResult->SetBorderSize(0);
	le_FrResult->AddEntry(h_default_bw,"Default","l");
	le_FrResult->AddEntry(h_scan[massbin_Max][widthbin_Max],Form("Mass %0.3f Width %.3f ",f0_MassMeanLow+(double)massbin_Max*(f0_MassMeanHigh-f0_MassMeanLow)/(nbin_f0_MassScan-1), f0_WidthLow+(double)widthbin_Max*(f0_WidthHigh-f0_WidthLow)/(nbin_f0_WidthScan -1) ),"l");
	le_FrResult->AddEntry(h_scan[massbin_Min][widthbin_Min],Form("Mass %0.3f Width %.3f ",f0_MassMeanLow+(double)massbin_Min*(f0_MassMeanHigh-f0_MassMeanLow)/(nbin_f0_MassScan-1), f0_WidthLow+(double)widthbin_Min*(f0_WidthHigh-f0_WidthLow)/(nbin_f0_WidthScan -1) ),"l");

	le_FrResult->Draw("SAME");

	TLine *tl1= new TLine(DtktkResmassCutMean-DtktkResmassCutWidth,0,DtktkResmassCutMean-DtktkResmassCutWidth,MaxVal);
	// tl1->Draw("SAME");

	TLine *tl2= new TLine(DtktkResmassCutMean+DtktkResmassCutWidth,0,DtktkResmassCutMean+DtktkResmassCutWidth,MaxVal);
	// tl2->Draw("SAME");


	SavePlotDirs(c_FrResult,"f0ShapeScan",{"Systematics","f0Shape"});

	cout<<"FrMax ratio = "<<FrMax/Fr_default_bw<<" , FrMin ratio = "<<FrMin/Fr_default_bw<<endl;



	// gPad->BuildLegend();




	/*
		 TCanvas *c_test_bw2=new TCanvas("c_test_bw2","c_test_bw2",800,800);
		 c_test_bw2->Divide(2,2);
		 c_test_bw2->cd(1);
		 h_scan[0]->Draw();
		 c_test_bw2->cd(2);
		 h_scan[1]->Draw();
		 c_test_bw2->cd(3);
		 h_scan[3]->Draw();
		 c_test_bw2->cd(4);
		 h_scan[4]->Draw();
		 */



	return 1;

	// end bw function method

/*

	int nbin_compare=1400;
	double binLow_compare=0.97;
	double binHigh_compare=1.67; 
	double binwidth_default= (binHigh_compare- binLow_compare)/nbin_compare;

	cout<<"check here"<<endl;


	TH1D *h_default=new TH1D("h_default","h_default",nbin_compare,binLow_compare,binHigh_compare); h_default->Sumw2();

	double Max=f1_exp->Eval(MaxBinVal);
	cout<<"Max = "<<Max<<endl;

	//	TF1 *f1_pol1=new TF1("f1_pol1","[0]*(x-[1])");
	f1_pol1->SetParameter(0, Max/(MaxBinVal-f0_MassMin) );
	f1_pol1->SetParameter(1, f0_MassMin);

	f1_pol1->SetRange(f0_MassMin,MaxBinVal);
	f1_pol1->Draw("SAME");

	//	double binVal=0;
	double binHight=0;

	//DtktkResmassCutMean
	//DtktkResmassCutWidth

	for(int i=0; i<nbin_compare;i++){
		binVal=h_default->GetBinCenter(i+1);

		if(binVal> f0_MassMin && binVal<MaxBinVal){
			binHight=f1_pol1->Eval(binVal);	
			//	cout<<"i = "<<i<<" , binVal = "<<binVal<<" ,  binHight = "<<binHight<<endl;
		}else if(binVal >=MaxBinVal){
			binHight=f1_exp->Eval(binVal);
		}else{
			binHight=0;
		}
		cout<<"i = "<<i<<" , binVal = "<<binVal<<" ,  binHight = "<<binHight<<endl;
		h_default->SetBinContent(i+1,binHight);
	}


	// scale
	double MaxScale =(double)100/(double)70;
	double MinScale =(double)10/(double)70;
	double ScaleDiff=MaxScale-MinScale;

	double CurrentScale=MinScale;
	int nScale=10;
	double binwidth_current=binwidth_default;
	double binHigh_compare_current=binHigh_compare;
	double MaxBinVal_Current=MaxBinVal;
	int MaxBin_Current=MaxBin;


	TH1D *h_Scale[nScale+1];

	for(int i=0; i<nScale+1;i++) // running to nScale, so total Nscale+1 points
	{
		CurrentScale=MinScale+i*ScaleDiff/nScale;
		binwidth_current=binwidth_default*CurrentScale;
		binHigh_compare_current=binLow_compare+(binHigh_compare-binLow_compare)*CurrentScale;
		cout<<"CurrentScale = "<<CurrentScale<<" , binHigh_compare_current "<<binHigh_compare_current<<endl;

		h_Scale[i]=new TH1D(Form("h_Scale%i",i),Form("h_Scale%i",i),nbin_compare,binLow_compare,binHigh_compare_current); 
		MaxBin_Current=h_Scale[i]->GetXaxis()->FindBin(MaxBinVal_Current); // MaxBinVal_Current : mass peak position
		cout<<"MaxBin_Current = "<<MaxBin_Current<<endl;
		f1_pol1->SetParameter(0, Max/(MaxBinVal_Current-f0_MassMin) );

		for(int j=0; j<nbin_compare;j++){
			binVal=h_Scale[i]->GetBinCenter(i+1);
			if(binVal >f0_MassMin && (j+1) < MaxBin_Current){
				binHight=f1_pol1->Eval(binVal);
			}else if( (j+1)>=MaxBin_Current &&  MaxBin+ (j+1-MaxBin_Current ) < nbin_compare ){
				binHight=h_default->GetBinContent(MaxBin+ (j+1-MaxBin_Current) );
			}else{
				binHight=0;
			}
			h_Scale[i]->SetBinContent(j+1, binHight);
			//			cout<<"j = "<<j<<" , binHight= "<<binHight<<endl;
		}// end for j<nbin_compare 
		h_Scale[i]->Scale(1/h_Scale[i]->Integral());
		cout<<"fr = "<<h_Scale[i]->Integral(h_default->GetXaxis()->FindBin(DtktkResmassCutMean-DtktkResmassCutWidth), h_default->GetXaxis()->FindBin(DtktkResmassCutMean + DtktkResmassCutWidth) )<<endl;

	}// end for i<nScale+1


	TCanvas *c_Scale=new TCanvas("c_Scale","c_Scale",800,800);
	c_Scale->Divide(2,2);
	c_Scale->cd(1);
	h_Scale[0]->Draw();
	c_Scale->cd(2);
	h_Scale[1]->Draw();
	c_Scale->cd(3);
	h_Scale[8]->Draw();
	c_Scale->cd(4);
	h_Scale[9]->Draw();



	h_default->Scale(1/h_default->Integral());

	double default_fr=h_default->Integral(h_default->GetXaxis()->FindBin(DtktkResmassCutMean-DtktkResmassCutWidth), h_default->GetXaxis()->FindBin(DtktkResmassCutMean + DtktkResmassCutWidth) );

	TCanvas *c_default = new TCanvas("c_default","c_default",800,800);
	c_default->cd();
	h_default->Draw();






	cout<<"default_fr = "<<default_fr<<endl;
	cout<<"MaxBin = "<<MaxBin<<" , MaxBinVal = "<<MaxBinVal<<endl;



	return 0;
*/

}
