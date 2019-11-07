
	// 1. binned fit
	// 2. unbinned fit
  // 3. background function variation
	// 4. signal width variation
	// 5. variable cut scan
	// 6. saving result, both number and fuction for plot
  // 7. plot marco (in another code)

	// must use new version of root, root6.02 crash

// add MC study & pull distribution later

// prompt/nonprompt weight: https://root-forum.cern.ch/t/roofit-problem-in-replacing-the-deprecated-setweightvar/11384
// using ds1->append(*d2) to combine ROODS , https://root.cern.ch/root/html/tutorials/roofit/rf402_datahandling.C.html


#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"

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
#include <RooFitResult.h>
#include <RooWorkspace.h>
#include <RooMCStudy.h>
#include "RooHist.h"
#include "RooConstVar.h"



#include "TFitter.h"
#include "TFitResult.h"

using namespace RooFit;
using namespace std;


// TCanvas *ctest1=new TCanvas("ctest1","ctest1",800,800);
// TCanvas *c_roofit = new TCanvas("c_roofit","c_roofit",800,800);

	double DsDataFitRangeLow =1.91;
	double DsDataFitRangeHigh = 2.11;

	double shiftY=0.05;
	double shiftX=0.32;
	double oneshift=0.070;

	TCanvas *c_binfit_MC[500];
	TCanvas *c_binfit_Data[500];
	TCanvas *c_roofit_MC[500] ;
	TCanvas *c_roofit_Data[500];
	TCanvas *c_roofit_Data_pull[500];
	TCanvas *c_roofit_Data_withpull[500];

	TCanvas *c_roofit_MCstudy[500];

	TCanvas *c_roofit_Data_Dca[500][30];

	int count_c_binfit=0;
	int count_c_rootfit=0;

	double FloatWidthErr_Norm=0;
	double FloatWidthVal_Norm=0;
	double FloatWidth_DsMass_Norm=1.965;

Double_t fitf_bkg1(Double_t *x, Double_t *par)
{
   if (TMath::Abs(x[0]) > 1.93 && TMath::Abs(x[0]) < 2.01) {
      TF1::RejectPoint();
      return 0;
   }

   // Double_t arg = 0;
   // if (par[2]) arg = (x[0] - par[1])/par[2];
//   Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
		Double_t fitval=par[0]+par[1]*x[0];
   return fitval;
}

Double_t fitf_bkg2(Double_t *x, Double_t *par)
{
   if (TMath::Abs(x[0]) > 1.93 && TMath::Abs(x[0]) < 2.01) {
      TF1::RejectPoint();
      return 0;
   }

   // Double_t arg = 0;
   // if (par[2]) arg = (x[0] - par[1])/par[2];
//   Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
		Double_t fitval=par[0]+par[1]*x[0]+par[2]*(x[0]*x[0]*2-1);
   return fitval;
}


int TF1FitBkg(TH1F* h_DsMassData,  Int_t ibin_Dpt, Double_t DptLow=0 ,Double_t DptHigh=10, TString str_PbPb="pp"){

	//gStyle->SetOptFit(1111);

	TF1 *bkg1=new TF1("fitf_bkg1",fitf_bkg1,DsDataFitRangeLow,DsDataFitRangeHigh,2);
	bkg1->SetParameters(h_DsMassData->Integral(),0);
	bkg1->SetLineColor(2);
	bkg1->SetParNames("Constant","1st order");
	h_DsMassData->Fit("fitf_bkg1","ME I0 S");
	h_DsMassData->Fit("fitf_bkg1","ME I0 S");
	h_DsMassData->Fit("fitf_bkg1","ME I0 S");
	h_DsMassData->Fit("fitf_bkg1","ME I0 S");
	h_DsMassData->Fit("fitf_bkg1","ME I0 S");
	h_DsMassData->Fit("fitf_bkg1","ME I0 S");

	TF1 *bkg2=new TF1("fitf_bkg2",fitf_bkg2,DsDataFitRangeLow,DsDataFitRangeHigh,3);
	bkg2->SetParameters(h_DsMassData->Integral(),0);
	bkg2->SetParNames("Constant","1st order","2nd order");
	bkg2->SetLineColor(4);
	h_DsMassData->Fit("fitf_bkg2","ME I0 S");
	h_DsMassData->Fit("fitf_bkg2","ME I0 S");
	h_DsMassData->Fit("fitf_bkg2","ME I0 S");
	h_DsMassData->Fit("fitf_bkg2","ME I0 S");
	h_DsMassData->Fit("fitf_bkg2","ME I0 S");
	h_DsMassData->Fit("fitf_bkg2","ME I0 S");





	TCanvas *c_TF1FitBkg=new TCanvas("c_TF1FitBkg","c_TF1FitBkg",800,800);
	c_TF1FitBkg->cd();
	h_DsMassData->Draw();
	h_DsMassData->SetTitle("");
	h_DsMassData->GetXaxis()->SetTitle("m_{KK#pi}");
	h_DsMassData->GetXaxis()->CenterTitle();
	bkg1->SetRange(DsDataFitRangeLow,DsDataFitRangeHigh);
	bkg1->Draw("SAME");
	bkg2->SetRange(DsDataFitRangeLow,DsDataFitRangeHigh);
	bkg2->Draw("SAME");

	gStyle->SetOptStat(0);


	TLegend *le_TF1FitBkg=new TLegend(0.55,0.65,0.85,0.85);
	le_TF1FitBkg->SetBorderSize(0);
	le_TF1FitBkg->AddEntry((TObject*)0,Form("%s , %.0f < D_{S} p_{T} < %.0f",str_PbPb.Data(),DptLow,DptHigh),"");
	le_TF1FitBkg->AddEntry(bkg1,Form("1st, %.3e/%i", bkg1->GetChisquare(),bkg1->GetNDF()),"l");
  le_TF1FitBkg->AddEntry(bkg2,Form("2nd, %.3e/%i", bkg2->GetChisquare(),bkg2->GetNDF()),"l");
	le_TF1FitBkg->Draw("same");

	
	SavePlotDirs(c_TF1FitBkg,Form("%s_Data_bkgtestTF1_pt%.0fto%.0f",str_PbPb.Data(),DptLow,DptHigh),{"SignalFit","bkgTest",str_PbPb.Data()});


	return 0;

}


int bkgFitTest(TH1F* h_DsMassData,  Int_t ibin_Dpt, Double_t DptLow=0 ,Double_t DptHigh=10, TString str_PbPb="pp"){
// using roofit with binned fit test chi2

	RooRealVar x("x","x",DsDataFitRangeLow,DsDataFitRangeHigh) ;
	x.setRange("sb_lo",DsDataFitRangeLow,1.93);
	x.setRange("sb_hi",2.01,DsDataFitRangeHigh);


  RooRealVar Cheb1("Cheb1","Cheb1",0,-1,1); // no input is better than with input
  RooRealVar Cheb2("Cheb2","Cheb2",0,-1,1);
  RooRealVar Cheb3("Cheb3","Cheb3",0,-1,1);

	RooDataHist dh("dh","dh",x,Import(*h_DsMassData));

	RooChebychev *Chebyshev1=new RooChebychev("Chebyshev1","Chebyshev1",x,RooArgList(Cheb1));

	Chebyshev1->fitTo(dh,Range("sb_lo,sb_hi"));

  RooPlot *frame_1stBkg=x.frame(Title("1st order fit"));
  dh.plotOn(frame_1stBkg,DataError(RooAbsData::SumW2));
	Chebyshev1->plotOn(frame_1stBkg);


	RooChebychev *Chebyshev2=new RooChebychev("Chebyshev2","Chebyshev2",x,RooArgList(Cheb1,Cheb2));

	Chebyshev2->fitTo(dh,Range("sb_lo,sb_hi"));

  RooPlot *frame_2ndBkg=x.frame(Title("2nd order fit"));
  dh.plotOn(frame_2ndBkg,DataError(RooAbsData::SumW2));
	Chebyshev2->plotOn(frame_2ndBkg);


	TCanvas *c_bkgFitTest=new TCanvas("c_bkgFitTest","c_bkgFitTest",800,800);
	c_bkgFitTest->Divide(2,2);
	c_bkgFitTest->cd(1);
  frame_1stBkg->Draw();

	c_bkgFitTest->cd(2);
	frame_2ndBkg->Draw();

	SavePlotDirs(c_bkgFitTest,Form("%s_Data_bkgtest_pt%.0fto%.0f",str_PbPb.Data(),DptLow,DptHigh),{"SignalFit","bkgTest",str_PbPb.Data()});


	cout<<"str_PbPb = "<<str_PbPb.Data()<<endl;

	return 0;

}

TF1 *binnedFit(TH1F* h_DsMassDsMC, TH1F* h_DsMassData, TH1F* h_RawBinFitYield, Int_t ibin_Dpt, Double_t DptLow=0 ,Double_t DptHigh=10 ,TString str_PbPb="pp", Int_t doBkgVariation=0 ,Int_t BkgFunction=0, Int_t doSigVariation=0, double SigWidthErr=0 ,TString extraName=""){


  double tex_upperY=0.95;


  TLatex* texCmsPre = new TLatex(0.15,tex_upperY, "#scale[1.25]{CMS} Projection");
  texCmsPre->SetNDC();
  texCmsPre->SetTextAlign(12);
  texCmsPre->SetTextSize(0.035);
  texCmsPre->SetTextFont(42);

  TLatex* texCmsSim = new TLatex(0.15,tex_upperY, "#scale[1.25]{CMS} Simulations");
  texCmsSim->SetNDC();
  texCmsSim->SetTextAlign(12);
  texCmsSim->SetTextSize(0.035);
  texCmsSim->SetTextFont(42);


  // TLatex* texColPbPb = new TLatex(0.88,tex_upperY, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  TLatex* texColPbPb = new TLatex(0.95,tex_upperY, "0.2 nb^{-1} (5.02 TeV PbPb)");
  texColPbPb->SetNDC();
  texColPbPb->SetTextAlign(32);
  texColPbPb->SetTextSize(0.035);
  texColPbPb->SetTextFont(42);

  // TLatex* texColpp = new TLatex(0.88,tex_upperY, "pp #sqrt{s} = 5.02 TeV");
  TLatex* texColpp = new TLatex(0.95,tex_upperY, "650 pb^{-1} (5.02 TeV pp )");
  texColpp->SetNDC();
  texColpp->SetTextAlign(32);
  texColpp->SetTextSize(0.035);
  texColpp->SetTextFont(42);


	cout<<"check a1"<<endl;
	// double DsMass=1.96828; // from parameters.h
	double DsMassCandWidth=0.05;
	double DsMassCandFitMeanWidth=0.02;

	// double DsDataFitRangeLow =1.91;
	// double DsDataFitRangeHigh = 2.11;


	cout<<"check a2"<<endl;
	// MC fit
	// TF1 *f1_DsMC = new TF1("f1_DsMC","[0]*( [4]*TMath::Gaus(x,[1],[2])/( sqrt(2*TMath::Pi()) *[2])  + (1-[4])*TMath::Gaus(x,[1],[3])/(sqrt(2*TMath::Pi()) *[3]) )"); // kTRUE nomalize Gaus
	TF1 *f1_DsMC = new TF1("f1_DsMC","[0]*( [4]*TMath::Gaus(x,[1],[2]*(1+[5]))/( sqrt(2*TMath::Pi())*[2]*(1+[5]) ) + (1-[4])*TMath::Gaus(x,[1],[3]*(1+[5]))/(sqrt(2*TMath::Pi())*[3]*(1+[5])  ) )"); // kTRUE nomalize Gaus
	f1_DsMC->SetParameter(0, h_DsMassDsMC->GetBinWidth(1) * h_DsMassDsMC->Integral( h_DsMassDsMC->FindBin(DsMass-DsMassCandWidth) , h_DsMassDsMC->FindBin(DsMass-DsMassCandWidth) ) ) ;
	f1_DsMC->SetParameter(1,DsMass);   // DsMean
	f1_DsMC->SetParLimits(1,DsMass-DsMassCandFitMeanWidth,DsMass+DsMassCandFitMeanWidth);
	f1_DsMC->SetParameter(2,0.01); // defaut (wider) gaus1 sigma
	f1_DsMC->SetParLimits(2,0.0001,0.1);
	f1_DsMC->SetParameter(3,0.001); // default (narrow) gaus2 sigma 
	f1_DsMC->SetParLimits(3,0.00001,0.1);
	f1_DsMC->SetParameter(4,0.5);  // fraction of Gaus1 / All
	f1_DsMC->SetParLimits(4,0,1);
	f1_DsMC->SetParameter(5,0); // for data /mc discrepancy
	f1_DsMC->SetParLimits(5,2,2);
	f1_DsMC->FixParameter(5,0);	

	f1_DsMC->SetLineColor(kRed);
	// fitting option : Q : quite, N: do not save fit function to histogram , 0 : not draw, M : improve fit result by TMinuit, L : log likelihood method, I : use integral instead of bin center, E : better erros estimation by minos , S : fit result returned in TFitResultPtr
	// see https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html#the-th1fit-method

	cout<<"check a3"<<endl;

	h_DsMassDsMC->Fit("f1_DsMC","QN0","",    DsMass-DsMassCandWidth,DsMass+DsMassCandWidth);
	h_DsMassDsMC->Fit("f1_DsMC","QN0","",    f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);
	h_DsMassDsMC->Fit("f1_DsMC","M IN0","",    f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);
	h_DsMassDsMC->Fit("f1_DsMC","L QN0","",  f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);
	h_DsMassDsMC->Fit("f1_DsMC","LN0","",    f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);
	h_DsMassDsMC->Fit("f1_DsMC","L MN0","",  f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);
	h_DsMassDsMC->Fit("f1_DsMC","L MN0","",  f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);
	h_DsMassDsMC->Fit("f1_DsMC","L M IN0","",f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);
	h_DsMassDsMC->Fit("f1_DsMC","L M IN0","",f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);
	h_DsMassDsMC->Fit("f1_DsMC","L M IN0","",f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);
	h_DsMassDsMC->Fit("f1_DsMC","L ME I0 S","",f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth);

	//c_fit[ibin_Dpt]=new TCanvas(Form("c_fit_%i",ibin_Dpt),Form("c_fit_%i",ibin_Dpt),800,800);
	//c_fit[ibin_Dpt]->cd();
	//h_DsMassDsMC->Draw();


	cout<<"check a4"<<endl;

	double yieldDsMC= f1_DsMC->Integral(f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth)/ h_DsMassDsMC->GetBinWidth(1);
	double yieldDsMCErr = f1_DsMC->Integral(DsMass-DsMassCandWidth,DsMass+DsMassCandWidth)/h_DsMassDsMC->GetBinWidth(1) * f1_DsMC->GetParError(0)/f1_DsMC->GetParameter(0);

	// cout<<"parameter 0 = "<<f1_DsMC->GetParameter(0)<<" BinWidth = "<< h_DsMassDsMC->GetBinWidth(1) <<endl;
	// cout<<"yield DsMC = "<< yieldDsMC <<" +- "<<yieldDsMCErr<<endl;

	c_binfit_MC[count_c_binfit]=new TCanvas(Form("c_binfit_MC_%i",count_c_binfit),"c_binfit_MC",800,800);
	c_binfit_MC[count_c_binfit]->cd();
	h_DsMassDsMC->Draw();
//	f1_DsMC->Draw("SAME"); // not work
	TF1 *f1_MC=(TF1*)h_DsMassDsMC->GetFunction("f1_DsMC");
	f1_MC->Draw("SAME");

	shiftY=0;
	TLatex *tl_binfitMC =new TLatex();
	tl_binfitMC->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s_MC_binfit",str_PbPb.Data()));  shiftY-=oneshift;
	tl_binfitMC->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%.0f < pt < %.0f ",DptLow,DptHigh)); shiftY-=oneshift;
	if(extraName!="")
	{ 	
	tl_binfitMC->DrawLatexNDC(textposx,textposy+shiftY,extraName); shiftY-=oneshift;
	}
	tl_binfitMC->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Mean=%.4f, frGaus1=%.3f  ",f1_MC->GetParameter(1),f1_MC->GetParameter(4) ));  shiftY-=oneshift;
	tl_binfitMC->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Sigma1=%.4f, Sigma2=%.4f",f1_MC->GetParameter(2), f1_MC->GetParameter(3)));  shiftY-=oneshift;
	// SavePlotDirs(c_binfit_MC[count_c_binfit],Form("%s_MC_binfit_pt%.0fto%.0f_%s",str_PbPb.Data(),DptLow,DptHigh,extraName.Data() ),{"SignalFit",str_PbPb,"binfit"} );


	// TF1 *f1_DsMix = new TF1("f1_DsMix"," [0]*( [4]*TMath::Gaus(x,[1],[2]*(1+[5]))/( sqrt(2*TMath::Pi())*[2]*(1+[5]) ) + (1-[4])*TMath::Gaus(x,[1],[3]*(1+[5]))/(sqrt(2*TMath::Pi())*[3]*(1+[5])  ) )  +  ([6]+[7]*x+[8]*x*x+[9]*x*x*x)  ");    // overall fit function for Data fitting, Ds , Dp ,background

	TF1 *f1_DsSignal = (TF1*)f1_DsMC->Clone("f1_DsSignal"); 
	// TF1 *f1_DsBkg= new TF1("f1_DsBkg","TMath::Chebyshev3(x,[0],[1],[2],[3])"); not work
	// TF1 *f1_DsBkg= new TF1("f1_DsBkg","[0]*(1+[1]*x+[2]*(2*x*x-1)+[3]*(4*x*x*x-3*x))");
	TF1 *f1_DsBkg= new TF1("f1_DsBkg","[0]*(1+[1]*x+[2]*(2*x*x-1))");
	// TF1 *f1_DsBkg= new TF1("f1_DsBkg","[0]*(1+[1]*x)");
	TF1 *f1_DsBkg_clone=(TF1*)f1_DsBkg->Clone("f1_DsBkg_clone");  // stupid bug in root, must use clone function
	TF1 *f1_DsMix = new TF1("f1_DsMix","f1_DsSignal+f1_DsBkg");	


	// Gaus2 for Ds 
/*
	f1_DsMix->SetParameter(0,1000);
	f1_DsMix->SetParLimits(0,0,100000000);		
	f1_DsMix->SetParameter(1,f1_DsMC->GetParameter(1)); // mass mean
	f1_DsMix->SetParLimits(1,DsMass-DsMassCandFitMeanWidth,DsMass+DsMassCandFitMeanWidth);
	f1_DsMix->SetParameter(2,f1_DsMC->GetParameter(2)); 	// gauss1 sigma
	f1_DsMix->SetParLimits(2,0.00001,0.1);	
	f1_DsMix->SetParameter(3,f1_DsMC->GetParameter(3));	// gauss2 sigma
	f1_DsMix->SetParLimits(3,0.00001,0.1);	
	f1_DsMix->SetParameter(4,f1_DsMC->GetParameter(4));	// fraction gaus1/ all
	f1_DsMix->SetParLimits(4,0,1);	
*/

	f1_DsMix->SetParameter(1,f1_DsMC->GetParameter(1)); // mass mean
	f1_DsMix->SetParLimits(1,DsMass-DsMassCandFitMeanWidth,DsMass+DsMassCandFitMeanWidth);
	f1_DsMix->SetParameter(5,0); // width ratio difference for Data/ MC, data width = (1+5)*MC width 	
	f1_DsMix->SetParLimits(5,-1,1);

	// fixing parameter for initial fit to get background first
	f1_DsMix->FixParameter(1,f1_DsMC->GetParameter(1));  // Ds Mass mean
	f1_DsMix->FixParameter(2,f1_DsMC->GetParameter(2));  // Ds Sigma1
	f1_DsMix->FixParameter(3,f1_DsMC->GetParameter(3));  // Ds Sigma2
	f1_DsMix->FixParameter(4,f1_DsMC->GetParameter(4));  // Ds Gaus1/all yield ratio
	// f1_DsMix->FixParameter(5,0);  // Ds Sigma1/all yield ratio

	f1_DsMix->SetLineColor(kRed);

	// cheby3 for bkg
	f1_DsMix->SetParameter(6,1000);
	f1_DsMix->SetParameter(7,0);
	// f1_DsMix->SetParameter(8,0);
	// f1_DsMix->SetParameter(9,0);

	f1_DsMix->SetParameter(6,h_DsMassData->Integral()*h_DsMassData->GetBinWidth(1));
	//f1_DsMix->SetParLimits(6,0,1000000000000);

	h_DsMassData->Fit("f1_DsMix","QN0","",    DsDataFitRangeLow, DsDataFitRangeHigh );
	h_DsMassData->Fit("f1_DsMix","QN0","",    DsDataFitRangeLow, DsDataFitRangeHigh);
	h_DsMassData->Fit("f1_DsMix","L QN0","",  DsDataFitRangeLow, DsDataFitRangeHigh);
	h_DsMassData->Fit("f1_DsMix","LN0","",    DsDataFitRangeLow, DsDataFitRangeHigh);
	h_DsMassData->Fit("f1_DsMix","L MN0","",  DsDataFitRangeLow, DsDataFitRangeHigh);
	h_DsMassData->Fit("f1_DsMix","L MN0","",  DsDataFitRangeLow, DsDataFitRangeHigh);
	h_DsMassData->Fit("f1_DsMix","L M IN0","",DsDataFitRangeLow, DsDataFitRangeHigh);
	h_DsMassData->Fit("f1_DsMix","L M IN0","",DsDataFitRangeLow, DsDataFitRangeHigh);

	f1_DsMix->ReleaseParameter(5);
	f1_DsMix->ReleaseParameter(1);

	h_DsMassData->Fit("f1_DsMix","QN0","",    DsDataFitRangeLow, DsDataFitRangeHigh );
	h_DsMassData->Fit("f1_DsMix","QN0","",    DsDataFitRangeLow, DsDataFitRangeHigh);
	h_DsMassData->Fit("f1_DsMix","L QN0","",  DsDataFitRangeLow, DsDataFitRangeHigh);
	h_DsMassData->Fit("f1_DsMix","LN0","",    DsDataFitRangeLow, DsDataFitRangeHigh);
	h_DsMassData->Fit("f1_DsMix","L MN0","",  DsDataFitRangeLow, DsDataFitRangeHigh);
	h_DsMassData->Fit("f1_DsMix","L MN0","",  DsDataFitRangeLow, DsDataFitRangeHigh);
	h_DsMassData->Fit("f1_DsMix","L M IN0","",DsDataFitRangeLow, DsDataFitRangeHigh);
	h_DsMassData->Fit("f1_DsMix","L M IN0","",DsDataFitRangeLow, DsDataFitRangeHigh);
	h_DsMassData->Fit("f1_DsMix","L EMI S0","",DsDataFitRangeLow, DsDataFitRangeHigh);

// fit status code : https://sft.its.cern.ch/jira/browse/ROOT-4445
// 4: MIGRAD error
// 40: MINOS error
// 400: HESSE error
// 4000: IMPROVE error


	int fitStatus=1;
	int fitIsValid=0;
	TFitResultPtr fitResut;
	double fitPrecision=1.e-8;
	while(fitStatus>0 && fitStatus!=4000){
		TFitter::SetPrecision(fitPrecision);
		fitResut=h_DsMassData->Fit("f1_DsMix","L EMI S0","",DsDataFitRangeLow, DsDataFitRangeHigh);
		fitStatus=fitResut->Status();
		fitIsValid=fitResut->IsValid();
		cout<<"fit precision: "<<TFitter::GetPrecision()<<"   status: "<<fitStatus<<" isValid = "<< fitIsValid <<endl;
		if(fitStatus){
			fitPrecision *= 3;
		}
		if(fitPrecision>1.e-3) {break;}
	}

	cout<<"f1_DsMix->GetParameter(0) = "<<f1_DsMix->GetParameter(0)<<endl;
	cout<<"f1_DsSignal->GetParameter(0) = "<<f1_DsSignal->GetParameter(0)<<endl;

//	TF1 *f1_DsSignal = new TF1("f1_DsSignal","[0]*( [4]*TMath::Gaus(x,[1],[2]*(1+[5]))/( sqrt(2*TMath::Pi()) *[2]*(1+[5])  )   +(1-[4])*TMath::Gaus(x,[1],[3]*(1+[5]))/(sqrt(2*TMath::Pi())*[3]*(1+[5])  ) )"); // kTRUE nomalize Gaus
	f1_DsSignal->SetParameters(f1_DsMix->GetParameter(0), f1_DsMix->GetParameter(1) , f1_DsMix->GetParameter(2) , f1_DsMix->GetParameter(3) , f1_DsMix->GetParameter(4), f1_DsMix->GetParameter(5));
	f1_DsSignal->SetParError(0,f1_DsMix->GetParError(0));
	f1_DsSignal->SetParError(1,f1_DsMix->GetParError(1));
	f1_DsSignal->SetParError(2,f1_DsMix->GetParError(2));
	f1_DsSignal->SetParError(3,f1_DsMix->GetParError(3));
	f1_DsSignal->SetParError(4,f1_DsMix->GetParError(4));
	f1_DsSignal->SetParError(5,f1_DsMix->GetParError(5));

	f1_DsSignal->SetFillColor(kOrange-3);
	f1_DsSignal->SetFillStyle(3002);
	f1_DsSignal->SetLineColor(kOrange-3);
	f1_DsSignal->SetLineWidth(2);
	f1_DsSignal->SetLineStyle(2);
  f1_DsSignal->SetRange(DsDataFitRangeLow,DsDataFitRangeHigh); // must setrange before plot or get empty


	// double yieldDsSignal= f1_DsSignal->Integral(f1_DsSignal->GetParameter(1)-DsMassCandWidth,f1_DsSignal->GetParameter(1)+DsMassCandWidth)/ h_DsMassData->GetBinWidth(1);
	double yieldDsSignal= f1_DsSignal->Integral(DsDataFitRangeLow,DsDataFitRangeHigh)/ h_DsMassData->GetBinWidth(1);
	double yieldDsSignalErr = f1_DsSignal->Integral(DsMass-DsMassCandWidth,DsMass+DsMassCandWidth)/h_DsMassData->GetBinWidth(1) * f1_DsSignal->GetParError(0)/f1_DsSignal->GetParameter(0);

	double yieldDsSignalTry = f1_DsSignal->GetParameter(0)/h_DsMassData->GetBinWidth(1);
	double yieldDsSignalTryErr = f1_DsSignal->GetParError(0)/h_DsMassData->GetBinWidth(1);

	cout<<"parameter 0 = "<<f1_DsSignal->GetParameter(0)<<" BinWidth = "<< h_DsMassData->GetBinWidth(1) <<endl;
	cout<<"yield DsSignal = "<< yieldDsSignal <<" +- "<<yieldDsSignalErr<<endl;
	cout<<"yield DsSignal from Para 0 = "<< yieldDsSignalTry <<" +- "<<yieldDsSignalTryErr<<endl;

//	TF1 *f1_DsBkg = new TF1("f1_DsBkg"," [0]+[1]*x+[2]*x*x+[3]*x*x*x ");    // overall fit function for Data fitting, Ds , Dp ,background
//	f1_DsBkg->SetParameters(f1_DsMix->GetParameter(6), f1_DsMix->GetParameter(7) ,f1_DsMix->GetParameter(8),f1_DsMix->GetParameter(9));
	f1_DsBkg->SetParameters(f1_DsMix->GetParameter(6), f1_DsMix->GetParameter(7) ,f1_DsMix->GetParameter(8));
	// f1_DsBkg->SetParameters(f1_DsMix->GetParameter(6), f1_DsMix->GetParameter(7));
	f1_DsBkg->SetParError(0,f1_DsMix->GetParError(6));
	f1_DsBkg->SetParError(1,f1_DsMix->GetParError(7));
	f1_DsBkg->SetParError(2,f1_DsMix->GetParError(8));
//	f1_DsBkg->SetParError(3,f1_DsMix->GetParError(9));

	f1_DsBkg->SetLineColor(4);
	f1_DsBkg->SetLineStyle(2);
	f1_DsBkg->SetRange(DsDataFitRangeLow,DsDataFitRangeHigh);

	// double yieldBkg= f1_DsBkg->Integral(DsDataFitRangeLow,DsDataFitRangeHigh)/ h_DsMassData->GetBinWidth(1);
	double yieldBkg= f1_DsBkg->Integral(DsDataFitRangeLow,DsDataFitRangeHigh)/ h_DsMassData->GetBinWidth(1);

	cout<<"yield Bkg = " <<yieldBkg<<endl;

	// c_fit->cd(2);

	// h_DsMassData->Draw("e");
	// f1_DsMix->Draw("same");
	// f1_DsBkg->Draw("same");
	// f1_DsSignal->Draw("same");

	// c_fitData[ibin_Dpt]->SaveAs(Form("plots/DsFitData_pt%.0fto%.0f_pp.pdf",DptLow,DptHigh));

	cout<<"ibin_Dpt = " <<ibin_Dpt<<" , yieldDsSignal = "<<yieldDsSignal<<endl;
	cout<<"fit status = "<<fitStatus<<endl;

	// from fit generate psedo data
	

	h_RawBinFitYield->SetBinContent(1, yieldDsSignal);
	h_RawBinFitYield->SetBinError(1, yieldDsSignalErr);

	if(fitStatus ==4){
		h_RawBinFitYield->SetBinContent(1, 0);
		h_RawBinFitYield->SetBinError(1, 0);
	}

  c_binfit_Data[count_c_binfit]=new TCanvas(Form("c_binfit_Data_%i",count_c_binfit),"c_binfit_Data",800,600);
  c_binfit_Data[count_c_binfit]->cd();

  c_binfit_Data[count_c_rootfit]->cd();
  c_binfit_Data[count_c_rootfit]->SetFillColor(0);
  c_binfit_Data[count_c_rootfit]->SetBorderMode(0);
  c_binfit_Data[count_c_rootfit]->SetFrameFillStyle(0);
  c_binfit_Data[count_c_rootfit]->SetFrameBorderMode(0);
  c_binfit_Data[count_c_rootfit]->SetLeftMargin( c_Lmg );
  c_binfit_Data[count_c_rootfit]->SetRightMargin( c_Rmg );
  c_binfit_Data[count_c_rootfit]->SetTopMargin( c_Tmg );
  c_binfit_Data[count_c_rootfit]->SetBottomMargin( c_Bmg );
  c_binfit_Data[count_c_rootfit]->SetTickx(0);
  c_binfit_Data[count_c_rootfit]->SetTicky(0);

  h_DsMassData->SetTitle("");
  h_DsMassData->GetXaxis()->SetTitle("m_{KK#pi} (GeV/c^{2})");
  h_DsMassData->GetXaxis()->CenterTitle();
	h_DsMassData->GetYaxis()->SetTitle(Form("Events/(%.3f)",h_DsMassData->GetBinWidth(2)));
  h_DsMassData->GetYaxis()->CenterTitle();
  h_DsMassData->Draw();


	h_DsMassData->SetMaximum(h_DsMassData->GetMaximum()*1.4);
	h_DsMassData->SetMinimum(h_DsMassData->GetMinimum()*0.8);
  h_DsMassData->Draw();
	TF1 *myfun=h_DsMassData->GetFunction("f1_DsMix");
  myfun->Draw("same");
	f1_DsBkg->SetRange(DsDataFitRangeLow,DsDataFitRangeHigh);
	f1_DsBkg->Draw("same"); //must setrange before draw or get empty
//	f1_DsSignal->Draw("same");

	shiftY=0.10;
	shiftX=0.44;
	TLatex *tl_binfitData =new TLatex();
	tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("D_{S}^{+} + D_{S}^{-}"));  shiftY-=oneshift;
	tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%.0f < pt < %.0f GeV/c ",DptLow,DptHigh)); shiftY-=oneshift;
	tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("|y| < 1.0")); shiftY-=oneshift;
	if(extraName!="")
	{ 	
	tl_binfitData->DrawLatexNDC(textposx,textposy+shiftY,extraName); shiftY-=oneshift;
	}
	tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Raw Yield = %.0f #pm %.0f ",yieldDsSignal,yieldDsSignalErr ));  shiftY-=oneshift;
	// tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Float width = %.3f",f1_DsMix->GetParameter(5))); shiftY-=oneshift;

	double le_x1=0.64;
	double le_x2=0.86;
	double le_y1=0.42;
	double le_y2=0.63;

	TLegend *le_binfit=new TLegend(le_x1,le_y1,le_x2,le_y2);
	le_binfit->SetBorderSize(0);
	le_binfit->AddEntry(h_DsMassData,"Data","lp");
	le_binfit->AddEntry(myfun,"Signal+Background","l");
	le_binfit->AddEntry(f1_DsBkg,"Background","l");
	le_binfit->Draw("SAME");


	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);

	
	writeExtraText=true;
	extraText="Projection";
	CMS_lumi(c_binfit_Data[count_c_binfit], 12, 10);

/*
  texCmsPre->Draw("SAME");
  if(str_PbPb=="pp"){
  texColpp->Draw("SAME");
  }else{
  texColPbPb->Draw("SAME");
  }
*/

	SavePlotDirs(c_binfit_Data[count_c_binfit],Form("%s_Data_fit_pt%.0fto%.0f_%s",str_PbPb.Data(),DptLow,DptHigh,extraName.Data() ),{"HLProjection","SignalFit",str_PbPb,"binfit"} );


/*
	TH1D *h1S=new TH1D("h1S","h1S",50,DsDataFitRangeLow,DsDataFitRangeHigh); h1S->Sumw2();
	TH1D *h1B=new TH1D("h1B","h1B",50,DsDataFitRangeLow,DsDataFitRangeHigh); h1B->Sumw2();
		
	TH1D *h1SB=new TH1D("h1SB","h1SB",50,DsDataFitRangeLow,DsDataFitRangeHigh); h1SB->Sumw2();

	h1S->FillRandom("f1_DsSignal",yieldDsSignal*.09);
	h1B->FillRandom("f1_DsBkg",yieldBkg*0.8);

	h1SB->Add(h1S);
	h1SB->Add(h1B);

	h1SB->Fit("f1_DsMix","QN0","",    DsDataFitRangeLow, DsDataFitRangeHigh );
	h1SB->Fit("f1_DsMix","QN0","",    DsDataFitRangeLow, DsDataFitRangeHigh);
	h1SB->Fit("f1_DsMix","L QN0","",  DsDataFitRangeLow, DsDataFitRangeHigh);
	h1SB->Fit("f1_DsMix","LN0","",    DsDataFitRangeLow, DsDataFitRangeHigh);
	h1SB->Fit("f1_DsMix","L MN0","",  DsDataFitRangeLow, DsDataFitRangeHigh);
	h1SB->Fit("f1_DsMix","L MN0","",  DsDataFitRangeLow, DsDataFitRangeHigh);
	h1SB->Fit("f1_DsMix","L M IN0","",DsDataFitRangeLow, DsDataFitRangeHigh);
	h1SB->Fit("f1_DsMix","L M IN0","",DsDataFitRangeLow, DsDataFitRangeHigh);
	h1SB->Fit("f1_DsMix","L EMI S0","",DsDataFitRangeLow, DsDataFitRangeHigh);



	fitStatus=1;
  fitPrecision=1.e-8;
  while(fitStatus>0 && fitStatus!=4000){
    TFitter::SetPrecision(fitPrecision);
    fitResut=h1SB->Fit("f1_DsMix","L EMI S0","",DsDataFitRangeLow, DsDataFitRangeHigh);
    fitStatus=fitResut->Status();
    fitIsValid=fitResut->IsValid();
    cout<<"fit precision: "<<TFitter::GetPrecision()<<"   status: "<<fitStatus<<" isValid = "<< fitIsValid <<endl;
    if(fitStatus){
      fitPrecision *= 3;
    }
    if(fitPrecision>1.e-3) {break;}
  }




// fit status code : https://sft.its.cern.ch/jira/browse/ROOT-4445
// 4: MIGRAD error
// 40: MINOS error
// 400: HESSE error
// 4000: IMPROVE error

//	TF1 *f1_DsSignal = new TF1("f1_DsSignal","[0]*( [4]*TMath::Gaus(x,[1],[2]*(1+[5]))/( sqrt(2*TMath::Pi()) *[2]*(1+[5])  )   +(1-[4])*TMath::Gaus(x,[1],[3]*(1+[5]))/(sqrt(2*TMath::Pi())*[3]*(1+[5])  ) )"); // kTRUE nomalize Gaus
	f1_DsSignal->SetParameters(f1_DsMix->GetParameter(0), f1_DsMix->GetParameter(1) , f1_DsMix->GetParameter(2) , f1_DsMix->GetParameter(3) , f1_DsMix->GetParameter(4), f1_DsMix->GetParameter(5));
	f1_DsSignal->SetParError(0,f1_DsMix->GetParError(0));
	f1_DsSignal->SetParError(1,f1_DsMix->GetParError(1));
	f1_DsSignal->SetParError(2,f1_DsMix->GetParError(2));
	f1_DsSignal->SetParError(3,f1_DsMix->GetParError(3));
	f1_DsSignal->SetParError(4,f1_DsMix->GetParError(4));
	f1_DsSignal->SetParError(5,f1_DsMix->GetParError(5));

	f1_DsSignal->SetFillColor(kOrange-3);
	f1_DsSignal->SetFillStyle(3002);
	f1_DsSignal->SetLineColor(kOrange-3);
	f1_DsSignal->SetLineWidth(2);
	f1_DsSignal->SetLineStyle(2);
  f1_DsSignal->SetRange(DsDataFitRangeLow,DsDataFitRangeHigh); // must setrange before plot or get empty


	// double yieldDsSignal= f1_DsSignal->Integral(f1_DsSignal->GetParameter(1)-DsMassCandWidth,f1_DsSignal->GetParameter(1)+DsMassCandWidth)/ h_DsMassData->GetBinWidth(1);
	yieldDsSignal= f1_DsSignal->Integral(DsDataFitRangeLow,DsDataFitRangeHigh)/ h1SB->GetBinWidth(1);
	yieldDsSignalErr = f1_DsSignal->Integral(DsMass-DsMassCandWidth,DsMass+DsMassCandWidth)/h1SB->GetBinWidth(1) * f1_DsSignal->GetParError(0)/f1_DsSignal->GetParameter(0);

	yieldDsSignalTry = f1_DsSignal->GetParameter(0)/h1SB->GetBinWidth(1);
	yieldDsSignalTryErr = f1_DsSignal->GetParError(0)/h1SB->GetBinWidth(1);

	cout<<"parameter 0 = "<<f1_DsSignal->GetParameter(0)<<" BinWidth = "<< h_DsMassData->GetBinWidth(1) <<endl;
	cout<<"yield DsSignal = "<< yieldDsSignal <<" +- "<<yieldDsSignalErr<<endl;
	cout<<"yield DsSignal from Para 0 = "<< yieldDsSignalTry <<" +- "<<yieldDsSignalTryErr<<endl;

//	TF1 *f1_DsBkg = new TF1("f1_DsBkg"," [0]+[1]*x+[2]*x*x+[3]*x*x*x ");    // overall fit function for Data fitting, Ds , Dp ,background
//	f1_DsBkg->SetParameters(f1_DsMix->GetParameter(6), f1_DsMix->GetParameter(7) ,f1_DsMix->GetParameter(8),f1_DsMix->GetParameter(9));
	f1_DsBkg->SetParameters(f1_DsMix->GetParameter(6), f1_DsMix->GetParameter(7) ,f1_DsMix->GetParameter(8));
	// f1_DsBkg->SetParameters(f1_DsMix->GetParameter(6), f1_DsMix->GetParameter(7));
	f1_DsBkg->SetParError(0,f1_DsMix->GetParError(6));
	f1_DsBkg->SetParError(1,f1_DsMix->GetParError(7));
	f1_DsBkg->SetParError(2,f1_DsMix->GetParError(8));
//	f1_DsBkg->SetParError(3,f1_DsMix->GetParError(9));

	f1_DsBkg->SetLineColor(4);
	f1_DsBkg->SetLineStyle(2);
	f1_DsBkg->SetRange(DsDataFitRangeLow,DsDataFitRangeHigh);

	// double yieldBkg= f1_DsBkg->Integral(DsDataFitRangeLow,DsDataFitRangeHigh)/ h_DsMassData->GetBinWidth(1);
	yieldBkg= f1_DsBkg->Integral(DsDataFitRangeLow,DsDataFitRangeHigh)/ h1SB->GetBinWidth(1);

	cout<<"yield Bkg = " <<yieldBkg<<endl;

	TCanvas *c_projection = new TCanvas("c_projection","c_projection",c_wtopx,c_wtopy,c_W,c_H);
	c_projection->cd();
	SetCanvas(c_projection);

  h1SB->SetTitle("");
  h1SB->GetXaxis()->SetTitle("m_{KK#pi} (GeV/c^{2})");
  h1SB->GetXaxis()->CenterTitle();
	h1SB->GetYaxis()->SetTitle(Form("Events/(%.3f)",h1SB->GetBinWidth(2)));
  h1SB->GetYaxis()->CenterTitle();

	h1SB->SetMaximum(h1SB->GetMaximum()*1.4);
	h1SB->SetMinimum(h1SB->GetMinimum()*0.8);
  h1SB->Draw();
	TF1 *myfunNew=h1SB->GetFunction("f1_DsMix");
  myfunNew->Draw("same");
	f1_DsBkg->Draw("same"); //must setrange before draw or get empty
//	f1_DsSignal->Draw("same");

	shiftY=0.10;
	shiftX=0.44;
	// TLatex *tl_binfitData =new TLatex();
	tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("D_{S}^{+} + D_{S}^{-}"));  shiftY-=oneshift;
	tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%.0f < pt < %.0f GeV/c ",DptLow,DptHigh)); shiftY-=oneshift;
	tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("|y| < 1.0")); shiftY-=oneshift;
	if(extraName!="")
	{ 	
	tl_binfitData->DrawLatexNDC(textposx,textposy+shiftY,extraName); shiftY-=oneshift;
	}
	tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Raw Yield = %.0f #pm %.0f ",yieldDsSignal,yieldDsSignalErr ));  shiftY-=oneshift;
	// tl_binfitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Float width = %.3f",f1_DsMix->GetParameter(5))); shiftY-=oneshift;

	le_x1=0.64;
	le_x2=0.86;
	le_y1=0.42;
	le_y2=0.63;

	if(str_PbPb=="PbPb" && ibin_Dpt==4){
		  le_y1=0.17;
		  le_y2=0.38;
	}
	
	TLegend *le_binfitNew=new TLegend(le_x1,le_y1,le_x2,le_y2);
	le_binfitNew->SetBorderSize(0);
	le_binfitNew->AddEntry(h1SB,"Data","lp");
	le_binfitNew->AddEntry(myfunNew,"Signal+Background","l");
	le_binfitNew->AddEntry(f1_DsBkg,"Background","l");
	le_binfitNew->Draw("SAME");


	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);



  // texCmsPre->Draw("SAME");
  if(str_PbPb=="pp"){
  // texColpp->Draw("SAME");
	CMS_lumi(c_projection, 11, 10);
  }else{
  // texColPbPb->Draw("SAME");
	CMS_lumi(c_projection, 12, 10);
  }


	SavePlotDirs(c_projection,Form("%s_Data_fit_pt%.0fto%.0f_%s",str_PbPb.Data(),DptLow,DptHigh,extraName.Data() ),{"HLProjection","SignalFit",str_PbPb,"binfit"} );

	if(str_PbPb=="PbPb" && (ibin_Dpt<=3  )  ){	
	SavePlotDirs(c_binfit_Data[count_c_binfit],Form("%s_Data_fit_pt%.0fto%.0f_%s",str_PbPb.Data(),DptLow,DptHigh,extraName.Data() ),{"HLProjection","SignalFit",str_PbPb,"binfit"} );
}
*/

	  // return f1_DsMix; // this empty when draw

		count_c_binfit++;
	  return myfun; // this works

}

RooFitResult *unbinnedFit(TTree* t_Data,TTree* t_MC, TH1F* h_DsMassDsMC, TH1F* h_DsMassData, TH1F* h_RawRooFitYield, TH1F* h_YieldDcabin ,Int_t isPbPb=0, Int_t ibin_Dpt=1, Double_t DptLow=0 ,Double_t DptHigh=10 ,TString str_PbPb="pp", Int_t doBkgVariation=0 ,Int_t BkgFunction=0, Int_t doSigVariation=0,double FloatWidthVal=0 , double FloatWidthErr=0, double FloatWidth_DsMass=1.967, TString extraName="", Bool_t doDcaBinsFit= false , Bool_t doMCstudy= false, TString extraName2=""){

		// doSigVariation== 1: vary the float width , 2 : single gaussian



  double tex_upperY=0.95;

  TLatex* texCmsPre = new TLatex(0.15,tex_upperY, "#scale[1.25]{CMS} Projection");
  texCmsPre->SetNDC();
  texCmsPre->SetTextAlign(12);
  texCmsPre->SetTextSize(0.035);
  texCmsPre->SetTextFont(42);

  TLatex* texCmsSim = new TLatex(0.15,tex_upperY, "#scale[1.25]{CMS} Simulations");
  texCmsSim->SetNDC();
  texCmsSim->SetTextAlign(12);
  texCmsSim->SetTextSize(0.035);
  texCmsSim->SetTextFont(42);


  // TLatex* texColPbPb = new TLatex(0.88,tex_upperY, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  TLatex* texColPbPb = new TLatex(0.95,tex_upperY, "0.2 nb^{-1} (5.02 TeV PbPb)");
  texColPbPb->SetNDC();
  texColPbPb->SetTextAlign(32);
  texColPbPb->SetTextSize(0.035);
  texColPbPb->SetTextFont(42);

  // TLatex* texColpp = new TLatex(0.88,tex_upperY, "pp #sqrt{s} = 5.02 TeV");
  TLatex* texColpp = new TLatex(0.95,tex_upperY, "650 pb^{-1} (5.02 TeV pp )");
  texColpp->SetNDC();
  texColpp->SetTextAlign(32);
  texColpp->SetTextSize(0.035);
  texColpp->SetTextFont(42);

//  TLatex* texSimPbPb = new TLatex(0.88,tex_upperY, "530 #mub^{-1} (5.02 TeV PbPb)");
  TLatex* texSimPbPb = new TLatex(0.95,tex_upperY, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  texSimPbPb->SetNDC();
  texSimPbPb->SetTextAlign(32);
  texSimPbPb->SetTextSize(0.035);
  texSimPbPb->SetTextFont(42);

  TLatex* texSimpp = new TLatex(0.95,tex_upperY, "pp #sqrt{s} = 5.02 TeV");
  // TLatex* texSimpp = new TLatex(0.88,tex_upperY, "27.4 pb^{-1} (5.02 TeV pp )");
  texSimpp->SetNDC();
  texSimpp->SetTextAlign(32);
  texSimpp->SetTextSize(0.035);
  texSimpp->SetTextFont(42);






	RooWorkspace *ws = new RooWorkspace("workspace");

	if(doBkgVariation>=1){extraName+=Form("BkgFunVar%i",BkgFunction);}
	if(doSigVariation>=1){extraName+=Form("SigFunVar%i",doSigVariation);}

  vector<TString> savedirs;
  savedirs.push_back("SignalFit");
  savedirs.push_back(str_PbPb);
  savedirs.push_back("roofit");
  if(extraName!=""){ savedirs.push_back(extraName); }


	// load parameter later for saving doScan plot name 
	
	// process flow
	// 1.read binned fit result
	// 2.read data (input tree * historgram are already applied cut)
	// 3. do fit
	// 3a. maybe refit MC first ? 
	// 4. check fit

  // 1. read binned fit result
		
	// 1.1 read binned fit mc result

	Int_t nbin_DmassDraw=50;

	TF1 *f1_DsMC=(TF1*)h_DsMassDsMC->GetFunction("f1_DsMC");
	double N_sigV_MC=f1_DsMC->GetParameter(0)/h_DsMassDsMC->GetBinWidth(1);
	double DsMassMeanV_MC=f1_DsMC->GetParameter(1);
	double DsWidth1V_MC=f1_DsMC->GetParameter(2);
	double DsWidth2V_MC=f1_DsMC->GetParameter(3);
	double DsGaus1FrV_MC=f1_DsMC->GetParameter(4);

	// 1.2 read binned fit data result
  TF1 *f1_DsMix=(TF1*)h_DsMassData->GetFunction("f1_DsMix");
  double N_sigV = f1_DsMix->GetParameter(0)/h_DsMassData->GetBinWidth(1);
  double DsMassMeanV= f1_DsMix->GetParameter(1);
  double DsWidth1V= f1_DsMix->GetParameter(2);
  double DsWidth2V= f1_DsMix->GetParameter(3);
  double DsGaus1FrV= f1_DsMix->GetParameter(4);
  double DsFloatWidthV= f1_DsMix->GetParameter(5);
  double N_bkgV = f1_DsMix->GetParameter(6)/h_DsMassData->GetBinWidth(1);
  double Cheb1V= f1_DsMix->GetParameter(7);
  double Cheb2V= f1_DsMix->GetParameter(8);
  // double Cheb3V= f1_DsMix->GetParameter(9);

	cout<<"binned fit parameters "<<endl;
	cout<<"N_sigV = "<<N_sigV<<endl;
	cout<<"DsMassMeanV = "<<DsMassMeanV<<endl;
	cout<<"DsWidth1V = "<<DsWidth1V  <<endl;
	cout<<"DsWidth2V = "<<DsWidth2V  <<endl;
	cout<<"DsGaus1FrV = "<<DsGaus1FrV  <<endl;
	cout<<"DsFloatWidthV = "<<DsFloatWidthV  <<endl;
	cout<<"N_bkgV = "<<N_bkgV  <<endl;
	cout<<"Cheb1V = "<<Cheb1V  <<endl;
	// cout<<"Cheb2V = "<<Cheb2V  <<endl;
	// cout<<"Cheb3V = "<<Cheb3V  <<endl;

	// 2. read mc & data

  RooRealVar Dmass("Dmass","Dmass",DsDataFitRangeLow,DsDataFitRangeHigh);
  RooRealVar Ddca("Ddca","Ddca",0,0.1); // temp
	// RooRealVar BasicWeight("BasicWeight","BasicWeight",0,1e15);
	// RooDataSet RooDS_MC("RooDS_MC","RooDS_MC",RooArgSet(Dmass),WeightVar(BasicWeight),Import(*t_MC));
	RooRealVar D0DataWeight("D0DataWeight","D0DataWeight",0,1e15);
	RooDataSet RooDS_MC("RooDS_MC","RooDS_MC",RooArgSet(Dmass),WeightVar(D0DataWeight),Import(*t_MC));
	// RooDataSet RooDS_MC("RooDS_MC","RooDS_MC",RooArgSet(Dmass),Import(*t_MC));
	RooDS_MC.Print("v");

	RooDataSet RooDSAll("RooDSAll","RooDSAll",RooArgSet(Dmass,Ddca),Import(*t_Data));
	RooDataSet *RooDS=(RooDataSet*)RooDSAll.reduce(RooArgSet(Dmass));
	// RooDataSet RooDS("RooDS","RooDS",RooArgSet(Dmass),Import(*t_Data));
	RooDS->Print("v");
	
	cout<<"after load rooDs"<<endl;

	c_roofit_MC[count_c_rootfit] = new TCanvas(Form("c_roofit_MC_%i",count_c_rootfit),"c_roofit_MC",800,800);
	c_roofit_MC[count_c_rootfit]->cd();
	RooPlot* massframe_MC = Dmass.frame(Title("Dmass"));
	RooDS_MC.plotOn(massframe_MC);
	massframe_MC->Draw();	

 	c_roofit_Data_pull[count_c_rootfit] = new TCanvas(Form("c_roofit_Data_pull_%i",count_c_rootfit),"c_roofit_Data_pull",c_wtopx,c_wtopy,c_W,c_H);
 	c_roofit_Data_withpull[count_c_rootfit] = new TCanvas(Form("c_roofit_Data_withpull_%i",count_c_rootfit),"c_roofit_Data_withpull",c_wtopx,c_wtopy,c_W,800);
 	c_roofit_Data[count_c_rootfit] = new TCanvas(Form("c_roofit_Data_%i",count_c_rootfit),"c_roofit_Data",c_wtopx,c_wtopy,c_W,c_H);
	c_roofit_Data[count_c_rootfit]->cd();
	// RooPlot* massframe = Dmass.frame(Title("Dmass"));
	RooPlot* massframe=new RooPlot("massframe","Dmass",Dmass,DsDataFitRangeLow,DsDataFitRangeHigh,nbin_DmassDraw);
	RooDS->plotOn(massframe,DataError(RooAbsData::SumW2));
	massframe->Draw();

 // 3. fit

	// 3.1 fit mc
	//RooRealVar DsMassMean_MC("DsMassMean_MC","DsMassMean_MC",DsMassMeanV_MC,1.965,1.971);
	//RooRealVar DsWidth1_MC("DsWidth1_MC","DsWidth1_MC",DsWidth1V_MC,0.000001,0.1);
	//RooRealVar DsWidth2_MC("DsWidth2_MC","DsWidth2_MC",DsWidth2V_MC,0.000001,0.1);
	//RooRealVar DsGaus1Fr_MC("DsGaus1Fr_MC","DsGaus1Fr_MC",DsGaus1FrV_MC,0.0,1);

	RooRealVar DsMassMean_MC("DsMassMean_MC","DsMassMean_MC",1.9690,1.965,1.971);
	RooRealVar DsWidth1_MC("DsWidth1_MC","DsWidth1_MC",0.011,0.008,0.02);
	RooRealVar DsWidth2_MC("DsWidth2_MC","DsWidth2_MC",0.0060,0.001,0.01);
	RooRealVar DsGaus1Fr_MC("DsGaus1Fr_MC","DsGaus1Fr_MC",0.65,0.3,1);

	RooGaussian Gaus1_MC("Gaus1_MC","gauss(Dmass,DsMassMean_MC,DsWidth1_MC)",Dmass,DsMassMean_MC,DsWidth1_MC);
	RooGaussian Gaus2_MC("Gaus2_MC","gauss(Dmass,DsMassMean_MC,DsWidth2_MC)",Dmass,DsMassMean_MC,DsWidth2_MC);
	RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(Gaus1_MC,Gaus2_MC),DsGaus1Fr_MC);


	RooFitResult *fitresult_MC;
	
	if(doSigVariation<=1){
	SigPdf_MC.fitTo(RooDS_MC,NumCPU(10));
	fitresult_MC = SigPdf_MC.fitTo(RooDS_MC,Save(),Range(1.935,2.005),NumCPU(10));
	fitresult_MC->Print("v");
	SigPdf_MC.plotOn(massframe_MC,LineColor(2));
	}else if(doSigVariation==2){ 
	Gaus1_MC.fitTo(RooDS_MC,NumCPU(10));
	fitresult_MC = Gaus1_MC.fitTo(RooDS_MC,Save(),Range(1.935,2.005),NumCPU(10));
	fitresult_MC->Print("v");
	Gaus1_MC.plotOn(massframe_MC,LineColor(2));
	}

	gStyle->SetOptFit(1111);
	gStyle->SetOptStat(1111);

	c_roofit_MC[count_c_rootfit]->cd();	

	massframe_MC->SetTitle("");
	massframe_MC->GetXaxis()->SetTitle("m_{KK#pi} (GeV)");
	massframe_MC->GetXaxis()->CenterTitle();
	massframe_MC->GetYaxis()->CenterTitle();
	massframe_MC->Draw();



  shiftY=0.05;
  TLatex *tl_roofitMC =new TLatex();
  // tl_roofitMC->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s_MC_unbinfit",str_PbPb.Data()));  shiftY-=oneshift;
  tl_roofitMC->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%.0f < p_{T} < %.0f GeV",DptLow,DptHigh)); shiftY-=oneshift;
  if(extraName!="" || extraName2!="")
  {
  tl_roofitMC->DrawLatexNDC(textposx+shiftX,textposy+shiftY,extraName+extraName2); shiftY-=oneshift;
  }
	if(doSigVariation<=1){
  tl_roofitMC->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Mean=%.4f, frGaus1=%.3f  ",DsMassMean_MC.getValV() , DsGaus1Fr_MC.getValV() ));  shiftY-=oneshift;
  tl_roofitMC->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Sigma1=%.4f, Sigma2=%.4f", DsWidth1_MC.getValV(), DsWidth2_MC.getValV()  ));  shiftY-=oneshift;
	}else if(doSigVariation==2){
  tl_roofitMC->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Mean=%.4f, Sigma=%.3f  ",DsMassMean_MC.getValV() , DsWidth1_MC.getValV() ));  shiftY-=oneshift;
	}

	texCmsSim->Draw("SAME");
	if(str_PbPb=="pp"){
	texSimpp->Draw("SAME");
	}else{
	texSimPbPb->Draw("SAME");
	}

  SavePlotDirs(c_roofit_MC[count_c_rootfit],Form("%s_MC_unbinfit_pt%.0fto%.0f_%s%s",str_PbPb.Data(),DptLow,DptHigh,extraName.Data(),extraName2.Data() ),savedirs );
  // SavePlotDirs(c_roofit_MC[count_c_rootfit],Form("%s_MC_unbinfit_pt%.0fto%.0f_%s%s",str_PbPb.Data(),DptLow,DptHigh,extraName.Data(),extraName2.Data() ),{"SignalFit",str_PbPb,"roofit"} );
	
	// c_roofit_MC->SaveAs("roofit.pdf"); // this works

	// cout<<"NumSig="<<NumSig.getValV()<<endl;

	// 3.2 fit data

	DsMassMeanV=DsMassMean_MC.getValV();
	DsWidth1V=DsWidth1_MC.getValV();
	DsWidth2V=DsWidth2_MC.getValV();
	DsGaus1FrV=DsGaus1Fr_MC.getValV();

	RooRealVar DsMassMean("DsMassMean","DsMassMean",DsMassMeanV,1.965,1.971);
	RooRealVar DsWidth1("DsWidth1","DsWidth1",DsWidth1V,0.000001,0.1);
	RooRealVar DsWidth2("DsWidth2","DsWidth2",DsWidth2V,0.000001,0.1);
	RooRealVar DsGaus1Fr("DsGaus1Fr","DsGaus1Fr",DsGaus1FrV,0.0,1);
	RooRealVar DsFloatWidth("DsFloatWidth","DsFloatWidth",DsFloatWidthV,-0.6,0.6);
	// RooRealVar DsFloatWidth("DsFloatWidth","DsFloatWidth",DsFloatWidthV,DsFloatWidthV,DsFloatWidthV);

	DsWidth1.setConstant(kTRUE);
	DsWidth2.setConstant(kTRUE);
	DsGaus1Fr.setConstant(kTRUE);
	// DsFloatWidth.setConstant(kTRUE);
	if(doSigVariation==1){
		DsFloatWidth.setVal(FloatWidthVal+FloatWidthErr);
		DsFloatWidth.setConstant(kTRUE);
		DsMassMean.setVal(FloatWidth_DsMass);
		DsMassMean.setConstant(kTRUE);
	}


	RooFormulaVar scale_width1("scale width1","scaled width1","DsWidth1*(1+DsFloatWidth)",RooArgSet(DsWidth1,DsFloatWidth));
	RooFormulaVar scale_width2("scale width2","scaled width2","DsWidth2*(1+DsFloatWidth)",RooArgSet(DsWidth2,DsFloatWidth));
	RooGaussian Gaus1("Gaus1","gauss(Dmass,DsMassMean,scale_width1)",Dmass,DsMassMean,scale_width1);
	RooGaussian Gaus2("Gaus2","gauss(Dmass,DsMassMean,scale_width2)",Dmass,DsMassMean,scale_width2);
	RooAddPdf SigPdf("SigPdf","SigPdf",RooArgList(Gaus1,Gaus2),DsGaus1Fr);
	
 	RooRealVar Cheb1("Cheb1","Cheb1",0,-1,1); // no input is better than with input
	RooRealVar Cheb2("Cheb2","Cheb2",0,-1,1);
	RooRealVar Cheb3("Cheb3","Cheb3",0,-1,1);

	// RooRealVar Cheb1("Cheb1","Cheb1",Cheb1V,-1,1);
	// RooRealVar Cheb2("Cheb2","Cheb2",Cheb2V,-1,1);
	// RooRealVar Cheb3("Cheb3","Cheb3",Cheb3V,-1,1);


	// RooChebychev BkgPdf("BkgPdf","BkgPdf",Dmass,RooArgList(Cheb1,Cheb2,Cheb3));

	RooChebychev *BkgPdf;

	if(doBkgVariation==1 && (DptLow>=20 || DptLow<=3)  ){
		BkgPdf=new RooChebychev("BkgPdf","BkgPdf",Dmass,RooArgList(Cheb1));
	}else if(doBkgVariation==1){
		BkgPdf=new RooChebychev("BkgPdf","BkgPdf",Dmass,RooArgList(Cheb1,Cheb2,Cheb3));
	}else{	
		BkgPdf=new RooChebychev("BkgPdf","BkgPdf",Dmass,RooArgList(Cheb1,Cheb2));
	}

	// RooChebychev BkgPdf("BkgPdf","BkgPdf",Dmass,RooArgList(Cheb1));
	RooRealVar NumSig("NumSig","Number of Signal",N_sigV,-1e4,4e5);
	RooRealVar NumBkg("NumBkg","Number of Background",N_bkgV,0,1e8);

	RooAddPdf RooDsMixPdf("RooDsMixPdf","RooDsMixPdf",RooArgList(SigPdf,*BkgPdf),RooArgList(NumSig,NumBkg));	
	RooAddPdf RooDsMixPdf_1G("RooDsMixPdf_1G","RooDsMixPdf_1G",RooArgList(Gaus1,*BkgPdf),RooArgList(NumSig,NumBkg));	

	Dmass.setRange("range1",1.91,2.04);
	Dmass.setRange("Allrange",1.91,2.11);

	// if(str_PbPb=="pp"&& (ibin_Dpt==0 || ibin_Dpt==2 || ibin_Dpt==3) ){
		// RooDsMixPdf.fitTo(RooDS,Extended(kTRUE),Range("range1"));
	// }
  // if(str_PbPb=="PbPb"&& ibin_Dpt==2){
	  // RooDsMixPdf.fitTo(RooDS,Extended(kTRUE),Range("range1"));
	// }

RooFitResult *fitresult=NULL;
RooHist *hpull=NULL;
if(doSigVariation<=1){
	// RooFitResult *fitresult = RooDsMixPdf.fitTo(RooDS,Extended(kTRUE),Save(),Range(1.91,2.06));
	RooDsMixPdf.fitTo(*RooDS,Extended(kTRUE),NumCPU(20));
	RooDsMixPdf.fitTo(*RooDS,Extended(kTRUE),NumCPU(20));
	fitresult = RooDsMixPdf.fitTo(*RooDS,Extended(kTRUE),Range("Allrange"),Save(),NumCPU(20));
	// RooFitResult *fitresult = RooDsMixPdf.fitTo(RooDS,Extended(kTRUE),Range("range1"),Save());

	fitresult->Print("v");
	
	RooDsMixPdf.plotOn(massframe,LineColor(2));
	hpull=massframe->pullHist();
	// RooDsMixPdf.plotOn(massframe,Components(SigPdf),DrawOption("F"),LineColor(kOrange-3),FillColor(kOrange-3),FillStyle(3002),MoveToBack());
	RooDsMixPdf.plotOn(massframe,Components(*BkgPdf),LineColor(4),LineStyle(kDashed));

}else if(doSigVariation==2){
	DsFloatWidth.setVal(0);
	DsFloatWidth.setConstant(kTRUE);
	// DsWidth1.setConstant(kFALSE);
	RooDsMixPdf_1G.fitTo(*RooDS,Extended(kTRUE),NumCPU(20));
	RooDsMixPdf_1G.fitTo(*RooDS,Extended(kTRUE),NumCPU(20));
	fitresult = RooDsMixPdf_1G.fitTo(*RooDS,Extended(kTRUE),Range("Allrange"),Save(),NumCPU(20));

	fitresult->Print("v");
	
	RooDsMixPdf_1G.plotOn(massframe,LineColor(2));
	hpull=massframe->pullHist();
	RooDsMixPdf_1G.plotOn(massframe,Components(*BkgPdf),LineColor(4),LineStyle(kDashed));

}

	// RooHist* hpull=massframe->pullHist();
  RooPlot *massframe_pull=new RooPlot("massframe_pull","Dmass",Dmass,DsDataFitRangeLow,DsDataFitRangeHigh,nbin_DmassDraw);
	massframe_pull->addPlotable(hpull,"P");

  c_roofit_Data_pull[count_c_rootfit]->cd();
	c_roofit_Data_pull[count_c_rootfit]->cd();
  c_roofit_Data_pull[count_c_rootfit]->SetFillColor(0);
  c_roofit_Data_pull[count_c_rootfit]->SetBorderMode(0);
  c_roofit_Data_pull[count_c_rootfit]->SetFrameFillStyle(0);
  c_roofit_Data_pull[count_c_rootfit]->SetFrameBorderMode(0);
  c_roofit_Data_pull[count_c_rootfit]->SetLeftMargin( c_Lmg );
  c_roofit_Data_pull[count_c_rootfit]->SetRightMargin( c_Rmg );
  c_roofit_Data_pull[count_c_rootfit]->SetTopMargin( c_Tmg );
  c_roofit_Data_pull[count_c_rootfit]->SetBottomMargin( c_Bmg );
  c_roofit_Data_pull[count_c_rootfit]->SetTickx(0);
  c_roofit_Data_pull[count_c_rootfit]->SetTicky(0);


	TLine *Lm1=new TLine(DsDataFitRangeLow,-1,DsDataFitRangeHigh,-1);
	Lm1->SetLineColor(3);
	Lm1->SetLineStyle(7);
	Lm1->SetLineWidth(2);

	TLine *Lp1=new TLine(DsDataFitRangeLow,1,DsDataFitRangeHigh,1);
	Lp1->SetLineColor(3);
	Lp1->SetLineStyle(7);
	Lp1->SetLineWidth(2);

  massframe_pull->SetTitle("");
  massframe_pull->GetXaxis()->SetTitle("m_{KK#pi} (GeV)");
  massframe_pull->GetXaxis()->CenterTitle();
	massframe_pull->GetYaxis()->SetTitle("pull");
  massframe_pull->GetYaxis()->CenterTitle();
  massframe_pull->Draw();

	Lm1->Draw("same");
	Lp1->Draw("same");

	SavePlotDirs(c_roofit_Data_pull[count_c_rootfit],Form("%s_Data_unbinfit_pull_pt%.0fto%.0f_%s%s",str_PbPb.Data(),DptLow,DptHigh,extraName.Data(),extraName2.Data() ),savedirs );

	cout<<"massframe GetMaximum() = "<<massframe->GetMaximum()<<" , min = "<<massframe->GetMinimum()<<endl;
	cout<<"h_DsMassData Maximum() = "<<h_DsMassData->GetMaximum()<<" , min = "<<h_DsMassData->GetMinimum()<<" nbin = "<<h_DsMassData->GetNbinsX()<<endl;

	int nbin_hdata=h_DsMassData->GetNbinsX();
	double hmax=h_DsMassData->GetMaximum()*(double)nbin_hdata/(double)nbin_DmassDraw;
	double hmin=h_DsMassData->GetMinimum()*(double)nbin_hdata/(double)nbin_DmassDraw;
	double hdiff=hmax-hmin;
	double maxExtra=0.3;
	double minExtra=0.15;
	double newMax=hmax+maxExtra*hdiff;
	double newMin=hmin-minExtra*hdiff; if(newMin<0){ newMin=0 ; }
	cout<<"newMax = "<<newMax<<" , newMin = "<<newMin<<endl;

	massframe->SetMaximum(newMax);// must put here after all plotOn
	massframe->SetMinimum(newMin);// must put here after all plotOn
	// massframe->GetYaxis()->SetRangeUser(200,1000);

	c_roofit_Data[count_c_rootfit]->cd();
  c_roofit_Data[count_c_rootfit]->SetFillColor(0);
  c_roofit_Data[count_c_rootfit]->SetBorderMode(0);
  c_roofit_Data[count_c_rootfit]->SetFrameFillStyle(0);
  c_roofit_Data[count_c_rootfit]->SetFrameBorderMode(0);
  c_roofit_Data[count_c_rootfit]->SetLeftMargin( c_Lmg );
  c_roofit_Data[count_c_rootfit]->SetRightMargin( c_Rmg );
  c_roofit_Data[count_c_rootfit]->SetTopMargin( c_Tmg );
  c_roofit_Data[count_c_rootfit]->SetBottomMargin( c_Bmg );
  c_roofit_Data[count_c_rootfit]->SetTickx(0);
  c_roofit_Data[count_c_rootfit]->SetTicky(0);

	massframe->SetTitle("");
	massframe->GetXaxis()->SetTitle("m_{KK#pi} (GeV)");
	massframe->GetXaxis()->CenterTitle();
	massframe->GetYaxis()->CenterTitle();
	massframe->Draw();


	// cout<<"massframe GetMaximum() = "<<massframe->GetMaximum()<<" , min = "<<massframe->GetMinimum()<<endl;

	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
  shiftY=0.05;
	shiftX=0.35;
  TLatex *tl_roofitData =new TLatex();
  // tl_roofitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s_Data_unbinfit",str_PbPb.Data()));  shiftY-=oneshift;
  tl_roofitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%.0f < p_{T} < %.0f GeV",DptLow,DptHigh)); shiftY-=oneshift;
  if(extraName!=""|| extraName2!="")
  {
  tl_roofitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,extraName+extraName2); shiftY-=oneshift;
  }
  tl_roofitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Raw Yield = %.1f #pm %.1f",NumSig.getValV(),NumSig.getError() ));  shiftY-=oneshift;
	if(doSigVariation<=1){
	// tl_roofitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Float Width = %.3f #pm %.3f ", DsFloatWidth.getValV(), DsFloatWidth.getError() ));  shiftY-=oneshift;
	}
  tl_roofitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,"|y|<1.0"); shiftY-=oneshift;
	// TString *savedirs={"SignalFit",str_PbPb,"roofit"};
	// vector<TString> savedirs;
	// savedirs.push_back("SignalFit");
	// savedirs.push_back(str_PbPb);
	// savedirs.push_back("roofit");
	// if(extraName!=""){ savedirs.push_back(extraName); }
	// SavePlotDirs(c_roofit_Data[count_c_rootfit],Form("%s_Data_unbinfit_pt%.0fto%.0f_%s",str_PbPb.Data(),DptLow,DptHigh,extraName.Data() ),{"SignalFit",str_PbPb,"roofit"} );

	texCmsPre->Draw("SAME");
	if(str_PbPb=="pp"){
	texColpp->Draw("SAME");
	}else{
	texColPbPb->Draw("SAME");
	}	

	SavePlotDirs(c_roofit_Data[count_c_rootfit],Form("%s_Data_unbinfit_pt%.0fto%.0f_%s%s",str_PbPb.Data(),DptLow,DptHigh,extraName.Data(),extraName2.Data() ),savedirs );


	// same plot for data & pull

	c_roofit_Data_withpull[count_c_rootfit]->cd();
	c_roofit_Data_withpull[count_c_rootfit]->cd();
  c_roofit_Data_withpull[count_c_rootfit]->SetFillColor(0);
  c_roofit_Data_withpull[count_c_rootfit]->SetBorderMode(0);
  c_roofit_Data_withpull[count_c_rootfit]->SetFrameFillStyle(0);
  c_roofit_Data_withpull[count_c_rootfit]->SetFrameBorderMode(0);
  c_roofit_Data_withpull[count_c_rootfit]->SetLeftMargin( c_Lmg );
  c_roofit_Data_withpull[count_c_rootfit]->SetRightMargin( c_Rmg );
  c_roofit_Data_withpull[count_c_rootfit]->SetTopMargin( c_Tmg );
  c_roofit_Data_withpull[count_c_rootfit]->SetBottomMargin( c_Bmg );
  c_roofit_Data_withpull[count_c_rootfit]->SetTickx(0);
  c_roofit_Data_withpull[count_c_rootfit]->SetTicky(0);


	TPad *pad1 = new TPad("pad1","top pad",0.0,0.25,1.0,1.0);
	pad1->SetTopMargin(0.08);
	pad1->SetBottomMargin(0.0);
	pad1->SetRightMargin(c_Rmg);
	pad1->SetLeftMargin(c_Lmg);
	pad1->Draw();
	TPad *pad2 = new TPad("pad2","bottom pad",0.0,0.00,1.0,0.25);
	pad2->SetTopMargin(0.0);
	pad2->SetBottomMargin(0.30);
	pad2->SetRightMargin(c_Rmg);
	pad2->SetLeftMargin(c_Lmg);
	pad2->Draw();
	
	pad1->cd();

	
	massframe->GetYaxis()->SetTitleSize(0.05);
	massframe->GetYaxis()->SetTitleOffset(1.15);
	massframe->GetYaxis()->SetLabelSize(0.05);	
	massframe->Draw();
  shiftY=0.05;
	shiftX=0.35;
  tl_roofitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%.0f < p_{T} < %.0f GeV",DptLow,DptHigh)); shiftY-=oneshift;
  if(extraName!=""|| extraName2!="")
  {
  tl_roofitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,extraName+extraName2); shiftY-=oneshift;
  }
  tl_roofitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Raw Yield = %.1f #pm %.1f",NumSig.getValV(),NumSig.getError() ));  shiftY-=oneshift;
	if(doSigVariation<=1){
	}
  tl_roofitData->DrawLatexNDC(textposx+shiftX,textposy+shiftY,"|y|<1.0"); shiftY-=oneshift;

	texCmsPre->Draw("SAME");
	if(str_PbPb=="pp"){
	texColpp->Draw("SAME");
	}else{
	texColPbPb->Draw("SAME");
	}	

	pad2->cd();
	massframe_pull->GetXaxis()->SetTitleSize(0.14);
	massframe_pull->GetXaxis()->SetLabelSize(0.11);	
	massframe_pull->GetYaxis()->SetTitleSize(0.14);
	massframe_pull->GetYaxis()->SetTitleOffset(0.4);
	massframe_pull->GetYaxis()->SetLabelSize(0.11);	
  massframe_pull->Draw();

	Lm1->Draw("same");
	Lp1->Draw("same");

	SavePlotDirs(c_roofit_Data_withpull[count_c_rootfit],Form("%s_Data_unbinfit_withpull_pt%.0fto%.0f_%s%s",str_PbPb.Data(),DptLow,DptHigh,extraName.Data(),extraName2.Data() ),savedirs );





	
	// c_roofit_Data->SaveAs("roofit.pdf"); // this works

	cout<<"NumSig="<<NumSig.getValV()<<" +- "<<NumSig.getError()<<endl;
	cout<<"NumBkg="<<NumBkg.getValV()<<endl;

	double NumSigVal=NumSig.getValV(); double NumSigErr=NumSig.getError();

	h_RawRooFitYield->SetBinContent(1,NumSig.getValV());
	h_RawRooFitYield->SetBinError(1,NumSig.getError());

 if(doBkgVariation==0 && doSigVariation ==0 )
	{ 
		FloatWidthVal_Norm=DsFloatWidth.getValV();   
		FloatWidthErr_Norm=DsFloatWidth.getError();  
		FloatWidth_DsMass_Norm=DsMassMean.getValV();
	} // save for SigFunVar fit

	// mcStudy might be very time consuming..
	if(doMCstudy){
		ws->import(DsMassMean);	
		ws->import(DsWidth1);	
		ws->import(DsWidth2);	
		ws->import(DsGaus1Fr);	
		ws->import(DsFloatWidth);	
		ws->import(Cheb1);
		ws->import(Cheb2);
		ws->import(scale_width1);	
		ws->import(scale_width2);	
		ws->import(NumSig);	
		ws->import(NumBkg);	
		ws->import(RooDsMixPdf);	
	

		TFile *ws_test= new TFile("ws_test.root","RECREATE");
 		ws->writeToFile("ws_test.root",kTRUE);

		RooMCStudy *mcstudy = new RooMCStudy(RooDsMixPdf,Dmass,Silence(),Extended(), FitOptions(Save(kTRUE),PrintEvalErrors(0),NumCPU(20)));
  // G e n e r a t e   a n d   f i t   e v e n t s
  // ---------------------------------------------

  // Generate and fit 1000 samples of Poisson(nExpected) events
  mcstudy->generateAndFit(100) ;



  // E x p l o r e   r e s u l t s   o f   s t u d y
  // ------------------------------------------------

  // Make plots of the distributions of mean, the error on mean and the pull of mean
  RooPlot* frame1 = mcstudy->plotParam(NumSig,Bins(40)) ;
  RooPlot* frame2 = mcstudy->plotError(NumSig,Bins(40)) ;
  RooPlot* frame3 = mcstudy->plotPull(NumSig,Bins(40),FitGauss(kTRUE)) ;

  // Plot distribution of minimized likelihood
  RooPlot* frame4 = mcstudy->plotNLL(Bins(40)) ;
/*
  // Make some histograms from the parameter dataset
  TH1* hh_cor_a0_s1f = mcstudy->fitParDataSet().createHistogram("hh",Cheb2,YVar(DsFloatWidth)) ;
  TH1* hh_cor_a0_a1  = mcstudy->fitParDataSet().createHistogram("hh",Cheb1,YVar(Cheb2)) ;

  // Access some of the saved fit results from individual toys
  TH2* corrHist000 = mcstudy->fitResult(0)->correlationHist("c000") ;
  TH2* corrHist127 = mcstudy->fitResult(27)->correlationHist("c27") ;
  TH2* corrHist953 = mcstudy->fitResult(53)->correlationHist("c53") ;
*/


  // Draw all plots on a canvas
  gStyle->SetPalette(1) ;
  gStyle->SetOptStat(0) ;
  c_roofit_MCstudy[count_c_rootfit] = new TCanvas(Form("c_roofit_MCstudy_%i",count_c_rootfit),Form("c_roofit_MCstudy_%i",count_c_rootfit),1000,1000) ;
  c_roofit_MCstudy[count_c_rootfit]->Divide(2,2) ;

	c_roofit_MCstudy[count_c_rootfit]->cd(1); gPad->SetLeftMargin(0.15);
	massframe->Draw();
  shiftY=0;
  // TLatex *tl_roofitData =new TLatex();
  tl_roofitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s_Data_unbinfit",str_PbPb.Data()));  shiftY-=oneshift;
  tl_roofitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%.0f < pt < %.0f ",DptLow,DptHigh)); shiftY-=oneshift;
  if(extraName!="")
  {
  tl_roofitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,extraName); shiftY-=oneshift;
  }
  tl_roofitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Raw Yield = %.1f #pm %.1f",NumSigVal,NumSigErr ));  shiftY-=oneshift;
  // tl_roofitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Float Width = %.3f #pm %.3f ", DsFloatWidth.getValV(), DsFloatWidth.getError() ));  shiftY-=oneshift;


  c_roofit_MCstudy[count_c_rootfit]->cd(2) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
  c_roofit_MCstudy[count_c_rootfit]->cd(3) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;
  c_roofit_MCstudy[count_c_rootfit]->cd(4) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;
  // c_roofit_MCstudy[count_c_rootfit]->cd(4) ; gPad->SetLeftMargin(0.15) ; frame4->GetYaxis()->SetTitleOffset(1.4) ; frame4->Draw() ;
//  c->cd(5) ; gPad->SetLeftMargin(0.15) ; hh_cor_a0_s1f->GetYaxis()->SetTitleOffset(1.4) ; hh_cor_a0_s1f->Draw("box") ;
//  c->cd(6) ; gPad->SetLeftMargin(0.15) ; hh_cor_a0_a1->GetYaxis()->SetTitleOffset(1.4) ; hh_cor_a0_a1->Draw("box") ;
//  c->cd(7) ; gPad->SetLeftMargin(0.15) ; corrHist000->GetYaxis()->SetTitleOffset(1.4) ; corrHist000->Draw("colz") ;
//  c->cd(8) ; gPad->SetLeftMargin(0.15) ; corrHist127->GetYaxis()->SetTitleOffset(1.4) ; corrHist127->Draw("colz") ;
//  c->cd(9) ; gPad->SetLeftMargin(0.15) ; corrHist953->GetYaxis()->SetTitleOffset(1.4) ; corrHist953->Draw("colz") ;
//	c->SaveAs("plots/mcstudy.png");
//	c->SaveAs("plots/mcstudy.pdf");
//	c->SaveAs("plots/mcstudy.C");

	// TString *savedirs={"SignalFit",str_PbPb,"roofit"};
	vector<TString> savedirs;
	savedirs.push_back("SignalFit");
	savedirs.push_back(str_PbPb);
	savedirs.push_back("roofit");
	if(extraName!=""){ savedirs.push_back(extraName); }
	// SavePlotDirs(c_roofit_Data[count_c_rootfit],Form("%s_Data_unbinfit_pt%.0fto%.0f_%s",str_PbPb.Data(),DptLow,DptHigh,extraName.Data() ),{"SignalFit",str_PbPb,"roofit"} );
	SavePlotDirs(c_roofit_MCstudy[count_c_rootfit],Form("%s_Data_unbinfit_pt%.0fto%.0f_MCStudy",str_PbPb.Data(),DptLow,DptHigh),{"SignalFit",str_PbPb,"roofit","MCStudy"} );


  // Make RooMCStudy object available on command line after
  // macro finishes
  gDirectory->Add(mcstudy) ;

	
	
	}// end mcStudy


 

	if(doDcaBinsFit){
		// RooDataSet RooDSDcaAll
		RooDataSet *RooDSDca[nbin_dca];
		RooFitResult *fitresultDca[nbin_dca];
		RooPlot* massframeDca[nbin_dca];

		for(int i=0; i<nbin_dca; i++){
			double dcaLow=bins_dca[i];
			double dcaHigh=bins_dca[i+1];
  		// RooDSDca[i]= new RooDataSet(Form("RooDSDca_%i",i),Form("RooDSDca_%i",i),RooArgSet(Dmass),Import(*t_Data),Cut(Form("Ddca>%f && Ddca<%f",dcaLow,dcaHigh)));
  		RooDSDca[i]=(RooDataSet*)RooDSAll.reduce(RooArgSet(Dmass),Form("Ddca>%f && Ddca<%f",dcaLow,dcaHigh));
	  	RooDSDca[i]->Print("v");

			DsFloatWidth.setConstant(kTRUE);
			DsMassMean.setConstant(kTRUE);
			RooDsMixPdf.fitTo(*RooDSDca[i],Extended(kTRUE),NumCPU(20));
			fitresultDca[i]=RooDsMixPdf.fitTo(*RooDSDca[i],Extended(kTRUE),NumCPU(20),Save());
  		fitresultDca[i]->Print("v");

			cout<<"idca bin = "<<i<<endl;
			cout<<"NumSig="<<NumSig.getValV()<<" +- "<<NumSig.getError()<<endl;
			cout<<"NumBkg="<<NumBkg.getValV()<<endl;
			// plot & save result
			
			h_YieldDcabin->SetBinContent(i+1, NumSig.getValV());
			h_YieldDcabin->SetBinError(i+1, NumSig.getError());


			massframeDca[i]=new RooPlot("massframe","Dmass",Dmass,DsDataFitRangeLow,DsDataFitRangeHigh,nbin_DmassDraw);
			RooDSDca[i]->plotOn(massframeDca[i]);

  		RooDsMixPdf.plotOn(massframeDca[i],LineColor(2));
	  	// RooDsMixPdf.plotOn(massframeDca[i],Components(SigPdf),DrawOption("F"),LineColor(kOrange-3),FillColor(kOrange-3),FillStyle(3002),MoveToBack());
		  RooDsMixPdf.plotOn(massframeDca[i],Components(*BkgPdf),LineColor(4),LineStyle(kDashed));

			c_roofit_Data_Dca[count_c_rootfit][i]=new TCanvas(Form("c_roofit_Data_Dca_%i_%i",count_c_rootfit,i), Form("c_roofit_Data_Dca_%i_%i",count_c_rootfit,i), 800,800 );
			c_roofit_Data_Dca[count_c_rootfit][i]->cd();	
			massframeDca[i]->Draw();

      shiftY=0;
      TLatex *tl_roofitDataDca =new TLatex();
      tl_roofitDataDca->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s_Data_unbinfit",str_PbPb.Data()));  shiftY-=oneshift;
      tl_roofitDataDca->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%.0f < pt < %.0f ",DptLow,DptHigh)); shiftY-=oneshift;
      tl_roofitDataDca->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%.4f < DCA < %.4f ",dcaLow,dcaHigh)); shiftY-=oneshift;
      if(extraName!="")
      {
      tl_roofitDataDca->DrawLatexNDC(textposx,textposy+shiftY,extraName); shiftY-=oneshift;
      }
      tl_roofitDataDca->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Raw Yield = %.1f #pm %.1f",NumSig.getValV(),NumSig.getError() ));  shiftY-=oneshift;
      tl_roofitDataDca->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Float Width = %.3f #pm %.3f ", DsFloatWidth.getValV(), DsFloatWidth.getError() ));  shiftY-=oneshift;


  		SavePlotDirs(c_roofit_Data_Dca[count_c_rootfit][i],Form("%s_Data_unbinfit_pt%.0fto%.0f_dcabin%i_%s",str_PbPb.Data(),DptLow,DptHigh,i,extraName.Data() ),{"SignalFit",str_PbPb,"roofit","Dcabins"} );

		
		}

	}




	count_c_rootfit++;

	return fitresult;
}


int SignalFit(int isPbPb=3, int ibin_Dpt=2,bool doMCstudySet=0, int doBkgVariationSet=0, int doSigVariationSet=0, int doCutScanSet=0){

	cout<<"isPbPb = "<<isPbPb<<" , ibin_Dpt = "<<ibin_Dpt<<" , doBkgVariationSet = "<<doBkgVariationSet<<" , doSigVariationSet = "<<doSigVariationSet<<" , doCutScanSet = "<<doCutScanSet<<" , doMCstudySet = "<<doMCstudySet<<endl;

	if(doMCstudySet){ // when doing the MC study, disable others 
	doBkgVariationSet=0;
	doSigVariationSet=0;
	doCutScanSet=0;
	}

	InitStyle();
	initParameter();
	setTDRStyle();

	int nbin_pt=nbin_pt_pp;
	double *bins_pt=bins_pt_pp;
 
 	double DptLow  = bins_pt[ibin_Dpt];
	double DptHigh = bins_pt[ibin_Dpt+1];

	TString dataName="./FitFile_pp_HLproject.root";
	TString mcName_Prompt="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output/MC_eff_pp_Prompt_phikkpi.root";
	TString mcName_NonPrompt="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output/MC_eff_pp_NonPrompt_phikkpi.root";
	TString outName=Form("output/FitResult_pp_pt%.0fto%.0f.root",DptLow,DptHigh);
	TString str_PbPb="pp";

	if(doMCstudySet){ // avoid overwritten other fit result
		outName=Form("output/MCstudy/FitResult_pp_pt%.0fto%.0f.root",DptLow,DptHigh);
	}

  double *Dchi2clMinScan_bins=Dchi2clMinScan_bins_pp;
	double *DalphaMaxScan_bins=DalphaMaxScan_bins_pp;
	double *DdlsMinScan_bins=DdlsMinScan_bins_pp;



	if(isPbPb==3){
		nbin_pt=nbin_pt_PbPb3;
		bins_pt=bins_pt_PbPb3;

 	  DptLow  = bins_pt[ibin_Dpt];
	  DptHigh = bins_pt[ibin_Dpt+1];

		dataName="./FitFile_PbPb3_HLproject.root";
		// dataName="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree/output/FitFile_PbPb3.root";
		mcName_Prompt="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output/MC_eff_PbPb3_Prompt_phikkpi.root";
		mcName_NonPrompt="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output/MC_eff_PbPb3_NonPrompt_phikkpi.root";
		outName=Form("output/FitResult_PbPb3_pt%.0fto%.0f.root",DptLow,DptHigh);
		str_PbPb="PbPb";


	if(doMCstudySet){
		outName=Form("output/MCstudy/FitResult_PbPb3_pt%.0fto%.0f.root",DptLow,DptHigh);
	}

    Dchi2clMinScan_bins=Dchi2clMinScan_bins_PbPb3;
		DalphaMaxScan_bins=DalphaMaxScan_bins_PbPb3;
		DdlsMinScan_bins=DdlsMinScan_bins_PbPb3;

	}
	// may need to consider shape contribution from nonprompt later, no significant difference for PNP, phi f0 channel

	TFile *f_data=TFile::Open(dataName.Data());
	TFile *f_mc_Prompt=TFile::Open(mcName_Prompt.Data());
	TFile *f_mc_NonPrompt=TFile::Open(mcName_NonPrompt.Data());

	TFile *f_out=TFile::Open(outName.Data(),"RECREATE");
	f_out->cd();

	TH1F *h_RawFitYield=new TH1F("h_RawFitYield","h_RawFitYield",1,0,1); h_RawFitYield->Sumw2();
	TH1F *h_RawBinFitYield=new TH1F("h_RawBinFitYield","h_RawBinFitYield",1,0,1); h_RawBinFitYield->Sumw2();
	TH1F *h_RawRooFitYield=new TH1F("h_RawRooFitYield","h_RawRooFitYield",1,0,1); h_RawRooFitYield->Sumw2();
  // TH1F *h_RawFitYield_Loose=new TH1F("h_RawFitYield_Loose","h_RawFitYield_Loose",nbin_pt,bins_pt);

	TH1F *h_YieldDcabin=new TH1F("h_YieldDcabin","h_YieldDcabin",nbin_dca,bins_dca); h_YieldDcabin->Sumw2();

	TH1F *h_DsMassDsMC_Prompt;
	TH1F *h_DsMassDsMC_NonPrompt;
	TH1F *h_DsMassData;
	TTree *t_DsMassData;
	TTree *t_DsMassDsMC_Prompt;
	TTree *t_DsMassDsMC_NonPrompt;

	TH1F *h_DsMassDsMCLoose;
	TH1F *h_DsMassDataLoose;
	TTree *t_DsMassDsMCLoose;
	TTree *t_DsMassDataLoose;


	cout<<"check 1"<<endl;
	// cout<<"mc histogram = "<<Form("h_DsMass_pt%.0fto%.0f",DptLow,DptHigh)<<endl;

	h_DsMassDsMC_Prompt=(TH1F*)f_mc_Prompt->Get(Form("h_DsMass_pt%.0fto%.0f",DptLow,DptHigh));
	h_DsMassDsMC_NonPrompt=(TH1F*)f_mc_NonPrompt->Get(Form("h_DsMass_pt%.0fto%.0f",DptLow,DptHigh));
	// h_DsMassDsMC_Prompt->Draw();
	cout<<"check 2"<<endl;
	t_DsMassDsMC_Prompt=(TTree*)f_mc_Prompt->Get(Form("t_DsMass_pt%.0fto%.0f",DptLow,DptHigh));
	t_DsMassDsMC_NonPrompt=(TTree*)f_mc_NonPrompt->Get(Form("t_DsMass_pt%.0fto%.0f",DptLow,DptHigh));
	cout<<"check 3"<<endl;
	h_DsMassData=(TH1F*)f_data->Get(Form("h_DsMassData_pt%.0fto%.0f",DptLow,DptHigh));
	h_DsMassData->Draw();
	cout<<"check 4"<<endl;
	// t_DsMassData=(TTree*)f_data->Get(Form("t_DsMassData_pt%.0fto%.0f",DptLow,DptHigh));


	// private MC doesn't have enough statistics, temporary use NonPrompt for PbPb , Dpt<6 bins , change back when official MC is ready
	TTree *t_DsMassDsMC;
	// if(isPbPb==3 && DptHigh<=6)
	if(DptHigh<=6)
	{
		t_DsMassDsMC=t_DsMassDsMC_NonPrompt;
	}else{
    t_DsMassDsMC=t_DsMassDsMC_Prompt;
	}

	cout<<"check 5"<<endl;
	// h_DsMassDsMCLoose=(TH1F*)f_mc->Get(Form("h_DsMassLoose_pt%.0fto%.0f",DptLow,DptHigh));
	// t_DsMassDsMCLoose=(TTree*)f_mc->Get(Form("t_DsMassLoose_pt%.0fto%.0f",DptLow,DptHigh));
	// h_DsMassDataLoose=(TH1F*)f_data->Get(Form("h_DsMassDataLoose_pt%.0fto%.0f",DptLow,DptHigh));
	// t_DsMassDataLoose=(TTree*)f_data->Get(Form("t_DsMassDataLoose_pt%.0fto%.0f",DptLow,DptHigh));

	h_DsMassData->Rebin(2);
	h_DsMassData->Rebin(2);
	// h_DsMassDataLoose->Rebin(2);

	if (ibin_Dpt>=5) {
		// h_DsMassData->Rebin(2);
		// h_DsMassDataLoose->Rebin(2);
	}
/*
	if(isPbPb){
		h_DsMassData->Rebin(2);
	}
*/

/*

  int tes1=TF1FitBkg( h_DsMassData, ibin_Dpt, DptLow , DptHigh, str_PbPb);
	int test=bkgFitTest(h_DsMassData, ibin_Dpt, DptLow , DptHigh, str_PbPb);
	cout<<"str_PbPb = "<<str_PbPb.Data()<<endl;
	return 1;
*/


	TF1 *f1_fitOutput=binnedFit(h_DsMassDsMC_NonPrompt,h_DsMassData,h_RawBinFitYield,ibin_Dpt,DptLow,DptHigh,str_PbPb);
	// TF1 *f1_fitOutput_test=binnedFit(h_DsMassDsMC_Prompt,h_DsMassData,h_RawBinFitYield,ibin_Dpt,DptLow,DptHigh,str_PbPb);
	cout<<"f1_fitOutput->GetParameter(0) ="<<f1_fitOutput->GetParameter(0)<<endl;


/*
// RooFitResult *unbinnedFit(TTree* t_Data,TTree* t_MC, TH1F* h_DsMassDsMC, TH1F* h_DsMassData, TH1F* h_RawRooFitYield, TH1F* h_YieldDcabin ,Int_t isPbPb=0, Int_t ibin_Dpt=1, Double_t DptLow=0 ,Double_t DptHigh=10 ,TString str_PbPb="pp", Int_t doBkgVariation=0 ,Int_t BkgFunction=0, Int_t doSigVariation=0,double FloatWidthVal=0 , double FloatWidthErr=0, double FloatWidth_DsMass=1.967, TString extraName="", Bool_t doDcaBinsFit= false , Bool_t doMCstudy= false, TString extraName2=""){
	int doBkgVariation=0; int BkgFunction=0; int doSigVariation=0; double FloatWidthVal=0; double FloatWidthErr=0; 
	double FloatWidth_DsMass=1.965; TString extraName=""; bool doDcaBinsFit=false; Bool_t doMCstudy=doMCstudySet; TString extraName2="";

	RooFitResult *Roofitresult = unbinnedFit(t_DsMassData, t_DsMassDsMC, h_DsMassDsMC_Prompt, h_DsMassData, h_RawRooFitYield, h_YieldDcabin, isPbPb , ibin_Dpt, DptLow, DptHigh, str_PbPb, doBkgVariation, BkgFunction, doSigVariation, FloatWidthVal, FloatWidthErr , FloatWidth_DsMass , extraName, doDcaBinsFit, doMCstudy);

	Roofitresult->Print("V");
	Roofitresult->SetName(Form("RoofitResult_pt%.0fto%.0f",DptLow,DptHigh));


	cout<<"FloatWidthVal_Norm = "<<FloatWidthVal_Norm<<endl;
	cout<<"FloatWidthErr_Norm = "<<FloatWidthErr_Norm<<endl;

	// do FiFun Variation

	if(doSigVariationSet){

	cout<<"\n\n ### doSigVariation start "<<endl;

	TH1F *h_RawRooFitYield_sigVarP=new TH1F("h_RawRooFitYield_sigVarP","h_RawRooFitYield_sigVarP",1,0,1); h_RawRooFitYield_sigVarP->Sumw2();
	TH1F *h_RawRooFitYield_sigVarN=new TH1F("h_RawRooFitYield_sigVarN","h_RawRooFitYield_sigVarN",1,0,1); h_RawRooFitYield_sigVarN->Sumw2();
	TH1F *h_RawRooFitYield_sigVar=new TH1F("h_RawRooFitYield_sigVar","h_RawRooFitYield_sigVar",1,0,1); h_RawRooFitYield_sigVarN->Sumw2();

	doBkgVariation=0; BkgFunction=0; doSigVariation=1; FloatWidthVal=FloatWidthVal_Norm; FloatWidthErr=FloatWidthErr_Norm; 
	FloatWidth_DsMass=FloatWidth_DsMass_Norm;	extraName="SigWidthP_"; doDcaBinsFit=false; doMCstudy=false;
	RooFitResult *Roofitresult_sigVarP = unbinnedFit(t_DsMassData, t_DsMassDsMC, h_DsMassDsMC_Prompt, h_DsMassData, h_RawRooFitYield_sigVarP, h_YieldDcabin, isPbPb , ibin_Dpt, DptLow, DptHigh, str_PbPb, doBkgVariation, BkgFunction, doSigVariation, FloatWidthVal, FloatWidthErr , FloatWidth_DsMass , extraName, doDcaBinsFit, doMCstudy);

	doBkgVariation=0; BkgFunction=0; doSigVariation=1; FloatWidthVal=FloatWidthVal_Norm; FloatWidthErr= -1*FloatWidthErr_Norm; 
	FloatWidth_DsMass=FloatWidth_DsMass_Norm; extraName="SigWidthN_"; doDcaBinsFit=false; doMCstudy=false;
	RooFitResult *Roofitresult_sigVarN = unbinnedFit(t_DsMassData, t_DsMassDsMC, h_DsMassDsMC_Prompt, h_DsMassData, h_RawRooFitYield_sigVarN, h_YieldDcabin, isPbPb , ibin_Dpt, DptLow, DptHigh, str_PbPb, doBkgVariation, BkgFunction, doSigVariation, FloatWidthVal, FloatWidthErr , FloatWidth_DsMass ,extraName, doDcaBinsFit, doMCstudy);

	doBkgVariation=0; BkgFunction=0; doSigVariation=2; FloatWidthVal=FloatWidthVal_Norm; FloatWidthErr= -1*FloatWidthErr_Norm; 
	FloatWidth_DsMass=FloatWidth_DsMass_Norm; extraName="SigGaus1_"; doDcaBinsFit=false; doMCstudy=false;
	RooFitResult *Roofitresult_sigVar = unbinnedFit(t_DsMassData, t_DsMassDsMC, h_DsMassDsMC_Prompt, h_DsMassData, h_RawRooFitYield_sigVar, h_YieldDcabin, isPbPb , ibin_Dpt, DptLow, DptHigh, str_PbPb, doBkgVariation, BkgFunction, doSigVariation, FloatWidthVal, FloatWidthErr , FloatWidth_DsMass ,extraName, doDcaBinsFit, doMCstudy);


	f_out->cd();
	h_RawRooFitYield_sigVarP->Write("",TObject::kOverwrite);
	h_RawRooFitYield_sigVarN->Write("",TObject::kOverwrite);
	h_RawRooFitYield_sigVar->Write("",TObject::kOverwrite);


	cout<<" ### doSigVariation done \n "<<endl;

}

	if(doBkgVariationSet){

	cout<<"\n\n ### doBkgVariation start "<<endl;

	TH1F *h_RawRooFitYield_bkgVar=new TH1F("h_RawRooFitYield_bkgVar","h_RawRooFitYield_bkgVar",1,0,1); h_RawRooFitYield_bkgVar->Sumw2();
	doBkgVariation=1; BkgFunction=1; doSigVariation=0; FloatWidthVal=0; FloatWidthErr=0; 
	FloatWidth_DsMass=FloatWidth_DsMass_Norm; extraName=""; doDcaBinsFit=false; doMCstudy=false;
	RooFitResult *Roofitresult_bkgVar = unbinnedFit(t_DsMassData, t_DsMassDsMC, h_DsMassDsMC_Prompt, h_DsMassData, h_RawRooFitYield_bkgVar, h_YieldDcabin, isPbPb , ibin_Dpt, DptLow, DptHigh, str_PbPb, doBkgVariation, BkgFunction, doSigVariation, FloatWidthVal, FloatWidthErr , FloatWidth_DsMass ,extraName, doDcaBinsFit, doMCstudy);

	f_out->cd();
	h_RawRooFitYield_bkgVar->Write("",TObject::kOverwrite);
	cout<<" ### doBkgVariation done \n "<<endl;

}


	////-- cut scan fit --////

	if(doCutScanSet){

	//// DalphaMaxScan ////

  cout<<"\n\n doCutScan start "<<endl;
  cout<<"\n\n do DalphaMaxScan start "<<endl;


  TH1F *h_DsMassData_DalphaMaxScan[nbin_DalphaMaxScan];
  TTree *t_DsMassData_DalphaMaxScan[nbin_DalphaMaxScan];
  TH1F *h_DsMassDsMC_Prompt_DalphaMaxScan[nbin_DalphaMaxScan];
  TTree *t_DsMassDsMC_Prompt_DalphaMaxScan[nbin_DalphaMaxScan];
  TH1F *h_DsMassDsMC_NonPrompt_DalphaMaxScan[nbin_DalphaMaxScan];
  TTree *t_DsMassDsMC_NonPrompt_DalphaMaxScan[nbin_DalphaMaxScan];

	TH1F *h_RawRooFitYield_DalphaMaxScan[nbin_DalphaMaxScan];
	TH1F *h_RawBinFitYield_DalphaMaxScan[nbin_DalphaMaxScan];
	TTree *t_DsMassDsMC_DalphaMaxScan[nbin_DalphaMaxScan];

  for(int ibin = 0; ibin<nbin_DalphaMaxScan; ibin++){

      h_DsMassDsMC_Prompt_DalphaMaxScan[ibin]=(TH1F*)f_mc_Prompt->Get(Form("h_DsMass_pt%.0fto%.0f_Dalpha_%i",DptLow,DptHigh,ibin));
			t_DsMassDsMC_Prompt_DalphaMaxScan[ibin]=(TTree*)f_mc_Prompt->Get(Form("t_DsMass_pt%.0fto%.0f_Dalpha_%i",DptLow,DptHigh,ibin));

      h_DsMassDsMC_NonPrompt_DalphaMaxScan[ibin]=(TH1F*)f_mc_NonPrompt->Get(Form("h_DsMass_pt%.0fto%.0f_Dalpha_%i",DptLow,DptHigh,ibin));
			t_DsMassDsMC_NonPrompt_DalphaMaxScan[ibin]=(TTree*)f_mc_NonPrompt->Get(Form("t_DsMass_pt%.0fto%.0f_Dalpha_%i",DptLow,DptHigh,ibin));

			h_DsMassData_DalphaMaxScan[ibin]=(TH1F*)f_data->Get(Form("h_DsMassData_pt%.0fto%.0f_Dalpha%.0f",DptLow,DptHigh,DalphaMaxScan_bins[ibin]*100));
			t_DsMassData_DalphaMaxScan[ibin]=(TTree*)f_data->Get(Form("t_DsMassData_pt%.0fto%.0f_Dalpha%.0f",DptLow,DptHigh,DalphaMaxScan_bins[ibin]*100));

			f_out->cd();
			h_RawRooFitYield_DalphaMaxScan[ibin]=new TH1F(Form("h_RawRooFitYield_DalphaMaxScan%.0f",DalphaMaxScan_bins[ibin]*100), Form("h_RawRooFitYield_DalphaMaxScan%.0f",DalphaMaxScan_bins[ibin]*100),1,0,1);
			h_RawBinFitYield_DalphaMaxScan[ibin]=new TH1F(Form("h_RawBinFitYield_DalphaMaxScan%.0f",DalphaMaxScan_bins[ibin]*100), Form("h_RawBinFitYield_DalphaMaxScan%.0f",DalphaMaxScan_bins[ibin]*100),1,0,1);
	// private MC doesn't have enough statistics, temporary use NonPrompt for PbPb , Dpt<6 bins , change back when official MC is ready
	// if(isPbPb==3 && DptHigh<=6)
	if(DptHigh<=6)
	{
		t_DsMassDsMC_DalphaMaxScan[ibin]=t_DsMassDsMC_NonPrompt_DalphaMaxScan[ibin];
	}else{
    t_DsMassDsMC_DalphaMaxScan[ibin]=t_DsMassDsMC_Prompt_DalphaMaxScan[ibin];
	}

	h_DsMassData_DalphaMaxScan[ibin]->Rebin(2);
	if (DptLow>=10) {
		h_DsMassData_DalphaMaxScan[ibin]->Rebin(2);
	}

// TF1 *binnedFit(TH1F* h_DsMassDsMC, TH1F* h_DsMassData, TH1F* h_RawBinFitYield, Int_t ibin_Dpt, Double_t DptLow=0 ,Double_t DptHigh=10 ,TString str_PbPb="pp", Int_t doBkgVariation=0 ,Int_t BkgFunction=0, Int_t doSigVariation=0, double SigWidthErr=0 ,TString extraName=""){

// RooFitResult *unbinnedFit(TTree* t_Data,TTree* t_MC, TH1F* h_DsMassDsMC, TH1F* h_DsMassData, TH1F* h_RawRooFitYield, TH1F* h_YieldDcabin ,Int_t isPbPb=0, Int_t ibin_Dpt=1, Double_t DptLow=0 ,Double_t DptHigh=10 ,TString str_PbPb="pp", Int_t doBkgVariation=0 ,Int_t BkgFunction=0, Int_t doSigVariation=0,double FloatWidthVal=0 , double FloatWidthErr=0, double FloatWidth_DsMass=1.967, TString extraName="", Bool_t doDcaBinsFit= false , Bool_t doMCstudy= false,TString extraName2="")

	doBkgVariation=0; BkgFunction=0; doSigVariation=0; FloatWidthVal=0; FloatWidthErr=0; 
	FloatWidth_DsMass=1.965;  extraName="DalphaScan"; doDcaBinsFit=false; doMCstudy=0; extraName2=Form("DalphaMaxScan%0.f",DalphaMaxScan_bins[ibin]*100);

	binnedFit(h_DsMassDsMC_Prompt_DalphaMaxScan[ibin],h_DsMassData_DalphaMaxScan[ibin],h_RawBinFitYield_DalphaMaxScan[ibin],ibin_Dpt,DptLow,DptHigh,str_PbPb, 0,0,0,0,extraName2); //// saving names has not been set.
	unbinnedFit(t_DsMassData_DalphaMaxScan[ibin], t_DsMassDsMC_DalphaMaxScan[ibin], h_DsMassDsMC_Prompt, h_DsMassData_DalphaMaxScan[ibin], h_RawRooFitYield_DalphaMaxScan[ibin], h_YieldDcabin, isPbPb , ibin_Dpt, DptLow, DptHigh, str_PbPb, doBkgVariation, BkgFunction, doSigVariation, FloatWidthVal, FloatWidthErr , FloatWidth_DsMass , extraName, doDcaBinsFit, doMCstudy, extraName2);
	// unbinnedFit(t_DsMassData_DalphaMaxScan[ibin], t_DsMassDsMC_DalphaMaxScan[ibin], h_DsMassDsMC_Prompt, h_DsMassData, h_RawRooFitYield_DalphaMaxScan[ibin], h_YieldDcabin, isPbPb , ibin_Dpt, DptLow, DptHigh, str_PbPb, doBkgVariation, BkgFunction, doSigVariation, FloatWidthVal, FloatWidthErr , FloatWidth_DsMass , extraName, doDcaBinsFit, doMCstudy, extraName2);

	/// temp use norm fit as init parameters
	// Roofitresult->Print("V");

	f_out->cd();
	h_RawRooFitYield_DalphaMaxScan[ibin]->Write("",TObject::kOverwrite);

	}

	//// Dchi2clMinScan ////


  cout<<"\n\n do Dchi2clMinScan start "<<endl;

  TH1F *h_DsMassData_Dchi2clMinScan[nbin_Dchi2clMinScan];
  TTree *t_DsMassData_Dchi2clMinScan[nbin_Dchi2clMinScan];
  TH1F *h_DsMassDsMC_Prompt_Dchi2clMinScan[nbin_Dchi2clMinScan];
  TTree *t_DsMassDsMC_Prompt_Dchi2clMinScan[nbin_Dchi2clMinScan];
  TH1F *h_DsMassDsMC_NonPrompt_Dchi2clMinScan[nbin_Dchi2clMinScan];
  TTree *t_DsMassDsMC_NonPrompt_Dchi2clMinScan[nbin_Dchi2clMinScan];

	TH1F *h_RawRooFitYield_Dchi2clMinScan[nbin_Dchi2clMinScan];
	TH1F *h_RawBinFitYield_Dchi2clMinScan[nbin_Dchi2clMinScan];
	TTree *t_DsMassDsMC_Dchi2clMinScan[nbin_Dchi2clMinScan];

  for(int ibin = 0; ibin<nbin_Dchi2clMinScan; ibin++){

      h_DsMassDsMC_Prompt_Dchi2clMinScan[ibin]=(TH1F*)f_mc_Prompt->Get(Form("h_DsMass_pt%.0fto%.0f_Dchi2cl_%i",DptLow,DptHigh,ibin));
			t_DsMassDsMC_Prompt_Dchi2clMinScan[ibin]=(TTree*)f_mc_Prompt->Get(Form("t_DsMass_pt%.0fto%.0f_Dchi2cl_%i",DptLow,DptHigh,ibin));

      h_DsMassDsMC_NonPrompt_Dchi2clMinScan[ibin]=(TH1F*)f_mc_NonPrompt->Get(Form("h_DsMass_pt%.0fto%.0f_Dchi2cl_%i",DptLow,DptHigh,ibin));
			t_DsMassDsMC_NonPrompt_Dchi2clMinScan[ibin]=(TTree*)f_mc_NonPrompt->Get(Form("t_DsMass_pt%.0fto%.0f_Dchi2cl_%i",DptLow,DptHigh,ibin));

			h_DsMassData_Dchi2clMinScan[ibin]=(TH1F*)f_data->Get(Form("h_DsMassData_pt%.0fto%.0f_Dchi2cl%.0f",DptLow,DptHigh,Dchi2clMinScan_bins[ibin]*100));
			t_DsMassData_Dchi2clMinScan[ibin]=(TTree*)f_data->Get(Form("t_DsMassData_pt%.0fto%.0f_Dchi2cl%.0f",DptLow,DptHigh,Dchi2clMinScan_bins[ibin]*100));

			f_out->cd();
			h_RawRooFitYield_Dchi2clMinScan[ibin]=new TH1F(Form("h_RawRooFitYield_Dchi2clMinScan%.0f",Dchi2clMinScan_bins[ibin]*100), Form("h_RawRooFitYield_Dchi2clMinScan%.0f",Dchi2clMinScan_bins[ibin]*100),1,0,1);
			h_RawBinFitYield_Dchi2clMinScan[ibin]=new TH1F(Form("h_RawBinFitYield_Dchi2clMinScan%.0f",Dchi2clMinScan_bins[ibin]*100), Form("h_RawBinFitYield_Dchi2clMinScan%.0f",Dchi2clMinScan_bins[ibin]*100),1,0,1);
	// private MC doesn't have enough statistics, temporary use NonPrompt for PbPb , Dpt<6 bins , change back when official MC is ready
	// if(isPbPb==3 && DptHigh<=6)
	if(DptHigh<=6)
	{
		t_DsMassDsMC_Dchi2clMinScan[ibin]=t_DsMassDsMC_NonPrompt_Dchi2clMinScan[ibin];
	}else{
    t_DsMassDsMC_Dchi2clMinScan[ibin]=t_DsMassDsMC_Prompt_Dchi2clMinScan[ibin];
	}

	h_DsMassData_Dchi2clMinScan[ibin]->Rebin(2);
	if (DptLow>=10) {
		h_DsMassData_Dchi2clMinScan[ibin]->Rebin(2);
	}

	// binnedFit(h_DsMassDsMC_Prompt,h_DsMassData,h_RawBinFitYield,ibin_Dpt,DptLow,DptHigh,str_PbPb); //// saving names has not been set.
	// TF1 *f1_fitOutput_test=binnedFit(h_DsMassDsMC_Prompt,h_DsMassData,h_RawBinFitYield,ibin_Dpt,DptLow,DptHigh,str_PbPb);
// RooFitResult *unbinnedFit(TTree* t_Data,TTree* t_MC, TH1F* h_DsMassDsMC, TH1F* h_DsMassData, TH1F* h_RawRooFitYield, TH1F* h_YieldDcabin ,Int_t isPbPb=0, Int_t ibin_Dpt=1, Double_t DptLow=0 ,Double_t DptHigh=10 ,TString str_PbPb="pp", Int_t doBkgVariation=0 ,Int_t BkgFunction=0, Int_t doSigVariation=0,double FloatWidthVal=0 , double FloatWidthErr=0, double FloatWidth_DsMass=1.967, TString extraName="", Bool_t doDcaBinsFit= false , Bool_t doMCstudy= false,TString extraName2="")
	doBkgVariation=0; BkgFunction=0; doSigVariation=0; FloatWidthVal=0; FloatWidthErr=0; 
	FloatWidth_DsMass=1.965;  extraName="Dchi2clScan"; doDcaBinsFit=false; doMCstudy=0; extraName2=Form("Dchi2clMinScan%0.f",Dchi2clMinScan_bins[ibin]*100);

	binnedFit(h_DsMassDsMC_Prompt_Dchi2clMinScan[ibin],h_DsMassData_Dchi2clMinScan[ibin],h_RawBinFitYield_Dchi2clMinScan[ibin],ibin_Dpt,DptLow,DptHigh,str_PbPb, 0,0,0,0,extraName2); //// saving names has not been set.
	unbinnedFit(t_DsMassData_Dchi2clMinScan[ibin], t_DsMassDsMC_Dchi2clMinScan[ibin], h_DsMassDsMC_Prompt, h_DsMassData_Dchi2clMinScan[ibin], h_RawRooFitYield_Dchi2clMinScan[ibin], h_YieldDcabin, isPbPb , ibin_Dpt, DptLow, DptHigh, str_PbPb, doBkgVariation, BkgFunction, doSigVariation, FloatWidthVal, FloatWidthErr , FloatWidth_DsMass , extraName, doDcaBinsFit, doMCstudy, extraName2);

	/// temp use norm fit as init parameters
	// Roofitresult->Print("V");

	f_out->cd();
	h_RawRooFitYield_Dchi2clMinScan[ibin]->Write("",TObject::kOverwrite);

	}


	//// DdlsMinScan ////

  cout<<"\n\n do DdlsMinScan start "<<endl;

  TH1F *h_DsMassData_DdlsMinScan[nbin_DdlsMinScan];
  TTree *t_DsMassData_DdlsMinScan[nbin_DdlsMinScan];
  TH1F *h_DsMassDsMC_Prompt_DdlsMinScan[nbin_DdlsMinScan];
  TTree *t_DsMassDsMC_Prompt_DdlsMinScan[nbin_DdlsMinScan];
  TH1F *h_DsMassDsMC_NonPrompt_DdlsMinScan[nbin_DdlsMinScan];
  TTree *t_DsMassDsMC_NonPrompt_DdlsMinScan[nbin_DdlsMinScan];

	TH1F *h_RawRooFitYield_DdlsMinScan[nbin_DdlsMinScan];
	TH1F *h_RawBinFitYield_DdlsMinScan[nbin_DdlsMinScan];
	TTree *t_DsMassDsMC_DdlsMinScan[nbin_DdlsMinScan];

  for(int ibin = 0; ibin<nbin_DdlsMinScan; ibin++){

      h_DsMassDsMC_Prompt_DdlsMinScan[ibin]=(TH1F*)f_mc_Prompt->Get(Form("h_DsMass_pt%.0fto%.0f_Ddls_%i",DptLow,DptHigh,ibin));
			t_DsMassDsMC_Prompt_DdlsMinScan[ibin]=(TTree*)f_mc_Prompt->Get(Form("t_DsMass_pt%.0fto%.0f_Ddls_%i",DptLow,DptHigh,ibin));

      h_DsMassDsMC_NonPrompt_DdlsMinScan[ibin]=(TH1F*)f_mc_NonPrompt->Get(Form("h_DsMass_pt%.0fto%.0f_Ddls_%i",DptLow,DptHigh,ibin));
			t_DsMassDsMC_NonPrompt_DdlsMinScan[ibin]=(TTree*)f_mc_NonPrompt->Get(Form("t_DsMass_pt%.0fto%.0f_Ddls_%i",DptLow,DptHigh,ibin));

			h_DsMassData_DdlsMinScan[ibin]=(TH1F*)f_data->Get(Form("h_DsMassData_pt%.0fto%.0f_Ddls%.0f",DptLow,DptHigh,DdlsMinScan_bins[ibin]*10));
			t_DsMassData_DdlsMinScan[ibin]=(TTree*)f_data->Get(Form("t_DsMassData_pt%.0fto%.0f_Ddls%.0f",DptLow,DptHigh,DdlsMinScan_bins[ibin]*10));

			f_out->cd();
			h_RawRooFitYield_DdlsMinScan[ibin]=new TH1F(Form("h_RawRooFitYield_DdlsMinScan%.0f",DdlsMinScan_bins[ibin]*100), Form("h_RawRooFitYield_DdlsMinScan%.0f",DdlsMinScan_bins[ibin]*100),1,0,1);
			h_RawBinFitYield_DdlsMinScan[ibin]=new TH1F(Form("h_RawBinFitYield_DdlsMinScan%.0f",DdlsMinScan_bins[ibin]*100), Form("h_RawBinFitYield_DdlsMinScan%.0f",DdlsMinScan_bins[ibin]*100),1,0,1);
	// private MC doesn't have enough statistics, temporary use NonPrompt for PbPb , Dpt<6 bins , change back when official MC is ready
	// if(isPbPb==3 && DptHigh<=6)
	if( DptHigh<=6)
	{
		t_DsMassDsMC_DdlsMinScan[ibin]=t_DsMassDsMC_NonPrompt_DdlsMinScan[ibin];
	}else{
    t_DsMassDsMC_DdlsMinScan[ibin]=t_DsMassDsMC_Prompt_DdlsMinScan[ibin];
	}

	h_DsMassData_DdlsMinScan[ibin]->Rebin(2);
	if (DptLow>=10) {
		h_DsMassData_DdlsMinScan[ibin]->Rebin(2);
	}

	// binnedFit(h_DsMassDsMC_Prompt,h_DsMassData,h_RawBinFitYield,ibin_Dpt,DptLow,DptHigh,str_PbPb); //// saving names has not been set.
	// TF1 *f1_fitOutput_test=binnedFit(h_DsMassDsMC_Prompt,h_DsMassData,h_RawBinFitYield,ibin_Dpt,DptLow,DptHigh,str_PbPb);
// RooFitResult *unbinnedFit(TTree* t_Data,TTree* t_MC, TH1F* h_DsMassDsMC, TH1F* h_DsMassData, TH1F* h_RawRooFitYield, TH1F* h_YieldDcabin ,Int_t isPbPb=0, Int_t ibin_Dpt=1, Double_t DptLow=0 ,Double_t DptHigh=10 ,TString str_PbPb="pp", Int_t doBkgVariation=0 ,Int_t BkgFunction=0, Int_t doSigVariation=0,double FloatWidthVal=0 , double FloatWidthErr=0, double FloatWidth_DsMass=1.967, TString extraName="", Bool_t doDcaBinsFit= false , Bool_t doMCstudy= false,TString extraName2="")
	doBkgVariation=0; BkgFunction=0; doSigVariation=0; FloatWidthVal=0; FloatWidthErr=0; 
	FloatWidth_DsMass=1.965;  extraName="DdlsScan"; doDcaBinsFit=false; doMCstudy=0; extraName2=Form("DdlsMinScan%0.f",DdlsMinScan_bins[ibin]*10);

	binnedFit(h_DsMassDsMC_Prompt_DdlsMinScan[ibin],h_DsMassData_DdlsMinScan[ibin],h_RawBinFitYield_DdlsMinScan[ibin],ibin_Dpt,DptLow,DptHigh,str_PbPb, 0,0,0,0,extraName2); //// saving names has not been set.
	unbinnedFit(t_DsMassData_DdlsMinScan[ibin], t_DsMassDsMC_DdlsMinScan[ibin], h_DsMassDsMC_Prompt, h_DsMassData_DdlsMinScan[ibin], h_RawRooFitYield_DdlsMinScan[ibin], h_YieldDcabin, isPbPb , ibin_Dpt, DptLow, DptHigh, str_PbPb, doBkgVariation, BkgFunction, doSigVariation, FloatWidthVal, FloatWidthErr , FloatWidth_DsMass , extraName, doDcaBinsFit, doMCstudy, extraName2);

	/// temp use norm fit as init parameters
	// Roofitresult->Print("V");

	f_out->cd();
	h_RawRooFitYield_DdlsMinScan[ibin]->Write("",TObject::kOverwrite);

	}


	// PhiMass Scan

	//// PhiMassScan ////

  cout<<"\n\n do PhiMassScan start "<<endl;

  TH1F *h_DsMassData_PhiMassScan[nbin_PhiMassScan];
  TTree *t_DsMassData_PhiMassScan[nbin_PhiMassScan];
  TH1F *h_DsMassDsMC_Prompt_PhiMassScan[nbin_PhiMassScan];
  TTree *t_DsMassDsMC_Prompt_PhiMassScan[nbin_PhiMassScan];
  TH1F *h_DsMassDsMC_NonPrompt_PhiMassScan[nbin_PhiMassScan];
  TTree *t_DsMassDsMC_NonPrompt_PhiMassScan[nbin_PhiMassScan];

	TH1F *h_RawRooFitYield_PhiMassScan[nbin_PhiMassScan];
	TH1F *h_RawBinFitYield_PhiMassScan[nbin_PhiMassScan];
	TTree *t_DsMassDsMC_PhiMassScan[nbin_PhiMassScan];

  for(int ibin = 0; ibin<nbin_PhiMassScan; ibin++){

      h_DsMassDsMC_Prompt_PhiMassScan[ibin]=(TH1F*)f_mc_Prompt->Get(Form("h_DsMass_pt%.0fto%.0f_PhiMass_%i",DptLow,DptHigh,ibin));
			t_DsMassDsMC_Prompt_PhiMassScan[ibin]=(TTree*)f_mc_Prompt->Get(Form("t_DsMass_pt%.0fto%.0f_PhiMass_%i",DptLow,DptHigh,ibin));

      h_DsMassDsMC_NonPrompt_PhiMassScan[ibin]=(TH1F*)f_mc_NonPrompt->Get(Form("h_DsMass_pt%.0fto%.0f_PhiMass_%i",DptLow,DptHigh,ibin));
			t_DsMassDsMC_NonPrompt_PhiMassScan[ibin]=(TTree*)f_mc_NonPrompt->Get(Form("t_DsMass_pt%.0fto%.0f_PhiMass_%i",DptLow,DptHigh,ibin));

			h_DsMassData_PhiMassScan[ibin]=(TH1F*)f_data->Get(Form("h_DsMassData_pt%.0fto%.0f_PhiMass_%i",DptLow,DptHigh,ibin));
			t_DsMassData_PhiMassScan[ibin]=(TTree*)f_data->Get(Form("t_DsMassData_pt%.0fto%.0f_PhiMass_%i",DptLow,DptHigh,ibin));

			f_out->cd();
			h_RawRooFitYield_PhiMassScan[ibin]=new TH1F(Form("h_RawRooFitYield_PhiMassScan_%i",ibin), Form("h_RawRooFitYield_PhiMassScan_%i",ibin),1,0,1);
			h_RawBinFitYield_PhiMassScan[ibin]=new TH1F(Form("h_RawBinFitYield_PhiMassScan_%i",ibin), Form("h_RawBinFitYield_PhiMassScan_%i",ibin),1,0,1);
	// private MC doesn't have enough statistics, temporary use NonPrompt for PbPb , Dpt<6 bins , change back when official MC is ready
	// if(isPbPb==3 && DptHigh<=6)
	if(DptHigh<=6)
	{
		t_DsMassDsMC_PhiMassScan[ibin]=t_DsMassDsMC_NonPrompt_PhiMassScan[ibin];
	}else{
    t_DsMassDsMC_PhiMassScan[ibin]=t_DsMassDsMC_Prompt_PhiMassScan[ibin];
	}

	h_DsMassData_PhiMassScan[ibin]->Rebin(2);
	if (DptLow>=10) {
		h_DsMassData_PhiMassScan[ibin]->Rebin(2);
	}

	// binnedFit(h_DsMassDsMC_Prompt,h_DsMassData,h_RawBinFitYield,ibin_Dpt,DptLow,DptHigh,str_PbPb); //// saving names has not been set.
	// TF1 *f1_fitOutput_test=binnedFit(h_DsMassDsMC_Prompt,h_DsMassData,h_RawBinFitYield,ibin_Dpt,DptLow,DptHigh,str_PbPb);
// RooFitResult *unbinnedFit(TTree* t_Data,TTree* t_MC, TH1F* h_DsMassDsMC, TH1F* h_DsMassData, TH1F* h_RawRooFitYield, TH1F* h_YieldDcabin ,Int_t isPbPb=0, Int_t ibin_Dpt=1, Double_t DptLow=0 ,Double_t DptHigh=10 ,TString str_PbPb="pp", Int_t doBkgVariation=0 ,Int_t BkgFunction=0, Int_t doSigVariation=0,double FloatWidthVal=0 , double FloatWidthErr=0, double FloatWidth_DsMass=1.967, TString extraName="", Bool_t doDcaBinsFit= false , Bool_t doMCstudy= false,TString extraName2="")
	doBkgVariation=0; BkgFunction=0; doSigVariation=0; FloatWidthVal=0; FloatWidthErr=0; 
	FloatWidth_DsMass=1.965;  extraName="PhiMassScan"; doDcaBinsFit=false; doMCstudy=0; extraName2=Form("PhiMassScan_%i",ibin);

	binnedFit(h_DsMassDsMC_Prompt_PhiMassScan[ibin],h_DsMassData_PhiMassScan[ibin],h_RawBinFitYield_PhiMassScan[ibin],ibin_Dpt,DptLow,DptHigh,str_PbPb, 0,0,0,0,extraName2); //// saving names has not been set.
	unbinnedFit(t_DsMassData_PhiMassScan[ibin], t_DsMassDsMC_PhiMassScan[ibin], h_DsMassDsMC_Prompt, h_DsMassData_PhiMassScan[ibin], h_RawRooFitYield_PhiMassScan[ibin], h_YieldDcabin, isPbPb , ibin_Dpt, DptLow, DptHigh, str_PbPb, doBkgVariation, BkgFunction, doSigVariation, FloatWidthVal, FloatWidthErr , FloatWidth_DsMass , extraName, doDcaBinsFit, doMCstudy, extraName2);

	/// temp use norm fit as init parameters
	// Roofitresult->Print("V");

	f_out->cd();
	h_RawRooFitYield_PhiMassScan[ibin]->Write("",TObject::kOverwrite);

	}





	cout<<" ### end , cut scan fit done"<<endl;

	///-- end cut scan fit --////
} // end if cut scan fit
*/
	f_out->cd();
	h_DsMassData->Write("",TObject::kOverwrite);
	f1_fitOutput->Write();// this is returned from fit function
	// myfunc->Write("myfunc",TObject::kOverwrite); // this is from function stored in histogram
	h_RawBinFitYield->Write("",TObject::kOverwrite);
	h_RawRooFitYield->Write("",TObject::kOverwrite);
	h_YieldDcabin->Write("",TObject::kOverwrite);
	// h_RawFitYield_Loose->Write("",TObject::kOverwrite);
	// Roofitresult->Write("",TObject::kOverwrite);

	f_out->Close();

	cout<<"ibin_Dpt = "<<ibin_Dpt<<" , DptLow = "<<DptLow<<" , DptHigh = "<<DptHigh<<endl;
	// cout<<"RooFit edm = "<<Roofitresult->edm()<<endl;

//	c_fitData[ibin_Dpt]->SaveAs(Form("plots/DsFitData_pt%.0fto%.0f_pp.pdf",DptLow,DptHigh));

/*  // test plot outside fit function
	TCanvas *ctest=new TCanvas("ctest","ctest",800,800);
	ctest->cd();
	h_DsMassData->Draw("same");
//	f1_fitOutput->Draw("same"); // still empty
	TF1 *myfunc=h_DsMassData->GetFunction("f1_DsMix");
	// myfunc->SetName(Form("f1_DsMix_pt%.0fto%.0f",DptLow,DptHigh)); // setname not work
	myfunc->Draw("same");
	ctest->SaveAs("testfitplot.pdf"); // this works
*/

	return 0;
}
