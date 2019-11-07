#include "MkkFit_para.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"

#include "TFitter.h"
#include "TFitResult.h"

TF1* fitDsMassMC(TH1D* h_DsMassDsMC);
TF1* fitDsMass(TH1D *hData, TF1 *f1MC,double &FitYield,double &FitYieldErr,double &BkgSubYield,double &BkgSubYieldErr,  TString s_ExtraName="",int fixSigShape=0, TF1 *f1_default=NULL);
//	gSystem->Exec("mkdir -p plots/fit");


  double DsMassCandWidth=0.05;
  double DsMassCandFitMeanWidth=0.02;
  double shiftY=0.05;
  double shiftX=0.32;
  double oneshift=0.070;

	double DptLowGL=5;
	double DptHighGL=40;
	TString s_ppPbPb="pp";
	TString s_type="";


int FitMkkpi(int isPbPb=0, double DptLow=6, double DptHigh=40){

	DptLowGL=DptLow;
	DptHighGL=DptHighGL;

	s_ppPbPb="pp";

	if(isPbPb==3)
	{
		s_ppPbPb="PbPb";
	}



  gSystem->Exec("mkdir -p plots/fit");

//	TString s_ppPbPb="pp";

	TString s_finName=Form("./his/%s_Dpt%.0fto%.0f.root",s_ppPbPb.Data(),DptLow,DptHigh);
	TString s_foutName=Form("./his/%s_Dpt%.0fto%.0f_mkkpiFit.root",s_ppPbPb.Data(),DptLow,DptHigh);

	s_type=Form("%s_Dpt%.0fto%.0f",s_ppPbPb.Data(),DptLow,DptHigh);

	TFile *fin=TFile::Open(s_finName.Data(),"READ");
	TFile *fout=TFile::Open(s_foutName.Data(),"RECREATE");

	TH1D *h_mkkpi_default_MCP_Phi=(TH1D*)fin->Get("h_mkkpi_default_MCP_Phi"); 
	TH1D *h_mkkpi_default_MCP_f0 =(TH1D*)fin->Get("h_mkkpi_default_MCP_f0");
	TH1D *h_mkk_MCP_Phi          =(TH1D*)fin->Get("h_mkk_MCP_Phi");
	TH1D *h_mkk_MCP_f0           =(TH1D*)fin->Get("h_mkk_MCP_f0");
	TH1D *h_mkk_MCP_Phi_fineBin  =(TH1D*)fin->Get("h_mkk_MCP_Phi_fineBin");
	TH1D *h_mkk_MCP_f0_fineBin   =(TH1D*)fin->Get("h_mkk_MCP_f0_fineBin");
	TH1D *h_mkkpi_default        =(TH1D*)fin->Get("h_mkkpi_default");

	TF1 *f1MC=fitDsMassMC(h_mkkpi_default_MCP_Phi);

	cout<<"f1MC->GetParameter(1) = "<<f1MC->GetParameter(1)<<endl;

	double FitYield;
	double FitYieldErr;
	double BkgSubYield;
	double BkgSubYieldErr;

	TF1 *f1Data=fitDsMass(h_mkkpi_default,f1MC, FitYield,FitYieldErr,BkgSubYield,BkgSubYieldErr);
	
	TH1D *h_mkkpi_mkkbins[nbin_mkk];
	TF1 *f1Data_mkkbins[nbin_mkk];

	TString extraName="";

	fout->cd();
	TH1D *h_NDs_mkkbin_Fit=new TH1D("h_NDs_mkkbin_Fit",";m_{KK}",nbin_mkk,bins_mkk); h_NDs_mkkbin_Fit->Sumw2();
	TH1D *h_NDs_mkkbin_BkgSub=new TH1D("h_NDs_mkkbin_BkgSub",";m_{KK}",nbin_mkk,bins_mkk); h_NDs_mkkbin_BkgSub->Sumw2();


	for(int i=0; i<nbin_mkk; i++)
	 // for(int i=0; i<1; i++)
	{
		h_mkkpi_mkkbins[i]=(TH1D*)fin->Get(Form("h_mkkpi_mkkbins_%i",i));
		extraName=Form("mkk%.0fto%.0f",bins_mkk[i]*1000,bins_mkk[i+1]*1000);
		f1Data_mkkbins[i]=fitDsMass(h_mkkpi_mkkbins[i],f1MC,FitYield,FitYieldErr,BkgSubYield,BkgSubYieldErr,extraName,1,f1Data);

		h_NDs_mkkbin_Fit->SetBinContent(i+1,FitYield);
		h_NDs_mkkbin_Fit->SetBinError(i+1,FitYieldErr);

		h_NDs_mkkbin_BkgSub->SetBinContent(i+1,BkgSubYield);
		h_NDs_mkkbin_BkgSub->SetBinError(i+1,BkgSubYieldErr);


	}

	// cout<<"FitYield = "<<FitYield<<endl;

		h_NDs_mkkbin_Fit->Write();
		h_NDs_mkkbin_BkgSub->Write();


	TLatex *tlatex=new TLatex();
	TCanvas *c_Nds_mkkbin_Fit=new TCanvas("c_Nds_mkkbin_Fit","c_Nds_mkkbin_Fit");
	c_Nds_mkkbin_Fit->cd();
	h_NDs_mkkbin_Fit->Draw();

	shiftY=0;
  tlatex->DrawLatexNDC(textposx+0.08,textposy+shiftY,Form("%s D_{S} %.0f < pt < %.0f ",s_ppPbPb.Data(),DptLowGL,DptHighGL)); shiftY-=oneshift;

	c_Nds_mkkbin_Fit->SaveAs(Form("./plots/fit/%s_Nds_mkkbin_Fit_Dpt%.0fto%.0f.png",s_type.Data(),DptLowGL,DptHighGL));

	TCanvas *c_Nds_mkkbin_BkgSub=new TCanvas("c_Nds_mkkbin_BkgSub","c_Nds_mkkbin_BkgSub");
	c_Nds_mkkbin_BkgSub->cd();
	h_NDs_mkkbin_BkgSub->Draw();

	shiftY=0;
  tlatex->DrawLatexNDC(textposx+0.08,textposy+shiftY,Form("%s D_{S} %.0f < pt < %.0f ",s_ppPbPb.Data(),DptLowGL,DptHighGL)); shiftY-=oneshift;

	c_Nds_mkkbin_BkgSub->SaveAs(Form("./plots/fit/%s_Nds_mkkbin_BkgSub_Dpt%.0fto%.0f.png",s_type.Data(),DptLowGL,DptHighGL));



	cout<<"s_ppPbPb = "<<s_ppPbPb.Data()<<" , isPbPb = "<<isPbPb<<endl;


	return 0;

}


TF1* fitDsMassMC(TH1D* h_DsMassDsMC){

  TF1 *f1_DsMC = new TF1("f1_DsMC","[0]*( [4]*TMath::Gaus(x,[1],[2]*(1+[5]))/( sqrt(2*TMath::Pi())*[2]*(1+[5]) ) + (1-[4])*TMath::Gaus(x,[1],[3]*(1+[5]))/(sqrt(2*TMath::Pi())*[3]*(1+[5])  ) )",DsMassFitLow,DsMassFitHigh); // kTRUE nomalize Gaus
  f1_DsMC->SetParameter(0, h_DsMassDsMC->GetBinWidth(1) * h_DsMassDsMC->Integral( h_DsMassDsMC->FindBin(DsMass-DsMassCandWidth) , h_DsMassDsMC->FindBin(DsMass+DsMassCandWidth) ) ) ;
  f1_DsMC->SetParameter(1,DsMass);   // DsMean
  f1_DsMC->SetParLimits(1,DsMass-DsMassCandFitMeanWidth,DsMass+DsMassCandFitMeanWidth);
  f1_DsMC->SetParameter(2,0.01); // defaut (wider) gaus1 sigma
  f1_DsMC->SetParLimits(2,0.0001,0.1);
  f1_DsMC->SetParameter(3,0.001); // default (narrow) gaus2 sigma
  f1_DsMC->SetParLimits(3,0.00001,0.015);
  f1_DsMC->SetParameter(4,0.5);  // fraction of Gaus1 / All
  f1_DsMC->SetParLimits(4,0,1);
  f1_DsMC->SetParameter(5,0); // for data /mc discrepancy
  f1_DsMC->SetParLimits(5,0,0);
  f1_DsMC->FixParameter(5,0);

  f1_DsMC->SetLineColor(kRed);

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

 double yieldDsMC= f1_DsMC->Integral(f1_DsMC->GetParameter(1)-DsMassCandWidth,f1_DsMC->GetParameter(1)+DsMassCandWidth)/ h_DsMassDsMC->GetBinWidth(1);
 double yieldDsMCErr = f1_DsMC->Integral(DsMass-DsMassCandWidth,DsMass+DsMassCandWidth)/h_DsMassDsMC->GetBinWidth(1) * f1_DsMC->GetParError(0)/f1_DsMC->GetParameter(0);

	TCanvas *c_MCfit=new TCanvas("c_MCfit","c_MCfit");
	c_MCfit->cd();
	h_DsMassDsMC->Draw();

  shiftY=0;
	TLatex *tl_binfitMC=new TLatex();
  tl_binfitMC->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s D_{S} %.0f < pt < %.0f ",s_ppPbPb.Data(),DptLowGL,DptHighGL)); shiftY-=oneshift;
 tl_binfitMC->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Mean=%.4f, frGaus1=%.3f  ",f1_DsMC->GetParameter(1),f1_DsMC->GetParameter(4) ));  shiftY-=oneshift;
  tl_binfitMC->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Sigma1=%.4f, Sigma2=%.4f",f1_DsMC->GetParameter(2), f1_DsMC->GetParameter(3)));  shiftY-=oneshift;



	f1_DsMC->Draw("same");


	c_MCfit->SaveAs(Form("./plots/fit/%s_MC_Mkkpi.png",s_type.Data()));

 return f1_DsMC;

}



TF1* fitDsMass(TH1D *h_DsMassData, TF1 *f1MC,double &FitYield,double &FitYieldErr,double &BkgSubYield,double &BkgSubYieldErr, TString s_ExtraName, int fixSigShape, TF1 *f1_default){

  // TF1* f = new TF1("fMass","[0]*([7]*([9]*TMath::Gaus(x,[1],[2]*(1+[11]))/(sqrt(2*3.1415927)*[2]*(1+[11]))+(1-[9])*TMath::Gaus(x,[1],[10]*(1+[11]))/(sqrt(2*3.1415927)*[10]*(1+[11])))+(1-[7]\
)*TMath::Gaus(x,[1],[8])/(sqrt(2*3.1415927)*[8]))+[3]+[4]*x+[5]*x*x+[6]*x*x*x", 1.7, 2.0);

	// FitYield=472;

	cout<<"check -1"<<endl;
  cout<<"f1MC->GetParameter(1) = "<<f1MC->GetParameter(1)<<endl;
	

  TF1 *f1_DsMC = new TF1("f1_DsMC","[0]*( [4]*TMath::Gaus(x,[1],[2]*(1+[5]))/( sqrt(2*TMath::Pi())*[2]*(1+[5]) ) + (1-[4])*TMath::Gaus(x,[1],[3]*(1+[5]))/(sqrt(2*TMath::Pi())*[3]*(1+[5])  ) )",DsMassFitLow,DsMassFitHigh); // kTRUE nomalize Gaus

  // f1_DsMC->SetParameter(0, h_DsMassDsMC->GetBinWidth(1) * h_DsMassDsMC->Integral( h_DsMassDsMC->FindBin(DsMass-DsMassCandWidth) , h_DsMassDsMC->FindBin(DsMass-DsMassCandWidth) ) ) ;
  f1_DsMC->SetParameter(1,f1MC->GetParameter(1));   // DsMean
  f1_DsMC->SetParLimits(1,DsMass-DsMassCandFitMeanWidth,DsMass+DsMassCandFitMeanWidth);
  f1_DsMC->SetParameter(2,f1MC->GetParameter(2)); // defaut (wider) gaus1 sigma
 f1_DsMC->SetParLimits(2,0.0001,0.1);
  f1_DsMC->SetParameter(3,f1MC->GetParameter(3)); // default (narrow) gaus2 sigma
  f1_DsMC->SetParLimits(3,0.00001,0.1);
  f1_DsMC->SetParameter(4,f1MC->GetParameter(4));  // fraction of Gaus1 / All
  f1_DsMC->SetParLimits(4,0,1);
  f1_DsMC->SetParameter(5,0); // for data /mc discrepancy
//  f1_DsMC->SetParLimits(5,0,0);
//  f1_DsMC->FixParameter(5,0);

  f1_DsMC->SetLineColor(kRed);


	cout<<"check 0"<<endl;



  TF1 *f1_DsSignal = (TF1*)f1_DsMC->Clone("f1_DsSignal");
  TF1 *f1_DsBkg= new TF1("f1_DsBkg","[0]*(1+[1]*x+[2]*(2*x*x-1))");
  TF1 *f1_DsBkg_clone=(TF1*)f1_DsBkg->Clone("f1_DsBkg_clone");  // stupid bug in root, must use clone function
  TF1 *f1_DsMix = new TF1("f1_DsMix","f1_DsSignal+f1_DsBkg");

	// f1_DsMix->SetParLimits(0,0,1000000);
  f1_DsMix->SetParameter(5,0); // width ratio difference for Data/ MC, data width = (1+5)*MC width
  f1_DsMix->SetParLimits(5,-1,1);
	f1_DsMix->SetParameter(0,0.1*h_DsMassData->Integral()*h_DsMassData->GetBinWidth(1));
  f1_DsMix->FixParameter(1,f1_DsMC->GetParameter(1));  // Ds Mass mean
  f1_DsMix->FixParameter(2,f1_DsMC->GetParameter(2));  // Ds Sigma1
  f1_DsMix->FixParameter(3,f1_DsMC->GetParameter(3));  // Ds Sigma2
  f1_DsMix->FixParameter(4,f1_DsMC->GetParameter(4));  // Ds Gaus1/all yield ratio
  f1_DsMix->SetLineColor(kRed);
//  f1_DsMix->SetParameter(6,1000);
  f1_DsMix->SetParameter(7,0);
  f1_DsMix->SetParameter(6,0.9*h_DsMassData->Integral()*h_DsMassData->GetBinWidth(1));


	if(fixSigShape==1 && f1_default!=NULL){
		cout<<"fixSigShape fit"<<endl;
		// cout<<"f1_default->GetParameter(7) = "<<f1_default->GetParameter(7)<<endl;
		// cout<<"f1_default->GetParameter(8) = "<<f1_default->GetParameter(8)<<endl;

		f1_DsMix->SetParameters(0.1*h_DsMassData->Integral()*h_DsMassData->GetBinWidth(1), f1_default->GetParameter(1), f1_default->GetParameter(2), f1_default->GetParameter(3), f1_default->GetParameter(4), f1_default->GetParameter(5), f1_default->GetParameter(6), f1_default->GetParameter(7), f1_default->GetParameter(8) );

		f1_DsMix->FixParameter(1,f1_default->GetParameter(1));  // Ds Mass mean
		f1_DsMix->FixParameter(2,f1_default->GetParameter(2));  
		f1_DsMix->FixParameter(3,f1_default->GetParameter(3));  
		f1_DsMix->FixParameter(4,f1_default->GetParameter(4));  
		f1_DsMix->FixParameter(5,f1_default->GetParameter(5));  

		f1_DsMix->FixParameter(7,f1_default->GetParameter(7));  
		f1_DsMix->FixParameter(8,f1_default->GetParameter(8));  

	}

	cout<<"check 1"<<endl;

  h_DsMassData->Fit("f1_DsMix","QN0","",    DsMassFitLow, DsMassFitHigh );
  h_DsMassData->Fit("f1_DsMix","QN0","",    DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L QN0","",  DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","LN0","",    DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L MN0","",  DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L MN0","",  DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L M IN0","",DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L M IN0","",DsMassFitLow, DsMassFitHigh);

	if(fixSigShape!=1){
  f1_DsMix->ReleaseParameter(5);
  f1_DsMix->ReleaseParameter(1);
	}

  h_DsMassData->Fit("f1_DsMix","QN0","",    DsMassFitLow, DsMassFitHigh );
  h_DsMassData->Fit("f1_DsMix","QN0","",    DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L QN0","",  DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","LN0","",    DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L MN0","",  DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L MN0","",  DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L M IN0","",DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L M I0","",DsMassFitLow, DsMassFitHigh);
// test release background

	f1_DsMix->ReleaseParameter(7);
	f1_DsMix->ReleaseParameter(8);
  h_DsMassData->Fit("f1_DsMix","L M IN0","",DsMassFitLow, DsMassFitHigh);
  h_DsMassData->Fit("f1_DsMix","L M I0","",DsMassFitLow, DsMassFitHigh);


	cout<<"check 2"<<endl;

	f1_DsMix->SetRange(DsMassFitLow,DsMassFitHigh);
	int fitStatus=1;
  int fitIsValid=0;
  TFitResultPtr fitResut;
  double fitPrecision=1.e-8;
  while(fitStatus>0 && fitStatus!=4000){
    TFitter::SetPrecision(fitPrecision);
    fitResut=h_DsMassData->Fit("f1_DsMix","L EMI S0","",DsMassFitLow, DsMassFitHigh);
    fitStatus=fitResut->Status();
    fitIsValid=fitResut->IsValid();
    cout<<"fit precision: "<<TFitter::GetPrecision()<<"   status: "<<fitStatus<<" isValid = "<< fitIsValid <<endl;
    if(fitStatus){
      fitPrecision *= 3;
    }
    if(fitPrecision>1.e-3) {break;}
  }
	fitResut=h_DsMassData->Fit("f1_DsMix","L EMI S0","",DsMassFitLow, DsMassFitHigh);
 cout<<"f1_DsMix->GetParameter(0) = "<<f1_DsMix->GetParameter(0)<<endl;
cout<<"f1_DsSignal->GetParameter(0) = "<<f1_DsSignal->GetParameter(0)<<endl;
  f1_DsSignal->SetParameters(f1_DsMix->GetParameter(0), f1_DsMix->GetParameter(1) , f1_DsMix->GetParameter(2) , f1_DsMix->GetParameter(3) , f1_DsMix->GetParameter(4), f1_DsMix->GetParameter(5));
  f1_DsSignal->SetParError(0,f1_DsMix->GetParError(0));
  f1_DsSignal->SetParError(1,f1_DsMix->GetParError(1));
  f1_DsSignal->SetParError(2,f1_DsMix->GetParError(2));
  f1_DsSignal->SetParError(3,f1_DsMix->GetParError(3));
  f1_DsSignal->SetParError(4,f1_DsMix->GetParError(4));
  f1_DsSignal->SetParError(5,f1_DsMix->GetParError(5));


	cout<<"check 3"<<endl;

  f1_DsSignal->SetFillColor(kOrange-3);
  f1_DsSignal->SetFillStyle(3002);
  f1_DsSignal->SetLineColor(kOrange-3);
  f1_DsSignal->SetLineWidth(2);
  f1_DsSignal->SetLineStyle(2);
  f1_DsSignal->SetRange(DsMassFitLow,DsMassFitHigh); // must setrange before plot or get empty
  double yieldDsSignal= f1_DsSignal->Integral(DsMassFitLow,DsMassFitHigh)/ h_DsMassData->GetBinWidth(1);
  double yieldDsSignalErr = f1_DsSignal->Integral(DsMass-DsMassCandWidth,DsMass+DsMassCandWidth)/h_DsMassData->GetBinWidth(1) * f1_DsSignal->GetParError(0)/f1_DsSignal->GetParameter(0);

  double yieldDsSignalTry = f1_DsSignal->GetParameter(0)/h_DsMassData->GetBinWidth(1);
  double yieldDsSignalTryErr = f1_DsSignal->GetParError(0)/h_DsMassData->GetBinWidth(1);


	cout<<"check 4"<<endl;

  cout<<"parameter 0 = "<<f1_DsSignal->GetParameter(0)<<" BinWidth = "<< h_DsMassData->GetBinWidth(1) <<endl;
  cout<<"yield DsSignal = "<< yieldDsSignal <<" +- "<<yieldDsSignalErr<<endl;
  cout<<"yield DsSignal from Para 0 = "<< yieldDsSignalTry <<" +- "<<yieldDsSignalTryErr<<endl;
  f1_DsBkg->SetParameters(f1_DsMix->GetParameter(6), f1_DsMix->GetParameter(7) ,f1_DsMix->GetParameter(8));
  // f1_DsBkg->SetParameters(f1_DsMix->GetParameter(6), f1_DsMix->GetParameter(7));
  f1_DsBkg->SetParError(0,f1_DsMix->GetParError(6));
  f1_DsBkg->SetParError(1,f1_DsMix->GetParError(7));
  f1_DsBkg->SetParError(2,f1_DsMix->GetParError(8));
//  f1_DsBkg->SetParError(3,f1_DsMix->GetParError(9));

  f1_DsBkg->SetLineColor(1);
  f1_DsBkg->SetLineStyle(2);
  f1_DsBkg->SetRange(DsMassFitLow,DsMassFitHigh);


  double covmat[3][3];
  for(int icov =0; icov<3; icov++){
    for(int jcov =0; jcov<3; jcov++){
      covmat[icov][jcov]=fitResut->CovMatrix(6+icov,6+jcov);
    }
  }

  double *covmatArr=covmat[0];
	

  // double yieldBkg= f1_DsBkg->Integral(DsMassFitLow,DsMassFitHigh)/ h_DsMassData->GetBinWidth(1);
  double yieldBkg= f1_DsBkg->Integral(DsMassFitLow,DsMassFitHigh)/ h_DsMassData->GetBinWidth(1);

  cout<<"yield Bkg = " <<yieldBkg<<endl;

	// bkg sub method

	

	
	double N2Sig=h_DsMassData->Integral(h_DsMassData->FindBin(DsMassMCMean-DsMass2SigWidth),h_DsMassData->FindBin(DsMassMCMean+DsMass2SigWidth));
	double binLow2SigVal=h_DsMassData->GetBinLowEdge(h_DsMassData->FindBin(DsMassMCMean-DsMass2SigWidth));
	double binHigh2SigVal=h_DsMassData->GetBinLowEdge(h_DsMassData->FindBin(DsMassMCMean+DsMass2SigWidth))+h_DsMassData->GetBinWidth(1);
	// double Nbkg2Sig=f1_DsBkg->Integral(DsMassMCMean-DsMass2SigWidth,DsMassMCMean+DsMass2SigWidth)/ h_DsMassData->GetBinWidth(1); // wrong, the edge need to match with th1 bin
	double Nbkg2Sig=f1_DsBkg->Integral(binLow2SigVal,binHigh2SigVal)/ h_DsMassData->GetBinWidth(1);
	double Nbkg2SigError=f1_DsBkg->IntegralError(binLow2SigVal,binHigh2SigVal,f1_DsBkg->GetParameters(), covmatArr )/ h_DsMassData->GetBinWidth(1);
	double NDsbkgSub=N2Sig-Nbkg2Sig;
	double NDsbkgSubError=sqrt(N2Sig+Nbkg2SigError*Nbkg2SigError);
	cout<<"N2Sig = "<<N2Sig<<" +- "<<sqrt(N2Sig)<<endl;
	cout<<"Nbkg2Sig = "<<Nbkg2Sig<<" +- "<<Nbkg2SigError<<endl;
	cout<<"NDsbkgSub = "<<NDsbkgSub<<" +- "<<NDsbkgSubError<<endl;




// #<{(|
	TCanvas *c_DataFit=new TCanvas("c_DataFit","c_DataFit");
	c_DataFit->cd();
	h_DsMassData->Draw();
	f1_DsMix->Draw("same");
	f1_DsBkg->Draw("same");
	f1_DsSignal->Draw("same");

  shiftY=0;
  TLatex *tl_binfitData =new TLatex();
//  tl_binfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s_Data_binfit",str_PbPb.Data()));  shiftY-=oneshift;
  tl_binfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("%s D_{S} %.0f < pt < %.0f ",s_ppPbPb.Data(),DptLowGL,DptHighGL)); shiftY-=oneshift;
  tl_binfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Fit Yield = %.1f #pm %.1f ",yieldDsSignal,yieldDsSignalErr ));  shiftY-=oneshift;
  tl_binfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("BkgSub Yield = %.1f #pm %.1f ",NDsbkgSub,NDsbkgSubError ));  shiftY-=oneshift;
  // tl_binfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,Form("Float width = %.3f",f1_DsMix->GetParameter(5))); shiftY-=oneshift;

  if(s_ExtraName!="")
  {
  tl_binfitData->DrawLatexNDC(textposx+0.28,textposy+shiftY,s_ExtraName); shiftY-=oneshift;
  }
  // SavePlotDirs(c_binfit_Data[count_c_binfit],Form("%s_Data_binfit_pt%.0fto%.0f_%s",str_PbPb.Data(),DptLow,DptHigh,extraName.Data() ),{"SignalFit",str_PbPb,"binfit"} );

  c_DataFit->SaveAs(Form("./plots/fit/%s_Data_Mkkpi_%s.png",s_type.Data(),s_ExtraName.Data()));


	FitYield=yieldDsSignal;
	FitYieldErr=yieldDsSignalErr;
	BkgSubYield=NDsbkgSub;
	BkgSubYieldErr=NDsbkgSubError;

	delete c_DataFit;
// |)}>#

	return f1_DsMix;
}


