#ifndef _ANA_BFEEDDOWN_UTI_H_
#define _ANA_BFEEDDOWN_UTI_H_

#include <iostream>
#include <iomanip>
#include <utility>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TBox.h>
#include <TCut.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TH3D.h>

// Remove error
void removeError(TH1F* h)
{
  for(int i=1;i<=h->GetNbinsX();i++)
    {
      h->SetBinError(i,0);
    }	
}

// divide by bin width
void divideBinWidth(TH1* h)
{
  h->Sumw2();
  for(int i=1;i<=h->GetNbinsX();i++)
    {
      Float_t val = h->GetBinContent(i);
      Float_t valErr = h->GetBinError(i);
      val/=h->GetBinWidth(i);
      valErr/=h->GetBinWidth(i);
      h->SetBinContent(i,val);
      h->SetBinError(i,valErr);
    }
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

void MutiplyBinWidth(TH1* h)
{
  h->Sumw2();
  for(int i=1;i<=h->GetNbinsX();i++)
    {
      Float_t val = h->GetBinContent(i);
      Float_t valErr = h->GetBinError(i);
      val*=h->GetBinWidth(i);
      valErr*=h->GetBinWidth(i);
      h->SetBinContent(i,val);
      h->SetBinError(i,valErr);
    }
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}


// make a histogram from TF1 function
TH1F* functionHist(TF1* f, TH1F* h, TString fHistname)
{
  TH1F* hF = (TH1F*)h->Clone(fHistname);
  for (int i=1;i<=h->GetNbinsX();i++)
    {
      Double_t var = f->Integral(h->GetBinLowEdge(i),h->GetBinLowEdge(i+1))/h->GetBinWidth(i);
      hF->SetBinContent(i,var);
      hF->SetBinError(i,0);
    }
  return hF;
}

TLegend* myLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
  TLegend* leg = new TLegend(x1,y1,x2,y2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  return leg; 
}

Double_t ErrorPro_aAPlusbB(Double_t AErr, Double_t BErr, Double_t aconst=1, Double_t bconst=1)
{ // f=aA+bB , if A,B are independent
  return sqrt(pow(AErr*aconst,2) + pow(BErr*bconst,2));
}

Double_t ErrorPro_aAPlusbB_Corr(Double_t AErr, Double_t BErr, Double_t aconst=1, Double_t bconst=1)
{ // f=aA+bB , if A,B are independent
  return sqrt(abs( pow(AErr*aconst,2) - pow(BErr*bconst,2) ) );
}



Double_t ErrorPro_AoverB(Double_t AErr, Double_t BErr, Double_t AVal, Double_t BVal) 
{ // f=A/B, if A,B are independent
	return AVal/BVal*sqrt(pow(AErr/AVal,2) + pow(BErr/BVal,2) ); 
}

Double_t ErrorPro_AoverB_Corr(Double_t AErr, Double_t BErr, Double_t AVal, Double_t BVal) 
{ // f=A/B, if A,B are totally correlated
	Double_t CErr=AErr;
	if(BErr<AErr) CErr=BErr; // the smaller one as subset

	return AVal/BVal*sqrt(pow(AErr/AVal,2) + pow(BErr/BVal,2) - 2/AVal/BVal*CErr*CErr  ); 
}

Double_t ErrorPro_AtimesB(Double_t AErr, Double_t BErr, Double_t AVal, Double_t BVal) 
{ // f=A/B, if A,B are independent
	return AVal*BVal*sqrt(pow(AErr/AVal,2) + pow(BErr/BVal,2) ); 
}

Double_t ErrorPro_AoverAplusB(Double_t AErr, Double_t BErr, Double_t AVal, Double_t BVal)
{ // f=A/(A+B), if  A,B are independent

	return AVal/(AVal+BVal)*sqrt( pow(AErr/AVal,2) + (AErr*AErr+BErr*BErr)/(pow(AVal+BVal,2)) - 2/(AVal*(AVal+BVal))*AErr*AErr );

}



#endif
