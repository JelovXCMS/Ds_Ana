#include "uti.h"

TString varname;
TString vartex;
void plotRatioDoubleratio(TString varname_="", TString vartex_="")
{
  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);
  gStyle->SetOptStat(0);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.8);

  varname=varname_;
  vartex=vartex_;
  
  TFile* infPP = new TFile(Form("outfDoubleratio/fPP_%s_DoubleRatio.root",varname.Data()));
  TH1D* hPP = (TH1D*)infPP->Get("hDoubleRatio");
  hPP->SetName(Form("hPP_%s",varname.Data()));
  TFile* infPbPb = new TFile(Form("outfDoubleratio/fPbPb_%s_DoubleRatio.root",varname.Data()));
  TH1D* hPbPb = (TH1D*)infPbPb->Get("hDoubleRatio");
  hPbPb->SetName(Form("hPbPb_%s",varname.Data()));
  TFile* infPPMB = new TFile(Form("outfDoubleratio/fPPMB_%s_DoubleRatio.root",varname.Data()));
  TH1D* hPPMB = (TH1D*)infPPMB->Get("hDoubleRatio");
  hPPMB->SetName(Form("hPPMB_%s",varname.Data()));
  TFile* infPbPbMB = new TFile(Form("outfDoubleratio/fPbPbMB_%s_DoubleRatio.root",varname.Data()));
  TH1D* hPbPbMB = (TH1D*)infPbPbMB->Get("hDoubleRatio");
  hPbPbMB->SetName(Form("hPbPbMB_%s",varname.Data()));

  TH1F* hRatio = (TH1F*)hPbPb->Clone("hRatio");
  hRatio->Divide(hPP);
  TH1F* hRatioMB = (TH1F*)hPbPbMB->Clone("hRatioMB");
  hRatioMB->Divide(hPPMB);

  Int_t varbins = hRatio->GetNbinsX();
  Float_t varstep = hRatio->GetBinWidth(1);
  Float_t varmin = hRatio->GetBinCenter(1)-0.5*varstep;
  Float_t varmax = hRatio->GetBinCenter(varbins)-0.5*varstep;
  Float_t aRatio[varbins],aRatioErr[varbins];
  Float_t aRatioMB[varbins],aRatioMBErr[varbins];
  Float_t aX[varbins],aZero[varbins];
  for(int i=0;i<varbins;i++)
    {
      aX[i] = hRatio->GetBinCenter(i+1)-0.5*varstep;
      aZero[i] = 0;
      aRatio[i] = hRatio->GetBinContent(i+1);
      aRatioErr[i] = hRatio->GetBinError(i+1);
      aRatioMB[i] = hRatioMB->GetBinContent(i+1);
      aRatioMBErr[i] = hRatioMB->GetBinError(i+1);
    }

  TGraphErrors* gRatio = new TGraphErrors(varbins,aX,aRatio,aZero,aRatioErr);
  TH2F* hemptyRatio = new TH2F("hemptyRatio","",20,varmin-0.5*varstep,varmax+0.5*varstep,10.,0.7,1.3);
  hemptyRatio->GetXaxis()->SetTitle(Form("%s",vartex.Data()));
  hemptyRatio->GetYaxis()->SetTitle("DoubleRatio^{PbPb} / DoubleRatio^{PP}");
  hemptyRatio->GetXaxis()->SetTitleOffset(1.);
  hemptyRatio->GetYaxis()->SetTitleOffset(.9);
  hemptyRatio->GetXaxis()->SetTitleSize(0.045);
  hemptyRatio->GetYaxis()->SetTitleSize(0.045);
  hemptyRatio->GetXaxis()->SetTitleFont(42);
  hemptyRatio->GetYaxis()->SetTitleFont(42);
  hemptyRatio->GetXaxis()->SetLabelFont(42);
  hemptyRatio->GetYaxis()->SetLabelFont(42);
  hemptyRatio->GetXaxis()->SetLabelSize(0.04);
  hemptyRatio->GetYaxis()->SetLabelSize(0.04);
  TCanvas* cRatio = new TCanvas("cRatio","",600,600);
  hemptyRatio->Draw();
  gRatio->Draw("psame");
  cRatio->SaveAs(Form("plotRatios/cRatio_%s.pdf",varname.Data()));

  TGraphErrors* gRatioMB = new TGraphErrors(varbins,aX,aRatioMB,aZero,aRatioMBErr);
  TH2F* hemptyRatioMB = new TH2F("hemptyRatioMB","",20,varmin-0.5*varstep,varmax+0.5*varstep,10.,0.7,1.3);
  hemptyRatioMB->GetXaxis()->SetTitle(Form("%s",vartex.Data()));
  hemptyRatioMB->GetYaxis()->SetTitle("DoubleRatio^{PbPb} / DoubleRatio^{PP}");
  hemptyRatioMB->GetXaxis()->SetTitleOffset(1.);
  hemptyRatioMB->GetYaxis()->SetTitleOffset(.9);
  hemptyRatioMB->GetXaxis()->SetTitleSize(0.045);
  hemptyRatioMB->GetYaxis()->SetTitleSize(0.045);
  hemptyRatioMB->GetXaxis()->SetTitleFont(42);
  hemptyRatioMB->GetYaxis()->SetTitleFont(42);
  hemptyRatioMB->GetXaxis()->SetLabelFont(42);
  hemptyRatioMB->GetYaxis()->SetLabelFont(42);
  hemptyRatioMB->GetXaxis()->SetLabelSize(0.04);
  hemptyRatioMB->GetYaxis()->SetLabelSize(0.04);
  TCanvas* cRatioMB = new TCanvas("cRatioMB","",600,600);
  hemptyRatioMB->Draw();
  gRatioMB->Draw("psame");
  cRatioMB->SaveAs(Form("plotRatios/cRatioMB_%s.pdf",varname.Data()));

}

int main(int argc, char *argv[])
{
  if(argc!=3)
    {
      std::cout << "Wrong number of inputs" << std::endl;
      return 1;
    }
  else
    {
      plotRatioDoubleratio(argv[1],argv[2]);
      return 0;
    }
}

