TSpline3* spLogDPt;
void divideBinWidth(TH1* h);

void calculateMeasuredDPtWeight(const char* label = "PromptPbPb")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas* c1 = new TCanvas();
  c1->SetLeftMargin(0.14);
  c1->SetBottomMargin(0.14);

  TFile* fInHP; 
  TFile* fInMB;
  if(strcmp(label, "PromptPbPb") == 0)
    {
      fInHP = new TFile("CrossSectionFONLLPbPb_HIN_16_001.root");
      fInMB = new TFile("CrossSectionFONLLPbPbMB_HIN_16_001.root");
    }
  else if(strcmp(label, "PromptPP") == 0)
    {
      fInHP = new TFile("CrossSectionFONLLPP_HIN_16_001.root");
      fInMB = new TFile("CrossSectionFONLLPPMB_HIN_16_001.root");
    }
  else
    {
      cout<<"wrong label: "<<label<<endl;
      return;
    }
  TH1D* hPtSigmaHP = fInHP->Get("hPtSigma");
  TH1D* hPtSigmaMB = fInMB->Get("hPtSigma");
  
  const int nPtBins = 14;
  double ptBins[nPtBins+1] = {2.,3.,4.,5.,6.,8.,10.,12.5,15.,20.,25.,30.,40.,60.,100.};
  TH1D* hPtSigma = new TH1D("hPtSigma", "hPtSigma", nPtBins, ptBins);
  hPtSigma->SetMarkerStyle(21);
  hPtSigma->SetMarkerColor(2);
  hPtSigma->SetLineColor(2);
  hPtSigma->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  
  for(int i=1; i<=hPtSigmaMB->GetXaxis()->GetNbins(); i++)
    {
      hPtSigma->SetBinContent(i, hPtSigmaMB->GetBinContent(i));
      hPtSigma->SetBinError(i, hPtSigmaMB->GetBinError(i));
    }
  for(int i=1; i<=hPtSigmaMB->GetXaxis()->GetNbins(); i++)
    {
      hPtSigma->SetBinContent(i+hPtSigmaMB->GetXaxis()->GetNbins(), hPtSigmaHP->GetBinContent(i));
      hPtSigma->SetBinError(i+hPtSigmaMB->GetXaxis()->GetNbins(), hPtSigmaHP->GetBinError(i));
    }

  c1->SetLogx();
  c1->SetLogy();

  TH1D* hStupidTGraph = new TH1D("hStupidTGraph", "", 100, 0, 100);
  hStupidTGraph->GetXaxis()->SetRangeUser(2, 100);
  hStupidTGraph->GetYaxis()->SetRangeUser(1, 1.e10);
  hStupidTGraph->SetXTitle(hPtSigma->GetXaxis()->GetTitle());
  hStupidTGraph->SetYTitle("dN/dp_{T} (GeV/c)^{-1}");
  hStupidTGraph->Draw();
  hPtSigma->Draw("psame");

  TFile* fFONLLPP = new TFile("FONLL_promptD_5TeV_y1.root");
  double norm=0.557;    //FF of D->D0
  TH1D* hPtFONLL = (TH1D*)fFONLLPP->Get("hdsdpt");
  hPtFONLL->Scale(norm);
  hPtFONLL->Scale(1./4);    // dN/dPt to dN per bin

  TH1D* hPtFONLLRebin = (TH1D*)hPtFONLL->Rebin(nPtBins,"hPtFONLLRebin",ptBins);
  divideBinWidth(hPtFONLLRebin);
  hPtFONLLRebin->SetMarkerStyle(24);
  hPtFONLLRebin->Draw("same");

  c1->SaveAs(Form("DptSpectrumMeasured_%s.pdf", label));

  c1->SetLogy(0);
  TH1D* hMeasuredFonllRatio = hPtSigma->Clone("hMeasuredFonllRatio");
  hMeasuredFonllRatio->Divide(hPtFONLLRebin);
  hMeasuredFonllRatio->GetYaxis()->SetRangeUser(0,4);
  hMeasuredFonllRatio->GetYaxis()->SetTitle("measured / FONLL");
  hMeasuredFonllRatio->Draw("");

  TF1* fPol3Log = new TF1("fPol3Log_measuredOverFonll", "[0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*log(x)*log(x)*log(x)", 0, 400);
  hMeasuredFonllRatio->Fit(fPol3Log);

  c1->SaveAs(Form("DptSpectrumMeasuredFonllRatio_%s.pdf", label));

  double* parameters = fPol3Log->GetParameters();
  
  cout<<endl;
  cout<<"measured Dpt weight"<<endl;
  cout<<Form("%1.5f+%1.5f*log(Dgenpt)+%1.5f*log(Dgenpt)*log(Dgenpt)+%1.5f*log(Dgenpt)*log(Dgenpt)*log(Dgenpt)", parameters[0], parameters[1], parameters[2], parameters[3])<<endl;
  cout<<Form("%1.5f+%1.5f*log(Gpt)+%1.5f*log(Gpt)*log(Gpt)+%1.5f*log(Gpt)*log(Gpt)*log(Gpt)",parameters[0], parameters[1], parameters[2], parameters[3])<<endl;
  cout<<endl;

  TFile* fOut = new TFile(Form("DPtMeasuredOverFonllWeight%s.root", label), "recreate");
  fPol3Log->Write();
  fOut->Close();
}

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
