void systematicsCutsSummary()
{
  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetTitleX(.0f);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const int nPtBins = 12;
  float ptBins[nPtBins+1] = {2.,3.,4.,5.,6.,7.,8.,10.,12,20.,40.,60.,100};

  float minAlphaVary = 0.16;
  float maxAlphaVary = 0.2;

  float minDLVaryPbPb = 4.;
  float maxDLVaryPbPb = 12.;
  float minChi2VaryPbPb = 0.05;
  float maxChi2VaryPbPb = 0.11;
  float minTrkDcaVaryPbPb =0.;
  float maxTrkDcaVaryPbPb =3.;

  float minDLVaryPbPbMB = 4.;
  float maxDLVaryPbPbMB = 8.;
  float minChi2VaryPbPbMB = 0.05;
  float maxChi2VaryPbPbMB = 0.11;
  float minTrkDcaVaryPbPbMB = 0.;
  float maxTrkDcaVaryPbPbMB = 4.5;

  float minDLVaryPbPbMBLowPt = 4.;
  float maxDLVaryPbPbMBLowPt = 8.;
  float minChi2VaryPbPbMBLowPt = 0.05;
  float maxChi2VaryPbPbMBLowPt = 0.11;
  float minTrkDcaVaryPbPbMBLowPt = 0.;
  float maxTrkDcaVaryPbPbMBLowPt = 6.;

  float minDLVaryPP = 4.;
  float maxDLVaryPP = 12.;
  float minChi2VaryPP = 0.05;
  float maxChi2VaryPP = 0.11;
  float minTrkDcaVaryPP = 0.;
  float maxTrkDcaVaryPP = 4.5;

  float minDLVaryPPMB = 4.;
  float maxDLVaryPPMB = 10.;
  float minChi2VaryPPMB = 0.05;
  float maxChi2VaryPPMB = 0.11;
  float minTrkDcaVaryPPMB = 0.;
  float maxTrkDcaVaryPPMB = 4.5;

  float minDLVaryPPMBLowPt = 4.;
  float maxDLVaryPPMBLowPt = 10.;
  float minChi2VaryPPMBLowPt = 0.05;
  float maxChi2VaryPPMBLowPt = 0.11;
  float minTrkDcaVaryPPMBLowPt = 0.;
  float maxTrkDcaVaryPPMBLowPt = 4.5;

  // cuts eff
  cout<<"cuts eff"<<endl;
  double maxRatio, minRatio;

  TH1D* hSysCutEffPbPb = new TH1D("hSysCutEffPbPb", ";p_{T} (GeV/c);relative systematics", nPtBins,ptBins);
  TH1D* hSysCutEffPP = new TH1D("hSysCutEffPP", ";p_{T} (GeV/c);relative systematics", nPtBins,ptBins);

  // PbPb
  TFile* fAlphaVaryRatioPbPb = new TFile("outfDoubleratio/fPbPb_Alpha_DoubleRatio.root");
  TH1F* hAlphaVaryRatioPbPb = (TH1F*)fAlphaVaryRatioPbPb->Get("hDoubleRatio");
  hAlphaVaryRatioPbPb->GetXaxis()->SetRangeUser(minAlphaVary, maxAlphaVary);
  maxRatio = hAlphaVaryRatioPbPb->GetMaximum();
  minRatio = hAlphaVaryRatioPbPb->GetMinimum();
  double maxDiffAlphaPbPb = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  TFile* fDecaylengthVaryRatioPbPb = new TFile("outfDoubleratio/fPbPb_Decaylength_DoubleRatio.root");
  TH1F* hDecaylengthVaryRatioPbPb = (TH1F*)fDecaylengthVaryRatioPbPb->Get("hDoubleRatio");
  hDecaylengthVaryRatioPbPb->GetXaxis()->SetRangeUser(minDLVaryPbPb, maxDLVaryPbPb);
  maxRatio = hDecaylengthVaryRatioPbPb->GetMaximum();
  minRatio = hDecaylengthVaryRatioPbPb->GetMinimum();
  double maxDiffDecaylengthPbPb = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  TFile* fChi2VaryRatioPbPb = new TFile("outfDoubleratio/fPbPb_Chi2_DoubleRatio.root");
  TH1F* hChi2VaryRatioPbPb = (TH1F*)fChi2VaryRatioPbPb->Get("hDoubleRatio");
  hChi2VaryRatioPbPb->GetXaxis()->SetRangeUser(minChi2VaryPbPb, maxChi2VaryPbPb);
  maxRatio = hChi2VaryRatioPbPb->GetMaximum();
  minRatio = hChi2VaryRatioPbPb->GetMinimum();
  double maxDiffChi2PbPb = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  TFile* fTrkDcaVaryRatioPbPb = new TFile("outfDoubleratio/fPbPb_TrkDca_DoubleRatio.root");
  TH1F* hTrkDcaVaryRatioPbPb = (TH1F*)fTrkDcaVaryRatioPbPb->Get("hDoubleRatio");
  hTrkDcaVaryRatioPbPb->GetXaxis()->SetRangeUser(minTrkDcaVaryPbPb, maxTrkDcaVaryPbPb);
  maxRatio = hTrkDcaVaryRatioPbPb->GetMaximum();
  minRatio = hTrkDcaVaryRatioPbPb->GetMinimum();
  double maxDiffTrkDcaPbPb = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  double diffCutsEffPbPb = sqrt(maxDiffAlphaPbPb*maxDiffAlphaPbPb + maxDiffDecaylengthPbPb*maxDiffDecaylengthPbPb + maxDiffChi2PbPb*maxDiffChi2PbPb +maxDiffTrkDcaPbPb*maxDiffTrkDcaPbPb);

  // PbPbMB
  TFile* fAlphaVaryRatioPbPbMB = new TFile("outfDoubleratio/fPbPbMB_Alpha_DoubleRatio.root");
  TH1F* hAlphaVaryRatioPbPbMB = (TH1F*)fAlphaVaryRatioPbPbMB->Get("hDoubleRatio");
  hAlphaVaryRatioPbPbMB->GetXaxis()->SetRangeUser(minAlphaVary, maxAlphaVary);
  maxRatio = hAlphaVaryRatioPbPbMB->GetMaximum();
  minRatio = hAlphaVaryRatioPbPbMB->GetMinimum();
  double maxDiffAlphaPbPbMB = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  TFile* fDecaylengthVaryRatioPbPbMB = new TFile("outfDoubleratio/fPbPbMB_Decaylength_DoubleRatio.root");
  TH1F* hDecaylengthVaryRatioPbPbMB = (TH1F*)fDecaylengthVaryRatioPbPbMB->Get("hDoubleRatio");
  hDecaylengthVaryRatioPbPbMB->GetXaxis()->SetRangeUser(minDLVaryPbPbMB, maxDLVaryPbPbMB);
  maxRatio = hDecaylengthVaryRatioPbPbMB->GetMaximum();
  minRatio = hDecaylengthVaryRatioPbPbMB->GetMinimum();
  double maxDiffDecaylengthPbPbMB = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  TFile* fChi2VaryRatioPbPbMB = new TFile("outfDoubleratio/fPbPbMB_Chi2_DoubleRatio.root");
  TH1F* hChi2VaryRatioPbPbMB = (TH1F*)fChi2VaryRatioPbPbMB->Get("hDoubleRatio");
  hChi2VaryRatioPbPbMB->GetXaxis()->SetRangeUser(minChi2VaryPbPbMB, maxChi2VaryPbPbMB);
  maxRatio = hChi2VaryRatioPbPbMB->GetMaximum();
  minRatio = hChi2VaryRatioPbPbMB->GetMinimum();
  double maxDiffChi2PbPbMB = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  TFile* fTrkDcaVaryRatioPbPbMB = new TFile("outfDoubleratio/fPbPbMB_TrkDca_DoubleRatio.root");
  TH1F* hTrkDcaVaryRatioPbPbMB = (TH1F*)fTrkDcaVaryRatioPbPbMB->Get("hDoubleRatio");
  hTrkDcaVaryRatioPbPbMB->GetXaxis()->SetRangeUser(minTrkDcaVaryPbPbMB, maxTrkDcaVaryPbPbMB);
  maxRatio = hTrkDcaVaryRatioPbPbMB->GetMaximum();
  minRatio = hTrkDcaVaryRatioPbPbMB->GetMinimum();
  double maxDiffTrkDcaPbPbMB = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  double diffCutsEffPbPbMB = sqrt(maxDiffAlphaPbPbMB*maxDiffAlphaPbPbMB + maxDiffDecaylengthPbPbMB*maxDiffDecaylengthPbPbMB + maxDiffChi2PbPbMB*maxDiffChi2PbPbMB + maxDiffTrkDcaPbPbMB*maxDiffTrkDcaPbPbMB);

  // PbPbMBLowPt
  TFile* fAlphaVaryRatioPbPbMBLowPt = new TFile("outfDoubleratio/fPbPbMBLowPt_Alpha_DoubleRatio.root");
  TH1F* hAlphaVaryRatioPbPbMBLowPt = (TH1F*)fAlphaVaryRatioPbPbMBLowPt->Get("hDoubleRatio");
  hAlphaVaryRatioPbPbMBLowPt->GetXaxis()->SetRangeUser(minAlphaVary, maxAlphaVary);
  maxRatio = hAlphaVaryRatioPbPbMBLowPt->GetMaximum();
  minRatio = hAlphaVaryRatioPbPbMBLowPt->GetMinimum();
  double maxDiffAlphaPbPbMBLowPt = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  TFile* fDecaylengthVaryRatioPbPbMBLowPt = new TFile("outfDoubleratio/fPbPbMBLowPt_Decaylength_DoubleRatio.root");
  TH1F* hDecaylengthVaryRatioPbPbMBLowPt = (TH1F*)fDecaylengthVaryRatioPbPbMBLowPt->Get("hDoubleRatio");
  hDecaylengthVaryRatioPbPbMBLowPt->GetXaxis()->SetRangeUser(minDLVaryPbPbMBLowPt, maxDLVaryPbPbMBLowPt);
  maxRatio = hDecaylengthVaryRatioPbPbMBLowPt->GetMaximum();
  minRatio = hDecaylengthVaryRatioPbPbMBLowPt->GetMinimum();
  double maxDiffDecaylengthPbPbMBLowPt = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  TFile* fChi2VaryRatioPbPbMBLowPt = new TFile("outfDoubleratio/fPbPbMBLowPt_Chi2_DoubleRatio.root");
  TH1F* hChi2VaryRatioPbPbMBLowPt = (TH1F*)fChi2VaryRatioPbPbMBLowPt->Get("hDoubleRatio");
  hChi2VaryRatioPbPbMBLowPt->GetXaxis()->SetRangeUser(minChi2VaryPbPbMBLowPt, maxChi2VaryPbPbMBLowPt);
  maxRatio = hChi2VaryRatioPbPbMBLowPt->GetMaximum();
  minRatio = hChi2VaryRatioPbPbMBLowPt->GetMinimum();
  double maxDiffChi2PbPbMBLowPt = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  TFile* fTrkDcaVaryRatioPbPbMBLowPt = new TFile("outfDoubleratio/fPbPbMBLowPt_TrkDca_DoubleRatio.root");
  TH1F* hTrkDcaVaryRatioPbPbMBLowPt = (TH1F*)fTrkDcaVaryRatioPbPbMBLowPt->Get("hDoubleRatio");
  hTrkDcaVaryRatioPbPbMBLowPt->GetXaxis()->SetRangeUser(minTrkDcaVaryPbPbMBLowPt, maxTrkDcaVaryPbPbMBLowPt);
  maxRatio = hTrkDcaVaryRatioPbPbMBLowPt->GetMaximum();
  minRatio = hTrkDcaVaryRatioPbPbMBLowPt->GetMinimum();
  double maxDiffTrkDcaPbPbMBLowPt = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  double diffCutsEffPbPbMBLowPt = sqrt(maxDiffAlphaPbPbMBLowPt*maxDiffAlphaPbPbMBLowPt + maxDiffDecaylengthPbPbMBLowPt*maxDiffDecaylengthPbPbMBLowPt + maxDiffChi2PbPbMBLowPt*maxDiffChi2PbPbMBLowPt + maxDiffTrkDcaPbPbMBLowPt*maxDiffTrkDcaPbPbMBLowPt);

  // pp
  TFile* fAlphaVaryRatioPP = new TFile("outfDoubleratio/fPP_Alpha_DoubleRatio.root");
  TH1F* hAlphaVaryRatioPP = (TH1F*)fAlphaVaryRatioPP->Get("hDoubleRatio");
  hAlphaVaryRatioPP->GetXaxis()->SetRangeUser(minAlphaVary, maxAlphaVary);
  maxRatio = hAlphaVaryRatioPP->GetMaximum();
  minRatio = hAlphaVaryRatioPP->GetMinimum();
  double maxDiffAlphaPP = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  TFile* fDecaylengthVaryRatioPP = new TFile("outfDoubleratio/fPP_Decaylength_DoubleRatio.root");
  TH1F* hDecaylengthVaryRatioPP = (TH1F*)fDecaylengthVaryRatioPP->Get("hDoubleRatio");
  hDecaylengthVaryRatioPP->GetXaxis()->SetRangeUser(minDLVaryPP, maxDLVaryPP);
  maxRatio = hDecaylengthVaryRatioPP->GetMaximum();
  minRatio = hDecaylengthVaryRatioPP->GetMinimum();
  double maxDiffDecaylengthPP = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  TFile* fChi2VaryRatioPP = new TFile("outfDoubleratio/fPP_Chi2_DoubleRatio.root");
  TH1F* hChi2VaryRatioPP = (TH1F*)fChi2VaryRatioPP->Get("hDoubleRatio");
  hChi2VaryRatioPP->GetXaxis()->SetRangeUser(minChi2VaryPP, maxChi2VaryPP);
  maxRatio = hChi2VaryRatioPP->GetMaximum();
  minRatio = hChi2VaryRatioPP->GetMinimum();
  double maxDiffChi2PP = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  TFile* fTrkDcaVaryRatioPP = new TFile("outfDoubleratio/fPP_TrkDca_DoubleRatio.root");
  TH1F* hTrkDcaVaryRatioPP = (TH1F*)fTrkDcaVaryRatioPP->Get("hDoubleRatio");
  hTrkDcaVaryRatioPP->GetXaxis()->SetRangeUser(minTrkDcaVaryPP, maxTrkDcaVaryPP);
  maxRatio = hTrkDcaVaryRatioPP->GetMaximum();
  minRatio = hTrkDcaVaryRatioPP->GetMinimum();
  double maxDiffTrkDcaPP = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  double diffCutsEffPP = sqrt(maxDiffAlphaPP*maxDiffAlphaPP + maxDiffDecaylengthPP*maxDiffDecaylengthPP + maxDiffChi2PP*maxDiffChi2PP + maxDiffTrkDcaPP*maxDiffTrkDcaPP);

  // ppMB
  TFile* fAlphaVaryRatioPPMB = new TFile("outfDoubleratio/fPPMB_Alpha_DoubleRatio.root");
  TH1F* hAlphaVaryRatioPPMB = (TH1F*)fAlphaVaryRatioPPMB->Get("hDoubleRatio");
  hAlphaVaryRatioPPMB->GetXaxis()->SetRangeUser(minAlphaVary, maxAlphaVary);
  maxRatio = hAlphaVaryRatioPPMB->GetMaximum();
  minRatio = hAlphaVaryRatioPPMB->GetMinimum();
  double maxDiffAlphaPPMB = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  TFile* fDecaylengthVaryRatioPPMB = new TFile("outfDoubleratio/fPPMB_Decaylength_DoubleRatio.root");
  TH1F* hDecaylengthVaryRatioPPMB = (TH1F*)fDecaylengthVaryRatioPPMB->Get("hDoubleRatio");
  hDecaylengthVaryRatioPPMB->GetXaxis()->SetRangeUser(minDLVaryPPMB, maxDLVaryPPMB);
  maxRatio = hDecaylengthVaryRatioPPMB->GetMaximum();
  minRatio = hDecaylengthVaryRatioPPMB->GetMinimum();
  double maxDiffDecaylengthPPMB = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  TFile* fChi2VaryRatioPPMB = new TFile("outfDoubleratio/fPPMB_Chi2_DoubleRatio.root");
  TH1F* hChi2VaryRatioPPMB = (TH1F*)fChi2VaryRatioPPMB->Get("hDoubleRatio");
  hChi2VaryRatioPPMB->GetXaxis()->SetRangeUser(minChi2VaryPPMB, maxChi2VaryPPMB);
  maxRatio = hChi2VaryRatioPPMB->GetMaximum();
  minRatio = hChi2VaryRatioPPMB->GetMinimum();
  double maxDiffChi2PPMB = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  TFile* fTrkDcaVaryRatioPPMB = new TFile("outfDoubleratio/fPPMB_TrkDca_DoubleRatio.root");
  TH1F* hTrkDcaVaryRatioPPMB = (TH1F*)fTrkDcaVaryRatioPPMB->Get("hDoubleRatio");
  hTrkDcaVaryRatioPPMB->GetXaxis()->SetRangeUser(minTrkDcaVaryPPMB, maxTrkDcaVaryPPMB);
  maxRatio = hTrkDcaVaryRatioPPMB->GetMaximum();
  minRatio = hTrkDcaVaryRatioPPMB->GetMinimum();
  double maxDiffTrkDcaPPMB = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  double diffCutsEffPPMB = sqrt(maxDiffAlphaPPMB*maxDiffAlphaPPMB + maxDiffDecaylengthPPMB*maxDiffDecaylengthPPMB + maxDiffChi2PPMB*maxDiffChi2PPMB + maxDiffTrkDcaPPMB*maxDiffTrkDcaPPMB);

  // ppMBLowPt
  TFile* fAlphaVaryRatioPPMBLowPt = new TFile("outfDoubleratio/fPPMBLowPt_Alpha_DoubleRatio.root");
  TH1F* hAlphaVaryRatioPPMBLowPt = (TH1F*)fAlphaVaryRatioPPMBLowPt->Get("hDoubleRatio");
  hAlphaVaryRatioPPMBLowPt->GetXaxis()->SetRangeUser(minAlphaVary, maxAlphaVary);
  maxRatio = hAlphaVaryRatioPPMBLowPt->GetMaximum();
  minRatio = hAlphaVaryRatioPPMBLowPt->GetMinimum();
  double maxDiffAlphaPPMBLowPt = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  TFile* fDecaylengthVaryRatioPPMBLowPt = new TFile("outfDoubleratio/fPPMBLowPt_Decaylength_DoubleRatio.root");
  TH1F* hDecaylengthVaryRatioPPMBLowPt = (TH1F*)fDecaylengthVaryRatioPPMBLowPt->Get("hDoubleRatio");
  hDecaylengthVaryRatioPPMBLowPt->GetXaxis()->SetRangeUser(minDLVaryPPMBLowPt, maxDLVaryPPMBLowPt);
  maxRatio = hDecaylengthVaryRatioPPMBLowPt->GetMaximum();
  minRatio = hDecaylengthVaryRatioPPMBLowPt->GetMinimum();
  double maxDiffDecaylengthPPMBLowPt = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  TFile* fChi2VaryRatioPPMBLowPt = new TFile("outfDoubleratio/fPPMBLowPt_Chi2_DoubleRatio.root");
  TH1F* hChi2VaryRatioPPMBLowPt = (TH1F*)fChi2VaryRatioPPMBLowPt->Get("hDoubleRatio");
  hChi2VaryRatioPPMBLowPt->GetXaxis()->SetRangeUser(minChi2VaryPPMBLowPt, maxChi2VaryPPMBLowPt);
  maxRatio = hChi2VaryRatioPPMBLowPt->GetMaximum();
  minRatio = hChi2VaryRatioPPMBLowPt->GetMinimum();
  double maxDiffChi2PPMBLowPt = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  TFile* fTrkDcaVaryRatioPPMBLowPt = new TFile("outfDoubleratio/fPPMBLowPt_TrkDca_DoubleRatio.root");
  TH1F* hTrkDcaVaryRatioPPMBLowPt = (TH1F*)fTrkDcaVaryRatioPPMBLowPt->Get("hDoubleRatio");
  hTrkDcaVaryRatioPPMBLowPt->GetXaxis()->SetRangeUser(minTrkDcaVaryPPMBLowPt, maxTrkDcaVaryPPMBLowPt);
  maxRatio = hTrkDcaVaryRatioPPMBLowPt->GetMaximum();
  minRatio = hTrkDcaVaryRatioPPMBLowPt->GetMinimum();
  double maxDiffTrkDcaPPMBLowPt = TMath::Max(fabs(maxRatio-1), fabs(minRatio-1));

  double diffCutsEffPPMBLowPt = sqrt(maxDiffAlphaPPMBLowPt*maxDiffAlphaPPMBLowPt + maxDiffDecaylengthPPMBLowPt*maxDiffDecaylengthPPMBLowPt + maxDiffChi2PPMBLowPt*maxDiffChi2PPMBLowPt + maxDiffTrkDcaPPMBLowPt*maxDiffTrkDcaPPMBLowPt);

  for(int i=1; i<=hSysCutEffPbPb->GetXaxis()->GetNbins(); i++)
    {
      hSysCutEffPbPb->SetBinContent(i, 0);
      hSysCutEffPP->SetBinContent(i, 0);
      if(hSysCutEffPbPb->GetXaxis()->GetBinCenter(i) < 6)
        {
          hSysCutEffPbPb->SetBinError(i, diffCutsEffPbPbMBLowPt);
          hSysCutEffPP->SetBinError(i, diffCutsEffPPMBLowPt);
        }
      else if(hSysCutEffPbPb->GetXaxis()->GetBinCenter(i) < 20)
	{
	  hSysCutEffPbPb->SetBinError(i, diffCutsEffPbPbMB);
	  hSysCutEffPP->SetBinError(i, diffCutsEffPPMB);
	}
      else
	{
	  hSysCutEffPbPb->SetBinError(i, diffCutsEffPbPb);
	  hSysCutEffPP->SetBinError(i, diffCutsEffPP);
	}
    }

  cout<<endl;
  cout<<"pp 2-6 & ";
  cout<<Form("%1.2f-%1.2f", minAlphaVary, maxAlphaVary)<<" & "<<Form("%1.1f", maxDiffAlphaPPMBLowPt*100)<<"$\\%$ & ";
  cout<<Form("%1.1f-%1.1f", minDLVaryPPMBLowPt, maxDLVaryPPMBLowPt)<<" & "<<Form("%1.1f", maxDiffDecaylengthPPMBLowPt*100)<<"$\\%$ & ";
  cout<<Form("%1.2f-%1.2f", minChi2VaryPPMBLowPt, maxChi2VaryPPMBLowPt)<<" & "<<Form("%1.1f", maxDiffChi2PPMBLowPt*100)<<"$\\%$ & ";
  cout<<Form("%1.1f-%1.1f", minTrkDcaVaryPPMBLowPt, maxTrkDcaVaryPPMBLowPt)<<" & "<<Form("%1.1f", maxDiffTrkDcaPPMBLowPt*100)<<"$\\%$ & ";
  cout<<Form("%1.1f", diffCutsEffPPMBLowPt*100)<<"$\\%$ \\\\"<<endl;
  cout<<"\\hline"<<endl;

  cout<<"pp 6-20 & ";
  cout<<Form("%1.2f-%1.2f", minAlphaVary, maxAlphaVary)<<" & "<<Form("%1.1f", maxDiffAlphaPPMB*100)<<"$\\%$ & ";
  cout<<Form("%1.1f-%1.1f", minDLVaryPPMB, maxDLVaryPPMB)<<" & "<<Form("%1.1f", maxDiffDecaylengthPPMB*100)<<"$\\%$ & ";
  cout<<Form("%1.2f-%1.2f", minChi2VaryPPMB, maxChi2VaryPPMB)<<" & "<<Form("%1.1f", maxDiffChi2PPMB*100)<<"$\\%$ & ";
  cout<<Form("%1.1f-%1.1f", minTrkDcaVaryPPMB, maxTrkDcaVaryPPMB)<<" & "<<Form("%1.1f", maxDiffTrkDcaPPMB*100)<<"$\\%$ & ";
  cout<<Form("%1.1f", diffCutsEffPPMB*100)<<"$\\%$ \\\\"<<endl;
  cout<<"\\hline"<<endl;

  cout<<"pp 20-100 & ";
  cout<<Form("%1.2f-%1.2f", minAlphaVary, maxAlphaVary)<<" & "<<Form("%1.1f", maxDiffAlphaPP*100)<<"$\\%$ & ";
  cout<<Form("%1.1f-%1.1f", minDLVaryPP, maxDLVaryPP)<<" & "<<Form("%1.1f", maxDiffDecaylengthPP*100)<<"$\\%$ & ";
  cout<<Form("%1.2f-%1.2f", minChi2VaryPP, maxChi2VaryPP)<<" & "<<Form("%1.1f", maxDiffChi2PP*100)<<"$\\%$ & ";
  cout<<Form("%1.1f-%1.1f", minTrkDcaVaryPP, maxTrkDcaVaryPP)<<" & "<<Form("%1.1f", maxDiffTrkDcaPP*100)<<"$\\%$ & ";
  cout<<Form("%1.1f", diffCutsEffPP*100)<<"$\\%$ \\\\"<<endl;
  cout<<"\\hline"<<endl;

  cout<<"PbPb 2-6 & ";
  cout<<Form("%1.2f-%1.2f", minAlphaVary, maxAlphaVary)<<" & "<<Form("%1.1f", maxDiffAlphaPbPbMBLowPt*100)<<"$\\%$ & ";
  cout<<Form("%1.1f-%1.1f", minDLVaryPbPbMBLowPt, maxDLVaryPbPbMBLowPt)<<" & "<<Form("%1.1f", maxDiffDecaylengthPbPbMBLowPt*100)<<"$\\%$ & ";
  cout<<Form("%1.2f-%1.2f", minChi2VaryPbPbMBLowPt, maxChi2VaryPbPbMBLowPt)<<" & "<<Form("%1.1f", maxDiffChi2PbPbMBLowPt*100)<<"$\\%$ & ";
  cout<<Form("%1.1f-%1.1f", minTrkDcaVaryPbPbMBLowPt, maxTrkDcaVaryPbPbMBLowPt)<<" & "<<Form("%1.1f", maxDiffTrkDcaPbPbMBLowPt*100)<<"$\\%$ & ";
  cout<<Form("%1.1f", diffCutsEffPbPbMBLowPt*100)<<"$\\%$ \\\\"<<endl;
  cout<<"\\hline"<<endl;

  cout<<"PbPb 6-20 & ";
  cout<<Form("%1.2f-%1.2f", minAlphaVary, maxAlphaVary)<<" & "<<Form("%1.1f", maxDiffAlphaPbPbMB*100)<<"$\\%$ & ";
  cout<<Form("%1.1f-%1.1f", minDLVaryPbPbMB, maxDLVaryPbPbMB)<<" & "<<Form("%1.1f", maxDiffDecaylengthPbPbMB*100)<<"$\\%$ & ";
  cout<<Form("%1.2f-%1.2f", minChi2VaryPbPbMB, maxChi2VaryPbPbMB)<<" & "<<Form("%1.1f", maxDiffChi2PbPbMB*100)<<"$\\%$ & ";
  cout<<Form("%1.1f-%1.1f", minTrkDcaVaryPbPbMB, maxTrkDcaVaryPbPbMB)<<" & "<<Form("%1.1f", maxDiffTrkDcaPbPbMB*100)<<"$\\%$ & ";
  cout<<Form("%1.1f", diffCutsEffPbPbMB*100)<<"$\\%$ \\\\"<<endl;
  cout<<"\\hline"<<endl;

  cout<<"PbPb 20-100 & ";
  cout<<Form("%1.2f-%1.2f", minAlphaVary, maxAlphaVary)<<" & "<<Form("%1.1f", maxDiffAlphaPbPb*100)<<"$\\%$ & ";
  cout<<Form("%1.1f-%1.1f", minDLVaryPbPb, maxDLVaryPbPb)<<" & "<<Form("%1.1f", maxDiffDecaylengthPbPb*100)<<"$\\%$ & ";
  cout<<Form("%1.2f-%1.2f", minChi2VaryPbPb, maxChi2VaryPbPb)<<" & "<<Form("%1.1f", maxDiffChi2PbPb*100)<<"$\\%$ & ";
  cout<<Form("%1.1f-%1.1f", minTrkDcaVaryPbPb, maxTrkDcaVaryPbPb)<<" & "<<Form("%1.1f", maxDiffTrkDcaPbPb*100)<<"$\\%$ & ";
  cout<<Form("%1.1f", diffCutsEffPbPb*100)<<"$\\%$ \\\\"<<endl;
  cout<<"\\hline"<<endl;


  /*
  cout<<"pp MB & "<<Form("%1.1f", maxDiffAlphaPPMB*100)<<"$\\%$ & "<<Form("%1.1f", maxDiffDecaylengthPPMB*100)<<"$\\%$ & "<<Form("%1.1f", maxDiffChi2PPMB*100)<<"$\\%$ & "<<Form("%1.1f", maxDiffTrkDcaPPMB*100)<<"$\\%$ & "<<Form("%1.1f", diffCutsEffPPMB*100)<<"$\\%$ \\\\"<<endl;
  cout<<"\\hline"<<endl;
  cout<<"pp MB low p_{T} & "<<Form("%1.1f", maxDiffAlphaPPMBLowPt*100)<<"$\\%$ & "<<Form("%1.1f", maxDiffDecaylengthPPMBLowPt*100)<<"$\\%$ & "<<Form("%1.1f", maxDiffChi2PPMBLowPt*100)<<"$\\%$ & "<<Form("%1.1f", maxDiffTrkDcaPPMBLowPt*100)<<"$\\%$ & "<<Form("%1.1f", diffCutsEffPPMBLowPt*100)<<"$\\%$ \\\\"<<endl;
  cout<<"\\hline"<<endl;
  cout<<"PbPb & "<<Form("%1.1f", maxDiffAlphaPbPb*100)<<"$\\%$ & "<<Form("%1.1f", maxDiffDecaylengthPbPb*100)<<"$\\%$ & "<<Form("%1.1f", maxDiffChi2PbPb*100)<<"$\\%$ & "<<Form("%1.1f", maxDiffTrkDcaPbPb*100)<<"$\\%$ & "<<Form("%1.1f", diffCutsEffPbPb*100)<<"$\\%$ \\\\"<<endl;
  cout<<"\\hline"<<endl;
  cout<<"PbPb MB & "<<Form("%1.1f", maxDiffAlphaPbPbMB*100)<<"$\\%$ & "<<Form("%1.1f", maxDiffDecaylengthPbPbMB*100)<<"$\\%$ & "<<Form("%1.1f", maxDiffChi2PbPbMB*100)<<"$\\%$ & "<<Form("%1.1f", maxDiffTrkDcaPbPbMB*100)<<"$\\%$ & "<<Form("%1.1f", diffCutsEffPbPbMB*100)<<"$\\%$ \\\\"<<endl;
  cout<<"\\hline"<<endl;
  cout<<"PbPb MB low p_{T} & "<<Form("%1.1f", maxDiffAlphaPbPbMBLowPt*100)<<"$\\%$ & "<<Form("%1.1f", maxDiffDecaylengthPbPbMBLowPt*100)<<"$\\%$ & "<<Form("%1.1f", maxDiffChi2PbPbMBLowPt*100)<<"$\\%$ & "<<Form("%1.1f", maxDiffTrkDcaPbPbMBLowPt*100)<<"$\\%$ & "<<Form("%1.1f", diffCutsEffPbPbMBLowPt*100)<<"$\\%$ \\\\"<<endl;
  cout<<"\\hline"<<endl;
  cout<<endl;
  */

  TCanvas* c1 = new TCanvas();
  
  TFile* fOut = new TFile("systematicsCuts.root", "recreate");
  fOut->WriteTObject(hSysCutEffPbPb);
  fOut->WriteTObject(hSysCutEffPP);
  fOut->Close();
}










