void calculateDPtSampleWeight(const char* label = "PromptPP")
{
  TFile* f1 = new TFile(Form("DGenPt%s0.root", label));
  float samplePtCuts[7] = {0, 2, 4, 10, 20, 40, 60};

  TCanvas* c1 = new TCanvas();
  c1->SetLogx();
  c1->SetLogy();
  
  TH1D* hDGenPt_sample[7];
  for(int i=0; i<7; i++)
    {
      hDGenPt_sample[i] = (TH1D*)f1->Get(Form("hDGenPt%s0_sampleDPt%1.0f", label, samplePtCuts[i]));
      if(i==0)
	hDGenPt_sample[i]->Draw();
      else
	hDGenPt_sample[i]->Draw("same");
    }

  c1->SaveAs(Form("DGenPt%s_samples.pdf", label));

  TH1D* hDGenPtWeighted = (TH1D*)hDGenPt_sample[0]->Clone(Form("hDGenPtWeighted%s", label));
  if(strcmp(label, "PromptPP")==0 || strcmp(label, "PromptPbPb")==0)
    hDGenPtWeighted->Scale(5.136e9/2.e6); //cross section with D filter (pb) / n events 
  if(strcmp(label, "NonPromptPP")==0 || strcmp(label, "NonPromptPbPb")==0)
    hDGenPtWeighted->Scale(1.72e8/2.e6); //cross section with D filter (pb) / n events

  float sampleWeights[7];
  int ptCutBins[8];
  for(int i=0; i<7; i++)
    ptCutBins[i] = hDGenPtWeighted->GetXaxis()->FindBin(samplePtCuts[i]+0.00001);
  ptCutBins[7] = hDGenPtWeighted->GetXaxis()->GetNbins()+1;

  for(int i=0; i<7; i++)
    {
      TH1D* hDGenPtTotal = (TH1D*)hDGenPt_sample[0]->Clone("hDGenPtTotal");
      for(int j=1; j<=i; j++)
	hDGenPtTotal->Add(hDGenPt_sample[j]);

      hDGenPtWeighted->Draw();
      hDGenPtTotal->Draw("same");

      sampleWeights[i] = hDGenPtWeighted->Integral(ptCutBins[i], ptCutBins[7]-1) / hDGenPtTotal->Integral(ptCutBins[i], ptCutBins[7]-1);
      cout<<"sampleWeight "<<i<<" : "<<sampleWeights[i]<<endl;

      delete hDGenPtWeighted;
      
      hDGenPtWeighted = (TH1D*)hDGenPtTotal->Clone(Form("hDGenPtWeighted%s", label));
      for(int j=0; j<=i; j++)
	{
	  int binStart = ptCutBins[j];
	  int binEnd = ptCutBins[j+1];
	  if(j==i) binEnd = ptCutBins[7];
	  for(int k=binStart; k<binEnd; k++)
	    {
	      hDGenPtWeighted->SetBinContent(k, hDGenPtWeighted->GetBinContent(k)*sampleWeights[j]);
	      hDGenPtWeighted->SetBinError(k, hDGenPtWeighted->GetBinError(k)*sampleWeights[j]);
	    }
	}

    }

  hDGenPtWeighted->Draw();
  c1->SaveAs(Form("DGenPt%sSamplePtWeighted.pdf", label));

  TH1D* hDGenPtWeight = (TH1D*)hDGenPtWeighted->Clone(Form("hDPtSampleWeight%s", label));
  for(int i=0; i<7; i++)
    for(int j=ptCutBins[i]; j<ptCutBins[i+1]; j++)
      {
	hDGenPtWeight->SetBinContent(j, sampleWeights[i]);
	hDGenPtWeight->SetBinError(j, 0);
      }
  hDGenPtWeight->SetBinContent(0, sampleWeights[0]);
  hDGenPtWeight->SetBinError(0, 0);
  hDGenPtWeight->SetBinContent(ptCutBins[7], sampleWeights[6]);
  hDGenPtWeight->SetBinError(ptCutBins[7], 0);

  TFile* fOut = new TFile(Form("DPtSampleWeight%s.root", label), "recreate");
  hDGenPtWeighted->Write();
  hDGenPtWeight->Write();
  fOut->Close();
  
}

