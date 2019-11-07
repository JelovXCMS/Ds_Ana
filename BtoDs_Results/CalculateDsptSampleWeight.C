#include "../include/parameters.h"
#include "../include/DsMinTreeLoad.h"
//calculate Dpt Sample weight based on cross section info 

	// PNPrompt : 0 prompt, 1 : Nonprompt

void CalculateDsptSampleWeight(int PNPrompt=1){


//	double bins_pt_Gen[]={0,1.8,3.8,5.7,9.5,19,999};
//	const int nbin_pt_Gen=sizeof(bins_pt_Gen)/sizeof(bins_pt_Gen[0]) -1;

//	double bins_pt_Norm[]={0,2,3,4,6,8,10,20,30,999};
//	const int nbin_pt_Norm=sizeof(bins_pt_Norm)/sizeof(bins_pt_Norm[0])-1;
	TString PNP_name="Prompt";
  TString bins_pt_Name[nbin_pt_Gen]={"0","1p8","3p8","5p7","9p5","19"};
  double *CS_X_FilterEff;

	double CS_X_FilterEff_Prompt[nbin_pt_Gen]={1.18E+09,4.66E+08,8.77E+07,2.10E+07,2.66E+06,1.57E+05};
	double CS_X_FilterEff_Nonprompt[nbin_pt_Gen]={5.68E+07,3.06E+07,9.13E+06,3.09E+06,5.65E+05,3.77E+04};

	if(PNPrompt==0)
	{
		CS_X_FilterEff=CS_X_FilterEff_Prompt;
	}

	if(PNPrompt==1) 
	{
		PNP_name="NonPrompt";   
		CS_X_FilterEff=CS_X_FilterEff_Nonprompt;
	}

	TFile *f_in[nbin_pt_Gen];
	TTree *t_gen[nbin_pt_Gen];

	TFile *f_out=new TFile(Form("./output/DsptSampleWeight_pp_%s_phi.root",PNP_name.Data()),"RECREATE");

	TH1D *h_Gpt_temp[nbin_pt_Gen];
	TH1D *h_Gpt=new TH1D("h_Gpt","h_Gpt",nbin_pt_Gen,bins_pt_Gen);
	h_Gpt->Sumw2();

	TH1D *h_DsptSampleWeight=new TH1D("h_DsptSampleWeight","h_DsptSampleWeight",nbin_pt_Gen,bins_pt_Gen);
	h_DsptSampleWeight->Sumw2();

  int GSignalTypeTrue=1;

  TCut cutGenTrue=Form("GSignalType==%i && GcollisionId==0 && TMath::Abs(Gy)<1",GSignalTypeTrue);
  TCut cutGenPNPrompt="GBAncestorpt<=0";
  if(PNPrompt==1){cutGenPNPrompt="GBAncestorpt>0";}

	for(int i=0; i<nbin_pt_Gen; i++){
		f_in[i]=TFile::Open(Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_%s_phikkpi_pt%s.root",PNP_name.Data(),bins_pt_Name[i].Data()));
	  t_gen[i]=(TTree*)f_in[i]->Get("ntGen");
		h_Gpt_temp[i]=new TH1D(Form("h_Gpt_temp_%i",i),Form("h_Gpt_temp_%i",i),nbin_pt_Gen,bins_pt_Gen); h_Gpt_temp[i]->Sumw2();
		t_gen[i]->Project(Form("h_Gpt_temp_%i",i),"Gpt",cutGenTrue&&cutGenPNPrompt);
		h_Gpt->Add(h_Gpt_temp[i]);
	}

	TCanvas *c2= new TCanvas();
	c2->cd();
	h_Gpt->Draw();

	// method 1 : by cross section*filter eff
	for(int i=0; i<nbin_pt_Gen; i++){
		if(i<nbin_pt_Gen-1){
			double CS_X_FilterEff_1=CS_X_FilterEff[i+1];

			cout<<"i = "<<i<<" ,bin1 = "<<h_Gpt->FindBin(bins_pt_Gen[i]+0.001)<<" , bin2 = "<<h_Gpt->FindBin(bins_pt_Gen[i+1]-0.001)<<endl;
			cout<<"Integral = "<<h_Gpt->Integral(h_Gpt->FindBin(bins_pt_Gen[i]+0.001),h_Gpt->FindBin(bins_pt_Gen[i+1]-0.001))<<endl;
			h_DsptSampleWeight->SetBinContent(i+1, (CS_X_FilterEff[i]-CS_X_FilterEff[i+1])/ (h_Gpt->Integral(h_Gpt->FindBin(bins_pt_Gen[i]+0.001),h_Gpt->FindBin(bins_pt_Gen[i+1]-0.001) )) );
		}else{
			h_DsptSampleWeight->SetBinContent(i+1, (CS_X_FilterEff[i]) / (h_Gpt->Integral(h_Gpt->FindBin(bins_pt_Gen[i]+0.001),h_Gpt->FindBin(bins_pt_Gen[i+1]-0.001) )) ); 
		}
	}

	TCanvas *c3= new TCanvas();
	c3->cd();
	h_DsptSampleWeight->Draw();

	// method 2 : (Hao) by previous samples weighted Sum

// Method 2 : iterative by previous samples reweighted sum
  TH1D *h_GptWeighted=(TH1D*)h_Gpt_temp[0]->Clone("h_GptWeighted");
  h_GptWeighted->Scale(1.72e8/2e6);  // unknown constant... should be propotional to CS*FilterEff*Lumi
  float sampleWeights[nbin_pt_Gen];
  int ptCutBins[nbin_pt_Gen+1];
  for(int i=0; i<nbin_pt_Gen; i++){
    ptCutBins[i]=h_GptWeighted->GetXaxis()->FindBin(bins_pt_Gen[i]+0.00001);
  }
  ptCutBins[nbin_pt_Gen]=h_GptWeighted->GetXaxis()->GetNbins()+1;
  for(int i=0; i<nbin_pt_Gen; i++){
    TH1D *h_GptTotal=(TH1D*)h_Gpt_temp[0]->Clone("h_GptTotal");
    for(int j=1; j<=i; j++){
      h_GptTotal->Add(h_Gpt_temp[j]);
    }
    sampleWeights[i] = h_GptWeighted->Integral(ptCutBins[i], ptCutBins[nbin_pt_Gen]-1) / h_GptTotal->Integral(ptCutBins[i], ptCutBins[nbin_pt_Gen]-1);
    cout<<"sampleWeights[i]"<<sampleWeights[i]<<endl;

    delete h_GptWeighted;
    h_GptWeighted=(TH1D*)h_GptTotal->Clone("h_GptWeighted");
    for(int j=0; j<=i; j++){
      int binStart= ptCutBins[j];
      int binEnd = ptCutBins[j+1];
      if(j==i) binEnd=ptCutBins[nbin_pt_Gen];
      for(int k =binStart; k<binEnd; k++){
        h_GptWeighted->SetBinContent(k, h_GptWeighted->GetBinContent(k)*sampleWeights[j]);
        h_GptWeighted->SetBinError(k, h_GptWeighted->GetBinError(k)*sampleWeights[j]);
      }
    }
  } // end for i<nbin_pt_Gen
  TH1D *h_DsptSampleWeight2= new TH1D("h_DsptSampleWeight2","h_DsptSampleWeight2",nbin_pt_Gen,bins_pt_Gen);
  h_DsptSampleWeight2->Sumw2();
  for(int i=0; i<7; i++){
    for(int j=ptCutBins[i]; j<ptCutBins[i+1]; j++)
      {
  h_DsptSampleWeight2->SetBinContent(j, sampleWeights[i]);
  h_DsptSampleWeight2->SetBinError(j, 0);
      }
  }
  h_DsptSampleWeight2->SetBinContent(0, sampleWeights[0]);
  h_DsptSampleWeight2->SetBinError(0, 0);
  h_DsptSampleWeight2->SetBinContent(ptCutBins[nbin_pt_Gen], sampleWeights[nbin_pt_Gen-1]);
  h_DsptSampleWeight2->SetBinError(ptCutBins[nbin_pt_Gen], 0);

  for(int i=0; i<nbin_pt_Gen;i++){
    cout<<"h_DsptSampleWeight2 bin i = "<<i<<" weight= "<<h_DsptSampleWeight2->GetBinContent(i+1)<<endl;
  }
  TCanvas *c8= new TCanvas();
  c8->cd();
  h_GptWeighted->Draw();


  TCanvas *c9=new TCanvas();
  c9->cd();
  h_DsptSampleWeight2->Draw();
// end method 2


	//check with unbiased sample

  DsMinTreeLoad DsGenCls[nbin_pt_Gen];
	// TH1D *h_Gpt_weighted=new TH1D("h_Gpt_weighted","h_Gpt_weighted",nbin_pt_Gen,bins_pt_Gen); h_Gpt_weighted->Sumw2();
	TH1D *h_Gpt_weighted=new TH1D("h_Gpt_weighted","h_Gpt_weighted",nbin_pt_Norm,bins_pt_Norm); h_Gpt_weighted->Sumw2();
	TH1D *h_Gpt_weighted2=new TH1D("h_Gpt_weighted2","h_Gpt_weighted2",nbin_pt_Norm,bins_pt_Norm); h_Gpt_weighted2->Sumw2();
	TH1D *h_Gpt_weightedFineBin=new TH1D("h_Gpt_weightedFineBin","h_Gpt_weightedFineBin",100,0,50); h_Gpt_weightedFineBin->Sumw2();
	TH1D *h_Gpt_weighted2FineBin=new TH1D("h_Gpt_weighted2FineBin","h_Gpt_weighted2FineBin",100,0,50); h_Gpt_weighted2FineBin->Sumw2();
	for(int i=0; i<nbin_pt_Gen;i++){
		DsGenCls[i].SetGenBranch(t_gen[i]);
		Long64_t nentries = t_gen[i]->GetEntries();
		for(int j=0;j<nentries; j++){
			t_gen[i]->GetEntry(j);
			for(int k=0; k<DsGenCls[i].Gsize; k++){
				if(DsGenCls[i].GSignalType[k]==GSignalTypeTrue && DsGenCls[i].GcollisionId[k]==0 && TMath::Abs(DsGenCls[i].Gy[k])<1){
					if( (PNPrompt==1 && DsGenCls[i].GBAncestorpt[k]>0 )||( PNPrompt==0 && DsGenCls[i].GBAncestorpt[k]<=0 ) ){
					h_Gpt_weighted->Fill(DsGenCls[i].Gpt[k],h_DsptSampleWeight->GetBinContent(h_DsptSampleWeight->FindBin(DsGenCls[i].Gpt[k])));
					h_Gpt_weighted2->Fill(DsGenCls[i].Gpt[k],h_DsptSampleWeight2->GetBinContent(h_DsptSampleWeight2->FindBin(DsGenCls[i].Gpt[k])));
					h_Gpt_weightedFineBin->Fill(DsGenCls[i].Gpt[k],h_DsptSampleWeight->GetBinContent(h_DsptSampleWeight->FindBin(DsGenCls[i].Gpt[k])));
					h_Gpt_weighted2FineBin->Fill(DsGenCls[i].Gpt[k],h_DsptSampleWeight->GetBinContent(h_DsptSampleWeight->FindBin(DsGenCls[i].Gpt[k])));
					}
				}
			}// end fo Gsize
		}// end for nentries
	}// end for nbin_pt_Gen

/*	
	TH1D *h_Gpt_ratio=(TH1D*)h_Gpt_temp[0]->Clone("h_Gpt_ratio");
	h_Gpt_ratio->SetTitle("h_Gpt_ratio");
	h_Gpt_ratio->Divide(h_Gpt_weighted);

	TCanvas *c4=new TCanvas();
	c4->cd();
	h_Gpt_ratio->Draw();

	TCanvas *c5=new TCanvas();
	c5->cd();
	h_Gpt_weighted->Draw();
*/
	TH1D *h_Gpt_ratio_method12=(TH1D*)h_Gpt_weighted->Clone("h_Gpt_ratio_method12");
	h_Gpt_ratio_method12->SetTitle("h_Gpt_ratio_method12");
	h_Gpt_ratio_method12->Divide(h_Gpt_weighted2);
	TCanvas *c10=new TCanvas();
	c10->cd();
	h_Gpt_ratio_method12->Draw();

/*
	TH1D *h_Gpt_ratio2=(TH1D*)h_Gpt_temp[0]->Clone("h_Gpt_ratio2");
	h_Gpt_ratio2->SetTitle("h_Gpt_ratio2");
	h_Gpt_ratio2->Divide(h_Gpt_weighted2);
	TCanvas *c11=new TCanvas();
	c11->cd();
	h_Gpt_ratio2->Draw();	
*/	


/* bad method
	double N_DsRef=0;
	double N_DsAdd=0;
	double Weight[nbin_pt_Gen];
	Weight[0]=1;
	h_DsptSampleWeight->SetBinContent(1,Weight[0]);
	
	// calculate relative weight to previous sample	
	for(int i=0; i<nbin_pt_Gen-1; i++){
		N_DsRef=h_Gpt[i]->Integral(i+2,nbin_pt_Gen);  // i+1 first bin, i+2 to scale bin
		N_DsAdd=h_Gpt[i+1]->Integral(i+2,nbin_pt_Gen);
		cout<<"i = "<<i<<" , N_DsRef = "<<N_DsRef<<" , N_DsAdd = "<<N_DsAdd<<endl;
		Weight[i+1]=N_DsRef/(N_DsRef+N_DsAdd)*Weight[i];
		h_DsptSampleWeight->SetBinContent(i+2,Weight[i+1]);
	}
*/

	for(int i=0; i<nbin_pt_Gen; i++){
		cout<<"h_DsptSampleWeight , bin "<<i<<" , weight = "<<h_DsptSampleWeight->GetBinContent(i+1)<<endl;;
	}


	f_out->cd();
	h_DsptSampleWeight->Write("",TObject::kOverwrite);
	h_DsptSampleWeight2->Write("",TObject::kOverwrite);
	h_Gpt_weighted->Write("",TObject::kOverwrite);
	h_Gpt_weighted2->Write("",TObject::kOverwrite);
	h_Gpt_weightedFineBin->Write("",TObject::kOverwrite);
	h_Gpt_weighted2FineBin->Write("",TObject::kOverwrite);
//	h_Gpt_ratio->Write("",TObject::kOverwrite);
	// h_Gpt_ratio->Write("",TObject::kOverwrite);	
	// for(int i=0; i<nbin_pt_Gen; i++){
		// h_Gpt[i]->Write("",TObject::kOverwrite);
	// }




}
