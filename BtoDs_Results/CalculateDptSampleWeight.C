#include "../include/parameters.h"
#include "../include/setBranches.h"
// calcualte weight for D0 sample

//calculate 

void CalculateDptSampleWeight(){

	// double bins_pt_GenD0[]={0,1.8,3.8,5.7,9.5,19,999};
//	double bins_pt_GenD0[]={0,2,4,10,20,40,60,999};
//	const int nbin_pt_GenD0=sizeof(bins_pt_GenD0)/sizeof(bins_pt_GenD0[0]) -1;

//  double bins_pt_Norm[]={0,1.8,3.8,5.7,9.5,19,999};
//  const int nbin_pt_Norm=sizeof(bins_pt_Norm)/sizeof(bins_pt_Norm[0]) -1;

	// these number from George : https://twiki.cern.ch/twiki/bin/view/CMS/MC_for_2015_pp5TeVSub
/*	double CS[nbin_pt_GenD0]={1.01E+3,1.01E+3,1.01E+3,48.5,2.31,0.174,3.46E-2};
	double FilterEff[nbin_pt_GenD0]={0.00014,0.000055,0.000016,0.000021,0.000025,0.000019,0.000013};
	double CS_X_FilterEff[nbin_pt_GenD0];
	for(int i=0; i<nbin_pt_GenD0; i++){
		CS_X_FilterEff[i]=CS[i]*FilterEff[i];
	}
*/
	// this number from Hao https://twiki.cern.ch/twiki/bin/view/CMS/Run2015HIHFMCRequest#D_meson_new_request_Mar_2017
	double CS_X_FilterEff[nbin_pt_GenD0]={1.72E+8,6.16E+07,1.92E+07,9.21E+05,7.40E+04,3.66E+03,3.81E+02};


	TChain *ch_in[nbin_pt_GenD0];
	for(int i=0; i<nbin_pt_GenD0 ; i++){
		ch_in[i]= new TChain("Dfinder/ntGen");
	}
		ch_in[0]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-0_pT-0_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt0_1/170426_204331/0000/*.root");
		ch_in[1]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-0_pT-2_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt2_1/170426_204802/0000/*.root");
		ch_in[2]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-0_pT-4_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt4_1/170426_204905/0000/*.root");
		ch_in[3]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-4_pT-10_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt10_1/170426_205009/0000/*.root");
		ch_in[4]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-10_pT-20_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt20_1/170426_205151/0000/*.root");
		ch_in[5]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-20_pT-40_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt40_1/170426_205317/0000/*.root");
		ch_in[6]->Add(" /mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-30_pT-60_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt60_1/170426_205430/0000/*.root");

	TFile *f_out=new TFile("./output/DptSampleWeight_pp_NonPrompt.root","RECREATE");
	TH1D *h_Gpt_temp[nbin_pt_GenD0];
	TH1D *h_Gpt=new TH1D("h_Gpt","h_Gpt",nbin_pt_GenD0,bins_pt_GenD0); h_Gpt->Sumw2();

	TH1D *h_DptSampleWeight=new TH1D("h_DptSampleWeight","h_DptSampleWeight",nbin_pt_GenD0,bins_pt_GenD0);
	h_DptSampleWeight->Sumw2();

  // int GSignalTypeTrue=1;
  int PNPrompt=1;

  TCut cutGenTrue=Form("(GisSignal==1 || GisSignal==2) && GcollisionId==0 && TMath::Abs(Gy)<1");
  TCut cutGenPNPrompt="GBAncestorpt<=0";
  if(PNPrompt==1){cutGenPNPrompt="GBAncestorpt>0";}

// Method 1 : by cross section X filter effiency
	for(int i=0; i<nbin_pt_GenD0; i++){
		// f_in[i]=TFile::Open(Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_Prompt_phikkpi_pt%s.root",bins_pt_Name[i].Data()));
	  // t_gen[i]=(TTree*)f_in[i]->Get("ntGen");
		h_Gpt_temp[i]=new TH1D(Form("h_Gpt_temp_%i",i),Form("h_Gpt_temp_%i",i),nbin_pt_GenD0,bins_pt_GenD0); h_Gpt_temp[i]->Sumw2();
		ch_in[i]->Project(Form("h_Gpt_temp_%i",i),"Gpt",cutGenTrue&&cutGenPNPrompt);
		h_Gpt->Add(h_Gpt_temp[i]);
	}

  TCanvas *c2= new TCanvas();
  c2->cd();
  h_Gpt->Draw();

	for(int i=0; i<nbin_pt_GenD0; i++){
    if(i<nbin_pt_GenD0-1){
      double CS_X_FilterEff_1=CS_X_FilterEff[i+1];
      cout<<"i = "<<i<<" ,bin1 = "<<h_Gpt->FindBin(bins_pt_GenD0[i]+0.001)<<" , bin2 = "<<h_Gpt->FindBin(bins_pt_GenD0[i+1]-0.001)<<endl;
      cout<<"Integral = "<<h_Gpt->Integral(h_Gpt->FindBin(bins_pt_GenD0[i]+0.001),h_Gpt->FindBin(bins_pt_GenD0[i+1]-0.001))<<endl;
      h_DptSampleWeight->SetBinContent(i+1, (CS_X_FilterEff[i]-CS_X_FilterEff[i+1])/ (h_Gpt->Integral(h_Gpt->FindBin(bins_pt_GenD0[i]+0.001),h_Gpt->FindBin(bins_pt_GenD0[i+1]-0.001) )) );
    }else{
      h_DptSampleWeight->SetBinContent(i+1, (CS_X_FilterEff[i]) / (h_Gpt->Integral(h_Gpt->FindBin(bins_pt_GenD0[i]+0.001),h_Gpt->FindBin(bins_pt_GenD0[i+1]-0.001) )) );
    }
  }
//// end of Method 1

// Method 2 : iterative by previous samples reweighted sum
	TH1D *h_GptWeighted=(TH1D*)h_Gpt_temp[0]->Clone("h_GptWeighted");
	h_GptWeighted->Scale(1.72e8/2e6);
	float sampleWeights[nbin_pt_GenD0];
	int ptCutBins[nbin_pt_GenD0+1];
	for(int i=0; i<nbin_pt_GenD0; i++){
		ptCutBins[i]=h_GptWeighted->GetXaxis()->FindBin(bins_pt_GenD0[i]+0.00001);
	}
  ptCutBins[7]=h_GptWeighted->GetXaxis()->GetNbins()+1;
	for(int i=0; i<nbin_pt_GenD0; i++){
		TH1D *h_GptTotal=(TH1D*)h_Gpt_temp[0]->Clone("h_GptTotal");
		for(int j=1; j<=i; j++){
			h_GptTotal->Add(h_Gpt_temp[j]);
		}
		sampleWeights[i] = h_GptWeighted->Integral(ptCutBins[i], ptCutBins[7]-1) / h_GptTotal->Integral(ptCutBins[i], ptCutBins[7]-1);
	  cout<<"sampleWeights[i]"<<sampleWeights[i]<<endl;

		delete h_GptWeighted;
		h_GptWeighted=(TH1D*)h_GptTotal->Clone("h_GptWeighted");
		for(int j=0; j<=i; j++){
			int binStart= ptCutBins[j];
			int binEnd = ptCutBins[j+1];
			if(j==i) binEnd=ptCutBins[7];
			for(int k =binStart; k<binEnd; k++){
				h_GptWeighted->SetBinContent(k, h_GptWeighted->GetBinContent(k)*sampleWeights[j]);
				h_GptWeighted->SetBinError(k, h_GptWeighted->GetBinError(k)*sampleWeights[j]);
			}
		}
	} // end for i<nbin_pt_GenD0

	TH1D *h_DptSampleWeight2= new TH1D("h_DptSampleWeight2","h_DptSampleWeight2",nbin_pt_GenD0,bins_pt_GenD0);
	h_DptSampleWeight2->Sumw2();
  for(int i=0; i<7; i++){
    for(int j=ptCutBins[i]; j<ptCutBins[i+1]; j++)
      {
  h_DptSampleWeight2->SetBinContent(j, sampleWeights[i]);
  h_DptSampleWeight2->SetBinError(j, 0);
      }
	}
  h_DptSampleWeight2->SetBinContent(0, sampleWeights[0]);
  h_DptSampleWeight2->SetBinError(0, 0);
  h_DptSampleWeight2->SetBinContent(ptCutBins[7], sampleWeights[6]);
  h_DptSampleWeight2->SetBinError(ptCutBins[7], 0);

	for(int i=0; i<nbin_pt_GenD0;i++){
		cout<<"h_DptSampleWeight2 bin i = "<<i<<" weight= "<<h_DptSampleWeight2->GetBinContent(i+1)<<endl;
	}



	TCanvas *c8= new TCanvas();
	c8->cd();
	h_GptWeighted->Draw();


	TCanvas *c9=new TCanvas();
	c9->cd();
	h_DptSampleWeight2->Draw();


  TCanvas *c3= new TCanvas();
  c3->cd();
  h_DptSampleWeight->Draw();

	TH1D *h_Gpt_weighted=new TH1D("h_Gpt_weighted","h_Gpt_weighted",nbin_pt_Norm,bins_pt_Norm); h_Gpt_weighted->Sumw2();
	TH1D *h_Gpt_weighted2=new TH1D("h_Gpt_weighted2","h_Gpt_weighted2",nbin_pt_Norm,bins_pt_Norm); h_Gpt_weighted2->Sumw2();
	TH1D *h_Gpt_weightedFineBin=new TH1D("h_Gpt_weightedFineBin","h_Gpt_weightedFineBin",100,0,50); h_Gpt_weightedFineBin->Sumw2();
	TH1D *h_Gpt_weighted2FineBin=new TH1D("h_Gpt_weighted2FineBin","h_Gpt_weighted2FineBin",100,0,50); h_Gpt_weighted2FineBin->Sumw2();

	for(int i=0; i<nbin_pt_GenD0; i++){
		SetGenBranches(ch_in[i],true,false);
		Long64_t nentries=ch_in[i]->GetEntries();
		for(int j=0; j<nentries; j++){
			ch_in[i]->GetEntry(j);
			for(int k=0; k<Gsize; k++){
				if( (GisSignal[k]==1||GisSignal[k]==2) && GcollisionId[k]==0 && TMath::Abs(Dy[k])<1 ){
			  	if( (PNPrompt==1 && GBAncestorpt[k]>0 )||( PNPrompt==0 && GBAncestorpt[k]<=0 ) ){
	          h_Gpt_weighted->Fill(Gpt[k],h_DptSampleWeight->GetBinContent(h_DptSampleWeight->FindBin(Gpt[k])));				
	          h_Gpt_weighted2->Fill(Gpt[k],h_DptSampleWeight2->GetBinContent(h_DptSampleWeight2->FindBin(Gpt[k])));				
	          h_Gpt_weightedFineBin->Fill(Gpt[k],h_DptSampleWeight->GetBinContent(h_DptSampleWeight->FindBin(Gpt[k])));				
	          h_Gpt_weighted2FineBin->Fill(Gpt[k],h_DptSampleWeight2->GetBinContent(h_DptSampleWeight2->FindBin(Gpt[k])));				
					}
				}
			}// end for Gsize
		}// end for nentries
	}// end for nbin_pt_GenD0

  TCanvas *c5=new TCanvas();
  c5->cd();
  h_Gpt_weighted->Draw();

	

/*
  TH1D *h_Gpt_ratio=new TH1D("h_Gpt_ratio","h_Gpt_ratio",nbin_pt_Norm,bins_pt_Norm); h_Gpt_ratio->Sumw2();
	ch_in[0]->Project("h_Gpt_ratio","Gpt",cutGenTrue&&cutGenPNPrompt);
  h_Gpt_ratio->Divide(h_Gpt_weighted);

  TCanvas *c4=new TCanvas();
  c4->cd();
  h_Gpt_ratio->Draw();


  TH1D *h_Gpt_ratio2=new TH1D("h_Gpt_ratio2","h_Gpt_ratio2",nbin_pt_Norm,bins_pt_Norm); h_Gpt_ratio->Sumw2();
	ch_in[0]->Project("h_Gpt_ratio2","Gpt",cutGenTrue&&cutGenPNPrompt);
  h_Gpt_ratio2->Divide(h_Gpt_weighted2);

  TCanvas *c10=new TCanvas();
  c10->cd();
  h_Gpt_ratio2->Draw();
*/


	TFile *f_D0_weightbyHao=TFile::Open("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/MCWeightMerge/DNtuple_NonPromptPP_mergeWeight.root");
	TTree *t_D0Gen_WeightbyHao=(TTree*)f_D0_weightbyHao->Get("ntGen");
	TH1D *h_Gpt_fromHaoWeight= new TH1D("h_Gpt_fromHaoWeight","h_Gpt_fromHaoWeight",nbin_pt_Norm,bins_pt_Norm);
	TH1D *h_Gpt_fromHaoWeightFineBin= new TH1D("h_Gpt_fromHaoWeightFineBin","h_Gpt_fromHaoWeightFineBin",100,0,50);

	t_D0Gen_WeightbyHao->Project("h_Gpt_fromHaoWeight","Gpt",(cutGenTrue && cutGenPNPrompt)*"GptSampleWeight");	
	t_D0Gen_WeightbyHao->Project("h_Gpt_fromHaoWeightFineBin","Gpt",(cutGenTrue && cutGenPNPrompt)*"GptSampleWeight");	

  TH1D *h_Gpt_ratio_toHao=new TH1D("h_Gpt_ratio_toHao","h_Gpt_ratio_toHao",nbin_pt_Norm,bins_pt_Norm); h_Gpt_ratio_toHao->Sumw2();
	ch_in[0]->Project("h_Gpt_ratio_toHao","Gpt",cutGenTrue&&cutGenPNPrompt);
  h_Gpt_ratio_toHao->Divide(h_Gpt_fromHaoWeight);

  TCanvas *c6=new TCanvas();
  c6->cd();
	h_Gpt_ratio_toHao->Draw();


	TH1D *h_Gpt_HaoDefault_repeatRatio=(TH1D*)h_Gpt_fromHaoWeight->Clone("h_Gpt_HaoDefault_repeatRatio");
	h_Gpt_HaoDefault_repeatRatio->SetTitle("h_Gpt_HaoDefault_repeatRatio");
	h_Gpt_HaoDefault_repeatRatio->Divide(h_Gpt_weighted2);

	TCanvas *c11=new TCanvas();
	c11->cd();
	h_Gpt_HaoDefault_repeatRatio->Draw();



/* bad method
	double N_DsRef=0;
	double N_DsAdd=0;
	double Weight[nbin_pt_GenD0];
	Weight[0]=1;
	h_DptSampleWeight->SetBinContent(1,Weight[0]);
	
	// calculate relative weight to previous sample	
	for(int i=0; i<nbin_pt_GenD0-1; i++){
		N_DsRef=h_Gpt[i]->Integral(i+2,nbin_pt_GenD0);  // i+1 first bin, i+2 to scale bin
		N_DsAdd=h_Gpt[i+1]->Integral(i+2,nbin_pt_GenD0);
		cout<<"i = "<<i<<" , N_DsRef = "<<N_DsRef<<" , N_DsAdd = "<<N_DsAdd<<endl;
		Weight[i+1]=N_DsRef/(N_DsRef+N_DsAdd)*Weight[i];
		h_DptSampleWeight->SetBinContent(i+2,Weight[i+1]);
	}
*/
	for(int i=0; i<nbin_pt_GenD0; i++){

		cout<<"h_DptSampleWeight , bin "<<i<<" , weight = "<<h_DptSampleWeight->GetBinContent(i+1)<<endl;;

	}


	f_out->cd();
	h_DptSampleWeight->Write("",TObject::kOverwrite);
	h_DptSampleWeight2->Write("",TObject::kOverwrite);
	h_Gpt_weighted->Write("",TObject::kOverwrite);	
	h_Gpt_weighted2->Write("",TObject::kOverwrite);	
	h_Gpt_weightedFineBin->Write("",TObject::kOverwrite);
	h_Gpt_fromHaoWeight->Write("",TObject::kOverwrite);
	h_Gpt_fromHaoWeightFineBin->Write("",TObject::kOverwrite);
  // h_Gpt_ratio->Write("",TObject::kOverwrite);	
	//for(int i=0; i<nbin_pt_GenD0; i++){
	//	h_Gpt[i]->Write("",TObject::kOverwrite);
	//}


}

