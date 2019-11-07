#include "../include/DsMinTreeLoad.h"
#include "../include/parameters.h"
#include "../include/uti.h"

#include <TH1.h>
#include <TH1D.h>
#include <TFile.h>

// stupid


void DsPtShape(){


	TFile *f_DsptWeight=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/CalculateWeight/DsptSampleWeight_pp_NonPrompt_phi.root");
	TH1D *h_DsptSampleWeight=(TH1D*)f_DsptWeight->Get("h_DsptSampleWeight");

	TChain *ch= new TChain("ntGen");
	ch->Add("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_NonPrompt_phikkpi_pt*.root");

	DsMinTreeLoad chGen;
	chGen.SetGenBranch(ch);

	Long64_t nentries=ch->GetEntries();
	double DsptWeight=0;

	// double nbin_pt=nbin_pt_pp;
	// double *bins_pt=bins_pt_pp;

	double bins_pt[]={0,2,3,4,6,8,10,20,30,999}; // using b->D binning
	const int nbin_pt=sizeof(bins_pt)/sizeof(bins_pt[0])-1;

	TFile *f_out=new TFile("NonPrompt_DsGpt_phikkpi_pp.root","RECREATE");
	
	TH1D *h_DsGpt=new TH1D("h_DsGpt","h_DsGpt",nbin_pt,bins_pt); h_DsGpt->Sumw2();
	TH1D *h_DsGpt_fineBin=new TH1D("h_DsGpt_fineBin","h_DsGpt_fineBin",200,0,50); h_DsGpt_fineBin->Sumw2();

	int GSignalTypeTrue=1;


	for(int i=0; i<nentries; i++){

		ch->GetEntry(i);
		cout<<"Gsize = "<<chGen.Gsize<<endl;
		for(int j=0; j<chGen.Gsize; j++){
			DsptWeight=h_DsptSampleWeight->GetBinContent(h_DsptSampleWeight->FindBin(chGen.Gpt[j]));
			// cout<<"Weight = "<<DsptWeight<<endl;
			if(DsptWeight==0){
				cout<<"DsptWeight = 0 , Gpt = "<<chGen.Gpt[j]<<" , Gpt Bin = "<<(h_DsptSampleWeight->FindBin(chGen.Gpt[j]))<<endl;
			};
			if(chGen.GSignalType[j]==1 && chGen.GcollisionId[j]==0 && TMath::Abs(chGen.Gy[j])<1){
				if(chGen.GBAncestorpt[j]>0){
					h_DsGpt->Fill(chGen.Gpt[j],DsptWeight);
					h_DsGpt_fineBin->Fill(chGen.Gpt[j],DsptWeight);
				}
			}
		}
	}

	f_out->cd();

	cout<<"h_DsGpt->Integral() = "<<h_DsGpt->Integral()<<endl;
	cout<<"h_DsGpt->GetBinContent(1) = "<<h_DsGpt->GetBinContent(1)<<endl;
	cout<<"h_DsGpt->GetBinContent(2) = "<<h_DsGpt->GetBinContent(2)<<endl;
	cout<<"h_DsGpt->GetBinContent(3) = "<<h_DsGpt->GetBinContent(3)<<endl;

	h_DsGpt->Scale(1/h_DsGpt->Integral());

	h_DsGpt->Write("",TObject::kOverwrite);
	h_DsGpt_fineBin->Write("",TObject::kOverwrite);


	// D0 part


	TFile *f_DptWeight=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/CalculateWeight/DptSampleWeight_pp_NonPrompt_phi.root");
	TH1D *h_DptSampleWeight=(TH1D*)f_DptWeight->Get("h_DptSampleWeight");

	TChain *ch_D0Gen= new TChain("Dfinder/ntGen");
    ch_D0Gen->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-0_pT-0_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt0_1/170426_204331/0000/*.root");
    ch_D0Gen->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-0_pT-2_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt2_1/170426_204802/0000/*.root");
    ch_D0Gen->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-0_pT-4_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt4_1/170426_204905/0000/*.root");
    ch_D0Gen->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-4_pT-10_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt10_1/170426_205009/0000/*.root");
    ch_D0Gen->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-10_pT-20_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt20_1/170426_205151/0000/*.root");

	TH1D *h_DGpt=new TH1D("h_DGpt","h_DGpt",nbin_pt,bins_pt); h_DGpt->Sumw2();
	TH1D *h_DGpt_fineBin=new TH1D("h_DGpt_fineBin","h_DGpt_fineBin",200,0,50); h_DGpt_fineBin->Sumw2(); 

	int MaxGsize=20000;
	int Gsize;
	float Gpt[MaxGsize];
	int GisSignal[MaxGsize];
	int GcollisionId[MaxGsize];
	float Gy[MaxGsize];
	float GBAncestorpt[MaxGsize];

  ch_D0Gen->SetBranchAddress("Gsize", &Gsize);
  ch_D0Gen->SetBranchAddress("Gy", Gy);
  ch_D0Gen->SetBranchAddress("Gpt", Gpt);
  ch_D0Gen->SetBranchAddress("GcollisionId", GcollisionId);
  ch_D0Gen->SetBranchAddress("GisSignal", GisSignal);
  ch_D0Gen->SetBranchAddress("GBAncestorpt", GBAncestorpt);

	nentries=ch_D0Gen->GetEntries();

	double DptWeight=0;

  for(int i=0; i<nentries; i++){

    ch_D0Gen->GetEntry(i);
    // cout<<"Gsize = "<<chGen.Gsize<<endl;
    for(int j=0; j<Gsize; j++){
      DptWeight=h_DptSampleWeight->GetBinContent(h_DptSampleWeight->FindBin(chGen.Gpt[j]));
      // cout<<"Weight = "<<DptWeight<<endl;
      if(DptWeight==0){
        cout<<"DptWeight = 0 , Gpt = "<<chGen.Gpt[j]<<" , Gpt Bin = "<<(h_DptSampleWeight->FindBin(chGen.Gpt[j]))<<endl;
      };
      if( (GisSignal[j]==1||GisSignal[j]==2) && GcollisionId[j]==0 && TMath::Abs(Gy[j])<1){
        if(GBAncestorpt[j]>0){
          h_DGpt->Fill(chGen.Gpt[j],DptWeight);
          h_DGpt_fineBin->Fill(chGen.Gpt[j],DptWeight);
        }
      }
    }
  }

/*
	TFile *f_D0=TFile::Open("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/MCWeightMerge/DNtuple_NonPromptPP_mergeWeight.root");
	TTree *t_D0Gen=(TTree*)f_D0->Get("ntGen");

	TH1D *h_DGpt=new TH1D("h_DGpt","h_DGpt",nbin_pt,bins_pt); h_DGpt->Sumw2();
	TH1D *h_DGpt_fineBin=new TH1D("h_DGpt_fineBin","h_DGpt_fineBin",200,0,50); h_DGpt_fineBin->Sumw2(); 

	t_D0Gen->Project("h_DGpt","Gpt","((GisSignal==1||GisSignal==2) &&GcollisionId==0 && abs(Gy)<1 && GBAncestorpt>0)*GptSampleWeight");
	t_D0Gen->Project("h_DGpt_fineBin","Gpt","((GisSignal==1||GisSignal==2) &&GcollisionId==0 && abs(Gy)<1 && GBAncestorpt>0)*GptSampleWeight");
*/

	h_DGpt->Scale(1/h_DGpt->Integral());
	f_out->cd();
	h_DGpt->Write("",TObject::kOverwrite);
	h_DGpt_fineBin->Write("",TObject::kOverwrite);

	// h_DsGpt->Scale(1/h_DsGpt->Integral());
	// h_DGpt->Scale(1/h_DGpt->Integral());

	TH1D *h_DsD0GptWeight=(TH1D*)h_DsGpt->Clone("h_DsD0GptWeight");
	h_DsD0GptWeight->Divide(h_DGpt);

	h_DsD0GptWeight->Write("",TObject::kOverwrite);


	// for test weight

	TH1D *h_DsGpt_fromMB=new TH1D("h_DsGpt_fromMB","h_DsGpt_fromMB",nbin_pt,bins_pt);
	TChain *ch_DsMB= new TChain("ntGen");
	ch_DsMB->Add("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_NonPrompt_phikkpi_pt0.root");
	ch_DsMB->Project("h_DsGpt_fromMB","Gpt","(GSignalType==1 &&GcollisionId==0 && abs(Gy)<1 && GBAncestorpt>0)");

	TH1D *h_DsGpt_MBWeightRatio=(TH1D*)h_DsGpt_fromMB->Clone("h_DsGpt_MBWeightRatio");
	h_DsGpt_MBWeightRatio->Divide(h_DsGpt);

	h_DsGpt_fromMB->Write("",TObject::kOverwrite);
	h_DsGpt_MBWeightRatio->Write("",TObject::kOverwrite);



	TH1D *h_DGpt_fromMB=new TH1D("h_DGpt_fromMB","h_DGpt_fromMB",nbin_pt,bins_pt);
	TChain *ch_DMB= new TChain("Dfinder/ntGen");
	ch_DMB->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-0_pT-0_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt0_1/170426_204331/0000/*.root");
	ch_DMB->Project("h_DGpt_fromMB","Gpt","( (GisSignal==1 || GisSignal==2) &&GcollisionId==0 && abs(Gy)<1 && GBAncestorpt>0)");

	TH1D *h_DGpt_MBWeightRatio=(TH1D*)h_DGpt_fromMB->Clone("h_DGpt_MBWeightRatio");
	h_DGpt_MBWeightRatio->Divide(h_DGpt);

	h_DGpt_fromMB->Write("",TObject::kOverwrite);
	h_DGpt_MBWeightRatio->Write("",TObject::kOverwrite);

	// so stupid just use first sample ratio




 


}
