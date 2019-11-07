// 1. read all files.. Easier start from a single merger file
// 2. calculate weight
// 3. create merge file
// 4. write to Merge file

#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/DsMinTreeLoad.h" 
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h" 
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h" 

#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"



#include <TString.h>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>
#include <TCut.h>

using namespace std;


int MC_MergeGenSamplePtWeight_D0(int isPbPb=0, int PNPrompt =0, TString MCFileList="", TString inMergeFile=""){


	InitStyle();

	TString str_isPbPb="pp";
	TString str_PNPrompt="Prompt";

	Int_t D0GenTrue=23333;

	if(isPbPb>0)	{
		str_isPbPb=Form("PbPb");
	}
	if(PNPrompt==1){
		str_PNPrompt="NonPrompt";
	}	

  // inMergeFile=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/%s_MC/DsMinTree_%s_MC_%s_%skkpi.root",str_isPbPb.Data(),str_isPbPb.Data(), str_PNPrompt.Data(),str_DsChannel.Data() );

	TFile *fout=new TFile(Form("./root_output_D0/D0_MC_GenSampleMerge_%s_%s.root",str_isPbPb.Data(),str_PNPrompt.Data()),"RECREATE");
	cout<<"write output file : "<<Form("./root_output_D0/D0_MC_GenSampleMerge_%s_%s.root",str_isPbPb.Data(),str_PNPrompt.Data())<<endl;

	Bool_t verbose=true;

	Bool_t REAL=false;

  TCut cutGenTrue=Form("(GisSignal==1 || GisSignal==2) && GcollisionId==0 && TMath::Abs(Gy)<1");
  TCut cutGenTrueGyAll=Form("(GisSignal==1 || GisSignal==2) && GcollisionId==0");
  TCut cutGenPNPrompt="GBAncestorpt<=0";
  if(PNPrompt==1){cutGenPNPrompt="GBAncestorpt>0";}

	// 1. read all files

	// TFile *f_in[nbin_pt_Gen];
	// TTree *t_gen[nbin_pt_Gen];
	// TTree *t_reco[nbin_pt_Gen];

  TChain *ch_in[nbin_pt_GenD0];
  for(int i=0; i<nbin_pt_GenD0 ; i++){
    ch_in[i]= new TChain("Dfinder/ntGen");
  }
	// prompt
	ch_in[0]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/PrmtD0_pThat-0_pT-0_pp_5p02-Pythia8/crab_Dfinder_ppMC_promptD0_pt0_1/170426_203300/0000/*.root");
	ch_in[1]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/PrmtD0_pThat-0_pT-2_pp_5p02-Pythia8/crab_Dfinder_ppMC_promptD0_pt2_1/170426_203421/0000/*.root");
	ch_in[2]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/PrmtD0_pThat-0_pT-4_pp_5p02-Pythia8/crab_Dfinder_ppMC_promptD0_pt4_1/170426_203538/0000/*.root");
	ch_in[3]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/PrmtD0_pThat-4_pT-10_pp_5p02-Pythia8/crab_Dfinder_ppMC_promptD0_pt10_1/170426_203644/0000/*.root");
	ch_in[4]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/PrmtD0_pThat-10_pT-20_pp_5p02-Pythia8/crab_Dfinder_ppMC_promptD0_pt20_1/170426_203755/0000/*.root");
	ch_in[5]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/PrmtD0_pThat-20_pT-40_pp_5p02-Pythia8/crab_Dfinder_ppMC_promptD0_pt40_1/170426_203924/0000/*.root");
	ch_in[6]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/PrmtD0_pThat-30_pT-60_pp_5p02-Pythia8/crab_Dfinder_ppMC_promptD0_pt60_1/170426_204135/0000/*.root");


	/* // nonprompt
    ch_in[0]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-0_pT-0_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt0_1/170426_204331/0000/*.root");
    ch_in[1]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-0_pT-2_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt2_1/170426_204802/0000/*.root");
    ch_in[2]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-0_pT-4_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt4_1/170426_204905/0000/*.root");
    ch_in[3]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-4_pT-10_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt10_1/170426_205009/0000/*.root");
    ch_in[4]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-10_pT-20_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt20_1/170426_205151/0000/*.root");
    ch_in[5]->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-20_pT-40_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt40_1/170426_205317/0000/*.root");
    ch_in[6]->Add(" /mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-30_pT-60_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt60_1/170426_205430/0000/*.root");
*/


	// fill each Gpt distribution 
  TH1D *h_Gpt_temp[nbin_pt_GenD0]; // Gpt for each sample
  TH1D *h_Gpt_temp_Gy1[nbin_pt_GenD0]; // Gpt for each sample
  TH1D *h_Gpt_temp_GyAll[nbin_pt_GenD0]; // Gpt for each sample

	cout<<"nbin_pt_GenD0 = "<<nbin_pt_GenD0<<endl;

	TCanvas *c_temp[nbin_pt_GenD0];

	// ch_in[0]->Draw("Gpt");

	for(int i =0; i<nbin_pt_GenD0; i++){
		cout<<"bins_pt_GenD0 = "<<bins_pt_GenD0[i]<<endl;
	}

  for(int i=0; i<nbin_pt_GenD0; i++){
		// ch_in[i]->Draw("Gpt",cutGenTrue&&cutGenPNPrompt);
		// return 1;
    h_Gpt_temp[i]=new TH1D(Form("h_Gpt_temp_%i",i),Form("h_Gpt_temp_%i",i),nbin_pt_GenD0,bins_pt_GenD0); h_Gpt_temp[i]->Sumw2();
    ch_in[i]->Project(Form("h_Gpt_temp_%i",i),"Gpt",cutGenTrue&&cutGenPNPrompt); // weight from pythia, account for bias2 option
    h_Gpt_temp_Gy1[i]=new TH1D(Form("h_Gpt_temp_Gy1_%i",i),Form("h_Gpt_temp_Gy1_%i",i),nbin_pt_Norm,bins_pt_Norm); h_Gpt_temp_Gy1[i]->Sumw2();
    ch_in[i]->Project(Form("h_Gpt_temp_Gy1_%i",i),"Gpt",cutGenTrue&&cutGenPNPrompt); // weight from pythia, account for bias2 option
    h_Gpt_temp_GyAll[i]=new TH1D(Form("h_Gpt_temp_GyAll_%i",i),Form("h_Gpt_temp_GyAll_%i",i),nbin_pt_Norm,bins_pt_Norm); h_Gpt_temp_GyAll[i]->Sumw2();
    ch_in[i]->Project(Form("h_Gpt_temp_GyAll_%i",i),"Gpt",cutGenTrueGyAll&&cutGenPNPrompt); // weight from pythia, account for bias2 option
		// ch_in[i]->Draw(Form("Gpt>>h_Gpt_temp_%i",i),cutGenTrue&&cutGenPNPrompt);	
	
		// c_temp[i]=new TCanvas(Form("c_temp_%i",i),Form("c_temp_%i",i),800,800);
		// c_temp[i]->cd();
		// h_Gpt_temp[i]->Draw();// checked ok	
		// break ;

  }

	// h_Gpt_temp[0]->Draw();
	// return 1;


	TH1D *h_GenSamplePtWeight=new TH1D("h_GenSamplePtWeight","h_GenSamplePtWeight",nbin_pt_GenD0,bins_pt_GenD0); h_GenSamplePtWeight->Sumw2();
	TH1D *h_GptWeighted_temp=(TH1D*)h_Gpt_temp[0]->Clone("h_GptWeighted_temp");
	TH1D *h_GptSumNext_temp;
	float GenSamplePtWeight[nbin_pt_GenD0]; // final weight
	GenSamplePtWeight[0]=1; // set first range pt weight to 1
	h_GenSamplePtWeight->SetBinContent(1,GenSamplePtWeight[0]);

	TCanvas *c_test[20];
	int c_testcount=0;

	TH1D *h_temp[nbin_pt_GenD0];

	for(int i=1; i<nbin_pt_GenD0; i++){
	  h_GptSumNext_temp=new TH1D("h_GptSumNext_temp","h_GptSumNext_temp",nbin_pt_GenD0,bins_pt_GenD0); h_GptSumNext_temp->Sumw2();
		for(int j=0; j<=i; j++){ // only compare & calculate for next sample
			h_GptSumNext_temp->Add(h_Gpt_temp[j]);	
			cout<<"i = "<<i<<" ,j = "<<j<<endl;	
		}
		// h_GptSumNext_temp->Draw();
		// GenSamplePtWeight[i]=h_GptWeighted_temp->Integral(h_GptWeighted_temp->GetXaxis()->FindBin(bins_pt_GenD0[i]+0.001),nbin_pt_GenD0+1)/h_GptSumNext_temp->Integral(h_GptWeighted_temp->GetXaxis()->FindBin(bins_pt_GenD0[i]+0.001),nbin_pt_GenD0+1) ;
		GenSamplePtWeight[i]=h_GptWeighted_temp->Integral(h_GptWeighted_temp->GetXaxis()->FindBin(bins_pt_GenD0[i]+0.001),h_GptWeighted_temp->GetXaxis()->FindBin(99))/h_GptSumNext_temp->Integral(h_GptSumNext_temp->GetXaxis()->FindBin(bins_pt_GenD0[i]+0.001),h_GptSumNext_temp->GetXaxis()->FindBin(99)) ;

		cout<<" \n\ni = "<<i<<endl;
		cout<<"h_GptWeighted_temp = "<<h_GptWeighted_temp->Integral(h_GptWeighted_temp->GetXaxis()->FindBin(bins_pt_GenD0[i]+0.001),nbin_pt_GenD0+1)<<" , h_GptSumNext_temp = "<<h_GptSumNext_temp->Integral(h_GptSumNext_temp->GetXaxis()->FindBin(bins_pt_GenD0[i]+0.001),nbin_pt_GenD0+1)<<endl;
		cout<<"GenSamplePtWeight[i] = "<<GenSamplePtWeight[i]<<endl;
		cout<<"binLow = "<<h_GptWeighted_temp->GetXaxis()->FindBin(bins_pt_GenD0[i]+0.001) << ", binHi = " <<nbin_pt_GenD0+1<<endl;
		cout<<"binLow Lowedge = "<<h_GptWeighted_temp->GetXaxis()->GetBinLowEdge( h_GptWeighted_temp->GetXaxis()->FindBin(bins_pt_GenD0[i]+0.001) )<<endl;
		h_GenSamplePtWeight->SetBinContent(i+1,GenSamplePtWeight[i]);

		c_test[c_testcount]=new TCanvas(Form("c_test%i",c_testcount), Form("c_test%i",c_testcount), 800,800);
		c_test[c_testcount]->cd();
		h_temp[i]=(TH1D*)h_GptWeighted_temp->Clone(Form("h_temp%i",i));
		h_temp[i]->Draw();
		c_testcount++;
		// c_test[c_testcount]->SaveAs(Form("./temp/h_GptWeighted_temp%i.png",c_testcount));


	  delete h_GptWeighted_temp;
		h_GptWeighted_temp=(TH1D*)h_GptSumNext_temp->Clone("h_GptWeighted_temp");
		for(int j=0; j<nbin_pt_GenD0; j++){
			if(j<=i){
				// TH1 first bin is 1 , not 0, need j+1
				h_GptWeighted_temp->SetBinContent(j+1,h_GptWeighted_temp->GetBinContent(j+1)*GenSamplePtWeight[j]);
			}else{  // apply weight to all after
				h_GptWeighted_temp->SetBinContent(j+1,h_GptWeighted_temp->GetBinContent(j+1)*GenSamplePtWeight[i]);
				// cout<<"jSampelweight j = "<< j<<endl;
			}
		} // end for j<nbin_pt_GenD0

		delete h_GptSumNext_temp;

	} // end for int i<nbin_pt_GenD0

	for(int i=0; i<nbin_pt_GenD0;i++){
		cout<<"GenSamplePtWeight i= "<<GenSamplePtWeight[i]<<endl;

	}



	// test weight alone
	// && build Gy

	fout->cd();

	TH1D *h_GptWeighted_Gy1=new TH1D("h_GptWeighted_Gy1","h_GptWeighted_Gy1",nbin_pt_Norm,bins_pt_Norm); h_GptWeighted_Gy1->Sumw2();
	TH1D *h_GptWeighted_GyAll=new TH1D("h_GptWeighted_GyAll","h_GptWeighted_GyAll",nbin_pt_Norm,bins_pt_Norm); h_GptWeighted_GyAll->Sumw2();
	for(int i=0; i<nbin_pt_GenD0; i++){
		h_GptWeighted_Gy1->Add(h_Gpt_temp_Gy1[i]);
		h_GptWeighted_GyAll->Add(h_Gpt_temp_GyAll[i]);
	}
	for(int i=0; i<nbin_pt_Norm; i++){
		h_GptWeighted_Gy1->SetBinContent(i+1, h_GptWeighted_Gy1->GetBinContent(i+1)*h_GenSamplePtWeight->GetBinContent(h_GenSamplePtWeight->GetXaxis()->FindBin(h_GptWeighted_Gy1->GetBinCenter(i+1))));
		h_GptWeighted_Gy1->SetBinError(i+1, h_GptWeighted_Gy1->GetBinError(i+1)*h_GenSamplePtWeight->GetBinContent(h_GenSamplePtWeight->GetXaxis()->FindBin(h_GptWeighted_Gy1->GetBinCenter(i+1))));

		h_GptWeighted_GyAll->SetBinContent(i+1, h_GptWeighted_GyAll->GetBinContent(i+1)*h_GenSamplePtWeight->GetBinContent(h_GenSamplePtWeight->GetXaxis()->FindBin(h_GptWeighted_Gy1->GetBinCenter(i+1))));
		h_GptWeighted_GyAll->SetBinError(i+1, h_GptWeighted_GyAll->GetBinError(i+1)*h_GenSamplePtWeight->GetBinContent(h_GenSamplePtWeight->GetXaxis()->FindBin(h_GptWeighted_Gy1->GetBinCenter(i+1))));
	}

	
	// GenSamplePtWeight result check
	h_Gpt_temp[0]->GetXaxis()->SetRangeUser(0,19);
	h_GptWeighted_temp->GetXaxis()->SetRangeUser(0,19);
	TCanvas *c_GenSamplePtWeight=Draw({h_Gpt_temp[0],h_GptWeighted_temp});
	c_GenSamplePtWeight->SaveAs(Form("./plots/D0GenSamplePtWeight_%s_%s.pdf",str_isPbPb.Data(),str_PNPrompt.Data()));
	TCanvas *c_GenSamplePtWeight_Compare=DrawCompare(h_Gpt_temp[0],h_GptWeighted_temp);
	c_GenSamplePtWeight_Compare->SaveAs(Form("./plots/D0GenSamplePtWeight_Compare_%s_%s.pdf",str_isPbPb.Data(),str_PNPrompt.Data()));
	
	// TCanvas *c_testWeightAlone=DrawCompare(h_GptWeighted_Gy1,h_GptWeighted_temp);
	// c_testWeightAlone->SaveAs(Form("./plots/D0GenSamplePtWeightTestAlone_Compare_%s_%s.pdf",str_isPbPb.Data(),str_PNPrompt.Data()));
	TCanvas *c_compareY1YAll=DrawCompare(h_GptWeighted_Gy1,h_GptWeighted_GyAll);
	c_compareY1YAll->SaveAs(Form("./plots/D0GenSamplePtGy1GyAll_Compare_%s_%s.pdf",str_isPbPb.Data(),str_PNPrompt.Data())); 

	h_GptWeighted_Gy1->Write();
	h_GptWeighted_GyAll->Write();


	return 1;
  // chekced ok

/*
	TH1D *h_GenSamplePtWeight_CheckRatio=(TH1D*)h_Gpt_temp[0]->Clone("h_GenSamplePtWeight_CheckRatio"); 
	h_GenSamplePtWeight_CheckRatio->Divide(h_GptWeighted_temp);
	// h_Gpt_temp[0]->SetAxisRange(0,19,"X");
	h_Gpt_temp[0]->GetXaxis()->SetRangeUser(0,40);
	h_GptWeighted_temp->GetXaxis()->SetRangeUser(0,40);
	h_GenSamplePtWeight_CheckRatio->GetXaxis()->SetRangeUser(0,40);
// h_GptWeighted_temp->SetAxisRange(0,19,"X");
	TCanvas *c_GptWeightedSample= new TCanvas("c_GptWeightedSample","c_GptWeightedSample",1200,600);
	c_GptWeightedSample->Divide(2,1);	
	c_GptWeightedSample->cd(1);
	h_GptWeighted_temp->Draw();
	c_GptWeightedSample->cd(2);
	h_GenSamplePtWeight_CheckRatio->Draw();
*/
 
	// 2. create merge file
/*
	// TFile *f_Merge=TFile::Open(inMergeFile.Data());
	TTree *t_gen_Merge=(TTree*)f_Merge->Get("ntGen");
	TTree *t_reco_Merge=(TTree*)f_Merge->Get("ntDs");
	DsMinTreeLoad DsTreeTest;
  DsTreeTest.SetGenBranch(t_gen_Merge,isPbPb); // an over lap issue with Gen & Reco tree
	DsTreeTest.SetBranch(t_reco_Merge,REAL,isPbPb);

  fout->cd();
	TTree *ntGen=t_gen_Merge->CloneTree(0);
	TTree *ntReco=t_reco_Merge->CloneTree(0);

	Float_t GptSampleWeight[MAX_XB];
	Float_t PbPbVzWeight;
	Float_t DgenptSampleWeight;

	// FONLL & D0 data weight
	Float_t GenFONLLWeight[MAX_XB];
	Float_t RecoFONLLWeight;


	Float_t GenFONLLRaaWeight[MAX_XB];
	Float_t RecoFONLLRaaWeight;

	Float_t GenD0DataWeight[MAX_XB];
	Float_t RecoD0DataWeight;





	TH1F *h_GptwithWeightFromTree=new TH1F("h_GptwithWeightFromTree","h_GptwithWeightFromTree",100,0,100); h_GptwithWeightFromTree->Sumw2();
	TH1F *h_GptnoWeightFromTree=new TH1F("h_GptnoWeightFromTree","h_GptnoWeightFromTree",100,0,100); h_GptnoWeightFromTree->Sumw2();


	TFile *f_FONLL=TFile::Open(Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_MergeWeight/FONLL/output/FONLLweightOverPythia_%s_%s_%s.root",str_isPbPb.Data(),str_PNPrompt.Data(),str_DsChannel.Data()),"READ");
	TH1D *h_Genpt_FONLLWeight;
	h_Genpt_FONLLWeight=(TH1D*)f_FONLL->Get("h_Genpt_FONLLWeight");
//	TH1D *h_Genpt_FONLLWeight_mt;
	if(PNPrompt==0){h_Genpt_FONLLWeight=(TH1D*)f_FONLL->Get("h_Genpt_FONLLWeight_mt"); } // Prompt use mt scaling weight


	TH1D *h_Genpt_FONLLRaaWeight=NULL;
	if(isPbPb){
	h_Genpt_FONLLRaaWeight=(TH1D*)f_FONLL->Get("h_Genpt_FONLLRaaWeight");
//	TH1D *h_Genpt_FONLLWeight_mt;
	if(PNPrompt==0){h_Genpt_FONLLRaaWeight=(TH1D*)f_FONLL->Get("h_Genpt_FONLLRaaWeight_mt"); } // Prompt use mt scaling weight
	}


	// TFile *f_DataWeight=TFile::Open( ,"READ");
	TF1 *f1_DataWeight=new TF1("f1_DataWeight","[0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*log(x)*log(x)*log(x)+ [4]*log(x)*log(x)*log(x)*log(x)");
	if(PNPrompt==0 && isPbPb==0 ){ f1_DataWeight->SetParameters(3.33929,-1.61349,0.438465,-0.0411468,0); }
	else if(PNPrompt==0 && isPbPb>0 ){ f1_DataWeight->SetParameters(4.21138,-4.0426,1.34612,-0.134656,0); }
	else if(PNPrompt==1 && isPbPb==0 ){ f1_DataWeight->SetParameters(0.55369,1.22995,-0.49766,0.05655,0); }
	else if(PNPrompt==1 && isPbPb>0 ){ f1_DataWeight->SetParameters(1.43138,0.37192,-0.71770,0.23269,-0.02181); }
	else { cout<<"wrong PNP or PbPb setting"<<endl;}

	TH1F *h_GptwithFONLLWeightFromTree=new TH1F("h_GptwithFONLLWeightFromTree","h_GptwithFONLLWeightFromTree",100,0,100); h_GptwithFONLLWeightFromTree->Sumw2();
	TH1F *h_GptwithFONLLRaaWeightFromTree=new TH1F("h_GptwithFONLLRaaWeightFromTree","h_GptwithFONLLRaaWeightFromTree",100,0,100); h_GptwithFONLLRaaWeightFromTree->Sumw2();
	TH1F *h_GptwithD0DataWeightFromTree=new TH1F("h_GptwithD0DataWeightFromTree","h_GptwithD0DataWeightFromTree",100,0,100); h_GptwithD0DataWeightFromTree->Sumw2();

	ntGen->Branch("GptSampleWeight",GptSampleWeight,"GptSampleWeight[Gsize]/F");
	ntReco->Branch("DgenptSampleWeight",&DgenptSampleWeight,"DgenptSampleWeight/F");

	ntGen->Branch("GenFONLLWeight",GenFONLLWeight,"GenFONLLWeight[Gsize]/F");
	ntReco->Branch("RecoFONLLWeight",&RecoFONLLWeight,"RecoFONLLWeight/F");

	ntGen->Branch("GenD0DataWeight",GenD0DataWeight,"GenD0DataWeight[Gsize]/F");
	ntReco->Branch("RecoD0DataWeight",&RecoD0DataWeight,"RecoD0DataWeight/F");

	if(isPbPb){
	ntGen->Branch("PbPbVzWeight",&PbPbVzWeight);
	ntReco->Branch("PbPbVzWeight",&PbPbVzWeight);

	ntGen->Branch("GenFONLLRaaWeight",GenFONLLRaaWeight,"GenFONLLRaaWeight[Gsize]/F");
	ntReco->Branch("RecoFONLLRaaWeight",&RecoFONLLRaaWeight,"RecoFONLLRaaWeight/F");
	}
	Long64_t totalN=0;
	Long64_t totalG=0;
	Long64_t totalG1=0;

	Long64_t nentries=0;

	// GenTree
	
	nentries=t_gen_Merge->GetEntries();
		for(int ien=0; ien<nentries; ien++){
			t_gen_Merge->GetEntry(ien);
			totalN++;

			if(isPbPb){
				PbPbVzWeight=f1_vzWeight->Eval(DsTreeTest.GenT_vz);
				// cout<<"DsTree[iGenSample].GenT_vz = "<<DsTree[iGenSample].GenT_vz<<" , PbPbVzWeight = "<<PbPbVzWeight<<endl;
			}

			for(int iGsize=0; iGsize<DsTreeTest.Gsize; iGsize++){
				GptSampleWeight[iGsize]=0;
				GptSampleWeight[iGsize]=h_GenSamplePtWeight->GetBinContent(h_GenSamplePtWeight->GetXaxis()->FindBin(DsTreeTest.Gpt[iGsize]));

				GenFONLLWeight[iGsize]=0;
				double Genpt=DsTreeTest.Gpt[iGsize];
				if(PNPrompt==1){ Genpt=DsTreeTest.GBAncestorpt[iGsize];  }
				GenFONLLWeight[iGsize]=h_Genpt_FONLLWeight->GetBinContent(h_Genpt_FONLLWeight->GetXaxis()->FindBin(Genpt) );

				if(isPbPb){
				GenFONLLRaaWeight[iGsize]=0;
				GenFONLLRaaWeight[iGsize]=h_Genpt_FONLLRaaWeight->GetBinContent(h_Genpt_FONLLRaaWeight->GetXaxis()->FindBin(Genpt) );
				}


				GenD0DataWeight[iGsize]=0;
				if( f1_DataWeight->Eval(Genpt) >=0){
				GenD0DataWeight[iGsize]=f1_DataWeight->Eval(Genpt);
				}

				// GptSampleWeight[iGsize]=h_GenSamplePtWeight->GetBinContent(h_GenSamplePtWeight->FindBin(DsTreeTest.Gpt[iGsize]));
				// cout<<"Gpt = "<<DsTree[iGenSample].Gpt[iGsize]<<" GptSampleWeight = "<<GptSampleWeight<<endl;
				// ntGen->Fill();	
				totalG++;

				if(DsTreeTest.GcollisionId[iGsize]==0 && DsTreeTest.GSignalType[iGsize]==GSignalTypeTrue && abs(DsTreeTest.Gy[iGsize])<=1 && ( (DsTreeTest.GBAncestorpt[iGsize]<=0 && PNPrompt==0 ) ||( DsTreeTest.GBAncestorpt[iGsize]>0 && PNPrompt==1  )  ) )
				{			
					// GptSampleWeight[iGsize]=h_GenSamplePtWeight->GetBinContent(h_GenSamplePtWeight->GetXaxis()->FindBin(DsTreeTest.Gpt[iGsize]));
					h_GptwithWeightFromTree->Fill(DsTreeTest.Gpt[iGsize],GptSampleWeight[iGsize]*DsTreeTest.GenT_weight);
					h_GptnoWeightFromTree->Fill(DsTreeTest.Gpt[iGsize]*DsTreeTest.GenT_weight);

					h_GptwithFONLLWeightFromTree->Fill(DsTreeTest.Gpt[iGsize],GptSampleWeight[iGsize]*DsTreeTest.GenT_weight*GenFONLLWeight[iGsize]);
					if(isPbPb){
					h_GptwithFONLLRaaWeightFromTree->Fill(DsTreeTest.Gpt[iGsize],GptSampleWeight[iGsize]*DsTreeTest.GenT_weight*GenFONLLRaaWeight[iGsize]);
					}
					h_GptwithD0DataWeightFromTree->Fill(DsTreeTest.Gpt[iGsize],GptSampleWeight[iGsize]*DsTreeTest.GenT_weight*GenFONLLWeight[iGsize]*GenD0DataWeight[iGsize]);

					totalG1++;
					}
				}
			ntGen->Fill();	
		}

	cout<<"totalN = "<<totalN<<endl;
	cout<<"totalG = "<<totalG<<endl;
	cout<<"totalG = "<<totalG1<<endl;


	TH1F *h_DgenptwithWeightFromTree=new TH1F("h_DgenptwithWeightFromTree","h_DgenptwithWeightFromTree",100,0,100); h_DgenptwithWeightFromTree->Sumw2();
	TH1F *h_DgenptnoWeightFromTree=new TH1F("h_DgenptnoWeightFromTree","h_DgenptnoWeightFromTree",100,0,100); h_DgenptnoWeightFromTree->Sumw2();



	TH1F *h_DgenptwithFONLLWeightFromTree=new TH1F("h_DgenptwithFONLLWeightFromTree","h_DgenptwithFONLLWeightFromTree",100,0,100); h_DgenptwithFONLLWeightFromTree->Sumw2();
	TH1F *h_DgenptwithFONLLRaaWeightFromTree=new TH1F("h_DgenptwithFONLLRaaWeightFromTree","h_DgenptwithFONLLRaaWeightFromTree",100,0,100); h_DgenptwithFONLLRaaWeightFromTree->Sumw2();
	TH1F *h_DgenptwithD0DataWeightFromTree=new TH1F("h_DgenptwithD0DataWeightFromTree","h_DgenptwithD0DataWeightFromTree",100,0,100); h_DgenptwithD0DataWeightFromTree->Sumw2();

	totalN=0;
	totalG=0;
	totalG1=0;

	// RecoTree // the simple tree structure is slow
	nentries=t_reco_Merge->GetEntries();
	for(int ien=0; ien<nentries; ien++){
		// cout<<"ien ="<<ien<<" nentries = "<<nentries<<endl;
		if(ien%100000==0) cout<<setw(7)<<ien<<" / "<<nentries<<endl;
		t_reco_Merge->GetEntry(ien);
		totalN++;
			if(isPbPb){
				PbPbVzWeight=f1_vzWeight->Eval(DsTreeTest.vz);
				// cout<<"DsTree[iGenSample].GenT_vz = "<<DsTree[iGenSample].GenT_vz<<" , PbPbVzWeight = "<<PbPbVzWeight<<endl;
			}

			DgenptSampleWeight=0;
			DgenptSampleWeight=h_GenSamplePtWeight->GetBinContent(h_GenSamplePtWeight->GetXaxis()->FindBin(DsTreeTest.Dgenpt));

			double Recopt=DsTreeTest.Dgenpt;
			if(PNPrompt==1) { Recopt=DsTreeTest.DgenBAncestorpt ; }

			RecoFONLLWeight=0;
			RecoFONLLWeight=h_Genpt_FONLLWeight->GetBinContent(h_Genpt_FONLLWeight->GetXaxis()->FindBin(Recopt));

			if(isPbPb){
			RecoFONLLRaaWeight=0;
			RecoFONLLRaaWeight=h_Genpt_FONLLRaaWeight->GetBinContent(h_Genpt_FONLLRaaWeight->GetXaxis()->FindBin(Recopt));
			}
				
			RecoD0DataWeight=0;
			RecoD0DataWeight=f1_DataWeight->Eval(Recopt);

			if(DsTreeTest.DgencollisionId==0 && DsTreeTest.DsGen==DsGenTrue &&( (DsTreeTest.DgenBAncestorpt<=0 && PNPrompt==0)||(DsTreeTest.DgenBAncestorpt>0 && PNPrompt==1) ) )  {
				totalG1++;
				h_DgenptwithWeightFromTree->Fill(DsTreeTest.Dgenpt,DgenptSampleWeight*DsTreeTest.weight);
				h_DgenptnoWeightFromTree->Fill(DsTreeTest.Dgenpt*DsTreeTest.weight);
				h_DgenptwithFONLLWeightFromTree->Fill(DsTreeTest.Dgenpt,DgenptSampleWeight*DsTreeTest.weight*RecoFONLLWeight);
				if(isPbPb){
				h_DgenptwithFONLLRaaWeightFromTree->Fill(DsTreeTest.Dgenpt,DgenptSampleWeight*DsTreeTest.weight*RecoFONLLRaaWeight);
				}

				h_DgenptwithD0DataWeightFromTree->Fill(DsTreeTest.Dgenpt,DgenptSampleWeight*DsTreeTest.weight*RecoFONLLWeight*RecoD0DataWeight);
			}
		
		ntReco->Fill();	
	} // end loop t_reco nentries

	cout<<"totalN = "<<totalN<<endl;
	cout<<"totalG = "<<totalG1<<endl;


	// cout<<"totalN = "<<totalN<<endl;
	// cout<<"totolG = "<<totalG<<endl;

	
	// reco tree
	
	// }


	TCanvas *c_result=new TCanvas("c_result","c_result",800,800);
	c_result->Divide(2,2);
	c_result->cd(1);
	gPad->SetLogx();
	gPad->SetLogy();
	h_GptwithWeightFromTree->Draw();
	c_result->cd(2);
	gPad->SetLogx();
	gPad->SetLogy();
	h_DgenptwithWeightFromTree->Draw();
	c_result->cd(3);
	gPad->SetLogx();
	gPad->SetLogy();
	h_GptwithFONLLWeightFromTree->Draw();
	c_result->cd(4);
	gPad->SetLogx();
	gPad->SetLogy();
	h_GptwithD0DataWeightFromTree->Draw();
	c_result->SaveAs(Form("plots/weightresultpt_%s_%s_%s.pdf",str_isPbPb.Data(),str_PNPrompt.Data(),str_DsChannel.Data()));


	fout->cd();
	h_GenSamplePtWeight->Write("",TObject::kOverwrite);
	ntGen->Write("",TObject::kOverwrite);
	h_GptWeighted_temp->Write("",TObject::kOverwrite);
	h_GptwithWeightFromTree->Write("",TObject::kOverwrite);
	h_GptnoWeightFromTree->Write("",TObject::kOverwrite);

	h_GptwithFONLLWeightFromTree->Write("",TObject::kOverwrite);
	h_GptwithD0DataWeightFromTree->Write("",TObject::kOverwrite);

	ntReco->Write("",TObject::kOverwrite);
  h_DgenptwithWeightFromTree->Write("",TObject::kOverwrite);
	h_DgenptnoWeightFromTree->Write("",TObject::kOverwrite);
  h_DgenptwithFONLLWeightFromTree->Write("",TObject::kOverwrite);
  h_DgenptwithD0DataWeightFromTree->Write("",TObject::kOverwrite);

	if(isPbPb){
	h_GptwithFONLLRaaWeightFromTree->Write("",TObject::kOverwrite);
  h_DgenptwithFONLLRaaWeightFromTree->Write("",TObject::kOverwrite);
	}
*/

	return 0; 

}


