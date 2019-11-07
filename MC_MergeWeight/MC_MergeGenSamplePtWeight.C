// 1. read all files.. Easier start from a single merger file
// 2. calculate weight
// 3. create merge file
// 4. write to Merge file

/*
// special handle for PbPb nornprompt MC weight, D0data weight is applied on FONLLRAA basis, so get new D0data weight=f1_data*FONLLRAAweight/FONLLWeight

h_Genpt_FONLLRaaWeight_mt is wrong , use h_Genpt_FONLLRaaWeight (no "_mt" ) for all
// use gr_D0Data_Weight for D0Data instead of f1_D0Data for prompt , mt_scaling is correct applied in gr but not in f1


*/





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


int MC_MergeGenSamplePtWeight(int isPbPb=0, int PNPrompt =0, int DsChannel=0, TString MCFileList="", TString inMergeFile=""){


	// InitStyle();

	TString str_isPbPb="pp";
	TString str_PNPrompt="Prompt";
	TString str_DsChannel="phi";

	Int_t DsGenTrue=23333;

	if(isPbPb>0)	{
		str_isPbPb=Form("PbPb");
	}
	if(PNPrompt==1){
		str_PNPrompt="NonPrompt";
	}	
	if(DsChannel==1){
		str_DsChannel="f0980";	
		DsGenTrue=24433;
	}
	if(DsChannel==2){
		str_DsChannel="kstar892";
		DsGenTrue=25544;
	}

	// MCFileList=Form("./MC_List/DsMinTree_PbPb_MC_Prompt_phikkpi.lis",);
	if(MCFileList==""){
	MCFileList=Form("./MC_List/DsMinTree_%s_MC_%s_%skkpi.lis",str_isPbPb.Data(),str_PNPrompt.Data(),str_DsChannel.Data() );}
	if(inMergeFile==""){
  // inMergeFile=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/%s_MC/DsMinTree_%s_MC_%s_%skkpi.root",str_isPbPb.Data(),str_isPbPb.Data(), str_PNPrompt.Data(),str_DsChannel.Data() ); }
  inMergeFile=Form("/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/%s_offMC/DsMinTree_%s_offMC_%s_%s.root",str_isPbPb.Data(),str_isPbPb.Data(), str_PNPrompt.Data(),str_DsChannel.Data() ); }

	if(!TFile::Open(inMergeFile.Data()))
	{
		cout<<" fail to open merge file input , exit"<<endl; 
		return 1;
	}

	TFile *f_Merge=TFile::Open(inMergeFile.Data());
	TFile *fout=new TFile(Form("/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/DsMinTree_MC_GenSampleMerge_%s_%s_%s.root",str_isPbPb.Data(),str_PNPrompt.Data(),str_DsChannel.Data()),"RECREATE");
	// TFile *fout=new TFile(Form("./root_output/DsMinTree_MC_GenSampleMerge_%s_%s_%s.root",str_isPbPb.Data(),str_PNPrompt.Data(),str_DsChannel.Data()),"RECREATE");
	// cout<<"write output file : "<<Form("./root_output/Ds_MC_GenSampleMerge_%s_%s_%s.root",str_isPbPb.Data(),str_PNPrompt.Data(),str_DsChannel.Data())<<endl;
	cout<<"write output file : "<<Form("/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/Ds_MC_GenSampleMerge_%s_%s_%s.root",str_isPbPb.Data(),str_PNPrompt.Data(),str_DsChannel.Data())<<endl;


	Bool_t verbose=true;

	initParameter();
	Bool_t REAL=false;

  int GSignalTypeTrue=1;
	if(DsChannel==1){GSignalTypeTrue=2;}; // f0980
	if(DsChannel==2){GSignalTypeTrue=3;}; // kstar892
  TCut cutGenTrue=Form("GSignalType==%i && GcollisionId==0 && TMath::Abs(Gy)<1",GSignalTypeTrue);
  TCut cutGenPNPrompt="GBAncestorpt<=0";
  if(PNPrompt==1){cutGenPNPrompt="GBAncestorpt>0";}


	// 1. read all files

	TFile *f_in[nbin_pt_Gen];
	TTree *t_gen[nbin_pt_Gen];
	TTree *t_reco[nbin_pt_Gen];
	DsMinTreeLoad DsTree[nbin_pt_Gen]; 

	TString MCFile_Array[nbin_pt_Gen]	;

	fstream inFileList(MCFileList.Data());
	if(!inFileList.is_open()){
		cout<<"MCFileList not exist"<<endl;
		return 1;
	}
	int ibin_pt_Gen=0;
	char buffer[300];
	while(inFileList>>MCFile_Array[ibin_pt_Gen] && ibin_pt_Gen<nbin_pt_Gen){ // somehow it would continue read at the end of file... 
	// while(!inFileList.eof()){
	// while(inFileList.good()){
		// inFileList.getline(buffer,sizeof(buffer));
		cout<<"reading MCFile i ="<<ibin_pt_Gen<<" from "<<MCFile_Array[ibin_pt_Gen].Data()<<endl;
		// cout<<"reading MCFile i ="<<ibin_pt_Gen<<" from "<<buffer<<endl;
		if(!TFile::Open(MCFile_Array[ibin_pt_Gen].Data()))
		{
			cout<<" fail to open MC input file, exit"<<endl; 
			return 2;
		}
		f_in[ibin_pt_Gen]=TFile::Open(MCFile_Array[ibin_pt_Gen].Data());
		// f_in[ibin_pt_Gen]=TFile::Open(buffer);
		t_gen[ibin_pt_Gen]=(TTree*)f_in[ibin_pt_Gen]->Get("ntGen");
		t_reco[ibin_pt_Gen]=(TTree*)f_in[ibin_pt_Gen]->Get("ntDs");
		DsTree[ibin_pt_Gen].SetGenBranch(t_gen[ibin_pt_Gen],isPbPb); // an over lap issue with Gen & Reco tree
		DsTree[ibin_pt_Gen].SetBranch(t_reco[ibin_pt_Gen],REAL,isPbPb); // an over lap issue with Gen & Reco tree
		// t_gen[ibin_pt_Gen]->Print();
		// t_reco[ibin_pt_Gen]->Print();
		ibin_pt_Gen++;
	 	cout<<"check i = "<<ibin_pt_Gen<<endl;
	}

	cout<<"Number of file load (ibin_pt_Gen, should be 6) = "<<ibin_pt_Gen<<endl;
	if(ibin_pt_Gen!=nbin_pt_Gen){
		cout<<"file number import incorrect , exit"<<endl;
		return 1;
	}

	// 2. Calculate weight, always apply pythia weight in any weight caculation

	// 2.1 vz weight function for PbPb MC , should be independet to other weight

	TF1 *f1_vzWeight = new TF1("f1_vzWeight","gaus(0)/(gaus(3))",-30.,30.);
  f1_vzWeight->SetParameters(0.08,0.44,5.12,0.08,3.25,5.23); 

	// 2.2 Calculate Dgenpt sample merge weight , by recursive using previous sample weight sum 

	// fill each Gpt distribution 
  TH1D *h_Gpt_temp[nbin_pt_Gen]; // Gpt for each sample

	cout<<"nbin_pt_Gen = "<<nbin_pt_Gen<<endl;

	TCanvas *c_temp[nbin_pt_Gen];

  for(int i=0; i<nbin_pt_Gen; i++){
    h_Gpt_temp[i]=new TH1D(Form("h_Gpt_temp_%i",i),Form("h_Gpt_temp_%i",i),nbin_pt_Gen,bins_pt_Gen); h_Gpt_temp[i]->Sumw2();
    t_gen[i]->Project(Form("h_Gpt_temp_%i",i),"Gpt",(TCut)(cutGenTrue&&cutGenPNPrompt)*"weight"); // weight from pythia, account for bias2 option
/*
		c_temp[i]=new TCanvas(Form("c_temp%i",i), Form("c_temp%i",i) , 800,800);
		c_temp[i]->cd();
	//	h_Gpt_temp[i]->Draw();
		t_gen[i]->Draw("Gpt",(TCut)(cutGenTrue&&cutGenPNPrompt)*"weight");
*/
  }



	TH1D *h_GenSamplePtWeight=new TH1D("h_GenSamplePtWeight","h_GenSamplePtWeight",nbin_pt_Gen,bins_pt_Gen); h_GenSamplePtWeight->Sumw2();
	TH1D *h_GptWeighted_temp=(TH1D*)h_Gpt_temp[0]->Clone("h_GptWeighted_temp");
	TH1D *h_GptSumNext_temp;
	float GenSamplePtWeight[nbin_pt_Gen]; // final weight
	GenSamplePtWeight[0]=1; // set first range pt weight to 1
	h_GenSamplePtWeight->SetBinContent(1,GenSamplePtWeight[0]);

	TCanvas *c_test[20];
	int c_testcount=0;

	TH1D *h_temp[nbin_pt_Gen];

	for(int i=1; i<nbin_pt_Gen; i++){
	  h_GptSumNext_temp=new TH1D("h_GptSumNext_temp","h_GptSumNext_temp",nbin_pt_Gen,bins_pt_Gen); h_GptSumNext_temp->Sumw2();
		for(int j=0; j<=i; j++){ // only compare & calculate for next sample
			h_GptSumNext_temp->Add(h_Gpt_temp[j]);	
			cout<<"i = "<<i<<" ,j = "<<j<<endl;	
		}
		// h_GptSumNext_temp->Draw();
		// GenSamplePtWeight[i]=h_GptWeighted_temp->Integral(h_GptWeighted_temp->GetXaxis()->FindBin(bins_pt_Gen[i]+0.001),nbin_pt_Gen+1)/h_GptSumNext_temp->Integral(h_GptWeighted_temp->GetXaxis()->FindBin(bins_pt_Gen[i]+0.001),nbin_pt_Gen+1) ;
		GenSamplePtWeight[i]=h_GptWeighted_temp->Integral(h_GptWeighted_temp->GetXaxis()->FindBin(bins_pt_Gen[i]+0.001),h_GptWeighted_temp->GetXaxis()->FindBin(99))/h_GptSumNext_temp->Integral(h_GptSumNext_temp->GetXaxis()->FindBin(bins_pt_Gen[i]+0.001),h_GptSumNext_temp->GetXaxis()->FindBin(99)) ;

		cout<<" \n\ni = "<<i<<endl;
		cout<<"h_GptWeighted_temp = "<<h_GptWeighted_temp->Integral(h_GptWeighted_temp->GetXaxis()->FindBin(bins_pt_Gen[i]+0.001),nbin_pt_Gen+1)<<" , h_GptSumNext_temp = "<<h_GptSumNext_temp->Integral(h_GptSumNext_temp->GetXaxis()->FindBin(bins_pt_Gen[i]+0.001),nbin_pt_Gen+1)<<endl;
		cout<<"GenSamplePtWeight[i] = "<<GenSamplePtWeight[i]<<endl;
		cout<<"binLow = "<<h_GptWeighted_temp->GetXaxis()->FindBin(bins_pt_Gen[i]+0.001) << ", binHi = " <<nbin_pt_Gen+1<<endl;
		cout<<"binLow Lowedge = "<<h_GptWeighted_temp->GetXaxis()->GetBinLowEdge( h_GptWeighted_temp->GetXaxis()->FindBin(bins_pt_Gen[i]+0.001) )<<endl;
		h_GenSamplePtWeight->SetBinContent(i+1,GenSamplePtWeight[i]);

		c_test[c_testcount]=new TCanvas(Form("c_test%i",c_testcount), Form("c_test%i",c_testcount), 800,800);
		c_test[c_testcount]->cd();
		h_temp[i]=(TH1D*)h_GptWeighted_temp->Clone(Form("h_temp%i",i));
		h_temp[i]->Draw();
		c_testcount++;
		// c_test[c_testcount]->SaveAs(Form("./temp/h_GptWeighted_temp%i.png",c_testcount));


	  delete h_GptWeighted_temp;
		h_GptWeighted_temp=(TH1D*)h_GptSumNext_temp->Clone("h_GptWeighted_temp");
		for(int j=0; j<nbin_pt_Gen; j++){
			if(j<=i){
				// TH1 first bin is 1 , not 0, need j+1
				h_GptWeighted_temp->SetBinContent(j+1,h_GptWeighted_temp->GetBinContent(j+1)*GenSamplePtWeight[j]);
			}else{  // apply weight to all after
				h_GptWeighted_temp->SetBinContent(j+1,h_GptWeighted_temp->GetBinContent(j+1)*GenSamplePtWeight[i]);
				// cout<<"jSampelweight j = "<< j<<endl;
			}
		} // end for j<nbin_pt_Gen

		delete h_GptSumNext_temp;

	} // end for int i<nbin_pt_Gen

	for(int i=0; i<nbin_pt_Gen;i++){
		cout<<"GenSamplePtWeight i= "<<GenSamplePtWeight[i]<<endl;

	}
	
	// GenSamplePtWeight result check
	h_Gpt_temp[0]->GetXaxis()->SetRangeUser(0,19);
	h_GptWeighted_temp->GetXaxis()->SetRangeUser(0,19);
	TCanvas *c_GenSamplePtWeight=Draw({h_Gpt_temp[0],h_GptWeighted_temp});
	c_GenSamplePtWeight->SaveAs(Form("./plots/GenSamplePtWeight_%s_%s_%s.pdf",str_isPbPb.Data(),str_PNPrompt.Data(),str_DsChannel.Data()));
	TCanvas *c_GenSamplePtWeight_Compare=DrawCompare(h_Gpt_temp[0],h_GptWeighted_temp);
	c_GenSamplePtWeight_Compare->SaveAs(Form("./plots/GenSamplePtWeight_Compare_%s_%s_%s.pdf",str_isPbPb.Data(),str_PNPrompt.Data(),str_DsChannel.Data()));
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

cout<<__LINE__<<endl;

	// TFile *f_Merge=TFile::Open(inMergeFile.Data());
	TTree *t_gen_Merge=(TTree*)f_Merge->Get("ntGen");
	TTree *t_reco_Merge=(TTree*)f_Merge->Get("ntDs");
	DsMinTreeLoad DsTreeTest;
cout<<__LINE__<<endl;

  DsTreeTest.SetGenBranch(t_gen_Merge,isPbPb); // an over lap issue with Gen & Reco tree
	DsTreeTest.SetBranch(t_reco_Merge,REAL,isPbPb,0);

cout<<__LINE__<<endl;

  fout->cd();
	TTree *ntGen=t_gen_Merge->CloneTree(0);
	TTree *ntReco=t_reco_Merge->CloneTree(0);

	Float_t GptSampleWeight[MAX_XB];
	Float_t PbPbVzWeight;
	Float_t DgenptSampleWeight;

	// FONLL & D0 data weight
	// Double_t GenFONLLWeight[MAX_XB];
	Float_t GenFONLLWeight[MAX_XB];
	Float_t RecoFONLLWeight;


	// Double_t GenFONLLRaaWeight[MAX_XB];
	Float_t RecoFONLLRaaWeight;
	Float_t GenFONLLRaaWeight[MAX_XB];
	// Float_t RecoFONLLRaaWeight;
	Float_t GenD0DataWeight[MAX_XB];
	// Double_t GenD0DataWeight[MAX_XB];
	Float_t GenDsDataWeight[MAX_XB];
	Float_t RecoD0DataWeight;
	Float_t RecoDsDataWeight;




/* no need for set binning
	double binFirst=0;
	double binLast=50;
	double binWidth=0.25;
	if(PNPrompt==1){binLast=100;}

  double current_bin=binFirst;
  vector<Double_t> v_bins;
  v_bins.push_back(0);
  while(current_bin<binLast){
    if(current_bin==10 || current_bin == 20  || current_bin == 40 || current_bin == 60)
    { binwidth*=2 ; }
    current_bin+=binwidth;
    v_bins.push_back(current_bin);
    cout<<"current_bin = "<<current_bin<<endl;
  }
  int nbin=v_bins.size();
  double bins[nbin];
  double bins_mtD0toDs[nbin-3];
  cout<<"nbin = "<<nbin<<endl;
  for(int i=0; i<nbin; i++){
      bins[i]=v_bins.at(i);
      cout<<"i = "<<i<<" , bins = "<< bins[i]<<endl;
      if(i>=3){
      bins_mtD0toDs[i-3]=sqrt(bins[i]*bins[i]+D0Mass*D0Mass-DsMass*DsMass);
      cout<<" bins_mtD0toDs = "<<bins_mtD0toDs[i-3]<<endl;
      }
  }
*/

cout<<__LINE__<<endl;

	TFile *f_FONLL=TFile::Open(Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_MergeWeight/FONLL/output/FONLLweightOverPythia_%s_%s_%s.root",str_isPbPb.Data(),str_PNPrompt.Data(),str_DsChannel.Data()),"READ");
	TH1D *h_Genpt_FONLLWeight;
	h_Genpt_FONLLWeight=(TH1D*)f_FONLL->Get("h_Genpt_FONLLWeight");
//	TH1D *h_Genpt_FONLLWeight_mt;
	// if(PNPrompt==0){h_Genpt_FONLLWeight=(TH1D*)f_FONLL->Get("h_Genpt_FONLLWeight_mt"); } // Prompt use mt scaling weight

cout<<__LINE__<<endl;

	TH1D *h_Genpt_FONLLRaaWeight=NULL;
	if(isPbPb){
	h_Genpt_FONLLRaaWeight=(TH1D*)f_FONLL->Get("h_Genpt_FONLLRaaWeight");
//	TH1D *h_Genpt_FONLLWeight_mt;
	// if(PNPrompt==0){h_Genpt_FONLLRaaWeight=(TH1D*)f_FONLL->Get("h_Genpt_FONLLRaaWeight_mt"); } // Prompt use mt scaling weight
	}

	TGraph *gr_D0Data_Weight=(TGraph*)f_FONLL->Get(Form("gr_D0Data_%s_weight",str_isPbPb.Data()));

cout<<__LINE__<<endl;

	// gr_D0Data_Weight->Draw();
	// return 1;
	

	// TFile *f_DataWeight=TFile::Open( ,"READ");
	TF1 *f1_DataWeight=new TF1("f1_DataWeight","[0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*log(x)*log(x)*log(x)+ [4]*log(x)*log(x)*log(x)*log(x)");
	TF1 *f1_DsDataWeight=new TF1("f1_DsDataWeight","([0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*log(x)*log(x)*log(x))*([4]+[5]*log(x)+[6]*log(x)*log(x)+[7]*log(x)*log(x)*log(x))*([8]+[9]*log(x)+[10]*log(x)*log(x))");
		f1_DsDataWeight->SetParameters(1,0,0,0,1,0,0,0,1,0,0);
	if(PNPrompt==0 && isPbPb==0 ){ 
		f1_DataWeight->SetParameters(3.33929,-1.61349,0.438465,-0.0411468,0); 
		f1_DsDataWeight->SetParameters(0.735902,0.0848251,-0.00260375,2.8312e-5,0.839458,0.0478619,0.00096107,1.96177e-06,0.899473,0.0268296,-5.44993e-05);
	}
	else if(PNPrompt==0 && isPbPb>0 ){ 
		f1_DataWeight->SetParameters(4.21138,-4.0426,1.34612,-0.134656,0); 
		f1_DsDataWeight->SetParameters(0.976754,0.0043164,0,0,0.978561,0.00406606,0,0,1,0,0);
	}
	else if(PNPrompt==1 && isPbPb==0 ){ f1_DataWeight->SetParameters(0.55369,1.22995,-0.49766,0.05655,0); }
	else if(PNPrompt==1 && isPbPb>0 ){ f1_DataWeight->SetParameters(1.43138,0.37192,-0.71770,0.23269,-0.02181); }
	else { cout<<"wrong PNP or PbPb setting"<<endl;}

	double binLow=2;
	double binHigh=50;
	int nbin=48;

	TH1D *h_GptwithWeightFromTree=new TH1D("h_GptwithWeightFromTree","h_GptwithWeightFromTree",nbin,binLow,binHigh); h_GptwithWeightFromTree->Sumw2();
	TH1D *h_GptnoWeightFromTree=new TH1D("h_GptnoWeightFromTree","h_GptnoWeightFromTree",nbin,binLow,binHigh); h_GptnoWeightFromTree->Sumw2();
	TH1D *h_GptwithFONLLWeightFromTree=new TH1D("h_GptwithFONLLWeightFromTree","h_GptwithFONLLWeightFromTree",nbin,binLow,binHigh); h_GptwithFONLLWeightFromTree->Sumw2();
	TH1D *h_GptwithFONLLRaaWeightFromTree=new TH1D("h_GptwithFONLLRaaWeightFromTree","h_GptwithFONLLRaaWeightFromTree",nbin,binLow,binHigh); h_GptwithFONLLRaaWeightFromTree->Sumw2();
	TH1D *h_GptwithD0DataWeightFromTree=new TH1D("h_GptwithD0DataWeightFromTree","h_GptwithD0DataWeightFromTree",nbin,binLow,binHigh); h_GptwithD0DataWeightFromTree->Sumw2();
	TH1D *h_GptwithDsDataWeightFromTree=new TH1D("h_GptwithDsDataWeightFromTree","h_GptwithDsDataWeightFromTree",nbin,binLow,binHigh); h_GptwithDsDataWeightFromTree->Sumw2();

cout<<__LINE__<<endl;

	ntGen->Branch("GptSampleWeight",GptSampleWeight,"GptSampleWeight[Gsize]/F");
	ntReco->Branch("DgenptSampleWeight",&DgenptSampleWeight,"DgenptSampleWeight/F");

	ntGen->Branch("GenFONLLWeight",GenFONLLWeight,"GenFONLLWeight[Gsize]/F");
	ntReco->Branch("RecoFONLLWeight",&RecoFONLLWeight,"RecoFONLLWeight/F");

	ntGen->Branch("GenD0DataWeight",GenD0DataWeight,"GenD0DataWeight[Gsize]/F");
	ntGen->Branch("GenDsDataWeight",GenDsDataWeight,"GenDsDataWeight[Gsize]/F");
	ntReco->Branch("RecoD0DataWeight",&RecoD0DataWeight,"RecoD0DataWeight/F");
	ntReco->Branch("RecoDsDataWeight",&RecoDsDataWeight,"RecoDsDataWeight/F");

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
	
cout<<__LINE__<<endl;

	nentries=t_gen_Merge->GetEntries();


		for(int ien=0; ien<nentries; ien++){
			t_gen_Merge->GetEntry(ien);
			totalN++;
			if(ien%500000==0) cout<<setw(7)<<ien<<" / "<<nentries<<endl;
			
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
				if(PNPrompt==0){

					// GenD0DataWeight[iGsize]=f1_DataWeight->Eval(Genpt);
					GenD0DataWeight[iGsize]=gr_D0Data_Weight->Eval(Genpt);

				}else if( f1_DataWeight->Eval(Genpt) >=0 && f1_DataWeight->Eval(Genpt) <=1000000000){
					GenD0DataWeight[iGsize]=f1_DataWeight->Eval(Genpt);
					if(isPbPb){
						Float_t RaaFactor=(GenFONLLRaaWeight[iGsize]/ GenFONLLWeight[iGsize]);
						if(RaaFactor>0.001 && RaaFactor<1000){
						GenD0DataWeight[iGsize]=f1_DataWeight->Eval(Genpt)*RaaFactor; // nonprompt PbPb D0 start from FONLL*Raa instead of FONLL, make formula consistent for latter use D0dataweighTotal= D0Dataweight*FONLL
						}
						if(GenD0DataWeight[iGsize] <=0 || GenD0DataWeight[iGsize] >= 10000){
							cout<<"GenD0DataWeight[iGsize] = "<<GenD0DataWeight[iGsize]<<endl;
							GenD0DataWeight[iGsize]=0.0001;
						}

					}
				}

				GenDsDataWeight[iGsize]=0;
				if( f1_DsDataWeight->Eval(Genpt) >=0 && f1_DsDataWeight->Eval(Genpt) <=1000000000){
				GenDsDataWeight[iGsize]=f1_DsDataWeight->Eval(Genpt);
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

					// h_GptwithD0DataWeightFromTree->Fill(DsTreeTest.Gpt[iGsize],GptSampleWeight[iGsize]*DsTreeTest.GenT_weight*GenFONLLWeight[iGsize]*GenD0DataWeight[iGsize]);
					h_GptwithD0DataWeightFromTree->Fill(DsTreeTest.Gpt[iGsize], GenD0DataWeight[iGsize]*GptSampleWeight[iGsize]*DsTreeTest.GenT_weight*GenFONLLWeight[iGsize]);
				// cout<<"GenD0DataWeight = "<<GenD0DataWeight[iGsize]<<" total = "<<GptSampleWeight[iGsize]*DsTreeTest.GenT_weight*GenFONLLWeight[iGsize]*GenD0DataWeight[iGsize]<<" before D0 = "<<GptSampleWeight[iGsize]*DsTreeTest.GenT_weight*GenFONLLWeight[iGsize]<<endl;

					h_GptwithDsDataWeightFromTree->Fill(DsTreeTest.Gpt[iGsize],GptSampleWeight[iGsize]*DsTreeTest.GenT_weight*GenFONLLWeight[iGsize]*GenD0DataWeight[iGsize]*GenDsDataWeight[iGsize]);

					totalG1++;
					} // if is signal
				}
			ntGen->Fill();	
		}

	cout<<"totalN = "<<totalN<<endl;
	cout<<"totalG = "<<totalG<<endl;
	cout<<"totalG = "<<totalG1<<endl;

cout<<__LINE__<<endl;

	TH1D *h_DgenptwithWeightFromTree=new TH1D("h_DgenptwithWeightFromTree","h_DgenptwithWeightFromTree",100,0,100); h_DgenptwithWeightFromTree->Sumw2();
	TH1D *h_DgenptnoWeightFromTree=new TH1D("h_DgenptnoWeightFromTree","h_DgenptnoWeightFromTree",100,0,100); h_DgenptnoWeightFromTree->Sumw2();



	TH1D *h_DgenptwithFONLLWeightFromTree=new TH1D("h_DgenptwithFONLLWeightFromTree","h_DgenptwithFONLLWeightFromTree",100,0,100); h_DgenptwithFONLLWeightFromTree->Sumw2();
	TH1D *h_DgenptwithFONLLRaaWeightFromTree=new TH1D("h_DgenptwithFONLLRaaWeightFromTree","h_DgenptwithFONLLRaaWeightFromTree",100,0,100); h_DgenptwithFONLLRaaWeightFromTree->Sumw2();
	TH1D *h_DgenptwithD0DataWeightFromTree=new TH1D("h_DgenptwithD0DataWeightFromTree","h_DgenptwithD0DataWeightFromTree",100,0,100); h_DgenptwithD0DataWeightFromTree->Sumw2();
	TH1D *h_DgenptwithDsDataWeightFromTree=new TH1D("h_DgenptwithDsDataWeightFromTree","h_DgenptwithDsDataWeightFromTree",100,0,100); h_DgenptwithDsDataWeightFromTree->Sumw2();

	totalN=0;
	totalG=0;
	totalG1=0;

cout<<__LINE__<<endl;

// adding DdlErr weight
	TFile *f_fit[4];
	TFile *f_fit2[4];
	TH1D *h_ratio[4];
	TH1D *h_ratio2[4];

	const int nbin_ratio=4;
	double bins_ratio[nbin_ratio+1]={2,4,6,10,40};

	double NorFactor[4];
  NorFactor[0]=1.78111/1.82891;
  NorFactor[1]=2.13961/2.17648;
  NorFactor[2]=1.31111/1.2707;
  NorFactor[3]=3.52482/3.45746;

	if(PNPrompt==1){
  NorFactor[0]=1.963222/2.00657;
  NorFactor[1]=2.96439/3.09694;
  NorFactor[2]=2.82868/2.78666;
  NorFactor[3]=1.08609/1.08427;
	}


cout<<__LINE__<<endl;

	TString Str_PbPb="pp";

/*
	for(int i =0; i<nbin_ratio; i++){
		f_fit[i]=TFile::Open(Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/ARC_HW/sPlot_varCompare_v2/fitout/%s_DdlErr_Dpt%.0fto%.0f.root",Str_PbPb.Data(),bins_ratio[i]*100,bins_ratio[i+1]*100));
    h_ratio[i]=(TH1D*)f_fit[i]->Get(Form("h_DdlErr_ratioSmooth"));
		f_fit2[i]=TFile::Open(Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/ARC_HW/sPlot_varCompare_v2/fitout_DdlErrWeight/%s_Ddl_Dpt%.0fto%.0f.root",Str_PbPb.Data(),bins_ratio[i]*100,bins_ratio[i+1]*100));
    h_ratio2[i]=(TH1D*)f_fit2[i]->Get(Form("h_Ddl_ratioSmooth"));
		h_ratio2[i]->Draw();
	}
*/
  double DdlErrMin=0.002;
  double DdlErrMax=0.1;
  double DdlMin=0.00;
  double DdlMax=1;
  double DptLowest=2;
  double DptHighest=40;

	Float_t DdlErrWeight;

	ntReco->Branch("DdlErrWeight",&DdlErrWeight,"DdlErrWeight/F");

cout<<__LINE__<<endl;
	// RecoTree // the simple tree structure is slow
	nentries=t_reco_Merge->GetEntries();
	for(int ien=0; ien<nentries; ien++){
		// cout<<"ien ="<<ien<<" nentries = "<<nentries<<endl;
		if(ien%100000==0) cout<<setw(7)<<ien<<" / "<<nentries<<endl;
		t_reco_Merge->GetEntry(ien);
		totalN++;
			
			DdlErrWeight=1;
/* //DdlErrWeight part
			if(DsTreeTest.Dpt>DptLowest && DsTreeTest.Dpt < DptHighest && DsTreeTest.DdlErr>DdlErrMin && DsTreeTest.DdlErr<DdlErrMax && DsTreeTest.Ddl>DdlMin && DsTreeTest.Ddl< DdlMax){
			for(int j=0; j<nbin_ratio; j++){
				if(DsTreeTest.Dpt>bins_ratio[j] && DsTreeTest.Dpt<bins_ratio[j+1]){
					DdlErrWeight=h_ratio[j]->GetBinContent(h_ratio[j]->FindBin(DsTreeTest.DdlErr))*NorFactor[j] * h_ratio2[j]->GetBinContent(h_ratio2[j]->FindBin(DsTreeTest.Ddl));
					if(DdlErrWeight<0.01){DdlErrWeight=0.01;}
					break;
				}
			}// end for j<nbin_ratio
			}// end if
*/

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
			RecoFONLLWeight*=1000;

			if(isPbPb){
			RecoFONLLRaaWeight=0;
			RecoFONLLRaaWeight=h_Genpt_FONLLRaaWeight->GetBinContent(h_Genpt_FONLLRaaWeight->GetXaxis()->FindBin(Recopt));
			RecoFONLLRaaWeight*=1000;

			}
				
			RecoD0DataWeight=0;
			if(PNPrompt==0){
					RecoD0DataWeight=gr_D0Data_Weight->Eval(Recopt);
				// RecoD0DataWeight=f1_DataWeight->Eval(Recopt); // pp nonprompt
			}else if( f1_DataWeight->Eval(Recopt) >=0 && f1_DataWeight->Eval(Recopt)<=100000000){
				if(isPbPb==0){
				RecoD0DataWeight=f1_DataWeight->Eval(Recopt); // pp nonprompt
				}else{

				// RecoD0DataWeight=f1_DataWeight->Eval(Recopt); // pp nonprompt
// #<{(|

				RecoD0DataWeight=f1_DataWeight->Eval(Recopt); // pp nonprompt

					Float_t RaaFactor=RecoFONLLRaaWeight/RecoFONLLWeight;
					if(RaaFactor>0.001&& RaaFactor<1000){
						RecoD0DataWeight*=RaaFactor;
					}else{
						// cout<<"RaaFactor = "<<RaaFactor<<endl;
					}
					// RecoD0DataWeight=RaaFactor*(Float_t)f1_DataWeight->Eval(Recopt); // PbPb nonprompt
					// cout<<"RecoD0DataWeight = "<<RecoD0DataWeight<<" , f1_DataWeight = "<<f1_DataWeight->Eval(Recopt)<<" ,RaaFactor = "<<RaaFactor<<endl;
					if(RecoD0DataWeight<0.00001 || RecoD0DataWeight >1000){
						RecoD0DataWeight=f1_DataWeight->Eval(Recopt);
					}
// |)}>#

						// GenD0DataWeight[iGsize]=f1_DataWeight->Eval(Genpt)*GenFONLLRaaWeight[iGsize]/ GenFONLLWeight[iGsize]; // nonprompt PbPb D0 start from FONLL*Raa instead of FONLL, make formula consistent for latter use D0dataweighTotal= D0Dataweight*FONLL
				}
			}



			RecoDsDataWeight=0;
        if( f1_DsDataWeight->Eval(Recopt) >=0 && f1_DsDataWeight->Eval(Recopt)<=100000000){
			RecoDsDataWeight=f1_DsDataWeight->Eval(Recopt);
			}


			if(DsTreeTest.DgencollisionId==0 && DsTreeTest.DsGen==DsGenTrue &&( (DsTreeTest.DgenBAncestorpt<=0 && PNPrompt==0)||(DsTreeTest.DgenBAncestorpt>0 && PNPrompt==1) ) )  {
				totalG1++;
				h_DgenptwithWeightFromTree->Fill(DsTreeTest.Dgenpt,DgenptSampleWeight*DsTreeTest.weight);
				h_DgenptnoWeightFromTree->Fill(DsTreeTest.Dgenpt*DsTreeTest.weight);
				h_DgenptwithFONLLWeightFromTree->Fill(DsTreeTest.Dgenpt,DgenptSampleWeight*DsTreeTest.weight*RecoFONLLWeight);
				if(isPbPb){
				h_DgenptwithFONLLRaaWeightFromTree->Fill(DsTreeTest.Dgenpt,DgenptSampleWeight*DsTreeTest.weight*RecoFONLLRaaWeight);
				}

				h_DgenptwithD0DataWeightFromTree->Fill(DsTreeTest.Dgenpt,DgenptSampleWeight*DsTreeTest.weight*RecoFONLLWeight*RecoD0DataWeight);
				h_DgenptwithDsDataWeightFromTree->Fill(DsTreeTest.Dgenpt,DgenptSampleWeight*DsTreeTest.weight*RecoFONLLWeight*RecoD0DataWeight*RecoDsDataWeight);

			}
		
		ntReco->Fill();	
	} // end loop t_reco nentries

	cout<<"totalN = "<<totalN<<endl;
	cout<<"totalG = "<<totalG1<<endl;


/*
	for(int iGenSample=0; iGenSample < nbin_pt_Gen; iGenSample++){

	// gen tree
		nentries=t_gen[iGenSample]->GetEntries();
		for(int ien=0; ien<nentries; ien++){
			t_gen[iGenSample]->GetEntry(ien);
			totalN++;

			if(isPbPb){
				PbPbVzWeight=f1_vzWeight->Eval(DsTree[iGenSample].GenT_vz);
				// cout<<"DsTree[iGenSample].GenT_vz = "<<DsTree[iGenSample].GenT_vz<<" , PbPbVzWeight = "<<PbPbVzWeight<<endl;
			}

			for(int iGsize=0; iGsize<DsTree[iGenSample].Gsize; iGsize++){
				GptSampleWeight[iGsize]=h_GenSamplePtWeight->GetBinContent(h_GenSamplePtWeight->FindBin(DsTree[iGenSample].Gpt[iGsize]));
				// cout<<"Gpt = "<<DsTree[iGenSample].Gpt[iGsize]<<" GptSampleWeight = "<<GptSampleWeight<<endl;
				// ntGen->Fill();	
				totalG++;
			}
			ntGen->Fill();	
		}
*/
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
	// h_GptwithDsDataWeightFromTree->Draw();
	c_result->SaveAs(Form("plots/weightresultpt_%s_%s_%s.pdf",str_isPbPb.Data(),str_PNPrompt.Data(),str_DsChannel.Data()));


	fout->cd();
	h_GenSamplePtWeight->Write("",TObject::kOverwrite);
	ntGen->Write("",TObject::kOverwrite);
	h_GptWeighted_temp->Write("",TObject::kOverwrite);
	h_GptwithWeightFromTree->Write("",TObject::kOverwrite);
	h_GptnoWeightFromTree->Write("",TObject::kOverwrite);

	h_GptwithFONLLWeightFromTree->Write("",TObject::kOverwrite);
	// h_GptwithD0DataWeightFromTree->Write("",TObject::kOverwrite);
	h_GptwithDsDataWeightFromTree->Write("",TObject::kOverwrite);

	ntReco->Write("",TObject::kOverwrite);
  h_DgenptwithWeightFromTree->Write("",TObject::kOverwrite);
	h_DgenptnoWeightFromTree->Write("",TObject::kOverwrite);
  h_DgenptwithFONLLWeightFromTree->Write("",TObject::kOverwrite);
  h_DgenptwithD0DataWeightFromTree->Write("",TObject::kOverwrite);
  h_DgenptwithDsDataWeightFromTree->Write("",TObject::kOverwrite);

	if(isPbPb){
	h_GptwithFONLLRaaWeightFromTree->Write("",TObject::kOverwrite);
  h_DgenptwithFONLLRaaWeightFromTree->Write("",TObject::kOverwrite);
	}

	h_GptwithD0DataWeightFromTree->Write("",TObject::kOverwrite);

	cout<<"MCFileList = "<<MCFileList<<endl;
	cout<<"inMergeFile = "<<inMergeFile<<endl;

	fout->Close();

	return 0; 

}


