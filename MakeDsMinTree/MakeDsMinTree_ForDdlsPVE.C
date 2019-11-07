using namespace std;
#include "../include/setBranches.h"
#include <memory>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <list>
#include <iomanip>
#include <cmath>

#include <TLorentzVector.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include <TString.h>
#include <TNtuple.h>
#include <TVector3.h>
#include <TRandom3.h>

// adding is PbPb30-100 

Bool_t istest = false;
bool fillZeroCandEvt = true;
Bool_t detailMode=false;
int MakeDsMinTree_ForDdlsPVE(TString infile="", TString outfile="", Bool_t REAL=false, Int_t isPbPb=1, Int_t startEntries=0, Int_t endEntries=-1, Bool_t skim=true, Bool_t gskim=true, Bool_t checkMatching=true, Bool_t iseos=false, Bool_t SkimHLTtree=false)
{
	// isPbPb=1 for centrality 30-100 trigger only, isPbPb=2 for GJ and trackonly 0-100

	TRandom3 Rdm;
	Rdm.SetSeed(0);
	

	if(istest)
	{
		//    infile="/home/peng43/work/Project/Ds_PbPb/CMSSW/DsFinder/TestSample/Dsfinder_pp_mc_detail_phikkpi_n128.root";
		//		infile="/mnt/hadoop/store/user/chengchi/Dsfinder_f0980kkpi_18220/pp_MC/Ds_f0980kkpi_pp/MC_pp_Ds_f0980kkpi_pt4/180225_150245/0000/Dsfinder_pp_mc_1.root";
		//    infile="/mnt/hadoop/store/user/chengchi/Ds_phikkpi_18220/pp_MC/Ds_phikkpi_pp/MC_pp_Ds_phikkpi_pt4/180225_003140/0000/finder_pp_mc_1.root";
	  infile="/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_Data/MB1/Dntuple_finder_pp_91.root_414.root";
		outfile="test_DsMinTree.root";
		REAL=false;
		isPbPb=false;
		skim=false;
		checkMatching=true;
		iseos=false;

		cout<<"in test mode"<<endl;

	}
	cout<<endl;
	if(REAL) cout<<"--- Processing - REAL DATA";
	else {cout<<"--- Processing - MC ";
		// skim=false;
	}
	if(isPbPb) cout<<" - PbPb";
	else cout<<" - pp";
	cout<<"skim = "<<skim<<endl;
	TString ifname;
	if(iseos) ifname = Form("root://eoscms.cern.ch//eos/cms%s",infile.Data());
	else ifname = infile;
	cout<<"ifname = "<<ifname<<endl;
	//  if(FfromMnt) ifname= Form("root://xrootd.rcac.purdue.edu//store",infile.Data());
	if (!TFile::Open(ifname))   { cout << " fail to open file" << endl; return 0;}

	TFile* fin = TFile::Open(ifname);
	TTree* t_Hlt = (TTree*)fin->Get("ntHlt"); // setbranch not built up yet
	TTree* t_Skim = (TTree*)fin->Get("ntSkim");
	TTree* t_Hi = (TTree*)fin->Get("ntHi");
	TTree* t_Ds =(TTree*)fin->Get("ntDPhikkpi");


  TTree *t_Gen= (TTree*)fin->Get("ntGen");
  SetGenBranches(t_Gen,true); // true for select only


	SetRecoBranches(t_Ds, REAL, detailMode);
	SetHltBranches(t_Hlt, isPbPb); // not build up yet
	SetSkimBranches(t_Skim, isPbPb);
	SetHiBranches(t_Hi,REAL, isPbPb);

	Long64_t nentries = t_Ds->GetEntries();
	if(endEntries>nentries || endEntries == -1) endEntries = nentries;
	TFile* fout = TFile::Open(Form("%s", outfile.Data()),"recreate");
	cout<<"--- Building trees"<<endl;

	fout->cd();
	TTree* nt_Ds = new TTree("ntDs","");

  TTree *ntGen=t_Gen->CloneTree(0);

	int nbin_Ddls=990;
	double binLow_Ddls=1;
	double binHigh_Ddls=100;

	TH1D *h_Ddls=new TH1D("h_Ddls","Decay Length Significance",nbin_Ddls,binLow_Ddls,binHigh_Ddls); h_Ddls->Sumw2();
	TH1D *h_Ddls_Gaus=new TH1D("h_Ddls_Gaus","Decay Length Significance",nbin_Ddls,binLow_Ddls,binHigh_Ddls); h_Ddls_Gaus->Sumw2();


	// add for PVerror study
	Float_t         PVxE;
	Float_t         PVyE;
	Float_t         PVzE;
	int             Trksize;

	t_Ds->SetBranchAddress("Trksize",&Trksize);
	t_Ds->SetBranchAddress("PVxE",&PVxE);
	t_Ds->SetBranchAddress("PVyE",&PVyE);
	t_Ds->SetBranchAddress("PVzE",&PVzE);


	// declare needed variables 
	//	const double MAX_XB=20000;
	Float_t         dmass;
	Float_t         dpt;
	Float_t         dchi2cl;  
	Float_t         dalpha;  
	Float_t         ddls;   
	Float_t         ddlsPPVE;   
	Float_t         ddlsMPVE;   
	Float_t         ddca;   
	Float_t         dtrk1Pt; 
	Float_t         dtrk2Pt; 
	Float_t         dtrk3Pt; 
	// Float_t         dtrk1PtErr; 
	 // Float_t         dtrk2PtErr; 
	// Float_t         dtrk3PtErr; 
	Float_t         dtktkResmass; 
	Float_t         dtktkRes_chi2cl; 
	Float_t         dtktkRes_svpvDistanceToSV;   
	Float_t         dtktkRes_dlsToSV;  

	Int_t           dhiBin;
	Float_t         dpthat;
	Float_t         dweight;
	Float_t         dEvtvx;
	Float_t         dEvtvy;
	Float_t         dEvtvz;
	Float_t         dNpart;
	Float_t         dNcoll;


	// Int_t           dgen;   
	Int_t           dsGen;  
	Float_t         dgenpt; 
	Int_t           dgencollisionId;  
	Float_t         dgenBAncestorpt;  
	Int_t           dgenBAncestorpdgId; 
	Int_t           dgenfromgenPV;                 



	// branches
	nt_Ds->Branch("hiBin",&dhiBin);
	nt_Ds->Branch("Dmass",                      &dmass                    );    
	nt_Ds->Branch("Dpt",                        &dpt                      );   
	nt_Ds->Branch("Dchi2cl",                    &dchi2cl                  );   
	nt_Ds->Branch("Dalpha",                     &dalpha                   );   
	nt_Ds->Branch("Ddls",                       &ddls                     );   
	nt_Ds->Branch("DdlsPPVE",                   &ddlsPPVE                 );   
	nt_Ds->Branch("DdlsMPVE",                   &ddlsMPVE                 );   
	nt_Ds->Branch("Ddca",                       &ddca                     );   
	nt_Ds->Branch("Dtrk1Pt",                    &dtrk1Pt                  );   
	nt_Ds->Branch("Dtrk2Pt",                    &dtrk2Pt                  );   
	nt_Ds->Branch("Dtrk3Pt",                    &dtrk3Pt                  );   
	// nt_Ds->Branch("Dtrk1PtErr",                 &dtrk1PtErr               );   
	// nt_Ds->Branch("Dtrk2PtErr",                 &dtrk2PtErr               );   
	// nt_Ds->Branch("Dtrk3PtErr",                 &dtrk3PtErr               );   
	nt_Ds->Branch("DtktkResmass",               &dtktkResmass             );   
	nt_Ds->Branch("DtktkRes_chi2cl",            &dtktkRes_chi2cl          );   
	nt_Ds->Branch("DtktkRes_svpvDistanceToSV",  &dtktkRes_svpvDistanceToSV);   
	nt_Ds->Branch("DtktkRes_dlsToSV",    &dtktkRes_dlsToSV  );   

  nt_Ds->Branch("DsGen",&dsGen); // add for using TMVA 
	nt_Ds->Branch("DgencollisionId",&dgencollisionId);
	nt_Ds->Branch("DgenBAncestorpt",&dgenBAncestorpt);

	nt_Ds->Branch("vx",&dEvtvx);
	nt_Ds->Branch("vy",&dEvtvy);
	nt_Ds->Branch("vz",&dEvtvz);


	if(!REAL){
		// Dgen info
		nt_Ds->Branch("pthat",&dpthat);
	  nt_Ds->Branch("weight",&dweight);
//		nt_Ds->Branch("DsGen",&dsGen);
		nt_Ds->Branch("Dgenpt",&dgenpt);
		// nt_Ds->Branch("DgencollisionId",&dgencollisionId);
//		nt_Ds->Branch("DgenBAncestorpt",&dgenBAncestorpt);
		nt_Ds->Branch("DgenBAncestorpdgId",&dgenBAncestorpdgId);
		nt_Ds->Branch("DgenfromgenPV",&dgenfromgenPV);


		ntGen->Branch("pthat",&dpthat);
		ntGen->Branch("hiBin",&dhiBin);
		ntGen->Branch("weight",&dweight);

	  ntGen->Branch("vx",&dEvtvx);
	  ntGen->Branch("vy",&dEvtvy);
	  ntGen->Branch("vz",&dEvtvz);


	if(isPbPb){
		nt_Ds->Branch("Npart",&dNpart);
		nt_Ds->Branch("Ncoll",&dNcoll);

		ntGen->Branch("Npart",&dNpart);
		ntGen->Branch("Ncoll",&dNcoll);

		}

	}


	/*  full tree structure
	// declare needed variables 
	const double MAX_XB=20000;
	Int_t           size;
	Int_t           dhiBin;
	Float_t         dmass[MAX_XB];  
	Float_t         dpt[MAX_XB];   
	Float_t         dchi2cl[MAX_XB];  
	Float_t         dalpha[MAX_XB];  
	Float_t         ddls[MAX_XB];   
	Float_t         ddca[MAX_XB];   
	Float_t         dtrk1Pt[MAX_XB]; 
	Float_t         dtrk2Pt[MAX_XB]; 
	Float_t         dtrk3Pt[MAX_XB]; 
	Float_t         dtrk1PtErr[MAX_XB]; 
	Float_t         dtrk2PtErr[MAX_XB]; 
	Float_t         dtrk3PtErr[MAX_XB]; 
	Float_t         dtktkResmass[MAX_XB]; 
	Float_t         dtktkRes_chi2cl[MAX_XB]; 
	Float_t         dtktkRes_svpvDistanceToSV[MAX_XB];   
	Float_t         dtktkRes_svpvDisErrToSV[MAX_XB];  

	Float_t         dpthat;
	Int_t           dgen[MAX_XB];   //[Dsize]
	Int_t           dsGen[MAX_XB];   //[Dsize]
	Float_t         dgenpt[MAX_XB];   //[Dsize]
	Int_t           dgencollisionId[MAX_XB];   //[Dsize]
	Float_t         dgenBAncestorpt[MAX_XB];   //[Dsize]
	Int_t           dgenBAncestorpdgId[MAX_XB];   //[Dsize]
	Int_t           dgenfromgenPV[MAX_XB];   //[Dsize]

	// branches
	nt->Branch("size",&size);
	nt->Branch("hiBin",&dhiBin);
	nt->Branch("Dmass",                      dmass                       ,"dmass[size]/F");                      
	nt->Branch("Dpt",                        dpt                         ,"dpt[size]/F");   
	nt->Branch("Dchi2cl",                    dchi2cl                     ,"dchi2cl[size]/F");  
	nt->Branch("Dalpha",                     dalpha                      ,"dalpha[size]/F");  
	nt->Branch("Ddls",                       ddls                        ,"ddls[size]/F");   
	nt->Branch("Ddca",                       ddca                        ,"ddca[size]/F");   
	nt->Branch("Dtrk1Pt",                    dtrk1Pt                     ,"dtrk1Pt[size]/F"); 
	nt->Branch("Dtrk2Pt",                    dtrk2Pt                     ,"dtrk2Pt[size]/F"); 
	nt->Branch("Dtrk3Pt",                    dtrk3Pt                     ,"dtrk3Pt[size]/F"); 
	nt->Branch("Dtrk1PtErr",                 dtrk1PtErr                  ,"dtrk1PtErr[size]/F"); 
	nt->Branch("Dtrk2PtErr",                 dtrk2PtErr                  ,"dtrk2PtErr[size]/F"); 
	nt->Branch("Dtrk3PtErr",                 dtrk3PtErr                  ,"dtrk3PtErr[size]/F"); 
	nt->Branch("DtktkResmass",               dtktkResmass                ,"dtktkResmass[size]/F"); 
	nt->Branch("DtktkRes_chi2cl",            dtktkRes_chi2cl             ,"dtktkRes_chi2cl[size]/F"); 
	nt->Branch("DtktkRes_svpvDistanceToSV",  dtktkRes_svpvDistanceToSV   ,"dtktkRes_svpvDistanceToSV[size]/F");
	nt->Branch("DtktkRes_svpvDisErrToSV",    dtktkRes_svpvDisErrToSV     ,"dtktkRes_svpvDisErrToSV[size]/F");  

	if(!REAL){
	// Dgen info
	nt->Branch("pthat",&dpthat);
	nt->Branch("DsGen", dsGen,"dsGen[size]/I");
	nt->Branch("Dgenpt", dgenpt,"dgenpt[size]/F");
	nt->Branch("DgencollisionId", dgencollisionId,"dgencollisionId[size]/I");
	nt->Branch("DgenBAncestorpt", dgenBAncestorpt,"dgenBAncestorpt[size]/F");
	nt->Branch("DgenBAncestorpdgId", dgenBAncestorpdgId,"dgenBAncestorpdgId[size]/I");
	nt->Branch("DgenfromgenPV",dgenfromgenPV, "dgenfromgenPV[size]/I");
	}
	*/

	cout<<"--- Building trees finished"<<endl;
	cout<<"--- Check the number of events for three trees"<<endl;
	cout<<t_Ds->GetEntries()<<" "<<t_Hlt->GetEntries()<<" "<<t_Hi->GetEntries();
	cout<<" "<<t_Skim->GetEntries()<<endl;
	cout<<endl;

	cout<<"--- Processing events"<<endl;

	// skim cuts // minimum cuts, try to keep same with Dfiner
	float trkPtCut=0.7;
	float trkEtaCut=1.5;
	float trkPtSnCut=0.3;

	float DsvpvDisSnCut=1.5; // min : pp 2.5 for Dpt>5, 4 for Dpt<5 , (need to rerun? 
	float DyCut=1;
	float Dchi2clCut=0.02; // min 0.05 in DsFinder
	float DalphaCut=0.25;

	// float DtktkRes_chi2clCut=0.02;
	// float DtktkResmassCutMean=1.01946;
	// float DtktkResmassCutWidth=0.009;
	float DptCut=2;
	if(isPbPb){DptCut=6;}

	// possiblely Dtktk_svpvtoSvSn
	// Dca cut or DcaSn (for prompt)

	int Dcount=0;
	int Dcount_trigger=0;

//	endEntries=50000;

	for(int i=startEntries;i<endEntries;i++)
	{
		t_Ds->GetEntry(i);
		t_Hlt->GetEntry(i);
		t_Skim->GetEntry(i);
		t_Hi->GetEntry(i);

		// cout<<"Trksize = "<<Trksize<<endl;
		// cout<<"PVxE = "<<PVxE<<endl;
		// cout<<"PVyE = "<<PVyE<<endl;
		// cout<<"PVzE = "<<PVzE<<endl;

		dhiBin=hiBin;
		dEvtvx=vx;
		dEvtvy=vy;
		dEvtvz=vz;

		dNpart=0;
		dNcoll=0;

		if(!REAL)	{
			dpthat=pthat;
			dweight=weight;
		}
		if(isPbPb){
			dNpart=Npart;
			dNcoll=Ncoll;
		}


		

	// fill Gentree before apply any cut
		if(!REAL){
			t_Gen->GetEntry(i);
			ntGen->Fill();
		}

		if(i%5000==0) cout<<setw(7)<<i<<" / "<<endEntries<<endl;

		// if some file has no event... hadd might have problem , need to use another marco to merge 
		// basic event selection
		// if(skim){
			// event selection
			if(abs(PVz)>15) continue;
			if(!isPbPb ){
				if(!pBeamScrapingFilter) continue;
				if(!pPAprimaryVertexFilter) continue;
			}
			if(isPbPb){
				if(!pclusterCompatibilityFilter) continue;
				if(!pprimaryVertexFilter) continue;
				if(!phfCoincFilter3) continue;  // 12 also auto pass if 3pass
			}
			// trigger selection , only apply to Data
/*
			cout<<"HLT_L1MinimumBiasHF1OR_part0_v1 = "<< HLT_L1MinimumBiasHF1OR_part0_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part1_v1 = "<< HLT_L1MinimumBiasHF1OR_part1_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part2_v1 = "<< HLT_L1MinimumBiasHF1OR_part2_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part3_v1 = "<< HLT_L1MinimumBiasHF1OR_part3_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part4_v1 = "<< HLT_L1MinimumBiasHF1OR_part4_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part5_v1 = "<< HLT_L1MinimumBiasHF1OR_part5_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part6_v1 = "<< HLT_L1MinimumBiasHF1OR_part6_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part7_v1 = "<< HLT_L1MinimumBiasHF1OR_part7_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part8_v1 = "<< HLT_L1MinimumBiasHF1OR_part8_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part9_v1 = "<< HLT_L1MinimumBiasHF1OR_part9_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part10_v1 = "<< HLT_L1MinimumBiasHF1OR_part10_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part11_v1 = "<< HLT_L1MinimumBiasHF1OR_part11_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part12_v1 = "<< HLT_L1MinimumBiasHF1OR_part12_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part13_v1 = "<< HLT_L1MinimumBiasHF1OR_part13_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part14_v1 = "<< HLT_L1MinimumBiasHF1OR_part14_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part15_v1 = "<< HLT_L1MinimumBiasHF1OR_part15_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part16_v1 = "<< HLT_L1MinimumBiasHF1OR_part16_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part17_v1 = "<< HLT_L1MinimumBiasHF1OR_part17_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part18_v1 = "<< HLT_L1MinimumBiasHF1OR_part18_v1<<endl;
			cout<<"HLT_L1MinimumBiasHF1OR_part19_v1 = "<< HLT_L1MinimumBiasHF1OR_part19_v1<<endl;
*/

			if(REAL){
				if(!isPbPb){
					if(!( HLT_L1MinimumBiasHF1OR_part1_v1  || HLT_L1MinimumBiasHF1OR_part2_v1  || HLT_L1MinimumBiasHF1OR_part3_v1   
								|| HLT_L1MinimumBiasHF1OR_part4_v1  || HLT_L1MinimumBiasHF1OR_part5_v1  || HLT_L1MinimumBiasHF1OR_part6_v1
								|| HLT_L1MinimumBiasHF1OR_part7_v1  || HLT_L1MinimumBiasHF1OR_part8_v1  || HLT_L1MinimumBiasHF1OR_part9_v1
								|| HLT_L1MinimumBiasHF1OR_part10_v1 || HLT_L1MinimumBiasHF1OR_part11_v1 || HLT_L1MinimumBiasHF1OR_part12_v1
								|| HLT_L1MinimumBiasHF1OR_part13_v1 || HLT_L1MinimumBiasHF1OR_part14_v1 || HLT_L1MinimumBiasHF1OR_part15_v1
								|| HLT_L1MinimumBiasHF1OR_part16_v1 || HLT_L1MinimumBiasHF1OR_part17_v1 || HLT_L1MinimumBiasHF1OR_part18_v1
								|| HLT_L1MinimumBiasHF1OR_part19_v1) ) continue;

				}

				if(isPbPb){
					if(isPbPb==1){
						if(! (HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part1_v1 || HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part2_v1
									|| HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part3_v1)) continue;
					}else{
						if(! (HLT_HIL1MinimumBiasHF2AND_part1_v1 || HLT_HIL1MinimumBiasHF2AND_part2_v1 || HLT_HIL1MinimumBiasHF2AND_part3_v1
									|| HLT_HIL1MinimumBiasHF2AND_part4_v1 || HLT_HIL1MinimumBiasHF2AND_part5_v1 || HLT_HIL1MinimumBiasHF2AND_part6_v1
									|| HLT_HIL1MinimumBiasHF2AND_part7_v1 || HLT_HIL1MinimumBiasHF2AND_part8_v1 || HLT_HIL1MinimumBiasHF2AND_part9_v1
									|| HLT_HIL1MinimumBiasHF2AND_part10_v1 || HLT_HIL1MinimumBiasHF2AND_part11_v1 ) ) continue;
					}
				}
			} // end REAL trigger

		// } // end skim

		Dcount_trigger++;
		// cout<<"check 2"<<endl;


		dhiBin=hiBin;
		if(!REAL)	{dpthat=pthat;}


		for(int j=0; j<Dsize ; j++){
			// if(skim){
				if(!Dtrk1highPurity[j] || !Dtrk2highPurity[j] || !Dtrk3highPurity[j]) continue;
				if(Dtrk1Pt[j]<trkPtCut || Dtrk2Pt[j]<trkPtCut || Dtrk3Pt[j]<trkPtCut) continue;
				if(abs(Dtrk1Eta[j])>trkEtaCut || abs(Dtrk1Eta[j])>trkEtaCut || abs(Dtrk1Eta[j])>trkEtaCut) continue;
				if(Dtrk1PtErr[j]/Dtrk1Pt[j] > trkPtSnCut || Dtrk2PtErr[j]/Dtrk2Pt[j] > trkPtSnCut || Dtrk3PtErr[j]/Dtrk3Pt[j] > trkPtSnCut) continue;
				if((DsvpvDistance[j]/DsvpvDisErr[j])< DsvpvDisSnCut) continue; 
				 if(abs(Dy[j])>DyCut) continue;
				if(Dchi2cl[j]<Dchi2clCut) continue;
				if(Dalpha[j]>DalphaCut) continue;
				if(Dpt[j]<DptCut) continue;

				// if(!REAL){
					// if(DgencollisionId!=0) continue;
				// }		

				//# add phi cut later
				// if(DtktkRes_chi2cl[i]<DtktkRes_chi2clCut) continue;
				// if(DtktkResmass[i]<DtktkResmassCutMean-DtktkResmassCutWidth || DtktkResmass[i]>DtktkResmassCutMean+DtktkResmassCutWidth)continue;
			// } // end skim cut

			dmass=Dmass[j];
			dpt=Dpt[j];
			dchi2cl=Dchi2cl[j];  
			dalpha=Dalpha[j];
			ddls=DsvpvDistance[j]/DsvpvDisErr[j];
	
			double PVE=sqrt(PVxE*PVxE+PVyE*PVyE+PVzE*PVzE);		
			ddlsPPVE=(DsvpvDistance[j]+PVE) /DsvpvDisErr[j];
			ddlsMPVE=(DsvpvDistance[j]-PVE) /DsvpvDisErr[j];



			ddca=Ddca[j]; //
			dtrk1Pt=Dtrk1Pt[j];
			dtrk2Pt=Dtrk2Pt[j];
			dtrk3Pt=Dtrk3Pt[j];
			// dtrk1PtErr                  =Dtrk1PtErr[j];
			// dtrk2PtErr                  =Dtrk2PtErr[j];
			// dtrk3PtErr                  =Dtrk3PtErr[j]; 
			dtktkResmass                =DtktkResmass[j];
			dtktkRes_chi2cl             =DtktkRes_chi2cl[j]; 
			dtktkRes_svpvDistanceToSV   =DtktkRes_svpvDistanceToSV[j];
			dtktkRes_dlsToSV            =dtktkRes_svpvDistanceToSV/DtktkRes_svpvDisErrToSV[j]; 

			dsGen=-1;
			dgenBAncestorpt=-1;
			dgencollisionId=-1;


			if(!REAL){
				// dgen              =Dgen[j];   
				dsGen             =DsGen[j];  
				dgenpt            =Dgenpt[j]; 
				dgencollisionId   =DgencollisionId[j];  
				dgenBAncestorpt   =DgenBAncestorpt[j];  
				dgenBAncestorpdgId=DgenBAncestorpdgId[j]; 
				dgenfromgenPV     =DgenfromgenPV[j];                
			}

			double GausW=PVE;
			double GausM=Rdm.Gaus(0,GausW);

			if(dpt>8 && dpt<10 && dchi2cl>0.03 && dalpha<0.12 && dsGen==23333&& dgencollisionId==0 && dgenBAncestorpt<=0){
			h_Ddls->Fill(DsvpvDistance[j]/DsvpvDisErr[j],weight);
			// cout<<"GausM = "<<GausM<<endl;
			// double GausBinLow=h_Ddls_Gaus->FindBin()	
			h_Ddls_Gaus->Fill((DsvpvDistance[j]+GausM)/DsvpvDisErr[j],weight);
			}
			// cout<<"PVE = "<<PVE<<" , DsvpvDistance[j] = "<< DsvpvDistance[j] <<"DsvpvDisErr[j"<< DsvpvDisErr[j]<<" , GausM =  "<<GausM <<endl;


			nt_Ds->Fill();
			Dcount++;
	
		} // end for j<Dsize

	} // for entries

	cout<<"Dcount = "<<Dcount<<endl;
  cout<<"Dcount_trigger= "<<Dcount_trigger<<endl;

	if(Dcount_trigger>0 ){
  cout<<"Dcount_trigger found "<<Dcount_trigger<<endl;
 }

	if(Dcount>0 ){
  cout<<"Dcount found "<<Dcount<<endl;
 }

	fout->cd();
	nt_Ds->Write("",TObject::kOverwrite); // save only the new version of this tree
	if(!REAL){
	ntGen->Write("",TObject::kOverwrite);}
	//fout->Write();


	h_Ddls->Write("",TObject::kOverwrite);
	h_Ddls_Gaus->Write("",TObject::kOverwrite);

	fout->Close();

	return 1;

} // end Make

int main(int argc, char *argv[])
{
	if(argc==3)
	{
		MakeDsMinTree_ForDdlsPVE(argv[1], argv[2]);
	}
	else if(argc==5)
	{
		cout<<"real = "<<atoi(argv[3])<<" ,isPbPb = "<<atoi(argv[4])<<endl;
		MakeDsMinTree_ForDdlsPVE(argv[1], argv[2],atoi(argv[3]), atoi(argv[4]));
	}
	else
	{
		std::cout << "Usage: loop.exe <input_collection> <output_file> <isReal> <isPbPb>" << std::endl;
		return 1;
	}

	return 0;
}

