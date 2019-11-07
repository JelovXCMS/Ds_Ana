// V2 : using tree structure, ie Dsize, D[iDsize]

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


// adding is PbPb30-100 

Bool_t istest = false;
bool fillZeroCandEvt = true;
Bool_t detailMode=false;
int MakeDsMinTreeV2(TString infile="", TString outfile="", Bool_t REAL=false, Int_t isPbPb=1, Int_t startEntries=0, Int_t endEntries=-1, Bool_t skim=true, Bool_t gskim=true, Bool_t checkMatching=true, Bool_t iseos=false, Bool_t SkimHLTtree=false)
{
	// isPbPb=1 for centrality 30-100 trigger only, isPbPb=2 for GJ and trackonly 0-100

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


	// declare needed variables 
	//	const double MAX_XB=20000;
	Int_t           dsize;
	Float_t         dmass[MAX_XB];
	Float_t         dpt[MAX_XB];
	Float_t         dchi2cl[MAX_XB];  
	Float_t         dalpha[MAX_XB];  
	Float_t         ddls[MAX_XB];   
	Float_t         ddca[MAX_XB];   
	Float_t         dtrk1Pt[MAX_XB]; 
	Float_t         dtrk2Pt[MAX_XB]; 
	Float_t         dtrk3Pt[MAX_XB]; 
	// Float_t         dtrk1PtErr; 
	// Float_t         dtrk2PtErr; 
	// Float_t         dtrk3PtErr; 
	Float_t         dtktkResmass[MAX_XB]; 
	Float_t         dtktkRes_chi2cl[MAX_XB]; 
	Float_t         dtktkRes_svpvDistanceToSV[MAX_XB];   
	Float_t         dtktkRes_dlsToSV[MAX_XB];  

	Int_t           dhiBin;
	Float_t         dpthat;
	Float_t         dweight;
	Float_t         dEvtvx;
	Float_t         dEvtvy;
	Float_t         dEvtvz;
	Float_t         dNpart;
	Float_t         dNcoll;


	// Int_t           dgen;   
	Int_t           dsGen[MAX_XB];  
	Float_t         dgenpt[MAX_XB]; 
	Int_t           dgencollisionId[MAX_XB];  
	Float_t         dgenBAncestorpt[MAX_XB];  
	Int_t           dgenBAncestorpdgId[MAX_XB]; 
	Int_t           dgenfromgenPV[MAX_XB];                 



	// branches
	nt_Ds->Branch("hiBin",&dhiBin);

	nt_Ds->Branch("vx",&dEvtvx);
	nt_Ds->Branch("vy",&dEvtvy);
	nt_Ds->Branch("vz",&dEvtvz);
	
	nt_Ds->Branch("Dsize",&dsize);

	nt_Ds->Branch("Dmass",                      dmass      ,"Dmass[Dsize]/F");    
	nt_Ds->Branch("Dpt",                        dpt        ,"Dpt[Dsize]/F"              );   
	nt_Ds->Branch("Dchi2cl",                    dchi2cl      ,"Dchi2cl[Dsize]/F"            );   
	nt_Ds->Branch("Dalpha",                     dalpha       ,"Dalpha[Dsize]/F"             );   
	nt_Ds->Branch("Ddls",                       ddls         ,"Ddls[Dsize]/F"               );   
	nt_Ds->Branch("Ddca",                       ddca         ,"Ddca[Dsize]/F"               );   
	nt_Ds->Branch("Dtrk1Pt",                    dtrk1Pt      ,"Dtrk1Pt[Dsize]/F"            );   
	nt_Ds->Branch("Dtrk2Pt",                    dtrk2Pt      ,"Dtrk2Pt[Dsize]/F"            );   
	nt_Ds->Branch("Dtrk3Pt",                    dtrk3Pt      ,"Dtrk3Pt[Dsize]/F"            );   
	// nt_Ds->Branch("Dtrk1PtErr",                 &dtrk1PtErr               );   
	// nt_Ds->Branch("Dtrk2PtErr",                 &dtrk2PtErr               );   
	// nt_Ds->Branch("Dtrk3PtErr",                 &dtrk3PtErr               );   
	nt_Ds->Branch("DtktkResmass",               dtktkResmass ,"DtktkResmass[Dsize]/F"          );   
	nt_Ds->Branch("DtktkRes_chi2cl",            dtktkRes_chi2cl  ,"DtktkRes_chi2cl[Dsize]/F"         );   
	nt_Ds->Branch("DtktkRes_svpvDistanceToSV",  dtktkRes_svpvDistanceToSV ,"DtktkRes_svpvDistanceToSV[Dsize]/F");   
	nt_Ds->Branch("DtktkRes_dlsToSV",           dtktkRes_dlsToSV ,"DtktkResmass[Dsize]/F" );   

  nt_Ds->Branch("DsGen",                      dsGen , "DsGen[Dsize]/I"); // add for using TMVA 
	nt_Ds->Branch("DgencollisionId",            dgencollisionId,"DgencollisionId[Dsize]/I");
	nt_Ds->Branch("DgenBAncestorpt",            dgenBAncestorpt,"DgenBAncestorpt[Dsize]/F");

	if(!REAL){
		// Dgen info
		nt_Ds->Branch("pthat",&dpthat);
	  nt_Ds->Branch("weight",&dweight);
//		nt_Ds->Branch("DsGen",&dsGen);
		ntGen->Branch("pthat",&dpthat);
		ntGen->Branch("hiBin",&dhiBin);
		ntGen->Branch("weight",&dweight);

	  ntGen->Branch("vx",&dEvtvx);
	  ntGen->Branch("vy",&dEvtvy);
	  ntGen->Branch("vz",&dEvtvz);





		nt_Ds->Branch("Dgenpt",dgenpt,"Dgenpt[Dsize]/F");
		// nt_Ds->Branch("DgencollisionId",&dgencollisionId);
//		nt_Ds->Branch("DgenBAncestorpt",&dgenBAncestorpt);
		nt_Ds->Branch("DgenBAncestorpdgId",dgenBAncestorpdgId,"DgenBAncestorpdgId[Dsize]/I");
		nt_Ds->Branch("DgenfromgenPV",dgenfromgenPV,"DgenfromgenPV[Dsize]/I");

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

		dhiBin=hiBin;
		dEvtvx=vx;
		dEvtvy=vy;
		dEvtvz=vz;

		dNpart=0;
		dNcoll=0;
		dsize=0;

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

			dmass[dsize]=Dmass[j];
			dpt[dsize]=Dpt[j];
			dchi2cl[dsize]=Dchi2cl[j];  
			dalpha[dsize]=Dalpha[j];
			ddls[dsize]=DsvpvDistance[j]/DsvpvDisErr[j];
			ddca[dsize]=Ddca[j]; //
			dtrk1Pt[dsize]=Dtrk1Pt[j];
			dtrk2Pt[dsize]=Dtrk2Pt[j];
			dtrk3Pt[dsize]=Dtrk3Pt[j];
			// dtrk1PtErr                  =Dtrk1PtErr[j];
			// dtrk2PtErr                  =Dtrk2PtErr[j];
			// dtrk3PtErr                  =Dtrk3PtErr[j]; 
			dtktkResmass[dsize]                =DtktkResmass[j];
			dtktkRes_chi2cl[dsize]             =DtktkRes_chi2cl[j]; 
			dtktkRes_svpvDistanceToSV[dsize]   =DtktkRes_svpvDistanceToSV[j];
			dtktkRes_dlsToSV[dsize]            =dtktkRes_svpvDistanceToSV[j]/DtktkRes_svpvDisErrToSV[j]; 

			dsGen[dsize]=-1;
			dgenBAncestorpt[dsize]=-1;
			dgencollisionId[dsize]=-1;


			if(!REAL){
				// dgen              =Dgen[j];   
				dsGen[dsize]             =DsGen[j];  
				dgenpt[dsize]            =Dgenpt[j]; 
				dgencollisionId[dsize]   =DgencollisionId[j];  
				dgenBAncestorpt[dsize]   =DgenBAncestorpt[j];  
				dgenBAncestorpdgId[dsize]=DgenBAncestorpdgId[j]; 
				dgenfromgenPV[dsize]     =DgenfromgenPV[j];                
			}

			// nt_Ds->Fill();
			Dcount++;
			dsize++;	
		} // end for j<Dsize
		// if(dsize>0){
		nt_Ds->Fill();
		// }

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
	fout->Close();

	return 1;

} // end Make

int main(int argc, char *argv[])
{
	if(argc==3)
	{
		MakeDsMinTreeV2(argv[1], argv[2]);
	}
	else if(argc==5)
	{
		cout<<"real = "<<atoi(argv[3])<<" ,isPbPb = "<<atoi(argv[4])<<endl;
		MakeDsMinTreeV2(argv[1], argv[2],atoi(argv[3]), atoi(argv[4]));
	}
	else
	{
		std::cout << "Usage: loop.exe <input_collection> <output_file> <isReal> <isPbPb>" << std::endl;
		return 1;
	}

	return 0;
}

