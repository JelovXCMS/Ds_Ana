// need to count MB trigged event only

#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <iostream>
#include <fstream>
#include <string>
#include <TString.h>

using namespace std;

int CountNEvents(string FileList="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/input_list/Ds_PbPb_Data_HIMB10_TrkOnly.lis"){

	if (FileList=="")
	{
		cout<<"FileList is empty , exit "<<endl;
		return -1;
	}

	Long64_t NTotal=0;
	Long64_t NTotal_noPvz=0;
	Long64_t NTotal_Trigger=0;

	ifstream infile(FileList.c_str());
	std::string filename;
	
	int ifile=0;
	while(infile>>filename){
		ifile++;
		if(!TFile::Open(filename.c_str())) {
			cout<<"can not open file : "<< filename.c_str()<<" exit or skip"<<endl;
			return -2;
		}
		TFile *fin = TFile::Open(filename.c_str());
		TTree *t_Hlt=(TTree*)fin->Get("ntHlt");
		TTree *t_Reco=(TTree*)fin->Get("ntDPhikkpi");
		TTree *t_Hi=(TTree*)fin->Get("ntHi");
		TTree *t_Skim=(TTree*)fin->Get("ntSkim");
		Long64_t NLocal=0;
		Long64_t NLocal_noPvz=0;
		Long64_t NLocal_Trigger=0;
		// NLocal=t_Hi->GetEntries();
		// NLocal=t_Hi->GetEntriesFast();
		t_Hlt->AddFriend(t_Reco);
		t_Hlt->AddFriend(t_Hi);
		t_Hlt->AddFriend(t_Skim);
		NLocal=t_Hlt->GetEntries("abs(Pvz)<15 && pclusterCompatibilityFilter && pprimaryVertexFilter && phfCoincFilter3 &&                                             ( HLT_HIL1MinimumBiasHF2AND_part1_v1 || HLT_HIL1MinimumBiasHF2AND_part2_v1 || HLT_HIL1MinimumBiasHF2AND_part3_v1 || HLT_HIL1MinimumBiasHF2AND_part4_v1       || HLT_HIL1MinimumBiasHF2AND_part5_v1 || HLT_HIL1MinimumBiasHF2AND_part6_v1 || HLT_HIL1MinimumBiasHF2AND_part7_v1 || HLT_HIL1MinimumBiasHF2AND_part8_v1       || HLT_HIL1MinimumBiasHF2AND_part9_v1 || HLT_HIL1MinimumBiasHF2AND_part10_v1 || HLT_HIL1MinimumBiasHF2AND_part11_v1) ");
		NLocal_noPvz=t_Hlt->GetEntries("pclusterCompatibilityFilter && pprimaryVertexFilter && phfCoincFilter3 &&                                             ( HLT_HIL1MinimumBiasHF2AND_part1_v1 || HLT_HIL1MinimumBiasHF2AND_part2_v1 || HLT_HIL1MinimumBiasHF2AND_part3_v1 || HLT_HIL1MinimumBiasHF2AND_part4_v1       || HLT_HIL1MinimumBiasHF2AND_part5_v1 || HLT_HIL1MinimumBiasHF2AND_part6_v1 || HLT_HIL1MinimumBiasHF2AND_part7_v1 || HLT_HIL1MinimumBiasHF2AND_part8_v1       || HLT_HIL1MinimumBiasHF2AND_part9_v1 || HLT_HIL1MinimumBiasHF2AND_part10_v1 || HLT_HIL1MinimumBiasHF2AND_part11_v1) ");
		// NLocal_Trigger=t_Hlt->GetEntries("                ( HLT_HIL1MinimumBiasHF2AND_part1_v1 || HLT_HIL1MinimumBiasHF2AND_part2_v1 || HLT_HIL1MinimumBiasHF2AND_part3_v1 || HLT_HIL1MinimumBiasHF2AND_part4_v1       || HLT_HIL1MinimumBiasHF2AND_part5_v1 || HLT_HIL1MinimumBiasHF2AND_part6_v1 || HLT_HIL1MinimumBiasHF2AND_part7_v1 || HLT_HIL1MinimumBiasHF2AND_part8_v1       || HLT_HIL1MinimumBiasHF2AND_part9_v1 || HLT_HIL1MinimumBiasHF2AND_part10_v1 || HLT_HIL1MinimumBiasHF2AND_part11_v1) ");

		NTotal+=NLocal;
		NTotal_noPvz+=NLocal_noPvz;
		// NTotal_Trigger+=NLocal_Trigger;

		if(ifile%1000==0){
			cout<<"ifile = "<<ifile<<endl;
		}

		// cout<<"NLocal = "<<NLocal<<" , NTotal = "<<NTotal<<endl;

		// must delete , memory issue
		delete t_Hlt;
		delete t_Reco;
		delete t_Hi;
		delete t_Skim;
		delete fin;
	}

	cout<<"finish FileList = "<<FileList.c_str()<<" , Total Files : "<< ifile <<" , NTotal = "<<NTotal<<" , NTotal_Trigger = "<<NTotal_Trigger<<" , NTotal_noPvz = "<<NTotal_noPvz<<endl;


	return 0;

}
