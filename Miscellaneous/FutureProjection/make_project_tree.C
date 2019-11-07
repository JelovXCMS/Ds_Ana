#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"

#include <TString.h>
#include <TCanvas.h>
#include "TBranch.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLatex.h"
#include <TH2F.h>
#include <TCut.h>
#include <TStyle.h>

#include <math.h>

#include <string>
#include <TChain.h>
#include <TFile.h>

#include "TH1.h"
#include "TROOT.h"
#include "TPad.h"
#include "TRandom.h"
#include "THStack.h"
#include "TH2.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TAxis.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <TCut.h>

	double lumi_pp_HL=650; // (pb-1 INT);
	double lumi_PbPb_HL=0.2*1000 ; // (nb-1->mub-1 MB);

	double lumi_pp_15=27.4; // (pb-1)
	double lumi_PbPb_15=43.9; // (mub-1   MB)


int make_project_tree(int isPbPb=3){

	cout<<"isPbPb = "<<isPbPb<<endl;

	double lumi_pp_ratio=(lumi_pp_HL/lumi_pp_15);
	double lumi_PbPb_ratio=(lumi_PbPb_HL/lumi_PbPb_15);
	double lumi_ratio=lumi_pp_ratio;

	int nbin_pt=nbin_pt_pp;
	double *bins_pt=bins_pt_pp;

	cout<<"lumi_pp_ratio = "<<lumi_pp_ratio<<" , lumi_PbPb_ratio = "<<lumi_PbPb_ratio<<endl;

	TString fin_name="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree/output/FitFile_pp.root";
	TString fout_name="FitFile_pp_HLproject.root";


	if(isPbPb==3){
		fin_name="/scratch/halstead/p/peng43/Ds_phikkpi/BuildFitFile_FromDsMinTree/DsFitFile_PbPb3_Data.root";
		fout_name="FitFile_PbPb3_HLproject.root";
		lumi_ratio=lumi_PbPb_ratio;
		nbin_pt=nbin_pt_PbPb3;
		bins_pt=bins_pt_PbPb3;
	}

	TFile *fin=TFile::Open(fin_name.Data(),"READ");

	cout<<"check 1"<<endl;

	TTree *tree_arr[nbin_pt];
	TH1F *th1f_arr[nbin_pt];

	TFile *fout=TFile::Open(fout_name.Data(),"recreate");

	cout<<"check 2"<<endl;
	TTree *tree_arr_new[nbin_pt];
	TH1F *th1f_arr_new[nbin_pt];


 	Float_t Dmass;


	for(int i=0; i<nbin_pt; i++){
		
	cout<<"check 3"<<endl;
		tree_arr[i]=(TTree*)fin->Get(Form("t_DsMassData_pt%.0fto%.0f",bins_pt[i],bins_pt[i+1]));

	cout<<"check 4"<<endl;
		th1f_arr[i]=(TH1F*)fin->Get(Form("h_DsMassData_pt%.0fto%.0f",bins_pt[i],bins_pt[i+1]));

	cout<<"check 5"<<endl;
		 // th1f_arr[i]->Draw();

	cout<<"check 6"<<endl;
		// th1f_arr_new[i]=(TH1F*)th1f_arr[i]->Clone(Form("h_DsMassData_pt%.0fto%.0f_new",bins_pt[i],bins_pt[i+1]));
		th1f_arr_new[i]=(TH1F*)th1f_arr[i]->Clone();
		
		for(int jbin=1; jbin<=th1f_arr[i]->GetNbinsX();jbin++){
			// cout<<"bin : "<<jbin<<" , content = "<<th1f_arr_new[i]->GetBinContent(jbin)<<endl;
			// cout<<"bin : "<<jbin<<" , content = "<<th1f_arr[i]->GetBinContent(jbin)<<endl;
			th1f_arr_new[i]->SetBinContent(jbin, th1f_arr[i]->GetBinContent(jbin)*lumi_ratio );
			th1f_arr_new[i]->SetBinError(jbin,th1f_arr[i]->GetBinError(jbin)*sqrt(lumi_ratio) );
      // cout<<"bin : "<<jbin<<" , content = "<<th1f_arr_new[i]->GetBinContent(jbin)<<endl;
		}

		fout->cd();
		th1f_arr_new[i]->Write("",TObject::kOverwrite);

/*		tree_arr[i]->Draw("Dmass");

		tree_arr[i]->SetBranchAddress("Dmass",&Dmass);

		fout->cd();
		tree_arr_new[i]=new TTree(Form("t_DsMassData_pt%.0fto%.0f",bins_pt[i],bins_pt[i+1]), Form("t_DsMassData_pt%.0fto%.0f",bins_pt[i],bins_pt[i+1]));

		tree_arr_new[i]->Branch("Dmass",&Dmass,"Dmass/F");

		Long64_t nentries=tree_arr[i]->GetEntries();

		for(int ien=0; ien<nentries; ien++){
	
			tree_arr[i]->GetEntry(ien);
			
			tree_arr_new[i]->Fill();
		// cout<<Dmass<<endl;	

		}


		fout->cd();
		// tree_arr_new[i]->Write();
*/
//		fout->Close();


		 // return 0;

	}


		fout->Close();




	return 0;


}
