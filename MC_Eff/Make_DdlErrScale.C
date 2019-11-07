#include "../include/uti.h"
#include "../include/DsMinTreeLoad.h"
#include "../include/parameters.h"


#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"


#include <TString.h>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>
#include <TLatex.h>


int Make_DdlErrScale(int isPbPb=0){

	TString s_ppPbPb="pp";
  int nbin_pt=nbin_pt_pp;
	double *bins_pt=bins_pt_pp;
  int startbin=0;
	if(isPbPb){
	s_ppPbPb="PbPb";
  nbin_pt=nbin_pt_PbPb3;
	bins_pt=bins_pt_PbPb3;
  startbin=2;
	}

	TFile *f_Scale[nbin_pt];
	for(int i=startbin; i<nbin_pt; i++){
		f_Scale[i]=TFile::Open();
	}


	TFile *f_DsTree=TFile::Open("/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/DsMinTree_MC_GenSampleMerge_pp_Prompt_phi.root");


}
