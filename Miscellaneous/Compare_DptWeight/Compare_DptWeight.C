#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/DsMinTreeLoad.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"

#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"



int Compare_DptWeight(int isPbPb=0, int PNPrompt = 0, int DsChannel=0){

	InitStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  TString MCFile;
  if(isPbPb==0 && PNPrompt==0 && DsChannel==0){MCFile=MCFile_Prompt_phikkpi_pp_Merge;  }
  if(isPbPb==0 && PNPrompt==1 && DsChannel==0){MCFile=MCFile_NonPrompt_phikkpi_pp_Merge;  }
  if(isPbPb==0 && PNPrompt==0 && DsChannel==1){MCFile=MCFile_Prompt_f0980kkpi_pp_Merge;  }
  if(isPbPb==0 && PNPrompt==1 && DsChannel==1){MCFile=MCFile_NonPrompt_f0980kkpi_pp_Merge;  }

  if(isPbPb && PNPrompt==0 && DsChannel==0){MCFile=MCFile_Prompt_phikkpi_PbPb_Merge;  }
  if(isPbPb && PNPrompt==1 && DsChannel==0){MCFile=MCFile_NonPrompt_phikkpi_PbPb_Merge;  }
  if(isPbPb && PNPrompt==0 && DsChannel==1){MCFile=MCFile_Prompt_f0980kkpi_PbPb_Merge;  }
  if(isPbPb && PNPrompt==1 && DsChannel==1){MCFile=MCFile_NonPrompt_f0980kkpi_PbPb_Merge;  }

	TString S_PNP="Prompt";
	if(PNPrompt==1){S_PNP="NonPrompt";}

  TString inputFile=MCFile;

  TFile *fin=TFile::Open(inputFile.Data());
	TTree *ntDs=(TTree*)fin->Get("ntDs");
	TTree *ntGen=(TTree*)fin->Get("ntGen");


	int nbin_Dchi2cl=40;
	double binLow_Dchi2cl=0.0;
	double binHigh_Dchi2cl=1.0;

	double DptLow=20;
	double DptHigh=40;

	TH1D *h_Dchi2cl_Weighted= new TH1D("h_Dchi2cl_Weighted","h_Dchi2cl_Weighted",nbin_Dchi2cl,binLow_Dchi2cl,binHigh_Dchi2cl);
	h_Dchi2cl_Weighted->Sumw2();
	ntDs->Draw("Dchi2cl>>h_Dchi2cl_Weighted",Form("(Dpt>%f&&Dpt<%f&&DsGen==23333)*weight*DgenptSampleWeight*RecoFONLLWeight*RecoD0DataWeight",DptLow,DptHigh));
	h_Dchi2cl_Weighted->Scale(1/h_Dchi2cl_Weighted->Integral());
	
	TH1D *h_Dchi2cl_NoWeight= new TH1D("h_Dchi2cl_NoWeight","h_Dchi2cl_NoWeight",nbin_Dchi2cl,binLow_Dchi2cl,binHigh_Dchi2cl);
	h_Dchi2cl_NoWeight->Sumw2();
	ntDs->Draw("Dchi2cl>>h_Dchi2cl_NoWeight",Form("(Dpt>%f&Dpt<%f&&DsGen==23333)",DptLow,DptHigh));	
	h_Dchi2cl_NoWeight->Scale(1/h_Dchi2cl_NoWeight->Integral());


	TCanvas *c2=DrawCompare(h_Dchi2cl_Weighted,h_Dchi2cl_NoWeight,"Vertex Probability");
	c2->cd();
	TLatex *tl = new TLatex();
	tl->DrawLatexNDC(0.25,0.5,Form("Dpt %.0f-%.0f GeV",DptLow,DptHigh));
	tl->DrawLatexNDC(0.25,0.43,Form("pp %s MC",S_PNP.Data()));
	c2->SaveAs(Form("Dpt_%.0f%.0f_%s_pp_Dchi2cl.pdf",DptLow,DptHigh,S_PNP.Data()));


	return 0;

}
