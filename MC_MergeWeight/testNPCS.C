#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"

void testNPCS(int isPbPb=0, int PNPrompt =1, int DsChannel=0){

  TString str_isPbPb="pp";
  TString str_PNPrompt="Prompt";
  TString str_DsChannel="phi";
  Int_t DsGenTrue=23333;

  if(isPbPb>0)  {
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


	TFile *f_pp=TFile::Open("./root_output/DsMinTree_MC_GenSampleMerge_pp_NonPrompt_phi.root","READ");
	TTree *t_gen_pp=(TTree*)f_pp->Get("ntGen");
	TH1D *h_genFONLL_pp = new TH1D("h_genFONLL_pp","h_genFONLL_pp",nbin_pt_pp,bins_pt_pp); h_genFONLL_pp->Sumw2();
	t_gen_pp->Project("h_genFONLL_pp","Gpt","( GSignalType==1 && GBAncestorpt >0 && abs(Gy)<1 && GcollisionId==0 )*(GptSampleWeight*GenFONLLWeight*weight) ");
	// h_genFONLL_pp->Draw();

	for(int i=0; i<nbin_pt_pp;i++){
		cout<<"bin i = "<<i<<" , h_genFONLL_pp = "<<h_genFONLL_pp->GetBinContent(i+1)<<endl;
	}

	TFile *f_PbPb=TFile::Open("./root_output/DsMinTree_MC_GenSampleMerge_PbPb_NonPrompt_phi.root","READ");
	TTree *t_gen_PbPb=(TTree*)f_PbPb->Get("ntGen");
	TH1D *h_genFONLL_PbPb = new TH1D("h_genFONLL_PbPb","h_genFONLL_PbPb",nbin_pt_PbPb3,bins_pt_PbPb3); h_genFONLL_PbPb->Sumw2();
	t_gen_PbPb->Project("h_genFONLL_PbPb","Gpt","( GSignalType==1 && GBAncestorpt >0 && abs(Gy)<1 && GcollisionId==0 )*(GptSampleWeight*GenFONLLWeight*weight*PbPbVzWeight*Ncoll) ");
	h_genFONLL_PbPb->Draw();

	for(int i=0; i<nbin_pt_PbPb3;i++){
		cout<<"bin i = "<<i<<" , h_genFONLL_PbPb*TAA = "<<h_genFONLL_PbPb->GetBinContent(i+1)* TAA0to100 <<endl;
	}






}
