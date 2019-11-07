#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"

int plotEffFineBin(int PNP=0){

	gStyle->SetOptStat(0);
	gSystem->Exec("mkdir -p plots_EffFine");

	TString MCFile=MCFile_Prompt_phikkpi_pp_Merge;
	TString s_cutBpt="DgenBAncestorpt<=0";
  TCut cutGenPNPrompt="GBAncestorpt<=0";
	TString s_PNP="Prompt";
	if(PNP==1){
		MCFile=MCFile_NonPrompt_phikkpi_pp_Merge;
		s_cutBpt="DgenBAncestorpt>0";
		cutGenPNPrompt="GBAncestorpt>0";
		s_PNP="NonPrompt";
	}



  TFile *f_Prompt=TFile::Open(MCFile.Data());
	TTree *nt_Prompt=(TTree*)f_Prompt->Get("ntDs");
	TTree *ntGen_Prompt=(TTree*)f_Prompt->Get("ntGen");

	int nbin=72;
	double binLow=2;
	double binHigh=20;

	TH1D *h_Prompt_Entry=new TH1D("h_Prompt_Entry",";p_{T}",nbin,binLow,binHigh);

	// nt_Prompt->Project("h_Prompt_Entry","Dpt","weight*DgenptSampleWeight*RecoD0DataWeight*(Ddls>2.5 && DsGen==23333 && DgencollisionId==0 &&Dalpha<0.12 && Dchi2cl>0.04 && DgenBAncestorpt<=0 )");
	nt_Prompt->Project("h_Prompt_Entry","Dpt",Form("weight*DgenptSampleWeight*RecoD0DataWeight*RecoFONLLWeight*(Ddls>2.5 && DsGen==23333 && DgencollisionId==0 &&Dalpha<0.12 && Dchi2cl>0.04 && %s)",s_cutBpt.Data() ) );

	// h_Prompt_Entry->Draw();

	int GSignalTypeTrue=1;
  TCut cutGenTrue=Form("GSignalType==%i && GcollisionId==0 && TMath::Abs(Gy)<1",GSignalTypeTrue);
  TString GenWeight_D0Data="weight*GptSampleWeight*GenFONLLWeight*GenD0DataWeight";


	TH1D *h_Prompt_Gen=new TH1D("h_Prompt_Gen",";p_{T}",nbin,binLow,binHigh);
	// ntGen_Prompt->Project("h_Prompt_Gen","Gpt",(TCut)((cutGenTrue && cutGenPNPrompt)*GenWeight_D0Data.Data()));
	ntGen_Prompt->Project("h_Prompt_Gen","Gpt",(TCut)((cutGenTrue && cutGenPNPrompt)*GenWeight_D0Data.Data()));

	// ntGen_Prompt->Draw("Gpt");
	// h_Prompt_Gen->Draw();

	TCanvas *c_eff=new TCanvas("c_eff","c_eff");
	c_eff->cd();

	TH1D *h_Prompt_Eff=(TH1D*)h_Prompt_Entry->Clone("h_Prompt_Eff");
	h_Prompt_Eff->GetYaxis()->SetTitle("Efficiency");
	h_Prompt_Eff->Divide(h_Prompt_Gen);
	h_Prompt_Eff->Draw();

	TLatex *tla=new TLatex();
	tla->DrawLatexNDC(textposx,textposy,Form("pp %s #phi",s_PNP.Data()) );


	c_eff->SaveAs(Form("plots_EffFine/EffFineBin_%s.png",s_PNP.Data()));
	

	return 0;

}
