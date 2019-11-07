#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/DsMinTreeLoad.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"

#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"

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

using namespace std;


int plot_weight(int isPbPb=0, int PNPrompt =0, int DsChannel=0, TString MCFileList="./MC_List/DsMinTree_PbPb_MC_Prompt_phikkpi.lis", TString inMergeFile="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/PbPb_MC/DsMinTree_PbPb_MC_Prompt_phikkpi.root"){

	TLatex *tla=new TLatex();

	InitStyle();
	setTDRStyle();

  TString str_isPbPb="pp";
  TString str_PNPrompt="Prompt";
  TString str_DsChannel="phi";

	TString str_DsChannel_plot="#phi#pi";

  Int_t DsGenTrue=23333;

  if(isPbPb>0)  {
    str_isPbPb=Form("PbPb");
  }
  if(PNPrompt==1){
    str_PNPrompt="NonPrompt";
  }
  if(DsChannel==1){
    str_DsChannel="f0980";
	  str_DsChannel_plot="f0#pi";
    DsGenTrue=24433;
  }
  if(DsChannel==2){
    str_DsChannel="kstar892";
	  str_DsChannel_plot="K^{*}K";
    DsGenTrue=25544;
  }

  // inMergeFile=Form("./root_output/DsMinTree_MC_GenSampleMerge_%s_%s_%s.root",str_isPbPb.Data(), str_PNPrompt.Data(),str_DsChannel.Data());
  inMergeFile=Form("/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree//DsMinTree_MC_GenSampleMerge_%s_%s_%s.root",str_isPbPb.Data(), str_PNPrompt.Data(),str_DsChannel.Data());

	cout<<"inMergeFile = "<<inMergeFile<<endl;
cout<<__LINE__<<endl;

  if(!TFile::Open(inMergeFile.Data()))
  {
    cout<<" fail to open merge file input , exit"<<endl;
    return 1;
  }

cout<<__LINE__<<endl;
  TFile *f_Merge=TFile::Open(inMergeFile.Data());

	TTree *ntGen=(TTree*)f_Merge->Get("ntGen");
cout<<__LINE__<<endl;

	if(isPbPb){
cout<<__LINE__<<endl;

		gStyle->SetOptStat(0);
		TString GenWeight="weight*GptSampleWeight*PbPbVzWeight*Ncoll*GenFONLLWeight*GenD0DataWeight";
		TH1D *h_hibin= new TH1D("h_hibin","h_hibin",40,0,200); h_hibin->Sumw2();
cout<<__LINE__<<endl;

		TCanvas *c_ncoll=new TCanvas("c_ncoll","c_ncoll",800,600);
		SetCanvas(c_ncoll);
		c_ncoll->cd();
		if(PNPrompt==0){
		ntGen->Draw("hiBin>>h_hibin","(GSignalType==1 && GcollisionId==0 && GBAncestorpt<=0) *(weight*GptSampleWeight*PbPbVzWeight*Ncoll*GenFONLLWeight*GenD0DataWeight)");}
		else if(PNPrompt==1){
		ntGen->Draw("hiBin>>h_hibin","(GSignalType==1&&GcollisionId==0&&GBAncestorpt>0)*(weight*GptSampleWeight*PbPbVzWeight*Ncoll*GenFONLLWeight*GenD0DataWeight)");}
		h_hibin->SetTitle("");
		h_hibin->GetXaxis()->SetTitle("hiBin");

cout<<__LINE__<<endl;

		SavePlotDirs(c_ncoll,Form("WeightResultVshiBin_%s_%s_%s",str_isPbPb.Data(),str_PNPrompt.Data(),str_DsChannel.Data()),{"MC","weight"});

		// return 1;
	}

cout<<__LINE__<<endl;
	
	TH1D *h_GenSamplePtWeight=(TH1D*)f_Merge->Get("h_GenSamplePtWeight");


	TCanvas *c_GenWeight=new TCanvas("c_GenWeight","c_GenWeight",c_wtopx,c_wtopy,c_W,c_H);
	c_GenWeight->cd();
	SetCanvas(c_GenWeight);
/*
  c_roofit_Data_pull[count_c_rootfit]->SetFillColor(0);
  c_roofit_Data_pull[count_c_rootfit]->SetBorderMode(0);
  c_roofit_Data_pull[count_c_rootfit]->SetFrameFillStyle(0);
  c_roofit_Data_pull[count_c_rootfit]->SetFrameBorderMode(0);
  c_roofit_Data_pull[count_c_rootfit]->SetLeftMargin( c_Lmg );
  c_roofit_Data_pull[count_c_rootfit]->SetRightMargin( c_Rmg );
  c_roofit_Data_pull[count_c_rootfit]->SetTopMargin( c_Tmg );
  c_roofit_Data_pull[count_c_rootfit]->SetBottomMargin( c_Bmg );
  c_roofit_Data_pull[count_c_rootfit]->SetTickx(0);
  c_roofit_Data_pull[count_c_rootfit]->SetTicky(0);
*/
	gPad->SetLogy();
	gPad->SetLogx();
	gStyle->SetTextFont(42);
	h_GenSamplePtWeight->GetXaxis()->SetRangeUser(2,50);
	h_GenSamplePtWeight->GetXaxis()->SetTitle("D_{S}^{#pm} gen p_{T} (GeV/c)");
	h_GenSamplePtWeight->GetXaxis()->SetTitleSize(0.05);
	h_GenSamplePtWeight->GetXaxis()->SetTitleOffset(1.0);
	h_GenSamplePtWeight->GetYaxis()->SetTitle("MC sample enhancement weight");

	h_GenSamplePtWeight->Draw();
	tla->DrawLatexNDC(textposx+0.30,textposy+0.05,Form("%s : %s D_{S}^{#pm}(%s)", str_isPbPb.Data(), str_PNPrompt.Data(), str_DsChannel_plot.Data()));

	SavePlotDirs(c_GenWeight,Form("GenSampleWeight_%s_%s_%s",str_isPbPb.Data(),str_PNPrompt.Data(),str_DsChannel.Data()),{"MC","weight"});

cout<<__LINE__<<endl;

	TH1D *h_GptwithWeightFromTree=(TH1D*)f_Merge->Get("h_GptwithWeightFromTree");
	TH1D *h_GptwithFONLLWeightFromTree=(TH1D*)f_Merge->Get("h_GptwithFONLLWeightFromTree");
	if(isPbPb){
		h_GptwithFONLLWeightFromTree=(TH1D*)f_Merge->Get("h_GptwithFONLLRaaWeightFromTree");
	}
	TH1D *h_GptwithD0DataWeightFromTree=(TH1D*)f_Merge->Get("h_GptwithD0DataWeightFromTree");

	TH1D *h_GptwithWeightFromTree_norm=(TH1D*)h_GptwithWeightFromTree->Clone();
	TH1D *h_GptwithFONLLWeightFromTree_norm=(TH1D*)h_GptwithFONLLWeightFromTree->Clone();
	TH1D *h_GptwithD0DataWeightFromTree_norm=(TH1D*)h_GptwithD0DataWeightFromTree->Clone();


/*
	TH1D *h_GptwithWeightFromTree_norm=new TH1D("h_GptwithWeightFromTree_norm","h_GptwithWeightFromTree_norm",48,2,50); h_GptwithWeightFromTree_norm->Sumw2();
	TH1D *h_GptwithFONLLWeightFromTree_norm=new TH1D("h_GptwithFONLLWeightFromTree_norm","h_GptwithFONLLWeightFromTree_norm",48,2,50); h_GptwithFONLLWeightFromTree_norm->Sumw2();
	TH1D *h_GptwithD0DataWeightFromTree_norm=new TH1D("h_GptwithD0DataWeightFromTree_norm","h_GptwithD0DataWeightFromTree_norm",48,2,50); h_GptwithD0DataWeightFromTree_norm->Sumw2();
*/
  int GSignalTypeTrue=1;
  if(DsChannel==1){GSignalTypeTrue=2;}; // f0980
  if(DsChannel==2){GSignalTypeTrue=3;}; // kstar892
  TCut cutGenTrue=Form("GSignalType==%i && GcollisionId==0 && TMath::Abs(Gy)<1",GSignalTypeTrue);
  TCut cutGenPNPrompt="GBAncestorpt<=0";
  if(PNPrompt==1){cutGenPNPrompt="GBAncestorpt>0";}
	TString PythiaWeight="(weight*GptSampleWeight*PbPbVzWeight*Ncoll)";
	TString FONLLWeight="(weight*GptSampleWeight*PbPbVzWeight*Ncoll*GenFONLLRaaWeight)";
	TString D0DataWeight="(weight*GptSampleWeight*PbPbVzWeight*Ncoll*GenFONLLWeight*GenD0DataWeight)";
	if(isPbPb==0){
	 PythiaWeight="(weight*GptSampleWeight)";
	 FONLLWeight="(weight*GptSampleWeight*GenFONLLWeight)";
	 D0DataWeight="(weight*GptSampleWeight*GenFONLLWeight*GenD0DataWeight)";
	}

	cout<<"cut  = "<< (TCut)(  cutGenTrue && cutGenPNPrompt ) * PythiaWeight<<endl;
	cout<<"cut D0 = "<< (TCut)(  cutGenTrue && cutGenPNPrompt)  * D0DataWeight<<endl;
/*
    ntGen->Draw("Gpt>>h_GptwithWeightFromTree_norm", (TCut)(  cutGenTrue && cutGenPNPrompt)  * PythiaWeight) ;
    ntGen->Draw("Gpt>>h_GptwithFONLLWeightFromTree_norm", (TCut)(  cutGenTrue && cutGenPNPrompt)  * FONLLWeight) ;
    ntGen->Draw("Gpt>>h_GptwithD0DataWeightFromTree_norm", (TCut)(  cutGenTrue && cutGenPNPrompt)  * D0DataWeight) ;
*/
cout<<__LINE__<<endl;

/*
	for(int i=1; i<=h_GptwithWeightFromTree->GetNbinsX(); i++){
		if(h_GptwithWeightFromTree->GetBinCenter(i)>=1.9 && h_GptwithWeightFromTree->GetBinCenter(i)<=43){

			h_GptwithWeightFromTree_norm->SetBinContent(i,h_GptwithWeightFromTree->GetBinContent(i));
			h_GptwithFONLLWeightFromTree_norm->SetBinContent(i,h_GptwithFONLLWeightFromTree->GetBinContent(i));
			h_GptwithD0DataWeightFromTree_norm->SetBinContent(i,h_GptwithD0DataWeightFromTree->GetBinContent(i));

			h_GptwithWeightFromTree_norm->SetBinError(i,h_GptwithWeightFromTree->GetBinError(i));
			h_GptwithFONLLWeightFromTree_norm->SetBinError(i,h_GptwithFONLLWeightFromTree->GetBinError(i));
			h_GptwithD0DataWeightFromTree_norm->SetBinError(i,h_GptwithD0DataWeightFromTree->GetBinError(i));
		}
	} // end for int

*/
			h_GptwithWeightFromTree_norm->Scale(1/h_GptwithWeightFromTree_norm->Integral());
			h_GptwithFONLLWeightFromTree_norm->Scale(1/h_GptwithFONLLWeightFromTree_norm->Integral());
			h_GptwithD0DataWeightFromTree_norm->Scale(1/h_GptwithD0DataWeightFromTree_norm->Integral());

	TCanvas *c_dis=new TCanvas("c_dis","c_dis",c_wtopx,c_wtopy,c_W,c_H);
	c_dis->cd();
	SetCanvas(c_dis);	
	gPad->SetLogy();
	// gPad->SetLogx();
	h_GptwithWeightFromTree_norm->GetXaxis()->SetTitle("D_{s}^{#pm} p_{T} (GeV/c)");
	h_GptwithWeightFromTree_norm->GetXaxis()->SetTitleSize(0.045);
	h_GptwithWeightFromTree_norm->GetXaxis()->SetTitleOffset(1.2);
	h_GptwithWeightFromTree_norm->SetMaximum(3);
	h_GptwithWeightFromTree_norm->Draw("same");
  h_GptwithFONLLWeightFromTree_norm->SetMarkerColor(2);	
	h_GptwithFONLLWeightFromTree_norm->SetLineColor(2);
	h_GptwithFONLLWeightFromTree_norm->Draw("same");
	h_GptwithD0DataWeightFromTree_norm->SetMarkerColor(4);
	h_GptwithD0DataWeightFromTree_norm->SetLineColor(4);
	h_GptwithD0DataWeightFromTree_norm->Draw("same");

	tla->DrawLatexNDC(textposx+0.30,textposy+0.05,Form("%s : %s D_{S}^{#pm}(%s)", str_isPbPb.Data(), str_PNPrompt.Data(), str_DsChannel_plot.Data()));


	TString str_FONLL="FONLL";
	if(isPbPb>0){
		if(PNPrompt==0)	{ str_FONLL="FONLL + PHSD";  }	
		else if(PNPrompt==1)	{ str_FONLL="FONLL + TAMU";  }	
	}

	TLegend *le=new TLegend(0.6,0.5,0.85,0.75);
	le->SetBorderSize(0);
	le->AddEntry(h_GptwithWeightFromTree_norm,"PYTHIA","lp");
	le->AddEntry(h_GptwithFONLLWeightFromTree_norm,str_FONLL.Data(),"lp");
	le->AddEntry(h_GptwithD0DataWeightFromTree_norm,"Data","lp");

	le->Draw("same");

	SavePlotDirs(c_dis,Form("WeightResult_%s_%s_%s",str_isPbPb.Data(),str_PNPrompt.Data(),str_DsChannel.Data()),{"MC","weight"});

	// Distribution ratio spectra

	TCanvas *c_GptRatio_FONLLtoData=new TCanvas("c_GptRatio_FONLLtoData","c_GptRatio_FONLLtoData");
	c_GptRatio_FONLLtoData->cd();

	TH1D *h_GptRatio_FONLLtoData=(TH1D*)h_GptwithFONLLWeightFromTree_norm->Clone("h_GptRatio_FONLLtoData");
	h_GptRatio_FONLLtoData->Divide(h_GptwithD0DataWeightFromTree_norm);

	h_GptRatio_FONLLtoData->SetTitle("");
	h_GptRatio_FONLLtoData->GetXaxis()->SetTitle("D_{S} p_{T} (GeV/c)");
	h_GptRatio_FONLLtoData->GetYaxis()->SetTitle("ratio");
	h_GptRatio_FONLLtoData->GetXaxis()->SetRangeUser(2,20);
	h_GptRatio_FONLLtoData->Draw();

	// TLatex *tla=new TLatex();
	tla->DrawLatexNDC(textposx,textposy-0.3,"MC p_{T} Shape Alternative/Default");
	
	
	c_GptRatio_FONLLtoData->SaveAs(Form("plots/ShapeRatio_%s_%s.png",str_isPbPb.Data(),str_PNPrompt.Data()));

	

	

	



	return 0;

}
