#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TGraph.h>
#include <iostream>
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"


using namespace std;


void PythiaRef(){

	// using bins_pt_Norm && nbin_pt_Norm from parameters bins_pt_Norm[]={0,2,3,4,5,6,8,10,20,40,100} 

	InitStyle();
	double TotalCS=1.01e12; // pb
	double FilterEff_D0=0.00509;
	double FilterEff_Ds=0.00118;

	TFile *f_out=new TFile("./output/PythiaRef.root","recreate");
/*
	TFile *f_D0=TFile::Open("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/MCWeightMerge/DNtuple_PromptPP_mergeWeight.root","READ");
	TTree *ntGen_D0=(TTree*)f_D0->Get("ntGen");

	f_out->cd();
	TH1D *h_Genpt_D0=new TH1D("h_Genpt_D0","h_Genpt_D0",nbin_pt_Norm,bins_pt_Norm); h_Genpt_D0->Sumw2();
	TH1D *h_Genpt_yall_D0=new TH1D("h_Genpt_yall_D0","h_Genpt_yall_D0",nbin_pt_Norm,bins_pt_Norm); h_Genpt_yall_D0->Sumw2();
	ntGen_D0->Project("h_Genpt_D0","Gpt","((GisSignal==1 || GisSignal==2) && GcollisionId==0 && TMath::Abs(Gy)<1 && GBAncestorpt<=0 ) * GptSampleWeight ");
	ntGen_D0->Project("h_Genpt_yall_D0","Gpt","((GisSignal==1 || GisSignal==2) && GcollisionId==0 && GBAncestorpt<=0 ) * GptSampleWeight ");
*/
	TFile *f_D0=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_MergeWeight/root_output_D0/D0_MC_GenSampleMerge_pp_Prompt.root");
	TH1D *h_Genpt_D0=(TH1D*)f_D0->Get("h_GptWeighted_Gy1");
	TH1D *h_Genpt_yall_D0=(TH1D*)f_D0->Get("h_GptWeighted_GyAll");

	double integral_y1_D0=h_Genpt_D0->Integral();
	double integral_yall_D0=h_Genpt_yall_D0->Integral();
	
	cout<<"D0 y1 : "	<<integral_y1_D0<<" , integral_yall = "<<integral_yall_D0<<endl;

	h_Genpt_D0->Scale(TotalCS*FilterEff_D0/integral_yall_D0/2);
	divideBinWidth(h_Genpt_D0);
	// Ds

	// TFile *f_Ds=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_MergeWeight/root_output/DsMinTree_MC_GenSampleMerge_pp_Prompt_phi.root","READ");
	TFile *f_Ds=TFile::Open("/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/DsMinTree_MC_GenSampleMerge_pp_Prompt_phi.root","READ");
	TTree *ntGen_Ds=(TTree*)f_Ds->Get("ntGen");

	f_out->cd();
	TH1D *h_Genpt_Ds=new TH1D("h_Genpt_Ds","h_Genpt_Ds",nbin_pt_Norm,bins_pt_Norm); h_Genpt_Ds->Sumw2();
	TH1D *h_Genpt_yall_Ds=new TH1D("h_Genpt_yall_Ds","h_Genpt_yall_Ds",nbin_pt_Norm,bins_pt_Norm); h_Genpt_yall_Ds->Sumw2();
	ntGen_Ds->Project("h_Genpt_Ds","Gpt","(GSignalType==1 && GcollisionId==0 && TMath::Abs(Gy)<1 && GBAncestorpt<=0 ) * GptSampleWeight*weight ");
	ntGen_Ds->Project("h_Genpt_yall_Ds","Gpt","( GSignalType==1 && GcollisionId==0 && GBAncestorpt<=0 ) * GptSampleWeight*weight ");

	double integral_y1_Ds=h_Genpt_Ds->Integral();
	double integral_yall_Ds=h_Genpt_yall_Ds->Integral();
	
	cout<<"Ds y1 : "	<<integral_y1_Ds<<" , integral_yall = "<<integral_yall_Ds<<endl;

	h_Genpt_Ds->Scale(TotalCS*FilterEff_Ds/integral_yall_Ds/2);
	divideBinWidth(h_Genpt_Ds);

	TH1D *h_Genpt_Ds_fine=new TH1D("h_Genpt_Ds_fine","h_Genpt_Ds_fine",40,0,40); h_Genpt_Ds_fine->Sumw2();
	TH1D *h_Genpt_yall_Ds_fine=new TH1D("h_Genpt_yall_Ds_fine","h_Genpt_yall_Ds_fine",40,0,40); h_Genpt_yall_Ds_fine->Sumw2();
	ntGen_Ds->Project("h_Genpt_Ds_fine","Gpt","(GSignalType==1 && GcollisionId==0 && TMath::Abs(Gy)<1 && GBAncestorpt<=0 ) * GptSampleWeight*weight ");
	ntGen_Ds->Project("h_Genpt_yall_Ds_fine","Gpt","( GSignalType==1 && GcollisionId==0 && GBAncestorpt<=0 ) * GptSampleWeight*weight ");

	integral_y1_Ds=h_Genpt_Ds_fine->Integral();
	integral_yall_Ds=h_Genpt_yall_Ds_fine->Integral();
	
	cout<<"Ds y1 : "	<<integral_y1_Ds<<" , integral_yall = "<<integral_yall_Ds<<endl;

	h_Genpt_Ds_fine->Scale(TotalCS*FilterEff_Ds/integral_yall_Ds/2);
	divideBinWidth(h_Genpt_Ds_fine);


	TCanvas *c_Pythia=new TCanvas("c_Pythia","c_Pythia",800,800);
	c_Pythia->cd();
	gPad->SetLogy();
	h_Genpt_D0->GetXaxis()->SetRangeUser(2,40);
	h_Genpt_Ds->GetXaxis()->SetRangeUser(2,40);
	h_Genpt_D0->Draw();
	h_Genpt_Ds->SetLineColor(2);
	h_Genpt_Ds->SetMarkerColor(2);
	h_Genpt_Ds->Draw("SAME");


	TH1D *h_DsOverD0_Pythia=(TH1D*)h_Genpt_Ds->Clone("h_DsOverD0_Pythia");
	h_DsOverD0_Pythia->Divide(h_Genpt_D0);

	TCanvas *c_DsOverD0_Pythia=new TCanvas("c_DsOverD0_Pythia","c_DsOverD0_Pythia",800,800);
	c_DsOverD0_Pythia->cd();
	// h_DsOverD0_Pythia->GetXaxis()->SetRangeUser(0,40);
	h_DsOverD0_Pythia->Draw();




	f_out->cd();
	h_Genpt_D0->Write("h_Genpt_D0");
	h_Genpt_Ds->Write("h_Genpt_Ds");
	h_DsOverD0_Pythia->Write("h_DsOverD0_Pythia");
	h_Genpt_Ds_fine->Write("h_Genpt_Ds_fine");

	TGraph *gr_Genpt_Ds=new TGraph();
	gr_Genpt_Ds->SetName("gr_Genpt_Ds");
	for(int i =1; i<h_Genpt_Ds_fine->GetNbinsX()+1;i++){
		cout<<"bincenter = "<<h_Genpt_Ds_fine->GetBinCenter(i)<<endl;
		gr_Genpt_Ds->SetPoint(i,h_Genpt_Ds_fine->GetBinCenter(i),h_Genpt_Ds_fine->GetBinContent(i))	;
	}
		gr_Genpt_Ds->Write();

	// new DsOverD0
	TFile *f_CSAll_ph0=TFile::Open("./TheoryPredict/CS_CharmAll_ph0.root","READ");
	TFile *f_CSAll_ph20=TFile::Open("./TheoryPredict/CS_CharmAll_ph20.root","READ");



	TH1D *hGpt_D0_ph0=(TH1D*)f_CSAll_ph0->Get("hGpt_D0");
	TH1D *hGpt_Ds_ph0=(TH1D*)f_CSAll_ph0->Get("hGpt_Ds");

	TH1D *hGpt_D0_ph20=(TH1D*)f_CSAll_ph20->Get("hGpt_D0");
	TH1D *hGpt_Ds_ph20=(TH1D*)f_CSAll_ph20->Get("hGpt_Ds");

	hGpt_Ds_ph0->Divide(hGpt_D0_ph0);
	hGpt_Ds_ph20->Divide(hGpt_D0_ph20);

	TGraph *gr_DsD0=new TGraph();
	gr_DsD0->SetName("gr_DsD0");
	for(int i =1; i<hGpt_Ds_ph0->GetNbinsX()+1;i++){
		if(hGpt_Ds_ph0->GetBinCenter(i)<20){
		gr_DsD0->SetPoint(i,hGpt_Ds_ph0->GetBinCenter(i),hGpt_Ds_ph0->GetBinContent(i));
		}else{
		cout<<"using ph20 ,hGpt_Ds_ph20->GetBinCenter() = "<<hGpt_Ds_ph20->GetBinCenter(i)<<endl;
		gr_DsD0->SetPoint(i,hGpt_Ds_ph20->GetBinCenter(i),hGpt_Ds_ph20->GetBinContent(i));
		}

	}

	f_out->cd();
	gr_DsD0->Write();

	// compare the y<1 difference , since the weight is applied for Gy<1 entry... 
	// the default Gen filter is not y<1 , the filter efficiency should be larger than the actual value for Gy<1 cut


}
