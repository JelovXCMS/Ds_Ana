// need to use original pythia B->Ds / B->D pt shape as additional weighting factor

// Goal calculate CS of NonPrompt Ds & dN/dpt for PbPb , from Hao's BtoD result (and maybe FONLL), change it to my ana_bins
// Inputs : 
//		1. BtoD input from Hao, (and errors), may also need FONLL input
//    2. FR_B & BtoD (and errors) -> DtoDsScale
//    3. BtoD shape & BtoDs shape difference (different weight procedures)



#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <iostream>
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"



using namespace std;


// bool verbose=true;


TH1D *h_applyWeight_fun(TH1D *h_ori , TH1D* h_weight, TString WeightName="", Int_t WeightMethod=1 ,const int nbin=nbin_pt_pp, double *bins=bins_pt_pp)
{
	if(WeightName==""){
		WeightName="weighted";	
	}
	TString HisName=Form("%s_%s",h_ori->GetName(),WeightName.Data());
	TH1D *h_out=new TH1D(HisName.Data(),HisName.Data(),nbin,bins); h_out->Sumw2();
	Double_t err=0;
	for(int i=0; i<nbin; i++){
		int bin_weight=h_weight->FindBin(h_ori->GetXaxis()->GetBinCenter(i+1));
		cout<<"bin_weight = "<<bin_weight<<endl;
		// h_out->SetBinContent(i+1, h_ori->GetBinContent(i+1)*h_weight->GetBinContent(i+1));
		// err=ErrorPro_AtimesB(h_ori->GetBinError(i+1), h_weight->GetBinError(i+1),  h_ori->GetBinContent(i+1), h_weight->GetBinContent(i+1));
		h_out->SetBinContent(i+1, h_ori->GetBinContent(i+1)*h_weight->GetBinContent(bin_weight));
		err=ErrorPro_AtimesB(h_ori->GetBinError(i+1), h_weight->GetBinError(bin_weight),  h_ori->GetBinContent(i+1), h_weight->GetBinContent(bin_weight));
		h_out->SetBinError(i+1, err);
		if(verbose){ cout<<"weight fun bin "<<i+1<<" , bincontent = "<<h_out->GetBinContent(i+1)<<" ,err = "<<h_out->GetBinError(i+1)<<endl;}
	}

	return h_out;
}


TH1D *h_tranformToAnaBin_fun(TH1D *h_input,TString HisName="" ,const int nbin=nbin_pt_pp, double *bins=bins_pt_pp, int doMutiplyWidth=1){

	if(HisName==""){
		HisName=Form("%s_AnaBin",h_input->GetName());
	}
	TH1D *h_input_temp=(TH1D*)h_input->Clone(Form("%s_temp_fun",h_input->GetName()));
	if(doMutiplyWidth == 1){
	MutiplyBinWidth(h_input_temp); // why this does not work?
	}

	if(verbose){
	cout<<"check MutiplyBinWidth"<<endl;
	for(int i =1; i<h_input_temp->GetNbinsX();i++){
		cout<<"bin "<<i<<" , ori = "<<h_input->GetBinContent(i)<<" , mutiply = "<<h_input_temp->GetBinContent(i)<<endl;
	}
	}
	
	TH1D *h_output=new TH1D(HisName.Data(),HisName.Data(),nbin,bins); h_output->Sumw2(); 
	double binIntegralError=0;
	double binIntegral=0;

	double binLow=0;
	double binHigh=0;

	for(int i=0; i<nbin; i++){
		binLow=bins[i];
		binHigh=bins[i+1];

		binIntegral=0;
		binIntegralError=0;

		binIntegral= h_input_temp->IntegralAndError(h_input_temp->FindBin(binLow+0.1),h_input_temp->FindBin(binHigh-0.1), binIntegralError);
		h_output->SetBinContent(i+1,binIntegral);	
		h_output->SetBinError(i+1,binIntegralError);	

		if(verbose){
		cout<<"ibin Low = "<< h_input_temp->FindBin(binLow+0.1)<<" , ibin High = "<<h_input_temp->FindBin(binHigh-0.1)<<endl;
		cout<<"binIntegral = "<< binIntegral<<endl;
		}

		}

	if(doMutiplyWidth==1){
	divideBinWidth(h_output);
	}

	return h_output;

}


int DtoDs(Int_t isPbPb=0){

	verbose=true;
	InitStyle();
	initParameter();	

	const int nbin_pt=nbin_pt_pp;
	double *bins_pt=bins_pt_pp;

	// double bins_pt_old[]={2,3,4,5,6,7,8,10,12,20,40,60};
	// const int nbins_pt_old=sizeof(bins_pt_old)/sizeof(bins_pt_old[0]);

	// constant
	////-- 1. Bfraction & BR D ,Ds

	double Fr_B0=0.389; // value from https://arxiv.org/pdf/1306.3663.pdf LHCb 7TeV pp 
	double Fr_Bp=0.381; 
	double Fr_Bs=0.105;

	double Fr_B0_Z=0.412; // value from pdg , Z->bb, 1.96TeV
	double Fr_Bp_Z=0.412; 
	double Fr_Bs_Z=0.088;

	double Fr_B0_p=0.34; // value from pdg , Z->bb, 1.96TeV
	double Fr_Bp_p=0.34; 
	double Fr_Bs_p=0.101;

	double Fr_B0_Err=0;
	double Fr_Bp_Err=0;
	double Fr_Bs_Err=0;

	// need to add other Fr B for systematic here

	////-- 2. Br to D(s) from PDG
	double BR_B0toD=0.555;
	double BR_BptoD=0.876;
	double BR_BstoD=0.00213; 

	double BR_B0toDs=0.103;
	double BR_BptoDs=0.09;
	double BR_BstoDs=0.93; 

	double delta_BR_B0toD=0.043;
	double delta_BR_BptoD=0.047;
	double delta_BR_BstoD=0.00043; 

	double delta_BR_B0toDs=0.021;
	double delta_BR_BptoDs=0.018;
	double delta_BR_BstoDs=0.25; 


	double deltaP_BR_B0toD=0.043;
	double deltaP_BR_BptoD=0.047;
	double deltaP_BR_BstoD=0.00043; 

	double deltaP_BR_B0toDs=0.021;
	double deltaP_BR_BptoDs=0.018;
	double deltaP_BR_BstoDs=0.25; 

	double deltaM_BR_B0toD=0.043;
	double deltaM_BR_BptoD=0.047;
	double deltaM_BR_BstoD=0.00043; 

	double deltaM_BR_B0toDs=0.018;
	double deltaM_BR_BptoDs=0.0162;
	double deltaM_BR_BstoDs=0.25; 



  double BR_B0toD_ErrRel=delta_BR_B0toD/BR_B0toD;
  double BR_BptoD_ErrRel=delta_BR_BptoD/BR_BptoD;
  double BR_BstoD_ErrRel=delta_BR_BstoD/BR_BstoD;

  double BR_B0toDs_ErrRel=delta_BR_B0toDs/BR_B0toDs;
  double BR_BptoDs_ErrRel=delta_BR_BptoDs/BR_BptoDs;
  double BR_BstoDs_ErrRel=delta_BR_BstoDs/BR_BstoDs;


	double BR_B0toDmax=BR_B0toD+deltaP_BR_B0toD;
	double BR_BptoDmax=BR_BptoD+deltaP_BR_BptoD;
	double BR_BstoDmax=BR_BstoD+deltaP_BR_BstoD; 

	double BR_B0toDsmax=BR_B0toDs+deltaP_BR_B0toDs;
	double BR_BptoDsmax=BR_BptoDs+deltaP_BR_BptoDs;
	double BR_BstoDsmax=BR_BstoDs+deltaP_BR_BstoDs; 


	double BR_B0toDmin=BR_B0toD-deltaM_BR_B0toD;
	double BR_BptoDmin=BR_BptoD-deltaM_BR_BptoD;
	double BR_BstoDmin=BR_BstoD-deltaM_BR_BstoD; 

	double BR_B0toDsmin=BR_B0toDs-deltaM_BR_B0toDs;
	double BR_BptoDsmin=BR_BptoDs-deltaM_BR_BptoDs;
	double BR_BstoDsmin=BR_BstoDs-deltaM_BR_BstoDs; 




	// add BR Err for systematic here


	double DtoDsScale=(Fr_B0*BR_B0toDs+Fr_Bp*BR_BptoDs+Fr_Bs*BR_BstoDs)/(Fr_B0*BR_B0toD+Fr_Bp*BR_BptoD+Fr_Bs*BR_BstoD); 

	cout<<"DtoDsScale = " <<DtoDsScale<<endl;

	double NumErr=sqrt(pow(Fr_B0*delta_BR_B0toDs,2)+pow(Fr_Bp*delta_BR_BptoDs,2)+pow(Fr_Bs*delta_BR_BstoDs,2));
	double DenoErr=sqrt(pow(Fr_B0*delta_BR_B0toD,2)+pow(Fr_Bp*delta_BR_BptoD,2)+pow(Fr_Bs*delta_BR_BstoD,2));
	double DtoDsScale_BrErrRel=sqrt( pow(NumErr/Fr_B0*BR_B0toDs+Fr_Bp*BR_BptoDs+Fr_Bs*BR_BstoDs,2) + pow( DenoErr/Fr_B0*BR_B0toD+Fr_Bp*BR_BptoD+Fr_Bs*BR_BstoD,2) );

	cout<<"DtoDsScale_BrErrRel = "<<DtoDsScale_BrErrRel<<endl;


	double DtoDsScale_Fr_Z=(Fr_B0_Z*BR_B0toDs+Fr_Bp_Z*BR_BptoDs+Fr_Bs_Z*BR_BstoDs)/(Fr_B0_Z*BR_B0toD+Fr_Bp_Z*BR_BptoD+Fr_Bs_Z*BR_BstoD);
	double DtoDsScale_Fr_p=(Fr_B0_p*BR_B0toDs+Fr_Bp_p*BR_BptoDs+Fr_Bs_p*BR_BstoDs)/(Fr_B0_p*BR_B0toD+Fr_Bp_p*BR_BptoD+Fr_Bs_p*BR_BstoD);

	cout<<"DtoDsScale_Fr_Z = "<<DtoDsScale_Fr_Z<<" ,DtoDsScale_Fr_p = "<<DtoDsScale_Fr_p<<endl;

	double DtoDsScaleMax=(Fr_B0_p*BR_B0toDsmax+Fr_Bp_p*BR_BptoDsmax+Fr_Bs_p*BR_BstoDsmax)/(Fr_B0_p*BR_B0toDmin+Fr_Bp_p*BR_BptoDmin+Fr_Bs_p*BR_BstoDmin);
	double DtoDsScaleMin=(Fr_B0_Z*BR_B0toDsmin+Fr_Bp_Z*BR_BptoDsmin+Fr_Bs_Z*BR_BstoDsmin)/(Fr_B0_Z*BR_B0toDmax+Fr_Bp_Z*BR_BptoDmax+Fr_Bs_Z*BR_BstoDmax);

	cout<<"DtoDsScaleMax = "<<DtoDsScaleMax<<endl;
	cout<<"DtoDsScaleMin = "<<DtoDsScaleMin<<endl;

	double DtoDsScale_BR_B0toDsmax= (Fr_B0*BR_B0toDsmax+Fr_Bp*BR_BptoDs+Fr_Bs*BR_BstoDs)/(Fr_B0*BR_B0toD+Fr_Bp*BR_BptoD+Fr_Bs*BR_BstoD);
	double DtoDsScale_BR_B0toDsmin= (Fr_B0*BR_B0toDsmin+Fr_Bp*BR_BptoDs+Fr_Bs*BR_BstoDs)/(Fr_B0*BR_B0toD+Fr_Bp*BR_BptoD+Fr_Bs*BR_BstoD);

	double DtoDsScale_BR_BptoDsmax= (Fr_B0*BR_B0toDs+Fr_Bp*BR_BptoDsmax+Fr_Bs*BR_BstoDs)/(Fr_B0*BR_B0toD+Fr_Bp*BR_BptoD+Fr_Bs*BR_BstoD);
	double DtoDsScale_BR_BptoDsmin= (Fr_B0*BR_B0toDs+Fr_Bp*BR_BptoDsmin+Fr_Bs*BR_BstoDs)/(Fr_B0*BR_B0toD+Fr_Bp*BR_BptoD+Fr_Bs*BR_BstoD);

	double DtoDsScale_BR_BstoDsmax= (Fr_B0*BR_B0toDs+Fr_Bp*BR_BptoDs+Fr_Bs*BR_BstoDsmax)/(Fr_B0*BR_B0toD+Fr_Bp*BR_BptoD+Fr_Bs*BR_BstoD);
	double DtoDsScale_BR_BstoDsmin= (Fr_B0*BR_B0toDs+Fr_Bp*BR_BptoDs+Fr_Bs*BR_BstoDsmin)/(Fr_B0*BR_B0toD+Fr_Bp*BR_BptoD+Fr_Bs*BR_BstoD);

	double DtoDsScale_BR_B0toDmax= (Fr_B0*BR_B0toDs+Fr_Bp*BR_BptoDs+Fr_Bs*BR_BstoDs)/(Fr_B0*BR_B0toDmax+Fr_Bp*BR_BptoD+Fr_Bs*BR_BstoD);
	double DtoDsScale_BR_B0toDmin= (Fr_B0*BR_B0toDs+Fr_Bp*BR_BptoDs+Fr_Bs*BR_BstoDs)/(Fr_B0*BR_B0toDmin+Fr_Bp*BR_BptoD+Fr_Bs*BR_BstoD);

	double DtoDsScale_BR_BptoDmax= (Fr_B0*BR_B0toDs+Fr_Bp*BR_BptoDs+Fr_Bs*BR_BstoDs)/(Fr_B0*BR_B0toD+Fr_Bp*BR_BptoDmax+Fr_Bs*BR_BstoD);
	double DtoDsScale_BR_BptoDmin= (Fr_B0*BR_B0toDs+Fr_Bp*BR_BptoDs+Fr_Bs*BR_BstoDs)/(Fr_B0*BR_B0toD+Fr_Bp*BR_BptoDmin+Fr_Bs*BR_BstoD);

	double DtoDsScale_BR_BstoDmax= (Fr_B0*BR_B0toDs+Fr_Bp*BR_BptoDs+Fr_Bs*BR_BstoDs)/(Fr_B0*BR_B0toD+Fr_Bp*BR_BptoD+Fr_Bs*BR_BstoDmax);
	double DtoDsScale_BR_BstoDmin= (Fr_B0*BR_B0toDs+Fr_Bp*BR_BptoDs+Fr_Bs*BR_BstoDs)/(Fr_B0*BR_B0toD+Fr_Bp*BR_BptoD+Fr_Bs*BR_BstoDmin);

	// return 1;
	// max/min DtoDsScale scan from variation of 1. FR , 2. BR_B

//return 1;
	cout<<"check1"<<endl;
////-- 3.  reading files

		TFile *f_BtoDs=TFile::Open("./output/BtoDs.root","RECREATE");

		// 3.0 FONLL prediction 

		// TFile *f_DsMC_FONLL_pp=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_MergeWeight/root_output/DsMinTree_MC_GenSampleMerge_pp_NonPrompt_phi.root","READ");
		// TFile *f_DsMC_FONLL_PbPb=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_MergeWeight/root_output/DsMinTree_MC_GenSampleMerge_PbPb_NonPrompt_phi.root","READ");


		TFile *f_DsMC_FONLL_pp=TFile::Open("/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/DsMinTree_MC_GenSampleMerge_pp_NonPrompt_phi.root","READ");
		TFile *f_DsMC_FONLL_PbPb=TFile::Open("/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/DsMinTree_MC_GenSampleMerge_PbPb_NonPrompt_phi.root","READ");

		TTree *t_DsMC_FONLL_pp=(TTree*)f_DsMC_FONLL_pp->Get("ntGen");
		TTree *t_DsMC_FONLL_PbPb=(TTree*)f_DsMC_FONLL_PbPb->Get("ntGen");

		f_BtoDs->cd();	

		TH1D *hBtoDsCrossSectionPP_FONLL= new TH1D("hBtoDsCrossSectionPP_FONLL","hBtoDsCrossSectionPP_FONLL",nbin_pt_pp,bins_pt_pp); hBtoDsCrossSectionPP_FONLL->Sumw2();
	  TH1D *hBtoDsdNdPtPbPb_FONLL = new TH1D("hBtoDsdNdPtPbPb_FONLL","hBtoDsdNdPtPbPb_FONLL",nbin_pt_PbPb3,bins_pt_PbPb3); hBtoDsdNdPtPbPb_FONLL->Sumw2();

		t_DsMC_FONLL_pp->Project("hBtoDsCrossSectionPP_FONLL","Gpt","(GSignalType==1 && GcollisionId==0 && TMath::Abs(Gy)<1 && GBAncestorpt>0)*(weight*GptSampleWeight*GenFONLLWeight)");
		t_DsMC_FONLL_PbPb->Project("hBtoDsdNdPtPbPb_FONLL","Gpt","(GSignalType==1 && GcollisionId==0 && TMath::Abs(Gy)<1 && GBAncestorpt>0)*(weight*GptSampleWeight*PbPbVzWeight*Ncoll*GenFONLLRaaWeight)");

	f_BtoDs->cd();
	// need to divide by binwidth
	
	divideBinWidth(hBtoDsCrossSectionPP_FONLL);
	divideBinWidth(hBtoDsdNdPtPbPb_FONLL);

	hBtoDsCrossSectionPP_FONLL->Write("",TObject::kOverwrite);
	hBtoDsdNdPtPbPb_FONLL->Write("",TObject::kOverwrite);

	// continue here	


		// 3.1 From Data B->D

	TFile *f_BtoD_CS=TFile::Open("./ori_root/BtoDCrossSection.root","READ");
	TFile *f_BtoD_PbPb=TFile::Open("./ori_root/BtoDdNdPt.root","READ");

	f_BtoDs->cd();	

	TH1D *hPromptDCrossSectionPP = (TH1D*)f_BtoD_CS->Get("hPromptDCrossSectionPP");
	TH1D *hBtoDCrossSectionPP = (TH1D*)f_BtoD_CS->Get("hBtoDCrossSectionPP");
	TH1D *hSysBtoDCrossSectionPP =(TH1D*)f_BtoD_CS->Get("hSysBtoDCrossSectionPP");

	TH1D *hPromptDdNdPtPbPb=(TH1D*)f_BtoD_PbPb->Get("hPromptDdNdPtPbPb");
	TH1D *hBtoDdNdPtPbPb=(TH1D*)f_BtoD_PbPb->Get("hBtoDdNdPtPbPb");
	TH1D *hSysBtoDdNdPtPbPb=(TH1D*)f_BtoD_PbPb->Get("hSysBtoDdNdPtPbPb");

	////-- 4. apply DtoDs scale

	TH1D *hBtoDsCrossSectionPP = (TH1D*)hBtoDCrossSectionPP->Clone("hBtoDsCrossSectionPP");
	hBtoDsCrossSectionPP->Scale(DtoDsScale);
	TH1D *hSysBtoDsCrossSectionPP =(TH1D*)hSysBtoDCrossSectionPP->Clone("hSysBtoDsCrossSectionPP");
	hSysBtoDsCrossSectionPP->Scale(DtoDsScale);
	// TH1D *hBtoDsCrossSectionPP_totalErr=(TH1D*)hBtoDsCrossSectionPP->Clone("hBtoDsCrossSectionPP_totalErr");

	TH1D *hBtoDsdNdPtPbPb=(TH1D*)hBtoDdNdPtPbPb->Clone("hBtoDsdNdPtPbPb");
	hBtoDsdNdPtPbPb->Scale(DtoDsScale);
	TH1D *hSysBtoDsdNdPtPbPb=(TH1D*)hSysBtoDdNdPtPbPb->Clone("hSysBtoDsdNdPtPbPb");
	hSysBtoDsdNdPtPbPb->Scale(DtoDsScale);

	cout<<"check2"<<endl;
	// test DtoDsScale
	plotylog=true;
	Draw({hBtoDCrossSectionPP,hBtoDsCrossSectionPP});
	Draw({hBtoDdNdPtPbPb,hBtoDsdNdPtPbPb});

	////-- 5. change to ana bin

	TH1D *hBtoDsCrossSectionPP_AnaBin=h_tranformToAnaBin_fun(hBtoDsCrossSectionPP);
	TH1D *hSysBtoDsCrossSectionPP_AnaBin=h_tranformToAnaBin_fun(hSysBtoDsCrossSectionPP);
	TH1D *hBtoDCrossSectionPP_AnaBin=h_tranformToAnaBin_fun(hBtoDCrossSectionPP);
	TH1D *hSysBtoDCrossSectionPP_AnaBin=h_tranformToAnaBin_fun(hSysBtoDCrossSectionPP);
	TH1D *hPromptDCrossSectionPP_AnaBin=h_tranformToAnaBin_fun(hPromptDCrossSectionPP);

	TH1D *hBtoDsdNdPtPbPb_AnaBin=h_tranformToAnaBin_fun(hBtoDsdNdPtPbPb,"",nbin_pt_PbPb3,bins_pt_PbPb3);
	TH1D *hBtoDdNdPtPbPb_AnaBin=h_tranformToAnaBin_fun(hBtoDdNdPtPbPb,"",nbin_pt_PbPb3,bins_pt_PbPb3);
	TH1D *hSysBtoDsdNdPtPbPb_AnaBin=h_tranformToAnaBin_fun(hSysBtoDsdNdPtPbPb,"",nbin_pt_PbPb3,bins_pt_PbPb3);
	TH1D *hSysBtoDdNdPtPbPb_AnaBin=h_tranformToAnaBin_fun(hSysBtoDdNdPtPbPb,"",nbin_pt_PbPb3,bins_pt_PbPb3);
	TH1D *hPromptDdNdPtPbPb_Anabin=h_tranformToAnaBin_fun(hPromptDdNdPtPbPb,"",nbin_pt_PbPb3,bins_pt_PbPb3);

	Draw({hBtoDsCrossSectionPP_AnaBin,hBtoDCrossSectionPP_AnaBin});
	Draw({hBtoDsdNdPtPbPb_AnaBin,hBtoDdNdPtPbPb_AnaBin});

	plotylog=false;

	TCanvas *c_pp_AnaBintest=new TCanvas("c_pp_AnaBintest","c_pp_AnaBintest",800,800);
	c_pp_AnaBintest->cd();
	hBtoDsCrossSectionPP->GetXaxis()->SetRangeUser(1,40);
	hBtoDsCrossSectionPP->SetLineColor(1);
	hBtoDsCrossSectionPP->SetMarkerColor(1);
	hBtoDsCrossSectionPP->Draw("SAME");
	hBtoDsCrossSectionPP_AnaBin->SetLineColor(2);
	hBtoDsCrossSectionPP_AnaBin->SetMarkerColor(2);
	hBtoDsCrossSectionPP_AnaBin->Draw("SAME");
	gPad->SetLogy();
	gPad->SetLogx();
	gPad->BuildLegend();
	c_pp_AnaBintest->SaveAs("plots/pdf/pp_AnaBintest.pdf");

	TCanvas *c_PbPb_AnaBintest=new TCanvas("c_PbPb_AnaBintest","c_PbPb_AnaBintest",800,800);
	c_PbPb_AnaBintest->cd();
	hBtoDsdNdPtPbPb->GetXaxis()->SetRangeUser(1,40);
	hBtoDsdNdPtPbPb->SetLineColor(1);
	hBtoDsdNdPtPbPb->SetMarkerColor(1);
	hBtoDsdNdPtPbPb->Draw("SAME");
	hBtoDsdNdPtPbPb_AnaBin->SetLineColor(2);
	hBtoDsdNdPtPbPb_AnaBin->SetMarkerColor(2);
	hBtoDsdNdPtPbPb_AnaBin->Draw("SAME");
	gPad->SetLogy();
	gPad->SetLogx();
	gPad->BuildLegend();
	c_PbPb_AnaBintest->SaveAs("plots/pdf/PbPb_AnaBintest.pdf");


	////-- 6. DtoDs shape reweight
	////  Method 1. directly reweight from DtoDs difference
	////  Method 2. using some pre-defined 


	// TFile *f_DtoDs_Weight=TFile::Open("NonPrompt_DsGpt_phikkpi_pp.root");


	//// the weight calculation might need to put into a seperate code?
	TFile *f_D0=TFile::Open("./output/DptSampleWeight_pp_NonPrompt.root");
  TH1D *h_D0Gpt_weighted=(TH1D*)f_D0->Get("h_Gpt_fromHaoWeight");
  h_D0Gpt_weighted->Scale(1/h_D0Gpt_weighted->Integral());
  TH1D *h_D0Gpt_weightedFineBin=(TH1D*)f_D0->Get("h_Gpt_fromHaoWeightFineBin");
  h_D0Gpt_weightedFineBin->Scale(1/h_D0Gpt_weightedFineBin->Integral());

/*
  TFile *f_Ds=TFile::Open("./output/DsptSampleWeight_pp_NonPrompt_phi.root");
  TH1D *h_DsGpt_weighted=(TH1D*)f_Ds->Get("h_Gpt_weighted2");
  h_DsGpt_weighted->Scale(1/h_DsGpt_weighted->Integral());
  TH1D *h_DsGpt_weightedFineBin=(TH1D*)f_Ds->Get("h_Gpt_weighted2FineBin");
  h_DsGpt_weightedFineBin->Scale(1/h_DsGpt_weightedFineBin->Integral());
*/
	// using weighted merge MC sample

	TH1D *h_DsGpt_weighted=new TH1D("h_DsGpt_weighted","h_DsGpt_weighted",nbin_pt_Norm,bins_pt_Norm); h_DsGpt_weighted->Sumw2();
	TH1D *h_DsGpt_weightedFineBin=new TH1D("h_DsGpt_weightedFineBin","h_DsGpt_weightedFineBin",100,0,50); h_DsGpt_weightedFineBin->Sumw2();
	t_DsMC_FONLL_pp->Project("h_DsGpt_weighted","Gpt","(GSignalType==1 && GcollisionId==0 && TMath::Abs(Gy)<1 && GBAncestorpt>0)*(weight*GptSampleWeight)");
	h_DsGpt_weighted->Scale(1/h_DsGpt_weighted->Integral());


	t_DsMC_FONLL_pp->Project("h_DsGpt_weightedFineBin","Gpt","(GSignalType==1 && GcollisionId==0 && TMath::Abs(Gy)<1 && GBAncestorpt>0)*(weight*GptSampleWeight)");
	h_DsGpt_weightedFineBin->Scale(1/h_DsGpt_weightedFineBin->Integral());

  TH1D *h_DsD0Ratio_pythia=(TH1D*)h_DsGpt_weighted->Clone("h_DsD0Ratio_pythia"); // this is cover full range pt bin, start from 0
  h_DsD0Ratio_pythia->SetTitle("h_DsD0Ratio_pythia");
  h_DsD0Ratio_pythia->Divide(h_D0Gpt_weighted);


	TCanvas *c_DsD0Ratio_pythia=new TCanvas("c_DsD0Ratio_pythia","c_DsD0Ratio_pythia",800,800);
	c_DsD0Ratio_pythia->cd();
	h_DsD0Ratio_pythia->Draw();
	c_DsD0Ratio_pythia->SaveAs("plots/pdf/c_DsD0Ratio_pythia.pdf");

	TH1D *h_D0toDs_pythiaWeight_AnaBin=h_tranformToAnaBin_fun(h_DsD0Ratio_pythia,"",nbin_pt_pp,bins_pt_pp,0);

	gStyle->SetOptStat(0);

	TCanvas *c_D0DsSpectra= new TCanvas("c_D0DsSpectra","c_D0DsSpectra",800,800);
	c_D0DsSpectra->cd();
	gPad->SetLogx();
	gPad->SetLogy();
	h_DsGpt_weightedFineBin->SetTitle("");
	h_DsGpt_weightedFineBin->GetXaxis()->SetTitle("p_{T} (GeV)");
	h_DsGpt_weightedFineBin->GetXaxis()->CenterTitle();
	h_DsGpt_weightedFineBin->SetLineColor(2);
	h_DsGpt_weightedFineBin->SetMarkerColor(2);
	h_DsGpt_weightedFineBin->Draw("SAME");
	h_D0Gpt_weightedFineBin->SetLineColor(4);
	h_D0Gpt_weightedFineBin->SetMarkerColor(4);
	h_D0Gpt_weightedFineBin->Draw("SAME");

	TLegend *le_D0Ds=new TLegend(0.2,0.2,0.45,0.45,NULL,"brNDC");
	le_D0Ds->SetBorderSize(0);
	le_D0Ds->SetTextSize(0.04);
	le_D0Ds->AddEntry((TObject*)0,"PYTHIA 8","");
	le_D0Ds->AddEntry(h_DsGpt_weightedFineBin,"nonprompt D_{S}^{#pm}","lp");
	le_D0Ds->AddEntry(h_D0Gpt_weightedFineBin,"nonprompt D^{0}","lp");
	le_D0Ds->Draw("SAME");

	SavePlotDirs(c_D0DsSpectra,"D0DsSpectra",{"BtoDs"});

	TLatex *tl=new TLatex();

	TCanvas *c_weight=new TCanvas("c_weight","c_weight",800,800);
	c_weight->cd();
	gPad->SetLogx();
	h_D0toDs_pythiaWeight_AnaBin->SetTitle("");
	h_D0toDs_pythiaWeight_AnaBin->GetXaxis()->SetTitle("p_{T} (GeV)");
	h_D0toDs_pythiaWeight_AnaBin->GetXaxis()->CenterTitle();
	h_D0toDs_pythiaWeight_AnaBin->GetYaxis()->SetTitle("nonprompt Ds/ nonprompt D0 Ratio");
	h_D0toDs_pythiaWeight_AnaBin->GetYaxis()->CenterTitle();
	h_D0toDs_pythiaWeight_AnaBin->SetLineColor(2);
	h_D0toDs_pythiaWeight_AnaBin->SetMarkerColor(2);
	h_D0toDs_pythiaWeight_AnaBin->Draw();

	tl->DrawLatexNDC(0.6,0.3,"PYTHIA 8");

	c_weight->SaveAs("plots/pdf/c_weight.pdf");

	SavePlotDirs(c_weight,"DtoDsScale",{"BtoDs"});

	cout<<"check3"<<endl;


// TH1D *h_applyWeight_fun(TH1D *h_ori , TH1D* h_weight, TString WeightName="", Int_t WeightMethod=1 ,const int nbin=nbin_pt_pp, double *bin=bins_pt_pp)

	TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight=h_applyWeight_fun(hBtoDsCrossSectionPP_AnaBin, h_D0toDs_pythiaWeight_AnaBin,"pythiaWeight");
	TH1D *hSysBtoDsCrossSectionPP_AnaBin_pythiaWeight=h_applyWeight_fun(hSysBtoDsCrossSectionPP_AnaBin, h_D0toDs_pythiaWeight_AnaBin,"pythiaWeight");

// need to check PbPb binning 
	TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight=h_applyWeight_fun(hBtoDsdNdPtPbPb_AnaBin, h_D0toDs_pythiaWeight_AnaBin,"pythiaWeight",1,nbin_pt_PbPb3,bins_pt_PbPb3);
	TH1D *hSysBtoDsdNdPtPbPb_AnaBin_pythiaWeight=h_applyWeight_fun(hSysBtoDsdNdPtPbPb_AnaBin, h_D0toDs_pythiaWeight_AnaBin,"pythiaWeight",1,nbin_pt_PbPb3,bins_pt_PbPb3);

	TCanvas *c_testweight_pp=new TCanvas("c_testweight_pp","c_testweight_pp",800,800);
	c_testweight_pp->cd();
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Draw("PLC PMC"); // only work after root 6.09
	hBtoDsCrossSectionPP_AnaBin->Draw("SAME PLC PMC");
	gPad->BuildLegend();
	gPad->SetLogy();
	gPad->SetLogx();
	c_testweight_pp->SaveAs("plots/pdf/testweight_pp.pdf");

	TCanvas *c_testweight_PbPb=new TCanvas("c_testweight_PbPb","c_testweight_PbPb",800,800);
	c_testweight_PbPb->cd();
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Draw("PLC PMC"); // only work after root 6.09
	hBtoDsCrossSectionPP_AnaBin->Draw("SAME PLC PMC");
	gPad->BuildLegend();
	gPad->SetLogy();
	gPad->SetLogx();
	c_testweight_PbPb->SaveAs("plots/pdf/testweight_PbPb.pdf");


// add stat error and syst error to one
	TH1D *hTotalErrBtoDsCrossSectionPP_AnaBin_pythiaWeight=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hTotalErrBtoDsCrossSectionPP_AnaBin_pythiaWeight");
	TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight_Errup=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_Errup");
	TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight_Errdown=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_Errdown");

	for(int i=0; i<nbin_pt_pp;i++){
		hTotalErrBtoDsCrossSectionPP_AnaBin_pythiaWeight->SetBinError(i+1, sqrt(pow(hBtoDsCrossSectionPP_AnaBin_pythiaWeight->GetBinError(i+1),2) + pow(hSysBtoDsCrossSectionPP_AnaBin_pythiaWeight->GetBinError(i+1),2)  ));
		hBtoDsCrossSectionPP_AnaBin_pythiaWeight_Errup->SetBinContent(i+1, hBtoDsCrossSectionPP_AnaBin_pythiaWeight_Errup->GetBinContent(i+1)+ hTotalErrBtoDsCrossSectionPP_AnaBin_pythiaWeight->GetBinError(i+1));
		hBtoDsCrossSectionPP_AnaBin_pythiaWeight_Errdown->SetBinContent(i+1, hBtoDsCrossSectionPP_AnaBin_pythiaWeight_Errdown->GetBinContent(i+1) - hTotalErrBtoDsCrossSectionPP_AnaBin_pythiaWeight->GetBinError(i+1));

	}

	TH1D *hTotalErrBtoDsdNdPtPbPb_AnaBin_pythiaWeight=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hTotalErrBtoDsdNdPtPbPb_AnaBin_pythiaWeight");
	TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Errup=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Errup");
	TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Errdown=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Errdown");
	for(int i=0; i<nbin_pt_PbPb3;i++){
		hTotalErrBtoDsdNdPtPbPb_AnaBin_pythiaWeight->SetBinError(i+1, sqrt(pow(hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->GetBinError(i+1),2) + pow(hSysBtoDsdNdPtPbPb_AnaBin_pythiaWeight->GetBinError(i+1),2)  ))	;
		hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Errup->SetBinContent(i+1, hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Errup->GetBinContent(i+1)+ hTotalErrBtoDsdNdPtPbPb_AnaBin_pythiaWeight->GetBinError(i+1));
		hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Errdown->SetBinContent(i+1, hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Errup->GetBinContent(i+1) -  hTotalErrBtoDsdNdPtPbPb_AnaBin_pythiaWeight->GetBinError(i+1));

	}

	TCanvas *c_totalErr_pp =new TCanvas("c_totalErr_pp","c_totalErr_pp",800,800);
	c_totalErr_pp->cd();
	gPad->SetLogy();
	gPad->SetLogx();
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Draw("PLC PMC");
	hSysBtoDsCrossSectionPP_AnaBin_pythiaWeight->Draw("SAME PLC PMC");
	hTotalErrBtoDsCrossSectionPP_AnaBin_pythiaWeight->Draw("SAME PLC PMC");
	gPad->BuildLegend();
	c_totalErr_pp->SaveAs("plots/pdf/totalErr_pp.pdf");


	TCanvas *c_totalErr_PbPb =new TCanvas("c_totalErr_PbPb","c_totalErr_PbPb",800,800);
	c_totalErr_PbPb->cd();
	gPad->SetLogy();
	gPad->SetLogx();
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Draw("PLC PMC");
	hSysBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Draw("SAME PLC PMC");
	hTotalErrBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Draw("SAME PLC PMC");
	gPad->BuildLegend();
	c_totalErr_PbPb->SaveAs("plots/pdf/totalErr_PbPb.pdf");

	TCanvas *c_BtoDs_pp=new TCanvas("c_BtoDs_pp","c_BtoDs_pp",800,800);
	c_BtoDs_pp->cd();
	gPad->SetLogx();
	gPad->SetLogy();
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight->SetTitle("");
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight->GetXaxis()->SetTitle("p_{T} (GeV)");
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight->GetXaxis()->CenterTitle();
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight->GetYaxis()->SetTitle("nonprompt D_{S} d#sigma/dp_{T} (pb GeV^{-1})");
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight->GetYaxis()->CenterTitle();
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Draw("SAME");

	TLegend *le_BtoDs_pp=new TLegend(0.68,0.7,0.85,0.85,NULL,"brNDC");
	le_BtoDs_pp->SetBorderSize(0);
	le_BtoDs_pp->SetTextSize(0.04);
	le_BtoDs_pp->AddEntry((TObject*)0,"pp","");
	le_BtoDs_pp->Draw("SAME");


	SavePlotDirs(c_BtoDs_pp,"BtoDs_pp",{"BtoDs"});

	TCanvas *c_BtoDs_PbPb=new TCanvas("c_BtoDs_PbPb","c_BtoDs_PbPb",800,800);
	c_BtoDs_PbPb->cd();
	gPad->SetLogx();
	gPad->SetLogy();
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->SetTitle("");
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->GetXaxis()->SetTitle("p_{T} (GeV)");
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->GetXaxis()->CenterTitle();
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->GetYaxis()->SetTitle("nonprompt D_{S} dN/dp_{T}");
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->GetYaxis()->CenterTitle();
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Draw("SAME");

	TLegend *le_BtoDs_PbPb=new TLegend(0.68,0.7,0.85,0.85,NULL,"brNDC");
	le_BtoDs_PbPb->SetBorderSize(0);
	le_BtoDs_PbPb->SetTextSize(0.04);
	le_BtoDs_PbPb->AddEntry((TObject*)0,"PbPb","");
	le_BtoDs_PbPb->Draw("SAME");


	SavePlotDirs(c_BtoDs_PbPb,"BtoDs_PbPb",{"BtoDs"});




	////-- last writing to file

	f_BtoDs->cd();
	hPromptDCrossSectionPP_AnaBin->Write("",TObject::kOverwrite);
	hBtoDCrossSectionPP_AnaBin->Write("",TObject::kOverwrite);

	hPromptDdNdPtPbPb_Anabin->Write("",TObject::kOverwrite);
	hBtoDdNdPtPbPb_AnaBin->Write("",TObject::kOverwrite);

	// hBtoDsCrossSectionPP->Write("",TObject::kOverwrite);
	// hSysBtoDsCrossSectionPP->Write("",TObject::kOverwrite);

	hBtoDsCrossSectionPP_AnaBin->Write("",TObject::kOverwrite);
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Write("",TObject::kOverwrite);
	hSysBtoDsCrossSectionPP_AnaBin_pythiaWeight->Write("",TObject::kOverwrite);
	hSysBtoDsCrossSectionPP_AnaBin->Write("",TObject::kOverwrite);
	hTotalErrBtoDsCrossSectionPP_AnaBin_pythiaWeight->Write("",TObject::kOverwrite);

	hBtoDsdNdPtPbPb_AnaBin->Write("",TObject::kOverwrite);
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Write("",TObject::kOverwrite);
	hSysBtoDsdNdPtPbPb_AnaBin->Write("",TObject::kOverwrite);
	hSysBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Write("",TObject::kOverwrite);
	hTotalErrBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Write("",TObject::kOverwrite);


	////-- for systematics

	TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Raaup=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Raaup");
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Raaup->SetTitle("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Raaup");
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Raaup->Scale(1.5);
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Raaup->Write("",TObject::kOverwrite);

	TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Raadown=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Raadown");
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Raadown->SetTitle("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Raadown");
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Raadown->Scale(1/1.5);
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Raadown->Write("",TObject::kOverwrite);


	// FONLL shape with Data CS/dNdpt

	double  bins_pt_DataNorm[]={2,3,4,5,6,8,10,20,40,100};
	const int nbin_pt_DataNorm=sizeof(bins_pt_DataNorm)/sizeof(bins_pt_DataNorm[0])-1;

	TH1D *hBtoDsdNdPtPbPb_NormBin=h_tranformToAnaBin_fun(hBtoDsdNdPtPbPb,"NormBin",nbin_pt_DataNorm,bins_pt_DataNorm);
	TH1D *hBtoDsCrossSectionPP_NormBin=h_tranformToAnaBin_fun(hBtoDsCrossSectionPP,"NormBin",nbin_pt_DataNorm,bins_pt_DataNorm);
	MutiplyBinWidth(hBtoDsdNdPtPbPb_NormBin);
	MutiplyBinWidth(hBtoDsCrossSectionPP_NormBin);

	TH1D *hBtoDsCrossSectionPP_FONLL_NormBin= new TH1D("hBtoDsCrossSectionPP_FONLL_NormBin","hBtoDsCrossSectionPP_FONLL_NormBin",nbin_pt_DataNorm,bins_pt_DataNorm); hBtoDsCrossSectionPP_FONLL_NormBin->Sumw2();
	TH1D *hBtoDsdNdPtPbPb_FONLL_NormBin = new TH1D("hBtoDsdNdPtPbPb_FONLL_NormBin","hBtoDsdNdPtPbPb_FONLL_NormBin",nbin_pt_DataNorm,bins_pt_DataNorm); hBtoDsdNdPtPbPb_FONLL_NormBin->Sumw2();

	t_DsMC_FONLL_pp->Project("hBtoDsCrossSectionPP_FONLL_NormBin","Gpt","(GSignalType==1 && GcollisionId==0 && TMath::Abs(Gy)<1 && GBAncestorpt>0)*(weight*GptSampleWeight*GenFONLLWeight)");
	t_DsMC_FONLL_PbPb->Project("hBtoDsdNdPtPbPb_FONLL_NormBin","Gpt","(GSignalType==1 && GcollisionId==0 && TMath::Abs(Gy)<1 && GBAncestorpt>0)*(weight*GptSampleWeight*PbPbVzWeight*Ncoll*GenFONLLRaaWeight)");

	hBtoDsdNdPtPbPb_NormBin->Scale(1/hBtoDsdNdPtPbPb_NormBin->Integral());
	hBtoDsCrossSectionPP_NormBin->Scale(1/hBtoDsCrossSectionPP_NormBin->Integral());

	hBtoDsCrossSectionPP_FONLL_NormBin->Scale(1/hBtoDsCrossSectionPP_FONLL_NormBin->Integral());
	hBtoDsdNdPtPbPb_FONLL_NormBin->Scale(1/hBtoDsdNdPtPbPb_FONLL_NormBin->Integral());

	TCanvas *c_afterscale=new TCanvas("c_afterscale","c_afterscale",800,800);
	c_afterscale->Divide(2,2);
	c_afterscale->cd(1);
	hBtoDsdNdPtPbPb_NormBin->Draw();
	c_afterscale->cd(2);
	hBtoDsdNdPtPbPb_FONLL_NormBin->Draw();
  c_afterscale->cd(3);
	hBtoDsCrossSectionPP_NormBin->Draw();
	c_afterscale->cd(4);
	hBtoDsCrossSectionPP_FONLL_NormBin->Draw();


	TH1D *hBtoDsCrossSectionPP_FONLL_Ratio=(TH1D*)hBtoDsCrossSectionPP_FONLL_NormBin->Clone("hBtoDsCrossSectionPP_FONLL_Ratio");
	hBtoDsCrossSectionPP_FONLL_Ratio->SetTitle("hBtoDsCrossSectionPP_FONLL_Ratio");
	hBtoDsCrossSectionPP_FONLL_Ratio->Divide(hBtoDsCrossSectionPP_NormBin);

	TH1D *hBtoDsdNdPtPbPb_FONLL_Ratio=(TH1D*)hBtoDsdNdPtPbPb_FONLL_NormBin->Clone("hBtoDsdNdPtPbPb_FONLL_Ratio");
	hBtoDsdNdPtPbPb_FONLL_Ratio->SetTitle("hBtoDsdNdPtPbPb_FONLL_Ratio");
	hBtoDsdNdPtPbPb_FONLL_Ratio->Divide(hBtoDsdNdPtPbPb_FONLL_Ratio,hBtoDsdNdPtPbPb_NormBin);

/*
	for(int i =0; i<nbin_pt_Norm; i++){
		cout<<"bin "<<i<<" , hBtoDsdNdPtPbPb_FONLL_NormBin = "<<hBtoDsdNdPtPbPb_FONLL_NormBin->GetBinContent(i+1)<<" , hBtoDsdNdPtPbPb_NormBin = "<<hBtoDsdNdPtPbPb_NormBin->GetBinContent(i+1)<<endl;
		hBtoDsdNdPtPbPb_FONLL_Ratio->SetBinContent(i+1, hBtoDsdNdPtPbPb_FONLL_NormBin->GetBinContent(i+1)/hBtoDsdNdPtPbPb_NormBin->GetBinContent(i+1));
	}
	*/

	TCanvas *c_FONLLShape_Ratio=new TCanvas("c_FONLLShape_Ratio","c_FONLLShape_Ratio",1200,600);
	c_FONLLShape_Ratio->Divide(2,1);
	c_FONLLShape_Ratio->cd(1);
	hBtoDsCrossSectionPP_FONLL_Ratio->Draw();
	c_FONLLShape_Ratio->cd(2);
	hBtoDsdNdPtPbPb_FONLL_Ratio->Draw();

	TH1D *hBtoDsCrossSectionPP_FONLL_Ratio_AnaBin=h_tranformToAnaBin_fun(hBtoDsCrossSectionPP_FONLL_Ratio,"",nbin_pt_pp,bins_pt_pp,0);
	TH1D *hBtoDsdNdPtPbPb_FONLL_Ratio_AnaBin=h_tranformToAnaBin_fun(hBtoDsdNdPtPbPb_FONLL_Ratio,"",nbin_pt_PbPb3,bins_pt_PbPb3,0);

  TCanvas *c_FONLLShape_Ratio2=new TCanvas("c_FONLLShape_Ratio2","c_FONLLShape_Ratio2",1200,600);
  c_FONLLShape_Ratio2->Divide(2,1);
  c_FONLLShape_Ratio2->cd(1);
  hBtoDsCrossSectionPP_FONLL_Ratio_AnaBin->Draw();
  c_FONLLShape_Ratio2->cd(2);
  hBtoDsdNdPtPbPb_FONLL_Ratio_AnaBin->Draw();


	TH1D *hBtoDsCrossSectionPP_FONLLShape=h_applyWeight_fun(hBtoDsCrossSectionPP_AnaBin,hBtoDsCrossSectionPP_FONLL_Ratio_AnaBin,"FONLLShape");
	TH1D *hBtoDsdNdPtPbPb_FONLLShape=h_applyWeight_fun(hBtoDsdNdPtPbPb_AnaBin,hBtoDsdNdPtPbPb_FONLL_Ratio_AnaBin,"FONLLShape",1,nbin_pt_PbPb3,bins_pt_PbPb3);

	f_BtoDs->cd();
	hBtoDsCrossSectionPP_FONLLShape->Write("",TObject::kOverwrite);
	hBtoDsdNdPtPbPb_FONLLShape->Write("",TObject::kOverwrite);



	// end FONLL shape with Data CS/dNdpt

	f_BtoDs->cd();
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight_Errup->Write("",TObject::kOverwrite);
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight_Errdown->Write("",TObject::kOverwrite);

	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Errup->Write("",TObject::kOverwrite);
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Errdown->Write("",TObject::kOverwrite);

		//DtoDs scale	
	TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBRup=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBRup");
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBRup->Scale((1+DtoDsScale_BrErrRel));
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBRup->Write("",TObject::kOverwrite);

	TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBRdown=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBRdown");
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBRdown->Scale((1-DtoDsScale_BrErrRel));
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBRdown->Write("",TObject::kOverwrite);

	TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBRup=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBRup");
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBRup->Scale((1+DtoDsScale_BrErrRel));
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBRup->Write("",TObject::kOverwrite);

	TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBRdown=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBRdown");
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBRdown->Scale((1-DtoDsScale_BrErrRel));
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBRdown->Write("",TObject::kOverwrite);


	// Br up & down seperately
 
  TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDmax=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDmax");
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDmax->Scale((DtoDsScale_BR_B0toDmax/DtoDsScale));
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDmax->Write("",TObject::kOverwrite);

  TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDmin=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDmin");
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDmin->Scale((DtoDsScale_BR_B0toDmin/DtoDsScale));
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDmin->Write("",TObject::kOverwrite);

  TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDmax=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDmax");
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDmax->Scale((DtoDsScale_BR_BptoDmax/DtoDsScale));
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDmax->Write("",TObject::kOverwrite);

  TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDmin=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDmin");
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDmin->Scale((DtoDsScale_BR_BptoDmin/DtoDsScale));
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDmin->Write("",TObject::kOverwrite);

  TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDmax=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDmax");
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDmax->Scale((DtoDsScale_BR_BstoDmax/DtoDsScale));
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDmax->Write("",TObject::kOverwrite);

  TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDmin=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDmin");
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDmin->Scale((DtoDsScale_BR_BstoDmin/DtoDsScale));
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDmin->Write("",TObject::kOverwrite);


  TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDsmax=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDsmax");
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDsmax->Scale((DtoDsScale_BR_B0toDsmax/DtoDsScale));
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDsmax->Write("",TObject::kOverwrite);

  TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDsmin=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDsmin");
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDsmin->Scale((DtoDsScale_BR_B0toDsmin/DtoDsScale));
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDsmin->Write("",TObject::kOverwrite);

  TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDsmax=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDsmax");
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDsmax->Scale((DtoDsScale_BR_BptoDsmax/DtoDsScale));
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDsmax->Write("",TObject::kOverwrite);

  TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDsmin=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDsmin");
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDsmin->Scale((DtoDsScale_BR_BptoDsmin/DtoDsScale));
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDsmin->Write("",TObject::kOverwrite);

  TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDsmax=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDsmax");
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDsmax->Scale((DtoDsScale_BR_BstoDsmax/DtoDsScale));
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDsmax->Write("",TObject::kOverwrite);

  TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDsmin=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDsmin");
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDsmin->Scale((DtoDsScale_BR_BstoDsmin/DtoDsScale));
  hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDsmin->Write("",TObject::kOverwrite);


			// pbPb

 
  TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDmax=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDmax");
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDmax->Scale((DtoDsScale_BR_B0toDmax/DtoDsScale));
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDmax->Write("",TObject::kOverwrite);

  TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDmin=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDmin");
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDmin->Scale((DtoDsScale_BR_B0toDmin/DtoDsScale));
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDmin->Write("",TObject::kOverwrite);

  TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDmax=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDmax");
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDmax->Scale((DtoDsScale_BR_BptoDmax/DtoDsScale));
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDmax->Write("",TObject::kOverwrite);

  TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDmin=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDmin");
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDmin->Scale((DtoDsScale_BR_BptoDmin/DtoDsScale));
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDmin->Write("",TObject::kOverwrite);

  TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDmax=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDmax");
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDmax->Scale((DtoDsScale_BR_BstoDmax/DtoDsScale));
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDmax->Write("",TObject::kOverwrite);

  TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDmin=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDmin");
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDmin->Scale((DtoDsScale_BR_BstoDmin/DtoDsScale));
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDmin->Write("",TObject::kOverwrite);


  TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDsmax=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDsmax");
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDsmax->Scale((DtoDsScale_BR_B0toDsmax/DtoDsScale));
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDsmax->Write("",TObject::kOverwrite);

  TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDsmin=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDsmin");
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDsmin->Scale((DtoDsScale_BR_B0toDsmin/DtoDsScale));
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDsmin->Write("",TObject::kOverwrite);

  TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDsmax=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDsmax");
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDsmax->Scale((DtoDsScale_BR_BptoDsmax/DtoDsScale));
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDsmax->Write("",TObject::kOverwrite);

  TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDsmin=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDsmin");
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDsmin->Scale((DtoDsScale_BR_BptoDsmin/DtoDsScale));
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDsmin->Write("",TObject::kOverwrite);

  TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDsmax=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDsmax");
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDsmax->Scale((DtoDsScale_BR_BstoDsmax/DtoDsScale));
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDsmax->Write("",TObject::kOverwrite);

  TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDsmin=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDsmin");
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDsmin->Scale((DtoDsScale_BR_BstoDsmin/DtoDsScale));
  hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDsmin->Write("",TObject::kOverwrite);


	// 


	TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleFrZ=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleFrZ");
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleFrZ->Scale(DtoDsScale_Fr_Z/DtoDsScale);
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleFrZ->Write("",TObject::kOverwrite);

	TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleFrp=(TH1D*)hBtoDsCrossSectionPP_AnaBin_pythiaWeight->Clone("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleFrp");
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleFrp->Scale(DtoDsScale_Fr_p/DtoDsScale);
	hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleFrp->Write("",TObject::kOverwrite);

	TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleFrZ=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleFrZ");
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleFrZ->Scale(DtoDsScale_Fr_Z/DtoDsScale);
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleFrZ->Write("",TObject::kOverwrite);

	TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleFrp=(TH1D*)hBtoDsdNdPtPbPb_AnaBin_pythiaWeight->Clone("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleFrp");
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleFrp->Scale(DtoDsScale_Fr_p/DtoDsScale);
	hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleFrp->Write("",TObject::kOverwrite);






	cout<<"DtoDsScale = " <<DtoDsScale<<endl;


	return 0;


}
