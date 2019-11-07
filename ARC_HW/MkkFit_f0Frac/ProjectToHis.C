
#include "MkkFit_para.h"

int ProjectToHis(int isPbPb=0, double DptLow=6, double DptHigh=40, double DVtxPCut=0.25, double DdlsCut=4.5){

	double TrkPtCut=1;

	TString s_ppPbPb="pp";

	// int bins_phi=

  gSystem->Exec("mkdir -p his");
  gSystem->Exec("mkdir -p plots");
  gSystem->Exec("mkdir -p plots/Ori");

	TString s_finName="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_Data/DsMinTree_pp_Data_MBAll.root";

	if(isPbPb){
		s_finName="/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/GJ/Ds_PbPb_Data_HIMB2_GJ.root";
		s_ppPbPb="PbPb";
	}

	TString s_foutName=Form("./his/%s_Dpt%.0fto%.0f.root",s_ppPbPb.Data(),DptLow,DptHigh);

	TString s_fMC_Prompt_Phi=MCFile_Prompt_phikkpi_pp_Merge;
	TString s_fMC_Prompt_f0=MCFile_Prompt_f0980kkpi_pp_Merge;
	
	if(isPbPb){
		s_fMC_Prompt_Phi=MCFile_Prompt_phikkpi_PbPb_Merge;
		s_fMC_Prompt_f0=MCFile_Prompt_f0980kkpi_PbPb_Merge;
	}
	

	// TString s_ppPbPb="pp";

	TFile *fin=TFile::Open(s_finName.Data(),"read");
	TTree *ntDs=(TTree*)fin->Get("ntDs");

	TFile *f_MCP_Phi=TFile::Open(s_fMC_Prompt_Phi.Data());
	TTree *ntDs_MCP_Phi=(TTree*)f_MCP_Phi->Get("ntDs");
	TFile *f_MCP_f0=TFile::Open(s_fMC_Prompt_f0.Data());
	TTree *ntDs_MCP_f0=(TTree*)f_MCP_f0->Get("ntDs");

	TFile *fout=TFile::Open(s_foutName.Data(),"recreate");
	fout->cd();

	TH1D *h_mkkpi_default_MCP_Phi=new TH1D("h_mkkpi_default_MCP_Phi","; m_{KK#pi}",nbin_DsMass*2, DsMassFitLow,DsMassFitHigh);
	h_mkkpi_default_MCP_Phi->Sumw2();

	ntDs_MCP_Phi->Project(Form("h_mkkpi_default_MCP_Phi"),"Dmass",Form("weight*DgenptSampleWeight*RecoD0DataWeight*(Dmass>%f && Dmass<%f && Dpt>%f && Dpt<%f && Dalpha<0.12 && Dchi2cl>%f && Ddls >%f && DtktkResmass>%f && DtktkResmass<%f && Dtrk1Pt>%f && Dtrk2Pt>%f && Dtrk3Pt >%f )",DsMassFitLow,DsMassFitHigh, DptLow, DptHigh, DVtxPCut, DdlsCut, DtktkResmassCutMean-DtktkResmassCutWidth, DtktkResmassCutMean+DtktkResmassCutWidth , TrkPtCut,TrkPtCut,TrkPtCut) );

	h_mkkpi_default_MCP_Phi->Write();

	TCanvas *ctest=new TCanvas("ctest","ctest");
  ctest->cd();
	h_mkkpi_default_MCP_Phi->Draw();
	ctest->SaveAs(Form("plots/Ori/OriHis_MCP_Phi_%s_pt%.0fto%.0f_mkkDefault.png",s_ppPbPb.Data(),DptLow,DptHigh));

	TH1D *h_mkkpi_default_MCP_f0=new TH1D("h_mkkpi_default_MCP_f0","; m_{KK#pi}",nbin_DsMass*2, DsMassFitLow,DsMassFitHigh);
	h_mkkpi_default_MCP_f0->Sumw2();

	ntDs_MCP_f0->Project(Form("h_mkkpi_default_MCP_f0"),"Dmass",Form("weight*DgenptSampleWeight*RecoD0DataWeight*(Dmass>%f && Dmass<%f && Dpt>%f && Dpt<%f && Dalpha<0.12 && Dchi2cl>%f && Ddls >%f && DtktkResmass>%f && DtktkResmass<%f && Dtrk1Pt>%f && Dtrk2Pt>%f && Dtrk3Pt >%f )",DsMassFitLow,DsMassFitHigh, DptLow, DptHigh, DVtxPCut, DdlsCut, DtktkResmassCutMean-DtktkResmassCutWidth, DtktkResmassCutMean+DtktkResmassCutWidth , TrkPtCut,TrkPtCut,TrkPtCut) );

	h_mkkpi_default_MCP_f0->Write();

  ctest->cd();
	h_mkkpi_default_MCP_f0->Draw();
	ctest->SaveAs(Form("plots/Ori/OriHis_MCP_f0_%s_pt%.0fto%.0f_mkkDefault.png",s_ppPbPb.Data(),DptLow,DptHigh));

	// mkk his for template fit

	TH1D *h_mkk_MCP_Phi=new TH1D("h_mkk_MCP_Phi",";m_{KK}",nbin_mkk,bins_mkk);
	h_mkk_MCP_Phi->Sumw2();
	ntDs_MCP_Phi->Project("h_mkk_MCP_Phi","DtktkResmass",Form("weight*DgenptSampleWeight*RecoD0DataWeight*(Dmass>%f && Dmass<%f && Dpt>%f && Dpt<%f && Dalpha<0.12 && Dchi2cl>%f && Ddls >%f && DtktkResmass>%f && DtktkResmass<%f && DsGen==23333 && Dtrk1Pt>%f && Dtrk2Pt>%f && Dtrk3Pt >%f )",DsMassMCMean-DsMass2SigWidth , DsMassMCMean+DsMass2SigWidth, DptLow, DptHigh, DVtxPCut, DdlsCut, bins_mkk[0], bins_mkk[nbin_mkk] , TrkPtCut,TrkPtCut,TrkPtCut) );
	h_mkk_MCP_Phi->Write();
	
	ctest->cd();
	h_mkk_MCP_Phi->Draw();
	ctest->SaveAs(Form("plots/Ori/OriHis_Mkk_MCP_Phi_%s_pt%.0fto%.0f.png",s_ppPbPb.Data(),DptLow,DptHigh));

	TH1D *h_mkk_MCP_f0=new TH1D("h_mkk_MCP_f0",";m_{KK}",nbin_mkk,bins_mkk);
	h_mkk_MCP_f0->Sumw2();
	ntDs_MCP_f0->Project("h_mkk_MCP_f0","DtktkResmass",Form("weight*DgenptSampleWeight*RecoD0DataWeight*(Dmass>%f && Dmass<%f && Dpt>%f && Dpt<%f && Dalpha<0.12 && Dchi2cl>%f && Ddls >%f && DtktkResmass>%f && DtktkResmass<%f && DsGen==24433 && Dtrk1Pt>%f && Dtrk2Pt>%f && Dtrk3Pt >%f )",DsMassMCMean-DsMass2SigWidth , DsMassMCMean+DsMass2SigWidth, DptLow, DptHigh, DVtxPCut, DdlsCut, bins_mkk[0], bins_mkk[nbin_mkk] , TrkPtCut,TrkPtCut,TrkPtCut) );
	h_mkk_MCP_f0->Write();
	
	ctest->cd();
	h_mkk_MCP_f0->Draw();
	ctest->SaveAs(Form("plots/Ori/OriHis_Mkk_MCP_f0_%s_pt%.0fto%.0f.png",s_ppPbPb.Data(),DptLow,DptHigh));


	TH1D *h_mkk_MCP_Phi_fineBin=new TH1D("h_mkk_MCP_Phi_fineBin",";m_{KK}",nbin_mkk*2,bin_mkkLow,bin_mkkHigh);
	h_mkk_MCP_Phi_fineBin->Sumw2();
	ntDs_MCP_Phi->Project("h_mkk_MCP_Phi_fineBin","DtktkResmass",Form("weight*DgenptSampleWeight*RecoD0DataWeight*(Dmass>%f && Dmass<%f && Dpt>%f && Dpt<%f && Dalpha<0.12 && Dchi2cl>%f && Ddls >%f && DtktkResmass>%f && DtktkResmass<%f && DsGen==23333 && Dtrk1Pt>%f && Dtrk2Pt>%f && Dtrk3Pt >%f )",DsMassMCMean-DsMass2SigWidth , DsMassMCMean+DsMass2SigWidth, DptLow, DptHigh, DVtxPCut, DdlsCut, bins_mkk[0], bins_mkk[nbin_mkk] , TrkPtCut,TrkPtCut,TrkPtCut) );
	h_mkk_MCP_Phi_fineBin->Write();
	
	ctest->cd();
	h_mkk_MCP_Phi_fineBin->Draw();
	ctest->SaveAs(Form("plots/Ori/OriHis_Mkk_MCP_Phi_fineBin_%s_pt%.0fto%.0f.png",s_ppPbPb.Data(),DptLow,DptHigh));

	TH1D *h_mkk_MCP_f0_fineBin=new TH1D("h_mkk_MCP_f0_fineBin",";m_{KK}",nbin_mkk*2,bin_mkkLow,bin_mkkHigh);
	h_mkk_MCP_f0_fineBin->Sumw2();
	ntDs_MCP_f0->Project("h_mkk_MCP_f0_fineBin","DtktkResmass",Form("weight*DgenptSampleWeight*RecoD0DataWeight*(Dmass>%f && Dmass<%f && Dpt>%f && Dpt<%f && Dalpha<0.12 && Dchi2cl>%f && Ddls >%f && DtktkResmass>%f && DtktkResmass<%f && DsGen==24433 && Dtrk1Pt>%f && Dtrk2Pt>%f && Dtrk3Pt >%f )",DsMassMCMean-DsMass2SigWidth , DsMassMCMean+DsMass2SigWidth, DptLow, DptHigh, DVtxPCut, DdlsCut, bins_mkk[0], bins_mkk[nbin_mkk] , TrkPtCut,TrkPtCut,TrkPtCut) );
	h_mkk_MCP_f0_fineBin->Write();
	
	ctest->cd();
	h_mkk_MCP_f0_fineBin->Draw();
	ctest->SaveAs(Form("plots/Ori/OriHis_Mkk_MCP_f0_fineBin_%s_pt%.0fto%.0f.png",s_ppPbPb.Data(),DptLow,DptHigh));




	TH1D *h_mkkpi_mkkbins[nbin_mkk];
	TH1D *h_mkkpi_default=new TH1D("h_mkkpi_default","; m_{KK#pi}",nbin_DsMass, DsMassFitLow,DsMassFitHigh);
	h_mkkpi_default->Sumw2();

	ntDs->Project(Form("h_mkkpi_default"),"Dmass",Form("Dmass>%f && Dmass<%f && Dpt>%f && Dpt<%f && Dalpha<0.12 && Dchi2cl>%f && Ddls >%f && DtktkResmass>%f && DtktkResmass<%f && Dtrk1Pt>%f && Dtrk2Pt>%f && Dtrk3Pt >%f  ",DsMassFitLow,DsMassFitHigh, DptLow, DptHigh, DVtxPCut, DdlsCut, DtktkResmassCutMean-DtktkResmassCutWidth, DtktkResmassCutMean+DtktkResmassCutWidth , TrkPtCut,TrkPtCut,TrkPtCut) );

	h_mkkpi_default->Write();

	// TCanvas *ctest=new TCanvas("ctest","ctest");
  ctest->cd();
	h_mkkpi_default->Draw();
	ctest->SaveAs(Form("plots/Ori/OriHis_%s_pt%.0fto%.0f_mkkDefault.png",s_ppPbPb.Data(),DptLow,DptHigh));

	delete ctest;	

	for(int i=0; i<nbin_mkk ; i++){

		cout<<"project bin : "<<i<<endl;
		h_mkkpi_mkkbins[i]=new TH1D(Form("h_mkkpi_mkkbins_%i",i),"; m_{KK#pi}", nbin_DsMass , DsMassFitLow, DsMassFitHigh );
		h_mkkpi_mkkbins[i]->Sumw2();
		ntDs->Project(Form("h_mkkpi_mkkbins_%i",i),"Dmass",Form("Dmass>%f && Dmass<%f && Dpt>%f && Dpt<%f && Dalpha<0.12 && Dchi2cl>%f && Ddls>%f && DtktkResmass>%f && DtktkResmass<%f && Dtrk1Pt>%f && Dtrk2Pt>%f && Dtrk3Pt >%f ",DsMassFitLow,DsMassFitHigh, DptLow, DptHigh, DVtxPCut,DdlsCut , bins_mkk[i], bins_mkk[i+1] , TrkPtCut,TrkPtCut,TrkPtCut) );
		TCanvas *ctest=new TCanvas("ctest","ctest");
		ctest->cd();
		h_mkkpi_mkkbins[i]->Draw();

		// ctest->SaveAs(Form("plots/Ori/OriHis_%i.pdf",i));
		ctest->SaveAs(Form("plots/Ori/OriHis_%s_pt%.0fto%.0f_mkk%.0fto%.0f.png",s_ppPbPb.Data(),DptLow,DptHigh,bins_mkk[i]*1000,bins_mkk[i+1]*1000));

		h_mkkpi_mkkbins[i]->Write();

		delete ctest;

	}

	fout->Close();



	return 0;
}
