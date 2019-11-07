#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/DsMinTreeLoad.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"

#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"

#include "TChain.h"

int Compare_Trk(int isPbPb=3, int PNPrompt = 0, int DsChannel=0, double DptLow=0, double DptHigh=10, double TrkPtLow=1, double TrkPtHigh=10){

	// InitStyle();

	TString DataFolder="/mnt/hadoop/store/group/hi/chengchi/Ds_phikkpi_18220/pp_data_resub/MinimumBias13/Data_pp_Ds_phikkpi_MB13/180224_215555/0000/";
	TString MCPFolder="/mnt/hadoop/store/group/hi/chengchi/Dsfinder_officialMC/pp_MC_0817/PrmtDs-PhiPi_pThat-0_DsPT-0_pp_5p02_Pythia8/Dfinder_offMC_pp_Ds_Prmt_PhiPi_pT0_112611/181126_161711/0000/"; 
  TString MCNPFolder="/mnt/hadoop/store/group/hi/chengchi/Dsfinder_officialMC/pp_MC_0817/NonPrDs-PhiPi_pThat-0_DsPT-0_pp_5p02_Pythia8/Dfinder_offMC_pp_Ds_NonPr_PhiPi_pT0_112611/181126_161749/0000/";

	TString s_PbPb="pp";



	if(isPbPb==3){

	DataFolder="/mnt/hadoop/store/group/hi/chengchi/Dsfinder_phikkpi_180305/PbPb_data/MB_GJ/HIMinimumBias3/Data_PbPb_Ds_phikkpi_HiMB3_GJpart3/180316_184608/0006/";
	MCPFolder="/mnt/hadoop/store/group/hi/chengchi/Dsfinder_officialMC/PbPb_MC/PrmtDs-PhiPi_pThat-0_DsPT-0_HydjetCymbMB_5p02_Pythia8/Dfinder_offMC_PbPb_Ds_Prmt_PhiPi_pT0_112611/181126_162002/0000/"; 
  MCNPFolder="/mnt/hadoop/store/group/hi/chengchi/Dsfinder_officialMC/PbPb_MC/NonPrDs-PhiPi_pThat-0_DsPT-0_HydjetCymbMB_5p02_Pythia8/Dfinder_offMC_PbPb_Ds_NonPr_PhiPi_pT0_112611/181126_162049/0000/";
	
	s_PbPb="PbPb3";

	}


	int useHighPtSample=0;
	if(useHighPtSample==1){
		DptLow=10;
		DptHigh=20;
	}


	TChain *t_Data=new TChain("Dfinder/root");
	TChain *t_MCP=new TChain("Dfinder/root");
	TChain *t_MCNP=new TChain("Dfinder/root");

	t_Data->Add(Form("%s*11.root",DataFolder.Data()));
	t_MCP->Add(Form("%s*1.root",MCPFolder.Data()));
	t_MCNP->Add(Form("%s*1.root",MCNPFolder.Data()));

	TFile *fout=new TFile(Form("fout_%s.root",s_PbPb.Data()),"recreate");
	fout->cd();

	int nbin=50;
	double binLow=0;
	double binHigh=0.035;

	TH1D *h_MCP_d0Err= new TH1D("h_MCP_d0Err","h_MCP_d0Err",nbin,binLow,binHigh); 
	h_MCP_d0Err->Sumw2();
	t_MCP->Draw("TrackInfo.d0error>>h_MCP_d0Err",Form("TrackInfo.d0error<%f&&TrackInfo.highPurity==1&&TrackInfo.ptErr/TrackInfo.pt<0.3&&TrackInfo.pt>%f",binHigh,TrkPtLow));
	// h_MCP_d0Err->Scale(1/h_MCP_d0Err->Integral());
	h_MCP_d0Err->SetLineColor(1);

	TCanvas *c3=new TCanvas("c3");
	c3->cd();

	TH1D *h_MCNP_d0Err= new TH1D("h_MCNP_d0Err","h_MCNP_d0Err",nbin,binLow,binHigh); 
	h_MCNP_d0Err->Sumw2();
	t_MCNP->Draw("TrackInfo.d0error>>h_MCNP_d0Err",Form("TrackInfo.d0error<%f&&TrackInfo.highPurity==1&&TrackInfo.ptErr/TrackInfo.pt<0.3&&TrackInfo.pt>%f",binHigh,TrkPtLow));
	// h_MCNP_d0Err->Scale(1/h_MCNP_d0Err->Integral());
	h_MCNP_d0Err->SetLineColor(4);


	TCanvas *c4=new TCanvas("c4");
	c4->cd();

	TH1D *h_Data_d0Err= new TH1D("h_Data_d0Err","h_Data_d0Err",nbin,binLow,binHigh); 
	h_Data_d0Err->Sumw2();
	t_Data->Draw("TrackInfo.d0error>>h_Data_d0Err",Form("TrackInfo.d0error<%f&&TrackInfo.highPurity==1&&TrackInfo.ptErr/TrackInfo.pt<0.3&&TrackInfo.pt>%f",binHigh,TrkPtLow));
	// h_Data_d0Err->Scale(1/h_Data_d0Err->Integral());
	h_Data_d0Err->SetLineColor(2);


	h_MCP_d0Err->Write();
	h_MCNP_d0Err->Write();
	h_Data_d0Err->Write();
	
	h_MCP_d0Err->Scale(1/h_MCP_d0Err->Integral());
	h_MCNP_d0Err->Scale(1/h_MCNP_d0Err->Integral());
	h_Data_d0Err->Scale(1/h_Data_d0Err->Integral());


	TCanvas *c2=new TCanvas("c2","c2",600,600);
	c2->cd();
	h_MCP_d0Err->GetXaxis()->SetTitle("Track d0Error");
	h_MCP_d0Err->SetTitle("");
	h_MCP_d0Err->Draw("same");
	h_MCNP_d0Err->Draw("same");
	h_Data_d0Err->Draw("same");

	TLegend *le=new TLegend(0.6,0.6,0.85,0.85);
	le->SetBorderSize(0);
	le->AddEntry(h_MCP_d0Err,"MCP","l")	;	
	le->AddEntry(h_MCNP_d0Err,"MCNP","l")	;	
	le->AddEntry(h_Data_d0Err,"Data","l")	;	
	le->Draw("same");

	c2->SaveAs(Form("TrackD0Error_%s.png",s_PbPb.Data()));


	return 0;

}

