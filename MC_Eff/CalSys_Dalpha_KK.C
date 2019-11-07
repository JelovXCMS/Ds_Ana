#include "../include/uti.h"
#include "../include/DsMinTreeLoad.h"
#include "../include/parameters.h"


#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting_simple.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"


#include <TString.h>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>
#include <TLatex.h>

int CalSys_Dalpha_KK(Int_t isPbPb=0,Int_t DsChannel=0){

  InitStyle();
  double tex_upperY=0.93;

  int nbin_pt=nbin_pt_pp;
  if(isPbPb>0){
    nbin_pt=nbin_pt_PbPb3;
  }


  TString str_PbPb="pp";
  if(isPbPb==3){
    str_PbPb="PbPb3";
  }
  TString str_DsChannel="phikkpi";
  // TString str_DsChannel_plot="D_{S}^{#pm} #rightarrow #phi #pi^{#pm}";
  TString str_DsChannel_plot="D_{S}^{+} #rightarrow #phi#pi^{+} #rightarrow K^{+}K^{-}#pi^{+}";
  if(DsChannel==1){
    str_DsChannel="f0kkpi";
    str_DsChannel_plot="D_{S}^{#pm} #rightarrow f0 #pi^{#pm}";
  }

  // TString fname_NP=Form("./output%s/MC_eff_%s_NonPrompt_%s.root",s_CutSet.Data(), str_PbPb.Data(),str_DsChannel.Data());
  TString fname_MCP_pp=Form("./output%s/MC_eff_pp_Prompt_%s.root",s_CutSet.Data(),str_DsChannel.Data());
  TString fname_MCP_PbPb=Form("./output%s/MC_eff_PbPb3_Prompt_%s.root",s_CutSet.Data(),str_DsChannel.Data());


	gSystem->Exec(Form("mkdir -p ./output%s/sys", s_CutSet.Data()) );	
  TString fname_Sys=Form("./output%s/sys/MC_DalphaKK_Sys_%s.root",s_CutSet.Data(),str_DsChannel.Data());

	TFile *fout=TFile::Open(fname_Sys,"recreate");

	TFile *fpp=TFile::Open(fname_MCP_pp.Data());
	TH1D *h_pp_Norm=(TH1D*)fpp->Get("h_RecoNorm");
	h_pp_Norm->SetLineColor(1);

	TH1D *h_pp_Norm_DalphaWt=(TH1D*)fpp->Get("h_RecoNorm_DalphaWt");
	TH1D *h_pp_Norm_DdlsWt=(TH1D*)fpp->Get("h_RecoNorm_DdlsWt");
	TH1D *h_pp_Norm_KKScl=(TH1D*)fpp->Get("h_RecoNorm_KKScl");

	h_pp_Norm_DalphaWt->SetLineColor(2);
	h_pp_Norm_DdlsWt->SetLineColor(kMagenta+3);
	h_pp_Norm_KKScl->SetLineColor(4);

	TCanvas *c_pp=new TCanvas("c_pp");
	c_pp->cd();
	h_pp_Norm->Draw();	
	h_pp_Norm_KKScl->Draw("same");	
	h_pp_Norm_DalphaWt->Draw("same");	
	h_pp_Norm_DdlsWt->Draw("same");	


	TH1D *h_pp_Dalpha_sys=new TH1D("h_pp_Dalpha_sys","h_pp_Dalpha_sys",nbin_pt_pp,bins_pt_pp);
	TH1D *h_pp_Ddls_sys=new TH1D("h_pp_Ddls_sys","h_pp_Ddls_sys",nbin_pt_pp,bins_pt_pp);
	TH1D *h_pp_KKScl_sys=new TH1D("h_pp_KKScl_sys","h_pp_KKScl_sys",nbin_pt_pp,bins_pt_pp);

	for(int i=1; i<=nbin_pt_pp;i++){
		double sys_Dalpha=fabs(h_pp_Norm->GetBinContent(i)-h_pp_Norm_DalphaWt->GetBinContent(i)) / h_pp_Norm->GetBinContent(i);
		h_pp_Dalpha_sys->SetBinContent(i, sys_Dalpha);
		cout<<"bin i "<<endl;
		cout<<"Dalpha : "<<sys_Dalpha*100<<endl;

		double sys_KKScl=fabs(h_pp_Norm->GetBinContent(i)-h_pp_Norm_KKScl->GetBinContent(i)) / h_pp_Norm->GetBinContent(i);
		h_pp_KKScl_sys->SetBinContent(i, sys_KKScl);
		cout<<"KK : "<<sys_KKScl*100<<endl;

		double sys_Ddls=fabs(h_pp_Norm->GetBinContent(i)-h_pp_Norm_DdlsWt->GetBinContent(i)) / h_pp_Norm->GetBinContent(i);
		h_pp_Ddls_sys->SetBinContent(i, sys_Ddls);
		cout<<"Ddls : "<<sys_Ddls*100<<endl;


	}

	TCanvas *c_pp_sys=new TCanvas("c_pp_sys");
	c_pp_sys->cd();
	h_pp_Dalpha_sys->Draw();
	h_pp_Ddls_sys->Draw("same");
	h_pp_KKScl_sys->Draw("same");

	fout->cd();
	h_pp_Dalpha_sys->Write();
	h_pp_Ddls_sys->Write();
	h_pp_KKScl_sys->Write();

	cout<<"pp done"<<endl;
 	cout<<"\n\n start PbPb"<<endl;


	TFile *fPbPb=TFile::Open(fname_MCP_PbPb.Data());
	TH1D *h_PbPb_Norm=(TH1D*)fPbPb->Get("h_RecoNorm");
	h_PbPb_Norm->SetLineColor(1);

	TH1D *h_PbPb_Norm_DalphaWt=(TH1D*)fPbPb->Get("h_RecoNorm_DalphaWt");
	TH1D *h_PbPb_Norm_DdlsWt=(TH1D*)fPbPb->Get("h_RecoNorm_DdlsWt");
	TH1D *h_PbPb_Norm_KKScl=(TH1D*)fPbPb->Get("h_RecoNorm_KKScl");

	h_PbPb_Norm_DalphaWt->SetLineColor(2);
	h_PbPb_Norm_DdlsWt->SetLineColor(kMagenta+3);
	h_PbPb_Norm_KKScl->SetLineColor(4);

	TCanvas *c_PbPb=new TCanvas("c_PbPb");
	c_PbPb->cd();
	h_PbPb_Norm->Draw();	
	h_PbPb_Norm_KKScl->Draw("same");	
	h_PbPb_Norm_DalphaWt->Draw("same");	
	h_PbPb_Norm_DdlsWt->Draw("same");	


	TH1D *h_PbPb_Dalpha_sys=new TH1D("h_PbPb_Dalpha_sys","h_PbPb_Dalpha_sys",nbin_pt_PbPb3,bins_pt_PbPb3);
	TH1D *h_PbPb_Ddls_sys=new TH1D("h_PbPb_Ddls_sys","h_PbPb_Ddls_sys",nbin_pt_PbPb3,bins_pt_PbPb3);
	TH1D *h_PbPb_KKScl_sys=new TH1D("h_PbPb_KKScl_sys","h_PbPb_KKScl_sys",nbin_pt_PbPb3,bins_pt_PbPb3);

	for(int i=3; i<=nbin_pt_PbPb3;i++){
		double sys_Dalpha=fabs(h_PbPb_Norm->GetBinContent(i)-h_PbPb_Norm_DalphaWt->GetBinContent(i)) / h_PbPb_Norm->GetBinContent(i);
		h_PbPb_Dalpha_sys->SetBinContent(i, sys_Dalpha);
		cout<<"bin i "<<endl;
		cout<<"Dalpha : "<<sys_Dalpha*100<<endl;

		double sys_KKScl=fabs(h_PbPb_Norm->GetBinContent(i)-h_PbPb_Norm_KKScl->GetBinContent(i)) / h_PbPb_Norm->GetBinContent(i);
		h_PbPb_KKScl_sys->SetBinContent(i, sys_KKScl);
		cout<<"KK : "<<sys_KKScl*100<<endl;

  	double sys_Ddls=fabs(h_PbPb_Norm->GetBinContent(i)-h_PbPb_Norm_DdlsWt->GetBinContent(i)) / h_PbPb_Norm->GetBinContent(i);
		h_PbPb_Ddls_sys->SetBinContent(i, sys_Ddls);
		cout<<"Ddls : "<<sys_Ddls*100<<endl;

	}

	TCanvas *c_PbPb_sys=new TCanvas("c_PbPb_sys");
	c_PbPb_sys->cd();
	h_PbPb_Dalpha_sys->Draw();
	h_PbPb_Ddls_sys->Draw("same");
	h_PbPb_KKScl_sys->Draw("same");

	fout->cd();
	h_PbPb_Dalpha_sys->Write();
	h_PbPb_Ddls_sys->Write();
	h_PbPb_KKScl_sys->Write();






	return 1;

}
