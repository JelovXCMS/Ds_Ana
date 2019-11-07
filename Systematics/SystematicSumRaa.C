#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"


#include <TString.h>
#include <TCanvas.h>
#include "TBranch.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLatex.h"
#include <TH2F.h>
#include <TCut.h>
#include <TStyle.h>
#include <TLine.h>


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

#include "TFitter.h"
#include "TFitResult.h"

// using namespace RooFit;
using namespace std;


int SystematicSumRaa(){
	int verbose_sys=1;


  double tex_upperY=0.95;

  TLatex* texCmsPre = new TLatex(0.15,tex_upperY, "#scale[1.25]{CMS} Preliminary");
  texCmsPre->SetNDC();
  texCmsPre->SetTextAlign(12);
  texCmsPre->SetTextSize(0.035);
  texCmsPre->SetTextFont(42);

  TLatex* texCmsSim = new TLatex(0.15,tex_upperY, "#scale[1.25]{CMS} Simulations");
  texCmsSim->SetNDC();
  texCmsSim->SetTextAlign(12);
  texCmsSim->SetTextSize(0.035);
  texCmsSim->SetTextFont(42);


  // TLatex* texColPbPb = new TLatex(0.88,tex_upperY, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  TLatex* texColPbPb = new TLatex(0.95,tex_upperY, "44 #mub^{-1} (5.02 TeV PbPb)");
  texColPbPb->SetNDC();
  texColPbPb->SetTextAlign(32);
  texColPbPb->SetTextSize(0.035);
  texColPbPb->SetTextFont(42);

  // TLatex* texColpp = new TLatex(0.88,tex_upperY, "pp #sqrt{s} = 5.02 TeV");
  TLatex* texColpp = new TLatex(0.95,tex_upperY, "38 nb^{-1} (5.02 TeV pp )");
  texColpp->SetNDC();
  texColpp->SetTextAlign(32);
  texColpp->SetTextSize(0.035);
  texColpp->SetTextFont(42);

  TLatex* texppPbPb = new TLatex(0.95, tex_upperY, "38 pb^{-1} (5.02 TeV pp) + 44 #mub^{-1} (5.02 TeV PbPb)");
  texppPbPb->SetNDC();
  texppPbPb->SetTextAlign(32);
  texppPbPb->SetTextSize(0.035);
  texppPbPb->SetTextFont(42);


	int startbin=0;

	  // if(isPbPb==3){startbin=2;}
  InitStyle();
  initParameter();
  setTDRStyle();
	

  gStyle->SetOptStat(0);
  TString str_PbPb="pp";
  TString str_PbPbtext="pp";


	cout<<"check 1"<<endl;

	TFile *fout=TFile::Open(Form("./output/SysSumRaa.root"),"RECREATE");
	fout->cd();

	TFile *fpp=TFile::Open("./output/SysSum_pp.root","read");
	TH1D *h_Total_SysRel_pp=(TH1D*)fpp->Get("h_Total_SysRel");
	// TH1D *h_partial_forRaa_SysRel_pp=(TH1D*)fpp->Get("h_partial_forRaa_SysRel");
	TH1D *h_DsOverD0_SysRel_pp=(TH1D*)fpp->Get("h_DsOverD0_SysRel");
	//TH1D *h_PromptDs_BtoDs_SysRel_pp=(TH1D*)fpp->Get("h_PromptDs_BtoDs_SysRel");
	TH1D *h_PromptDs_BtoDs_SysRel_pp=(TH1D*)fpp->Get("h_PromptDs_BtoDs_SysRel_ForRaa");
	TH1D *h_PromptDs_BtoDsCancel_SysRel_pp=(TH1D*)fpp->Get("h_PromptDs_BtoDs_SysRel_ForRaaCancel");
	TH1D *h_PromptDs_CutScanAll_SysRel_pp=(TH1D*)fpp->Get("h_PromptDs_CutScanAll_SysRel");
	TH1D *h_DdlsScale_SysRel_pp=(TH1D*)fpp->Get("h_DdlsScale_SysRel");
	TH1D *h_RawRooFitYield_pdfVar_RelErr_pp=(TH1D*)fpp->Get("h_RawRooFitYield_pdfVar_RelErr");
	TH1D *h_PromptDs_MCShape_SysRel_pp=(TH1D*)fpp->Get("h_PromptDs_MCShape_SysRel");
	

	TFile *fPbPb3=TFile::Open("./output/SysSum_PbPb3.root","read");
	TH1D *h_Total_SysRel_PbPb3=(TH1D*)fPbPb3->Get("h_Total_SysRel");
	// TH1D *h_partial_forRaa_SysRel_PbPb3=(TH1D*)fpp->Get("h_partial_forRaa_SysRel");
	TH1D *h_DsOverD0_SysRel_PbPb3=(TH1D*)fPbPb3->Get("h_DsOverD0_SysRel");
	TH1D *h_PromptDs_BtoDs_SysRel_PbPb3=(TH1D*)fPbPb3->Get("h_PromptDs_BtoDs_SysRel_ForRaa");
	TH1D *h_PromptDs_BtoDsCancel_SysRel_PbPb3=(TH1D*)fPbPb3->Get("h_PromptDs_BtoDs_SysRel_ForRaaCancel");
	// TH1D *h_PromptDs_BtoDs_SysRel_PbPb3=(TH1D*)fPbPb3->Get("h_PromptDs_BtoDs_SysRel");
	TH1D *h_DdlsScale_SysRel_PbPb3=(TH1D*)fPbPb3->Get("h_DdlsScale_SysRel");

	TH1D *h_PromptDs_CutScanAll_SysRel_PbPb3=(TH1D*)fPbPb3->Get("h_PromptDs_CutScanAll_SysRel");
	TH1D *h_RawRooFitYield_pdfVar_RelErr_PbPb3=(TH1D*)fPbPb3->Get("h_RawRooFitYield_pdfVar_RelErr");
	TH1D *h_PromptDs_MCShape_SysRel_PbPb3=(TH1D*)fPbPb3->Get("h_PromptDs_MCShape_SysRel");
	

  TFile *fin_pp=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/CrossSection_dNdpt/output/PromptDsCrossSectionPP_FixShape.root","OPEN");
  TFile *fin_PbPb3=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/CrossSection_dNdpt/output/PromptDsdNdptPbPb3_FixShape.root","OPEN");
  TH1D *h_PromptDs_CS_pp=(TH1D*)fin_pp->Get("h_PromptDs_pp");
  TH1D *h_PromptDs_dNdpt_PbPb3=(TH1D*)fin_PbPb3->Get("h_PromptDs_PbPb3");


	fout->cd();
	TH1D *h_Raa_SysRel=new TH1D("h_Raa_SysRel","h_Raa_SysRel",nbin_pt_PbPb3,bins_pt_PbPb3);

	TH1D *h_Raa_SysRel_Trk=new TH1D("h_Raa_SysRel_Trk","h_Raa_SysRel_Trk",nbin_pt_PbPb3,bins_pt_PbPb3);
	TH1D *h_Raa_SysRel_Cut=new TH1D("h_Raa_SysRel_Cut","h_Raa_SysRel_Cut",nbin_pt_PbPb3,bins_pt_PbPb3);
	TH1D *h_Raa_SysRel_DlE=new TH1D("h_Raa_SysRel_DlE","h_Raa_SysRel_DlE",nbin_pt_PbPb3,bins_pt_PbPb3);
	TH1D *h_Raa_SysRel_Pdf=new TH1D("h_Raa_SysRel_Pdf","h_Raa_SysRel_Pdf",nbin_pt_PbPb3,bins_pt_PbPb3);
	TH1D *h_Raa_SysRel_PtS=new TH1D("h_Raa_SysRel_PtS","h_Raa_SysRel_PtS",nbin_pt_PbPb3,bins_pt_PbPb3);
	TH1D *h_Raa_SysRel_BtD=new TH1D("h_Raa_SysRel_BtD","h_Raa_SysRel_BtD",nbin_pt_PbPb3,bins_pt_PbPb3);
	TH1D *h_Raa_SysRel_Phi=new TH1D("h_Raa_SysRel_Phi","h_Raa_SysRel_Phi",nbin_pt_PbPb3,bins_pt_PbPb3);

	h_Total_SysRel_pp->SetTitle("h_Total_SysRel_pp");
	h_DsOverD0_SysRel_pp->SetTitle("h_DsOverD0_SysRel_pp");
	
	h_Total_SysRel_PbPb3->SetName("h_Total_SysRel_PbPb3");
	h_DsOverD0_SysRel_PbPb3->SetName("h_DsOverD0_SysRel_PbPb3");

	int j=0; // pp bin, i : PbPb bin
	double BtoDs_Sys=0;
	double CutScan_Sys=0;
	double DdlErr_Sys=0;
	double pdfVar_Sys=0;
	double MCShape_Sys=0;
	double Raa_Sys=0;
	double Tracking_Sys=0.18;
	for(int i=2; i<nbin_pt_PbPb3; i++){
		j=i+2;

		double dsigma_pp=h_PromptDs_CS_pp->GetBinContent(j+1);
		double dsigma_PbPb=h_PromptDs_dNdpt_PbPb3->GetBinContent(i+1)/TAA0to100;
		double errBtoDs_pp=h_PromptDs_BtoDs_SysRel_pp->GetBinContent(j+1);
		double errBtoDs_PbPb=h_PromptDs_BtoDs_SysRel_PbPb3->GetBinContent(i+1);
		double errBtoDsCancel_pp=h_PromptDs_BtoDsCancel_SysRel_pp->GetBinContent(j+1);
		double errBtoDsCancel_PbPb=h_PromptDs_BtoDsCancel_SysRel_PbPb3->GetBinContent(i+1);

		// BtoDs_Sys=sqrt( pow(errBtoDs_pp,2) + pow( errBtoDs_PbPb,2) - 2/dsigma_PbPb*dsigma_pp*dsigma_pp*dsigma_pp  );
		// BtoDs_Sys=sqrt( pow(errBtoDs_pp,2) + pow( errBtoDs_PbPb,2) - 2.0*errBtoDs_pp*errBtoDs_pp*dsigma_PbPb/dsigma_pp  );
		BtoDs_Sys=sqrt( pow(errBtoDs_pp,2) + pow( errBtoDs_PbPb,2) + abs( pow(errBtoDsCancel_PbPb,2)-pow( errBtoDsCancel_pp ,2 ) ) );
		cout<<setprecision(4)<<std::fixed;
		
//		cout<<"errBtoDs_pp = "<<errBtoDs_pp <<" , errBtoDs_PbPb = "<<errBtoDs_PbPb<<" , sqrt(last) =  "<< sqrt(2.0*errBtoDs_pp*errBtoDs_pp*dsigma_PbPb/dsigma_pp)	<<" ,alter = "<<sqrt(2.0*errBtoDs_pp*errBtoDs_PbPb*dsigma_pp/dsigma_PbPb)<<endl;
		// cout<<"dsigma_pp = "<<dsigma_pp<<"dsigma_PbPb = "<<dsigma_PbPb<<endl;
		// cout<<"errBtoDs_pp^2 = "<<errBtoDs_pp*errBtoDs_pp<<" , errBtoDs_PbPb^2 = "<<errBtoDs_PbPb*errBtoDs_PbPb<<" , last = "<<2.0*errBtoDs_pp*errBtoDs_pp*(dsigma_PbPb/dsigma_pp)<<" , alter = "<<2.0*errBtoDs_PbPb*errBtoDs_PbPb*(dsigma_pp/dsigma_PbPb)<<endl; 


		// BtoDs_Sys=h_PromptDs_BtoDs_SysRel_pp->GetBinContent(j+1) > h_PromptDs_BtoDs_SysRel_PbPb3->GetBinContent(i+1) ? h_PromptDs_BtoDs_SysRel_pp->GetBinContent(j+1) : h_PromptDs_BtoDs_SysRel_PbPb3->GetBinContent(i+1) ;
		// BtoDs_Sys=sqrt(pow(h_PromptDs_BtoDs_SysRel_pp->GetBinContent(j+1),2) + pow(h_PromptDs_BtoDs_SysRel_PbPb3->GetBinContent(i+1),2));
		CutScan_Sys=sqrt(pow(h_PromptDs_CutScanAll_SysRel_pp->GetBinContent(j+1),2) + pow(h_PromptDs_CutScanAll_SysRel_PbPb3->GetBinContent(i+1),2));
		DdlErr_Sys=sqrt(pow(h_DdlsScale_SysRel_pp->GetBinContent(j+1),2) + pow(h_DdlsScale_SysRel_PbPb3->GetBinContent(i+1),2));
		pdfVar_Sys=sqrt(pow(h_RawRooFitYield_pdfVar_RelErr_pp->GetBinContent(j+1),2) + pow(h_RawRooFitYield_pdfVar_RelErr_PbPb3->GetBinContent(i+1),2));
		MCShape_Sys=sqrt(pow(h_PromptDs_MCShape_SysRel_pp->GetBinContent(j+1),2) +pow(h_PromptDs_MCShape_SysRel_PbPb3->GetBinContent(i+1),2) );


		// double PhiRatio_RelErr=sqrt(0.014/0.954*0.014/0.954 + 0.009/0.896*0.009/0.896);	
		double PhiRatio_RelErr=sqrt(PhiRatioErr_pp*PhiRatioErr_pp+PhiRatioErr_PbPb*PhiRatioErr_PbPb);	

		// Raa_Sys=sqrt(pow(BtoDs_Sys,2) +  pow(Tracking_Sys,2) + pow(CutScan_Sys,2) + pow(pdfVar_Sys,2) + pow(MCShape_Sys,2));
		// Raa_Sys=sqrt(pow(BtoDs_Sys,2) +  pow(Tracking_Sys,2) + pow(CutScan_Sys,2) + pow(pdfVar_Sys,2) + pow(MCShape_Sys,2) +  pow(PhiRatio_RelErr,2) + pow(DdlErr_Sys,2));

		// remove PhiRatio_RelErr
		Raa_Sys=sqrt(pow(BtoDs_Sys,2) +  pow(Tracking_Sys,2) + pow(CutScan_Sys,2) + pow(pdfVar_Sys,2) + pow(MCShape_Sys,2) + pow(DdlErr_Sys,2));

		h_Raa_SysRel->SetBinContent(i+1,Raa_Sys);

		h_Raa_SysRel_Trk->SetBinContent(i+1,0.18);
		h_Raa_SysRel_Cut->SetBinContent(i+1,CutScan_Sys);
		h_Raa_SysRel_DlE->SetBinContent(i+1,DdlErr_Sys);
		h_Raa_SysRel_Pdf->SetBinContent(i+1,pdfVar_Sys);
		h_Raa_SysRel_PtS->SetBinContent(i+1,MCShape_Sys);
		h_Raa_SysRel_BtD->SetBinContent(i+1,BtoDs_Sys);
		h_Raa_SysRel_Phi->SetBinContent(i+1,PhiRatio_RelErr);

		cout<<"\n\nibin = "<<i<<" , Raa_Sys = "<<setprecision(1)<<std::fixed<<Raa_Sys*100<<" %"<<endl;
		cout<<"CutScan_Sys = "<<setprecision(1)<<std::fixed<<CutScan_Sys*100<<"% , pdfVar_Sys = "<<setprecision(1)<<std::fixed<<pdfVar_Sys*100<<"% ,MCShape_Sys = "<<setprecision(1)<<std::fixed<<MCShape_Sys*100<<"% , BtoDs_Sys = "<<setprecision(1)<<std::fixed<<BtoDs_Sys*100<<"%"<<" DdlErr_Sys = "<<DdlErr_Sys*100<<"% "<<endl;

		if(verbose_sys==1){
			cout<<"CutScan_Sys_pp = "<<h_PromptDs_CutScanAll_SysRel_pp->GetBinContent(j+1)*100<<", CutScan_Sys_PbPb =  "<<h_PromptDs_CutScanAll_SysRel_PbPb3->GetBinContent(i+1)*100<<endl;	
			cout<<"pdf_Sys_pp = "<<h_RawRooFitYield_pdfVar_RelErr_pp->GetBinContent(j+1)*100<<", pdf_Sys_PbPb =  "<<h_RawRooFitYield_pdfVar_RelErr_PbPb3->GetBinContent(i+1)*100<<endl;	
		}

	}

   cout<<"\n\n --------------------\n Raa syst. "<<endl;
   cout<<"Selection efficiency ";
  for(int i=2; i<nbin_pt_PbPb3; i++){
    cout<<" & "<<setprecision(1)<<setw(4)<<h_Raa_SysRel_Cut->GetBinContent(i+1)*100;
  }
    cout<<" \\\\ \\hline"<<endl;
   cout<<"Signal extraction ";
  for(int i=2; i<nbin_pt_PbPb3; i++){
    cout<<" & "<<setprecision(1)<<setw(4)<<h_Raa_SysRel_Pdf->GetBinContent(i+1)*100;
  }
    cout<<" \\\\ \\hline"<<endl;
   cout<<"MC $p_{T}$ shape";
  for(int i=2; i<nbin_pt_PbPb3; i++){
    cout<<" & "<<setprecision(1)<<setw(4)<<h_Raa_SysRel_PtS->GetBinContent(i+1)*100;
  }
    cout<<" \\\\ \\hline"<<endl;
   cout<<"MC Decay Lenth Tune";
  for(int i=2; i<nbin_pt_PbPb3; i++){
    cout<<" & "<<setprecision(1)<<setw(4)<<h_Raa_SysRel_DlE->GetBinContent(i+1)*100;
  }
    cout<<" \\\\ \\hline"<<endl;
 
   cout<<"Non-prompt \\Ds";
  for(int i=2; i<nbin_pt_PbPb3; i++){
    cout<<" & "<<setprecision(1)<<setw(4)<<h_Raa_SysRel_BtD->GetBinContent(i+1)*100;
  }
    cout<<" \\\\ \\hline"<<endl;

/*
		double PhiRatio_RelErr=sqrt(PhiRatioErr_pp*PhiRatioErr_pp+PhiRatioErr_PbPb*PhiRatioErr_PbPb);	
   cout<<"$\\Ds$ $\\to$ $\\phi$ $\\pi^{\\pm}$  ratio ";
    cout<<" & \\multicolumn{4}{c|}{"<<setprecision(1)<<PhiRatio_RelErr*100<<"}";
    cout<<" \\\\ \\hline"<<endl;
*/

   cout<<"Total bin by bin";
  for(int i=2; i<nbin_pt_PbPb3; i++){
    cout<<" & "<<setprecision(1)<<setw(4)<<h_Raa_SysRel->GetBinContent(i+1)*100;
  }
    cout<<" \\\\ \\hline"<<endl;








	TCanvas *c_SysRaa=new TCanvas("c_SysRaa","c_SysRaa",c_wtopx,c_wtopy,c_W,c_H);
	c_SysRaa->cd();
	SetCanvas(c_SysRaa);
	gPad->SetLogx();

	h_Raa_SysRel->GetXaxis()->SetRangeUser(6,40);
	h_Raa_SysRel->SetMaximum(1.0);
	h_Raa_SysRel->SetMinimum(0);
  h_Raa_SysRel->SetFillColor(kGray);
  h_Raa_SysRel->SetFillStyle(1001);
  h_Raa_SysRel->SetTitle("");
  h_Raa_SysRel->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_Raa_SysRel->GetXaxis()->CenterTitle();
  h_Raa_SysRel->GetYaxis()->SetTitle("uncertainty");
  h_Raa_SysRel->GetYaxis()->CenterTitle();

	h_Raa_SysRel->Draw("SAME E2");

	h_Raa_SysRel_Trk->SetLineColor(kViolet+2);
	h_Raa_SysRel_Trk->Draw("SAME");
	h_Raa_SysRel_Cut->SetLineColor(kOrange+2);
	h_Raa_SysRel_Cut->Draw("SAME");
	h_Raa_SysRel_DlE->SetLineColor(kOrange-6);
	h_Raa_SysRel_DlE->Draw("SAME");
	h_Raa_SysRel_Pdf->SetLineColor(kBlue+1);
	h_Raa_SysRel_Pdf->Draw("SAME");
	h_Raa_SysRel_PtS->SetLineColor(kGreen+4);
	h_Raa_SysRel_PtS->Draw("SAME");
	h_Raa_SysRel_BtD->SetLineColor(kRed+2);
	h_Raa_SysRel_BtD->Draw("SAME");

//	h_Raa_SysRel_Phi->SetLineColor(kYellow+3);
//	h_Raa_SysRel_Phi->Draw("SAME");

	texCmsPre->Draw("SAME");
	texppPbPb->Draw("SAME");

	TLegend *le=new TLegend(0.5,0.59,0.85,0.89);
	le->SetBorderSize(0);

	le->AddEntry(h_Raa_SysRel_Trk,"Tracking Efficiency","l");
	le->AddEntry(h_Raa_SysRel_Cut,"Selection efficiency","l");
	le->AddEntry(h_Raa_SysRel_DlE,"Decay Lenth Errr Scale","l");
	le->AddEntry(h_Raa_SysRel_Pdf,"Signal extraction","l");
	le->AddEntry(h_Raa_SysRel_PtS,"MC p_{T} shape","l");
	le->AddEntry(h_Raa_SysRel_BtD,"Nonprompt D_{S}","l");
	// le->AddEntry(h_Raa_SysRel_Phi,"#phi#pi Channel Ratio","l");
	le->AddEntry(h_Raa_SysRel,"Total for R_{AA}","lf");

	le->Draw("SAME");

	SavePlotDirs(c_SysRaa,"SystRaa",{"Systematics"});

	// h_Lumi_SysRel->Draw("SAME");
	// h_f0_SysRel->Draw("SAME");

	fout->cd();
	h_Total_SysRel_pp->Write("h_Total_SysRel_pp",TObject::kOverwrite);
	h_DsOverD0_SysRel_pp->Write("h_DsOverD0_SysRel_pp",TObject::kOverwrite);
	h_Total_SysRel_PbPb3->Write("h_Total_SysRel_PbPb3",TObject::kOverwrite);
	h_DsOverD0_SysRel_PbPb3->Write("h_DsOverD0_SysRel_PbPb3",TObject::kOverwrite);
	h_Raa_SysRel->Write("",TObject::kOverwrite);

	h_Raa_SysRel_Trk->Write("",TObject::kOverwrite);
	h_Raa_SysRel_Cut->Write("",TObject::kOverwrite);
	h_Raa_SysRel_DlE->Write("",TObject::kOverwrite);
	h_Raa_SysRel_Pdf->Write("",TObject::kOverwrite);
	h_Raa_SysRel_PtS->Write("",TObject::kOverwrite);
	h_Raa_SysRel_BtD->Write("",TObject::kOverwrite);
	h_Raa_SysRel_Phi->Write("",TObject::kOverwrite);

	fout->Close();

	return 0;

}
