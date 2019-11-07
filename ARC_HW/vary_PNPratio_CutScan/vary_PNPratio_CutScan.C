#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"

// #include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"
// #include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
// #include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"


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
#include <TGraph.h>
#include <TGraphErrors.h>

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
#include "TFitResultPtr.h"

double textposx=0.2;
double textposy=0.80;

double shiftY=0;
double shiftX=0.28;
double oneshift=0.075;


int vary_PNPratio_CutScan(int isPbPb=0, Int_t startbin=0, Int_t start_var=0){
  if(isPbPb==3){startbin=2;}
   // if(isPbPb==3){startbin=4;}
  if(isPbPb==0){startbin=0;}
  // StartBin=startbin;

	gSystem->Exec("mkdir -p plots");
	gSystem->Exec("mkdir -p plots/Summary");

	TLatex *tltx=new TLatex();
  TFitResultPtr fitR;
  double x0[1]={0};
  double x0err[1];

	TLine *tli=new TLine(0,0,0,0);
	tli->SetLineColor(2);
	tli->SetLineStyle(6);

  gStyle->SetOptStat(0);

  TString str_PbPb="pp";
	TString sstr_ppPbPb="PP";
  TString str_PbPbtext="pp";
  int nbin_pt=nbin_pt_pp;
  double *bins_pt=bins_pt_pp;

  double LumiNevt=LumiSum;
  TString str_eff_Prompt_phi=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_pp_Prompt_phikkpi.root",s_CutSet.Data());
  TString str_eff_NonPrompt_phi=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_pp_NonPrompt_phikkpi.root",s_CutSet.Data());
  TString str_eff_Prompt_f0=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_pp_Prompt_f0kkpi.root",s_CutSet.Data());
  TString str_eff_NonPrompt_f0=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_pp_NonPrompt_f0kkpi.root",s_CutSet.Data());

	TString fcs_name="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/CrossSection_dNdpt/output/PromptDsCrossSectionPP_FixShape.root";

  if(isPbPb==3){
    LumiNevt=NevtPbPb3;
    nbin_pt=nbin_pt_PbPb3;
    bins_pt=bins_pt_PbPb3;
    str_eff_Prompt_phi=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_PbPb3_Prompt_phikkpi.root",s_CutSet.Data());
    str_eff_NonPrompt_phi=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_PbPb3_NonPrompt_phikkpi.root",s_CutSet.Data());
    str_eff_Prompt_f0=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_PbPb3_Prompt_f0kkpi.root",s_CutSet.Data());
    str_eff_NonPrompt_f0=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_Eff/output%s/MC_eff_PbPb3_NonPrompt_f0kkpi.root",s_CutSet.Data());
     str_PbPb="PbPb3";
     sstr_ppPbPb="PbPb3";
		str_PbPbtext="PbPb";
			fcs_name="/home/peng43/work/Project/Ds_PbPb//Ds_PbPbAna/CrossSection_dNdpt/output/PromptDsdNdptPbPb3_FixShape.root";
  }

  // reading fit raw yield


	TFile *f_cs=TFile::Open(fcs_name.Data());
	TH1D  *h_CSfr_PromptDs=(TH1D*)f_cs->Get(Form("h_CSfr_PromptDs_%s",str_PbPb.Data()));

	// h_CSfr_PromptDs->Draw();


   cout<<__LINE__<<endl;
  TString s_f_raw=Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/SignalFit/output%s/RawFitYield_FixShapeTrkPtScan_%s.root",s_CutSet.Data(),str_PbPb.Data());

  TFile *f_raw=TFile::Open(s_f_raw.Data(),"read");
  TH1F *h_RawFitYield=(TH1F*)f_raw->Get("h_RawBinFitYield");
  cout<<__LINE__<<endl;
 // for Efficiency histogram
  TFile *f_Prompt_phikkpi= TFile::Open(str_eff_Prompt_phi.Data());
  TFile *f_NonPrompt_phikkpi= TFile::Open(str_eff_NonPrompt_phi.Data());
  TFile *f_Prompt_f0kkpi= TFile::Open(str_eff_Prompt_f0.Data());
  TFile *f_NonPrompt_f0kkpi= TFile::Open(str_eff_NonPrompt_f0.Data());

  TString str_effH="h_RecoNormEff";

  TH1D *h_Eff_Prompt_phikkpi=(TH1D*)f_Prompt_phikkpi->Get(str_effH.Data());
  TH1D *h_Eff_NonPrompt_phikkpi=(TH1D*)f_NonPrompt_phikkpi->Get(str_effH.Data());
  TH1D *h_Eff_Prompt_f0kkpi=(TH1D*)f_Prompt_f0kkpi->Get(str_effH.Data());
  TH1D *h_Eff_NonPrompt_f0kkpi=(TH1D*)f_NonPrompt_f0kkpi->Get(str_effH.Data());

  TH1D *h_Eff_Prompt_AllBR=(TH1D*)h_Eff_Prompt_phikkpi->Clone("h_Eff_Prompt_AllBR");
  h_Eff_Prompt_AllBR->Sumw2();
  h_Eff_Prompt_AllBR->Add(h_Eff_Prompt_phikkpi,h_Eff_Prompt_f0kkpi,BRphi,BRf0);

  TH1D *h_Eff_NonPrompt_AllBR=(TH1D*)h_Eff_NonPrompt_phikkpi->Clone("h_Eff_NonPrompt_AllBR");
  h_Eff_NonPrompt_AllBR->Sumw2();
  h_Eff_NonPrompt_AllBR->Add(h_Eff_NonPrompt_phikkpi,h_Eff_NonPrompt_f0kkpi,BRphi,BRf0);

  // read efficiency for cut scan for Ddls
  TH1D *h_Eff_Prompt_phikkpi_DdlsMinScan[nbin_DdlsMinScan];
  TH1D *h_Eff_Prompt_f0kkpi_DdlsMinScan[nbin_DdlsMinScan];
  TH1D *h_Eff_NonPrompt_phikkpi_DdlsMinScan[nbin_DdlsMinScan];
   TH1D *h_Eff_NonPrompt_f0kkpi_DdlsMinScan[nbin_DdlsMinScan];
  TH1D *h_Eff_Prompt_AllBR_DdlsMinScan[nbin_DdlsMinScan];
  TH1D *h_Eff_NonPrompt_AllBR_DdlsMinScan[nbin_DdlsMinScan];
  str_effH="h_RecoNormEff_Ddls";
  TH1F *h_RawFitYield_DdlsMinScan[nbin_DdlsMinScan];
 
  for(int i=0; i<nbin_DdlsMinScan; i++){
    h_Eff_Prompt_phikkpi_DdlsMinScan[i]=(TH1D*)f_Prompt_phikkpi->Get(Form("%s_%i",str_effH.Data(),i));
    h_Eff_Prompt_f0kkpi_DdlsMinScan[i]=(TH1D*)f_Prompt_f0kkpi->Get(Form("%s_%i",str_effH.Data(),i));
    h_Eff_NonPrompt_phikkpi_DdlsMinScan[i]=(TH1D*)f_NonPrompt_phikkpi->Get(Form("%s_%i",str_effH.Data(),i));
    h_Eff_NonPrompt_f0kkpi_DdlsMinScan[i]=(TH1D*)f_NonPrompt_f0kkpi->Get(Form("%s_%i",str_effH.Data(),i));
    h_Eff_Prompt_AllBR_DdlsMinScan[i]=(TH1D*)h_Eff_Prompt_phikkpi_DdlsMinScan[i]->Clone(Form("h_Eff_Prompt_AllBR_DdlsMinScan_%i",i));
    h_Eff_Prompt_AllBR_DdlsMinScan[i]->Sumw2();
    h_Eff_Prompt_AllBR_DdlsMinScan[i]->Add( h_Eff_Prompt_phikkpi_DdlsMinScan[i],h_Eff_Prompt_f0kkpi_DdlsMinScan[i],BRphi,BRf0);
    h_Eff_NonPrompt_AllBR_DdlsMinScan[i]=(TH1D*)h_Eff_NonPrompt_phikkpi_DdlsMinScan[i]->Clone(Form("h_Eff_NonPrompt_AllBR_DdlsMinScan_%i",i));
    h_Eff_NonPrompt_AllBR_DdlsMinScan[i]->Sumw2();
    h_Eff_NonPrompt_AllBR_DdlsMinScan[i]->Add( h_Eff_NonPrompt_phikkpi_DdlsMinScan[i],h_Eff_NonPrompt_f0kkpi_DdlsMinScan[i],BRphi,BRf0);
    h_RawFitYield_DdlsMinScan[i]=(TH1F*)f_raw->Get(Form("h_RawRooFitYield_DdlsMinScan_%i",i));

  }
  // end efficiency for cut scan for Ddls

  
	// reading original cut scan 

	TFile *f_cutScan=TFile::Open(Form("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/Systematics/output/CutScanSys_FixShape_%s.root",str_PbPb.Data()));

	TH1D *h_PromptDs_DdlsMinScan_ratio_pt[nbin_pt];

	TGraphErrors *gr_PromptDs_DdlsMinScan_Syst[nbin_pt];
	TGraphErrors *gr_PromptDs_DdlsMinScan_Syst_Def[nbin_pt];
	TGraphErrors *gr_PromptDs_DdlsMinScan_Slope[nbin_pt];
	// TGraph *gr_Best_PromptFr=new TGraph();

	TH1D *h_Best_PromptFr=new TH1D("h_Best_PromptFr",";p_{T};Prompt Fraction",nbin_pt,bins_pt); h_Best_PromptFr->Sumw2();

	double fr_Low=0.7;
	double fr_High=1.0;
	double fr_step=0.01;
	double fr_best=0;
	double sys_best=100;
	const int nstep=(fr_High-fr_Low)/fr_step + 1;

	TH1D *h_PromptDs_DdlsMinScan_vary_PNP[nstep][nbin_pt];

	int nbin_var=nbin_DdlsMinScan;

	TF1 *f1_pol1=new TF1("f1_pol1","[0]+x*[1]");
	f1_pol1->SetLineColor(2);
	f1_pol1->SetLineStyle(2);
	f1_pol1->SetParameter(0,1);
	f1_pol1->SetParameter(1,0);
	f1_pol1->SetRange(0,8);

	double his_high=2.0;
	double his_low=0.4;

	TCanvas *c_def=new TCanvas("c_def","c_def",800,800);

	TCanvas *c_var=new TCanvas("c_var","c_var",800,800);

	TCanvas *c_sysSum=new TCanvas("c_sysSum","c_sysSum",800,800);

	// nbin_pt=1; // for test

	for(int i=startbin; i<nbin_pt; i++){

		gr_PromptDs_DdlsMinScan_Syst[i]=new TGraphErrors();
		gr_PromptDs_DdlsMinScan_Syst_Def[i]=new TGraphErrors();
		gr_PromptDs_DdlsMinScan_Syst_Def[i]->SetLineColor(2);
		gr_PromptDs_DdlsMinScan_Syst_Def[i]->SetMarkerStyle(23);
		gr_PromptDs_DdlsMinScan_Syst_Def[i]->SetMarkerColor(2);
		gr_PromptDs_DdlsMinScan_Slope[i]=new TGraphErrors();

		double fr_def=h_CSfr_PromptDs->GetBinContent(i+1);
		cout<<"fr_def = "<<fr_def<<endl;
		h_PromptDs_DdlsMinScan_ratio_pt[i]=(TH1D*)f_cutScan->Get(Form("h_PromptDs_DdlsMinScan_ratio_pt%.0fto%.0f",bins_pt[i],bins_pt[i+1]));
		h_PromptDs_DdlsMinScan_ratio_pt[i]->GetXaxis()->SetTitle("Decay Length Significance");
		h_PromptDs_DdlsMinScan_ratio_pt[i]->GetXaxis()->SetTitleOffset(0.9);
		h_PromptDs_DdlsMinScan_ratio_pt[i]->GetYaxis()->SetTitleOffset(0.9);
		h_PromptDs_DdlsMinScan_ratio_pt[i]->GetYaxis()->SetTitle("#sigma_{var} / #sigma_{default}");

		h_PromptDs_DdlsMinScan_ratio_pt[i]->Fit("f1_pol1","EMIS0");
		h_PromptDs_DdlsMinScan_ratio_pt[i]->Fit("f1_pol1","EMIS0");
		fitR=h_PromptDs_DdlsMinScan_ratio_pt[i]->Fit("f1_pol1","EMIS0");
    fitR->GetConfidenceIntervals(1,1,1,x0,x0err,0.683,false);
	
		c_def->cd();
		h_PromptDs_DdlsMinScan_ratio_pt[i]->SetTitle("");
		h_PromptDs_DdlsMinScan_ratio_pt[i]->SetMaximum(his_high);
		h_PromptDs_DdlsMinScan_ratio_pt[i]->SetMinimum(his_low);
		h_PromptDs_DdlsMinScan_ratio_pt[i]->Draw();
		f1_pol1->Draw("same");

		double Systematic= abs(f1_pol1->GetParameter(0)-1);
		double SysError = x0err[0];

		double Slope= f1_pol1->GetParameter(1);
	  double SlopeErr= f1_pol1->GetParError(1);

		// tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Systematic %.1f%% #pm %.1f%%", Systematic*100, SysError*100));
		// tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Systematic %.1f%% #pm %.1f%%", Systematic*100, SysError*100));

		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s  %.0f < p_{T} < %.0f ", str_PbPb.Data(), bins_pt[i],bins_pt[i+1] )); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("(Default) Prompt fr :%.2f ", fr_def)); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Systematic : %.1f%% #pm %.1f%%", Systematic*100, SysError*100)); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Slope : %.2e #pm %.2e", Slope, SlopeErr));
		shiftY=0;

		c_def->SaveAs(Form("./plots/%s_Ddls_pt%.0fto%.0f_defaultFR.png",str_PbPb.Data(),bins_pt[i],bins_pt[i+1]));

		gr_PromptDs_DdlsMinScan_Syst_Def[i]->SetPoint(1, fr_def, Systematic );
		gr_PromptDs_DdlsMinScan_Syst_Def[i]->SetPointError(1, 0, SysError );


		double PromptEff_def=h_Eff_Prompt_AllBR->GetBinContent(i+1);
		double NonPromptEff_def=h_Eff_NonPrompt_AllBR->GetBinContent(i+1);
 
		// return 1;
		fr_best=0;
		sys_best=100;
		
		int start_var=0;
		if(isPbPb==3 && i<=2){
			start_var=1;
		}

		for(int k =0; k<nstep; k++){	
			h_PromptDs_DdlsMinScan_vary_PNP[k][i]=(TH1D*)h_PromptDs_DdlsMinScan_ratio_pt[i]->Clone(Form("h_PromptDs_DdlsMinScan_vary_PNP%i_pt%.0fto%.0f",80+k,bins_pt[i],bins_pt[i+1]));
			h_PromptDs_DdlsMinScan_vary_PNP[k][i]->GetXaxis()->SetTitle("Decay Length Significance");
			double fr_new=fr_Low+k*fr_step;

		for(int j=start_var; j<nbin_var; j++){
			double PromptEff=h_Eff_Prompt_AllBR_DdlsMinScan[j]->GetBinContent(i+1);
			double NonPromptEff=h_Eff_NonPrompt_AllBR_DdlsMinScan[j]->GetBinContent(i+1);
			double var_ratio=h_PromptDs_DdlsMinScan_ratio_pt[i]->GetBinContent(j+2);		
	
			cout<<"PromptEff = "<<PromptEff<<" NonPromptEff = "<<NonPromptEff<<endl;
			cout<<"var_ratio = "<<var_ratio<<endl;

			// for(int k=0; k<nstep; k++){

			double newRatio= var_ratio * ( PromptEff_def + (1-fr_new)/fr_new*NonPromptEff_def ) / ( PromptEff+(1-fr_new)/fr_new*NonPromptEff ) / ( PromptEff_def + (1-fr_def)/fr_def*NonPromptEff_def ) * ( PromptEff+(1-fr_def)/fr_def*NonPromptEff  );
			double newRatioErr= newRatio* h_PromptDs_DdlsMinScan_ratio_pt[i]->GetBinError(j+2)  /var_ratio;


			cout<<"new_ratio = "<<newRatio<<endl;

			h_PromptDs_DdlsMinScan_vary_PNP[k][i]->SetBinContent(j+2, newRatio);
			h_PromptDs_DdlsMinScan_vary_PNP[k][i]->SetBinError(j+2, newRatioErr);

			// }

			} // end for j varbin
		h_PromptDs_DdlsMinScan_vary_PNP[k][i]->Fit("f1_pol1","EMIS0");
		h_PromptDs_DdlsMinScan_vary_PNP[k][i]->Fit("f1_pol1","EMIS0");
		fitR=h_PromptDs_DdlsMinScan_vary_PNP[k][i]->Fit("f1_pol1","EMIS0");
    fitR->GetConfidenceIntervals(1,1,1,x0,x0err,0.683,false);
	
		c_var->cd();
		h_PromptDs_DdlsMinScan_vary_PNP[k][i]->SetMaximum(his_high);
		h_PromptDs_DdlsMinScan_vary_PNP[k][i]->SetMinimum(his_low);
		h_PromptDs_DdlsMinScan_vary_PNP[k][i]->GetXaxis()->SetTitle("Decay Length Significance");
		h_PromptDs_DdlsMinScan_vary_PNP[k][i]->SetTitle("");
		// h_PromptDs_DdlsMinScan_vary_PNP[k][i]->SetXTitle("11");
		h_PromptDs_DdlsMinScan_vary_PNP[k][i]->Draw();
		f1_pol1->Draw("same");

		Systematic= abs(f1_pol1->GetParameter(0)-1);
		SysError = x0err[0];

		Slope= f1_pol1->GetParameter(1);
	  SlopeErr= f1_pol1->GetParError(1);

		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s  %.0f < p_{T} < %.0f ", str_PbPb.Data(), bins_pt[i],bins_pt[i+1] )); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Prompt fr :%.2f ", fr_new)); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Systematic : %.1f%% #pm %.1f%%", Systematic*100, SysError*100)); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Slope : %.2e #pm %.2e", Slope, SlopeErr));

		shiftY=0;

		c_var->SaveAs(Form("./plots/%s_Ddls_pt%.0fto%.0f_PNPFR%.0f.png",str_PbPb.Data(),bins_pt[i],bins_pt[i+1],fr_new*100));

		gr_PromptDs_DdlsMinScan_Syst[i]->SetPoint(k, fr_new, Systematic );
		gr_PromptDs_DdlsMinScan_Syst[i]->SetPointError(k, 0, SysError );

		gr_PromptDs_DdlsMinScan_Slope[i]->SetPoint(k, fr_new, Slope );
		gr_PromptDs_DdlsMinScan_Slope[i]->SetPointError(k, 0, SlopeErr );

			if(Systematic<sys_best){
				fr_best=fr_new;	
				sys_best=Systematic;
			}

		} // end for k frPNP variation

		// gr_Best_PromptFr->SetPoint(i,(bins_pt[i+1]+bins_pt[i]) )
		h_Best_PromptFr->SetBinContent(i+1, fr_best);
		h_Best_PromptFr->SetBinError(i+1, 0.01);

		double gr_Sys_high=0.8;
		double gr_Sys_low=-0.25;

		c_sysSum->cd();
		gr_PromptDs_DdlsMinScan_Syst[i]->SetMarkerStyle(21);
		gr_PromptDs_DdlsMinScan_Syst[i]->SetMinimum(gr_Sys_low);
		gr_PromptDs_DdlsMinScan_Syst[i]->SetMaximum(gr_Sys_high);
		gr_PromptDs_DdlsMinScan_Syst[i]->GetXaxis()->SetTitle("Prompt fraction");
		gr_PromptDs_DdlsMinScan_Syst[i]->GetYaxis()->SetTitle("Systematics");
		gr_PromptDs_DdlsMinScan_Syst[i]->GetYaxis()->SetTitleOffset(1.2);
		gr_PromptDs_DdlsMinScan_Syst[i]->Draw("ap");

		gr_PromptDs_DdlsMinScan_Syst_Def[i]->Draw("p");
/*
		tli->SetX1(fr_def);
		tli->SetX2(fr_def);
		tli->SetY1(gr_Sys_high);
		tli->SetY2(gr_Sys_low);
		tli->Draw("same");
*/
		tltx->DrawLatexNDC(textposx+shiftX-0.3,textposy+shiftY,Form("%s  %.0f < p_{T} < %.0f , #sigma_{Decay Length}", str_PbPb.Data(), bins_pt[i],bins_pt[i+1] )); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX-0.3,textposy+shiftY,Form("default Prompt fraction : %.3f ", fr_def )); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX-0.3,textposy+shiftY,Form("Best Prompt fraction : %.2f ", fr_best )); shiftY-=oneshift;
		shiftY=0;

		c_sysSum->SaveAs(Form("./plots/Summary/%s_Ddls_pt%.0fto%.0f_varPNPFR_syst.png",str_PbPb.Data(),bins_pt[i],bins_pt[i+1]));


	} // end for i pt bin	


   	// h_PromptDs_DdlsMinScan_ratio_pt[0]->Draw();                                                    


	// h_PromptDs_DdlsMinScan_vary_PNP[0][1]->Draw();
	// h_PromptDs_DdlsMinScan_vary_PNP[0][1]->SetLineColor(2);
	// h_PromptDs_DdlsMinScan_vary_PNP[20][1]->Draw("same");
/*
	TCanvas *ctest=new TCanvas("ctest","ctest");
	ctest->cd();

	// gr_PromptDs_DdlsMinScan_Syst[0]->Draw("ap");
	// gr_PromptDs_DdlsMinScan_Slope[0]->SetMaximum(1.2);
	// gr_PromptDs_DdlsMinScan_Slope[0]->SetMinimum(-0.2);
	// gr_PromptDs_DdlsMinScan_Slope[0]->SetMarkerStyle(22);
	
	// gr_PromptDs_DdlsMinScan_Slope[0]->Draw("ap");
	gr_PromptDs_DdlsMinScan_Syst[0]->Draw("ap");
*/
	cout<<"nstep = "<<nstep<<endl;

	// test

	TCanvas *c_PromptFr=new TCanvas("c_PromptFr","c_PromptFr",800,800);
	c_PromptFr->cd();
	h_CSfr_PromptDs->SetMarkerStyle(21);
	h_CSfr_PromptDs->SetMaximum(1.1);
	h_CSfr_PromptDs->SetMinimum(0);
	h_CSfr_PromptDs->GetXaxis()->SetTitle("D_{S} p_{T} (GeV/c)");
	h_CSfr_PromptDs->GetYaxis()->SetTitle("Prompt Fraction");
	h_CSfr_PromptDs->GetYaxis()->SetTitleOffset(1.1);
	h_CSfr_PromptDs->SetTitle("");
	h_CSfr_PromptDs->Draw();
	h_Best_PromptFr->SetLineColor(2);
	h_Best_PromptFr->SetMarkerColor(2);
	h_Best_PromptFr->SetMarkerStyle(21);
	h_Best_PromptFr->Draw("same");

	TLegend *le_PromptFr=new TLegend(0.25,0.2,0.6,0.5) ;
	le_PromptFr->SetBorderSize(0);
	le_PromptFr->AddEntry((TObject*)0,Form("%s D_{S}",str_PbPbtext.Data()), "");
	le_PromptFr->AddEntry(h_CSfr_PromptDs,"Default","pl");
	le_PromptFr->AddEntry(h_Best_PromptFr,"Best FR_{prompt} for Ddls","pl");
	le_PromptFr->Draw("same");

	c_PromptFr->SaveAs(Form("./plots/Summary/%s_Ddls_PromptFr.png",str_PbPb.Data()) ) ;



	gSystem->Exec("mkdir -p output");
	TFile *fout=new TFile(Form("./output/%s.root",str_PbPb.Data()),"RECREATE");
	fout->cd();
	h_Best_PromptFr->Write();
	h_CSfr_PromptDs->Write();
	for(int i =startbin; i<nbin_pt; i++){
			h_PromptDs_DdlsMinScan_ratio_pt[i]->Write();
			gr_PromptDs_DdlsMinScan_Syst[i]->Write();
			gr_PromptDs_DdlsMinScan_Slope[i]->Write();
		for(int k=0; k<nstep; k++){

			h_PromptDs_DdlsMinScan_vary_PNP[k][i]->Write();
		}

	}

	return 0;

}
