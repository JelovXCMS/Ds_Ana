#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting_simple.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/CMS_lumi.C"
#include "varCompare_para.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TFitter.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"

int SystSum(int isPbPb=0){

	TString s_ppPbPb="pp";
	TString s_ppPbPb3="pp";
	int nbin_pt=nbin_pt_pp;
	double *bins_pt=bins_pt_pp;
	int startbin=0;
	if(isPbPb){
		s_ppPbPb="PbPb";
		s_ppPbPb3="PbPb3";
		nbin_pt=nbin_pt_PbPb3;
		bins_pt=bins_pt_PbPb3;
		startbin=2;
	}

	TFile *fout=new TFile(Form("./RatioOut_BestScale/CutScanSys_FixShape_%s.root",s_ppPbPb3.Data()),"recreate");

	fout->cd();	

	TString s_var="Ddls";
	cout<<"\n"<<s_var<<endl;
	double Ddls_sys[nbin_pt];
	double Ddls_sys_Scale[nbin_pt];
	double Ddls_sys_noScale[nbin_pt];
	TFile *fins[nbin_pt];
	TH1D *h_Syst[nbin_pt];
	TH1D *h_DdlsScale_Sys[nbin_pt];
	TH1D *h_DdlsScale_SysNoScale[nbin_pt];
	for(int i =startbin; i<nbin_pt; i++){
			fins[i]=TFile::Open(Form("./RatioOut_BestScale/%s/%s_%s_Dpt%.0fto%.0fFixShape.root",s_var.Data(), s_ppPbPb.Data(), s_var.Data(),bins_pt[i],bins_pt[i+1] ));
			h_Syst[i]=(TH1D*)fins[i]->Get("h_Syst");
			cout<<"bin i = "<<i<<" , sys  = "<<h_Syst[i]->GetBinContent(1)*100<<endl;
				Ddls_sys[i]=h_Syst[i]->GetBinContent(1);
			h_DdlsScale_Sys[i]=(TH1D*)fins[i]->Get("h_DdlsScale_Sys");
			h_DdlsScale_SysNoScale[i]=(TH1D*)fins[i]->Get("h_DdlsScale_SysNoScale");
				Ddls_sys_Scale[i]=h_DdlsScale_Sys[i]->GetBinContent(1)	;		
				Ddls_sys_noScale[i]=h_DdlsScale_SysNoScale[i]->GetBinContent(1)	;		
			// cout<<"Ddls_sys_Scale[i] = "<<Ddls_sys_Scale[i]<<" , Ddls_sys_noScale[i] = "<<Ddls_sys_noScale[i]<<endl;	
	}

// Approval_Em method Dalpha
	TFile *f_Dalpha_KK=TFile::Open(Form("../../MC_Eff/output%s/sys/MC_DalphaKK_Sys_phikkpi.root",s_CutSet.Data()));
	TH1D *h_pp_Ddls=(TH1D*)f_Dalpha_KK->Get("h_pp_Ddls_sys");
	TH1D *h_PbPb_Ddls=(TH1D*)f_Dalpha_KK->Get("h_PbPb_Ddls_sys");
	TH1D *h_Ddls=h_pp_Ddls;
	if(isPbPb) h_Ddls=h_PbPb_Ddls;
	for(int i =startbin; i<nbin_pt; i++){
				
			Ddls_sys[i]=h_Ddls->GetBinContent(i+1);
			cout<<"bin i = "<<i<<" , sys  = "<<h_Ddls->GetBinContent(i+1)*100<<endl;
	}




	s_var="Dalpha";
	cout<<"\n"<<s_var<<endl;
	double Dalpha_sys[nbin_pt];
	// TFile *fins[nbin_pt];
	// TH1D *h_Syst[nbin_pt];

/* old method
	for(int i =startbin; i<nbin_pt; i++){
			fins[i]=TFile::Open(Form("./RatioOut_BestScale/%s/%s_%s_Dpt%.0fto%.0fFixShape.root",s_var.Data(), s_ppPbPb.Data(), s_var.Data(),bins_pt[i],bins_pt[i+1] ));
			h_Syst[i]=(TH1D*)fins[i]->Get("h_Syst");
			cout<<"bin i = "<<i<<" , sys  = "<<h_Syst[i]->GetBinContent(1)*100<<endl;
				Dalpha_sys[i]=h_Syst[i]->GetBinContent(1);
	}*/
// Approval_Em method Dalpha
	// TFile *f_Dalpha_KK=TFile::Open(Form("../../MC_Eff/output%s/sys/MC_DalphaKK_Sys_phikkpi.root",s_CutSet.Data()));
	TH1D *h_pp_Dalpha=(TH1D*)f_Dalpha_KK->Get("h_pp_Dalpha_sys");
	TH1D *h_PbPb_Dalpha=(TH1D*)f_Dalpha_KK->Get("h_PbPb_Dalpha_sys");
	TH1D *h_Dalpha=h_pp_Dalpha;
	if(isPbPb) h_Dalpha=h_PbPb_Dalpha;
	for(int i =startbin; i<nbin_pt; i++){
				
			Dalpha_sys[i]=h_Dalpha->GetBinContent(i+1);
			cout<<"bin i = "<<i<<" , sys  = "<<h_Dalpha->GetBinContent(i+1)*100<<endl;
	}

	// return 1;

	s_var="Dchi2cl";
	cout<<"\n"<<s_var<<endl;
	double Dchi2cl_sys[nbin_pt];
	// TFile *fins[nbin_pt];
	// TH1D *h_Syst[nbin_pt];
	for(int i =startbin; i<nbin_pt; i++){
			fins[i]=TFile::Open(Form("./RatioOut_BestScale/%s/%s_%s_Dpt%.0fto%.0fFixShape.root",s_var.Data(), s_ppPbPb.Data(), s_var.Data(),bins_pt[i],bins_pt[i+1] ));
			h_Syst[i]=(TH1D*)fins[i]->Get("h_Syst");
			cout<<"bin i = "<<i<<" , sys  = "<<h_Syst[i]->GetBinContent(1)*100<<endl;
				Dchi2cl_sys[i]=h_Syst[i]->GetBinContent(1);
	}

	s_var="PhiMass";
	cout<<"\n"<<s_var<<endl;
	double PhiMass_sys[nbin_pt];
	// TFile *fins[nbin_pt];
	// TH1D *h_Syst[nbin_pt];
	/* // old Phimass for approval freeze
	for(int i =startbin; i<nbin_pt; i++){
			fins[i]=TFile::Open(Form("./RatioOut_BestScale/%s/%s_%s_Dpt%.0fto%.0fFixShape.root",s_var.Data(), s_ppPbPb.Data(), s_var.Data(),bins_pt[i],bins_pt[i+1] ));
			h_Syst[i]=(TH1D*)fins[i]->Get("h_Syst");
			cout<<"bin i = "<<i<<" , sys  = "<<h_Syst[i]->GetBinContent(1)*100<<endl;
				PhiMass_sys[i]=h_Syst[i]->GetBinContent(1);
	}
*/
// new Approval_Em method Phimass
	TH1D *h_pp_KKScl=(TH1D*)f_Dalpha_KK->Get("h_pp_KKScl_sys");
	TH1D *h_PbPb_KKScl=(TH1D*)f_Dalpha_KK->Get("h_PbPb_KKScl_sys");
	TH1D *h_KKScl=h_pp_KKScl;
	if(isPbPb) h_KKScl=h_PbPb_KKScl;
	for(int i =startbin; i<nbin_pt; i++){
				
			PhiMass_sys[i]=h_KKScl->GetBinContent(i+1);
			cout<<"bin i = "<<i<<" , sys  = "<<h_KKScl->GetBinContent(i+1)*100<<endl;
	}



// cal Sum
	TH1D *h_PromptDs_CutScanAll_SysRel=new TH1D("h_PromptDs_CutScanAll_SysRel","",nbin_pt,bins_pt);
	TH1D *h_DdlErrScale_SysRel=new TH1D("h_DdlsScale_SysRel","",nbin_pt,bins_pt);
	TH1D *h_Ddls_SysRel=new TH1D("h_Ddls_SysRel","",nbin_pt,bins_pt);
	TH1D *h_Dalpha_SysRel=new TH1D("h_Dalpha_SysRel","",nbin_pt,bins_pt);
	TH1D *h_Dchi2cl_SysRel=new TH1D("h_Dchi2cl_SysRel","",nbin_pt,bins_pt);
	TH1D *h_PhiMass_SysRel=new TH1D("h_PhiMass_SysRel","",nbin_pt,bins_pt);

	cout<<"\n Sum"<<endl;

	double SysSum[nbin_pt];
	for(int i =startbin; i<nbin_pt; i++){
		SysSum[i]=sqrt(PhiMass_sys[i]*PhiMass_sys[i] + Dalpha_sys[i]*Dalpha_sys[i] + Dchi2cl_sys[i]*Dchi2cl_sys[i] +  Ddls_sys[i]*Ddls_sys[i]);
		h_PromptDs_CutScanAll_SysRel->SetBinContent(i+1,SysSum[i]);
		h_DdlErrScale_SysRel->SetBinContent(i+1,Ddls_sys_Scale[i]);
		h_Ddls_SysRel->SetBinContent(i+1,Ddls_sys[i]);
		h_Dalpha_SysRel->SetBinContent(i+1,Dalpha_sys[i]);
		h_Dchi2cl_SysRel->SetBinContent(i+1,Dchi2cl_sys[i]);
		h_PhiMass_SysRel->SetBinContent(i+1,PhiMass_sys[i]);
			cout<<"bin i = "<<i<<" , sys  = "<<SysSum[i]*100<<endl;
	}

	cout<<"\n DdlErr Scale"<<endl;
	for(int i =startbin; i<nbin_pt; i++){
			cout<<"bin i = "<<i<<" , sys  = "<<Ddls_sys_Scale[i]*100<<endl;

	}
	cout<<"\n DdlErr 0Scale"<<endl;
	for(int i =startbin; i<nbin_pt; i++){
			cout<<"bin i = "<<i<<" , sys  = "<<Ddls_sys_noScale[i]*100<<endl;

	}

	TCanvas *c_Sys=new TCanvas("c_Sys","",800,800);
	c_Sys->cd();
	gStyle->SetOptStat(0);

		h_PromptDs_CutScanAll_SysRel->SetMaximum(0.4);
		h_PromptDs_CutScanAll_SysRel->SetLineColor(1);
		h_PromptDs_CutScanAll_SysRel->GetXaxis()->SetTitle("p_T (GeV/c)");
		h_PromptDs_CutScanAll_SysRel->GetYaxis()->SetTitle("Relative Systematics");
		h_PromptDs_CutScanAll_SysRel->GetYaxis()->SetTitleOffset(1.3);
		h_PromptDs_CutScanAll_SysRel->Draw();
		h_DdlErrScale_SysRel->Draw("same");
		h_Ddls_SysRel->SetLineColor(kRed+1);
		h_Ddls_SysRel->Draw("same");
		h_Dalpha_SysRel->SetLineColor(kViolet+1);
		h_Dalpha_SysRel->Draw("same");
		h_Dchi2cl_SysRel->SetLineColor(kGreen+2);
		h_Dchi2cl_SysRel->Draw("same");
		h_PhiMass_SysRel->SetLineColor(kYellow+1);
		h_PhiMass_SysRel->Draw("same");

	TLegend *le=new TLegend(0.4,0.45,0.8,0.85);
	le->SetBorderSize(0);
	le->AddEntry((TObject*)0,Form("%s",s_ppPbPb.Data()),"");
	le->AddEntry(h_Ddls_SysRel,"Decay Length significance","l");
	le->AddEntry(h_Dalpha_SysRel,"Pointing Angle","l");
	le->AddEntry(h_Dchi2cl_SysRel,"Vertex Probability","l");
	le->AddEntry(h_PhiMass_SysRel,"Phi Mass Cut","l");
	le->AddEntry(h_DdlErrScale_SysRel,"Decay Length Error Scale","l");
	le->Draw("same");

	gPad->SetLogx();

	fout->cd();
	h_PromptDs_CutScanAll_SysRel->Write();
	h_DdlErrScale_SysRel->Write();

	c_Sys->SaveAs(Form("./RatioOut_BestScale/%s_sys.png",s_ppPbPb.Data()));

		cout<<"\n\n------------\n";
    cout<<"Pointing Angle";
	for(int i =startbin; i<nbin_pt; i++){
			cout<<" &  "<<setprecision(1)<<setw(4)<<std::fixed<<Dalpha_sys[i]*100;
	}
    cout<<" \\\\ \\hline"<<endl;
    cout<<"Vertex Probability";
	for(int i =startbin; i<nbin_pt; i++){
			cout<<" &  "<<setprecision(1)<<setw(4)<<std::fixed<<Dchi2cl_sys[i]*100;
	}
    cout<<" \\\\ \\hline"<<endl;
    cout<<"Decay length significance";
	for(int i =startbin; i<nbin_pt; i++){
			cout<<" &  "<<setprecision(1)<<setw(4)<<std::fixed<<Ddls_sys[i]*100;
	}
    cout<<" \\\\ \\hline"<<endl;
    cout<<"$\\phi$ mass window";
	for(int i =startbin; i<nbin_pt; i++){
			cout<<" &  "<<setprecision(1)<<setw(4)<<std::fixed<<PhiMass_sys[i]*100;
	}
    cout<<" \\\\ \\hline"<<endl;





	return 1;
}
