#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>

void Alice_Ds(){

	TFile *Alice_result=TFile::Open("./output/Alice_result.root","recreate");

	const int nbin_DsD0_0to10=4;
	double bins_pt_0to10[nbin_DsD0_0to10]={5,7,10,14};
	double bins_pt_0to10_err[nbin_DsD0_0to10]={0.5,0.5,0.5,0.5};
	double DsD0_0to10[nbin_DsD0_0to10]={0.38,0.417,0.329,0.336};
	double DsD0_0to10_Errhi[nbin_DsD0_0to10]={0.468,0.498,0.388,0.402};
	double DsD0_0to10_Errlo[nbin_DsD0_0to10]={0.281,0.321,0.259,0.251};

	for(int i=0; i<nbin_DsD0_0to10;i++){
		DsD0_0to10_Errhi[i]=DsD0_0to10_Errhi[i]-DsD0_0to10[i];
		DsD0_0to10_Errlo[i]=DsD0_0to10[i]-DsD0_0to10_Errlo[i];
	}

	TGraphAsymmErrors *gr_DsD0_0to10=new TGraphAsymmErrors(nbin_DsD0_0to10,bins_pt_0to10,DsD0_0to10,bins_pt_0to10_err,bins_pt_0to10_err,DsD0_0to10_Errlo,DsD0_0to10_Errhi);

	const int nbin_DsD0_pp=4;
	double bins_pt_pp[nbin_DsD0_pp]={3,5,7,9};
	double bins_pt_pp_err[nbin_DsD0_pp]={0.5,0.5,0.5,0.5};
	double DsD0_pp[nbin_DsD0_pp]={0.181,0.211,0.226,0.240};
	double DsD0_pp_Errhi[nbin_DsD0_pp]={0.204,0.240,0.251,0.266};
	double DsD0_pp_Errlo[nbin_DsD0_pp]={0.156,0.185,0.196,0.207};

	for(int i=0; i<nbin_DsD0_pp;i++){
		DsD0_pp_Errhi[i]=DsD0_pp_Errhi[i]-DsD0_pp[i];
		DsD0_pp_Errlo[i]=DsD0_pp[i]-DsD0_pp_Errlo[i];
	}

	TGraphAsymmErrors *gr_DsD0_pp=new TGraphAsymmErrors(nbin_DsD0_pp,bins_pt_pp,DsD0_pp,bins_pt_pp_err,bins_pt_pp_err,DsD0_pp_Errlo,DsD0_pp_Errhi);
	
	TCanvas *c_DsD0=new TCanvas("c_DsD0","c_DsD0",800,600);
	c_DsD0->cd();
	gr_DsD0_0to10->SetFillStyle(3005);
	gr_DsD0_0to10->SetFillColor(6);
	gr_DsD0_0to10->Draw("AP2*");
	gr_DsD0_pp->SetFillStyle(3005);
	gr_DsD0_pp->SetFillColor(4);
	gr_DsD0_pp->Draw("p2same");


	const int nbin_Raa=4;
	double bins_Raa_pt[nbin_Raa]={5,7,10,14};
	double bins_Raa_pt_err[nbin_Raa]={0.5,0.5,0.5,0.5};
	double Raa[nbin_Raa]={0.573,0.390,0.291,0.285};
	double Raa_Errhi[nbin_Raa]={0.724,0.475,0.357,0.376};
	double Raa_Errlo[nbin_Raa]={0.409,0.285,0.220,0.148};

	for(int i=0; i<nbin_Raa; i++){
		Raa_Errhi[i]=Raa_Errhi[i]-Raa[i];
		Raa_Errlo[i]=Raa[i]-Raa_Errlo[i];
	}


	TGraphAsymmErrors *gr_Raa=new TGraphAsymmErrors(nbin_Raa,bins_Raa_pt,Raa,bins_Raa_pt_err,bins_Raa_pt_err,Raa_Errlo,Raa_Errhi);


	TCanvas *c_Raa=new TCanvas("c_Raa","c_Raa",800,600);
	c_Raa->cd();
	gr_Raa->Draw("AP2*");

	Alice_result->cd();
	gr_DsD0_0to10->Write("gr_DsD0_0to10");	
	gr_DsD0_pp->Write("gr_DsD0_pp");
	gr_Raa->Write("gr_Raa");


}
