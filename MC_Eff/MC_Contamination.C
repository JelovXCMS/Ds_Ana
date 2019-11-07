#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "vector"

#include <iostream>
#include <fstream>
#include <TString.h>

#include "TBranch.h"
#include "TTree.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"

#include "TRegexp.h"
#include "TNtuple.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TTreeReader.h"
#include "TTimeStamp.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TF1.h"
#include <numeric>
#include "TMath.h"
#include "TLegend.h"
#include "TLatex.h"
#include <TH2D.h>
#include <TH2.h>
#include <TH1D.h>
#include <TH1.h>
#include <TCut.h>
#include <TF1.h>

#include <TStyle.h>

#include <TProfile.h>

#include <TCanvas.h>


#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"


void SetHistStyle(TH1* h,int color = 1, int markerstyle =20, double markersize = 0.9 , int linestyle=1){
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(markerstyle);
  h->SetMarkerSize(markersize);
  h->SetLineStyle(linestyle);
}


void MC_Contamination(){


	InitStyle();

	double tex_upperY=0.93;

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


  TLatex* texColPbPb = new TLatex(0.88,tex_upperY, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
  texColPbPb->SetNDC();
  texColPbPb->SetTextAlign(32);
  texColPbPb->SetTextSize(0.035);
  texColPbPb->SetTextFont(42);

  TLatex* texColpp = new TLatex(0.88,tex_upperY, "pp #sqrt{s} = 5.02 TeV");
  texColpp->SetNDC();
  texColpp->SetTextAlign(32);
  texColpp->SetTextSize(0.035);
  texColpp->SetTextFont(42);



	TFile *f_phi=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/pp_MC_DsMinTree/DsMinTree_pp_MC_Ds_phikkpi_pt4.root");
	TTree *t_phi=(TTree*)f_phi->Get("ntDs");

	TFile *f_f0kkpi=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/pp_MC_DsMinTree/DsMinTree_pp_MC_Ds_f0980kkpi_pt4.root");
	TTree *t_f0kkpi=(TTree*)f_f0kkpi->Get("ntDs");


	TFile *f_f0m1010kkpi=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/pp_MC_DsMinTree/DsMinTree_pp_MC_Ds_f0m1010kkpi_pt4.root");
	TTree *t_f0m1010kkpi=(TTree*)f_f0m1010kkpi->Get("ntDs");

	TFile *f_f0m990w100kkpi=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/pp_MC_DsMinTree/DsMinTree_pp_MC_Ds_f0m990w100kkpi_pt4.root");
	TTree *t_f0m990w100kkpi=(TTree*)f_f0m990w100kkpi->Get("ntDs");

	TFile *f_f0m980w10kkpi=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/pp_MC_DsMinTree/DsMinTree_pp_MC_Ds_f0m980w10kkpi_pt4.root");
	TTree *t_f0m980w10kkpi=(TTree*)f_f0m980w10kkpi->Get("ntDs");
/*
	TFile *f_f0pipipi=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/pp_MC_DsMinTree/DsMinTree_pp_MC_Ds_f0980pipipi_pt4.root");
	TTree *t_f0pipipi=(TTree*)f_f0pipipi->Get("ntDs");
*/

	TFile *f_kstar=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/pp_MC_DsMinTree/DsMinTree_pp_MC_Ds_kstar892kkpi_pt4.root");
	TTree *t_kstar=(TTree*)f_kstar->Get("ntDs");

// Dplus not important

	double nGen_phi=31968;  // count number of ntGen, meet GSignalType
	double nGen_f0kkpi=32264;
	double nGen_f0m1010kkpi=16161;  // not 100% coorect for f0mX
	double nGen_f0m990w100kkpi=16382;
	double nGen_f0m980w10kkpi=16369;
	double nGen_kstar=32116;

	double BR_phi=2.27;
	double BR_f0kkpi=1.15;
	double BR_kstar=2.61; 

	TCut cut_Dpt="Dpt>4 && Dpt<10";
	TCut cut_GenTrue_phi="DsGen==23333";
	TCut cut_GenTrue_f0kkpi="DsGen==24433";
	TCut cut_GenTrue_kstar="DsGen==25544 || DsGen==25555";

	TCut cut_GenTrue_phi_wrong="DsGen==23344";
	TCut cut_GenTrue_f0kkpi_wrong="DsGen==24477";

  float DtktkResmassCutMean=1.01946;
  float DtktkResmassCutWidth=0.009;


	TCut cut_phiMass=Form("DtktkResmass>%f && DtktkResmass<%f ", DtktkResmassCutMean-DtktkResmassCutWidth, DtktkResmassCutMean+DtktkResmassCutWidth);
	TCut cut_phicl2="DtktkRes_chi2cl>0.05";


	TFile *fout=new TFile("pp_MC_contimination.root","RECREATE");

	int nbin_tktkResMass=60;
	double binLo_tktkResMas=0.98;
	double binHi_tktkResMass=1.09;

	seth(nbin_tktkResMass,binLo_tktkResMas,binHi_tktkResMass);
	auto h_phi_tktkResMass=geth("h_phi_tktkResMass");
	auto h_f0kkpi_tktkResMass=geth("h_f0kkpi_tktkResMass");
	auto h_f0m1010kkpi_tktkResMass=geth("h_f0m1010kkpi_tktkResMass");
	auto h_f0m990w100kkpi_tktkResMass=geth("h_f0m990w100kkpi_tktkResMass");
	auto h_f0m980w10kkpi_tktkResMass=geth("h_f0m980w10kkpi_tktkResMass");
	auto h_kstar_tktkResMass=geth("h_kstar_tktkResMass");

	auto h_phi_tktkResMass_wrong=geth("h_phi_tktkResMass_wrong");
	auto h_f0kkpi_tktkResMass_wrong=geth("h_f0kkpi_tktkResMass_wrong");
	auto h_phif0_tktkResMass_wrong=geth("h_fhif0_tktkResMass_wrong");
	

	t_phi->Project("h_phi_tktkResMass","DtktkResmass", cut_Dpt && cut_phicl2 && cut_GenTrue_phi);
	t_f0kkpi->Project("h_f0kkpi_tktkResMass","DtktkResmass", cut_Dpt && cut_phicl2 && cut_GenTrue_f0kkpi);
	t_f0m1010kkpi->Project("h_f0m1010kkpi_tktkResMass","DtktkResmass", cut_Dpt && cut_phicl2 && cut_GenTrue_f0kkpi);
	t_f0m990w100kkpi->Project("h_f0m990w100kkpi_tktkResMass","DtktkResmass", cut_Dpt && cut_phicl2 && cut_GenTrue_f0kkpi);
	t_f0m980w10kkpi->Project("h_f0m980w10kkpi_tktkResMass","DtktkResmass", cut_Dpt && cut_phicl2 && cut_GenTrue_f0kkpi);
	t_kstar->Project("h_kstar_tktkResMass","DtktkResmass", cut_Dpt && cut_phicl2 && cut_GenTrue_kstar);

	t_phi->Project("h_phi_tktkResMass_wrong","DtktkResmass", cut_Dpt && cut_phicl2 && cut_GenTrue_phi_wrong);
	t_f0kkpi->Project("h_f0kkpi_tktkResMass_wrong","DtktkResmass", cut_Dpt && cut_phicl2 && cut_GenTrue_f0kkpi_wrong);


	cout<<"h_phi_tktkResMass_wrong->Integral() = "<<h_phi_tktkResMass_wrong->Integral()<<endl;
	cout<<"h_f0kkpi_tktkResMass_wrong->Integral() = "<<h_f0kkpi_tktkResMass_wrong->Integral()<<endl;

	h_phi_tktkResMass->Scale(BR_phi/nGen_phi);
	h_f0kkpi_tktkResMass->Scale(BR_f0kkpi/nGen_f0kkpi);
	h_f0m1010kkpi_tktkResMass->Scale(BR_f0kkpi/nGen_f0m1010kkpi);
	h_f0m990w100kkpi_tktkResMass->Scale(BR_f0kkpi/nGen_f0m990w100kkpi);
	h_f0m980w10kkpi_tktkResMass->Scale(BR_f0kkpi/nGen_f0m980w10kkpi);
	h_kstar_tktkResMass->Scale(BR_kstar/nGen_kstar);


	h_phi_tktkResMass_wrong->Scale(BR_phi/nGen_phi);
	h_f0kkpi_tktkResMass_wrong->Scale(BR_f0kkpi/nGen_f0kkpi);

	h_phif0_tktkResMass_wrong->Add(h_phi_tktkResMass_wrong);
	h_phif0_tktkResMass_wrong->Add(h_f0kkpi_tktkResMass_wrong);


	// continue here 2018 6.15


	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TCanvas *c_phi=new TCanvas("c_phi","c_phi",900,900);
	c_phi->cd();
	h_phi_tktkResMass->GetXaxis()->SetTitle("m_{KK} (GeV)");
	h_phi_tktkResMass->GetXaxis()->CenterTitle();
	h_phi_tktkResMass->SetLineColor(2);
	h_phi_tktkResMass->SetTitle("");
	h_phi_tktkResMass->SetMarkerColor(2);
	h_phi_tktkResMass->Draw("SAME");
	h_f0kkpi_tktkResMass->SetLineColor(4);
	h_f0kkpi_tktkResMass->SetMarkerColor(4);
	h_f0kkpi_tktkResMass->Draw("SAME");
	h_kstar_tktkResMass->SetLineColor(6);
	h_kstar_tktkResMass->SetMarkerColor(6);
	h_kstar_tktkResMass->Draw("SAME");
	h_phif0_tktkResMass_wrong->SetLineColor(7);
	h_phif0_tktkResMass_wrong->SetMarkerColor(7);
	h_phif0_tktkResMass_wrong->Draw("SAME");

	texCmsSim->Draw("SAME");
	texColpp->Draw("SAME");

	TLegend *le_phi=new TLegend(0.55,0.55,0.85,0.85,NULL,"brNDC");
	le_phi->SetBorderSize(0);	
	le_phi->SetTextFont(42);
	// le_phi->SetTextSize(0.04);
	le_phi->SetFillStyle(0);
	// le_phi->AddEntry((TObject*)0,"4 GeV < p_{T} < 10 GeV ","");
	le_phi->AddEntry(h_phi_tktkResMass,"#phi(->KK) #pi","lp");
	le_phi->AddEntry(h_f0kkpi_tktkResMass,"f0(->KK) #pi","lp");
	le_phi->AddEntry(h_kstar_tktkResMass,"K^{*}(892)K ","lp");
	le_phi->AddEntry(h_phif0_tktkResMass_wrong,"wrong mass assignment","lp");

	le_phi->Draw("SAME");

	SavePlotDirs(c_phi,"mkk",{"MC"});

	/*
	h_phi_tktkResMass_wrong->SetLineColor(6);
	h_phi_tktkResMass_wrong->SetMarkerColor(6);
	h_phi_tktkResMass_wrong->Draw("SAME");
	h_f0kkpi_tktkResMass_wrong->SetLineColor(7);
	h_f0kkpi_tktkResMass_wrong->SetMarkerColor(7);
	h_f0kkpi_tktkResMass_wrong->Draw("SAME");
	*/

	// return;


/*
	Draw({h_phi_tktkResMass,h_f0kkpi_tktkResMass,h_kstar_tktkResMass});
	Draw({h_phi_tktkResMass,h_f0kkpi_tktkResMass,h_f0m1010kkpi_tktkResMass,h_f0m990w100kkpi_tktkResMass,h_f0m980w10kkpi_tktkResMass,h_kstar_tktkResMass});

	Draw({h_f0kkpi_tktkResMass,h_f0m1010kkpi_tktkResMass,h_f0m990w100kkpi_tktkResMass,h_f0m980w10kkpi_tktkResMass});
*/

  int nbin_DsMass=60;
  double binLo_DsMass=1.91;
  double binHi_DsMass=2.11;

	double DsMassMin=1.96828;
	double DsMassCandWidth=0.03;
	double DsMassCandLo=DsMassMin-DsMassCandWidth;
	double DsMassCandHi=DsMassMin+DsMassCandWidth;

  seth(nbin_DsMass,binLo_DsMass,binHi_DsMass);
  auto h_phi_DsMass=geth("h_phi_DsMass");
  auto h_f0kkpi_DsMass=geth("h_f0kkpi_DsMass");
  auto h_f0m1010kkpi_DsMass=geth("h_f0m1010kkpi_DsMass");
  auto h_f0m990w100kkpi_DsMass=geth("h_f0m990w100kkpi_DsMass");
  auto h_f0m980w10kkpi_DsMass=geth("h_f0m980w10kkpi_DsMass");
  auto h_kstar_DsMass=geth("h_kstar_DsMass");

  auto h_phi_DsMass_wrong=geth("h_phi_DsMass_wrong");
  auto h_f0kkpi_DsMass_wrong=geth("h_f0kkpi_DsMass_wrong");
  auto h_phif0_DsMass_wrong=geth("h_phif0_DsMass_wrong");


  t_phi->Project("h_phi_DsMass","Dmass", cut_Dpt && cut_phicl2 && cut_GenTrue_phi && cut_phiMass);
  t_f0kkpi->Project("h_f0kkpi_DsMass","Dmass", cut_Dpt && cut_phicl2 && cut_GenTrue_f0kkpi && cut_phiMass);
  t_f0m1010kkpi->Project("h_f0m1010kkpi_DsMass","Dmass", cut_Dpt && cut_phicl2 && cut_GenTrue_f0kkpi && cut_phiMass);
  t_f0m990w100kkpi->Project("h_f0m990w100kkpi_DsMass","Dmass", cut_Dpt && cut_phicl2 && cut_GenTrue_f0kkpi && cut_phiMass);
  t_f0m980w10kkpi->Project("h_f0m980w10kkpi_DsMass","Dmass", cut_Dpt && cut_phicl2 && cut_GenTrue_f0kkpi && cut_phiMass);
  t_kstar->Project("h_kstar_DsMass","Dmass", cut_Dpt && cut_phicl2 && cut_GenTrue_kstar && cut_phiMass);

  t_phi->Project("h_phi_DsMass_wrong","Dmass", cut_Dpt && cut_phicl2 && cut_GenTrue_phi_wrong && cut_phiMass);
  t_f0kkpi->Project("h_f0kkpi_DsMass_wrong","Dmass", cut_Dpt && cut_phicl2 && cut_GenTrue_f0kkpi_wrong && cut_phiMass);

  h_phi_DsMass->Scale(BR_phi/nGen_phi);
  h_f0kkpi_DsMass->Scale(BR_f0kkpi/nGen_f0kkpi);
  h_f0m1010kkpi_DsMass->Scale(BR_f0kkpi/nGen_f0m1010kkpi);
  h_f0m990w100kkpi_DsMass->Scale(BR_f0kkpi/nGen_f0m990w100kkpi);
  h_f0m980w10kkpi_DsMass->Scale(BR_f0kkpi/nGen_f0m980w10kkpi);
  h_kstar_DsMass->Scale(BR_kstar/nGen_kstar);

  h_phi_DsMass_wrong->Scale(BR_phi/nGen_phi);
  h_f0kkpi_DsMass_wrong->Scale(BR_f0kkpi/nGen_f0kkpi);
	h_phif0_DsMass_wrong->Add(h_phi_DsMass_wrong);
	h_phif0_DsMass_wrong->Add(h_f0kkpi_DsMass_wrong);


	TCanvas *c_Ds=new TCanvas("c_Ds","c_Ds",900,900);
	c_Ds->cd();
	h_phi_DsMass->GetXaxis()->SetTitle("m_{KK#pi} (GeV)");
	h_phi_DsMass->GetXaxis()->CenterTitle();
	h_phi_DsMass->SetLineColor(2);
	h_phi_DsMass->SetTitle("");
	h_phi_DsMass->SetMarkerColor(2);
	h_phi_DsMass->Draw("SAME");

	h_f0kkpi_DsMass->SetLineColor(4);
	h_f0kkpi_DsMass->SetMarkerColor(4);
	h_f0kkpi_DsMass->Draw("SAME");
	h_kstar_DsMass->SetLineColor(6);
	h_kstar_DsMass->SetMarkerColor(6);
	h_kstar_DsMass->Draw("SAME");
	h_phif0_DsMass_wrong->SetLineColor(7);
	h_phif0_DsMass_wrong->SetMarkerColor(7);
	h_phif0_DsMass_wrong->Draw("SAME");

	texCmsSim->Draw("SAME");
	texColpp->Draw("SAME");

	TLegend *le_Ds=new TLegend(0.55,0.55,0.85,0.85,NULL,"brNDC");
	le_Ds->SetBorderSize(0);	
	le_Ds->SetTextFont(42);
	// le_Ds->SetTextSize(0.04);
	le_Ds->SetFillStyle(0);
	// le_Ds->AddEntry((TObject*)0,"4 GeV < p_{T} < 10 GeV ","");
	le_Ds->AddEntry(h_phi_DsMass,"#phi(->KK) #pi","lp");
	le_Ds->AddEntry(h_f0kkpi_DsMass,"f0(->KK) #pi","lp");
	le_Ds->AddEntry(h_kstar_DsMass,"K^{*}(892)K ","lp");
	le_Ds->AddEntry(h_phif0_DsMass_wrong,"wrong mass assignment","lp");

	le_Ds->Draw("SAME");

	SavePlotDirs(c_Ds,"mkkpi",{"MC"});






	double n_phi_Ds=h_phi_DsMass->Integral(h_phi_DsMass->FindBin(DsMassCandLo),h_phi_DsMass->FindBin(DsMassCandHi));
	double n_f0kkpi_Ds=h_f0kkpi_DsMass->Integral(h_f0kkpi_DsMass->FindBin(DsMassCandLo),h_f0kkpi_DsMass->FindBin(DsMassCandHi));
	double n_f0m1010kkpi_Ds=h_f0m1010kkpi_DsMass->Integral(h_f0m1010kkpi_DsMass->FindBin(DsMassCandLo),h_f0m1010kkpi_DsMass->FindBin(DsMassCandHi));
	double n_f0m990w100kkpi_Ds=h_f0m990w100kkpi_DsMass->Integral(h_f0m990w100kkpi_DsMass->FindBin(DsMassCandLo),h_f0m990w100kkpi_DsMass->FindBin(DsMassCandHi));
	double n_f0m980w10kkpi_Ds=h_f0m980w10kkpi_DsMass->Integral(h_f0m980w10kkpi_DsMass->FindBin(DsMassCandLo),h_f0m980w10kkpi_DsMass->FindBin(DsMassCandHi));
	double n_kstar_Ds=h_kstar_DsMass->Integral(h_kstar_DsMass->FindBin(DsMassCandLo),h_kstar_DsMass->FindBin(DsMassCandHi));

	cout<<"n_phi_Ds = "<<n_phi_Ds<<endl;
	cout<<"n_f0kkpi_Ds = "<<n_f0kkpi_Ds<<endl;
	cout<<"n_kstar_Ds = "<<n_kstar_Ds<<endl;

	cout<<"n_phi_Ds / all = "<< n_phi_Ds/ (n_phi_Ds+n_f0kkpi_Ds+n_kstar_Ds)<<endl;		
	cout<<"n_f0kkpi_Ds / all = "<< n_f0kkpi_Ds/ (n_phi_Ds+n_f0kkpi_Ds+n_kstar_Ds)<<endl;		
	cout<<"n_f0m1010kkpi_Ds / all = "<< n_f0m1010kkpi_Ds/ (n_phi_Ds+n_f0m1010kkpi_Ds+n_kstar_Ds)<<endl;		
	cout<<"n_f0m990w100kkpi_Ds / all = "<< n_f0m990w100kkpi_Ds/ (n_phi_Ds+n_f0m990w100kkpi_Ds+n_kstar_Ds)<<endl;		
	cout<<"n_f0m980w10kkpi_Ds / all = "<< n_f0m980w10kkpi_Ds/ (n_phi_Ds+n_f0m980w10kkpi_Ds+n_kstar_Ds)<<endl;		

//  Draw({h_phi_DsMass,h_f0kkpi_DsMass,h_kstar_DsMass});

//  Draw({h_phi_DsMass,h_f0kkpi_DsMass,h_f0m1010kkpi_DsMass,h_f0m990w100kkpi_DsMass,h_f0m980w10kkpi_DsMass,h_kstar_DsMass});
//  Draw({h_f0kkpi_DsMass,h_f0m1010kkpi_DsMass,h_f0m990w100kkpi_DsMass,h_f0m980w10kkpi_DsMass});





	WriteAllHists();

}
