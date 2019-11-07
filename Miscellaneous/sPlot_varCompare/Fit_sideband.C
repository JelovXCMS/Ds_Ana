#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/tdrstyle.C"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"


#include <TString.h>
#include <TCanvas.h>
#include "TBranch.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLatex.h"
#include <TH2F.h>
#include <TCut.h>
#include <TStyle.h>

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

#include <RooFit.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooHistPdf.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooGlobalFunc.h>
#include <RooGaussian.h>
#include <RooPolynomial.h>
#include <RooGenericPdf.h>
#include <RooFormulaVar.h>
#include <RooChebychev.h>
#include <RooExponential.h>
#include <RooFitResult.h>
#include <RooWorkspace.h>
#include <RooMCStudy.h>
#include "RooHist.h"
#include "RooConstVar.h"
#include "RooMsgService.h"

#include "RooStats/SPlot.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"

#include "varCompare_para.h"

using namespace RooStats;


using namespace RooFit;
using namespace std;
  double DsDataFitRangeLow =1.91;
  double DsDataFitRangeHigh = 2.11;

  double shiftY=0.05;
  double shiftX=0.32;
  double oneshift=0.070;

  TCanvas *c_binfit_MC[500];
  TCanvas *c_binfit_Data[500];
  TCanvas *c_roofit_MC[500] ;
  TCanvas *c_roofit_Data[500];
  TCanvas *c_roofit_Data_pull[500];
  TCanvas *c_roofit_Data_withpull[500];

  TCanvas *c_roofit_MCstudy[500];
  TCanvas *c_roofit_GoF[500];

  TCanvas *c_roofit_Data_Dca[500][30];

  int count_c_binfit=0;
  int count_c_rootfit=0;

  double FloatWidthErr_Norm=0;
  double FloatWidthVal_Norm=0;
  double FloatWidth_DsMass_Norm=1.965;



int Fit_sideband(int isPbPb=0, TString var="Dalpha"){

	gStyle->SetOptStat(0);
	TString Str_PbPb="pp";
	if(isPbPb){Str_PbPb="PbPb";}


  TString dataName="./output/pp_fitTest.root";
  TString mcName_Prompt="./output/MC_pp_Prompt_phikkpi.root";
  TString mcName_NonPrompt="./output/MC_pp_NonPrompt_phikkpi.root";
  TString mcName_NonPrompt_f0="./output/MC_pp_NonPrompt_f0kkpi.root";

	double Dpt_Low=Dpt_Low_pp;
  double Dpt_High=Dpt_Hight_pp;
	if(isPbPb){
		Dpt_Low=Dpt_Low_PbPb;
		Dpt_High=Dpt_Hight_PbPb;

		dataName="./output/PbPb_fitTest.root";
	 	mcName_Prompt="./output/MC_PbPb_Prompt_phikkpi.root";
	 	mcName_NonPrompt="./output/MC_PbPb_NonPrompt_phikkpi.root";
	 	mcName_NonPrompt_f0="./output/MC_PbPb_NonPrompt_f0kkpi.root";

	}

  TFile *f_data=TFile::Open(dataName.Data());
  TFile *f_mc_Prompt=TFile::Open(mcName_Prompt.Data());
  TFile *f_mc_NonPrompt=TFile::Open(mcName_NonPrompt.Data());
  TFile *f_mc_NonPrompt_f0=TFile::Open(mcName_NonPrompt_f0.Data());


  TTree *t_DsMassDsMC_Prompt=(TTree*)f_mc_Prompt->Get(Form("t_for%s",var.Data()));
  TTree *t_DsMassDsMC_NonPrompt=(TTree*)f_mc_NonPrompt->Get(Form("t_for%s",var.Data()));
  TTree *t_DsMassData=(TTree*)f_data->Get(Form("t_for%s",var.Data()));

	TFile *f_out=TFile::Open(Form("./compare_root/%s_%s.root",Str_PbPb.Data(), var.Data()),"recreate");

	
	// import mc & fit
  RooRealVar Dmass("Dmass","Dmass",DsDataFitRangeLow,DsDataFitRangeHigh);
  // RooRealVar Ddca("Ddca","Ddca",0,0.1); // temp
  // RooRealVar Ddls("Ddls","Ddls",0,200); // temp
  RooRealVar D0DataWeight("D0DataWeight","D0DataWeight",0,1e15);
  RooDataSet RooDS_MC("RooDS_MC","RooDS_MC",RooArgSet(Dmass,D0DataWeight),WeightVar(D0DataWeight),Import(*t_DsMassDsMC_Prompt));
  RooDS_MC.Print("v");

  RooRealVar DsMassMean_MC("DsMassMean_MC","DsMassMean_MC",1.9690,1.965,1.971);
  RooRealVar DsWidth1_MC("DsWidth1_MC","DsWidth1_MC",0.011,0.008,0.02);
  RooRealVar DsWidth2_MC("DsWidth2_MC","DsWidth2_MC",0.0060,0.001,0.01);
  RooRealVar DsGaus1Fr_MC("DsGaus1Fr_MC","DsGaus1Fr_MC",0.65,0.3,1);

  RooGaussian Gaus1_MC("Gaus1_MC","gauss(Dmass,DsMassMean_MC,DsWidth1_MC)",Dmass,DsMassMean_MC,DsWidth1_MC);
  RooGaussian Gaus2_MC("Gaus2_MC","gauss(Dmass,DsMassMean_MC,DsWidth2_MC)",Dmass,DsMassMean_MC,DsWidth2_MC);
  RooAddPdf SigPdf_MC("SigPdf_MC","SigPdf_MC",RooArgList(Gaus1_MC,Gaus2_MC),DsGaus1Fr_MC);
  SigPdf_MC.fitTo(RooDS_MC,NumCPU(10));
  SigPdf_MC.fitTo(RooDS_MC,NumCPU(10));



	// import data
  //RooDataSet RooDSAll("RooDSAll","RooDSAll",RooArgSet(Dmass,Ddls),Import(*t_DsMassData));
  //RooDataSet *RooDS=(RooDataSet*)RooDSAll.reduce(RooArgSet(Dmass));
	
	RooDataSet *RooDS= new RooDataSet("RooDS","RooDS",RooArgSet(Dmass),Import(*t_DsMassData));

  double DsMassMeanV=DsMassMean_MC.getValV();
  double DsWidth1V=DsWidth1_MC.getValV();
  double DsWidth2V=DsWidth2_MC.getValV();
  double DsGaus1FrV=DsGaus1Fr_MC.getValV();

  RooRealVar DsMassMean("DsMassMean","DsMassMean",DsMassMeanV,1.965,1.971);
  RooRealVar DsWidth1("DsWidth1","DsWidth1",DsWidth1V,0.000001,0.1);
  RooRealVar DsWidth2("DsWidth2","DsWidth2",DsWidth2V,0.000001,0.1);
  RooRealVar DsGaus1Fr("DsGaus1Fr","DsGaus1Fr",DsGaus1FrV,0.0,1);
  RooRealVar DsFloatWidth("DsFloatWidth","DsFloatWidth",0,-0.6,0.6);
  // RooRealVar DsFloatWidth("DsFloatWidth","DsFloatWidth",DsFloatWidthV,DsFloatWidthV,DsFloatWidthV);

  DsWidth1.setConstant(kTRUE);
  DsWidth2.setConstant(kTRUE);
  DsGaus1Fr.setConstant(kTRUE);
  RooFormulaVar scale_width1("scale width1","scaled width1","DsWidth1*(1+DsFloatWidth)",RooArgSet(DsWidth1,DsFloatWidth));
  RooFormulaVar scale_width2("scale width2","scaled width2","DsWidth2*(1+DsFloatWidth)",RooArgSet(DsWidth2,DsFloatWidth));
  RooGaussian Gaus1("Gaus1","gauss(Dmass,DsMassMean,scale_width1)",Dmass,DsMassMean,scale_width1);
  RooGaussian Gaus2("Gaus2","gauss(Dmass,DsMassMean,scale_width2)",Dmass,DsMassMean,scale_width2);
  RooAddPdf SigPdf("SigPdf","SigPdf",RooArgList(Gaus1,Gaus2),DsGaus1Fr);
  RooRealVar Cheb1("Cheb1","Cheb1",0,-1,1); // no input is better than with input
  RooRealVar Cheb2("Cheb2","Cheb2",0,-1,1);
  RooRealVar Cheb3("Cheb3","Cheb3",0,-1,1);
  RooChebychev *BkgPdf;
  BkgPdf=new RooChebychev("BkgPdf","BkgPdf",Dmass,RooArgList(Cheb1,Cheb2));
  RooRealVar NumSig("NumSig","Number of Signal",1000,-1e4,4e8);
  RooRealVar NumBkg("NumBkg","Number of Background",10000,0,1e10);

  RooAddPdf *RooDsMixPdf= new RooAddPdf("RooDsMixPdf","RooDsMixPdf",RooArgList(SigPdf,*BkgPdf),RooArgList(NumSig,NumBkg));
  RooDsMixPdf->fitTo(*RooDS,Extended(kTRUE),NumCPU(20));
  RooDsMixPdf->fitTo(*RooDS,Extended(kTRUE),NumCPU(20));
  RooDsMixPdf->fitTo(*RooDS,Extended(kTRUE),NumCPU(20));

	TCanvas *c_test=new TCanvas("c_test","c_test",800,600);
	c_test->cd();

  RooPlot* massframe=new RooPlot("massframe","Dmass",Dmass,DsDataFitRangeLow,DsDataFitRangeHigh,50);
  // RooPlot* massframe=new RooPlot("massframe","Dmass",Dmass,DsDataFitRangeLow,DsDataFitRangeHigh,nbin_DmassDraw);
  RooDS->plotOn(massframe,DataError(RooAbsData::SumW2));
  RooDsMixPdf->plotOn(massframe,LineColor(2));
  massframe->Draw();

/*
	cout<<" \n\nRooDsMixPdf = "<< RooDsMixPdf->getVal()<<endl;

	RooArgSet nset(Dmass);
	cout<<"RooDsMixPdf[Dmass] = "<< RooDsMixPdf->getVal(&nset)<<endl;

	RooAbsReal* igx = RooDsMixPdf->createIntegral(Dmass) ;
	cout << "gx_Int[x] = " << igx->getVal() << endl ;

	Dmass.setRange("signal",1.9490,1.9890);

	RooAbsReal* igx_sig = RooDsMixPdf->createIntegral(Dmass,NormSet(Dmass),Range("signal")) ;
	cout << "gx_Int[x|signal]_Norm[x] = " << igx_sig->getVal() << endl ;	
*/
  cout<<" \n\nSigPdf = "<< SigPdf.getVal()<<endl;

	RooArgSet nset(Dmass);
	cout<<"SigPdf[Dmass] = "<< SigPdf.getVal(&nset)<<endl;

	RooAbsReal* igx = SigPdf.createIntegral(Dmass) ;
	cout << "gx_Int[x] = " << igx->getVal() << endl ;

	Dmass.setRange("signal",1.9490,1.9890);

	RooAbsReal* igx_sig = SigPdf.createIntegral(Dmass,NormSet(Dmass),Range("signal")) ;
	cout << "Sig_Int[x|signal]_Norm[x] = " << igx_sig->getVal() << endl ;	

	RooAbsReal* igx_bkg = BkgPdf->createIntegral(Dmass,NormSet(Dmass),Range("signal")) ;
	cout << "Bkg_Int[x|signal]_Norm[x] = " << igx_bkg->getVal() << endl ;	


	cout<<" \n------------------ \n"<<endl;

	double N_Sig_2sig=NumSig.getValV()* igx_sig->getVal();
	double N_Bkg_2sig=NumBkg.getValV()* igx_bkg->getVal();

	cout<<"N_Sig_2sig = "<<N_Sig_2sig<<endl;
	cout<<"N_Bkg_2sig = "<<N_Bkg_2sig<<endl;

// end fitting
	// return 1;

	// start to plot var by sideband subtraction

	double *bins_var=bins_Ddls_pp;
	int nbin_var= nbin_Ddls_pp;

	if(var=="Dalpha"){
		bins_var=bins_Dalpha_pp;
		nbin_var=nbin_Dalpha_pp;
	}else if(var=="Dchi2cl"){
		bins_var=bins_Dchi2cl_pp;
		nbin_var=nbin_Dchi2cl_pp;
	}

	if(isPbPb){
		bins_var=bins_Ddls_PbPb;
		nbin_var= nbin_Ddls_PbPb;
		if(var=="Dalpha"){
		bins_var=bins_Dalpha_PbPb;
		nbin_var=nbin_Dalpha_PbPb;
		}else if(var=="Dchi2cl"){
		bins_var=bins_Dchi2cl_PbPb;
		nbin_var=nbin_Dchi2cl_PbPb;
		}
	}


	TH1D *h_var_sideband= new TH1D(Form("h_%s_sideband",var.Data()),Form("h_%s_sideband",var.Data()),nbin_var, bins_var); h_var_sideband->Sumw2();
	t_DsMassData->Project(Form("h_%s_sideband",var.Data()),var.Data(),"Dmass<1.93 || Dmass>2.01");

	TH1D *h_var_Data2sig= new TH1D(Form("h_%s_Data2sig",var.Data()),Form("h_%s_Data2sig",var.Data()),nbin_var, bins_var); h_var_Data2sig->Sumw2();
	t_DsMassData->Project(Form("h_%s_Data2sig",var.Data()),var.Data(),"Dmass>1.949 && Dmass<1.989");
	
	h_var_sideband->Scale(N_Bkg_2sig/h_var_sideband->Integral());
	
	h_var_Data2sig->Add(h_var_sideband,-1);
	h_var_Data2sig->Scale(1/h_var_Data2sig->Integral());

	// MC histogram 
	TH1D *h_var_PromptMC =new TH1D(Form("h_%s_PromptMC",var.Data()), Form("h_%s_PromptMC",var.Data()), nbin_var, bins_var); h_var_PromptMC->Sumw2();
	t_DsMassDsMC_Prompt->Project(Form("h_%s_PromptMC",var.Data()),var.Data(),"(Dmass>1.949 && Dmass<1.989)*D0DataWeight");

	h_var_PromptMC->Scale(1/h_var_PromptMC->Integral());

	TH1D *h_var_NonPromptMC =new TH1D(Form("h_%s_NonPromptMC",var.Data()), Form("h_%s_NonPromptMC",var.Data()), nbin_var, bins_var); h_var_NonPromptMC->Sumw2();
	t_DsMassDsMC_NonPrompt->Project(Form("h_%s_NonPromptMC",var.Data()),var.Data(),"(Dmass>1.949 && Dmass<1.989)*D0DataWeight");

  h_var_NonPromptMC->Scale(1/h_var_NonPromptMC->Integral());	



	double fr_withcut_prompt=0.85; // just for get quick result

	// extract fr_prompt from fit and BtoDs estimation
	bool use_fr_fromBtoD=true;
	if(use_fr_fromBtoD){
		double N_yield=NumSig.getValV();
		double LumiNevt=LumiSum;
		if(isPbPb){LumiNevt=NevtPbPb3;}


		TFile *f_NonPromptDs=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BtoDs_Results/output/BtoDs.root");
	  // TH1D *hBtoDsCrossSectionPP_AnaBin=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin");
	  TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight");
	
	  // TH1D *hBtoDsdNdPtPbPb_AnaBin=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin");
	  TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight");
	
	  TH1D *hBtoDs_AnaBin_pythiaWeight=hBtoDsCrossSectionPP_AnaBin_pythiaWeight;  // this is differential cross section
		if(isPbPb){ hBtoDs_AnaBin_pythiaWeight=hBtoDsdNdPtPbPb_AnaBin_pythiaWeight; }

		MutiplyBinWidth(hBtoDs_AnaBin_pythiaWeight); 
		double CS_Integral=hBtoDs_AnaBin_pythiaWeight->Integral(hBtoDs_AnaBin_pythiaWeight->GetXaxis()->FindBin(Dpt_Low+0.00001), hBtoDs_AnaBin_pythiaWeight->GetXaxis()->FindBin(Dpt_High-0.00001) );	

		cout<<"CS_Integral = "<<CS_Integral<<endl;

		TH1D *h_eff_nonprompt_phi=(TH1D*)f_mc_NonPrompt->Get(Form("h_RecoEff_for%s",var.Data()));
		TH1D *h_eff_nonprompt_f0=(TH1D*)f_mc_NonPrompt_f0->Get(Form("h_RecoEff_for%s",var.Data()));

		
	
		double NonPrompt_phi_Eff=h_eff_nonprompt_phi->GetBinContent(1);
		double NonPrompt_f0_Eff=h_eff_nonprompt_f0->GetBinContent(1);

		cout<<"NonPrompt_phi_Eff = "<<NonPrompt_phi_Eff<<endl;
		cout<<"NonPrompt_f0_Eff = "<<NonPrompt_f0_Eff<<endl;
	
		double N_NonPrompt_yield=CS_Integral*(2*LumiNevt)*( BRphi*NonPrompt_phi_Eff+ BRf0*NonPrompt_f0_Eff );

		cout<<"LumiNevt = "<<LumiNevt<<endl;

		fr_withcut_prompt=1-(N_NonPrompt_yield/N_yield);

		cout<<"fr_withcut_prompt = "<<fr_withcut_prompt<<" , N_NonPrompt_yield = "<<N_NonPrompt_yield<<" , N_yield = "<<N_yield<<endl;


	}





	TH1D *h_var_MixMC= new TH1D("h_var_MixMC","h_var_MixMC",nbin_var, bins_var); h_var_MixMC->Sumw2();

	h_var_MixMC->Add(h_var_PromptMC,h_var_NonPromptMC,fr_withcut_prompt,1-fr_withcut_prompt);


	cout<<"start drawing"<<endl;

	TCanvas *c_Ddls_sideband = new TCanvas("c_Ddls_sideband","c_Ddls_sideband",800,600);
	c_Ddls_sideband->Divide(2,2);
	c_Ddls_sideband->cd(1);
	h_var_Data2sig->Draw();
	c_Ddls_sideband->cd(2);
//	h_var_sideband->Draw();
	h_var_PromptMC->Draw();
	c_Ddls_sideband->cd(3);
	h_var_NonPromptMC->Draw();
	c_Ddls_sideband->cd(4);

	h_var_MixMC->SetMaximum(0.3);
	h_var_MixMC->SetMinimum(0.0);
	h_var_MixMC->Draw();
	h_var_Data2sig->SetLineColor(2);
	h_var_Data2sig->Draw("same");

	c_Ddls_sideband->SaveAs(Form("./plots/%s_%s.png",Str_PbPb.Data(),var.Data()));

	InitStyle();
	setTDRStyle();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TCanvas *c_var_compare = new TCanvas("c_var_compare","c_var_compare",c_wtopx,c_wtopy,c_W,c_H);
	SetCanvas(c_var_compare);
	c_var_compare->cd();

  h_var_MixMC->SetMaximum(0.3);
  h_var_MixMC->SetMinimum(0.0);
	h_var_MixMC->SetTitle("");
	h_var_MixMC->GetXaxis()->SetTitle(var.Data());
	h_var_MixMC->GetYaxis()->SetTitle("Normalized Distribution");
  h_var_MixMC->Draw();
  h_var_Data2sig->SetLineColor(2);
	h_var_Data2sig->SetMarkerColor(2);
  h_var_Data2sig->Draw("same");

	TLegend *le = new TLegend(0.6,0.62,0.88,0.85);
	le->SetBorderSize(0);
	le->AddEntry((TObject*)0,Form("%s",Str_PbPb.Data()),"");
	le->AddEntry(h_var_MixMC,"MC","l");
	le->AddEntry(h_var_Data2sig,"data sideband sub.","l");
	le->Draw("same");

	SavePlotDirs(c_var_compare,Form("%s_%s_var_MCDataCompare",Str_PbPb.Data(),var.Data()),{"Miscellaneous","var_MCDataCompare"});



	f_out->cd();
	h_var_MixMC->Write();
	h_var_Data2sig->Write();


	// cumulative for all above, considering the cut is > certain value

	TH1D *h_var_MixMC_Cum= new TH1D("h_var_MixMC_Cum","h_var_MixMC_Cum",nbin_var, bins_var); h_var_MixMC_Cum->Sumw2();
	TH1D *h_var_Data2sig_Cum= new TH1D("h_var_Data2sig_Cum","h_var_Data2sig_Cum",nbin_var, bins_var); h_var_Data2sig_Cum->Sumw2();
	double cum_MixMC=0;
	double cum_Data=0; 

	for(int i=1; i<=nbin_var; i++){
		for(int j=i; j<=nbin_var; j++){
			cum_MixMC+=h_var_MixMC->GetBinContent(j);
			cum_Data+=h_var_Data2sig->GetBinContent(j);			

		} // end for j
		h_var_MixMC_Cum->SetBinContent(i,cum_MixMC);		
		h_var_Data2sig_Cum->SetBinContent(i,cum_Data);		

	}// end for i 


	TH1D *h_var_Ratio_Cum= new TH1D("h_var_Ratio_Cum","h_var_Ratio_Cum",nbin_var, bins_var); h_var_Ratio_Cum->Sumw2();
	h_var_Ratio_Cum->Divide(h_var_Data2sig_Cum,h_var_MixMC_Cum);

	h_var_MixMC_Cum->Write();
	h_var_Data2sig_Cum->Write();
	h_var_Ratio_Cum->Write();



	return 1;



/* not work ,for Splot.
	RooMsgService::instance().setSilentMode(true);


  RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",*RooDS, RooDsMixPdf, RooArgList(NumSig,NumBkg) );

   std::cout << std::endl <<  "Yield of Sig is "
               << NumSig.getVal() << ".  From sWeights it is "
               << sData->GetYieldFromSWeight("NumSig") << std::endl;
   std::cout << "Yield of Bkg is "
               << NumBkg.getVal() << ".  From sWeights it is "
               << sData->GetYieldFromSWeight("NumBkg") << std::endl
               << std::endl;

   std::cout << "import new dataset with sWeights" << std::endl;
//   ws->import(*data, Rename("dataWithSWeights"));
//	RooDataSet *data_wts=(RooDataSet)
  RooDsMixPdf->fitTo(*RooDS,Extended(kTRUE),NumCPU(20));
*/

	// plot the distribution
/*
	TCanvas *c_Ddls= new TCanvas("c_Ddls","c_Ddls",800,600);
	c_Ddls->cd();
   RooDataSet * dataw_s = new RooDataSet(RooDS->GetName(),RooDS->GetTitle(),RooDS,*RooDS->get(),0,"NumSig_sw") ;
	cout<<"check 1"<<endl;

   RooPlot* frame2 = Ddls.frame() ;
	cout<<"check 2"<<endl;
   dataw_s->plotOn(frame2, DataError(RooAbsData::SumW2) ) ;
	cout<<"check 3"<<endl;
   frame2->SetTitle("Ddls distribution for Sig");
   frame2->Draw() ;
*/




	return 0;

}



