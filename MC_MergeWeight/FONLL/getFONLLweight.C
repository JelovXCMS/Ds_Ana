/*
// remove some part of dotestmtscling_pythia , keep in older version (approval freeze)
// implment D0 data weight from HIN-16-016


// compare pythia from merged sample and FONLL to give FONLL based weight




// 2019. 10. 4 input changed to TGraph , and use liner interpolations for weight

*/



#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/DsMinTreeLoad.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"

#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"



#include <TString.h>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>
#include <TCut.h>

using namespace std;


TCanvas *can[500];
int can_count=0;

int getFONLLweight(int isPbPb=3, int PNPrompt =0, int DsChannel=0){

	bool dotestmtscling_pythia=0; // just for test

	double D0Mass=1.86483;

	double Fr_B0=0.389; // value from https://arxiv.org/pdf/1306.3663.pdf LHCb 7TeV pp
  double Fr_Bp=0.381;
  double Fr_Bs=0.105;

  double BR_B0toDs=0.103;
  double BR_BptoDs=0.09;
  double BR_BstoDs=0.93;

	double Fr_BtoDs=Fr_B0*BR_B0toDs+Fr_Bp*BR_BptoDs+Fr_Bs*BR_BstoDs;
	cout<<"Fr_BtoDs = "<<Fr_BtoDs<<endl;

	InitStyle();
	initParameter();

  TString str_isPbPb="pp";
	TString str_PNPrompt="Prompt";
	TString str_DsChannel="phi";

  Int_t DsGenTrue=23333;

  int GSignalTypeTrue=1;

  if(DsChannel==1){GSignalTypeTrue=2;}; // f0980
  if(DsChannel==2){GSignalTypeTrue=3;}; // kstar892
  TCut cutGenTrue=Form("GSignalType==%i && GcollisionId==0 && TMath::Abs(Gy)<1",GSignalTypeTrue);
  TCut cutGenPNPrompt="GBAncestorpt<=0";
  if(PNPrompt==1){cutGenPNPrompt="GBAncestorpt>0";}

  TString GenWeight="weight*GptSampleWeight";
  TString RecoWeight="weight*DgenptSampleWeight";
  if(isPbPb){
    GenWeight="weight*GptSampleWeight*PbPbVzWeight*Ncoll";
    // GenWeight="weight*GptSampleWeight";
    RecoWeight="weight*DgenptSampleWeight*PbPbVzWeight*Ncoll";
  }

	double binFirst=0;
	double binLast=50;
	double binwidth=0.25;

  if(isPbPb>0)  {
    str_isPbPb=Form("PbPb");
  }
  if(PNPrompt==1){
    str_PNPrompt="NonPrompt";
		binLast=150;
  }
  if(DsChannel==1){
    str_DsChannel="f0980";
    DsGenTrue=24433;
  }
  if(DsChannel==2){
    str_DsChannel="kstar892";
    DsGenTrue=25544;
  }

	double current_bin=binFirst;
	vector<Double_t> v_bins;
	v_bins.push_back(0);
	while(current_bin<binLast){
		if(current_bin==10 || current_bin == 20  || current_bin == 40 || current_bin == 60) { binwidth*=2 ; }
		current_bin+=binwidth;
		v_bins.push_back(current_bin);	
		cout<<"current_bin = "<<current_bin<<endl;
	}


	int nbin=v_bins.size();
	double bins[nbin];
	double bins_mtD0toDs[nbin-3];
	cout<<"nbin = "<<nbin<<endl;
	for(int i=0; i<nbin; i++){
			bins[i]=v_bins.at(i);
			cout<<"i = "<<i<<" , D0 ptbins = "<< bins[i]<<"\t";
			if(i>=3){
			bins_mtD0toDs[i-3]=sqrt(bins[i]*bins[i]+D0Mass*D0Mass-DsMass*DsMass);
			cout<<" bins_mtD0toDs = "<<bins_mtD0toDs[i-3]<<endl;
			}
	}


  TCut cutRecoTrue=Form("DsGen==%i && DgencollisionId==0 ",DsGenTrue);
  TCut cutRecoPNPrompt="DgenBAncestorpt<=0";
	if(PNPrompt==1){cutRecoPNPrompt="DgenBAncestorpt>0" ;}


	cout<<"running "<<str_isPbPb<<" , "<<str_PNPrompt<<" , "<<str_DsChannel<<endl;

	// read FONLL inputs

	TFile *f_FONLL=TFile::Open("FONLL.root","READ");
	TH1D *h_FONLL_D0=(TH1D*)f_FONLL->Get("h_FONLL_PromptD0");
	TH1D *h_FONLL_Ds_byD0MTs=(TH1D*)f_FONLL->Get("h_FONLL_PromptDs_byD0MTs");
	TH1D *h_FONLL_B=(TH1D*)f_FONLL->Get("h_FONLL_B");
	TH1D *h_FONLL_B_toDs=(TH1D*)f_FONLL->Get("h_FONLL_B_toDs");

	TH1D *h_DsRaa_PHSD=(TH1D*)f_FONLL->Get("h_DsRaa_PHSD");
	TH1D *h_BRaa_TAMU=(TH1D*)f_FONLL->Get("h_BRaa_TAMU");

	TGraph *gr_DsRaa=(TGraph*)f_FONLL->Get("gr_DsRaa_PHSD");

	TGraph *gr_FONLL_D0=(TGraph*)f_FONLL->Get("gr_FONLL_PromptD0");
	TGraph *gr_FONLL_Ds_byD0MTs=(TGraph*)f_FONLL->Get("gr_FONLL_PromptDs_byD0MTs");
	TGraph *gr_FONLL_B=(TGraph*)f_FONLL->Get("gr_FONLL_B");
	TGraph *gr_FONLL_B_toDs=(TGraph*)f_FONLL->Get("gr_FONLL_B_toDs");

	TGraph *gr_DsRaa_PHSD=(TGraph*)f_FONLL->Get("gr_DsRaa_PHSD");
	TGraph *gr_BRaa_TAMU=(TGraph*)f_FONLL->Get("gr_BRaa_TAMU");
	TGraphErrors *grErr_BRaa_TAMU=(TGraphErrors*)f_FONLL->Get("grErr_BRaa_TAMU");

	TGraph *gr_D0Data_pp_weight=(TGraph*)f_FONLL->Get("gr_D0Data_pp_weight");
	TGraph *gr_D0Data_PbPb_weight=(TGraph*)f_FONLL->Get("gr_D0Data_PbPb_weight");

	// from HIN-16-016 result
  TF1 *f1_DataWeight=new TF1("f1_DataWeight","[0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*log(x)*log(x)*log(x)+ [4]*log(x)*log(x)*log(x)*log(x)");
if(PNPrompt==0 && isPbPb==0 ){ f1_DataWeight->SetParameters(3.33929,-1.61349,0.438465,-0.0411468,0); }  
else if(PNPrompt==0 && isPbPb>0 ){f1_DataWeight->SetParameters(4.21138,-4.0426,1.34612,-0.134656,0); }	 
else if(PNPrompt==1 && isPbPb==0 ){ f1_DataWeight->SetParameters(0.55369,1.22995,-0.49766,0.05655,0); }
else if(PNPrompt==1 && isPbPb>0 ){ f1_DataWeight->SetParameters(1.43138,0.37192,-0.71770,0.23269,-0.02181); }

	cout<<"read FONLL root done "<<endl;

	TH1D *h_Raa=NULL;
	TGraph *gr_Raa=nullptr;
	TH1D *h_FONLL;
	TGraph *gr_FONLL=nullptr;

	if(PNPrompt==0){
		h_FONLL=h_FONLL_Ds_byD0MTs;
		// h_FONLL=h_FONLL_D0;
		h_Raa=h_DsRaa_PHSD;

		gr_FONLL=gr_FONLL_Ds_byD0MTs;
		gr_Raa=gr_DsRaa_PHSD;
	}else if(PNPrompt==1){
		h_FONLL=h_FONLL_B_toDs;
//		h_FONLL->Scale(Fr_BtoDs);
		h_Raa=h_BRaa_TAMU;
		
		gr_FONLL=gr_FONLL_B_toDs;
		gr_Raa=gr_BRaa_TAMU;

	}else{
		cout<<"PNP unknow : "<<PNPrompt<<endl;
		return 1;
	}

	if(isPbPb){
		h_FONLL->Scale(TAA0to100);
		// h_FONLL_Ds_byD0MTs->Scale(TAA0to100);

		for(int i=0;i<gr_FONLL->GetN();i++){
			gr_FONLL->GetY()[i]*=TAA0to100;
		}	
		for(int i=0;i<gr_FONLL_Ds_byD0MTs->GetN();i++){
			gr_FONLL_Ds_byD0MTs->GetY()[i]*=TAA0to100;
		}	

	}

	can[can_count]=new TCanvas(Form("can%i",can_count),Form("can%i",can_count),800,800);
	can[can_count]->cd();
	h_FONLL->Draw();
	// h_FONLL_Ds_byD0MTs->SetLineColor(2);
	// h_FONLL_Ds_byD0MTs->Draw("same");
	

	can_count++;
// return 1;	
	
////--  mt scaling for Prompt ---////

  int nbin_mt=180;
	double bin_Ds[nbin_mt+1];
	double bin_D0[nbin_mt+1];
	//   double DsMass=1.96828;
	TH1D *h_PromptD0_pp_Gpt_mtscaling=NULL;


	// TFile *f_pythia=TFile::Open(Form("../root_output/DsMinTree_MC_GenSampleMerge_%s_%s_%s.root",str_isPbPb.Data(),str_PNPrompt.Data(),str_DsChannel.Data()),"READ");
	TFile *f_pythia=TFile::Open(Form("/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/DsMinTree_MC_GenSampleMerge_%s_%s_%s.root",str_isPbPb.Data(),str_PNPrompt.Data(),str_DsChannel.Data()),"READ");


	TTree *t_gen=(TTree*)f_pythia->Get("ntGen");
	// TTree *t_reco=(TTree*)f_pythia->Get("ntDs");

cout<<__LINE__<<endl;

	t_gen->Print();
	// t_reco->Print();

cout<<__LINE__<<endl;

	TFile *f_out=TFile::Open(Form("./output/FONLLweightOverPythia_%s_%s_%s.root", str_isPbPb.Data(),str_PNPrompt.Data(),str_DsChannel.Data()),"RECREATE");
	TH1D *h_Genpt_pythia=new TH1D("h_Genpt_pythia","h_Genpt_pythia",nbin-1,bins); h_Genpt_pythia->Sumw2();
	TH1D *h_Genpt_pythia_mtbin=new TH1D("h_Genpt_pythia_mtbin","h_Genpt_pythia_mtbin",nbin-4,bins_mtD0toDs); h_Genpt_pythia_mtbin->Sumw2();


cout<<__LINE__<<endl;

	TString var_pt="Gpt";
	if(PNPrompt==1){
		var_pt="GBAncestorpt";
	}
	
	// TCut cutGenAll=(TCut)(cutGenTrue&&cutGenPNPrompt)*GenWeight.Data();
	// t_gen->Draw(Form("%s>>h_Genpt_pythia",var_pt.Data()), (TCut)(cutGenTrue && cutGenPNPrompt)*GenWeight.Data() );
	t_gen->Project("h_Genpt_pythia",var_pt.Data(), (TCut)(cutGenTrue && cutGenPNPrompt)*GenWeight.Data() );
	divideBinWidth(h_Genpt_pythia);

	// return 1;

	t_gen->Project("h_Genpt_pythia_mtbin",var_pt.Data(), (TCut)(cutGenTrue && cutGenPNPrompt)*GenWeight.Data() );
	divideBinWidth(h_Genpt_pythia_mtbin);


	can[can_count]=new TCanvas(Form("can%i",can_count),Form("can%i",can_count),800,800);
	can[can_count]->cd();
	h_Genpt_pythia->Draw();
	can_count++;

	cout<<"project pythia default Gpt done"<<endl;


	//-- test mt scaling with pythia D0,Ds  continue
	if (dotestmtscling_pythia){

	TH1D *h_Genpt_pythia_binmt=new TH1D("h_Genpt_pythia_binmt","h_Genpt_pythia_binmt",nbin_mt,bin_Ds); h_Genpt_pythia_binmt->Sumw2();
	t_gen->Project("h_Genpt_pythia_binmt",var_pt.Data(), (TCut)(cutGenTrue && cutGenPNPrompt)*GenWeight.Data() );
	divideBinWidth(h_Genpt_pythia_binmt);

	can[can_count]=new TCanvas(Form("can%i",can_count),Form("can%i",can_count),800,800);
	can[can_count]->cd();
	TH1D *h_D0Ds_Pythia_mtratio=(TH1D*)h_PromptD0_pp_Gpt_mtscaling->Clone("h_D0Ds_Pythia_mtratio");
	h_D0Ds_Pythia_mtratio->Divide(h_Genpt_pythia_binmt);
	h_D0Ds_Pythia_mtratio->Draw();
	can_count++;

	}
	//-- end test mt scaling 



	// Calculate Weight
	double FONLLWeight=0;
	double FONLLWeight_mt=0;
	double FONLLRaaWeight=0;
	double FONLLRaaWeight_mt=0;
	TH1D *h_Genpt_FONLLWeight= new TH1D("h_Genpt_FONLLWeight","h_Genpt_FONLLWeight",nbin-1,bins); h_Genpt_FONLLWeight->Sumw2();
	TH1D *h_Genpt_FONLLWeight_mt= new TH1D("h_Genpt_FONLLWeight_mt","h_Genpt_FONLLWeight_mt",nbin-4,bins_mtD0toDs); h_Genpt_FONLLWeight_mt->Sumw2();
	TH1D *h_Genpt_FONLLRaaWeight= new TH1D("h_Genpt_FONLLRaaWeight","h_Genpt_FONLLRaaWeight",nbin-1,bins); h_Genpt_FONLLRaaWeight->Sumw2();
	TH1D *h_Genpt_FONLLRaaWeight_mt= new TH1D("h_Genpt_FONLLRaaWeight_mt","h_Genpt_FONLLRaaWeight_mt",nbin-1,bins); h_Genpt_FONLLRaaWeight_mt->Sumw2();

	for(int ibin=0; ibin<nbin-1 ; ibin++){


		// FONLLWeight=h_FONLL->GetBinContent( h_FONLL->FindBin(h_Genpt_pythia->GetBinCenter(ibin+1)) ) / h_Genpt_pythia->GetBinContent(ibin+1);
		FONLLWeight=h_FONLL->GetBinContent( h_FONLL->FindBin(h_Genpt_FONLLWeight->GetBinCenter(ibin+1)  ) ) / h_Genpt_pythia->GetBinContent( h_Genpt_pythia->FindBin( h_Genpt_FONLLWeight->GetBinCenter(ibin+1)  ));
		h_Genpt_FONLLWeight->SetBinContent( ibin+1 , FONLLWeight );	
		h_Genpt_FONLLWeight->SetBinError(ibin+1 , 0 );	


		FONLLRaaWeight=h_FONLL->GetBinContent( h_FONLL->FindBin(h_Genpt_FONLLRaaWeight->GetBinCenter(ibin+1)) ) / h_Genpt_pythia->GetBinContent( h_Genpt_pythia->FindBin( h_Genpt_FONLLWeight->GetBinCenter(ibin+1) )) *gr_Raa->Eval(h_Genpt_FONLLRaaWeight->GetBinCenter(ibin+1),0,"" )  ;

/*
		FONLLRaaWeight=h_FONLL->GetBinContent( h_FONLL->FindBin(h_Genpt_pythia->GetBinCenter(ibin+1)) ) / h_Genpt_pythia->GetBinContent(ibin+1) *h_Raa->GetBinContent(h_Raa->FindBin(h_Genpt_pythia->GetBinCenter(ibin+1)))  ;

		if(PNPrompt==0){
		FONLLRaaWeight=h_FONLL->GetBinContent( h_FONLL->FindBin(h_Genpt_pythia->GetBinCenter(ibin+1)) ) / h_Genpt_pythia->GetBinContent(ibin+1) *gr_DsRaa->Eval(h_Genpt_pythia->GetBinCenter(ibin+1),0,"" )  ;
		}
*/

		h_Genpt_FONLLRaaWeight->SetBinContent( ibin+1 , FONLLRaaWeight );	
		h_Genpt_FONLLRaaWeight->SetBinError(ibin+1 , 0 );	

		cout<<"weight : "<<FONLLWeight<<", h_Genpt_FONLLWeight center = "<<h_Genpt_FONLLWeight->GetBinCenter(ibin+1)<< ", h_Genpt_pythia bin center = "<<h_Genpt_pythia->GetBinCenter( h_Genpt_pythia->FindBin(h_Genpt_FONLLWeight->GetBinCenter(ibin+1)) )<<" , h_FONLL bincenter = "<< h_FONLL->GetBinCenter( h_FONLL->FindBin(h_Genpt_FONLLWeight->GetBinCenter(ibin+1)  ) )<<endl;
		// cout<<"weight : "<<FONLLWeight<<", h_FONLL bin center = "<<h_FONLL->GetBinCenter( h_FONLL->FindBin(h_Genpt_pythia->GetBinCenter(ibin+1)) )<< ", h_Genpt_pythia bin center = "<<h_Genpt_pythia->GetBinCenter(ibin+1)<<endl;
		// cout<<"h_FONLL = "<< h_FONLL->GetBinContent( h_FONLL->FindBin(h_Genpt_pythia->GetBinCenter(ibin+1)) ) << " , h_Genpt_pythia = "<<h_Genpt_pythia->GetBinContent(ibin+1)<<endl;


		if(ibin<nbin-4){   // this is wrong

		// FONLLWeight=h_FONLL->GetBinContent( h_FONLL->FindBin(h_Genpt_FONLLWeight->GetBinCenter(ibin+1)  ) ) / h_Genpt_pythia->GetBinContent( h_Genpt_pythia->FindBin( h_Genpt_FONLLWeight->GetBinCenter(ibin+1)  ));
		// FONLLWeight_mt=h_FONLL_Ds_byD0MTs->GetBinContent( h_FONLL_Ds_byD0MTs->FindBin(h_Genpt_pythia_mtbin->GetBinCenter(ibin+1)) ) / h_Genpt_pythia_mtbin->GetBinContent(ibin+1);
		FONLLWeight_mt=h_FONLL_Ds_byD0MTs->GetBinContent( h_FONLL_Ds_byD0MTs->FindBin(h_Genpt_FONLLWeight_mt->GetBinCenter(ibin+1)) ) / h_Genpt_pythia_mtbin->GetBinContent(  h_Genpt_pythia_mtbin->FindBin( h_Genpt_FONLLWeight_mt->GetBinCenter(ibin+1) ) );
		h_Genpt_FONLLWeight_mt->SetBinContent( ibin+1 , FONLLWeight_mt );	
		h_Genpt_FONLLWeight_mt->SetBinError(ibin+1 , 0 );	


		FONLLRaaWeight_mt=h_FONLL_Ds_byD0MTs->GetBinContent( h_FONLL_Ds_byD0MTs->FindBin(h_Genpt_FONLLRaaWeight_mt->GetBinCenter(ibin+1)) ) / h_Genpt_pythia_mtbin->GetBinContent( h_Genpt_pythia_mtbin->FindBin( h_Genpt_FONLLRaaWeight_mt->GetBinCenter(ibin+1) )  ) *gr_Raa->Eval(h_Genpt_FONLLRaaWeight_mt->GetBinCenter(ibin+1),0,"" )  ;
/*
		FONLLRaaWeight_mt=h_FONLL_Ds_byD0MTs->GetBinContent( h_FONLL_Ds_byD0MTs->FindBin(h_Genpt_pythia_mtbin->GetBinCenter(ibin+1)) ) / h_Genpt_pythia_mtbin->GetBinContent(ibin+1)* h_Raa->GetBinContent(h_Raa->FindBin(h_Genpt_pythia_mtbin->GetBinCenter(ibin+1))) ;
		if(PNPrompt==0){
		FONLLRaaWeight_mt=h_FONLL_Ds_byD0MTs->GetBinContent( h_FONLL_Ds_byD0MTs->FindBin(h_Genpt_pythia_mtbin->GetBinCenter(ibin+1)) ) / h_Genpt_pythia_mtbin->GetBinContent(ibin+1)*gr_DsRaa->Eval(h_Genpt_pythia_mtbin->GetBinCenter(ibin+1),0,"" );
	}
*/

		h_Genpt_FONLLRaaWeight_mt->SetBinContent( ibin+1 , FONLLRaaWeight_mt );	
		h_Genpt_FONLLRaaWeight_mt->SetBinError(ibin+1 , 0 );	


		// cout<<"\n mt weight : "<<FONLLWeight_mt<<", h_FONLL bin center = "<<h_FONLL_Ds_byD0MTs->GetBinCenter( h_FONLL_Ds_byD0MTs->FindBin(h_Genpt_pythia_mtbin->GetBinCenter(ibin+1)) )<< ", h_Genpt_pythia bin center = "<<h_Genpt_pythia_mtbin->GetBinCenter(ibin+1)<<endl;
		// cout<<"mt h_FONLL = "<< h_FONLL_Ds_byD0MTs->GetBinContent( h_FONLL_Ds_byD0MTs->FindBin(h_Genpt_pythia->GetBinCenter(ibin+1)) ) << " , h_Genpt_pythia = "<<h_Genpt_pythia_mtbin->GetBinContent(ibin+1)<<endl;

			
		} // end if ibin<nbin-4 for mtbin

	}

	// return 1;

	can[can_count]=new TCanvas(Form("can%i",can_count),Form("can%i",can_count),800,800);
	can[can_count]->cd();
	h_Genpt_FONLLWeight->Draw(); 
	can_count++;


	// return 1;


	can[can_count]=new TCanvas(Form("can%i",can_count),Form("can%i",can_count),800,800);
	can[can_count]->cd();
	h_Genpt_FONLLWeight_mt->Draw(); 
	can_count++;

	TGraph *gr_Genpt_FONLLWeight = new TGraph(); gr_Genpt_FONLLWeight->SetName("gr_Genpt_FONLLWeight"); 
	TGraph *gr_Genpt_FONLLWeight_mt = new TGraph(); gr_Genpt_FONLLWeight_mt->SetName("gr_Genpt_FONLLWeight_mt");
	TGraph *gr_Genpt_FONLLRaaWeight = new TGraph(); gr_Genpt_FONLLRaaWeight->SetName("gr_Genpt_FONLLRaaWeight");
	TGraph *gr_Genpt_FONLLRaaWeight_mt = new TGraph(); gr_Genpt_FONLLRaaWeight_mt->SetName("gr_Genpt_FONLLRaaWeight_mt");

	int i=1;

	// while(binNow<binHigh)
	for(auto binpt:v_bins){
		double pythiaBinContent=h_Genpt_pythia->GetBinContent( h_Genpt_pythia->FindBin(binpt));
		double RaaValue=gr_Raa->Eval(binpt);
		double FONLLValue=gr_FONLL->Eval(binpt);
//		double D0DataRelative=gr_D0Data_pp_weight

    FONLLWeight=FONLLValue / pythiaBinContent;
		FONLLRaaWeight=FONLLValue/ pythiaBinContent * RaaValue; 

		gr_Genpt_FONLLWeight->SetPoint(i, binpt , FONLLWeight);
		gr_Genpt_FONLLWeight_mt->SetPoint(i, binpt , FONLLWeight); 
		gr_Genpt_FONLLRaaWeight->SetPoint(i, binpt , FONLLRaaWeight); 
		gr_Genpt_FONLLRaaWeight_mt->SetPoint(i, binpt , FONLLRaaWeight);  

		++i;
	}

	// smooth histogram 


	f_out->cd();
	h_Genpt_FONLLWeight->Write("",TObject::kOverwrite);
	
	if(PNPrompt==0){
	h_Genpt_FONLLWeight_mt->Write("",TObject::kOverwrite);
	}

	if(isPbPb){
	  h_Genpt_FONLLRaaWeight->Write("",TObject::kOverwrite);
	  if(PNPrompt==0){
 	 h_Genpt_FONLLRaaWeight_mt->Write("",TObject::kOverwrite);
  }
		
	}




	gr_D0Data_pp_weight->Write("",TObject::kOverwrite);
	gr_D0Data_PbPb_weight->Write("",TObject::kOverwrite);


  gr_Genpt_FONLLWeight->Write("",TObject::kOverwrite);
  gr_Genpt_FONLLWeight_mt->Write("",TObject::kOverwrite);
  gr_Genpt_FONLLRaaWeight->Write("",TObject::kOverwrite);
  gr_Genpt_FONLLRaaWeight_mt->Write("",TObject::kOverwrite);

/*
  TF1 *f1_test=new TF1("f1_test","[0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*log(x)*log(x)*log(x)+ [4]*log(x)*log(x)*log(x)*log(x)");
  // TF1 *f1_test=new TF1("f1_test","[0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*log(x)*log(x)*log(x)");
//  TF1 *f1_test=new TF1("f1_test","[0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*log(x)*log(x)*log(x)");
  // TF1 *f1_test=new TF1("f1_test","[0]+[1]*x+[2]*x*x+[3]*x*x*x");
	f1_test->SetRange(1,50);
	f1_test->SetLineColor(2);
	gr_Genpt_FONLLWeight_mt->Fit("f1_test","LEMS0FI");
	gr_Genpt_FONLLWeight_mt->Fit("f1_test","LEMS0FI");
	gr_Genpt_FONLLWeight_mt->Fit("f1_test","LEMS0FI");
	gr_Genpt_FONLLWeight_mt->Fit("f1_test","LEMS0FI");
	gr_Genpt_FONLLWeight_mt->Fit("f1_test","LEMS0FI");
	gr_Genpt_FONLLWeight_mt->Fit("f1_test","LEMS0FI");

	TCanvas *c_fit=new TCanvas("c_fit","c_fit");
	c_fit->cd();
	gr_Genpt_FONLLWeight_mt->Draw("ap");
	f1_test->Draw("same");
*/
	// apply additional factor for BR to estimate CS

	cout<<"finish "<<str_isPbPb<<" , "<<str_PNPrompt<<" , "<<str_DsChannel<<endl;

	return 0;

}// end int getFONLLweight

