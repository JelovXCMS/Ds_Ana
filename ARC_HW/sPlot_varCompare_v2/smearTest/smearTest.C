#include "../varCompare_para.h"


int smearTest(int isPbPb=0, TString var_compare="DdlErr", TString var_cut="Dpt" , double var_cutLow=6, double var_cutHigh=10, double bins_var_Low=0, double bins_var_High=1, double bins_var_DrawLow=0, double bins_var_DrawHigh=0.1 ){

	TString Str_PbPb="pp";
	TString Str_PbPb3="pp";
	if(isPbPb){
		Str_PbPb="PbPb";
		Str_PbPb3="PbPb3";
	}	

  TString S_dataName=Form("../rootF/%s_fitFile.root",Str_PbPb3.Data());
  TString S_MCP=Form("../rootF/%sMC_phiPrompt_fitFile.root",Str_PbPb3.Data());
  TString S_MCNP=Form("../rootF/%sMC_phiNonPrompt_fitFile.root",Str_PbPb3.Data());
	TString S_fit=Form("../fitout/%s_%s_%s%.0fto%.0f.root",Str_PbPb.Data(),var_compare.Data(),var_cut.Data(),var_cutLow*100,var_cutHigh*100);


	TFile *f_data=TFile::Open(S_dataName);
	TFile *f_MCP=TFile::Open(S_MCP);
	TFile *f_MCNP=TFile::Open(S_MCNP);
	TFile *f_fit=TFile::Open(S_fit);

	TH1D *h_fr_PromptMC=(TH1D*)f_data->Get("h_fr_PromptMC");
	TH1D *h_DdlErr_MixMC=(TH1D*)f_data->Get("h_DdlErr_MixMC");

	TH1D *h_DdlErr_MixMC_test=(TH1D*)h_DdlErr_MixMC->Clone("h_DdlErr_MixMC_test");
	double content;
	double x, w ,z ;
	const UInt_t nbinx=h_DdlErr_MixMC->GetNbinsX();

`

	cout<<"test"<<endl;



	return 0;

}

int main(int argc, char*argv[]){

  if(argc==6){
    smearTest(atoi(argv[1]), argv[2], argv[3], atof(argv[4]), atof(argv[5]) );
  }else if(argc==10){
    smearTest(atoi(argv[1]), argv[2], argv[3], atof(argv[4]), atof(argv[5]) , atof(argv[6]) , atof(argv[7]) , atof(argv[8]) , atof(argv[9]));
  }else{
    smearTest();
    cout<<"wrong number of input parameters , need 6\n int isPbPb=0, TString var_compare=DdxyzErr, TString var_cut=Dtrk1Pt , double var_cutLow=0.75, double var_cutHigh=1.25"<<endl;
    return 1;
  }
//int Fit_sideband(int isPbPb=0, TString var_compare="DdxyzErr", TString var_cut="Dtrk1Pt" , double var_cutLow=0.75, double var_cutHigh=1.25){


  return 0;
}

