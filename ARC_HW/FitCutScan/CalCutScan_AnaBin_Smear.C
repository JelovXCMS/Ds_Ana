// double bins_var[]={3.5,4.0,4.5,5.0};
// const int nbin_var=4;

// double bins_var[]={1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0};
// const int nbin_var=11;
// double bins_var[]={1.5,2.5,3.5,4.5,5.5};
// const int nbin_var=5;


// double bins_var[]={0.08,0.11,0.14,0.17,0.2};
// const int nbin_var=5;

// double bins_var[]={1.5,2,2.5,3,3.5,4,4.5,5,5.5};
// const int nbin_var=9;
// double bins_var[]={3.25 ,3.5,3.75,4.0,4.25,4.5,4.75,5,5.25};
// const int nbin_var=9;

double PNPRatio_Glo=0.85;
double PNPRatioUp_Glo=0.00;
double PNPRatioDown_Glo=0.00;


// double bins_var[]={0.02,0.07,0.12,0.17,0.22 ,0.27};
// const int nbin_var=6;
// scan_val=(0.02 0.07 0.12 0.17 0.22)

// double bins_var[]={0.2,0.25,0.3,0.35};
// const int nbin_var=4;

// double bins_var[]={0.05,0.1,0.15,0.20,0.25,0.3,0.35,0.4,0.45};
// double bins_var[]={0.20,0.25,0.3,0.35,0.4};
// const int nbin_var=9;
// const int nbin_var=5;

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

// double textposx=0.2;
// double textposy=0.77;


Float_t SmearRelFactorArr[]={0.01,0.05,0.1,0.2,0.3};
Float_t SmearAbsFactorArr[]={0.0001,0.001,0.005,0.02,0.03,0.05,0.08,0.1,0.15};
Float_t ScaleErrFactorArr[]={0.97,1,1.05,1.15,1.29};
// Float_t SmearFactorArr[]={0.01,0.05};
const int nSmrRelF=sizeof(SmearRelFactorArr)/sizeof(SmearRelFactorArr[0]);
const int nSmrAbsF=sizeof(SmearAbsFactorArr)/sizeof(SmearAbsFactorArr[0]);
const int nSclErrF=sizeof(ScaleErrFactorArr)/sizeof(ScaleErrFactorArr[0]);
Float_t Ddls_SmrRelF[nSmrRelF];
Float_t DdlErr_SmrRelF[nSmrRelF];
Float_t Ddls_SmrAbsF[nSmrAbsF];
Float_t DdlErr_SmrAbsF[nSmrAbsF];
Float_t Ddls_SclErrF[nSclErrF];
Float_t DdlErr_SclErrF[nSclErrF];





double shiftY=0.06;
double shiftX=0.3;
double oneshift=0.075;

int reWeight=0;
TString s_reWeight="DdxyzErrWeight";


 // TH1D *h_Eff_fun(isPbPb,f_MCP,DptLow,DptHigh,Dalpha_cut, Dchi2cl_cut, Ddls_cut, DtktkResmass_cut){

// }
/*
	 TFile *MakeSmearFile(TFile *f_MCP,TString fname,int nSmear=2,double SmearFactor=0.1){
	 TRandom3 Rdm;
	 Rdm.SetSeed(0);
	 TTree *t_ori=(TTree*)f_MCP->Get("t_fit");

	 TFile *f_temp=new TFile(fname.Data(),"recreate");





	 }
	 */


int smearTree(TFile *,TFile *, int , double ,double );

TH1D *h_CSCal_fun(int isPbPb,TH1D *hRawYield, TH1D* hBtoDsCSdNdpt_temp,TFile *f_mc_Prompt,TFile *f_mc_NonPrompt,double DptLow,double DptHigh,double Dalpha_cut,double Dchi2cl_cut,double Ddls_cut,double DtktkResmass_cut,int fixPNPRatio,double PNPRatio,double PhiRatio, double &Eff_MCP, double &Eff_MCPerr, double &Eff_MCNP, double &Eff_MCNPerr, double &PromptYield, double &PromptYieldErr, TString s_Ddls="Ddls", TString s_SMrWt="1"){
	// hRawYield[i] is integral in pt ranges,
	// hBtoDsCSdNdpt_temp is 1/dpt,
	// calculate integral first, then get 1/dpt in the end

	double LumiNevt=LumiSum;
	if(isPbPb){
		LumiNevt=NevtPbPb3;
	}

	// get Efficiency

	//  TString DataCuts=Form("Dpt>%f && Dpt<%f && Dalpha < %f && Dchi2cl > %f && Ddls >%f && %s > %f && %s <%f", Dpt_Low, Dpt_High, Dalpha_cut, Dchi2cl_cut, Ddls_cut, var_cut.Data(), var_cutLow,var_cut.Data(), var_cutHigh); 
	TString DataCuts=Form("Dpt>%f && Dpt<%f && Dalpha<%f && Dchi2cl>%f && %s>%f && DtktkResmass>%f && DtktkResmass<%f", DptLow,DptHigh,Dalpha_cut,Dchi2cl_cut,s_Ddls.Data(),Ddls_cut, DtktkResmassCutMean-DtktkResmass_cut, DtktkResmassCutMean+DtktkResmass_cut);
	TString S_reWeight="1";
	if(reWeight==1){
		S_reWeight=s_reWeight;
	}

	TTree *t_MCP=(TTree*)f_mc_Prompt->Get("t_fit");
	TTree *t_MCNP=(TTree*)f_mc_NonPrompt->Get("t_fit");

	TH1D *h_reco_cut_MCP=new TH1D("h_reco_cut_MCP","h_reco_cut_MCP",40,0,40); h_reco_cut_MCP->Sumw2();
	t_MCP->Project("h_reco_cut_MCP","Dpt",(TCut)Form("TotalWeight*%s*%s*(%s)",S_reWeight.Data(),s_SMrWt.Data(),DataCuts.Data() ));


	TH1D *h_reco_cut_MCNP=new TH1D("h_reco_cut_MCNP","h_reco_cut_MCNP",40,0,40); h_reco_cut_MCNP->Sumw2();
	t_MCNP->Project("h_reco_cut_MCNP","Dpt",(TCut)Form("TotalWeight*%s*%s*(%s)",S_reWeight.Data(),s_SMrWt.Data(),DataCuts.Data() ));

	TH1D *hGen_MCP=(TH1D*)f_mc_Prompt->Get("hGen_pt");
	TH1D *hGen_MCNP=(TH1D*)f_mc_NonPrompt->Get("hGen_pt");

	double Nerr_reco_MCP;
	double Nerr_reco_MCNP;
	double Nerr_Gen_MCP;
	double Nerr_Gen_MCNP;

	double N_reco_MCP=h_reco_cut_MCP->IntegralAndError( h_reco_cut_MCP->FindBin(DptLow+0.001), h_reco_cut_MCP->FindBin(DptHigh-0.001) , Nerr_reco_MCP);
	double N_reco_MCNP=h_reco_cut_MCNP->IntegralAndError( h_reco_cut_MCNP->FindBin(DptLow+0.001), h_reco_cut_MCNP->FindBin(DptHigh-0.001) , Nerr_reco_MCNP);
	double N_Gen_MCP=hGen_MCP->IntegralAndError(hGen_MCP->FindBin(DptLow+0.001), hGen_MCP->FindBin(DptHigh-0.001) , Nerr_Gen_MCP);
	double N_Gen_MCNP=hGen_MCNP->IntegralAndError(hGen_MCNP->FindBin(DptLow+0.001), hGen_MCNP->FindBin(DptHigh-0.001) , Nerr_Gen_MCNP);

	// cout<<"N_reco_MCP = "<<N_reco_MCP<<endl;

	// double Eff_MCP=N_reco_MCP/N_Gen_MCP;
	// double Eff_MCNP=N_reco_MCNP/N_Gen_MCNP;

	Eff_MCP=N_reco_MCP/N_Gen_MCP;
	Eff_MCNP=N_reco_MCNP/N_Gen_MCNP;

	Eff_MCPerr=Eff_MCP * sqrt(Nerr_reco_MCP/N_reco_MCP * Nerr_reco_MCP/N_reco_MCP + Nerr_Gen_MCP/N_Gen_MCP *Nerr_Gen_MCP/N_Gen_MCP );
	Eff_MCNPerr=Eff_MCNP * sqrt(Nerr_reco_MCNP/N_reco_MCNP * Nerr_reco_MCNP/N_reco_MCNP + Nerr_Gen_MCNP/N_Gen_MCNP *Nerr_Gen_MCNP/N_Gen_MCNP );


	cout<<"Eff_MCP = "<<Eff_MCP<<" +- "<<Eff_MCPerr<<endl;
	cout<<"Eff_MCNP = "<<Eff_MCNP<<" +- "<<Eff_MCNPerr<<endl;


	TH1D *hBtoDsCSdNdpt_Integral=(TH1D*)hBtoDsCSdNdpt_temp->Clone("hBtoDsCSdNdpt_Integral");
	MutiplyBinWidth(hBtoDsCSdNdpt_Integral);

	double CS_BtoDs=hBtoDsCSdNdpt_Integral->Integral(hBtoDsCSdNdpt_Integral->FindBin(DptLow+0.001), hBtoDsCSdNdpt_Integral->FindBin(DptHigh-0.001))/(DptHigh-DptLow);
	double N_yield=hRawYield->GetBinContent(1)/(DptHigh-DptLow);
	double CS_PromptDs = (N_yield*PhiRatio/(2*LumiNevt) - (CS_BtoDs* (BRphi*Eff_MCNP ) ) ) / ( BRphi*Eff_MCP );
	cout<<"CS_PromptDs = "<<CS_PromptDs<<endl;

	double fr_Prompt=CS_PromptDs/(CS_PromptDs+CS_BtoDs);
	cout<<"fr_Prompt = "<<fr_Prompt<<endl;


	double N_yield_Up=(hRawYield->GetBinContent(1)+hRawYield->GetBinError(1))/(DptHigh-DptLow);
	double CS_PromptDs_Up = (N_yield_Up*PhiRatio/(2*LumiNevt) - (CS_BtoDs* (BRphi*Eff_MCNP ) ) ) / ( BRphi*Eff_MCP );
	double fr_Prompt_Up=CS_PromptDs_Up/(CS_PromptDs_Up+CS_BtoDs);

	double N_yield_Down=(hRawYield->GetBinContent(1)-hRawYield->GetBinError(1))/(DptHigh-DptLow);
	double CS_PromptDs_Down = (N_yield_Down*PhiRatio/(2*LumiNevt) - (CS_BtoDs* (BRphi*Eff_MCNP ) ) ) / ( BRphi*Eff_MCP );
	double fr_Prompt_Down=CS_PromptDs_Down/(CS_PromptDs_Down+CS_BtoDs);

	// double frErr_Prompt=fr_Prompt*hRawYield->GetBinError(1)/hRawYield->GetBinContent(1);


	// continued here 20190401
	if(!fixPNPRatio){
		PNPRatio_Glo=fr_Prompt;
		PNPRatioUp_Glo=fr_Prompt_Up;
		PNPRatioDown_Glo=fr_Prompt_Down;

		cout<<"PNPRatio_Glo = "<<PNPRatio_Glo<<endl;
		cout<<"fr_Prompt = "<<fr_Prompt<<" , fr_Prompt_Up = "<<fr_Prompt_Up<<" , fr_Prompt_Down = "<<fr_Prompt_Down<<endl;
		cout<<"RawYield Err = "<< hRawYield->GetBinError(1)<<"hRawYield = "<<hRawYield->GetBinContent(1)<<endl;

	}

	if(fixPNPRatio){
		CS_PromptDs=(N_yield*PhiRatio/(2*LumiNevt)) / ( ( BRphi*Eff_MCP ) + (1.0-PNPRatio)/(PNPRatio)*(BRphi*Eff_MCNP ) );
		cout<<"CS_PromptDs new = "<<CS_PromptDs<<" , PNPRatio = "<<PNPRatio<<endl;
	}

	PromptYield= CS_PromptDs*BRphi*Eff_MCP*2*LumiNevt;
	PromptYieldErr= PromptYield * hRawYield->GetBinError(1)/ hRawYield->GetBinContent(1);	


	TH1D *h_CSresult=new TH1D("h_CSresult","h_CSresult",1,0,1);
	h_CSresult->SetBinContent(1,CS_PromptDs);
	h_CSresult->SetBinError(1,CS_PromptDs*hRawYield->GetBinError(1)/hRawYield->GetBinContent(1));

	delete h_reco_cut_MCP;
	delete h_reco_cut_MCNP;
	delete hGen_MCP;
	delete hGen_MCNP;
	delete t_MCP;
	delete t_MCNP;

	return h_CSresult;
}


int CalCutScan_AnaBin_Smear(int isPbPb=0, int useAnaBin=1, int ibin_Dpt=7, double DptLow=20, double DptHigh=40, TString var_scan="Dalpha", int FixShape=1, int useReweight=0, TString ReweightTree="DdxyzErrWeight", int doPFrScan=0, int fixPFr=1, double PfrVal=1.0, int doSmear=1, int nSmear=10){

	// InitStyle();
	initParameter();
	setTDRStyle();
	InitStyle();

	double PhiMass=DtktkResmassCutMean;

	double *DalphaMax_bins=DalphaMax_bins_pp;
	double *Dchi2clMin_bins=Dchi2clMin_bins_pp;
	double *DdlsMin_bins=DdlsMin_bins_pp;
	double *bins_pt=bins_pt_pp;
	int nbin_pt=nbin_pt_pp;


	double Dalpha_cut=0.12;
	// double Dchi2cl_cut=0.1;
	// double Dchi2cl_cut=0.04;
	double Dchi2cl_cut=0.25;
	double Ddls_cut=2.5;
	double DtktkResmass_cut=DtktkResmassCutWidth;

	double Dalpha_cut_Def=0.12;
	double Dchi2cl_cut_Def=0.25;
	double Ddls_cut_Def=2.5;
	double DtktkResmass_cut_Def=DtktkResmassCutWidth;


	if(isPbPb==0 && DptHigh<=6){
		Dchi2cl_cut=0.1;
		Ddls_cut=3.5;
	}
	if(isPbPb==0 && DptLow>=6){
		Dchi2cl_cut=0.03;
		Ddls_cut=2.4;
	}
	if(isPbPb==3){
		Dchi2cl_cut=0.25;
		Ddls_cut=4.8;
	}
	// if(isPbPb==3 && DptHigh<=10){
	// Dchi2cl_cut=0.25;
	// Ddls_cut=4.8;
	// }
	// if(isPbPb==3 && DptLow>=10){
	// Dchi2cl_cut=0.08;
	// Ddls_cut=4.7;
	// }

	// import Default Cut
	// double DtktkResmass_cut=DtktkResmassCutWidth;  // 0.009



	if(isPbPb==3){
		nbin_pt=nbin_pt_PbPb3;
		DalphaMax_bins=DalphaMax_bins_PbPb3;
		Dchi2clMin_bins=Dchi2clMin_bins_PbPb3;
		DdlsMin_bins=DdlsMin_bins_PbPb3;
		bins_pt=bins_pt_PbPb3;
	}
	if(useAnaBin){
		Dalpha_cut_Def=DalphaMax_bins[ibin_Dpt];
		Dalpha_cut=DalphaMax_bins[ibin_Dpt];
		Dchi2cl_cut_Def=Dchi2clMin_bins[ibin_Dpt];
		Dchi2cl_cut=Dchi2clMin_bins[ibin_Dpt];
		Ddls_cut_Def=DdlsMin_bins[ibin_Dpt];
		Ddls_cut=DdlsMin_bins[ibin_Dpt];
		DptLow=bins_pt[ibin_Dpt];
		DptHigh=bins_pt[ibin_Dpt+1];
	}




	double bins_var_Dalpha[]={0.08,0.11,0.14,0.17,0.2};
	const int nbin_var_Dalpha=5;

	// need specail handle
	double bins_var_PhiMass[]={0.008,0.0085,0.009,0.0095,0.01,0.0105,0.0110};
	const int nbin_var_PhiMass=7;

	double bins_var_Ddls_pp[]={1.5, 2.5,3.5,4.5,5.5};
	const int nbin_var_Ddls_pp=5;

	double bins_var_Ddls_PbPb[]={2.5,3.5,4.5,5.5,6.5};
	const int nbin_var_Ddls_PbPb=5;


	double bins_var_Ddls_PbPb_lowpt[]={3.5,4.5,5.5,6.5};
	const int nbin_var_Ddls_PbPb_lowpt=4;

	double bins_var_Dchi2cl_pp[]={0.02,0.12,0.22,0.32,0.42};
	const int nbin_var_Dchi2cl_pp=5;

	double bins_var_Dchi2cl_PbPb[]={0.05,0.15,0.25,0.35,0.45};
	const int nbin_var_Dchi2cl_PbPb=5;

	double *bins_var=bins_var_Dalpha;
	int nbin_var=nbin_var_Dalpha;


	if(var_scan=="Ddls"){
		bins_var=bins_var_Ddls_pp;
		nbin_var=nbin_var_Ddls_pp;
		if(isPbPb && DptHigh<=10){
			bins_var=bins_var_Ddls_PbPb_lowpt;
			nbin_var=nbin_var_Ddls_PbPb_lowpt;
		}else if(isPbPb && DptLow>=10){
			bins_var=bins_var_Ddls_PbPb;
			nbin_var=nbin_var_Ddls_PbPb;
		}
	}
	else if(var_scan=="Dchi2cl"){
		bins_var=bins_var_Dchi2cl_pp;
		nbin_var=nbin_var_Dchi2cl_pp;
		if(isPbPb ){
			bins_var=bins_var_Dchi2cl_PbPb;
			nbin_var=nbin_var_Dchi2cl_PbPb;
		}	
	}
	else if(var_scan=="Dalpha"){
		bins_var=bins_var_Dalpha;
		nbin_var=nbin_var_Dalpha;
	}
	else if(var_scan=="PhiMass"){
		bins_var=bins_var_PhiMass;
		nbin_var=nbin_var_PhiMass;
	}else{
		cout<<"unkonw var for scan : "<<var_scan<<endl;
		return 2;
	}



	TLatex *tltx=new TLatex();

	TString s_FixShape="";
	if(FixShape){
		s_FixShape="FixShape";
	}

	TString s_PbPb="pp";
	TString s_PbPb3="pp";
	if(isPbPb){
		s_PbPb="PbPb";
		s_PbPb3="PbPb3";
	}


	TString s_reWeightFname="";
	if(useReweight){
		reWeight=1;
		s_reWeight=ReweightTree;
	}
	if(reWeight==1){
		s_reWeightFname="_"+s_reWeight;
		cout<<"s_reWeightFname = "<<s_reWeightFname<<endl;
	}

	// TString mcName_Prompt=Form("/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_cutScan/%sMC_phiPrompt_fitFile%s.root",s_PbPb3.Data(),s_reWeightFname.Data());
	TString mcName_Prompt=Form("/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_cutScan/%sMC_phiPrompt_fitFile%s.root",s_PbPb3.Data(),s_reWeightFname.Data());
	TString mcName_NonPrompt=Form("/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_cutScan/%sMC_phiNonPrompt_fitFile%s.root",s_PbPb3.Data(),s_reWeightFname.Data());

	TFile *f_MCP=TFile::Open(mcName_Prompt.Data());
	TFile *f_MCNP=TFile::Open(mcName_NonPrompt.Data());

	// TTree *t_MCP=(TTree*)f_MCP->Get("t_fit");
	// t_MCP->Draw("Dmass");


	TFile *f_fit[nbin_var];
	TFile *f_fit_Default;
	TH1D *h_RawYield[nbin_var];
	TH1D *h_RawYield_Default;
	int index_best=0;
	double bestErr=1;
	int fixPNPRatio=0;
	double PNPRatio=0.9;
	double PhiRatio=0.95;	
	if(isPbPb){
		PhiRatio=0.91;
	}

	cout<<"isPbPb = "<<isPbPb<<endl;
	cout<<"useAnaBin = "<<useAnaBin<<endl;
	cout<<"ibin_Dpt = "<<ibin_Dpt<<endl;
	cout<<"DptLow = "<<DptLow<<endl;
	cout<<"DptHigh = "<<DptHigh<<endl;
	cout<<"var_scan = "<<var_scan<<endl;
	cout<<"FixShape = "<<FixShape<<endl;
	cout<<"useReweight = "<<useReweight<<endl;
	cout<<"ReweightTree = "<<ReweightTree<<endl;
	cout<<"doPFrScan = "<<doPFrScan<<endl;
	cout<<"fixPFr = "<<fixPFr<<endl;
	cout<<"PfrVal = "<<PfrVal<<endl;
	cout<<"PhiRatio = "<<PhiRatio<<endl;

	// end 0



	// 1. read FitYield and Get The Best
	f_fit_Default=TFile::Open(Form("./fitout/%s_%sDpt%.0fto%.0f_Default.root", s_PbPb.Data(),s_FixShape.Data(),DptLow,DptHigh ));
	h_RawYield_Default=(TH1D*)f_fit_Default->Get("h_RawYield");

	for(int i=0;i<nbin_var; i++){
		f_fit[i]=TFile::Open(Form("./fitout/%s/%s_%sDpt%.0fto%.0f_%s%.0f.root",var_scan.Data(), s_PbPb.Data(),s_FixShape.Data(),DptLow,DptHigh, var_scan.Data(),bins_var[i]*100 ));
		h_RawYield[i]=(TH1D*)f_fit[i]->Get("h_RawYield");

		// h_RawYield[i]->Draw();
		if(h_RawYield[i]->GetBinError(1)/h_RawYield[i]->GetBinContent(1) <bestErr){
			index_best=i;
			bestErr=h_RawYield[i]->GetBinError(1)/h_RawYield[i]->GetBinContent(1);
		}
		cout<<"i  = "<<i<<" , RelErr = "<<h_RawYield[i]->GetBinError(1)/h_RawYield[i]->GetBinContent(1)<<endl;

	}
	cout<<"index_best = "<<index_best<<endl;
	cout<<"bestErr = "<<bestErr<<endl;
	cout<<"bestVar = "<<bins_var[index_best]<<endl;

	// index_best=4;

	// 2. use Best to Calculate PromptFr

	TH1D *h_CS[nbin_var];
	TH1D *h_CS_Default;

	// TH1D *h_EffMCP[nbin_var];
	// TH1D *h_EffMCP_Default;

	// TH1D *h_EffMCP[nbin_var];
	// TH1D *h_EffMCP_Default;

	double Eff_MCP_Default;
	double EffErr_MCP_Default;  
	double Eff_MCNP_Default;
	double EffErr_MCNP_Default;

	double Eff_MCP[nbin_var];
	double EffErr_MCP[nbin_var];  
	double Eff_MCNP[nbin_var];
	double EffErr_MCNP[nbin_var];

	double RawYield_Default;
	double RawYieldErr_Default;

	double PromptYield_Default;
	double PromptYieldErr_Default;

	// double RawYield[nbin_var];
	// double RawYieldErr[nbin_var];

	double RawYield;
	double RawYieldErr;

	double PromptYield[nbin_var];
	double PromptYieldErr[nbin_var];

	TString s_fixPFr="";
	if(fixPFr){
		fixPNPRatio=1;
		PNPRatio=PfrVal;
		PNPRatio_Glo=PfrVal;
		s_fixPFr=Form("_fixPFr%.0f",PfrVal*100);
	}

	// c_SingleRatio->SaveAs(Form("./Ratioplots/%s%s/%s_%s_SingleRatio_Dpt%.0fto%.0f%s.png", var_scan.Data(),s_reWeightFname.Data() , s_PbPb.Data(),var_scan.Data(),DptLow,DptHigh,s_FixShape.Data()));
	gSystem->Exec(Form("mkdir -p RatioOut/%s%s/", var_scan.Data(),s_reWeightFname.Data()));
	// TString foutName=Form("./RatioOut/%s%s/%s_%s%s_Dpt%.0fto%.0f%s.root", var_scan.Data(),s_reWeightFname.Data() , s_PbPb.Data(),var_scan.Data(),s_fixPFr.Data(),DptLow,DptHigh,s_FixShape.Data() );
	TString s_PFrScan="";
	if(doPFrScan){
		s_PFrScan="_PFrScan";
	}
	TString foutName=Form("./RatioOut/%s%s/%s_%s%s_Dpt%.0fto%.0f%s%s.root", var_scan.Data(),s_reWeightFname.Data() , s_PbPb.Data(),var_scan.Data(),s_fixPFr.Data(),DptLow,DptHigh,s_FixShape.Data(),s_PFrScan.Data() );
	TFile *fout=new TFile(foutName.Data(),"recreate");	


	TFile *f_NonPromptDs=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BtoDs_Results/output/BtoDs.root");
	TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight");
	TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight");
	TH1D *hBtoDs_AnaBin_pythiaWeight=hBtoDsCrossSectionPP_AnaBin_pythiaWeight;  // this is differential cross section
	if(isPbPb){ hBtoDs_AnaBin_pythiaWeight=hBtoDsdNdPtPbPb_AnaBin_pythiaWeight; }


	TString s_XTitle="Decay Length Significance";
	double DefaultCutVal=0;

	if(var_scan=="Ddls"){
		Ddls_cut=bins_var[index_best];
		DefaultCutVal=DdlsMin_bins[ibin_Dpt];
		s_XTitle="Decay Length Significance";
	}
	else if(var_scan=="Dchi2cl"){
		Dchi2cl_cut=bins_var[index_best];
		DefaultCutVal=Dchi2clMin_bins[ibin_Dpt];
		s_XTitle="Vertex Probability";
	}
	else if(var_scan=="Dalpha"){
		Dalpha_cut=bins_var[index_best];
		DefaultCutVal=DalphaMax_bins[ibin_Dpt];
		s_XTitle="Pointing Angle";
	}
	else if(var_scan=="PhiMass"){
		DtktkResmass_cut=bins_var[index_best];
		DefaultCutVal=DtktkResmassCutWidth;
		s_XTitle="#phi Mass window";
	}else{
		cout<<"unkonw var for scan : "<<var_scan<<endl;
		return 2;
	}





	if(doSmear){

		// create smeared tree
		TString mcName_Prompt_new=Form("/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_cutScan/%sMC_phiPrompt_fitFile%s_smear.root",s_PbPb3.Data(),s_reWeightFname.Data());
		TString mcName_NonPrompt_new=Form("/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_cutScan/%sMC_phiNonPrompt_fitFile%s_smear.root",s_PbPb3.Data(),s_reWeightFname.Data());
		TFile *f_mcp_new=new TFile(mcName_Prompt_new.Data(),"recreate");
		smearTree(f_mcp_new,f_MCP,nSmear,DptLow,DptHigh);
		TFile *f_mcnp_new=new TFile(mcName_NonPrompt_new.Data(),"recreate");
		smearTree(f_mcnp_new,f_MCNP,nSmear,DptLow,DptHigh);


		TGraphErrors *gr_DdlsScaleScan=new TGraphErrors();
		gr_DdlsScaleScan->SetName("gr_DdlsScaleScan");

		double Scale_Best=0;
		double Sys_Best=100;
		double Scale_Best_Pfr=0;

		// for(int m=0;m<nSmrAbsF;m++)
		for(int m=0;m<nSclErrF;m++)
		{

		// t_new->Branch(Form("Ddls_Abs%.0fem5",SmearAbsFactorArr[i]*1e5),&Ddls_SmrAbsF[i]);
		// run the result
		fixPNPRatio=0;
		if(fixPFr){
			fixPNPRatio=1;
			PNPRatio_Glo=PfrVal;
		}	

		// t_new->Branch(Form("Ddls_Scl%.0fem3",ScaleErrFactorArr[i]*1e3),&Ddls_SclErrF[i]);
		TString s_Ddls=Form("Ddls_Scl%.0fem3",ScaleErrFactorArr[m]*1e3);
		TString s_Ddls_text=Form("Ddls Err. Scale : %.2f",ScaleErrFactorArr[m]);
		// TString s_Ddls=Form("Ddls_Abs%.0fem5",SmearAbsFactorArr[m]*1e5);
		// TString s_Ddls=Form("Ddls_Rel%.0fem5",SmearRelFactorArr[m]*1e5);
		TString s_SMrWt="SMrWt";
		h_CS_Default=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_AnaBin_pythiaWeight ,f_mcp_new, f_mcnp_new, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default, s_Ddls,s_SMrWt);
		RawYield_Default=h_RawYield_Default->GetBinContent(1);
		RawYieldErr_Default=h_RawYield_Default->GetBinError(1);

		 cout<<"PNPRatio_Glo = "<<PNPRatio_Glo<<endl;

		fixPNPRatio=1;
		// PNPRatio_Glo=1;
		for(int i=0; i<nbin_var; i++){
			if(i==index_best && !useAnaBin)continue;

			if(var_scan=="Ddls"){
				Ddls_cut=bins_var[i];
			}
			else if(var_scan=="Dchi2cl"){
				Dchi2cl_cut=bins_var[i];
			}
			else if(var_scan=="Dalpha"){
				Dalpha_cut=bins_var[i];
			 }
			else if(var_scan=="PhiMass"){
				DtktkResmass_cut=bins_var[i];
			}else{
				cout<<"unkonw var for scan : "<<var_scan<<endl;
				return 2;
			}

			h_CS[i]=h_CSCal_fun(isPbPb,h_RawYield[i],hBtoDs_AnaBin_pythiaWeight ,f_mcp_new, f_mcnp_new, DptLow,DptHigh,Dalpha_cut,Dchi2cl_cut,Ddls_cut,DtktkResmass_cut, fixPNPRatio, PNPRatio_Glo, PhiRatio, Eff_MCP[i],EffErr_MCP[i],Eff_MCNP[i],EffErr_MCNP[i],PromptYield[i],PromptYieldErr[i],s_Ddls,s_SMrWt);
		}

		TGraphErrors *grMCP=new TGraphErrors();
		TGraphErrors *grMCNP=new TGraphErrors();
		TGraphErrors *grDataAll=new TGraphErrors();
		TGraphErrors *grDataPrompt=new TGraphErrors();

		TGraphErrors *gr_CS_Default=new TGraphErrors();
		gr_CS_Default->SetName("gr_CS_Default");


		TGraphErrors *gr_CS=new TGraphErrors();
		gr_CS->SetName("gr_CS");

		TGraphErrors *gr=new TGraphErrors();
		gr->SetName("gr_DoubleRatio");
		grMCP->SetName("grMCP");
		grMCNP->SetName("grMCNP");
		grDataAll->SetName("grDataAll");
		grDataPrompt->SetName("grDataPrompt");

		// double BestCS=h_CS[index_best]->GetBinContent(1);
		// double BestCSRelErr=h_CS[index_best]->GetBinError(1)/h_CS[index_best]->GetBinContent(1);

		double CS_Default=h_CS_Default->GetBinContent(1);
		 double CSRelErr_Default=h_CS_Default->GetBinError(1)/h_CS_Default->GetBinContent(1);
		double CSErr_Default=h_CS_Default->GetBinError(1);


		gr_CS_Default->SetPoint(0,DefaultCutVal,CS_Default);
		gr_CS_Default->SetPointError(0,0,CSRelErr_Default*CS_Default);

		double Syst=0;
		double SystError=0;

		for(int i=0; i<nbin_var; i++){
			// if(i==index_best)continue;
			//    if(var_scan=="Ddls"){
			//      Ddls_cut=bins_var[i];
			//    }
			double CS=h_CS[i]->GetBinContent(1);
			double CSRelErr=h_CS[i]->GetBinError(1)/h_CS[i]->GetBinContent(1);
			double CSErr=h_CS[i]->GetBinError(1);
			// double CSRatio=CS/BestCS;
			// double CSRatioErr=sqrt(abs(CSRelErr*CSRelErr-BestCSRelErr*BestCSRelErr))*CSRatio;

			double CSRatio=CS/CS_Default;
			double CSRatioErr=sqrt(abs(CSErr*CSErr-CSErr_Default*CSErr_Default))/CS_Default;
			// double CSRatioErr=sqrt(abs(CSRelErr*CSRelErr-CSRelErr_Default*CSRelErr_Default))*CSRatio;
			// double CSRatioErr=CSRelErr*CSRatio;

			if(abs(CSRatio-1) > Syst ){
				Syst=abs(CSRatio-1);
				SystError=CSRelErr;	
			}

			gr_CS->SetPoint(i,bins_var[i],CS);
			gr_CS->SetPointError(i,0,CSRelErr*CS);

			gr->SetPoint(i,bins_var[i] ,CSRatio);
			gr->SetPointError(i,0.,CSRatioErr);		
 /*
			double EffRatio_MCP=Eff_MCP[i]/Eff_MCP_Default;
			double EffRatioErr_MCP=EffRatio_MCP * sqrt(abs(EffErr_MCP[i]/Eff_MCP[i]*EffErr_MCP[i]/Eff_MCP[i] - EffErr_MCP_Default/Eff_MCP_Default* EffErr_MCP_Default/Eff_MCP_Default ) );

			grMCP->SetPoint(i,bins_var[i], EffRatio_MCP);
			grMCP->SetPointError(i,0., EffRatioErr_MCP);

			double EffRatio_MCNP=Eff_MCNP[i]/Eff_MCNP_Default;
			double EffRatioErr_MCNP=EffRatio_MCNP * sqrt(abs(EffErr_MCNP[i]/Eff_MCNP[i]*EffErr_MCNP[i]/Eff_MCNP[i] - EffErr_MCNP_Default/Eff_MCNP_Default* EffErr_MCNP_Default/Eff_MCNP_Default ) );

			grMCNP->SetPoint(i,bins_var[i], EffRatio_MCNP);
			grMCNP->SetPointError(i,0., EffRatioErr_MCNP);

			RawYield=h_RawYield[i]->GetBinContent(1);
			RawYieldErr=h_RawYield[i]->GetBinError(1);

			double RawYield_Ratio=RawYield/RawYield_Default;
			double RawYield_RatioErr=RawYield_Ratio * sqrt( abs( RawYieldErr/RawYield*RawYieldErr/RawYield - RawYieldErr_Default/RawYield_Default* RawYieldErr_Default/RawYield_Default ) );

			grDataAll->SetPoint(i, bins_var[i], RawYield_Ratio);
			grDataAll->SetPointError(i, 0., RawYield_RatioErr);


			double PromptYield_Ratio=PromptYield[i]/PromptYield_Default;
			double PromptYield_RatioErr=PromptYield_Ratio * sqrt( abs( PromptYieldErr[i]/PromptYield[i]*PromptYieldErr[i]/PromptYield[i] - PromptYieldErr_Default/PromptYield_Default* PromptYieldErr_Default/PromptYield_Default ) );

			grDataPrompt->SetPoint(i, bins_var[i], PromptYield_Ratio);
			grDataPrompt->SetPointError(i, 0., PromptYield_RatioErr);
*/

		}

		 // gr->GetXaxis()->SetRangeUser(0,4.5);
		gStyle->SetOptStat(0);
		int nbinTemp=1;
		double binsTempLow=0;
		double binsTempHigh=5.5;
		if(isPbPb){binsTempHigh=5.5;}
		if(var_scan=="Dchi2cl"){
			binsTempHigh=0.3;
			if(isPbPb){
				binsTempHigh=0.5;
			}
		}

		if(var_scan=="PhiMass"){
			binsTempLow=bins_var[0] - 0.5*(bins_var[1]- bins_var[0]) ;
		}

		// test binHigh 
		binsTempHigh=bins_var[nbin_var-1] + 0.5*(bins_var[nbin_var-1]- bins_var[nbin_var-2]) ;

/*
		TCanvas *c_SingleRatio= new TCanvas("c_SingleRatio","c_SingleRatio",800,800);
		c_SingleRatio->cd();

		TH1D *htemp2=new TH1D("htemp2",Form(";%s; Varied cut/Default cut",s_XTitle.Data()),nbinTemp,binsTempLow,binsTempHigh);
		htemp2->SetMaximum(2.5);
		htemp2->SetMinimum(0);
		htemp2->Draw();
		grMCP->SetLineColor(2);
		grMCP->SetMarkerColor(2);
		grMCP->Draw("p");
		grMCNP->SetLineColor(6);
		grMCNP->SetMarkerColor(6);
		grMCNP->Draw("p");
		grDataAll->SetLineColor(1);
		 grDataAll->SetMarkerColor(1);
		grDataAll->Draw("p");
		grDataPrompt->SetLineColor(4);
		grDataPrompt->SetMarkerColor(4);
		grDataPrompt->Draw("p");

		TLegend *le_SingleRatio=new TLegend(0.5,0.55,0.8,0.8);
		le_SingleRatio->SetBorderSize(0);
		le_SingleRatio->AddEntry(grMCP,"MC Prompt D_{S}^{+}","p");
		le_SingleRatio->AddEntry(grMCNP,"MC NonPrompt D_{S}^{+}","p");
		 le_SingleRatio->AddEntry(grDataAll,"Data D_{S}^{+} Yield","p");
		le_SingleRatio->AddEntry(grDataPrompt,"Data Prompt D_{S}^{+} Yield","p");
		le_SingleRatio->Draw("same");
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s : %.0f < D_{S} p_{T} < %.0f",s_PbPb.Data(),DptLow,DptHigh));

		gSystem->Exec(Form("mkdir -p Ratioplots/%s%s/",var_scan.Data(),s_reWeightFname.Data()));
		c_SingleRatio->SaveAs(Form("./Ratioplots/%s%s/%s_%s_SingleRatio%s_Dpt%.0fto%.0f%s.png", var_scan.Data(),s_reWeightFname.Data() , s_PbPb.Data(),var_scan.Data(),s_fixPFr.Data(),DptLow,DptHigh,s_FixShape.Data()));
*/
		TCanvas *c_out=new TCanvas("c_out","c_out",800,800);
		c_out->cd();

		TH1D *htemp=new TH1D("htemp",Form(";%s; #sigma_{Varied cut}/#sigma_{Default cut}",s_XTitle.Data()),nbinTemp,binsTempLow,binsTempHigh);
		htemp->SetBinContent(1,1);
		htemp->SetBinError(1,h_CS[index_best]->GetBinError(1)/h_CS[index_best]->GetBinContent(1));
		htemp->SetFillStyle(3001);
		htemp->SetFillColor(kGray);
		htemp->SetMarkerSize(0);
		htemp->SetMarkerColor(kGray);
		htemp->SetLineWidth(0);
		htemp->SetMaximum(2);
		htemp->SetMinimum(0);
		htemp->Draw("E2");
		gr->SetMarkerStyle(21);
		 gr->SetMarkerSize(0.4);
		gr->Draw("p");

		TLine *tl=new TLine(DefaultCutVal,0,DefaultCutVal,2);
		tl->SetLineColor(4);
		tl->SetLineStyle(6);
		tl->Draw("same");

		// grMCP->Draw("p");
		// grMCNP->Draw("p");

		// TF1 *f1=new TF1("f1","[0]+x*[1]");
		TFitResultPtr fitR;
		double x0[1]={0};
		double x0err[1];
		TF1 *f1=new TF1("f1","[0]+(x-[2])*[1]");
		f1->SetLineColor(2);
		f1->SetParameter(0,1);
		// f1->FixParameter(0,1);
		f1->SetParameter(1,0);
		f1->SetParameter(2,DefaultCutVal);
		// f1->FixParameter(2,bins_var[index_best]);
		f1->SetRange(0,binsTempHigh);

		gr->Fit("f1","LEMS0F");
		gr->Fit("f1","LEMS0F");
		fitR=gr->Fit("f1","LEMIS0F");
		// fitR=gr->Fit("f1","EMIS0");
		fitR->GetConfidenceIntervals(1,1,1,x0,x0err,0.683,false);

		TLine *tl1=new TLine(binsTempLow,1,binsTempHigh,1);
		tl1->SetLineColor(2);
		if(var_scan=="Ddls" || var_scan=="Dchi2cl"){	
			f1->Draw("same");
			Syst=(abs(f1->Eval(0)-1));
			SystError=x0err[0];
		 }else{
			tl1->Draw("same");
		}

		if(Syst<Sys_Best){
			Sys_Best=Syst;
			Scale_Best=ScaleErrFactorArr[m];	
			Scale_Best_Pfr=PNPRatio_Glo;
		}



		gStyle->SetOptFit(0);


		shiftY=0.06;

		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s : %.0f < D_{S} p_{T} < %.0f",s_PbPb.Data(),DptLow,DptHigh)); shiftY-=oneshift;

		if(var_scan=="Ddls" || var_scan=="Dchi2cl"){
			tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Systematics : %.1f%% #pm %.1f%%",Syst*100, SystError*100  )); shiftY-=oneshift;
			tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Slope : %.3f #pm %.3f",f1->GetParameter(1), f1->GetParError(1)) ); shiftY-=oneshift;
		}else{
			tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Systematics : %.1f%% ",Syst*100 )); shiftY-=oneshift;
		}
		if(fixPFr){
			tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Fix Prompt Ratio =  %.2f",PfrVal )); shiftY-=oneshift;
		}else{
			tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Default Prompt Ratio =  %.2f",PNPRatio_Glo )); shiftY-=oneshift;
		}
			tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s ",s_Ddls_text.Data() )); shiftY-=oneshift;



		shiftY=0.06;
		gSystem->Exec(Form("mkdir -p Ratioplots_smear/%s%s/",var_scan.Data(),s_reWeightFname.Data()));
		c_out->SaveAs(Form("./Ratioplots_smear/%s%s/%s_%s%s_Dpt%.0fto%.0f%s_%s.png",var_scan.Data(),s_reWeightFname.Data(),s_PbPb.Data(),var_scan.Data(),s_fixPFr.Data(),DptLow,DptHigh,s_FixShape.Data(), s_Ddls.Data() ));



		gr_DdlsScaleScan->SetPoint(m,ScaleErrFactorArr[m],Syst);
		gr_DdlsScaleScan->SetPointError(m,0.0,SystError);
		


	} // end for nSmrAbsF

	TCanvas *c_DdlsScaleScan=new TCanvas("c_DdlsScaleScan","",600,600);
	c_DdlsScaleScan->cd();
	gr_DdlsScaleScan->GetXaxis()->SetTitle("DdlErr Scale");
	gr_DdlsScaleScan->GetYaxis()->SetTitle("Syst.");
	gr_DdlsScaleScan->SetMaximum(0.8);
	gr_DdlsScaleScan->SetMinimum(-0.2);
	gr_DdlsScaleScan->Draw("ap");

	shiftY=0.06;
	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s : %.0f < D_{S} p_{T} < %.0f",s_PbPb.Data(),DptLow,DptHigh)); shiftY-=oneshift;
	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Best Scale : %.2f",Scale_Best)); shiftY-=oneshift;


	shiftY=0.06;

	c_DdlsScaleScan->SaveAs(Form("./Ratioplots_smear/%s%s/%s_%s%s_Dpt%.0fto%.0f%s_DdlErrScaleScan.png",var_scan.Data(),s_reWeightFname.Data(),s_PbPb.Data(),var_scan.Data(),s_fixPFr.Data(),DptLow,DptHigh,s_FixShape.Data()));

	fout->cd();
	gr_DdlsScaleScan->Write("",TObject::kOverwrite);

	cout<<"\nsmearing result "<<endl;
	cout<<"Best Scale : "<<Scale_Best<<endl;
	cout<<"Best Scale Pfr : "<<Scale_Best_Pfr<<endl;

	TGraphErrors *gr_ScaleBest=new TGraphErrors();
	gr_ScaleBest->SetName("gr_ScaleBest");
	gr_ScaleBest->SetPoint(0,(DptLow+DptHigh)/2,Scale_Best);
	gr_ScaleBest->SetPointError(0,(DptLow+DptHigh)/2-DptLow,0);
	gr_ScaleBest->Write();

	TGraphErrors *gr_ScaleBestPfr=new TGraphErrors();
	gr_ScaleBestPfr->SetName("gr_ScaleBestPfr");
	gr_ScaleBestPfr->SetPoint(0,(DptLow+DptHigh)/2,Scale_Best_Pfr);
	gr_ScaleBestPfr->SetPointError(0,(DptLow+DptHigh)/2-DptLow,0);
	gr_ScaleBestPfr->Write();



	} // end doSmear

///////////////////////////////
//// end do Smear /////////////
///////////////////////////////

	// return 1;

	/////////////////////////
	/// for PFr syst error //
	/////////////////////////

	// do the BtoD sys estimation before default, since we need PNPRatio_Glo for later systematics

	// FONLL shape with Data CS

	// if()

	fixPNPRatio=0;
	TH1D *hBtoDsCrossSectionPP_FONLLShape=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_FONLLShape");
	TH1D *hBtoDsdNdPtPbPb_FONLLShape=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_FONLLShape");
	TH1D *hBtoDs_FONLLShape=hBtoDsCrossSectionPP_FONLLShape;
	if(isPbPb){
		hBtoDs_FONLLShape=hBtoDsdNdPtPbPb_FONLLShape;
	}

	TH1D *h_PromptDs_BtoDsFONLLShape=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_FONLLShape ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);	
	double FONLL_PNPRatio=PNPRatio_Glo;

	cout<<"FONLL PNPRatio_Glo = "<<PNPRatio_Glo<<endl;

	// BtoD total err up & down

	TH1D *hBtoDsCrossSectionPP_Errup=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_Errup");
	TH1D *hBtoDsdNdPtPbPb_Errup=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Errup");
	TH1D *hBtoDs_Errup=hBtoDsCrossSectionPP_Errup;
	if(isPbPb){
		hBtoDs_Errup=hBtoDsdNdPtPbPb_Errup;
	}

	TH1D *h_PromptDs_Errup=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_Errup ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);
	double BtoDErrup_PNPRatio=PNPRatio_Glo;

	TH1D *hBtoDsCrossSectionPP_Errdown=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_Errdown");
	TH1D *hBtoDsdNdPtPbPb_Errdown=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_Errdown");
	TH1D *hBtoDs_Errdown=hBtoDsCrossSectionPP_Errdown;
	if(isPbPb){
		hBtoDs_Errdown=hBtoDsdNdPtPbPb_Errdown;
	}

	TH1D *h_PromptDs_Errdown=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_Errdown ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);
	double BtoDErrdown_PNPRatio=PNPRatio_Glo;


	// Scale Br up & down seperately 


	TH1D *hBtoDsCrossSectionPP_ScaleBR_B0toDmax=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDmax");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_B0toDmax=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDmax");
	TH1D *hBtoDs_ScaleBR_B0toDmax=hBtoDsCrossSectionPP_ScaleBR_B0toDmax;
	if(isPbPb){
		hBtoDs_ScaleBR_B0toDmax=hBtoDsdNdPtPbPb_ScaleBR_B0toDmax;
	}

	TH1D *h_PromptDs_ScaleBR_B0toDmax=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_ScaleBR_B0toDmax ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);
	double B0toDmax_PNPRatio=PNPRatio_Glo;


	TH1D *hBtoDsCrossSectionPP_ScaleBR_B0toDmin=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDmin");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_B0toDmin=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDmin");
	TH1D *hBtoDs_ScaleBR_B0toDmin=hBtoDsCrossSectionPP_ScaleBR_B0toDmin;
	if(isPbPb){
		hBtoDs_ScaleBR_B0toDmin=hBtoDsdNdPtPbPb_ScaleBR_B0toDmin;
	}

	TH1D *h_PromptDs_ScaleBR_B0toDmin=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_ScaleBR_B0toDmin ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);
	double B0toDmin_PNPRatio=PNPRatio_Glo;

	//

	TH1D *hBtoDsCrossSectionPP_ScaleBR_BptoDmax=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDmax");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_BptoDmax=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDmax");
	TH1D *hBtoDs_ScaleBR_BptoDmax=hBtoDsCrossSectionPP_ScaleBR_BptoDmax;
	if(isPbPb){
		hBtoDs_ScaleBR_BptoDmax=hBtoDsdNdPtPbPb_ScaleBR_BptoDmax;
	}
	TH1D *h_PromptDs_ScaleBR_BptoDmax=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_ScaleBR_BptoDmax ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);
	double BptoDmax_PNPRatio=PNPRatio_Glo;


	TH1D *hBtoDsCrossSectionPP_ScaleBR_BptoDmin=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDmin");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_BptoDmin=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDmin");
	TH1D *hBtoDs_ScaleBR_BptoDmin=hBtoDsCrossSectionPP_ScaleBR_BptoDmin;
	if(isPbPb){
		hBtoDs_ScaleBR_BptoDmin=hBtoDsdNdPtPbPb_ScaleBR_BptoDmin;
	}

	TH1D *h_PromptDs_ScaleBR_BptoDmin=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_ScaleBR_BptoDmin ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);
	double BptoDmin_PNPRatio=PNPRatio_Glo;


	//

	TH1D *hBtoDsCrossSectionPP_ScaleBR_BstoDmax=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDmax");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_BstoDmax=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDmax");
	TH1D *hBtoDs_ScaleBR_BstoDmax=hBtoDsCrossSectionPP_ScaleBR_BstoDmax;
	if(isPbPb){
		hBtoDs_ScaleBR_BstoDmax=hBtoDsdNdPtPbPb_ScaleBR_BstoDmax;
	}
	TH1D *h_PromptDs_ScaleBR_BstoDmax=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_ScaleBR_BstoDmax ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);
	double BstoDmax_PNPRatio=PNPRatio_Glo;


	TH1D *hBtoDsCrossSectionPP_ScaleBR_BstoDmin=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDmin");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_BstoDmin=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDmin");
	TH1D *hBtoDs_ScaleBR_BstoDmin=hBtoDsCrossSectionPP_ScaleBR_BstoDmin;
	if(isPbPb){
		hBtoDs_ScaleBR_BstoDmin=hBtoDsdNdPtPbPb_ScaleBR_BstoDmin;
	}
	TH1D *h_PromptDs_ScaleBR_BstoDmin=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_ScaleBR_BstoDmin ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);
	double BstoDmin_PNPRatio=PNPRatio_Glo;

	/// BtoDs

	TH1D *hBtoDsCrossSectionPP_ScaleBR_B0toDsmax=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDsmax");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_B0toDsmax=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDsmax");
	TH1D *hBtoDs_ScaleBR_B0toDsmax=hBtoDsCrossSectionPP_ScaleBR_B0toDsmax;
	if(isPbPb){
		hBtoDs_ScaleBR_B0toDsmax=hBtoDsdNdPtPbPb_ScaleBR_B0toDsmax;
	}

	TH1D *h_PromptDs_ScaleBR_B0toDsmax=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_ScaleBR_B0toDsmax ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);
	double B0toDsmax_PNPRatio=PNPRatio_Glo;


	TH1D *hBtoDsCrossSectionPP_ScaleBR_B0toDsmin=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_B0toDsmin");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_B0toDsmin=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_B0toDsmin");
	TH1D *hBtoDs_ScaleBR_B0toDsmin=hBtoDsCrossSectionPP_ScaleBR_B0toDsmin;
	if(isPbPb){
		hBtoDs_ScaleBR_B0toDsmin=hBtoDsdNdPtPbPb_ScaleBR_B0toDsmin;
	}

	TH1D *h_PromptDs_ScaleBR_B0toDsmin=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_ScaleBR_B0toDsmin ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);
	double B0toDsmin_PNPRatio=PNPRatio_Glo;

	//


	TH1D *hBtoDsCrossSectionPP_ScaleBR_BptoDsmax=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDsmax");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_BptoDsmax=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDsmax");
	TH1D *hBtoDs_ScaleBR_BptoDsmax=hBtoDsCrossSectionPP_ScaleBR_BptoDsmax;
	if(isPbPb){
		hBtoDs_ScaleBR_BptoDsmax=hBtoDsdNdPtPbPb_ScaleBR_BptoDsmax;
	}

	TH1D *h_PromptDs_ScaleBR_BptoDsmax=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_ScaleBR_BptoDsmax ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);
	double BptoDsmax_PNPRatio=PNPRatio_Glo;



	TH1D *hBtoDsCrossSectionPP_ScaleBR_BptoDsmin=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BptoDsmin");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_BptoDsmin=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BptoDsmin");
	TH1D *hBtoDs_ScaleBR_BptoDsmin=hBtoDsCrossSectionPP_ScaleBR_BptoDsmin;
	if(isPbPb){
		hBtoDs_ScaleBR_BptoDsmin=hBtoDsdNdPtPbPb_ScaleBR_BptoDsmin;
	}

	TH1D *h_PromptDs_ScaleBR_BptoDsmin=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_ScaleBR_BptoDsmin ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);
	double BptoDsmin_PNPRatio=PNPRatio_Glo;


	//


	TH1D *hBtoDsCrossSectionPP_ScaleBR_BstoDsmax=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDsmax");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_BstoDsmax=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDsmax");
	TH1D *hBtoDs_ScaleBR_BstoDsmax=hBtoDsCrossSectionPP_ScaleBR_BstoDsmax;
	if(isPbPb){
		hBtoDs_ScaleBR_BstoDsmax=hBtoDsdNdPtPbPb_ScaleBR_BstoDsmax;
	}

	TH1D *h_PromptDs_ScaleBR_BstoDsmax=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_ScaleBR_BstoDsmax ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);
	double BstoDsmax_PNPRatio=PNPRatio_Glo;



	TH1D *hBtoDsCrossSectionPP_ScaleBR_BstoDsmin=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleBR_BstoDsmin");
	TH1D *hBtoDsdNdPtPbPb_ScaleBR_BstoDsmin=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleBR_BstoDsmin");
	TH1D *hBtoDs_ScaleBR_BstoDsmin=hBtoDsCrossSectionPP_ScaleBR_BstoDsmin;
	if(isPbPb){
		hBtoDs_ScaleBR_BstoDsmin=hBtoDsdNdPtPbPb_ScaleBR_BstoDsmin;
	}

	TH1D *h_PromptDs_ScaleBR_BstoDsmin=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_ScaleBR_BstoDsmin ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);
	double BstoDsmin_PNPRatio=PNPRatio_Glo;



	// Scale Fr_Z & Fr_p

	TH1D *hBtoDsCrossSectionPP_ScaleFrZ=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleFrZ");
	TH1D *hBtoDsdNdPtPbPb_ScaleFrZ=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleFrZ");
	TH1D *hBtoDs_ScaleFrZ=hBtoDsCrossSectionPP_ScaleFrZ;
	if(isPbPb){
		hBtoDs_ScaleFrZ=hBtoDsdNdPtPbPb_ScaleFrZ;
	}

	TH1D *h_PromptDs_ScaleFrZ=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_ScaleFrZ ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);
	double BtoDsScaleFrZ_PNPRatio=PNPRatio_Glo;


	TH1D *hBtoDsCrossSectionPP_ScaleFrp=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight_ScaleFrp");
	TH1D *hBtoDsdNdPtPbPb_ScaleFrp=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight_ScaleFrp");
	TH1D *hBtoDs_ScaleFrp=hBtoDsCrossSectionPP_ScaleFrp;
	if(isPbPb){
		hBtoDs_ScaleFrp=hBtoDsdNdPtPbPb_ScaleFrp;
	}

	TH1D *h_PromptDs_ScaleFrp=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_ScaleFrp ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);
	double BtoDsScaleFrp_PNPRatio=PNPRatio_Glo;




	/////////////////////////
	/// for PFr syst error //
	/////////////////////////
	//// -- end BtoDs syst --////



	// do the nominal after BtoD Syst

	fixPNPRatio=0;
	h_CS_Default=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_AnaBin_pythiaWeight ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PNPRatio, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);
	RawYield_Default=h_RawYield_Default->GetBinContent(1);
	RawYieldErr_Default=h_RawYield_Default->GetBinError(1);

	// return 1;

	// h_EffMCP_Default=h_Eff_fun(isPbPb,f_MCP,DptLow,DptHigh,Dalpha_cut, Dchi2cl_cut, Ddls_cut, DtktkResmass_cut);

	// h_CS[index_best]=h_CSCal_fun(isPbPb,h_RawYield[index_best],hBtoDs_AnaBin_pythiaWeight ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut,Dchi2cl_cut,Ddls_cut,DtktkResmass_cut, fixPNPRatio, PNPRatio, PhiRatio);
	// fixPNPRatio=1;
	// h_CS[index_best]=h_CSCal_fun(isPbPb,h_RawYield[index_best],hBtoDs_AnaBin_pythiaWeight ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut,Dchi2cl_cut,Ddls_cut,DtktkResmass_cut, fixPNPRatio, PNPRatio_Glo, PhiRatio);

	cout<<"PNPRatio_Glo = "<<PNPRatio_Glo<<endl;

	fixPNPRatio=1;
	for(int i=0; i<nbin_var; i++){
		if(i==index_best && !useAnaBin)continue;

		if(var_scan=="Ddls"){
			Ddls_cut=bins_var[i];
		}
		else if(var_scan=="Dchi2cl"){
			Dchi2cl_cut=bins_var[i];
		}
		else if(var_scan=="Dalpha"){
			Dalpha_cut=bins_var[i];
		}
		else if(var_scan=="PhiMass"){
			DtktkResmass_cut=bins_var[i];
		}else{
			cout<<"unkonw var for scan : "<<var_scan<<endl;
			return 2;
		}


		h_CS[i]=h_CSCal_fun(isPbPb,h_RawYield[i],hBtoDs_AnaBin_pythiaWeight ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut,Dchi2cl_cut,Ddls_cut,DtktkResmass_cut, fixPNPRatio, PNPRatio_Glo, PhiRatio, Eff_MCP[i],EffErr_MCP[i],Eff_MCNP[i],EffErr_MCNP[i],PromptYield[i],PromptYieldErr[i]);
	}

	TGraphErrors *grMCP=new TGraphErrors();
	TGraphErrors *grMCNP=new TGraphErrors();
	TGraphErrors *grDataAll=new TGraphErrors();
	TGraphErrors *grDataPrompt=new TGraphErrors();

	TGraphErrors *gr_CS_Default=new TGraphErrors();
	gr_CS_Default->SetName("gr_CS_Default");


	TGraphErrors *gr_CS=new TGraphErrors();
	gr_CS->SetName("gr_CS");

	TGraphErrors *gr=new TGraphErrors();
	gr->SetName("gr_DoubleRatio");
	grMCP->SetName("grMCP");
	grMCNP->SetName("grMCNP");
	grDataAll->SetName("grDataAll");
	grDataPrompt->SetName("grDataPrompt");

	double BestCS=h_CS[index_best]->GetBinContent(1);
	double BestCSRelErr=h_CS[index_best]->GetBinError(1)/h_CS[index_best]->GetBinContent(1);

	double CS_Default=h_CS_Default->GetBinContent(1);
	double CSRelErr_Default=h_CS_Default->GetBinError(1)/h_CS_Default->GetBinContent(1);
	double CSErr_Default=h_CS_Default->GetBinError(1);


	gr_CS_Default->SetPoint(0,DefaultCutVal,CS_Default);
	gr_CS_Default->SetPointError(0,0,CSRelErr_Default*CS_Default);

	double Syst=0;
	double SystError=0;

	for(int i=0; i<nbin_var; i++){
		// if(i==index_best)continue;
		//    if(var_scan=="Ddls"){
		//      Ddls_cut=bins_var[i];
		//    }
		double CS=h_CS[i]->GetBinContent(1);
		double CSRelErr=h_CS[i]->GetBinError(1)/h_CS[i]->GetBinContent(1);
		double CSErr=h_CS[i]->GetBinError(1);
		// double CSRatio=CS/BestCS;
		// double CSRatioErr=sqrt(abs(CSRelErr*CSRelErr-BestCSRelErr*BestCSRelErr))*CSRatio;

		double CSRatio=CS/CS_Default;
		double CSRatioErr=sqrt(abs(CSErr*CSErr-CSErr_Default*CSErr_Default))/CS_Default;
		// double CSRatioErr=sqrt(abs(CSRelErr*CSRelErr-CSRelErr_Default*CSRelErr_Default))*CSRatio;
		// double CSRatioErr=CSRelErr*CSRatio;

		if(abs(CSRatio-1) > Syst ){
			Syst=abs(CSRatio-1);
			SystError=CSRelErr;	
		}

		gr_CS->SetPoint(i,bins_var[i],CS);
		gr_CS->SetPointError(i,0,CSRelErr*CS);

		gr->SetPoint(i,bins_var[i] ,CSRatio);
		gr->SetPointError(i,0.,CSRatioErr);		

		double EffRatio_MCP=Eff_MCP[i]/Eff_MCP_Default;
		double EffRatioErr_MCP=EffRatio_MCP * sqrt(abs(EffErr_MCP[i]/Eff_MCP[i]*EffErr_MCP[i]/Eff_MCP[i] - EffErr_MCP_Default/Eff_MCP_Default* EffErr_MCP_Default/Eff_MCP_Default ) );

		grMCP->SetPoint(i,bins_var[i], EffRatio_MCP);
		grMCP->SetPointError(i,0., EffRatioErr_MCP);

		double EffRatio_MCNP=Eff_MCNP[i]/Eff_MCNP_Default;
		double EffRatioErr_MCNP=EffRatio_MCNP * sqrt(abs(EffErr_MCNP[i]/Eff_MCNP[i]*EffErr_MCNP[i]/Eff_MCNP[i] - EffErr_MCNP_Default/Eff_MCNP_Default* EffErr_MCNP_Default/Eff_MCNP_Default ) );

		grMCNP->SetPoint(i,bins_var[i], EffRatio_MCNP);
		grMCNP->SetPointError(i,0., EffRatioErr_MCNP);

		RawYield=h_RawYield[i]->GetBinContent(1);
		RawYieldErr=h_RawYield[i]->GetBinError(1);

		double RawYield_Ratio=RawYield/RawYield_Default;
		double RawYield_RatioErr=RawYield_Ratio * sqrt( abs( RawYieldErr/RawYield*RawYieldErr/RawYield - RawYieldErr_Default/RawYield_Default* RawYieldErr_Default/RawYield_Default ) );

		grDataAll->SetPoint(i, bins_var[i], RawYield_Ratio);
		grDataAll->SetPointError(i, 0., RawYield_RatioErr);



		double PromptYield_Ratio=PromptYield[i]/PromptYield_Default;
		double PromptYield_RatioErr=PromptYield_Ratio * sqrt( abs( PromptYieldErr[i]/PromptYield[i]*PromptYieldErr[i]/PromptYield[i] - PromptYieldErr_Default/PromptYield_Default* PromptYieldErr_Default/PromptYield_Default ) );

		grDataPrompt->SetPoint(i, bins_var[i], PromptYield_Ratio);
		grDataPrompt->SetPointError(i, 0., PromptYield_RatioErr);


	}
	// gr->GetXaxis()->SetRangeUser(0,4.5);
	gStyle->SetOptStat(0);
	int nbinTemp=1;
	double binsTempLow=0;
	double binsTempHigh=5.5;
	if(isPbPb){binsTempHigh=5.5;}
	if(var_scan=="Dchi2cl"){
		binsTempHigh=0.3;
		if(isPbPb){
			binsTempHigh=0.5;
		}
	}

	if(var_scan=="PhiMass"){
		binsTempLow=bins_var[0] - 0.5*(bins_var[1]- bins_var[0]) ;
	}

	// test binHigh 
	binsTempHigh=bins_var[nbin_var-1] + 0.5*(bins_var[nbin_var-1]- bins_var[nbin_var-2]) ;


	TCanvas *c_SingleRatio= new TCanvas("c_SingleRatio","c_SingleRatio",800,800);
	c_SingleRatio->cd();

	TH1D *htemp2=new TH1D("htemp2",Form(";%s; Varied cut/Default cut",s_XTitle.Data()),nbinTemp,binsTempLow,binsTempHigh);
	htemp2->SetMaximum(2.5);
	htemp2->SetMinimum(0);
	htemp2->Draw();
	grMCP->SetLineColor(2);
	grMCP->SetMarkerColor(2);
	grMCP->Draw("p");
	grMCNP->SetLineColor(6);
	grMCNP->SetMarkerColor(6);
	grMCNP->Draw("p");
	grDataAll->SetLineColor(1);
	grDataAll->SetMarkerColor(1);
	grDataAll->Draw("p");
	grDataPrompt->SetLineColor(4);
	grDataPrompt->SetMarkerColor(4);
	grDataPrompt->Draw("p");

	TLegend *le_SingleRatio=new TLegend(0.5,0.55,0.8,0.8);
	le_SingleRatio->SetBorderSize(0);
	le_SingleRatio->AddEntry(grMCP,"MC Prompt D_{S}^{+}","p");
	le_SingleRatio->AddEntry(grMCNP,"MC NonPrompt D_{S}^{+}","p");
	le_SingleRatio->AddEntry(grDataAll,"Data D_{S}^{+} Yield","p");
	le_SingleRatio->AddEntry(grDataPrompt,"Data Prompt D_{S}^{+} Yield","p");
	le_SingleRatio->Draw("same");
	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s : %.0f < D_{S} p_{T} < %.0f",s_PbPb.Data(),DptLow,DptHigh));

	gSystem->Exec(Form("mkdir -p Ratioplots/%s%s/",var_scan.Data(),s_reWeightFname.Data()));
	c_SingleRatio->SaveAs(Form("./Ratioplots/%s%s/%s_%s_SingleRatio%s_Dpt%.0fto%.0f%s.png", var_scan.Data(),s_reWeightFname.Data() , s_PbPb.Data(),var_scan.Data(),s_fixPFr.Data(),DptLow,DptHigh,s_FixShape.Data()));

	TCanvas *c_out=new TCanvas("c_out","c_out",800,800);
	c_out->cd();

	TH1D *htemp=new TH1D("htemp",Form(";%s; #sigma_{Varied cut}/#sigma_{Default cut}",s_XTitle.Data()),nbinTemp,binsTempLow,binsTempHigh);
	htemp->SetBinContent(1,1);
	htemp->SetBinError(1,h_CS[index_best]->GetBinError(1)/h_CS[index_best]->GetBinContent(1));
	htemp->SetFillStyle(3001);
	htemp->SetFillColor(kGray);
	htemp->SetMarkerSize(0);
	htemp->SetMarkerColor(kGray);
	htemp->SetLineWidth(0);
	htemp->SetMaximum(2);
	htemp->SetMinimum(0);
	htemp->Draw("E2");
	gr->SetMarkerStyle(21);
	gr->SetMarkerSize(0.4);
	gr->Draw("p");

	TLine *tl=new TLine(DefaultCutVal,0,DefaultCutVal,2);
	tl->SetLineColor(4);
	tl->SetLineStyle(6);
	tl->Draw("same");

	// grMCP->Draw("p");
	// grMCNP->Draw("p");

	// TF1 *f1=new TF1("f1","[0]+x*[1]");
	TFitResultPtr fitR;
	double x0[1]={0};
	double x0err[1];
	TF1 *f1=new TF1("f1","[0]+(x-[2])*[1]");
	f1->SetLineColor(2);
	f1->SetParameter(0,1);
	// f1->FixParameter(0,1);
	f1->SetParameter(1,0);
	f1->SetParameter(2,DefaultCutVal);
	// f1->FixParameter(2,bins_var[index_best]);
	f1->SetRange(0,binsTempHigh);

	gr->Fit("f1","LEMS0F");
	gr->Fit("f1","LEMS0F");
	fitR=gr->Fit("f1","LEMIS0F");
	// fitR=gr->Fit("f1","EMIS0");
	fitR->GetConfidenceIntervals(1,1,1,x0,x0err,0.683,false);

	TLine *tl1=new TLine(binsTempLow,1,binsTempHigh,1);
	tl1->SetLineColor(2);
	if(var_scan=="Ddls" || var_scan=="Dchi2cl"){	
		f1->Draw("same");
		Syst=(abs(f1->Eval(0)-1));
		SystError=x0err[0];
	}else{
		tl1->Draw("same");
	}




	gStyle->SetOptFit(0);

	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s : %.0f < D_{S} p_{T} < %.0f",s_PbPb.Data(),DptLow,DptHigh)); shiftY-=oneshift;

	if(var_scan=="Ddls" || var_scan=="Dchi2cl"){
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Systematics : %.1f%% #pm %.1f%%",Syst*100, SystError*100  )); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Slope : %.3f #pm %.3f",f1->GetParameter(1), f1->GetParError(1)) ); shiftY-=oneshift;
	}else{
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Systematics : %.1f%% ",Syst*100 )); shiftY-=oneshift;
	}
	if(fixPFr){
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Fix Prompt Ratio =  %.2f",PfrVal )); shiftY-=oneshift;
	}

	gSystem->Exec(Form("mkdir -p Ratioplots/%s%s/",var_scan.Data(),s_reWeightFname.Data()));
	c_out->SaveAs(Form("./Ratioplots/%s%s/%s_%s%s_Dpt%.0fto%.0f%s.png",var_scan.Data(),s_reWeightFname.Data(),s_PbPb.Data(),var_scan.Data(),s_fixPFr.Data(),DptLow,DptHigh,s_FixShape.Data()));

	// save sys to histogram 

	TH1D *h_Syst=new TH1D("h_Syst","",1,0,1);
	h_Syst->SetBinContent(1,Syst );
	h_Syst->SetBinError(1,SystError);

	fout->cd();
	h_Syst->Write();
	gr_CS->Write();
	gr_CS_Default->Write();
	gr->Write();
	grMCP->Write();
	grMCNP->Write();
	grDataAll->Write();
	grDataPrompt->Write();
	htemp->Write("Stat_error");
	// calculate sys

	double *VarCut=DdlsMin_bins_pp;
	if(isPbPb){
		VarCut=DdlsMin_bins_PbPb3;
	}
	if(var_scan=="Dchi2cl"){
		VarCut=Dchi2clMin_bins_pp;
		if(isPbPb){
			VarCut=Dchi2clMin_bins_PbPb3;
		}
	}
	double f1_at0=f1->Eval(0);

	ofstream f_sys;
	gSystem->Exec("mkdir -p Syst"); 
	f_sys.open(Form("./Syst/%s_%s_Dpt%.0fto%.0f.txt",s_PbPb.Data(),var_scan.Data(),DptLow,DptHigh));

	for(int i =0; i<nbin_pt; i++){
		cout<<"Dpt = "<<bins_pt[i]<<" to "<<bins_pt[i+1]<<endl;
		if(bins_pt[i]>=DptLow && bins_pt[i+1]<=DptHigh){
			// if(VarCut[i]<0.02){VarCut[i]=0.02;}
			cout<<"var cut = "<<VarCut[i]<<" f1val = "<<f1->Eval(VarCut[i])<<endl;
			cout<<"Syst. = "<<abs(f1_at0/f1->Eval(VarCut[i]) -1)*100<<" %"<<endl;
			f_sys<<"Dpt = "<<bins_pt[i]<<" to "<<bins_pt[i+1]<<endl;
			f_sys<<"Syst. = "<<abs(f1_at0/f1->Eval(VarCut[i]) -1)*100<<" %"<<endl;
		}

	}



	// TH1D *h_PFrDefault=new TH1D("h_PFrDefault","h_PFrDefault",1,0,1);
	// h_PFrDefault->SetBinContent(1,PNPRatio_Glo);
	// h_PFrDefault->SetBinError(1,PNPRatioErr_Glo);

	cout<<"PNPRatio_Glo = "<<PNPRatio_Glo<<" , PNPRatioUp_Glo = "<<PNPRatioUp_Glo<<" , PNPRatioDown_Glo = "<<PNPRatioDown_Glo<<endl;

	TGraphAsymmErrors *gr_PFrDefault=new TGraphAsymmErrors();
	gr_PFrDefault->SetName("gr_PFrDefault");
	gr_PFrDefault->SetPoint(0,(DptLow+DptHigh)/2,PNPRatio_Glo);
	gr_PFrDefault->SetPointError(0,(DptLow+DptHigh)/2-DptLow, DptHigh-(DptLow+DptHigh)/2,PNPRatio_Glo-PNPRatioDown_Glo, PNPRatioUp_Glo-PNPRatio_Glo);


	cout<<"gr PFr up= "<<gr_PFrDefault->GetErrorYhigh(0)<<" gr PFr down = "<<gr_PFrDefault->GetErrorYlow(0)<<endl;
	double pointx=0;
	double pointy=0;
	gr_PFrDefault->GetPoint(0,pointx,pointy);

	cout<<"gr PFr = "<<pointy<<" , Ds pt center= "<<pointx<<endl;


	// gr_PFrDefault->SetPointError(1,DptLow,DptHigh,PNPRatioUp_Glo-PNPRatioUp_Glo,PNPRatioDown_Glo-PNPRatio_Glo);

	// double CS_PromptDs=h_CS_Default->GetBinContent(1);
	// double CSErr_PromptDs=h_CS_Default->GetBinError(1);

	// TH1D *hBtoDsCSdNdpt_Integral_temp=(TH1D*)hBtoDs_AnaBin_pythiaWeight->Clone("hBtoDsCSdNdpt_Integral_temp");
	// MutiplyBinWidth(hBtoDsCSdNdpt_Integral_temp);
	// double CSErr_NonPromptDs=0;
	// double CS_NonPromptDs=hBtoDsCSdNdpt_Integral_temp->Integral(hBtoDsCSdNdpt_Integral_temp->FindBin(DptLow+0.001), hBtoDsCSdNdpt_Integral_temp->FindBin(DptHigh-0.001))/(DptHigh-DptLow);

	// double CS_NonPromptDs=hBtoDs_AnaBin_pythiaWeight->GetBinContent(ibin_Dpt+1);
	// double CSErr_NonPromptDs=hBtoDs_AnaBin_pythiaWeight->GetBinError(ibin_Dpt+1);;
	fout->cd();
	// h_PFrDefault->Write();	
	gr_PFrDefault->Write();
	gr_PFrDefault->Draw("ap");


	double xcenter=(DptLow+DptHigh)/2;
	double xlow=(DptLow+DptHigh)/2-DptLow;
	double xhigh=DptHigh-(DptLow+DptHigh)/2;
	// calculate Fr syst
	double FONLL_PNPRatio_RelErr=abs(FONLL_PNPRatio - PNPRatio_Glo)/PNPRatio_Glo;

	double BtoDErrup_PNPRatio_RelErr=abs(BtoDErrup_PNPRatio - PNPRatio_Glo)/PNPRatio_Glo; 
	double BtoDErrdown_PNPRatio_RelErr=abs(BtoDErrdown_PNPRatio - PNPRatio_Glo)/PNPRatio_Glo; 
	double BtoDErr_PNPRatio_RelErr= BtoDErrup_PNPRatio_RelErr >BtoDErrdown_PNPRatio_RelErr ? BtoDErrup_PNPRatio_RelErr : BtoDErrdown_PNPRatio_RelErr;

	double BtoDsScaleFrZ_PNPRatio_RelErr=abs(BtoDsScaleFrZ_PNPRatio - PNPRatio_Glo)/PNPRatio_Glo;
	double BtoDsScaleFrp_PNPRatio_RelErr=abs(BtoDsScaleFrp_PNPRatio - PNPRatio_Glo)/PNPRatio_Glo;
	double BtoDsScaleFr_PNPRatio_RelErr=BtoDsScaleFrZ_PNPRatio_RelErr > BtoDsScaleFrp_PNPRatio_RelErr ? BtoDsScaleFrZ_PNPRatio_RelErr : BtoDsScaleFrp_PNPRatio_RelErr;


	double B0toDmax_PNPRatio_RelErr=abs(B0toDmax_PNPRatio - PNPRatio_Glo)/PNPRatio_Glo;
	double B0toDmin_PNPRatio_RelErr=abs(B0toDmin_PNPRatio - PNPRatio_Glo)/PNPRatio_Glo;
	double B0toD_PNPRatio_RelErr=B0toDmax_PNPRatio_RelErr > B0toDmin_PNPRatio_RelErr ? B0toDmax_PNPRatio_RelErr : B0toDmin_PNPRatio_RelErr;
	double BptoDmax_PNPRatio_RelErr=abs(BptoDmax_PNPRatio - PNPRatio_Glo)/PNPRatio_Glo;
	double BptoDmin_PNPRatio_RelErr=abs(BptoDmin_PNPRatio - PNPRatio_Glo)/PNPRatio_Glo;
	double BptoD_PNPRatio_RelErr=BptoDmax_PNPRatio_RelErr > BptoDmin_PNPRatio_RelErr ? BptoDmax_PNPRatio_RelErr : BptoDmin_PNPRatio_RelErr;
	double BstoDmax_PNPRatio_RelErr=abs(BstoDmax_PNPRatio - PNPRatio_Glo)/PNPRatio_Glo;
	double BstoDmin_PNPRatio_RelErr=abs(BstoDmin_PNPRatio - PNPRatio_Glo)/PNPRatio_Glo;
	double BstoD_PNPRatio_RelErr=BstoDmax_PNPRatio_RelErr > BstoDmin_PNPRatio_RelErr ? BstoDmax_PNPRatio_RelErr : BstoDmin_PNPRatio_RelErr;

	double B0toDsmax_PNPRatio_RelErr=abs(B0toDsmax_PNPRatio - PNPRatio_Glo)/PNPRatio_Glo;
	double B0toDsmin_PNPRatio_RelErr=abs(B0toDsmin_PNPRatio - PNPRatio_Glo)/PNPRatio_Glo;
	double B0toDs_PNPRatio_RelErr=B0toDsmax_PNPRatio_RelErr > B0toDsmin_PNPRatio_RelErr ? B0toDsmax_PNPRatio_RelErr : B0toDsmin_PNPRatio_RelErr;
	double BptoDsmax_PNPRatio_RelErr=abs(BptoDsmax_PNPRatio - PNPRatio_Glo)/PNPRatio_Glo;
	double BptoDsmin_PNPRatio_RelErr=abs(BptoDsmin_PNPRatio - PNPRatio_Glo)/PNPRatio_Glo;
	double BptoDs_PNPRatio_RelErr=BptoDsmax_PNPRatio_RelErr > BptoDsmin_PNPRatio_RelErr ? BptoDsmax_PNPRatio_RelErr : BptoDsmin_PNPRatio_RelErr;
	double BstoDsmax_PNPRatio_RelErr=abs(BstoDsmax_PNPRatio - PNPRatio_Glo)/PNPRatio_Glo;
	double BstoDsmin_PNPRatio_RelErr=abs(BstoDsmin_PNPRatio - PNPRatio_Glo)/PNPRatio_Glo;
	double BstoDs_PNPRatio_RelErr=BstoDsmax_PNPRatio_RelErr > BstoDsmin_PNPRatio_RelErr ? BstoDsmax_PNPRatio_RelErr : BstoDsmin_PNPRatio_RelErr;

	double DtoDs_sysTotal_RelErr=sqrt(B0toD_PNPRatio_RelErr*B0toD_PNPRatio_RelErr + BptoD_PNPRatio_RelErr*BptoD_PNPRatio_RelErr + BstoD_PNPRatio_RelErr*BstoD_PNPRatio_RelErr + B0toDs_PNPRatio_RelErr*B0toDs_PNPRatio_RelErr + BptoDs_PNPRatio_RelErr*BptoDs_PNPRatio_RelErr + BstoDs_PNPRatio_RelErr*BstoDs_PNPRatio_RelErr + BtoDsScaleFr_PNPRatio_RelErr*BtoDsScaleFr_PNPRatio_RelErr + FONLL_PNPRatio_RelErr*FONLL_PNPRatio_RelErr + BtoDErr_PNPRatio_RelErr*BtoDErr_PNPRatio_RelErr);



	TGraphAsymmErrors *gr_PFrDefault_Sys=new TGraphAsymmErrors();
	gr_PFrDefault_Sys->SetName("gr_PFrDefault_Sys");
	gr_PFrDefault_Sys->SetPoint(0,(DptLow+DptHigh)/2,PNPRatio_Glo);

	double yerrUp=PNPRatio_Glo*DtoDs_sysTotal_RelErr;
	if((yerrUp+PNPRatio_Glo) > 1){
		cout<<"up err exceed 1 , change to 1"<<endl;
		yerrUp=1-PNPRatio_Glo;		
	}

	gr_PFrDefault_Sys->SetPointError(0,(DptLow+DptHigh)/2-DptLow, DptHigh-(DptLow+DptHigh)/2 ,  DtoDs_sysTotal_RelErr*PNPRatio_Glo, yerrUp);

	gr_PFrDefault_Sys->SetLineColor(2);
	gr_PFrDefault_Sys->SetFillColor(kRed-4);
	gr_PFrDefault_Sys->SetFillStyle(3005);
	gr_PFrDefault_Sys->Draw("apF2");
	gr_PFrDefault_Sys->Write();

	gr_PFrDefault->Draw("p");
	// sys

	cout<<"DtoDs_sysTotal_RelErr = "<<DtoDs_sysTotal_RelErr<<endl;

	// return 1;

	// doPFrScan

	if(doPFrScan && !fixPFr){
		int NpfScanMax=100;
		TH1D *h_CS_Default_PFrScan[NpfScanMax];
		TH1D *h_CS_PFrScan[NpfScanMax][nbin_var];

		TGraphErrors *gr_PFrScan[NpfScanMax];
		TGraphErrors *gr_PFrScan_Sys=new TGraphErrors();
		TGraphErrors *gr_PFrScan_Sys_Def=new TGraphErrors();
		gr_PFrScan_Sys_Def->SetLineColor(2);
		gr_PFrScan_Sys_Def->SetMarkerColor(2);
		gr_PFrScan_Sys_Def->SetMarkerStyle(23);

		gr_PFrScan_Sys_Def->SetPoint(1,PNPRatio_Glo,Syst);
		gr_PFrScan_Sys_Def->SetPointError(1,0,SystError);

		fixPNPRatio=1;
		double PfrValLow=0.7;
		double PfrValHigh=1.0;
		double PfrValstep=0.01;
		PfrVal=PfrValLow;
		int j=0;

		double PFr_Best=0;
		double Sys_Best=100;


		while(PfrVal<=PfrValHigh+0.00001 && j<=NpfScanMax){

			gr_PFrScan[j]=new TGraphErrors();


			h_CS_Default_PFrScan[j]=h_CSCal_fun(isPbPb,h_RawYield_Default,hBtoDs_AnaBin_pythiaWeight ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut_Def,Dchi2cl_cut_Def,Ddls_cut_Def,DtktkResmass_cut_Def, fixPNPRatio, PfrVal, PhiRatio, Eff_MCP_Default, EffErr_MCP_Default , Eff_MCNP_Default, EffErr_MCNP_Default, PromptYield_Default,PromptYieldErr_Default);

			cout<<"PfrVal = "<<PfrVal<<endl;
			cout<<"h_CS_Default = "<<h_CS_Default_PFrScan[j]->GetBinContent(1);

			for(int i=0; i<nbin_var; i++){
				if(i==index_best && !useAnaBin)continue;

				if(var_scan=="Ddls"){
					Ddls_cut=bins_var[i];
				}
				else if(var_scan=="Dchi2cl"){
					Dchi2cl_cut=bins_var[i];
				}
				else if(var_scan=="Dalpha"){
					Dalpha_cut=bins_var[i];
				}
				else if(var_scan=="PhiMass"){
					DtktkResmass_cut=bins_var[i];
				}else{
					cout<<"unkonw var for scan : "<<var_scan<<endl;
					return 2;
				}

				h_CS_PFrScan[j][i]=h_CSCal_fun(isPbPb,h_RawYield[i],hBtoDs_AnaBin_pythiaWeight ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut,Dchi2cl_cut,Ddls_cut,DtktkResmass_cut, fixPNPRatio, PfrVal, PhiRatio, Eff_MCP[i],EffErr_MCP[i],Eff_MCNP[i],EffErr_MCNP[i],PromptYield[i],PromptYieldErr[i]);
			}


			Syst=0;
			SystError=0;

			CS_Default=h_CS_Default_PFrScan[j]->GetBinContent(1);
			CSRelErr_Default=h_CS_Default_PFrScan[j]->GetBinError(1)/h_CS_Default_PFrScan[j]->GetBinContent(1);
			CSErr_Default=h_CS_Default_PFrScan[j]->GetBinError(1);

			for(int i=0; i<nbin_var; i++){
				double CS=h_CS_PFrScan[j][i]->GetBinContent(1);
				double CSRelErr=h_CS_PFrScan[j][i]->GetBinError(1)/h_CS_PFrScan[j][i]->GetBinContent(1);
				double CSErr=h_CS_PFrScan[j][i]->GetBinError(1);
				// double CSRatio=CS/BestCS;
				// double CSRatioErr=sqrt(abs(CSRelErr*CSRelErr-BestCSRelErr*BestCSRelErr))*CSRatio;

				double CSRatio=CS/CS_Default;
				double CSRatioErr=sqrt(abs(CSErr*CSErr-CSErr_Default*CSErr_Default))/CS_Default;
				// double CSRatioErr=sqrt(abs(CSRelErr*CSRelErr-CSRelErr_Default*CSRelErr_Default))*CSRatio;
				// double CSRatioErr=CSRelErr*CSRatio;

				// for Dalpha && Phimass cut scan
				if(abs(CSRatio-1) > Syst ){
					Syst=abs(CSRatio-1);
					SystError=CSRelErr;	
				}

				gr_PFrScan[j]->SetPoint(i,bins_var[i] ,CSRatio);
				gr_PFrScan[j]->SetPointError(i,0.,CSRatioErr);		

			} // end for i<nbin_var;

			TCanvas *c_out2=new TCanvas("c_out2","c_out2",800,800);
			c_out2->cd();

			TH1D *htemp3=new TH1D("htemp3",Form(";%s; #sigma_{Varied cut}/#sigma_{Default cut}",s_XTitle.Data()),nbinTemp,binsTempLow,binsTempHigh);
			htemp3->SetBinContent(1,1);
			htemp3->SetBinError(1,h_CS[index_best]->GetBinError(1)/h_CS[index_best]->GetBinContent(1));
			htemp3->SetFillStyle(3001);
			htemp3->SetFillColor(kGray);
			htemp3->SetMarkerSize(0);
			htemp3->SetMarkerColor(kGray);
			htemp3->SetLineWidth(0);
			htemp3->SetMaximum(2);
			htemp3->SetMinimum(0);
			htemp3->Draw("E2");
			gr_PFrScan[j]->SetMarkerStyle(21);
			gr_PFrScan[j]->SetMarkerSize(0.4);
			gr_PFrScan[j]->Draw("p");

			TLine *tl2=new TLine(DefaultCutVal,0,DefaultCutVal,2);
			tl2->SetLineColor(4);
			tl2->SetLineStyle(6);
			tl2->Draw("same");

			// grMCP->Draw("p");
			// grMCNP->Draw("p");

			// TF1 *f1=new TF1("f1","[0]+x*[1]");
			// TFitResultPtr fitR;

			f1->SetParameter(2,DefaultCutVal);
			f1->SetRange(0,binsTempHigh);

			gr_PFrScan[j]->Fit("f1","LEMS0F");
			gr_PFrScan[j]->Fit("f1","LEMS0F");
			fitR=gr_PFrScan[j]->Fit("f1","LEMIS0F");
			// fitR=gr->Fit("f1","EMIS0");
			fitR->GetConfidenceIntervals(1,1,1,x0,x0err,0.683,false);

			// TLine *tl1=new TLine(binsTempLow,1,binsTempHigh,1);
			// tl1->SetLineColor(2);
			if(var_scan=="Ddls" || var_scan=="Dchi2cl"){	
				f1->Draw("same");
				Syst=(abs(f1->Eval(0)-1));
				SystError=x0err[0];
			}else{
				tl1->Draw("same");
			}

			if(Syst<Sys_Best){
				Sys_Best=Syst;
				PFr_Best=PfrVal;
			}

			shiftY=0.06;
			tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s : %.0f < D_{S} p_{T} < %.0f",s_PbPb.Data(),DptLow,DptHigh)); shiftY-=oneshift;

			if(var_scan=="Ddls" || var_scan=="Dchi2cl"){
				tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Systematics : %.1f%% #pm %.1f%%",Syst*100, SystError*100  )); shiftY-=oneshift;
				tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Slope : %.3f #pm %.3f",f1->GetParameter(1), f1->GetParError(1)) ); shiftY-=oneshift;
			}else{
				tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Systematics : %.1f%% ",Syst*100 )); shiftY-=oneshift;
			}

			gSystem->Exec(Form("mkdir -p RatioScanPlots/%s%s/",var_scan.Data(),s_reWeightFname.Data()));
			c_out2->SaveAs(Form("./RatioScanPlots/%s%s/%s_%s_Dpt%.0fto%.0f%s_PFr%.0f.png",var_scan.Data(),s_reWeightFname.Data(),s_PbPb.Data(),var_scan.Data(),DptLow,DptHigh,s_FixShape.Data(), PfrVal*100));

			gr_PFrScan_Sys->SetPoint(j,PfrVal,Syst);
			gr_PFrScan_Sys->SetPointError(j,0,SystError);



			delete tl2;
			delete c_out2;


			PfrVal+=PfrValstep;
			j++;
		} // end PfrVal<=PfrValHigh

		TCanvas *c_PFrScan=new TCanvas("c_PFrScan","c_PFrScan",800,800);
		c_PFrScan->cd();
		gr_PFrScan_Sys->GetXaxis()->SetTitle("Prompt Fraction");
		gr_PFrScan_Sys->GetYaxis()->SetTitle("Systematics");
		gr_PFrScan_Sys->SetMaximum(0.8);
		gr_PFrScan_Sys->SetMinimum(-0.2);	
		gr_PFrScan_Sys->Draw("ap");
		gr_PFrScan_Sys_Def->Draw("p");


		shiftX=0.1;
		shiftY=0.06;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s : %.0f < D_{S} p_{T} < %.0f , %s",s_PbPb.Data(),DptLow,DptHigh, s_XTitle.Data() )); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Default Prompt fraction : %.3f",PNPRatio_Glo)); shiftY-=oneshift;
		tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("Best Prompt fraction : %.3f",PFr_Best)); shiftY-=oneshift;

		c_PFrScan->SaveAs(Form("./RatioScanPlots/%s%s/%s_%s_Dpt%.0fto%.0f%s_PFrScanSyst.png",var_scan.Data(),s_reWeightFname.Data(),s_PbPb.Data(),var_scan.Data(),DptLow,DptHigh,s_FixShape.Data()));

		fout->cd();
		gr_PFrScan_Sys->SetName("gr_PFrScan_Sys");
		gr_PFrScan_Sys->Write();

		TH1D *h_PFrBest=new TH1D("h_PFrBest","h_PFrBest",1,0,1);
		h_PFrBest->SetBinContent(1,PFr_Best);

		h_PFrBest->Write();

		TGraphErrors *gr_PFrBest=new TGraphErrors();
		gr_PFrBest->SetName("gr_PFrBest");
		gr_PFrBest->SetPoint(0,(DptLow+DptHigh)/2,PFr_Best );
		gr_PFrBest->SetPointError(0,(DptLow+DptHigh)/2-DptLow,0);
		gr_PFrBest->Write();



	} // end if doPFrScan



	/*
		 TTree *t_Ds_MCPrompt=(TTree*)f_mc_Prompt->Get(Form("t_fit"));
		 TTree *t_Ds_MCNonPrompt=(TTree*)f_mc_NonPrompt->Get(Form("t_fit"));

		 double N_yield=h_RawYield[index_best]->GetBinContent(1);
		 double LumiNevt=LumiSum;
		 if(isPbPb){LumiNevt=NevtPbPb3;}

		 TFile *f_NonPromptDs=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BtoDs_Results/output/BtoDs.root");
		 TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight");
		 TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight");
		 TH1D *hBtoDs_AnaBin_pythiaWeight=hBtoDsCrossSectionPP_AnaBin_pythiaWeight;  // this is differential cross section
		 if(isPbPb){ hBtoDs_AnaBin_pythiaWeight=hBtoDsdNdPtPbPb_AnaBin_pythiaWeight; }
		 MutiplyBinWidth(hBtoDs_AnaBin_pythiaWeight);
		 double CS_Integral=hBtoDs_AnaBin_pythiaWeight->Integral(hBtoDs_AnaBin_pythiaWeight->GetXaxis()->FindBin(Dpt_Low+0.00001), hBtoDs_AnaBin_pythiaWeight->GetXaxis()->FindBin(Dpt_High-0.00001) );
		 cout<<"CS_Integral = "<<CS_Integral<<endl;
		 */


	// int CalCutScan_AnaBin_temp(int isPbPb=0, int useAnaBin=7, int ibin_Dpt=7, double DptLow=10, double DptHigh=20, TString var_scan="Ddls", int FixShape=1, int useReweight=0, TString ReweightTree="DdxyzErrWeight", int doPFrScan=0, int fixPFr=0, double PfrVal=0.9){
	cout<<"isPbPb = "<<isPbPb<<endl;
	cout<<"useAnaBin = "<<useAnaBin<<endl;
	cout<<"ibin_Dpt = "<<ibin_Dpt<<endl;
	cout<<"DptLow = "<<DptLow<<endl;
	cout<<"DptHigh = "<<DptHigh<<endl;
	cout<<"var_scan = "<<var_scan<<endl;
	cout<<"FixShape = "<<FixShape<<endl;
	cout<<"useReweight = "<<useReweight<<endl;
	cout<<"ReweightTree = "<<ReweightTree<<endl;
	cout<<"doPFrScan = "<<doPFrScan<<endl;
	cout<<"fixPFr = "<<fixPFr<<endl;
	cout<<"PfrVal = "<<PfrVal<<endl;


	return 0;

}


int main(int argc , char*argv[]){

	// int CalCutScan_AnaBin_temp(int isPbPb=0, int useAnaBin=7, int ibin_Dpt=7, double DptLow=10, double DptHigh=20, TString var_scan="Ddls", int FixShape=1, int useReweight=0, TString ReweightTree="DdxyzErrWeight", int doPFrScan=0, int fixPFr=0, double PfrVal=0.9){
	// int CalCutScan_AnaBin_temp(int isPbPb=0, int useAnaBin=1, int ibin_Dpt=7, double DptLow=10, double DptHigh=20, TString var_scan="Dalpha", int FixShape=1, int useReweight=1, TString ReweightTree="DdxyzErrWeight", int doPFrScan=1, int fixPFr=0, double PfrVal=0.7){
// int CalCutScan_AnaBin_Smear(int isPbPb=0, int useAnaBin=1, int ibin_Dpt=7, double DptLow=20, double DptHigh=40, TString var_scan="Ddls", int FixShape=1, int useReweight=0, TString ReweightTree="DdxyzErrWeight", int doPFrScan=0, int fixPFr=0, double PfrVal=0.7, int doSmear=1, int nSmear=10){
	if(argc==13){

		CalCutScan_AnaBin_Smear( atoi(argv[1]) , atoi(argv[2]) , atoi(argv[3]), atof(argv[4]) , atof(argv[5]) , argv[6] , atoi(argv[7]) , atoi(argv[8]) , argv[9] , atoi(argv[10]) , atoi(argv[11]), atof(argv[12]) );
	}else{
		cout<<"wrong number of input"<<endl;
		return 1;
	}
	return 0;

}

int smearTree(TFile *f_mcp_new,TFile *f_mc, int nSmear=10, double DptLow=2,double DptHigh=40){

	TRandom3 Rdm;
	Rdm.SetSeed(0);

	TTree *t_mcp=(TTree*)f_mc->Get("t_fit");
	TH1D *hGen_pt=(TH1D*)f_mc->Get("hGen_pt");
	f_mcp_new->cd();
	hGen_pt->Write("",TObject::kOverwrite);

	Float_t Dmass;
	Float_t Dpt;
	Float_t Ddls;
	Float_t Ddl;
	Float_t DdlErr;
	Float_t Dalpha;
	Float_t Dchi2cl;
	Float_t DtktkResmass;
	Float_t TotalWeight;  
	Float_t DlxyBS;
	Float_t DlxyBSErr;
	Float_t DlxyBSs;
	Float_t Dalpha_BS_2D;

	t_mcp->SetBranchAddress("Dmass",&Dmass);
	t_mcp->SetBranchAddress("Ddls",&Ddls);
	t_mcp->SetBranchAddress("Dpt",&Dpt);
	t_mcp->SetBranchAddress("DdlErr",&DdlErr);
	t_mcp->SetBranchAddress("Ddl",&Ddl);
	t_mcp->SetBranchAddress("DtktkResmass",&DtktkResmass);
	t_mcp->SetBranchAddress("Dalpha",&Dalpha);
	t_mcp->SetBranchAddress("Dchi2cl",&Dchi2cl);
	t_mcp->SetBranchAddress("TotalWeight",&TotalWeight);

	t_mcp->SetBranchAddress("Dalpha_BS_2D",&Dalpha_BS_2D);
	t_mcp->SetBranchAddress("DlxyBS",&DlxyBS);
	t_mcp->SetBranchAddress("DlxyBSErr",&DlxyBSErr);
	t_mcp->SetBranchAddress("DlxyBSs",&DlxyBSs);

	Float_t SMrWt=(float)1.0/nSmear;

	// Float_t SmearRelFactorArr[]={0.01,0.05,0.1,0.2,0.3};
	// Float_t SmearAbsFactorArr[]={0.0001,0.001,0.02,0.05,0.1};
	// Float_t SmearFactorArr[]={0.01,0.05};
	// const int nSmrRelF=sizeof(SmearRelFactorArr)/sizeof(SmearRelFactorArr[0]);
	// const int nSmrAbsF=sizeof(SmearAbsFactorArr)/sizeof(SmearAbsFactorArr[0]);
	// Float_t Ddls_SmrRelF[nSmrRelF];
	// Float_t Ddls_SmrAbsF[nSmrAbsF];

	f_mcp_new->cd();
	TTree *t_new=new TTree("t_fit","t_fit");
	t_new->Branch("Dmass",&Dmass);
	t_new->Branch("Dpt",&Dpt);
	t_new->Branch("Ddls",&Ddls);
	t_new->Branch("Ddl",&Ddl);
	t_new->Branch("DdlErr",&DdlErr);
	t_new->Branch("Dalpha",&Dalpha);
	t_new->Branch("Dchi2cl",&Dchi2cl);
	t_new->Branch("DtktkResmass",&DtktkResmass);
	t_new->Branch("TotalWeight",&TotalWeight);
	t_new->Branch("SMrWt",&SMrWt);
	for(int i=0; i<nSmrRelF; i++){
		t_new->Branch(Form("Ddls_Rel%.0fem5",SmearRelFactorArr[i]*1e5),&Ddls_SmrRelF[i]);
		t_new->Branch(Form("DdlErr_Rel%.0fem5",SmearRelFactorArr[i]*1e5),&DdlErr_SmrRelF[i]);
	}
	for(int i=0; i<nSmrAbsF; i++){
		t_new->Branch(Form("Ddls_Abs%.0fem5",SmearAbsFactorArr[i]*1e5),&Ddls_SmrAbsF[i]);
		t_new->Branch(Form("DdlErr_Abs%.0fem5",SmearAbsFactorArr[i]*1e5),&DdlErr_SmrAbsF[i]);
	}
	for(int i=0; i<nSclErrF; i++){
		t_new->Branch(Form("Ddls_Scl%.0fem3",ScaleErrFactorArr[i]*1e3),&Ddls_SclErrF[i]);
		t_new->Branch(Form("DdlErr_Scl%.0fem5",ScaleErrFactorArr[i]*1e3),&DdlErr_SclErrF[i]);
	}


	Long64_t nentries=t_mcp->GetEntries();
	double smF=0.1;
	double DdlErrTemp=0;

	for(Long64_t i=0;i<nentries; i++)
		// for(Long64_t i=0;i<500; i++)
	{   
		if(i%200000==0) {cout<<setw(10)<<i<<" / "<<nentries<<endl;}
		t_mcp->GetEntry(i);

		if(Dpt<DptLow || Dpt>DptHigh) {
			continue;
		}


		for(int j=0; j<nSmear;j++){



			for(int k=0;k<nSmrRelF;k++){
				DdlErr_SmrRelF[k]=Rdm.Gaus(DdlErr,DdlErr*SmearRelFactorArr[k]);
				while(DdlErr_SmrRelF[k]<0){
	        DdlErr_SmrRelF[k]=Rdm.Gaus(DdlErr,DdlErr*SmearRelFactorArr[k]);					
				}
				Ddls_SmrRelF[k]=Ddl / DdlErr_SmrRelF[k];
				// Ddls_SmrRelF[k]=Ddl / Rdm.Gaus(DdlErr,DdlErr*SmearRelFactorArr[k]);
				// while(Ddls_SmrRelF[k]<0){
					// Ddls_SmrRelF[k]=Ddl / Rdm.Gaus(DdlErr,DdlErr*SmearRelFactorArr[k]);
				// }       
			}
			for(int k=0;k<nSmrAbsF;k++){
				DdlErr_SmrAbsF[k]=Rdm.Gaus(DdlErr,SmearAbsFactorArr[k]);
				while(DdlErr_SmrAbsF[k]<0){
	        DdlErr_SmrAbsF[k]=Rdm.Gaus(DdlErr,SmearAbsFactorArr[k]);				
				}
				Ddls_SmrAbsF[k]=Ddl / DdlErr_SmrAbsF[k];
				// Ddls_SmrAbsF[k]=Ddl / Rdm.Gaus(DdlErr,DdlErr*SmearAbsFactorArr[k]);
				// while(Ddls_SmrAbsF[k]<0){
					// Ddls_SmrAbsF[k]=Ddl / Rdm.Gaus(DdlErr,DdlErr*SmearAbsFactorArr[k]);
				// }       
			}
			for(int k=0;k<nSclErrF;k++){
				DdlErr_SclErrF[k]=DdlErr*ScaleErrFactorArr[k];
				Ddls_SclErrF[k]=Ddl / DdlErr_SclErrF[k];
				if(ScaleErrFactorArr[k]>1){
					Ddls_SclErrF[k]=Rdm.Gaus(Ddls_SclErrF[k],sqrt(1.0*1.0-1.0/ScaleErrFactorArr[k]*1.0/ScaleErrFactorArr[k]*1.0));
				}
				// Ddls_SmrAbsF[k]=Ddl / Rdm.Gaus(DdlErr,DdlErr*SmearAbsFactorArr[k]);
				// while(Ddls_SmrAbsF[k]<0){
					// Ddls_SmrAbsF[k]=Ddl / Rdm.Gaus(DdlErr,DdlErr*SmearAbsFactorArr[k]);
				// }       
			}

			t_new->Fill();
		}

	} 

	f_mcp_new->cd();
	t_new->Write("",TObject::kOverwrite);

	return 1; 


}

