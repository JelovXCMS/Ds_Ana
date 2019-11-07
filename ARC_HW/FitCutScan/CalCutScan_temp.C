// double bins_var[]={3.5,4.0,4.5,5.0};
// const int nbin_var=4;

// double bins_var[]={1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0};
// const int nbin_var=11;
double bins_var[]={1.5,2.5,3.5,4.5,5.5};
const int nbin_var=5;

// double bins_var[]={3.25 ,3.5,3.75,4.0,4.25,4.5,4.75,5,5.25};
// const int nbin_var=9;

double PNPRatio_Glo=0.85;

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
#include "varCompare_para.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"

  double textposx=0.2;
  double textposy=0.77;

  double shiftY=0.0;
  double shiftX=0.3;
  double oneshift=0.075;

	int reWeight=0;
	TString s_reWeight="DdxyzErrWeight";


TH1D *h_CSCal_fun(int isPbPb,TH1D *hRawYield, TH1D* hBtoDsCSdNdpt_temp,TFile *f_mc_Prompt,TFile *f_mc_NonPrompt,double DptLow,double DptHigh,double Dalpha_cut,double Dchi2cl_cut,double Ddls_cut,double DtktkResmass_cut,int fixPNPRatio,double PNPRatio,double PhiRatio){
	// hRawYield[i] is integral in pt ranges,
	// hBtoDsCSdNdpt_temp is 1/dpt,
	// calculate integral first, then get 1/dpt in the end

	double LumiNevt=LumiSum;
	if(isPbPb){
		LumiNevt=NevtPbPb3;
	}

	// get Efficiency

//  TString DataCuts=Form("Dpt>%f && Dpt<%f && Dalpha < %f && Dchi2cl > %f && Ddls >%f && %s > %f && %s <%f", Dpt_Low, Dpt_High, Dalpha_cut, Dchi2cl_cut, Ddls_cut, var_cut.Data(), var_cutLow,var_cut.Data(), var_cutHigh); 
 TString DataCuts=Form("Dpt>%f && Dpt<%f && Dalpha<%f && Dchi2cl>%f && Ddls>%f && DtktkResmass>%f && DtktkResmass<%f", DptLow,DptHigh,Dalpha_cut,Dchi2cl_cut,Ddls_cut, DtktkResmassCutMean-DtktkResmass_cut, DtktkResmassCutMean+DtktkResmass_cut);
 TString S_reWeight="1";
	if(reWeight==1){
	S_reWeight=s_reWeight;
	}

	TTree *t_MCP=(TTree*)f_mc_Prompt->Get("t_fit");
	TTree *t_MCNP=(TTree*)f_mc_NonPrompt->Get("t_fit");

	TH1D *h_reco_cut_MCP=new TH1D("h_reco_cut_MCP","h_reco_cut_MCP",40,0,40); h_reco_cut_MCP->Sumw2();
	t_MCP->Project("h_reco_cut_MCP","Dpt",(TCut)Form("TotalWeight*%s*(%s)",S_reWeight.Data(),DataCuts.Data() ));


	TH1D *h_reco_cut_MCNP=new TH1D("h_reco_cut_MCNP","h_reco_cut_MCNP",40,0,40); h_reco_cut_MCNP->Sumw2();
	t_MCNP->Project("h_reco_cut_MCNP","Dpt",(TCut)Form("TotalWeight*%s*(%s)",S_reWeight.Data(),DataCuts.Data() ));

	TH1D *hGen_MCP=(TH1D*)f_mc_Prompt->Get("hGen_pt");
	TH1D *hGen_MCNP=(TH1D*)f_mc_NonPrompt->Get("hGen_pt");

	double N_reco_MCP=h_reco_cut_MCP->Integral( h_reco_cut_MCP->FindBin(DptLow+0.001), h_reco_cut_MCP->FindBin(DptHigh-0.001) );
	double N_reco_MCNP=h_reco_cut_MCNP->Integral( h_reco_cut_MCNP->FindBin(DptLow+0.001), h_reco_cut_MCNP->FindBin(DptHigh-0.001) );
	double N_Gen_MCP=hGen_MCP->Integral(hGen_MCP->FindBin(DptLow+0.001), hGen_MCP->FindBin(DptHigh-0.001) );
	double N_Gen_MCNP=hGen_MCNP->Integral(hGen_MCNP->FindBin(DptLow+0.001), hGen_MCNP->FindBin(DptHigh-0.001) );

	// cout<<"N_reco_MCP = "<<N_reco_MCP<<endl;

	double Eff_MCP=N_reco_MCP/N_Gen_MCP;
	double Eff_MCNP=N_reco_MCNP/N_Gen_MCNP;

	cout<<"Eff_MCP = "<<Eff_MCP<<endl;
	cout<<"Eff_MCNP = "<<Eff_MCNP<<endl;


	TH1D *hBtoDsCSdNdpt_Integral=(TH1D*)hBtoDsCSdNdpt_temp->Clone("hBtoDsCSdNdpt_Integral");
	MutiplyBinWidth(hBtoDsCSdNdpt_Integral);

	double CS_BtoDs=hBtoDsCSdNdpt_Integral->Integral(hBtoDsCSdNdpt_Integral->FindBin(DptLow+0.001), hBtoDsCSdNdpt_Integral->FindBin(DptHigh-0.001))/(DptHigh-DptLow);
	double N_yield=hRawYield->GetBinContent(1)/(DptHigh-DptLow);

  double CS_PromptDs = (N_yield*PhiRatio/(2*LumiNevt) - (CS_BtoDs* (BRphi*Eff_MCNP ) ) ) / ( BRphi*Eff_MCP );
	cout<<"CS_PromptDs = "<<CS_PromptDs<<endl;
	
 	double fr_Prompt=CS_PromptDs/(CS_PromptDs+CS_BtoDs);
	cout<<"fr_Prompt = "<<fr_Prompt<<endl;

	if(!fixPNPRatio){
		PNPRatio_Glo=fr_Prompt;
	}

	if(fixPNPRatio){
    CS_PromptDs=(N_yield*PhiRatio/(2*LumiNevt)) / ( ( BRphi*Eff_MCP ) + (1.0-PNPRatio)/(PNPRatio)*(BRphi*Eff_MCNP ) );
	cout<<"CS_PromptDs new = "<<CS_PromptDs<<" , PNPRatio = "<<PNPRatio<<endl;
	}


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


int CalCutScan_temp(int isPbPb=0, double DptLow=10, double DptHigh=20, TString var_scan="Ddls", int FixShape=1){

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
	if(reWeight==1){
		s_reWeightFname="_"+s_reWeight;
		cout<<"s_reWeightFname = "<<s_reWeightFname<<endl;
	}

  TString mcName_Prompt=Form("/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_cutScan/%sMC_phiPrompt_fitFile%s.root",s_PbPb3.Data(),s_reWeightFname.Data());
  TString mcName_NonPrompt=Form("/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_cutScan/%sMC_phiNonPrompt_fitFile%s.root",s_PbPb3.Data(),s_reWeightFname.Data());

	TFile *f_MCP=TFile::Open(mcName_Prompt.Data());
	TFile *f_MCNP=TFile::Open(mcName_NonPrompt.Data());

	// TTree *t_MCP=(TTree*)f_MCP->Get("t_fit");
	// t_MCP->Draw("Dmass");
	

	TFile *f_fit[nbin_var];
	TH1D *h_RawYield[nbin_var];
	int index_best=0;
	double bestErr=1;
	
	double Dalpha_cut=0.12;
	// double Dchi2cl_cut=0.1;
	// double Dchi2cl_cut=0.04;
	double Dchi2cl_cut=0.25;
	double Ddls_cut=2.5;
	double DtktkResmass_cut=DtktkResmassCutWidth;
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


	int fixPNPRatio=0;
	double PNPRatio=0.9;
	double PhiRatio=0.95;	
	if(isPbPb){
	PhiRatio=0.91;
	}

	// 1. read FitYield and GetThe Best

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
 

    TFile *f_NonPromptDs=TFile::Open("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BtoDs_Results/output/BtoDs.root");
    TH1D *hBtoDsCrossSectionPP_AnaBin_pythiaWeight=(TH1D*)f_NonPromptDs->Get("hBtoDsCrossSectionPP_AnaBin_pythiaWeight");
     TH1D *hBtoDsdNdPtPbPb_AnaBin_pythiaWeight=(TH1D*)f_NonPromptDs->Get("hBtoDsdNdPtPbPb_AnaBin_pythiaWeight");
    TH1D *hBtoDs_AnaBin_pythiaWeight=hBtoDsCrossSectionPP_AnaBin_pythiaWeight;  // this is differential cross section
    if(isPbPb){ hBtoDs_AnaBin_pythiaWeight=hBtoDsdNdPtPbPb_AnaBin_pythiaWeight; }

	if(var_scan=="Ddls"){
		Ddls_cut=bins_var[index_best];
	}
	if(var_scan=="Dchi2cl"){
		Dchi2cl_cut=bins_var[index_best];
	}


	h_CS[index_best]=h_CSCal_fun(isPbPb,h_RawYield[index_best],hBtoDs_AnaBin_pythiaWeight ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut,Dchi2cl_cut,Ddls_cut,DtktkResmass_cut, fixPNPRatio, PNPRatio, PhiRatio);
	// fixPNPRatio=1;
	// h_CS[index_best]=h_CSCal_fun(isPbPb,h_RawYield[index_best],hBtoDs_AnaBin_pythiaWeight ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut,Dchi2cl_cut,Ddls_cut,DtktkResmass_cut, fixPNPRatio, PNPRatio_Glo, PhiRatio);

	cout<<"PNPRatio_Glo = "<<PNPRatio_Glo<<endl;

	fixPNPRatio=1;
	for(int i=0; i<nbin_var; i++){
		if(i==index_best)continue;

		if(var_scan=="Ddls"){
			Ddls_cut=bins_var[i];
		}
		if(var_scan=="Dchi2cl"){
			Dchi2cl_cut=bins_var[i];
		}

		h_CS[i]=h_CSCal_fun(isPbPb,h_RawYield[i],hBtoDs_AnaBin_pythiaWeight ,f_MCP, f_MCNP, DptLow,DptHigh,Dalpha_cut,Dchi2cl_cut,Ddls_cut,DtktkResmass_cut, fixPNPRatio, PNPRatio_Glo, PhiRatio);
	}

	TGraphErrors *gr=new TGraphErrors();
	double BestCS=h_CS[index_best]->GetBinContent(1);
	double BestCSRelErr=h_CS[index_best]->GetBinError(1)/h_CS[index_best]->GetBinContent(1);
	
	for(int i=0; i<nbin_var; i++){
		// if(i==index_best)continue;
//    if(var_scan=="Ddls"){
//      Ddls_cut=bins_var[i];
//    }
		double CS=h_CS[i]->GetBinContent(1);
		double CSRelErr=h_CS[i]->GetBinError(1)/h_CS[i]->GetBinContent(1);
		double CSRatio=CS/BestCS;
		double CSRatioErr=sqrt(abs(CSRelErr*CSRelErr-BestCSRelErr*BestCSRelErr))*CSRatio;
		// double CSRatioErr=CSRelErr*CSRatio;

		gr->SetPoint(i,bins_var[i] ,CSRatio);
		gr->SetPointError(i,0.,CSRatioErr);		

	}
	// gr->GetXaxis()->SetRangeUser(0,4.5);
	gStyle->SetOptStat(0);
	int nbinTemp=1;
	double binsTempLow=0;
	double binsTempHigh=4.5;
	if(isPbPb){binsTempHigh=5.5;}
	if(var_scan=="Dchi2cl"){
		binsTempHigh=0.3;
		if(isPbPb){
			binsTempHigh=0.5;
		}
	}

	TCanvas *c_out=new TCanvas("c_out","c_out");
	c_out->cd();

	TH1D *htemp=new TH1D("htemp",Form(";%s; #sigma_{Varied cut}/#sigma_{Best cut}",var_scan.Data()),nbinTemp,binsTempLow,binsTempHigh);
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

	// TF1 *f1=new TF1("f1","[0]+x*[1]");
	TF1 *f1=new TF1("f1","[0]+(x-[2])*[1]");
	f1->SetLineColor(2);
	f1->SetParameter(0,1);
	// f1->FixParameter(0,1);
	f1->SetParameter(1,0);
	f1->SetParameter(2,bins_var[index_best]);
	f1->FixParameter(2,bins_var[index_best]);
	f1->SetRange(0,binsTempHigh);

	gr->Fit("f1","LEMS0F");
	gr->Fit("f1","LEMS0F");
	
	f1->Draw("same");

	tltx->DrawLatexNDC(textposx+shiftX,textposy+shiftY,Form("%s : %.0f < D_{S} p_{T} < %.0f",s_PbPb.Data(),DptLow,DptHigh));


	gSystem->Exec(Form("mkdir -p plots/"));
	c_out->SaveAs(Form("./plots/%s_%s_Dpt%.0fto%.0f%s.png",s_PbPb.Data(),var_scan.Data(),DptLow,DptHigh,s_FixShape.Data()));

	// calculate sys
	int nbin_pt=nbin_pt_pp;
	double *bins_pt=bins_pt_pp;
	if(isPbPb){
		nbin_pt=nbin_pt_PbPb3;
		bins_pt=bins_pt_PbPb3;
	}


	initParameter();
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




	return 0;

}


int main(int argc , char*argv[]){

	if(argc==5){

			CalCutScan_temp( atoi(argv[1]) , atof(argv[2]) , atof(argv[3]) , argv[4]  );
	}else{
		cout<<"wrong number of input"<<endl;
		return 1;
	}
	return 0;

}

