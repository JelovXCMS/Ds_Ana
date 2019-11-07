
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


int smearTree(TFile *f_mcp_new,TFile *f_mc, int nSmear=20, double DptLow=2,double DptHigh=40, int isPbPb=0){


    double BestScale;

	double BestScales_pp[]={0.99,1.06,1.07,1.09,1.07,1.14,1.14,1.29};
	double BestScales_PbPb[]={1.0,1.0,1.0,1.0,1.08,1.16,1.14,1.17};
	double *BestScales=BestScales_pp;
	if(isPbPb){
			BestScales=BestScales_PbPb;
	}

  TRandom3 Rdm;
  Rdm.SetSeed(0);

  TTree *t_mcp=(TTree*)f_mc->Get("t_fit");
  // TH1D *hGen_pt=(TH1D*)f_mc->Get("hGen_pt");
  f_mcp_new->cd();
  // hGen_pt->Write("",TObject::kOverwrite);

	cout<<__LINE__<<endl;

  Float_t Dmass;
  Float_t Dpt;
  Float_t Ddls;
  Float_t Ddl;
  Float_t DdlErr;
  Float_t Dalpha;
  Float_t Dchi2cl;
  Float_t DtktkResmass;
  Float_t TotalWeight;

	Float_t Dtrk1Pt;
	Float_t Dtrk2Pt;
	Float_t Dtrk3Pt;

  Int_t   DsGen;
  Int_t   DgencollisionId;
  Float_t DgenBAncestorpt;
  Float_t RecoFONLLWeight;
  Float_t RecoD0DataWeight;
  Float_t RecoDsDataWeight;
  Float_t weight;
  Float_t DgenptSampleWeight;
  Float_t Ncoll;
  Float_t PbPbVzWeight;
  Float_t RecoFONLLRaaWeight ;

  Float_t DlxyBS;
  Float_t DlxyBSErr;
  Float_t DlxyBSs;
  Float_t Dalpha_BS_2D;

  t_mcp->SetBranchAddress("Dmass",&Dmass);
  t_mcp->SetBranchAddress("Ddls",&Ddls);
  t_mcp->SetBranchAddress("Dpt",&Dpt);
  t_mcp->SetBranchAddress("DdlErr",&DdlErr);
  t_mcp->SetBranchAddress("Ddl",&Ddl);
  // t_mcp->SetBranchAddress("DtktkResmass",&DtktkResmass);
  t_mcp->SetBranchAddress("Dalpha",&Dalpha);
  t_mcp->SetBranchAddress("Dchi2cl",&Dchi2cl);
  // t_mcp->SetBranchAddress("TotalWeight",&TotalWeight);
  t_mcp->SetBranchAddress("TotalWeight",&TotalWeight);

/*
t_mcp->SetBranchAddress("Dtrk1Pt"            ,&       Dtrk1Pt);
t_mcp->SetBranchAddress("Dtrk2Pt"            ,&       Dtrk2Pt);
t_mcp->SetBranchAddress("Dtrk3Pt"            ,&       Dtrk3Pt);
t_mcp->SetBranchAddress("DsGen"              ,&       DsGen);
t_mcp->SetBranchAddress("DgencollisionId"    ,&       DgencollisionId);
t_mcp->SetBranchAddress("DgenBAncestorpt"    ,&       DgenBAncestorpt);
t_mcp->SetBranchAddress("RecoFONLLWeight"    ,&       RecoFONLLWeight);
t_mcp->SetBranchAddress("RecoD0DataWeight"   ,&       RecoD0DataWeight);
t_mcp->SetBranchAddress("RecoDsDataWeight"   ,&       RecoDsDataWeight);
t_mcp->SetBranchAddress("weight"             ,&       weight);
t_mcp->SetBranchAddress("DgenptSampleWeight" ,&       DgenptSampleWeight);
t_mcp->SetBranchAddress("Ncoll"              ,&       Ncoll);
t_mcp->SetBranchAddress("PbPbVzWeight"       ,&       PbPbVzWeight);
t_mcp->SetBranchAddress("RecoFONLLRaaWeight" ,&      RecoFONLLRaaWeight) ;
*/


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
  // t_new->Branch("DtktkResmass",&DtktkResmass);
  t_new->Branch("TotalWeight",&TotalWeight);
  t_new->Branch("SMrWt",&SMrWt);

/*
t_new->Branch("Dtrk1Pt"            ,&       Dtrk1Pt);
t_new->Branch("Dtrk2Pt"            ,&       Dtrk2Pt);
t_new->Branch("Dtrk3Pt"            ,&       Dtrk3Pt);
t_new->Branch("DsGen"              ,&       DsGen);
t_new->Branch("DgencollisionId"    ,&       DgencollisionId);
t_new->Branch("DgenBAncestorpt"    ,&       DgenBAncestorpt);
t_new->Branch("RecoFONLLWeight"    ,&       RecoFONLLWeight);
t_new->Branch("RecoD0DataWeight"   ,&       RecoD0DataWeight);
t_new->Branch("RecoDsDataWeight"   ,&       RecoDsDataWeight);
t_new->Branch("weight"             ,&       weight);
t_new->Branch("DgenptSampleWeight" ,&       DgenptSampleWeight);
t_new->Branch("Ncoll"              ,&       Ncoll);
t_new->Branch("PbPbVzWeight"       ,&       PbPbVzWeight);
t_new->Branch("RecoFONLLRaaWeight" ,&      RecoFONLLRaaWeight) ;
*/



/*
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
*/
  float Ddls_BestScale;
  float Ddl_BestScale;
  float DdlErr_BestScale;

    t_new->Branch(Form("Ddls_BestScale"),&Ddls_BestScale);
    t_new->Branch(Form("Ddl_BestScale"),&Ddl_BestScale);
    t_new->Branch(Form("DdlErr_BestScale"),&DdlErr_BestScale);
/*
  float Ddls_BestScale[nScales];
  float Ddl_BestScale[nScales];
  float DdlErr_BestScale[nScales];


  for(int i=0; i<nScales; i++){
    t_new->Branch(Form("Ddls_Scl%.0fem3",BestScales[i]*1e3),&Ddls_BestScale[i]);
    t_new->Branch(Form("Ddl_Scl%.0fem3",BestScales[i]*1e3),&Ddl_BestScale[i]);
    t_new->Branch(Form("DdlErr_Scl%.0fem5",BestScales[i]*1e3),&DdlErr_BestScale[i]);
  }
*/


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
		if(isPbPb && Dpt<6) continue;
		if(Dalpha > 0.2) continue;
		// if(DgencollisionId!=0) continue;


			BestScale=1;
		if(Dpt<3){
			BestScale=BestScales[0];
		}else if(Dpt<4){
			BestScale=BestScales[1];	
		}else if(Dpt<5){
			BestScale=BestScales[2];	
		}else if(Dpt<6){
			BestScale=BestScales[3];	
		}else if(Dpt<8){
			BestScale=BestScales[4];	
		}else if(Dpt<10){
			BestScale=BestScales[5];	
		}else if(Dpt<20){
			BestScale=BestScales[6];	
		}else if(Dpt<40){
			BestScale=BestScales[7];	
		}

    for(int j=0; j<nSmear;j++){
/*
        for(int k=0;k<nScales;k++){
        DdlErr_BestScale[k]=DdlErr*BestScales[k];
        Ddls_BestScale[k]=Ddl / DdlErr_BestScale[k];
        Ddl_BestScale[k]=Ddl ;
        if(BestScales[k]>1){
          Ddls_BestScale[k]=Rdm.Gaus(Ddls_BestScale[k],sqrt(1.0*1.0-1.0/BestScales[k]*1.0/BestScales[k]*1.0));
          Ddl_BestScale[k]=Rdm.Gaus(Ddl_BestScale[k],DdlErr_BestScale[k]*sqrt(1.0*1.0-1.0/BestScales[k]*1.0/BestScales[k]*1.0));
        }
        }
*/
        DdlErr_BestScale=DdlErr*BestScale;
        Ddls_BestScale=Ddl / DdlErr_BestScale;
        Ddl_BestScale=Ddl ;
        if(BestScale>1){
          Ddls_BestScale=Rdm.Gaus(Ddls_BestScale,sqrt(1.0*1.0-1.0/BestScale*1.0/BestScale*1.0));
          Ddl_BestScale=Rdm.Gaus(Ddl_BestScale,DdlErr_BestScale*sqrt(1.0*1.0-1.0/BestScale*1.0/BestScale*1.0));
        }

      t_new->Fill();
    } // end for nsmear

  } // end nentries;

  f_mcp_new->cd();
  t_new->Write("",TObject::kOverwrite);

  return 0;


}








int Make_MC_DlsScale(int isPbPb=3, int nSmear=20){

	 double DptLow=2;
	 double DptHigh=40;

	 TString s_PbPb3="pp";
	
	if(isPbPb){
		DptLow=6;
		DptHigh=40;
		s_PbPb3="PbPb3";
	
	}
	  

    TString mcName_Prompt_old=Form("./rootF/%sMC_phiPrompt_fitFile.root",s_PbPb3.Data());
    TString mcName_Prompt_new=Form("./rootF/%sMC_phiPrompt_fitFile_DdlsScl.root",s_PbPb3.Data());
    // TString mcName_NonPrompt_new=Form("/scratch/halstead/p/peng43/Ds_phikkpi/FitFile_cutScan/%sMC_phiNonPrompt_fitFile%s_smear.root",s_PbPb3.Data(),s_reWeightFname.Data());

		TFile *f_MCP=TFile::Open(mcName_Prompt_old.Data());
    TFile *f_MCP_new=new TFile(mcName_Prompt_new.Data(),"recreate");
    smearTree(f_MCP_new,f_MCP,nSmear,DptLow,DptHigh,isPbPb);
    // TFile *f_mcnp_new=new TFile(mcName_NonPrompt_new.Data(),"recreate");
    // smearTree(f_mcnp_new,f_MCNP,nSmear,DptLow,DptHigh,isPbPb);
		// f_MCP_new->Write();
		TH1D *hMCP_Gen_pt=(TH1D*)f_MCP->Get("hGen_pt");	
		TH1D *hMCP_Gen_pt_D0Weight=(TH1D*)f_MCP->Get("hGen_pt_D0Weight");	
		TH1D *hMCP_Dpt=(TH1D*)f_MCP->Get("hDpt");	
		TH1D *hMCP_EffDpt=(TH1D*)f_MCP->Get("hEffDpt");	
		f_MCP_new->cd();
		hMCP_Gen_pt->Write();
		hMCP_Gen_pt_D0Weight->Write();
		hMCP_Dpt->Write();
		hMCP_EffDpt->Write();

		f_MCP_new->Close();
		
		cout<<"\nprompt smear done\n"<<endl;
		

    TString mcName_NonPrompt_old=Form("./rootF/%sMC_phiNonPrompt_fitFile.root",s_PbPb3.Data());
    TString mcName_NonPrompt_new=Form("./rootF/%sMC_phiNonPrompt_fitFile_DdlsScl.root",s_PbPb3.Data());

		TFile *f_MCNP=TFile::Open(mcName_NonPrompt_old.Data());
    TFile *f_MCNP_new=new TFile(mcName_NonPrompt_new.Data(),"recreate");
    smearTree(f_MCNP_new,f_MCNP,nSmear,DptLow,DptHigh,isPbPb);

		TH1D *hMCNP_Gen_pt=(TH1D*)f_MCNP->Get("hGen_pt");	
		TH1D *hMCNP_Gen_pt_D0Weight=(TH1D*)f_MCNP->Get("hGen_pt_D0Weight");	
		TH1D *hMCNP_Dpt=(TH1D*)f_MCNP->Get("hDpt");	
		TH1D *hMCNP_EffDpt=(TH1D*)f_MCNP->Get("hEffDpt");	
		f_MCNP_new->cd();
		hMCNP_Gen_pt->Write();
		hMCNP_Gen_pt_D0Weight->Write();
		hMCNP_Dpt->Write();
		hMCNP_EffDpt->Write();

	  f_MCNP_new->Close();

		cout<<"\nnonprompt smear done\n"<<endl;


// int smearTree(TFile *f_mcp_new,TFile *f_mc, int nSmear=20, double DptLow=2,double DptHigh=40, int isPbPb=0){



	return 1;

}





