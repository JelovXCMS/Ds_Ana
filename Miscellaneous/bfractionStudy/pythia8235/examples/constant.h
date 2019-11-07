#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TMath.h"
#include "TFile.h"

#define psi2_switch               2 
//0 for bbc, 1 for tpc, 2 for default
#define method_ab                 1
#define method_cd                 0
#define px_cutoff                 1
// 1-9 for 10-90% tail of all events

#define No_Centr                  9
#define No_PhiS                   1
#define NoAssocPt                 9

#define PX_read_switch            1
#define dphi_switch               1
#define phibinning                48//24;
#define etabinning                80//80;//48;

const float PtLow[NoAssocPt]= {0.15,0.5,1.0,1.5,2.0,2.5,3.0,4.0,6.0};
const float PtHigh[NoAssocPt]={0.5,1.0,1.5,2.0,2.5,3.0,4.0,6.0,10.0};

int cen_flag=-1;
int centr[No_Centr];
const int ArrayLength = 1000;
int nEvents = 100000;
int ntracks = 1000;


TH1F *dphi_far[NoAssocPt];
TH1F *dphi_close[NoAssocPt];
TH1F *mdphi_far[NoAssocPt];
TH1F *mdphi_close[NoAssocPt];
// for MC
TF1* tmp_eta;
TF1* tmp_phi;
TF1* tmp_pt;

//



TH1F *PX_1;
TH1F *PX_2;
TH1F *PX_3;
TH1F *PX_4;

TProfile *bbc_res_raw;
TProfile *bbc_res_sub;
TProfile *bbc_res_corr;
TProfile *EP_v2_L;
TProfile *RP_v2_L;
TProfile *EP_v2_H;
TProfile *RP_v2_H;
TProfile *Ratio_v2_L;
TProfile *Ratio_v2_H;

TProfile *Px_EP_v2_L;
TProfile *Px_RP_v2_L;
TProfile *Px_EP_v2_H;
TProfile *Px_RP_v2_H;
TProfile *Px_Ratio_v2_L;
TProfile *Px_Ratio_v2_H;
TProfile *v2_check;

TH2D *AMPT_histMdetGenMC[No_Centr];
TH1F *Phi_sp_folded[No_Centr];
TH1F *Phi_sp_rp[No_Centr];
TH1F *Phi_sp_raw[No_Centr];

TH1F *F_Phi_sp_folded[No_Centr];
TH1F *F_Phi_sp_rp[No_Centr];


TH1F *TrigTot;
int TotTrig;

float px_cut1;
float px_cut2;
float px_cut3;
float px_cut4;
float px_tmp1;
float px_tmp2;
float px_tmp3;
float px_tmp4;

float bbcpsi2;
float bbcres_sub;
float bbcres_full;
float tpcpsi2;
float psi2;
int TotTrig1;
int TotTrig2;
int TotTrig3;
int TotTrig4;

float Phi_T_M[ArrayLength];
float Pt_T_M[ArrayLength];
float Eta_T_M[ArrayLength];
float Phi_T_P[ArrayLength];
float Pt_T_P[ArrayLength];
float Eta_T_P[ArrayLength];



