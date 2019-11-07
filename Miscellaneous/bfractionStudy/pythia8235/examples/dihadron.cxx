#define dihadron_cxx
#include "dihadron.h"
#include "constant.h"
//#include "pid.h"
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector.h>
#include <iostream>
#include <fstream>
#include <TString.h>
#include <TRandom3.h>

//#define PI TMath::Pi()
using namespace std;
/*
double dilu(double *x, double *par){
	double result;
	result = par[0]*(exp(-0.5*(par[1]*par[1]))+((sqrt(0.5*pi))*par[1]*cos(2*x[0])*exp(-0.5*(par[1]*par[1])*sin(2*x[0])*sin(2*x[0])))*(1+erf((par[1]*cos(2*x[0]))/sqrt(2))));
	return result;
}

double PijCal(double chi){
	TF1 *Dphi=new TF1("Dphi",dilu,-0.5*pi,0.5*pi,2);
	Dphi->SetParameter(0, 2./(2*pi));// this parameter should be 1/2*PI for initial form.
	Dphi->SetParameter(1, chi);//0.217974
	double r = Dphi->GetRandom();
	return r;

}

double Eta_Dis(){
	TF1 *tmp_eta=new TF1("tmp_eta","1",-1, 1);
	double result=	tmp_eta->GetRandom();
	//result = 2*(result - 0.5);
	//double result = eta_gen->Rndm();
	//result = 2*(result - 0.5);
	return result;
}


double Phi_Dis(double v2){
	TF1 *tmp_phi=new TF1("tmp_phi","1+2*[0]*cos(2*x)",-pi, pi);
	tmp_phi->SetParameter(0, v2);
	double result=	tmp_phi->GetRandom();
	//result = 2*(result - 0.5);
	return result;
}

double Pt_Dis( double para){
	TF1 *tmp_pt=new TF1("tmp_pt","x*exp(-x/0.3)",0.05, 15);
	double result=	tmp_pt->GetRandom();
	//tmp_pt->SetParameter(0, para);
	//TF1 *tmp_pt=new TF1("tmp_pt","1",0, 1);
	//double rnd = tmp_pt->GetRandom();
	//double result = -0.3*log(1-rnd);


	return result;
}
*/



void dihadron::Loop()
{
	if (fChain == 0) return;
	TFile *infile;
	infile = new TFile("px.root");
	//if(!(PX_read_switch)){
	if(1){
		TH1F *read_px_1, *read_px_2, *read_px_3, *read_px_4;	
		read_px_1 = (TH1F*)infile->Get("PX_1");
		read_px_2 = (TH1F*)infile->Get("PX_2");
		read_px_3 = (TH1F*)infile->Get("PX_3");
		read_px_4 = (TH1F*)infile->Get("PX_4");
		cout<<"check read in procedure start"<<endl;

		float en1 = 0.1*px_cutoff*(read_px_1->GetEntries());
		int sum1=0;
		float x1=0;
		for(int k=1; k<10001; k++){
			sum1 += (read_px_1->GetBinContent(k));
			x1 = read_px_1->GetBinCenter(k);
			if(sum1 > en1) break;
		}
		px_cut1 = x1;

		float en2 = 0.1*px_cutoff*(read_px_2->GetEntries());
		int sum2=0;
		float x2=0;
		for(int k=1; k<10001; k++){
			sum2 += (read_px_2->GetBinContent(k));
			x2 = read_px_2->GetBinCenter(k);
			if(sum2 > en2) break;
		}
		px_cut2 = x2;

		float en3 = 0.1*px_cutoff*(read_px_3->GetEntries());
		int sum3=0;
		float x3=0;
		for(int k=1; k<10001; k++){
			sum3 += (read_px_3->GetBinContent(k));
			x3 = read_px_3->GetBinCenter(k);
			if(sum3 > en3) break;
		}
		px_cut3 = x3;

	float en4 = 0.1*px_cutoff*(read_px_4->GetEntries());
		int sum4=0;
		float x4=0;
		for(int k=1; k<10001; k++){
			sum4 += (read_px_4->GetBinContent(k));
			x4 = read_px_4->GetBinCenter(k);
			if(sum4 > en4) break;
		}
		px_cut4 = x4;

		cout<<"read in procedure done"<<endl;
		//================
	}
	infile->Close();


	TFile *hfile = new TFile("Test.root","RECREATE");
	if(dphi_switch){	
		for(int k=0; k<NoAssocPt; k++){
			dphi_far[k]= new TH1F(Form("dphi_far%d",k),"dphi far region",phibinning,-1.5,-1.5+(2*pi));
			dphi_close[k]= new TH1F(Form("dphi_close%d",k),"dphi close region",phibinning,-1.5,-1.5+(2*pi));
			dphi_far[k]->Sumw2();
			dphi_close[k]->Sumw2();
		}
	}

	if(PX_read_switch){
		PX_1=new TH1F("PX_1","",10000,-600,1);
		PX_2=new TH1F("PX_2","",10000,-600,1);
		PX_3=new TH1F("PX_3","",10000,-600,1);
		PX_4=new TH1F("PX_4","",10000,-600,1);
	}
	TrigTot = new TH1F("TrigTot","Number of Trig.Particles",1,0,1);
	TrigTot->Sumw2();



	TH1D *Nch  = new TH1D ("Nch","Nch",3000,0.,3000.);
	TH1D *hpt  = new TH1D ("hpt"," ",500,0.,20.);
	TH1D *hphi  = new TH1D ("hphi"," ",200,-2*pi,2*pi);
	v2_check= new TProfile("v2_check","",10, 0, 10, -100, 100);
	TH1D *heta = new TH1D ("heta"," ",1000,-5.,5.);


	Double_t nentries = Int_t(fChain->GetEntries());
	Long64_t nbytes = 0, nb = 0;

	gRandom->SetSeed(0);

	TF1 *tmp_eta=new TF1("tmp_eta","1",-1, 1);
	TF1 *tmp_phi=new TF1("tmp_phi","1+2*[0]*cos(2*x)",-pi, pi);
	tmp_phi->SetParameter(0, 0.05);
	TF1 *tmp_pt=new TF1("tmp_pt","x*exp(-x/0.3)",0.05, 15);

	TRandom3 eta_gen;
	TRandom3 phi_gen;
	TRandom3 pt_gen;
	float eta_buff[ntracks];
	float pt_buff[ntracks];
	float phi_buff[ntracks];
	for (Long64_t jentry=0; jentry< nEvents;jentry++) {
		TotTrig = 0;

		if(jentry%10 == 0){
			cout<<">>"<<(100*double(jentry))/(nEvents) <<"\% Events accomplished"<<endl;
		}
	
		for(Int_t i=0; i<ntracks; i++){//1st loop
			float eta = tmp_eta->GetRandom();
			float pt = tmp_pt->GetRandom();
			float phi = tmp_phi->GetRandom();
			float v2 = cos(2*phi);
			eta_buff[i] = eta;
			pt_buff[i] = pt;
			phi_buff[i] = phi;
			hpt->Fill(pt, 1/pt);
			heta->Fill(eta);
			hphi->Fill(phi);
			v2_check->Fill(0.5, v2);
		}
		//continue;
		//cout<<"chech"<<endl;

		for(int i = 0; i<ArrayLength; i++){//initialize the array
			Phi_T_M[i] = -999;
			Pt_T_M[i] = -999;
			Eta_T_M[i] = -999;
			Phi_T_P[i] = -999;
			Pt_T_P[i] = -999;
			Eta_T_P[i] = -999;
		}

		TotTrig1=0;
		TotTrig2=0;
		TotTrig3=0;
		TotTrig4=0;

		for(Int_t i=0; i<ntracks; i++){//2nd loop
			float pt = pt_buff[i];
			float eta = eta_buff[i]; 
			float phi = phi_buff[i];
			if(pt < 3 || pt > 10) continue;

			px_tmp1 = 0;
			px_tmp2 = 0;
			px_tmp3 = 0;
			px_tmp4 = 0;
			for(Int_t j=0; j<ntracks; j++){//for px calculation
				if(i == j) continue;
				float npt = pt_buff[j];
				float neta = eta_buff[j]; 
				float nphi = phi_buff[j];
				float dPhi= phi-nphi;
				if(cos(dPhi) >= 0) continue;
				if(fabs(neta) > 1) continue;
				if(npt < 0.2 || npt > 10) continue;
				float px_tmp= npt*cos(dPhi);
				if(neta < -0.5){
					px_tmp1 += px_tmp;
				}
				if(neta > -0.5 && neta < 0.){
					px_tmp2 += px_tmp;
				}
				if(neta > 0. && neta < 0.5){
					px_tmp3 += px_tmp;
				}
				if(neta > 0.5 ){
					px_tmp4 += px_tmp;
				}

			}
			float limit_1;
			float limit_2;
			float limit_3;
			float limit_4;
			//---------------------------------------------------------------------

			if(PX_read_switch){
				PX_1->Fill(px_tmp1);
				PX_2->Fill(px_tmp2);
				PX_3->Fill(px_tmp3);
				PX_4->Fill(px_tmp4);
				limit_1= 0;
				limit_2= 0;
				limit_3= 0;
				limit_4= 0;
			}		
			else{
				limit_1= px_cut1;
				limit_2= px_cut2;
				limit_3= px_cut3;
				limit_4= px_cut4;
			}

			if((px_tmp1 > limit_1) && (px_tmp4 > limit_4)) continue;
			if((px_tmp1 <= limit_1) && (px_tmp4 <= limit_4)) continue;

			TotTrig += 1;
			if((px_tmp1 <= limit_1) && (px_tmp4 >  limit_4)){//minus
				Phi_T_M[TotTrig1]=phi;
				Pt_T_M[TotTrig1]=pt;
				Eta_T_M[TotTrig1]=eta;
				TotTrig1++;
			}
			if(px_tmp4<= limit_4 && px_tmp1 > limit_1)//plus
			{
				Phi_T_P[TotTrig4]=phi;
				Pt_T_P[TotTrig4]=pt;
				Eta_T_P[TotTrig4]=eta;
				TotTrig4++;
			} 

		}//end of 2nd loop

		TrigTot->Fill(0.5,TotTrig);
		float delta_phi, delta_eta, delta_pt;
		float tPhi, tPt, tEta;
		for(int i=0; i<TotTrig1; i++){//same event correlation loop1;
			tPhi=Phi_T_M[i];
			tPt=Pt_T_M[i];
			tEta=Eta_T_M[i];
			for(int j=0; j<ntracks; j++){
				float a_pt = pt_buff[j];
				float a_eta = eta_buff[j]; 
				float a_phi = phi_buff[j];
				if(method_ab){if(fabs(a_eta) > 0.5) continue;}
				if(method_cd){if(fabs(a_eta) < 0.5) continue;}
				if(a_pt < 0.2 || a_pt > 10) continue;
				delta_phi = tPhi-a_phi;
				delta_eta = tEta-a_eta;
				delta_pt = tPt-a_pt;
				if(fabs(delta_phi)< 1e-6 && fabs(delta_eta) <1e-6 && fabs(delta_pt) <1e-6) continue;
				int pt_region= GetPtRegion(a_pt);
				if(delta_phi < -pi) delta_phi += 2*pi;//folding
				if(delta_phi > pi)  delta_phi -= 2*pi;
				if(delta_phi < -1.5)  delta_phi += 2*pi;
				if(a_eta > 0.){dphi_far[pt_region]->Fill(delta_phi);} 
				if(a_eta < 0.){dphi_close[pt_region]->Fill(delta_phi);}
			}
		}
		for(int i=0; i<TotTrig4; i++){//same event correlation loop4;
			tPhi=Phi_T_P[i];
			tPt=Pt_T_P[i];
			tEta=Eta_T_P[i];
			int s_region= GetPhiSRegion(tPhi, psi2, No_PhiS);
			for(int j=0; j<ntracks; j++){
				float a_pt = pt_buff[j];
				float a_eta = eta_buff[j]; 
				float a_phi = phi_buff[j];

				if(method_ab){if(fabs(a_eta) > 0.5) continue;}
				if(method_cd){if(fabs(a_eta) < 0.5) continue;}
				if(a_pt < 0.2 || a_pt > 10) continue;
				delta_phi = tPhi-a_phi;
				delta_eta = tEta-a_eta;
				delta_pt = tPt-a_pt;
				if(fabs(delta_phi)< 1e-6 && fabs(delta_eta) <1e-6 && fabs(delta_pt) <1e-6) continue;
				int pt_region= GetPtRegion(a_pt);
				if(delta_phi < -pi) delta_phi += 2*pi;//folding
				if(delta_phi > pi)  delta_phi -= 2*pi;
				if(delta_phi < -1.5)  delta_phi += 2*pi;
				if(a_eta < 0.){dphi_far[pt_region]->Fill(delta_phi);} 
				if(a_eta > 0.){dphi_close[pt_region]->Fill(delta_phi);}
			}
		}


	}//end of event loop;
	hfile->Write();
	hfile->Close();
}


int main(int argc, char **argv)
{
	dihadron t;
	t.Loop();

	return 0;
	// return **argc;
}
