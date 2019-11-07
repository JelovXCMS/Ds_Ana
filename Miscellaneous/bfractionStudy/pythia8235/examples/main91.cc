// main91.cc is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It studies the charged multiplicity distribution at the LHC.
// Modified by Rene Brun, Axel Naumann and Bernhard Meirose
// to use ROOT for histogramming.

// Stdlib header file for input and output.
#include <iostream>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
//#include "constant.h"
// ROOT, for histogramming.
#include "TH1.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"


#define NoAssocPt                 9

using namespace Pythia8;

int main(int argc, char* argv[]) {

	// Create the ROOT application environment.
	TApplication theApp("hist", &argc, argv);

	// Create Pythia instance and set it up to generate hard QCD processes
	// above pTHat = 20 GeV for pp collisions at 14 TeV.
	Pythia pythia;
	pythia.readString("HardQCD:all = on");
	pythia.readString("PhaseSpace:pTHatMin = 20.");
	//pythia.readString("PhaseSpace:pTHatMin = 0.");
	//pythia.readString("Beams:eCM = 14000.");
	pythia.readString("Beams:eCM = 200.");
	pythia.init();

	// Create file on which histogram(s) can be saved.
	TFile* outFile = new TFile("hist.root", "RECREATE");
	double PI = 3.14159265359;
	int phibinning = 48;
	const float PtLow[NoAssocPt]= {0.15,0.5,1.0,1.5,2.0,2.5,3.0,4.0,6.0};
	const float PtHigh[NoAssocPt]={0.5,1.0,1.5,2.0,2.5,3.0,4.0,6.0,10.0};
	const float ptLow[9]= {0.15,0.5,1.0,1.5,2.0,2.5,3.0,4.0,6.0};
	const float ptHigh[9]={0.5,1.0,1.5,2.0,2.5,3.0,4.0,6.0,10.0};

	double px_tmp1, px_tmp2, px_tmp3, px_tmp4;
	// Book histogram.
	TH1F *mult = new TH1F("mult","charged multiplicity", 100, -0.5, 799.5);
	TH1F *npt = new TH1F("npt","npt", 300, -0.5, 30);
	TH1F *neta = new TH1F("neta","neta", 100, -5, 5);
	TH1F *nphi = new TH1F("nphi","nphi", 100, -PI, PI);
	TH1F *dphi_far[NoAssocPt];
	TH1F *dphi_close[NoAssocPt];
	TH1F *PX_1, *PX_2, *PX_3, *PX_4;
	for(int k=0; k<NoAssocPt; k++){
		dphi_far[k]= new TH1F(Form("dphi_far%d",k),"dphi far region",phibinning,-1.5,-1.5+(2*PI));
		dphi_close[k]= new TH1F(Form("dphi_close%d",k),"dphi close region",phibinning,-1.5,-1.5+(2*PI));
		dphi_far[k]->Sumw2();
		dphi_close[k]->Sumw2();
	}
	TH1F *TrigTot;
	TrigTot = new TH1F("TrigTot","Number of Trig.Particles",1,0,1);
	TrigTot->Sumw2();

	PX_1=new TH1F("PX_1","",10000,-600,1);
	PX_2=new TH1F("PX_2","",10000,-600,1);
	PX_3=new TH1F("PX_3","",10000,-600,1);
	PX_4=new TH1F("PX_4","",10000,-600,1);


	// Begin event loop. Generate event; skip if generation aborted.
	for (int iEvent = 0; iEvent < 3000; ++iEvent) {
		if (!pythia.next()) continue;
		int nCharged = 0;
		int TotTrig =0;
		for (int i = 0; i < pythia.event.size(); ++i){//first particle loop starts here
			if (pythia.event[i].isFinal() && pythia.event[i].isCharged() &&(fabs(pythia.event[i].eta()) < 1.)){
				double tpt =  pythia.event[i].pT();
				double teta = pythia.event[i].eta();
				double tphi = pythia.event[i].phi();
				npt->Fill(tpt);
				neta->Fill(teta);
				nphi->Fill(tphi);
				++nCharged;


				if(tpt < 3 || tpt > 10) continue;
				px_tmp1 = 0;
				px_tmp2 = 0;
				px_tmp3 = 0;
				px_tmp4 = 0;
				for (int j = 0; j < pythia.event.size(); ++j){
					if (pythia.event[j].isFinal() && pythia.event[j].isCharged() &&(fabs(pythia.event[j].eta()) < 1.)){
						if(i == j) continue;
						double apt =  pythia.event[j].pT();
						double aeta = pythia.event[j].eta();
						double aphi = pythia.event[j].phi();
						double dphi = tphi - aphi;
						if(cos(dphi) >= 0) continue;
						if(apt < 0.2 || apt > 10) continue;
						double px_tmp= apt*cos(dphi);
						if(aeta < -0.5){
							px_tmp1 += px_tmp;
						}
						if(aeta > -0.5 && aeta < 0.){
							px_tmp2 += px_tmp;
						}
						if(aeta > 0. && aeta < 0.5){
							px_tmp3 += px_tmp;
						}
						if(aeta > 0.5 ){
							px_tmp4 += px_tmp;
						}
					}}
				PX_1->Fill(px_tmp1);
				PX_2->Fill(px_tmp2);
				PX_3->Fill(px_tmp3);
				PX_4->Fill(px_tmp4);
				int eta_flag = -999;
				double pxtail = -8.6;
				if((px_tmp1 > pxtail) && (px_tmp4 > pxtail)) continue;
				if((px_tmp1 <= pxtail) && (px_tmp4 <= pxtail)) continue;
				TotTrig += 1;
				if((px_tmp1 <= pxtail) && (px_tmp4 >  pxtail)){// minus
					eta_flag = -1;
				}
				if(px_tmp4<= pxtail && px_tmp1 > pxtail){//plus
					eta_flag = 1;
				}
				for (int k = 0; k < pythia.event.size(); ++k){
					if (pythia.event[k].isFinal() && pythia.event[k].isCharged() &&(fabs(pythia.event[k].eta()) < 1.)){
						if(i == k) continue;
						double apt =  pythia.event[k].pT();
						double aeta = pythia.event[k].eta();
						double aphi = pythia.event[k].phi();
						double delta_phi = tphi - aphi;
						if(apt < 0.2 || apt > 10) continue;
						int pT_region= -999;
						for(int i=0;i<9;i++)
						{
							if(apt<ptHigh[i]&&apt>=ptLow[i]) { pT_region = i; break;} 
						}

						int pt_region= pT_region;
						if(delta_phi < -PI) delta_phi += 2*PI;//folding
						if(delta_phi > PI)  delta_phi -= 2*PI;
						if(delta_phi < -1.5)  delta_phi += 2*PI;
						if(eta_flag == 1){
							if(aeta < 0.){dphi_far[pt_region]->Fill(delta_phi);} 
							if(aeta > 0.){dphi_close[pt_region]->Fill(delta_phi);}
						}
						if(eta_flag == -1){
							if(aeta > 0.){dphi_far[pt_region]->Fill(delta_phi);} 
							if(aeta < 0.){dphi_close[pt_region]->Fill(delta_phi);}
						}

					}}
			}}//first particle loop ends here

		TrigTot->Fill(0.5,TotTrig);
		mult->Fill( nCharged );


	}

	// Statistics on event generation.
	pythia.stat();

	// Show histogram. Possibility to close it.
	//mult->Draw();
	//std::cout << "\nDouble click on the histogram window to quit.\n";
	//gPad->WaitPrimitive();

	// Save histogram on file and close file.
	outFile->Write();
	outFile->Close();
	//mult->Write();
	//npt->Write();
	//neta->Write();
	//nphi->Write();
	//delete outFile;

	// Done.
	return 0;
}
