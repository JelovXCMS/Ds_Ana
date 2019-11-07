#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

void TheoryPredict(){

	TFile *f_TAMU_result =TFile::Open("./output/TheoryPredict.root","recreate");

	double TAMU_DsD0_pp_val[]={0.2714,0.2766,0.27698,0.28098,0.28079,0.27921,0.27635,0.27206,0.26777,0.26487,0.26425,0.26261,0.26181,0.2615,0.25996,0.26089,0.26041,0.25978,0.26139,0.26087,0.26068,0.26036,0.26166,0.26032,0.26194,0.26071,0.26035,0.26024,0.26241,0.2605,0.26065,0.26078,0.26338,0.25959,0.25807,0.25882,0.26005,0.25891,0.26427,0.26071,0.26209,0.2607,0.26061,0.2608,0.26046,0.25654,0.25678,0.2538,0.25956,0.25026,0.25543,0.25523,0.25336,0.25433,0.25689,0.25895,0.25783,0.24663,0.25289,0.25809,};
	const int n_TAMU_DsD0_pp=sizeof(TAMU_DsD0_pp_val)/sizeof(TAMU_DsD0_pp_val[0]);
	double TAMU_DsD0_pp_pt[n_TAMU_DsD0_pp];
	for(int i =0; i<n_TAMU_DsD0_pp; i++){
		TAMU_DsD0_pp_pt[i]=i*0.2+0.1;
		cout<<"TAMU_DsD0_pp_pt = "<<TAMU_DsD0_pp_pt[i]<<" : "<<TAMU_DsD0_pp_val[i]<<endl;
	}

	TGraph *gr_TAMU_DsD0_pp=new TGraph(n_TAMU_DsD0_pp,TAMU_DsD0_pp_pt,TAMU_DsD0_pp_val); 
	gr_TAMU_DsD0_pp->SetName("gr_TAMU_DsD0_pp");

	f_TAMU_result->cd();
	gr_TAMU_DsD0_pp->Write();

	// gr_DsD0_pp->Draw("");

	double TAMU_DsD0_PbPb020_withHadDif_val[n_TAMU_DsD0_pp]={
 0.37885, 	
 0.38518,	
 0.38922,	
 0.39418,	
 0.39997,	
 0.40768,	
 0.41562,	
 0.42259,	
 0.43263,	
 0.44044,	
 0.45056,	
 0.46176,	
 0.47175,	
 0.482	,  
 0.49045,	
 0.49851,	
 0.50341,	
 0.50608,	
 0.50876,	
 0.50797,	
 0.50621,	
 0.50333,	
 0.49725,	
 0.49159,	
 0.4843	,  
 0.47556,	
 0.46736,	
 0.45791,	
 0.44869,	
 0.43961,	
 0.43168,	
 0.42448,	
 0.4172	,  
 0.40828,	
 0.40177,	
 0.39361,	
 0.38707,	
 0.38156,	
 0.37551,	
 0.37002,	
 0.3639	,  
 0.35742,	
 0.35386,	
 0.34692,	
 0.34178,	
 0.33837,	
 0.3357	,  
 0.33186,	
 0.3309	,  
 0.32728,	
	0.32394,	
	0.3238	,
	0.31855	,
	0.31309	,
	0.31159	,
	0.31001	,
	0.30904	,
	0.30184	,
	0.29813	,
	0.2989	
	};

	double TAMU_DsD0_PbPb020_noHadDif_val[n_TAMU_DsD0_pp]={
0.38847, 
0.39025,
0.39175,
0.39457,
0.39911,
0.40445,
0.40983,
0.41658,
0.42674,
0.43725,
0.44918,
0.46163,
0.47266,
0.48394,
0.49217,
0.50046,
0.50528,
0.50696,
0.50775,
0.50466,
0.50023,
0.49599,
0.48746,
0.47941,
0.47051,
0.46113,
0.4539,
0.44489,
0.43679,
0.42859,
0.42045,
0.4112,
0.40335,
0.39471,
0.38961,
0.38404,
0.37731,
0.37392,
0.3687,
0.36268,
0.35637,
0.3504,
0.34415,
0.33769,
0.33696,
0.33266,
0.33189,
0.32931,
0.326,
0.3225,
0.31765,
0.31472,
0.31297,
0.30544,
0.30504,
0.30667,
0.30563,
0.30216,
0.2949,
0.29755
	};

	TGraph *gr_TAMU_DsD0_PbPb020_noHadDif=new TGraph(n_TAMU_DsD0_pp,TAMU_DsD0_pp_pt,TAMU_DsD0_PbPb020_noHadDif_val); 
	gr_TAMU_DsD0_PbPb020_noHadDif->SetName("gr_TAMU_DsD0_PbPb020_noHadDif");

	TGraph *gr_TAMU_DsD0_PbPb020_withHadDif=new TGraph(n_TAMU_DsD0_pp,TAMU_DsD0_pp_pt,TAMU_DsD0_PbPb020_withHadDif_val); 
	gr_TAMU_DsD0_PbPb020_withHadDif->SetName("gr_TAMU_DsD0_PbPb020_withHadDif");

	gr_TAMU_DsD0_PbPb020_noHadDif->Write();
	gr_TAMU_DsD0_PbPb020_withHadDif->Write();


	double CuJet_pt[]={
8.0  	,	
8.5 	,	
9.0 	,	
9.5	  ,	
10.0	,	
10.5	,	
11.0	,	
11.5	,	
12.0	,	
12.5	,	
13.0	,	
13.5	,	
14.0	,	
14.5	,	
15.0	,	
15.5	,	
16.0	,	
16.5	,	
17.0	,	
17.5	,	
18.0	,	
18.5	,	
19.0	,	
19.5	,	
20.0	,	
20.5	,	
21.0	,	
21.5	,	
22.0	,	
22.5	,	
23.0	,	
23.5	,	
24.0	,	
24.5	,	
25.0	,	
25.5	,	
26.0	,	
26.5	,	
27.0	,	
27.5	,	
28.0	,	
28.5	,	
29.0	,	
29.5	,	
30.0	,	
30.5	,	
31.0	,	
31.5	,	
32.0	,	
32.5	,	
33.0	,	
33.5	,	
34.0	,	
34.5	,	
35.0	,	
35.5	,	
36.0	,	
36.5	,	
37.0	,	
37.5	,	
38.0	,	
38.5	,	
39.0	,	
39.5	,	
40.0		
}; 

	int nCuJet=sizeof(CuJet_pt)/sizeof(CuJet_pt[0]);

	double CuJet_Raa_alpha1cm0p28[]={
0.240636, 
0.248390,
0.255914,
0.263216,
0.270305,
0.277188,
0.283871,
0.290363,
0.296669,
0.302797,
0.308754,
0.314545,
0.320176,
0.325654,
0.330985,
0.336173,
0.341223,
0.346143,
0.350935,
0.355606,
0.360159,
0.364599,
0.368931,
0.373159,
0.377287,
0.381319,
0.385258,
0.389109,
0.392874,
0.396557,
0.400162,
0.403691,
0.407147,
0.410534,
0.413854,
0.417110,
0.420304,
0.423439,
0.426516,
0.429539,
0.432510,
0.435430,
0.438301,
0.441126,
0.443905,
0.446642,
0.449337,
0.451992,
0.454608,
0.457187,
0.459730,
0.462238,
0.464713,
0.467156,
0.469567,
0.471948,
0.474300,
0.476624,
0.478920,
0.481189,
0.483433,
0.485651,
0.487845,
0.490014,
0.492160
};

	double CuJet_Raa_alpha0p8cm0p22[]={
0.312488, 
0.320627,
0.328540,
0.336232,
0.343713,
0.350987,
0.358062,
0.364944,
0.371639,
0.378153,
0.384492,
0.390661,
0.396666,
0.402513,
0.408206,
0.413750,
0.419151,
0.424414,
0.429542,
0.434540,
0.439413,
0.444165,
0.448799,
0.453321,
0.457733,
0.462040,
0.466245,
0.470351,
0.474363,
0.478282,
0.482113,
0.485859,
0.489521,
0.493104,
0.496610,
0.500042,
0.503402,
0.506692,
0.509916,
0.513075,
0.516171,
0.519207,
0.522186,
0.525107,
0.527975,
0.530790,
0.533555,
0.536270,
0.538938,
0.541560,
0.544138,
0.546673,
0.549167,
0.551620,
0.554035,
0.556412,
0.558752,
0.561057,
0.563328,
0.565565,
0.567770,
0.569944,
0.572087,
0.574200,
0.576285
};


	TGraph *gr_Raa_CuJet_alpha1cm0p28=new TGraph(nCuJet , CuJet_pt , CuJet_Raa_alpha1cm0p28); 
	gr_Raa_CuJet_alpha1cm0p28->SetName("gr_CuJet_Raa_alpha1cm0p28");
	gr_Raa_CuJet_alpha1cm0p28->Write();

	TGraph *gr_Raa_CuJet_alpha0p8cm0p22=new TGraph(nCuJet , CuJet_pt , CuJet_Raa_alpha0p8cm0p22); 
	gr_Raa_CuJet_alpha0p8cm0p22->SetName("gr_CuJet_Raa_alpha0p8cm0p22");
	gr_Raa_CuJet_alpha0p8cm0p22->Write();
	

	double PHSD_Raa_pt[]={ 0.5,1.5,2.5,3.5,4.5,5.5,7,9,12,16};
	double PHSD_Raa_val[]={ 
0.743591099, 
1.088786742,
1.202523859,
0.844312333,
0.557596527,
0.394009335,
0.283843785,
0.212354049,
0.235129642,
0.258787857
};
	int nPHSD_Raa=sizeof(PHSD_Raa_pt)/sizeof(PHSD_Raa_pt[0]);

	TGraph *gr_PHSD_Raa=new TGraph(nPHSD_Raa,PHSD_Raa_pt,PHSD_Raa_val);
	gr_PHSD_Raa->SetName("gr_PHSD_Raa");
	gr_PHSD_Raa->Write();


	double PHSD_DsD0_pp_pt[]={1,3,5,7,9,12,16};
	int nPHSD_DsD0=sizeof(PHSD_DsD0_pp_pt)/sizeof(PHSD_DsD0_pp_pt[0]);
	double PHSD_DsD0_pp_val[]={
0.173727832,	
0.190932735,	
0.200997246,	
0.203403323,	
0.212203434,	
0.218798775,	
0.218537283	
};

	double PHSD_DsD0_PbPb_val[]={
0.245833495,
0.254393963,
0.268494497,
0.265658358,
0.236222019,
0.21348136,
0.227276164
	};

	double PHSD_DsD0_PbPbpp[nPHSD_DsD0];
	for(int i =0; i<nPHSD_DsD0;i++){
		PHSD_DsD0_PbPbpp[i]=PHSD_DsD0_PbPb_val[i]/PHSD_DsD0_pp_val[i];
	}


	TGraph *gr_PHSD_DsD0_pp=new TGraph(nPHSD_DsD0,PHSD_DsD0_pp_pt,PHSD_DsD0_pp_val);
	gr_PHSD_DsD0_pp->SetName("gr_PHSD_DsD0_pp");
	gr_PHSD_DsD0_pp->Write();

	TGraph *gr_PHSD_DsD0_PbPb=new TGraph(nPHSD_DsD0,PHSD_DsD0_pp_pt,PHSD_DsD0_PbPb_val);
	gr_PHSD_DsD0_PbPb->SetName("gr_PHSD_DsD0_PbPb");
	gr_PHSD_DsD0_PbPb->Write();

	TGraph *gr_PHSD_DsD0_PbPbpp=new TGraph(nPHSD_DsD0,PHSD_DsD0_pp_pt,PHSD_DsD0_PbPbpp);
	gr_PHSD_DsD0_PbPbpp->SetName("gr_PHSD_DsD0_PbPbpp");
	gr_PHSD_DsD0_PbPbpp->Write();


	return;

}
	// return 0;
/*
	TFile *Alice_result=TFile::Open("./output/Alice_result.root","recreate");

	const int nbin_DsD0_0to10=4;
	double bins_pt_0to10[nbin_DsD0_0to10]={5,7,10,14};
	double bins_pt_0to10_err[nbin_DsD0_0to10]={0.5,0.5,0.5,0.5};
	double DsD0_0to10[nbin_DsD0_0to10]={0.38,0.417,0.329,0.336};
	double DsD0_0to10_Errhi[nbin_DsD0_0to10]={0.468,0.498,0.388,0.402};
	double DsD0_0to10_Errlo[nbin_DsD0_0to10]={0.281,0.321,0.259,0.251};

	for(int i=0; i<nbin_DsD0_0to10;i++){
		DsD0_0to10_Errhi[i]=DsD0_0to10_Errhi[i]-DsD0_0to10[i];
		DsD0_0to10_Errlo[i]=DsD0_0to10[i]-DsD0_0to10_Errlo[i];
	}

	TGraphAsymmErrors *gr_DsD0_0to10=new TGraphAsymmErrors(nbin_DsD0_0to10,bins_pt_0to10,DsD0_0to10,bins_pt_0to10_err,bins_pt_0to10_err,DsD0_0to10_Errlo,DsD0_0to10_Errhi);

	const int nbin_DsD0_pp=4;
	double bins_pt_pp[nbin_DsD0_pp]={3,5,7,9};
	double bins_pt_pp_err[nbin_DsD0_pp]={0.5,0.5,0.5,0.5};
	double DsD0_pp[nbin_DsD0_pp]={0.181,0.211,0.226,0.240};
	double DsD0_pp_Errhi[nbin_DsD0_pp]={0.204,0.240,0.251,0.266};
	double DsD0_pp_Errlo[nbin_DsD0_pp]={0.156,0.185,0.196,0.207};

	for(int i=0; i<nbin_DsD0_pp;i++){
		DsD0_pp_Errhi[i]=DsD0_pp_Errhi[i]-DsD0_pp[i];
		DsD0_pp_Errlo[i]=DsD0_pp[i]-DsD0_pp_Errlo[i];
	}

	TGraphAsymmErrors *gr_DsD0_pp=new TGraphAsymmErrors(nbin_DsD0_pp,bins_pt_pp,DsD0_pp,bins_pt_pp_err,bins_pt_pp_err,DsD0_pp_Errlo,DsD0_pp_Errhi);
	
	TCanvas *c_DsD0=new TCanvas("c_DsD0","c_DsD0",800,600);
	c_DsD0->cd();
	gr_DsD0_0to10->SetFillStyle(3005);
	gr_DsD0_0to10->SetFillColor(6);
	gr_DsD0_0to10->Draw("AP2*");
	gr_DsD0_pp->SetFillStyle(3005);
	gr_DsD0_pp->SetFillColor(4);
	gr_DsD0_pp->Draw("p2same");


	const int nbin_Raa=4;
	double bins_Raa_pt[nbin_Raa]={5,7,10,14};
	double bins_Raa_pt_err[nbin_Raa]={0.5,0.5,0.5,0.5};
	double Raa[nbin_Raa]={0.573,0.390,0.291,0.285};
	double Raa_Errhi[nbin_Raa]={0.724,0.475,0.357,0.376};
	double Raa_Errlo[nbin_Raa]={0.409,0.285,0.220,0.148};

	for(int i=0; i<nbin_Raa; i++){
		Raa_Errhi[i]=Raa_Errhi[i]-Raa[i];
		Raa_Errlo[i]=Raa[i]-Raa_Errlo[i];
	}


	TGraphAsymmErrors *gr_Raa=new TGraphAsymmErrors(nbin_Raa,bins_Raa_pt,Raa,bins_Raa_pt_err,bins_Raa_pt_err,Raa_Errlo,Raa_Errhi);


	TCanvas *c_Raa=new TCanvas("c_Raa","c_Raa",800,600);
	c_Raa->cd();
	gr_Raa->Draw("AP2*");

	Alice_result->cd();
	gr_DsD0_0to10->Write("gr_DsD0_0to10");	
	gr_DsD0_pp->Write("gr_DsD0_pp");
	gr_Raa->Write("gr_Raa");
*/

