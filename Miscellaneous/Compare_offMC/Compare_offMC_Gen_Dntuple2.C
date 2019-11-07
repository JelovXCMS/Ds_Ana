#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"



int c_count=0;
TCanvas *cn[500];
TH1F *h_ori[500];
TH1F *h_new[500];


TCanvas *Compare_his(TTree *t_ori,TTree *t_new, TString var, int nbin, double binLow, double binHigh, TCut cut= ""){

	c_count++;

	h_ori[c_count]=new TH1F(Form("h_Private_%i",c_count),Form("%s_Private",var.Data()),nbin,binLow,binHigh);
	h_new[c_count]=new TH1F(Form("h_Official_%i",c_count),Form("%s_Official",var.Data()),nbin,binLow,binHigh);

	h_ori[c_count]->SetTitle("Private");
	h_new[c_count]->SetTitle("Official");

	t_ori->Project(Form("h_Private_%i",c_count),var.Data(),cut*"weight");
	t_new->Project(Form("h_Official_%i",c_count),var.Data(),cut*"weight");

	h_ori[c_count]->Scale(1/h_ori[c_count]->Integral());
	h_new[c_count]->Scale(1/h_new[c_count]->Integral());

//	cn[c_count]=new TCanvas(Form("c%i",c_count),  Form("c%i",c_count),800,800);	
//	cn[c_count]->cd();

	cn[c_count]=Draw({h_ori[c_count],h_new[c_count]});
	DrawCompare(h_ori[c_count],h_new[c_count],var.Data());

/*
	h_ori[c_count]->SetTitle(var.Data());
	h_ori[c_count]->SetLineColor(4);
	h_ori[c_count]->Draw();
	h_new[c_count]->SetLineColor(2);
	h_new[c_count]->Draw("SAME");
*/
	return cn[c_count];

}



void Compare_offMC_Gen_Dntuple2(int PNP=1){

	InitStyle();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TFile *f_ori=TFile::Open("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_MC/Dntuple_PbPb_MC_NonPrompt22_f0980kkpi_pt19.root","READ");

	TTree *t_ntDs_ori=(TTree*)f_ori->Get("ntDPhikkpi");
	TTree *t_ntGen_ori=(TTree*)f_ori->Get("ntGen");
	TTree *t_ntHlt_ori=(TTree*)f_ori->Get("ntHlt");
	TTree *t_ntSkim_ori=(TTree*)f_ori->Get("ntSkim");
	TTree *t_ntHi_ori=(TTree*)f_ori->Get("ntHi");
//	TTree *t_ntRoot_ori=(TTree*)f_ori->Get("Dfinder/root");

	t_ntDs_ori->AddFriend(t_ntGen_ori);
	t_ntDs_ori->AddFriend(t_ntHlt_ori);
	t_ntDs_ori->AddFriend(t_ntSkim_ori);
	t_ntDs_ori->AddFriend(t_ntHi_ori);
//	t_ntDs_ori->AddFriend(t_ntRoot_ori);


	TFile *f_new=TFile::Open("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/offPbPb_MC/Ds_NonPrompt_f0kkpi_pt19.root","READ");

	TTree *t_ntDs_new=(TTree*)f_new->Get("ntDPhikkpi");
	TTree *t_ntGen_new=(TTree*)f_new->Get("ntGen");
	TTree *t_ntHlt_new=(TTree*)f_new->Get("ntHlt");
	TTree *t_ntSkim_new=(TTree*)f_new->Get("ntSkim");
	TTree *t_ntHi_new=(TTree*)f_new->Get("ntHi");
//	TTree *t_ntRoot_new=(TTree*)f_new->Get("Dfinder/root");

	t_ntDs_new->AddFriend(t_ntGen_new);
	t_ntDs_new->AddFriend(t_ntHlt_new);
	t_ntDs_new->AddFriend(t_ntSkim_new);
	t_ntDs_new->AddFriend(t_ntHi_new);
//	t_ntDs_new->AddFriend(t_ntRoot_new);


	TCanvas *c3=Compare_his(t_ntDs_ori,t_ntDs_new,"pthat",48,2,50);
	TCanvas *c4=Compare_his(t_ntDs_ori,t_ntDs_new,"Gpt",10,20,40);
	TCanvas *c5=Compare_his(t_ntDs_ori,t_ntDs_new,"Gtk1pt",20,0,20);
	TCanvas *c6=Compare_his(t_ntDs_ori,t_ntDs_new,"Gtk2pt",20,0,20);
	TCanvas *c7=Compare_his(t_ntDs_ori,t_ntDs_new,"Gtk3pt",20,0,20);
	TCanvas *c8=Compare_his(t_ntDs_ori,t_ntDs_new,"Dsize",20,0,20);
	TCanvas *c9=Compare_his(t_ntDs_ori,t_ntDs_new,"Dpt",20,20,40,"DsGen==24433");
	TCanvas *c10=Compare_his(t_ntDs_ori,t_ntDs_new,"Dalpha",20,0,0.25,"DsGen==24433");
	TCanvas *c11=Compare_his(t_ntDs_ori,t_ntDs_new,"Dchi2cl",20,0,0.25,"DsGen==24433");
	TCanvas *c12=Compare_his(t_ntDs_ori,t_ntDs_new,"DsvpvDistance/DsvpvDisErr",20,0,8,"DsGen==24433");
	TCanvas *c13=Compare_his(t_ntDs_ori,t_ntDs_new,"hiBin",100,0,200);

	// TCanvas *c13=Compare_his(t_ntDs_ori,t_ntDs_new,"EvtInfo.PVxE",40,0,0.0003);
	// TCanvas *c14=Compare_his(t_ntDs_ori,t_ntDs_new,"TrackInfo.size",20,0,20);
	// TCanvas *c15=Compare_his(t_ntDs_ori,t_ntDs_new,"TrackInfo.ptErr",20,0,0.08);
	// TCanvas *c16=Compare_his(t_ntDs_ori,t_ntDs_new,"DInfo.tktkRes_mass",20,0.98,1.19);
	// TCanvas *c17=Compare_his(t_ntDs_ori,t_ntDs_new,"DInfo.vtxXErr",20,0,0.0004);

// efficiency compare

	double Dchi2clMin=0.03;
	double DalphaMax=0.12;
	double DdlsMin=2.5;
	double TrkptAcc=0.7;
	double TrketaAcc=1.5;


  TCut cutDTrue="DsGen==24433 && DgencollisionId==0";
  TCut cutDPrompt="DgenBAncestorpt<=0";
  TCut cutDNonPrompt="DgenBAncestorpt>0";
  TCut cutDMass=Form("TMath::Abs(Dmass-%f)<%f",DsMass,DsMassRange);

  TCut cutDchi2cl=Form("Dchi2cl > %f", Dchi2clMin);
  TCut cutDalpha=Form("Dalpha <%f", DalphaMax);
  TCut cutDdls=Form("DsvpvDistance/DsvpvDisErr > %f", DdlsMin);

  TCut cutDNorm = cutDTrue && cutDMass && cutDchi2cl && cutDalpha && cutDdls;

  TCut cutGenTrue="GSignalType==2 && GcollisionId==0 && TMath::Abs(Gy)<1";
  TCut cutGenPrompt="GBAncestorpt<=0";
  TCut cutGenNonPrompt="GBAncestorpt>0";
//  TCut cutGenMass=Form("TMath::Abs(Gmass-%f)<%f",DsMass,DsMassRange); // not stored
  TCut cutGenAcc=Form("Gtk2pt >%f && GRestk1pt>%f && GRestk2pt>%f && TMath::Abs(Gtk2eta)<%f && TMath::Abs(GRestk1eta)<%f && TMath::Abs(GRestk2eta)<%f  ",TrkptAcc,TrkptAcc,TrkptAcc,TrketaAcc,TrketaAcc,TrketaAcc);  // must careful using tr2, Restrk1 & Restk2
	TCut cuthiBin="";

	TCut cutDPNP=cutDPrompt;
	TCut cutGenPNP=cutGenPrompt;
	if(PNP==1){ // nonprompt
	 cutDPNP=cutDNonPrompt;
	 cutGenPNP=cutGenNonPrompt;
	}

	double DptLow=20;
	double DptHigh=40;

	int DptBins=5;

  TCut cutDPt=Form("Dpt> %f && Dpt< %f",DptLow,DptHigh);
  TCut cutGenPt=Form("Gpt> %f && Gpt< %f",DptLow,DptHigh);

	
  TH1D *h_GenAlltemp_ori = new TH1D("h_GenAlltemp_ori","h_GenAlltemp_ori",DptBins,DptLow,DptHigh); h_GenAlltemp_ori->Sumw2();
  TH1D *h_GenAcctemp_ori = new TH1D("h_GenAcctemp_ori","h_GenAcctemp_ori",DptBins,DptLow,DptHigh); h_GenAcctemp_ori->Sumw2();
  TH1D *h_RecoNormtemp_ori = new TH1D("h_RecoNormtemp_ori","h_RecoNormtemp_ori",DptBins,DptLow,DptHigh); h_RecoNormtemp_ori->Sumw2();

  t_ntDs_ori->Project("h_GenAlltemp_ori","Gpt",cutGenTrue && cutGenPt && cutGenPNP && cutGenPt && cuthiBin);
  t_ntDs_ori->Project("h_GenAcctemp_ori","Gpt",cutGenTrue && cutGenAcc && cutGenPt && cutGenPNP && cutGenPt && cuthiBin);
  t_ntDs_ori->Project("h_RecoNormtemp_ori","Dpt",cutDNorm && cutDPt && cutDPNP && cutDPt && cuthiBin);

  TH1D *h_GenAccEfftemp_ori=(TH1D*)h_GenAcctemp_ori->Clone("h_GenAccEfftemp_ori");  h_GenAccEfftemp_ori->Sumw2();
  TH1D *h_RecoNormEfftemp_ori=(TH1D*)h_RecoNormtemp_ori->Clone("h_RecoNormEfftemp_ori"); h_RecoNormEfftemp_ori->Sumw2();

  h_GenAccEfftemp_ori->Divide(h_GenAccEfftemp_ori,h_GenAlltemp_ori,1,1,"B");
  h_RecoNormEfftemp_ori->Divide(h_RecoNormEfftemp_ori,h_GenAlltemp_ori,1,1,"B");

	
  TH1D *h_GenAlltemp_new = new TH1D("h_GenAlltemp_new","h_GenAlltemp_new",DptBins,DptLow,DptHigh); h_GenAlltemp_new->Sumw2();
  TH1D *h_GenAcctemp_new = new TH1D("h_GenAcctemp_new","h_GenAcctemp_new",DptBins,DptLow,DptHigh); h_GenAcctemp_new->Sumw2();
  TH1D *h_RecoNormtemp_new = new TH1D("h_RecoNormtemp_new","h_RecoNormtemp_new",DptBins,DptLow,DptHigh); h_RecoNormtemp_new->Sumw2();

  t_ntDs_new->Project("h_GenAlltemp_new","Gpt",(cutGenTrue && cutGenPt && cutGenPNP && cutGenPt && cuthiBin)*"weight");
  t_ntDs_new->Project("h_GenAcctemp_new","Gpt",(cutGenTrue && cutGenAcc && cutGenPt && cutGenPNP && cutGenPt && cuthiBin)*"weight");
  t_ntDs_new->Project("h_RecoNormtemp_new","Dpt",(cutDNorm && cutDPt && cutDPNP && cutDPt && cuthiBin)*"weight");

  TH1D *h_GenAccEfftemp_new=(TH1D*)h_GenAcctemp_new->Clone("h_GenAccEfftemp_new");  h_GenAccEfftemp_new->Sumw2();
  TH1D *h_RecoNormEfftemp_new=(TH1D*)h_RecoNormtemp_new->Clone("h_RecoNormEfftemp_new"); h_RecoNormEfftemp_new->Sumw2();

  h_GenAccEfftemp_new->Divide(h_GenAccEfftemp_new,h_GenAlltemp_new,1,1,"B");
  h_RecoNormEfftemp_new->Divide(h_RecoNormEfftemp_new,h_GenAlltemp_new,1,1,"B");


	Draw({h_GenAccEfftemp_ori,h_GenAccEfftemp_new});
	DrawCompare(h_GenAccEfftemp_ori,h_GenAccEfftemp_new);

	h_RecoNormEfftemp_ori->SetTitle("Private");
	h_RecoNormEfftemp_new->SetTitle("Official");

	Draw({h_RecoNormEfftemp_ori,h_RecoNormEfftemp_new});
	DrawCompare(h_RecoNormEfftemp_ori,h_RecoNormEfftemp_new,"Efficiency");


}
