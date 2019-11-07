#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/InitStyle.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/helpers/plotting.h"

int c_count=0;
TCanvas *cn[500];
TH1F *h_ori[500];
TH1F *h_new[500];


TCanvas *Compare_his(TTree *t_ori,TTree *t_new, TString var, int nbin, double binLow, double binHigh, TCut cut= ""){

	c_count++;

	h_ori[c_count]=new TH1F(Form("h_ori_%i",c_count),Form("%s_ori",var.Data()),nbin,binLow,binHigh);
	h_new[c_count]=new TH1F(Form("h_new_%i",c_count),Form("%s_new",var.Data()),nbin,binLow,binHigh);

	t_ori->Project(Form("h_ori_%i",c_count),var.Data(),cut);
	t_new->Project(Form("h_new_%i",c_count),var.Data(),cut*"weight");

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



void Compare_Bias22_Gen(){

	InitStyle();

	TFile *f_ori=TFile::Open("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_NonPrompt_phikkpi_pt5p7.root");

	TTree *t_ntDs_ori=(TTree*)f_ori->Get("ntDPhikkpi");
	TTree *t_ntGen_ori=(TTree*)f_ori->Get("ntGen");
	TTree *t_ntHlt_ori=(TTree*)f_ori->Get("ntHlt");
	TTree *t_ntSkim_ori=(TTree*)f_ori->Get("ntSkim");
	TTree *t_ntHi_ori=(TTree*)f_ori->Get("ntHi");

	t_ntDs_ori->AddFriend(t_ntGen_ori);
	t_ntDs_ori->AddFriend(t_ntHlt_ori);
	t_ntDs_ori->AddFriend(t_ntSkim_ori);
	t_ntDs_ori->AddFriend(t_ntHi_ori);


	TFile *f_new=TFile::Open("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_NonPrompt22_phikkpi_pt5p7.root");

	TTree *t_ntDs_new=(TTree*)f_new->Get("ntDPhikkpi");
	TTree *t_ntGen_new=(TTree*)f_new->Get("ntGen");
	TTree *t_ntHlt_new=(TTree*)f_new->Get("ntHlt");
	TTree *t_ntSkim_new=(TTree*)f_new->Get("ntSkim");
	TTree *t_ntHi_new=(TTree*)f_new->Get("ntHi");

	t_ntDs_new->AddFriend(t_ntGen_new);
	t_ntDs_new->AddFriend(t_ntHlt_new);
	t_ntDs_new->AddFriend(t_ntSkim_new);
	t_ntDs_new->AddFriend(t_ntHi_new);


	TCanvas *c3=Compare_his(t_ntDs_ori,t_ntDs_new,"pthat",48,2,50);
	TCanvas *c4=Compare_his(t_ntDs_ori,t_ntDs_new,"Gpt",20,6,20);
	TCanvas *c5=Compare_his(t_ntDs_ori,t_ntDs_new,"Gtk1pt",20,0,10);
	TCanvas *c6=Compare_his(t_ntDs_ori,t_ntDs_new,"Gtk2pt",20,0,10);
	TCanvas *c7=Compare_his(t_ntDs_ori,t_ntDs_new,"Gtk3pt",20,0,10);
	TCanvas *c8=Compare_his(t_ntDs_ori,t_ntDs_new,"Dsize",20,0,10);
	TCanvas *c9=Compare_his(t_ntDs_ori,t_ntDs_new,"Dpt",20,6,10,"DsGen==23333");
	TCanvas *c10=Compare_his(t_ntDs_ori,t_ntDs_new,"Dalpha",20,0,0.25,"DsGen==23333");
	TCanvas *c11=Compare_his(t_ntDs_ori,t_ntDs_new,"Dchi2cl",20,0,0.25,"DsGen==23333");
	TCanvas *c12=Compare_his(t_ntDs_ori,t_ntDs_new,"DsvpvDistance/DsvpvDisErr",20,0,8,"DsGen==23333");

	// TCanvas *c4=Compare_his(t_ntDs_ori,t_ntDs_new,"Gpt",48,2,50);



}
