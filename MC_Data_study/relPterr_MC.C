
int relPterr_MC(){

	gStyle->SetOptStat(0);

	TChain *ch_Ds=new TChain("ntDPhikkpi");
	TChain *ch_Hlt=new TChain("ntHlt");
	TChain *ch_Skim=new TChain("ntSkim");
	TChain *ch_Hi=new TChain("ntHi");
	
  ch_Ds->Add("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_MC/Ds_PbPb_MC_Prompt_phikkpi_pt9p5/*.root");
  ch_Hlt->Add("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_MC/Ds_PbPb_MC_Prompt_phikkpi_pt9p5/*.root");
  ch_Skim->Add("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_MC/Ds_PbPb_MC_Prompt_phikkpi_pt9p5/*.root");
  ch_Hi->Add("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_MC/Ds_PbPb_MC_Prompt_phikkpi_pt9p5/*.root");

	ch_Ds->AddFriend(ch_Hlt);
	ch_Ds->AddFriend(ch_Skim);
	ch_Ds->AddFriend(ch_Hi);

	double DptHi=20;
	double DptLow=10;
	
	TCut cut_event="pclusterCompatibilityFilter&&pprimaryVertexFilter&&phfCoincFilter3";
	TCut cut_trigger="HLT_HIL1MinimumBiasHF2AND_part1_v1 || HLT_HIL1MinimumBiasHF2AND_part2_v1 || HLT_HIL1MinimumBiasHF2AND_part3_v1 || HLT_HIL1MinimumBiasHF2AND_part4_v1 || HLT_HIL1MinimumBiasHF2AND_part5_v1 || HLT_HIL1MinimumBiasHF2AND_part6_v1 || HLT_HIL1MinimumBiasHF2AND_part7_v1 || HLT_HIL1MinimumBiasHF2AND_part8_v1 || HLT_HIL1MinimumBiasHF2AND_part9_v1 || HLT_HIL1MinimumBiasHF2AND_part10_v1 || HLT_HIL1MinimumBiasHF2AND_part11_v1";
	TCut cut_trk="Dtrk1highPurity && Dtrk2highPurity && Dtrk3highPurity && Dtrk1Pt>1 && Dtrk2Pt >1 && Dtrk3Pt >1 && abs(Dtrk1Eta)<1.5 && abs(Dtrk2Eta)<1.5 && abs(Dtrk3Eta)<1.5 && Dtrk1PtErr/Dtrk1Pt <0.1 && Dtrk2PtErr/Dtrk2Pt <0.1 && Dtrk3PtErr/Dtrk3Pt <0.1";
	TCut cut_mass="Dmass>1.90 && Dmass<2.11 && DtktkResmass > 1 && DtktkResmass < 1.2";

	TCut cut_Dpt=Form("Dpt> %f && Dpt< %f", DptLow, DptHi);

	TCut cut_Dtrue="DsGen==23333 && DgencollisionId==0 && DgenBAncestorpt<=0";

	TCut cutAll= cut_event && cut_trk && cut_Dpt && cut_mass && cut_Dtrue;


	cout<<"cutAll = "<<cutAll<<endl;
	
	int nbin_relpterr=100;
	double binlow_relpterr=0.005;
	double binhigh_relpterr=0.02;

	TFile *fout=new TFile("MC_trkrelpterr.root","recreate");

	fout->cd();

	TH1D *h_trk1relpterr= new TH1D("h_trk1relpterr","h_trk1relpterr",nbin_relpterr,binlow_relpterr,binhigh_relpterr); h_trk1relpterr->Sumw2();
	TH1D *h_trk2relpterr= new TH1D("h_trk2relpterr","h_trk2relpterr",nbin_relpterr,binlow_relpterr,binhigh_relpterr); h_trk2relpterr->Sumw2();
	TH1D *h_trk3relpterr= new TH1D("h_trk3relpterr","h_trk3relpterr",nbin_relpterr,binlow_relpterr,binhigh_relpterr); h_trk3relpterr->Sumw2();
	TH1D *h_trksrelpterr= new TH1D("h_trksrelpterr","h_trksrelpterr",nbin_relpterr,binlow_relpterr,binhigh_relpterr); h_trksrelpterr->Sumw2();


	// ch_Ds->Draw("Dtrk1PtErr/Dtrk1Pt",cutAll);
	// ch_Ds->Draw("Dtrk2PtErr/Dtrk2Pt",cutAll);
	// ch_Ds->Draw("Dtrk3PtErr/Dtrk3Pt",cutAll);

	ch_Ds->Project("h_trk1relpterr","Dtrk1PtErr/Dtrk1Pt",cutAll);
	ch_Ds->Project("h_trk2relpterr","Dtrk2PtErr/Dtrk2Pt",cutAll);
	ch_Ds->Project("h_trk3relpterr","Dtrk3PtErr/Dtrk3Pt",cutAll);
	
	h_trksrelpterr->Add(h_trk1relpterr);
	h_trksrelpterr->Add(h_trk2relpterr);
	h_trksrelpterr->Add(h_trk3relpterr);

	h_trk1relpterr->Write();
	h_trk2relpterr->Write();
	h_trk3relpterr->Write();
	h_trksrelpterr->Write();


	// compare Dchi2cl/ Ddls for different relTrkpterr distribtuion.

	int nbin_Dchi2cl=25;
	double binlow_Dchi2cl=0;
	double binhigh_Dchi2cl=1;

	// TCut cut_large="Dtrk1PtErr/Dtrk1Pt >0.009 && Dtrk2PtErr/Dtrk2Pt >0.009 && Dtrk3PtErr/Dtrk3Pt >0.009";
	// TCut cut_small="Dtrk1PtErr/Dtrk1Pt <0.009 && Dtrk2PtErr/Dtrk2Pt <0.009 && Dtrk3PtErr/Dtrk3Pt <0.009";


	TCut cut_large="Dtrk1PtErr/Dtrk1Pt >0.009 && Dtrk3PtErr/Dtrk3Pt >0.009 ";
	TCut cut_small="Dtrk1PtErr/Dtrk1Pt <0.009 && Dtrk3PtErr/Dtrk3Pt <0.009 ";

	TH1D *h_Dchi2cl_relpterr_large= new TH1D("h_Dchi2cl_relpterr_large","h_Dchi2cl_relpterr_large",nbin_Dchi2cl,binlow_Dchi2cl,binhigh_Dchi2cl); h_Dchi2cl_relpterr_large->Sumw2();
	TH1D *h_Dchi2cl_relpterr_small= new TH1D("h_Dchi2cl_relpterr_small","h_Dchi2cl_relpterr_small",nbin_Dchi2cl,binlow_Dchi2cl,binhigh_Dchi2cl); h_Dchi2cl_relpterr_small->Sumw2();

	ch_Ds->Project("h_Dchi2cl_relpterr_large","Dchi2cl",(TCut)cutAll&&cut_large);
	ch_Ds->Project("h_Dchi2cl_relpterr_small","Dchi2cl",(TCut)cutAll&&cut_small);

	h_Dchi2cl_relpterr_large->Scale(1/h_Dchi2cl_relpterr_large->Integral());
	h_Dchi2cl_relpterr_small->Scale(1/h_Dchi2cl_relpterr_small->Integral());

	h_Dchi2cl_relpterr_large->Write();
	h_Dchi2cl_relpterr_small->Write();

	
	TCanvas *c_1=new TCanvas("c_1","c_1",800,600);
	c_1->cd();
	h_Dchi2cl_relpterr_large->SetTitle("");
	h_Dchi2cl_relpterr_large->GetXaxis()->SetTitle("Vertex Probability");
	h_Dchi2cl_relpterr_large->Draw("same");
	h_Dchi2cl_relpterr_small->SetLineColor(2);
	h_Dchi2cl_relpterr_small->Draw("same");

	TLegend *le_Dchi2cl=new TLegend(0.7,0.7,0.85,0.85);
	le_Dchi2cl->SetBorderSize(0);
	le_Dchi2cl->AddEntry(h_Dchi2cl_relpterr_large,"rel pt err >0.009","l");
	le_Dchi2cl->AddEntry(h_Dchi2cl_relpterr_small,"rel pt err <0.009","l");
	le_Dchi2cl->Draw("same");

	//ddls

	int nbin_Ddls=25;
	double binlow_Ddls=0;
	double binhigh_Ddls=40;


	TH1D *h_Ddls_relpterr_large= new TH1D("h_Ddls_relpterr_large","h_Ddls_relpterr_large",nbin_Ddls,binlow_Ddls,binhigh_Ddls); h_Ddls_relpterr_large->Sumw2();
	TH1D *h_Ddls_relpterr_small= new TH1D("h_Ddls_relpterr_small","h_Ddls_relpterr_small",nbin_Ddls,binlow_Ddls,binhigh_Ddls); h_Ddls_relpterr_small->Sumw2();

	ch_Ds->Project("h_Ddls_relpterr_large","DsvpvDistance/DsvpvDisErr",(TCut)cutAll&&cut_large);
	ch_Ds->Project("h_Ddls_relpterr_small","DsvpvDistance/DsvpvDisErr",(TCut)cutAll&&cut_small);

	h_Ddls_relpterr_large->Scale(1/h_Ddls_relpterr_large->Integral());
	h_Ddls_relpterr_small->Scale(1/h_Ddls_relpterr_small->Integral());

	h_Ddls_relpterr_large->Write();
	h_Ddls_relpterr_small->Write();

	
	TCanvas *c_2=new TCanvas("c_2","c_2",800,600);
	c_2->cd();
	h_Ddls_relpterr_large->SetTitle("");
	h_Ddls_relpterr_large->GetXaxis()->SetTitle("Decay length significance");
	h_Ddls_relpterr_large->Draw("same");
	h_Ddls_relpterr_small->SetLineColor(2);
	h_Ddls_relpterr_small->Draw("same");

	TLegend *le_Ddls=new TLegend(0.7,0.7,0.85,0.85);
	le_Ddls->SetBorderSize(0);
	le_Ddls->AddEntry(h_Ddls_relpterr_large,"rel pt err >0.009","l");
	le_Ddls->AddEntry(h_Ddls_relpterr_small,"rel pt err <0.009","l");
	le_Ddls->Draw("same");




return 0;
}
