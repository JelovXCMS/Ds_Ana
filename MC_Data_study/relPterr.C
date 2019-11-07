
int relPterr(){

	TChain *ch_Ds=new TChain("ntDPhikkpi");
	TChain *ch_Hlt=new TChain("ntHlt");
	TChain *ch_Skim=new TChain("ntSkim");
	TChain *ch_Hi=new TChain("ntHi");
	
  ch_Ds->Add("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_Data/GJ/HIMB3_GJpart3/*5.root");
  ch_Hlt->Add("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_Data/GJ/HIMB3_GJpart3/*5.root");
  ch_Skim->Add("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_Data/GJ/HIMB3_GJpart3/*5.root");
  ch_Hi->Add("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_Data/GJ/HIMB3_GJpart3/*5.root");

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

	TCut cutAll= cut_event && cut_trigger && cut_trk && cut_Dpt && cut_mass;

	cout<<"cutAll = "<<cutAll<<endl;
	
	int nbin_relpterr=100;
	double binlow_relpterr=0.005;
	double binhigh_relpterr=0.02;

	TFile *fout=new TFile("f_trkrelpterr.root","recreate");

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



return 0;
}
