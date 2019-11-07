
int phi_f0_MCcompare(){

	gStyle->SetOptStat(0);

	TChain *phich_Ds=new TChain("ntDPhikkpi");
	TChain *phich_Hlt=new TChain("ntHlt");
	TChain *phich_Skim=new TChain("ntSkim");
	TChain *phich_Hi=new TChain("ntHi");
	
  phich_Ds->Add("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_MC/Ds_PbPb_MC_Prompt_phikkpi_pt9p5/*.root");
  phich_Hlt->Add("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_MC/Ds_PbPb_MC_Prompt_phikkpi_pt9p5/*.root");
  phich_Skim->Add("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_MC/Ds_PbPb_MC_Prompt_phikkpi_pt9p5/*.root");
  phich_Hi->Add("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_MC/Ds_PbPb_MC_Prompt_phikkpi_pt9p5/*.root");

	phich_Ds->AddFriend(phich_Hlt);
	phich_Ds->AddFriend(phich_Skim);
	phich_Ds->AddFriend(phich_Hi);

	TChain *f0ch_Ds=new TChain("ntDPhikkpi");
	TChain *f0ch_Hlt=new TChain("ntHlt");
	TChain *f0ch_Skim=new TChain("ntSkim");
	TChain *f0ch_Hi=new TChain("ntHi");
	
  f0ch_Ds->Add("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_MC/Ds_PbPb_MC_Prompt_f0980kkpi_pt9p5/*.root");
  f0ch_Hlt->Add("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_MC/Ds_PbPb_MC_Prompt_f0980kkpi_pt9p5/*.root");
  f0ch_Skim->Add("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_MC/Ds_PbPb_MC_Prompt_f0980kkpi_pt9p5/*.root");
  f0ch_Hi->Add("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_MC/Ds_PbPb_MC_Prompt_f0980kkpi_pt9p5/*.root");

	f0ch_Ds->AddFriend(f0ch_Hlt);
	f0ch_Ds->AddFriend(f0ch_Skim);
	f0ch_Ds->AddFriend(f0ch_Hi);


	double DptHi=20;
	double DptLow=10;
	
	TCut cut_event="pclusterCompatibilityFilter&&pprimaryVertexFilter&&phfCoincFilter3";
	TCut cut_trigger="HLT_HIL1MinimumBiasHF2AND_part1_v1 || HLT_HIL1MinimumBiasHF2AND_part2_v1 || HLT_HIL1MinimumBiasHF2AND_part3_v1 || HLT_HIL1MinimumBiasHF2AND_part4_v1 || HLT_HIL1MinimumBiasHF2AND_part5_v1 || HLT_HIL1MinimumBiasHF2AND_part6_v1 || HLT_HIL1MinimumBiasHF2AND_part7_v1 || HLT_HIL1MinimumBiasHF2AND_part8_v1 || HLT_HIL1MinimumBiasHF2AND_part9_v1 || HLT_HIL1MinimumBiasHF2AND_part10_v1 || HLT_HIL1MinimumBiasHF2AND_part11_v1";
	TCut cut_trk="Dtrk1highPurity && Dtrk2highPurity && Dtrk3highPurity && Dtrk1Pt>1 && Dtrk2Pt >1 && Dtrk3Pt >1 && abs(Dtrk1Eta)<1.5 && abs(Dtrk2Eta)<1.5 && abs(Dtrk3Eta)<1.5 && Dtrk1PtErr/Dtrk1Pt <0.1 && Dtrk2PtErr/Dtrk2Pt <0.1 && Dtrk3PtErr/Dtrk3Pt <0.1";
	TCut cut_mass="Dmass>1.90 && Dmass<2.11 && DtktkResmass > 1 && DtktkResmass < 1.2";

	TCut cut_Dpt=Form("Dpt> %f && Dpt< %f", DptLow, DptHi);

	TCut cut_Dtrue="DgencollisionId==0 && DgenBAncestorpt<=0";

	TCut cutAll= cut_event && cut_trk && cut_Dpt && cut_mass && cut_Dtrue;

	TCut cut_Dphi="DsGen==23333";
	TCut cut_Df0="DsGen==24433";

	cout<<"cutAll = "<<cutAll<<endl;
	
	int nbin_Dchi2cl=20;
	double binlow_Dchi2cl=0;
	double binhigh_Dchi2cl=1;

	// TCut cut_large="Dtrk1PtErr/Dtrk1Pt >0.009 && Dtrk2PtErr/Dtrk2Pt >0.009 && Dtrk3PtErr/Dtrk3Pt >0.009";
	// TCut cut_small="Dtrk1PtErr/Dtrk1Pt <0.009 && Dtrk2PtErr/Dtrk2Pt <0.009 && Dtrk3PtErr/Dtrk3Pt <0.009";


	TH1D *h_Dchi2cl_relpterr_phi= new TH1D("h_Dchi2cl_relpterr_phi","h_Dchi2cl_relpterr_phi",nbin_Dchi2cl,binlow_Dchi2cl,binhigh_Dchi2cl); h_Dchi2cl_relpterr_phi->Sumw2();
	TH1D *h_Dchi2cl_relpterr_f0= new TH1D("h_Dchi2cl_relpterr_f0","h_Dchi2cl_relpterr_f0",nbin_Dchi2cl,binlow_Dchi2cl,binhigh_Dchi2cl); h_Dchi2cl_relpterr_f0->Sumw2();

	phich_Ds->Project("h_Dchi2cl_relpterr_phi","Dchi2cl",(TCut) cutAll && cut_Dphi);
	f0ch_Ds->Project("h_Dchi2cl_relpterr_f0","Dchi2cl",(TCut)cutAll && cut_Df0);

	h_Dchi2cl_relpterr_phi->Scale(1/h_Dchi2cl_relpterr_phi->Integral());
	h_Dchi2cl_relpterr_f0->Scale(1/h_Dchi2cl_relpterr_f0->Integral());

	TCanvas *c_Dchi2cl=new TCanvas("c_Dchi2cl","c_Dchi2cl",800,600);
	c_Dchi2cl->cd();
	h_Dchi2cl_relpterr_phi->SetTitle("");
	h_Dchi2cl_relpterr_phi->GetXaxis()->SetTitle("Vertex Probability");
	h_Dchi2cl_relpterr_phi->Draw("same");
	h_Dchi2cl_relpterr_f0->SetLineColor(2);
	h_Dchi2cl_relpterr_f0->Draw("same");

	TLegend *le_Dchi2cl=new TLegend(0.7,0.2,0.85,0.35);
	le_Dchi2cl->SetBorderSize(0);
	le_Dchi2cl->AddEntry(h_Dchi2cl_relpterr_phi,"#phi channel","l");
	le_Dchi2cl->AddEntry(h_Dchi2cl_relpterr_f0,"f0 channel","l");
	le_Dchi2cl->Draw("same");

	TH1D *h_Dchi2cl_ratio=(TH1D*)h_Dchi2cl_relpterr_phi->Clone("h_Dchi2cl_ratio");
	h_Dchi2cl_ratio->Divide(h_Dchi2cl_relpterr_f0);

	TCanvas *c_Dchi2cl_ratio= new TCanvas("c_Dchi2cl_ratio","c_Dchi2cl_ratio",800,600);
	c_Dchi2cl_ratio->cd();
	h_Dchi2cl_ratio->Draw();


	int nbin_Ddls=20;
	double binlow_Ddls=2.5;
	double binhigh_Ddls=15;

	// TCut cut_large="Dtrk1PtErr/Dtrk1Pt >0.009 && Dtrk2PtErr/Dtrk2Pt >0.009 && Dtrk3PtErr/Dtrk3Pt >0.009";
	// TCut cut_small="Dtrk1PtErr/Dtrk1Pt <0.009 && Dtrk2PtErr/Dtrk2Pt <0.009 && Dtrk3PtErr/Dtrk3Pt <0.009";


	TH1D *h_Ddls_relpterr_phi= new TH1D("h_Ddls_relpterr_phi","h_Ddls_relpterr_phi",nbin_Ddls,binlow_Ddls,binhigh_Ddls); h_Ddls_relpterr_phi->Sumw2();
	TH1D *h_Ddls_relpterr_f0= new TH1D("h_Ddls_relpterr_f0","h_Ddls_relpterr_f0",nbin_Ddls,binlow_Ddls,binhigh_Ddls); h_Ddls_relpterr_f0->Sumw2();

	phich_Ds->Project("h_Ddls_relpterr_phi","DsvpvDistance/DsvpvDisErr",(TCut) cutAll && cut_Dphi);
	f0ch_Ds->Project("h_Ddls_relpterr_f0","DsvpvDistance/DsvpvDisErr",(TCut)cutAll && cut_Df0);

	h_Ddls_relpterr_phi->Scale(1/h_Ddls_relpterr_phi->Integral());
	h_Ddls_relpterr_f0->Scale(1/h_Ddls_relpterr_f0->Integral());

	TCanvas *c_Ddls=new TCanvas("c_Ddls","c_Ddls",800,600);
	c_Ddls->cd();
	h_Ddls_relpterr_phi->SetTitle("");
	h_Ddls_relpterr_phi->GetXaxis()->SetTitle("Decay Length Significance");
	h_Ddls_relpterr_phi->Draw("same");
	h_Ddls_relpterr_f0->SetLineColor(2);
	h_Ddls_relpterr_f0->Draw("same");

	TLegend *le_Ddls=new TLegend(0.7,0.7,0.85,0.85);
	le_Ddls->SetBorderSize(0);
	le_Ddls->AddEntry(h_Ddls_relpterr_phi,"#phi channel","l");
	le_Ddls->AddEntry(h_Ddls_relpterr_f0,"f0 channel","l");
	le_Ddls->Draw("same");


	TH1D *h_Ddls_ratio=(TH1D*)h_Ddls_relpterr_phi->Clone("h_Ddls_ratio");
	h_Ddls_ratio->Divide(h_Ddls_relpterr_f0);

	TCanvas *c_Ddls_ratio= new TCanvas("c_Ddls_ratio","c_Ddls_ratio",800,600);
	c_Ddls_ratio->cd();
	h_Ddls_ratio->Draw();




/*
	TCanvas *c_1=new TCanvas("c_1","c_1",800,600);
	c_1->cd();
	h_Dchi2cl_relpterr_large->Draw("same");
	h_Dchi2cl_relpterr_small->SetLineColor(2);
	h_Dchi2cl_relpterr_small->Draw("same");


	//ddls

	int nbin_Ddls=25;
	double binlow_Ddls=0;
	double binhigh_Ddls=40;


	TH1D *h_Ddls_relpterr_large= new TH1D("h_Ddls_relpterr_large","h_Ddls_relpterr_large",nbin_Ddls,binlow_Ddls,binhigh_Ddls); h_Ddls_relpterr_large->Sumw2();
	TH1D *h_Ddls_relpterr_small= new TH1D("h_Ddls_relpterr_small","h_Ddls_relpterr_small",nbin_Ddls,binlow_Ddls,binhigh_Ddls); h_Ddls_relpterr_small->Sumw2();

	phich_Ds->Project("h_Ddls_relpterr_large","DsvpvDistance/DsvpvDisErr",(TCut)cutAll&&cut_large);
	phich_Ds->Project("h_Ddls_relpterr_small","DsvpvDistance/DsvpvDisErr",(TCut)cutAll&&cut_small);

	h_Ddls_relpterr_large->Scale(1/h_Ddls_relpterr_large->Integral());
	h_Ddls_relpterr_small->Scale(1/h_Ddls_relpterr_small->Integral());

	h_Ddls_relpterr_large->Write();
	h_Ddls_relpterr_small->Write();

	
	TCanvas *c_2=new TCanvas("c_2","c_2",800,600);
	c_2->cd();
	h_Ddls_relpterr_large->Draw("same");
	h_Ddls_relpterr_small->SetLineColor(2);
	h_Ddls_relpterr_small->Draw("same");
*/





return 0;
}
