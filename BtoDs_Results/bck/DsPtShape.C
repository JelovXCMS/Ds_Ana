/// abndon, directly calculate in DtoDs

#include "../include/DsMinTreeLoad.h"
#include "../include/parameters.h"
#include "../include/uti.h"

#include <TH1.h>
#include <TH1D.h>
#include <TFile.h>


void DsPtShape(){



	double bins_pt[]={0,2,3,4,6,8,10,20,30,999}; // using b->D binning
	const int nbin_pt=sizeof(bins_pt)/sizeof(bins_pt[0])-1;

	TFile *f_out=new TFile("NonPrompt_DsGpt_phikkpi_pp.root","RECREATE");
	
	TH1D *h_DsGpt=new TH1D("h_DsGpt","h_DsGpt",nbin_pt,bins_pt); h_DsGpt->Sumw2();

	TH1D *h_DGpt=new TH1D("h_DGpt","h_DGpt",nbin_pt,bins_pt); h_DGpt->Sumw2();

	TChain *ch_DsMB= new TChain("ntGen");
	ch_DsMB->Add("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_MC/DsMinTree_pp_MC_Ds_NonPrompt_phikkpi_pt0.root");
	ch_DsMB->Project("h_DsGpt","Gpt","(GSignalType==1 &&GcollisionId==0 && abs(Gy)<1 && GBAncestorpt>0)");

	TChain *ch_DMB= new TChain("Dfinder/ntGen");
	ch_DMB->Add("/mnt/hadoop/store/user/hqiu/DfinderProduction_newNewMC/NonPrD0_pThat-0_pT-0_pp_5p02-Pythia8/crab_Dfinder_ppMC_nonPromptD0_pt0_1/170426_204331/0000/*.root");
	ch_DMB->Project("h_DGpt","Gpt","( (GisSignal==1 || GisSignal==2) &&GcollisionId==0 && abs(Gy)<1 && GBAncestorpt>0)");

	h_DsGpt->Scale(1/h_DsGpt->Integral());
	h_DGpt->Scale(1/h_DGpt->Integral());

	TH1D *h_D0toDs_pythiaWeight=(TH1D*)h_DsGpt->Clone("h_D0toDs_pythiaWeight");
	h_D0toDs_pythiaWeight->SetTitle("h_D0toDs_pythiaWeight");
	h_D0toDs_pythiaWeight->Divide(h_DGpt);

	f_out->cd();
	h_DsGpt->Write("",TObject::kOverwrite);
	h_DGpt->Write("",TObject::kOverwrite);
	h_D0toDs_pythiaWeight->Write("",TObject::kOverwrite);

	gStyle->SetOptStat(0);
	TCanvas *c_PythiaWeight=new TCanvas("c_PythiaWeight","c_PythiaWeight",800,1600);
	c_PythiaWeight->Divide(1,2);
	c_PythiaWeight->cd(1);
	gPad->SetLogy();
	h_DsGpt->GetXaxis()->SetRangeUser(0,30);
	h_DsGpt->SetLineColor(2);
	h_DsGpt->Draw();
	h_DGpt->SetLineColor(4);
	h_DGpt->Draw("SAME");

	TLegend *le_PythiaWeight= new TLegend(0.55,0.6,0.85,0.85);
	le_PythiaWeight->SetBorderSize(0);
	le_PythiaWeight->AddEntry(h_DsGpt,"Ds pT","l");
	le_PythiaWeight->AddEntry(h_DGpt,"D pT","l");
	le_PythiaWeight->Draw("SAME");

	c_PythiaWeight->cd(2);
	h_D0toDs_pythiaWeight->GetXaxis()->SetRangeUser(0,30);
	h_D0toDs_pythiaWeight->Draw();

	c_PythiaWeight->SaveAs("DtoDs_PythiaWeight.pdf");

}
