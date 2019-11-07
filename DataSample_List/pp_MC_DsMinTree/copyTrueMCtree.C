
void  copyTrueMCtree(){

	TFile *fin=TFile::Open("DsMinTree_pp_MC_Ds_phikkpi_pt4.root","read");
	TTree *t_in=(TTree*)fin->Get("ntDs");
	Long64_t nentries=t_in->GetEntries();
	int DsGen;
	t_in->SetBranchAddress("DsGen",&DsGen);

	TFile *fout=new TFile("DsMinTree_pp_MC_Ds_phikkpi_pt4_true.root","RECREATE");
	TTree *t_out=t_in->CloneTree(0);

	for(Long64_t i=0; i<nentries; i++){
		t_in->GetEntry(i);
		if(DsGen==23333) t_out->Fill();	
	}

	t_out->Print();
	t_out->AutoSave();

}
