
#ifndef DsMinTreeLoad_h
#define DsMinTreeLoad_h

#define MAX_XB        20000

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


// follow /home/peng43/work/Project/GspDiBjet2017/GspDiBjet2017/include/branch_class

class DsMinTreeLoad {

	public :
		TTree          *ntDs;   //!pointer to the analyzed TTree or TChain
	  TTree          *ntGen;

		// common
    Int_t           hiBin;
		Float_t         weight;
	  Float_t         pthat;
		Float_t         vx;
		Float_t         vy;
		Float_t         vz;
		
		Float_t         Ncoll;
		Float_t         Npart;

    Int_t           GenT_hiBin;
		Float_t         GenT_weight;
	  Float_t         GenT_pthat;
		Float_t         GenT_vx;
		Float_t         GenT_vy;
		Float_t         GenT_vz;

		Float_t         GenT_Ncoll;
		Float_t         GenT_Npart;


	// Reco Tree
		Int_t           DsGen;
  	Float_t         Dgenpt;
  	Int_t           DgencollisionId;
  	Float_t         DgenBAncestorpt;
  	Int_t           DgenBAncestorpdgId;
  	Int_t           DgenfromgenPV;
	
		Float_t         DgendecayvtxXtoPv;
		Float_t         DgendecayvtxYtoPv;
		Float_t         DgendecayvtxZtoPv;

		Float_t         Dmass;
		Float_t         Dpt;
		Float_t         Dchi2cl;
		Float_t         Dalpha;
		Float_t         Dalpha_BS_2D;
		Float_t         Ddls;
		
		Float_t         Ddls_BestScale;

		Float_t         Ddl;
		Float_t         DdlErr;
		Float_t         DlxyBS;
		Float_t         DlxyBSErr;
		Float_t         DlxyBSs;
		Float_t         Ddca;
		Float_t         Dtrk1Pt;
		Float_t         Dtrk2Pt;
		Float_t         Dtrk3Pt;
		Float_t         DtktkResmass;
		Float_t         DtktkRes_chi2cl;
	  Float_t         DtktkRes_svpvDistanceToSV;
    Float_t         DtktkRes_dlsToSV;

		Float_t         DvtxX;
		Float_t         DvtxY;
		Float_t         DvtxZ;
		Float_t         DvtxXErr;
		Float_t         DvtxYErr;
		Float_t         DvtxZErr;
		Float_t         Dd0Err;
		Float_t         DdxyzErr;	

		Float_t         PVxE;
		Float_t         PVyE;
		Float_t         PVzE;

		


	// Gen tree variables

Float_t         GPVx;
Float_t         GPVy;
Float_t         GPVz;
Int_t           Gsize;
Float_t         Gy[MAX_XB];   //[Gsize]
Float_t         Geta[MAX_XB];   //[Gsize]
Float_t         Gphi[MAX_XB];   //[Gsize]
Float_t         Gpt[MAX_XB];   //[Gsize]
Int_t           GpdgId[MAX_XB];   //[Gsize]
Int_t           GcollisionId[MAX_XB];   //[Gsize]
Int_t           GisSignal[MAX_XB];   //[Gsize]
Int_t           GSignalType[MAX_XB];   //[Gsize]
Float_t         GBAncestorpt[MAX_XB];   //[Gsize]
Int_t           GBAncestorpdgId[MAX_XB];   //[Gsize]
Float_t         Gtk1pt[MAX_XB];   //[Gsize]
Float_t         Gtk1eta[MAX_XB];   //[Gsize]
Float_t         Gtk1y[MAX_XB];   //[Gsize]
Float_t         Gtk1phi[MAX_XB];   //[Gsize]
Float_t         Gtk2pt[MAX_XB];   //[Gsize]
Float_t         Gtk2eta[MAX_XB];   //[Gsize]
Float_t         Gtk2y[MAX_XB];   //[Gsize]
Float_t         Gtk2phi[MAX_XB];   //[Gsize]

Float_t         GRestk1pt[MAX_XB];
Float_t         GRestk1eta[MAX_XB];
Float_t         GRestk1y[MAX_XB];
Float_t         GRestk1phi[MAX_XB];
Float_t         GRestk2pt[MAX_XB];
Float_t         GRestk2eta[MAX_XB];
Float_t         GRestk2y [MAX_XB];
Float_t         GRestk2phi[MAX_XB];




		DsMinTreeLoad()=default;
		virtual void SetBranch(TTree *tree, Bool_t REAL, Int_t isPbPb, Int_t readScale);
		virtual void SetGenBranch(TTree *tree, Int_t isPbPb);	


};

#endif

void DsMinTreeLoad::SetBranch(TTree *tree, Bool_t REAL=true, Int_t isPbPb=0, Int_t readScale=1)
{
	if (!tree) return;
	ntDs = tree;

	if(readScale){
		ntDs->SetBranchAddress("Ddls_BestScale"           ,  &Ddls_BestScale);
	}

	ntDs->SetBranchAddress("Dmass"          ,  &Dmass);
	ntDs->SetBranchAddress("Dpt"            ,  &Dpt);
	ntDs->SetBranchAddress("Dchi2cl"        ,  &Dchi2cl);
	ntDs->SetBranchAddress("Dalpha"         ,  &Dalpha);
	ntDs->SetBranchAddress("Dalpha_BS_2D"         ,  &Dalpha_BS_2D);
	ntDs->SetBranchAddress("Ddls"           ,  &Ddls);
	ntDs->SetBranchAddress("Ddl"            ,  &Ddl);
	ntDs->SetBranchAddress("DdlErr"         ,  &DdlErr);
	ntDs->SetBranchAddress("Ddl"            ,  &Ddl);
	ntDs->SetBranchAddress("DlxyBS"           ,  &DlxyBS);
	ntDs->SetBranchAddress("DlxyBSErr"           ,  &DlxyBSErr);
	ntDs->SetBranchAddress("DlxyBSs"           ,  &DlxyBSs);
	ntDs->SetBranchAddress("Dtrk1Pt"        ,  &Dtrk1Pt);
	ntDs->SetBranchAddress("Dtrk2Pt"        ,  &Dtrk2Pt);
	ntDs->SetBranchAddress("Dtrk3Pt"        ,  &Dtrk3Pt);
	ntDs->SetBranchAddress("DtktkResmass"   ,  &DtktkResmass);
	ntDs->SetBranchAddress("DtktkRes_chi2cl",  &DtktkRes_chi2cl);

  ntDs->SetBranchAddress("DtktkRes_svpvDistanceToSV",  &DtktkRes_svpvDistanceToSV);
  ntDs->SetBranchAddress("DtktkRes_dlsToSV",    &DtktkRes_dlsToSV  );

  ntDs->SetBranchAddress("DsGen"          ,  &DsGen);
  ntDs->SetBranchAddress("DgencollisionId",&DgencollisionId);
  ntDs->SetBranchAddress("DgenBAncestorpt",&DgenBAncestorpt);


  ntDs->SetBranchAddress("hiBin",&hiBin);
  ntDs->SetBranchAddress("vx",&vx);
  ntDs->SetBranchAddress("vy",&vy);
  ntDs->SetBranchAddress("vz",&vz);

  // ntDs->SetBranchAddress("Dd0Err",&Dd0Err);
  ntDs->SetBranchAddress("DdxyzErr",&DdxyzErr);

  ntDs->SetBranchAddress("DvtxX",&DvtxX);
  ntDs->SetBranchAddress("DvtxY",&DvtxY);
  ntDs->SetBranchAddress("DvtxZ",&DvtxZ);
  ntDs->SetBranchAddress("DvtxXErr",&DvtxXErr);
  ntDs->SetBranchAddress("DvtxYErr",&DvtxYErr);
  ntDs->SetBranchAddress("DvtxZErr",&DvtxZErr);

  ntDs->SetBranchAddress("PVxE",&PVxE);
  ntDs->SetBranchAddress("PVyE",&PVyE);
  ntDs->SetBranchAddress("PVzE",&PVzE);



	if(!REAL){
    ntDs->SetBranchAddress("pthat",&pthat);
	  ntDs->SetBranchAddress("weight",&weight);
    ntDs->SetBranchAddress("Dgenpt",&Dgenpt);
    ntDs->SetBranchAddress("DgenBAncestorpdgId",&DgenBAncestorpdgId);
    ntDs->SetBranchAddress("DgenfromgenPV",&DgenfromgenPV);
 
	  ntDs->SetBranchAddress("DgendecayvtxXtoPv",&DgendecayvtxXtoPv);
	  ntDs->SetBranchAddress("DgendecayvtxYtoPv",&DgendecayvtxYtoPv);
	  ntDs->SetBranchAddress("DgendecayvtxZtoPv",&DgendecayvtxZtoPv);


	  if(isPbPb){
		  ntDs->SetBranchAddress("Ncoll",&Ncoll);
	  	ntDs->SetBranchAddress("Npart",&Npart);
	  }
	//	ntDs->SetBranchAddress("DsGen"          ,  &DsGen);
  // ntDs->SetBranchAddress("DgencollisionId",&DgencollisionId);
  // ntDs->SetBranchAddress("DgenBAncestorpt",&DgenBAncestorpt);
	}

}

void DsMinTreeLoad::SetGenBranch(TTree *tree, Int_t isPbPb=0)
{
	if(!tree) return;
	ntGen = tree;

  ntGen->SetBranchAddress("GPVx", &GPVx);
  ntGen->SetBranchAddress("GPVy", &GPVy);
  ntGen->SetBranchAddress("GPVz", &GPVz);
  ntGen->SetBranchAddress("Gsize", &Gsize);
  ntGen->SetBranchAddress("Gy", Gy);
  ntGen->SetBranchAddress("Geta", Geta);
  ntGen->SetBranchAddress("Gphi", Gphi);
  ntGen->SetBranchAddress("Gpt", Gpt);
  ntGen->SetBranchAddress("GpdgId", GpdgId);
  ntGen->SetBranchAddress("GcollisionId", GcollisionId);
  ntGen->SetBranchAddress("GisSignal", GisSignal);
  ntGen->SetBranchAddress("GSignalType", GSignalType);
  ntGen->SetBranchAddress("GBAncestorpt", GBAncestorpt);
  ntGen->SetBranchAddress("GBAncestorpdgId", GBAncestorpdgId);
  ntGen->SetBranchAddress("Gtk1pt", Gtk1pt);
  ntGen->SetBranchAddress("Gtk1eta", Gtk1eta);
  // ntGen->SetBranchAddress("Gtk1y", Gtk1y);
  // ntGen->SetBranchAddress("Gtk1phi", Gtk1phi);
  ntGen->SetBranchAddress("Gtk2pt", Gtk2pt);
  ntGen->SetBranchAddress("Gtk2eta", Gtk2eta);
  // ntGen->SetBranchAddress("Gtk2y", Gtk2y);
  // ntGen->SetBranchAddress("Gtk2phi", Gtk2phi);

  ntGen->SetBranchAddress("GRestk1pt",    GRestk1pt    );
  ntGen->SetBranchAddress("GRestk1eta",   GRestk1eta);
  // ntGen->SetBranchAddress("GRestk1y",     GRestk1y);
  // ntGen->SetBranchAddress("GRestk1phi",   GRestk1phi);
  ntGen->SetBranchAddress("GRestk2pt",    GRestk2pt);
  ntGen->SetBranchAddress("GRestk2eta",   GRestk2eta);
  // ntGen->SetBranchAddress("GRestk2y",     GRestk2y);
  // ntGen->SetBranchAddress("GRestk2phi",   GRestk2phi);

	ntGen->SetBranchAddress("hiBin",&GenT_hiBin);
	ntGen->SetBranchAddress("weight",&GenT_weight);
	ntGen->SetBranchAddress("vx",&GenT_vx);
	ntGen->SetBranchAddress("vy",&GenT_vy);
	ntGen->SetBranchAddress("vz",&GenT_vz);

  // ntGen->SetBranchAddress("pthat",&GenT_pthat);
	if(isPbPb){
		ntGen->SetBranchAddress("Ncoll",&GenT_Ncoll);
		ntGen->SetBranchAddress("Npart",&GenT_Npart);
	}
	

}




