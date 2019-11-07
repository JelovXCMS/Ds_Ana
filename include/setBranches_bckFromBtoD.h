// add three variable for new sample
// add genVtxInfo 
// ntDkpi : DgenPtWeight
// ntHi : DPtSample
// ntGen : GptWeight

//  Bool_t hasDptFilteredWeight=false, Bool_t hasDgenVtx=false 
//Bool_t useDptFilteredWeight=false, Bool_t useDgenVtx=false
// Bool_t isPbPb
// Bool_t isData=false, Bool_t isPbPb=false, Bool_t isPthatWeight=false,  Bool_t hasDptFilteredWeight=false
// Bool_t isData=false, Bool_t isPbPb=false, Bool_t isPthatWeight=false
/*
void SetRecoBranches(TTree* nt, Bool_t hasDptFilteredWeight=false, Bool_t hasDgenVtx=false, Bool_t hasHardQcdWeight=false)
void RecoBranches(TTree* dnt,  Bool_t useDptFilteredWeight=false, Bool_t useDgenVtx=false,  Bool_t useFonllDptWeight=false)
void SetHltBranches(TTree* nt, Bool_t isPbPb)
void HltBranches(TTree* nt, Bool_t isPbPb)
void SetHiBranches(TTree* nt, Bool_t isData=false, Bool_t isPbPb=false, Bool_t isPthatWeight=false,  Bool_t hasDptFilteredWeight=false)
void HiBranches(TTree* nt, Bool_t isData=false, Bool_t isPbPb=false, Bool_t isPthatWeight=false, Bool_t hasDptFilteredWeight=false)
void SetSkimBranches(TTree* nt,  Bool_t isPbPb=false)
void SkimBranches(TTree* nt,  Bool_t isPbPb=false){
void SetGenBranches(TTree* nt, bool hasDptFilteredWeight=false, bool hasDgenVtx=false, Bool_t hasHardQcdWeight=false)
void GenBranches(TTree* nt,  bool useDptFilteredWeight=false, bool useDgenVtx=false, Bool_t useFonllDptWeight=false)
*/



#ifndef _ANA_BFEEDDOWN_SETBRANCHES_H_
#define _ANA_BFEEDDOWN_SETBRANCHES_H_

//#include "uti.h"
#define MAX_XB        20000
#define kMaxEvtPlanes 200

// DntupleBranch

// EvtInfo
Int_t     RunNo;
Int_t     EvtNo;
Int_t     LumiNo;
Int_t     Dsize;
Float_t   PVx;
Float_t   PVy;
Float_t   PVxE;
Float_t   PVyE;
Float_t   PVzE;
Float_t   PVz;
// Dinfo
Float_t   Dmass[MAX_XB];
Float_t   Dpt[MAX_XB];
Float_t   Deta[MAX_XB];
Float_t   Dphi[MAX_XB];
Float_t   Dy[MAX_XB];
Float_t   Dchi2ndf[MAX_XB];
Float_t   Dchi2cl[MAX_XB];
Float_t   Dalpha[MAX_XB];
Float_t   DsvpvDistance[MAX_XB];
Float_t   DsvpvDisErr[MAX_XB];
Float_t   DsvpvDistance_2D[MAX_XB];
Float_t   DsvpvDisErr_2D[MAX_XB];
Float_t   DlxyBS[MAX_XB];
Float_t   DlxyBSErr[MAX_XB];
Float_t   DvtxX[MAX_XB];
Float_t   DvtxY[MAX_XB];
Float_t   DvtxZ[MAX_XB];
Float_t   Dtrk1Pt[MAX_XB];
Float_t   Dtrk2Pt[MAX_XB];
Float_t   Dtrk1Eta[MAX_XB];
Float_t   Dtrk2Eta[MAX_XB];
Float_t   Dtrk1Phi[MAX_XB];
Float_t   Dtrk2Phi[MAX_XB];
Float_t   Dtrk1PtErr[MAX_XB];
Float_t   Dtrk2PtErr[MAX_XB];
Float_t   Dtrk1PixelHit[MAX_XB];
Float_t   Dtrk2PixelHit[MAX_XB];
Float_t   Dtrk1StripHit[MAX_XB];
Float_t   Dtrk2StripHit[MAX_XB];
Float_t   Dtrk1nStripLayer[MAX_XB];
Float_t   Dtrk2nStripLayer[MAX_XB];
Float_t   Dtrk1nPixelLayer[MAX_XB];
Float_t   Dtrk2nPixelLayer[MAX_XB];
Float_t   Dtrk1Chi2ndf[MAX_XB];
Float_t   Dtrk2Chi2ndf[MAX_XB];
Bool_t    Dtrk1highPurity[MAX_XB];
Bool_t    Dtrk2highPurity[MAX_XB];
//Dgen info
Int_t   Dgen[MAX_XB];
Int_t     DgennDa[MAX_XB];
Int_t     DgenIndex[MAX_XB];
Float_t   Dgenpt[MAX_XB];
Float_t   Dgeneta[MAX_XB];
Float_t   Dgenphi[MAX_XB];
Float_t   Dgeny[MAX_XB];
Int_t     DgenBAncestorpdgId[MAX_XB];
Float_t   DgenBAncestorpt[MAX_XB];
Int_t     DgencollisionId[MAX_XB];
Float_t 	DgenPtWeight[MAX_XB];

Float_t   DgenPtPtHatWeight[MAX_XB];

Float_t FonllDptWeight[MAX_XB];


//Dgen Vtx 
  float   DgenprodvtxX[MAX_XB];
  float   DgenprodvtxY[MAX_XB];
  float   DgenprodvtxZ[MAX_XB];
  float   DgendecayvtxX[MAX_XB];
  float   DgendecayvtxY[MAX_XB];
  float   DgendecayvtxZ[MAX_XB];
  int     DgenfromgenPV[MAX_XB];



   TBranch        *b_DsvpvDistance;   //!
   TBranch        *b_DsvpvDisErr;   //!
   TBranch        *b_DsvpvDistance_2D;   //!
   TBranch        *b_DsvpvDisErr_2D;   //!



void SetRecoBranches(TTree* nt, Bool_t hasDptFilteredWeight=false, Bool_t hasDgenVtx=false, Bool_t hasHardQcdWeight=false, Bool_t hasFonllDptWeight=false)
{

	nt->SetBranchAddress("RunNo",&RunNo);
	nt->SetBranchAddress("EvtNo",&EvtNo);
	nt->SetBranchAddress("LumiNo",&LumiNo);
	nt->SetBranchAddress("PVx",&PVx);
	nt->SetBranchAddress("PVy",&PVy);
	nt->SetBranchAddress("PVxE",&PVxE);
	nt->SetBranchAddress("PVyE",&PVyE);
	nt->SetBranchAddress("PVzE",&PVzE);	
	nt->SetBranchAddress("Dsize",&Dsize);
	nt->SetBranchAddress("PVz",&PVz);
	nt->SetBranchAddress("Dmass",Dmass);
	nt->SetBranchAddress("Dpt",Dpt);
	nt->SetBranchAddress("Deta",Deta);
	nt->SetBranchAddress("Dphi",Dphi);
	nt->SetBranchAddress("Dy",Dy);
	nt->SetBranchAddress("Dchi2ndf",Dchi2ndf);
	nt->SetBranchAddress("Dchi2cl",Dchi2cl);
	nt->SetBranchAddress("Dalpha",Dalpha);
/*
	nt->SetBranchAddress("DsvpvDistance",DsvpvDistance);
	nt->SetBranchAddress("DsvpvDisErr",DsvpvDisErr);
	nt->SetBranchAddress("DsvpvDistance_2D",DsvpvDistance_2D);
	nt->SetBranchAddress("DsvpvDisErr_2D",DsvpvDisErr_2D);
*/

  nt->SetBranchAddress("DsvpvDistance", DsvpvDistance, &b_DsvpvDistance);
  nt->SetBranchAddress("DsvpvDisErr", DsvpvDisErr, &b_DsvpvDisErr);
  nt->SetBranchAddress("DsvpvDistance_2D", DsvpvDistance_2D, &b_DsvpvDistance_2D);
  nt->SetBranchAddress("DsvpvDisErr_2D", DsvpvDisErr_2D, &b_DsvpvDisErr_2D);

	nt->SetBranchAddress("DlxyBS",DlxyBS);
	nt->SetBranchAddress("DlxyBSErr",DlxyBSErr);
	nt->SetBranchAddress("DvtxX",DvtxX);
	nt->SetBranchAddress("DvtxY",DvtxY);
  nt->SetBranchAddress("DvtxZ",DvtxZ);
	nt->SetBranchAddress("Dtrk1Pt",Dtrk1Pt);
	nt->SetBranchAddress("Dtrk2Pt",Dtrk2Pt);
	nt->SetBranchAddress("Dtrk1Eta",Dtrk1Eta);
	nt->SetBranchAddress("Dtrk2Eta",Dtrk2Eta);
	nt->SetBranchAddress("Dtrk1Phi",Dtrk1Phi);
	nt->SetBranchAddress("Dtrk2Phi",Dtrk2Phi);
	nt->SetBranchAddress("Dtrk1PtErr",Dtrk1PtErr);
	nt->SetBranchAddress("Dtrk2PtErr",Dtrk2PtErr);
	nt->SetBranchAddress("Dtrk1PixelHit",Dtrk1PixelHit);
	nt->SetBranchAddress("Dtrk2PixelHit",Dtrk2PixelHit);
	nt->SetBranchAddress("Dtrk1StripHit",Dtrk1StripHit);
	nt->SetBranchAddress("Dtrk2StripHit",Dtrk2StripHit);
	nt->SetBranchAddress("Dtrk1nStripLayer",Dtrk1nStripLayer);
	nt->SetBranchAddress("Dtrk2nStripLayer",Dtrk2nStripLayer);
	nt->SetBranchAddress("Dtrk1nPixelLayer",Dtrk1nPixelLayer);
	nt->SetBranchAddress("Dtrk2nPixelLayer",Dtrk2nPixelLayer);
	nt->SetBranchAddress("Dtrk1Chi2ndf",Dtrk1Chi2ndf);
	nt->SetBranchAddress("Dtrk2Chi2ndf",Dtrk2Chi2ndf);
	nt->SetBranchAddress("Dtrk1highPurity",Dtrk1highPurity);
	nt->SetBranchAddress("Dtrk2highPurity",Dtrk2highPurity);
	nt->SetBranchAddress("Dgen",Dgen);
	nt->SetBranchAddress("DgennDa",DgennDa);
	nt->SetBranchAddress("DgenIndex",DgenIndex);
	nt->SetBranchAddress("Dgenpt",Dgenpt);
	nt->SetBranchAddress("Dgeneta",Dgeneta);
	nt->SetBranchAddress("Dgenphi",Dgenphi);
	nt->SetBranchAddress("Dgeny",Dgeny);
	nt->SetBranchAddress("DgenBAncestorpdgId",DgenBAncestorpdgId);
	nt->SetBranchAddress("DgenBAncestorpt",DgenBAncestorpt);
	nt->SetBranchAddress("DgencollisionId",DgencollisionId);

	if(hasDptFilteredWeight){
  nt->SetBranchAddress("DgenPtWeight",DgenPtWeight);
	}
	if(hasDgenVtx){
    nt->SetBranchAddress("DgenprodvtxX",DgenprodvtxX);
    nt->SetBranchAddress("DgenprodvtxY",DgenprodvtxY);
    nt->SetBranchAddress("DgenprodvtxZ",DgenprodvtxZ);
    nt->SetBranchAddress("DgendecayvtxX",DgendecayvtxX);
    nt->SetBranchAddress("DgendecayvtxY",DgendecayvtxY);
    nt->SetBranchAddress("DgendecayvtxZ",DgendecayvtxZ);
    nt->SetBranchAddress("DgenfromgenPV",DgenfromgenPV);	
	}
	if(hasHardQcdWeight){
  nt->SetBranchAddress("DgenPtPtHatWeight",DgenPtPtHatWeight);
	}

  if(hasFonllDptWeight){
  nt->SetBranchAddress("FonllDptWeight",FonllDptWeight);
  }


}



void RecoBranches(TTree* dnt,  Bool_t useDptFilteredWeight=false, Bool_t useDgenVtx=false, Bool_t useFonllDptWeight=false)
{

	dnt->Branch("RunNo",&RunNo);
	dnt->Branch("EvtNo",&EvtNo);
	dnt->Branch("LumiNo",&LumiNo);
	dnt->Branch("PVx",&PVx);
	dnt->Branch("PVy",&PVy);
	dnt->Branch("PVxE",&PVxE);
	dnt->Branch("PVyE",&PVyE);
	dnt->Branch("PVzE",&PVzE);	
	dnt->Branch("Dsize",&Dsize);
	dnt->Branch("PVz",&PVz);

	dnt->Branch("Dmass",Dmass,"Dmass[Dsize]/F");
	dnt->Branch("Dpt",Dpt,"Dpt[Dsize]/F");
	dnt->Branch("Deta",Deta,"Deta[Dsize]/F");
	dnt->Branch("Dphi",Dphi,"Dphi[Dsize]/F");
	dnt->Branch("Dy",Dy,"Dy[Dsize]/F");
	dnt->Branch("DvtxX",DvtxX,"DvtxX[Dsize]/F");
	dnt->Branch("DvtxY",DvtxY,"DvtxY[Dsize]/F");
  dnt->Branch("DvtxZ",DvtxZ,"DvtxZ[Dsize]/F");
	// dnt->Branch("Dd0",Dd0,"Dd0[Dsize]/F");
	// dnt->Branch("Dd0Err",Dd0Err,"Dd0Err[Dsize]/F");
	// dnt->Branch("Ddxyz",Ddxyz,"Ddxyz[Dsize]/F");
	// dnt->Branch("DdxyzErr",DdxyzErr,"DdxyzErr[Dsize]/F");
	dnt->Branch("Dchi2ndf",Dchi2ndf,"Dchi2ndf[Dsize]/F");
	dnt->Branch("Dchi2cl",Dchi2cl,"Dchi2cl[Dsize]/F");
	// dnt->Branch("Ddtheta",Ddtheta,"Ddtheta[Dsize]/F");
	// dnt->Branch("Dlxy",Dlxy,"Dlxy[Dsize]/F");
	dnt->Branch("Dalpha",Dalpha,"Dalpha[Dsize]/F");
	dnt->Branch("DsvpvDistance",DsvpvDistance,"DsvpvDistance[Dsize]/F");
	dnt->Branch("DsvpvDisErr",DsvpvDisErr,"DsvpvDisErr[Dsize]/F");
	dnt->Branch("DsvpvDistance_2D",DsvpvDistance_2D,"DsvpvDistance_2D[Dsize]/F");
	dnt->Branch("DsvpvDisErr_2D",DsvpvDisErr_2D,"DsvpvDisErr_2D[Dsize]/F");
	// dnt->Branch("DtktkRes_chi2ndf",DtktkRes_chi2ndf,"DtktkRes_chi2ndf[Dsize]/F");
	// dnt->Branch("DtktkRes_chi2cl",DtktkRes_chi2cl,"DtktkRes_chi2cl[Dsize]/F");
	// dnt->Branch("DtktkRes_alpha",DtktkRes_alpha,"DtktkRes_alpha[Dsize]/F");
	// dnt->Branch("DtktkRes_svpvDistance",DtktkRes_svpvDistance,"DtktkRes_svpvDistance[Dsize]/F");
	// dnt->Branch("DtktkRes_svpvDisErr",DtktkRes_svpvDisErr,"DtktkRes_svpvDisErr[Dsize]/F");
	dnt->Branch("DlxyBS",DlxyBS,"DlxyBS[Dsize]/F");
	dnt->Branch("DlxyBSErr",DlxyBSErr,"DlxyBSErr[Dsize]/F");
	// dnt->Branch("DMaxDoca",DMaxDoca,"DMaxDoca[Dsize]/F");
	dnt->Branch("Dtrk1Pt",Dtrk1Pt,"Dtrk1Pt[Dsize]/F");
	dnt->Branch("Dtrk2Pt",Dtrk2Pt,"Dtrk2Pt[Dsize]/F");
	dnt->Branch("Dtrk1Eta",Dtrk1Eta,"Dtrk1Eta[Dsize]/F");
	dnt->Branch("Dtrk2Eta",Dtrk2Eta,"Dtrk2Eta[Dsize]/F");
	dnt->Branch("Dtrk1Phi",Dtrk1Phi,"Dtrk1Phi[Dsize]/F");
	dnt->Branch("Dtrk2Phi",Dtrk2Phi,"Dtrk2Phi[Dsize]/F");
	dnt->Branch("Dtrk1PtErr",Dtrk1PtErr,"Dtrk1PtErr[Dsize]/F");
	dnt->Branch("Dtrk2PtErr",Dtrk2PtErr,"Dtrk2PtErr[Dsize]/F");
	// dnt->Branch("Dtrk1Dxy",Dtrk1Dxy,"Dtrk1Dxy[Dsize]/F");
	// dnt->Branch("Dtrk2Dxy",Dtrk2Dxy,"Dtrk2Dxy[Dsize]/F");
	dnt->Branch("Dtrk1PixelHit",Dtrk1PixelHit,"Dtrk1PixelHit[Dsize]/F");
	dnt->Branch("Dtrk2PixelHit",Dtrk2PixelHit,"Dtrk2PixelHit[Dsize]/F");
	dnt->Branch("Dtrk1StripHit",Dtrk1StripHit,"Dtrk1StripHit[Dsize]/F");
	dnt->Branch("Dtrk2StripHit",Dtrk2StripHit,"Dtrk2StripHit[Dsize]/F");
	dnt->Branch("Dtrk1nStripLayer",Dtrk1nStripLayer,"Dtrk1nStripLayer[Dsize]/F");
	dnt->Branch("Dtrk2nStripLayer",Dtrk2nStripLayer,"Dtrk2nStripLayer[Dsize]/F");
	dnt->Branch("Dtrk1nPixelLayer",Dtrk1nPixelLayer,"Dtrk1nPixelLayer[Dsize]/F");
	dnt->Branch("Dtrk2nPixelLayer",Dtrk2nPixelLayer,"Dtrk2nPixelLayer[Dsize]/F");
	dnt->Branch("Dtrk1Chi2ndf",Dtrk1Chi2ndf,"Dtrk1Chi2ndf[Dsize]/F");
	dnt->Branch("Dtrk2Chi2ndf",Dtrk2Chi2ndf,"Dtrk2Chi2ndf[Dsize]/F");
	// dnt->Branch("Dtrk1MassHypo",Dtrk1MassHypo,"Dtrk1MassHypo[Dsize]/F");
	// dnt->Branch("Dtrk2MassHypo",Dtrk2MassHypo,"Dtrk2MassHypo[Dsize]/F");
	// dnt->Branch("Dtrk1Algo",Dtrk1Algo,"Dtrk1Algo[Dsize]/I");
	// dnt->Branch("Dtrk2Algo",Dtrk2Algo,"Dtrk2Algo[Dsize]/I");
	// dnt->Branch("Dtrk1originalAlgo",Dtrk1originalAlgo,"Dtrk1originalAlgo[Dsize]/I");
	// dnt->Branch("Dtrk2originalAlgo",Dtrk2originalAlgo,"Dtrk2originalAlgo[Dsize]/I");
	dnt->Branch("Dtrk1highPurity",Dtrk1highPurity,"Dtrk1highPurity[Dsize]/O");
	dnt->Branch("Dtrk2highPurity",Dtrk2highPurity,"Dtrk2highPurity[Dsize]/O");

	dnt->Branch("Dgen",Dgen,"Dgen[Dsize]/I");
	dnt->Branch("DgenIndex",DgenIndex,"DgenIndex[Dsize]/I");
	dnt->Branch("DgennDa",DgennDa,"DgennDa[Dsize]/I");
	dnt->Branch("Dgenpt",Dgenpt,"Dgenpt[Dsize]/F");
	dnt->Branch("Dgeneta",Dgeneta,"Dgeneta[Dsize]/F");
	dnt->Branch("Dgenphi",Dgenphi,"Dgenphi[Dsize]/F");
	dnt->Branch("Dgeny",Dgeny,"Dgeny[Dsize]/F");
	dnt->Branch("DgencollisionId",DgencollisionId,"DgencollisionId[Dsize]/I");
	dnt->Branch("DgenBAncestorpt",DgenBAncestorpt,"DgenBAncestorpt[Dsize]/F");
	dnt->Branch("DgenBAncestorpdgId",DgenBAncestorpdgId,"DgenBAncestorpdgId[Dsize]/I");

	if(useDptFilteredWeight){
	dnt->Branch("DgenPtWeight",DgenPtWeight,"DgenPtWeight[Dsize]/F");
  dnt->Branch("DgenPtPtHatWeight",DgenPtPtHatWeight,"DgenPtPtHatWeight[Dsize]/F");
	}
	if(useDgenVtx){
    dnt->Branch("DgenprodvtxX",DgenprodvtxX,"DgenprodvtxX[Dsize]/F");
    dnt->Branch("DgenprodvtxY",DgenprodvtxY,"DgenprodvtxY[Dsize]/F");
    dnt->Branch("DgenprodvtxZ",DgenprodvtxZ,"DgenprodvtxZ[Dsize]/F");
    dnt->Branch("DgendecayvtxX",DgendecayvtxX,"DgendecayvtxX[Dsize]/F");
    dnt->Branch("DgendecayvtxY",DgendecayvtxY,"DgendecayvtxY[Dsize]/F");
    dnt->Branch("DgendecayvtxZ",DgendecayvtxZ,"DgendecayvtxZ[Dsize]/F");
    dnt->Branch("DgenfromgenPV",DgenfromgenPV,"DgenfromgenPV[Dsize]/I");
	}
	if(useFonllDptWeight){
		dnt->Branch("FonllDptWeight",FonllDptWeight,"FonllDptWeight[Dsize]/F");
	}


}

// Hlt Tree


Int_t     HLT_HIL1MinimumBiasHF1AND_v1;
Int_t     HLT_HIL1MinimumBiasHF2AND_part1_v1;
Int_t     HLT_HIL1MinimumBiasHF2AND_part2_v1;
Int_t     HLT_HIL1MinimumBiasHF2AND_part3_v1;
Int_t     HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part1_v1;
Int_t     HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part2_v1;
Int_t     HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part3_v1;
Int_t     HLT_HIL1Centralityext50100MinimumumBiasHF2AND_v1;
void SetHltBranches(TTree* nt, Bool_t isPbPb)
{
	if( isPbPb){
		nt->SetBranchAddress("HLT_HIL1MinimumBiasHF1AND_v1",&HLT_HIL1MinimumBiasHF1AND_v1);
		nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part1_v1",&HLT_HIL1MinimumBiasHF2AND_part1_v1);
		nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part2_v1",&HLT_HIL1MinimumBiasHF2AND_part2_v1);
		nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part3_v1",&HLT_HIL1MinimumBiasHF2AND_part3_v1);
		nt->SetBranchAddress("HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part1_v1",&HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part1_v1);
		nt->SetBranchAddress("HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part2_v1",&HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part2_v1);
		nt->SetBranchAddress("HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part3_v1",&HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part3_v1);
		nt->SetBranchAddress("HLT_HIL1Centralityext50100MinimumumBiasHF2AND_v1",&HLT_HIL1Centralityext50100MinimumumBiasHF2AND_v1);
	}

}

void HltBranches(TTree* nt, Bool_t isPbPb)
{
  if( isPbPb){
    nt->Branch("HLT_HIL1MinimumBiasHF1AND_v1",&HLT_HIL1MinimumBiasHF1AND_v1);
    nt->Branch("HLT_HIL1MinimumBiasHF2AND_part1_v1",&HLT_HIL1MinimumBiasHF2AND_part1_v1);
    nt->Branch("HLT_HIL1MinimumBiasHF2AND_part2_v1",&HLT_HIL1MinimumBiasHF2AND_part2_v1);
    nt->Branch("HLT_HIL1MinimumBiasHF2AND_part3_v1",&HLT_HIL1MinimumBiasHF2AND_part3_v1);
    nt->Branch("HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part1_v1",&HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part1_v1);
    nt->Branch("HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part2_v1",&HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part2_v1);
    nt->Branch("HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part3_v1",&HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part3_v1);
    nt->Branch("HLT_HIL1Centralityext50100MinimumumBiasHF2AND_v1",&HLT_HIL1Centralityext50100MinimumumBiasHF2AND_v1);
  }

}




Int_t     hiBin;
Int_t     hiNevtPlane;
Float_t   hiEvtPlanes[kMaxEvtPlanes];
Float_t   hiEvtPlanesqx[kMaxEvtPlanes];
Float_t   hiEvtPlanesqy[kMaxEvtPlanes];
Float_t   pthatweight;
Float_t   pthat;
Int_t     DPtSample;

void SetHiBranches(TTree* nt, Bool_t isData=false, Bool_t isPbPb=false, Bool_t isPthatWeight=false,  Bool_t hasDptFilteredWeight=false)
{
	nt->SetBranchAddress("pthat",&pthat);
  if(hasDptFilteredWeight){
	nt->SetBranchAddress("DPtSample",&DPtSample);
  }

	if(isPbPb){
		nt->SetBranchAddress("hiBin",&hiBin);
		nt->SetBranchAddress("hiNevtPlane",&hiNevtPlane);
		nt->SetBranchAddress("hiEvtPlanes",hiEvtPlanes);
		// nt->SetBranchAddress("hiEvtPlanesqx",hiEvtPlanesqx);
		// nt->SetBranchAddress("hiEvtPlanesqy",hiEvtPlanesqy);
	}
	if(isPthatWeight)
	{
		nt->SetBranchAddress("pthatweight",&pthatweight);
	}
}

void HiBranches(TTree* nt, Bool_t isData=false, Bool_t isPbPb=false, Bool_t isPthatWeight=false)
{
  nt->Branch("pthat",&pthat);
	nt->Branch("DPtSample",&DPtSample);
  if(isPbPb){
    nt->Branch("hiBin",&hiBin);
    nt->Branch("hiNevtPlane",&hiNevtPlane);
    nt->Branch("hiEvtPlanes",hiEvtPlanes,"hiEvtPlanes[hiNevtPlane]/F");
//    nt->Branch("hiEvtPlanesqx",hiEvtPlanesqx,"hiEvtPlanesqx[hiNevtPlane]/F");
//    nt->Branch("hiEvtPlanesqy",hiEvtPlanesqy,"hiEvtPlanesqy[hiNevtPlane]/F");
//    if(!isData) nt->SetBranchAddress("pthatweight",&pthatweight);
  }
  if(isPthatWeight)
  {
    nt->Branch("pthatweight",&pthatweight);
  }
}



	 Int_t     pcollisionEventSelection;
	 Int_t     pprimaryVertexFilter;
	 Int_t     phfCoincFilter3;
	 Int_t     pclusterCompatibilityFilter;

Int_t     pPAprimaryVertexFilter;
Int_t     pBeamScrapingFilter;


void SetSkimBranches(TTree* nt, Bool_t isPbPb=false)
{

	if(isPbPb)
	{
//		 nt->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
		 nt->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
		 nt->SetBranchAddress("phfCoincFilter3",&phfCoincFilter3);
		 nt->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);
	}
	if(!isPbPb)
	{
	nt->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter);
	nt->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
	}
}

void SkimBranches(TTree* nt, Bool_t isPbPb=false){

	if(isPbPb)
	{
//	nt->Branch("pcollisionEventSelection",&pcollisionEventSelection);
	nt->Branch("pprimaryVertexFilter",&pprimaryVertexFilter);
	nt->Branch("phfCoincFilter3",&phfCoincFilter3);
	nt->Branch("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);	
	}

	
	if(!isPbPb)
	{
	nt->Branch("pPAprimaryVertexFilter",&pPAprimaryVertexFilter);
  nt->Branch("pBeamScrapingFilter",&pBeamScrapingFilter);
	}
}



// ntGenClass

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
Float_t         GPtWeight[MAX_XB];

Float_t 				GPtPtHatWeight[MAX_XB];
Float_t 			FonllGptWeight[MAX_XB];


// Gen production vertex and decay vertex
Int_t     GfromgenPV[MAX_XB];
Float_t   GprodvtxX[MAX_XB];//gen production vertex
Float_t   GprodvtxY[MAX_XB];
Float_t   GprodvtxZ[MAX_XB];
Float_t   GdecayvtxX[MAX_XB];//gen decay vertex
Float_t   GdecayvtxY[MAX_XB];
Float_t   GdecayvtxZ[MAX_XB];




void SetGenBranches(TTree* nt, Bool_t hasDptFilteredWeight=false, Bool_t hasDgenVtx=false, Bool_t hasHardQcdWeight=false, Bool_t hasFonllDptWeight=false)
{

  nt->SetBranchAddress("GPVx", &GPVx);
  nt->SetBranchAddress("GPVy", &GPVy);
  nt->SetBranchAddress("GPVz", &GPVz);
	nt->SetBranchAddress("Gsize", &Gsize);
	nt->SetBranchAddress("Gy", Gy);
	nt->SetBranchAddress("Geta", Geta);
	nt->SetBranchAddress("Gphi", Gphi);
	nt->SetBranchAddress("Gpt", Gpt);
	nt->SetBranchAddress("GpdgId", GpdgId);
	nt->SetBranchAddress("GcollisionId", GcollisionId);
	nt->SetBranchAddress("GisSignal", GisSignal);
	nt->SetBranchAddress("GBAncestorpt", GBAncestorpt);
	nt->SetBranchAddress("GBAncestorpdgId", GBAncestorpdgId);
	nt->SetBranchAddress("Gtk1pt", Gtk1pt);
	nt->SetBranchAddress("Gtk1eta", Gtk1eta);
	nt->SetBranchAddress("Gtk1y", Gtk1y);
	nt->SetBranchAddress("Gtk1phi", Gtk1phi);
	nt->SetBranchAddress("Gtk2pt", Gtk2pt);
	nt->SetBranchAddress("Gtk2eta", Gtk2eta);
	nt->SetBranchAddress("Gtk2y", Gtk2y);
	nt->SetBranchAddress("Gtk2phi", Gtk2phi);

	if(hasDptFilteredWeight){
    nt->SetBranchAddress("GPtWeight",GPtWeight);
	}
	if(hasDgenVtx){
    nt->SetBranchAddress("GfromgenPV",GfromgenPV);
    nt->SetBranchAddress("GprodvtxX",GprodvtxX);
    nt->SetBranchAddress("GprodvtxY",GprodvtxY);
    nt->SetBranchAddress("GprodvtxZ",GprodvtxZ);
    nt->SetBranchAddress("GdecayvtxX",GdecayvtxX);
    nt->SetBranchAddress("GdecayvtxY",GdecayvtxY);
    nt->SetBranchAddress("GdecayvtxZ",GdecayvtxZ);
	
	}
  if(hasHardQcdWeight){
  nt->SetBranchAddress("GPtPtHatWeight",GPtPtHatWeight);
  }

	if(hasFonllDptWeight){
	nt->SetBranchAddress("FonllGptWeight",FonllGptWeight);
	}


}

void GenBranches(TTree* nt,  bool useDptFilteredWeight=false, bool useDgenVtx=false, Bool_t useFonllDptWeight=false){

    nt->Branch("GPVx",&GPVx);
    nt->Branch("GPVy",&GPVy);
    nt->Branch("GPVz",&GPVz);
    nt->Branch("Gsize",&Gsize);
    nt->Branch("Gy",Gy,"Gy[Gsize]/F");
    nt->Branch("Geta",Geta,"Geta[Gsize]/F");
    nt->Branch("Gphi",Gphi,"Gphi[Gsize]/F");
    nt->Branch("Gpt",Gpt,"Gpt[Gsize]/F");
    nt->Branch("GpdgId",GpdgId,"GpdgId[Gsize]/I");
    nt->Branch("GcollisionId",GcollisionId,"GcollisionId[Gsize]/I");
    nt->Branch("GisSignal",GisSignal,"GisSignal[Gsize]/I");
    nt->Branch("GBAncestorpt",GBAncestorpt,"GBAncestorpt[Gsize]/F");
    nt->Branch("GBAncestorpdgId",GBAncestorpdgId,"GBAncestorpdgId[Gsize]/I");
    nt->Branch("Gtk1pt",Gtk1pt,"Gtk1pt[Gsize]/F");
    nt->Branch("Gtk1eta",Gtk1eta,"Gtk1eta[Gsize]/F");
    nt->Branch("Gtk1y",Gtk1y,"Gtk1y[Gsize]/F");
    nt->Branch("Gtk1phi",Gtk1phi,"Gtk1phi[Gsize]/F");
    nt->Branch("Gtk2pt",Gtk2pt,"Gtk2pt[Gsize]/F");
    nt->Branch("Gtk2eta",Gtk2eta,"Gtk2eta[Gsize]/F");
    nt->Branch("Gtk2y",Gtk2y,"Gtk2y[Gsize]/F");
    nt->Branch("Gtk2phi",Gtk2phi,"Gtk2phi[Gsize]/F");
		if(useDptFilteredWeight){
		nt->Branch("GPtWeight",GPtWeight,"GPtWeight[Gsize]/F");
		nt->Branch("GPtPtHatWeight",GPtPtHatWeight,"GPtPtHatWeight[Gsize]/F");
		}
		if(useDgenVtx){
    nt->Branch("GfromgenPV",GfromgenPV,"GfromgenPV[Gsize]/I");
    nt->Branch("GprodvtxX",GprodvtxX,"GprodvtxX[Gsize]/F");
    nt->Branch("GprodvtxY",GprodvtxY,"GprodvtxY[Gsize]/F");
    nt->Branch("GprodvtxZ",GprodvtxZ,"GprodvtxZ[Gsize]/F");
    nt->Branch("GdecayvtxX",GdecayvtxX,"GdecayvtxX[Gsize]/F");
    nt->Branch("GdecayvtxY",GdecayvtxY,"GdecayvtxY[Gsize]/F");
    nt->Branch("GdecayvtxZ",GdecayvtxZ,"GdecayvtxZ[Gsize]/F");

		}

  if(useFonllDptWeight){
    nt->Branch("FonllGptWeight",FonllGptWeight,"FonllGptWeight[Gsize]/F");
  }


}




#endif
