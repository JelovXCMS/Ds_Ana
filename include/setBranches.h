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


// renew 20180301 for Ds Ntuple , hlt need to be updated
#include <TBranch.h>
#include <TTree.h>

#ifndef _ANA_BFEEDDOWN_SETBRANCHES_H_
#define _ANA_BFEEDDOWN_SETBRANCHES_H_

//#include "uti.h"
#define MAX_XB        20000
#define kMaxEvtPlanes 200

// DntupleBranch

   Int_t           RunNo;
   Int_t           EvtNo;
   Int_t           LumiNo;
   Int_t           Dsize;
   Float_t         PVx;
   Float_t         PVy;
   Float_t         PVz;
   Float_t         PVnchi2;
   Float_t         BSx;
   Float_t         BSy;
   Float_t         BSz;
   Float_t         PVxE;
   Float_t         PVyE;
   Float_t         PVzE;
   Float_t         BSxErr;
   Float_t         BSyErr;
   Float_t         BSzErr;
   Float_t         BSdxdz;
   Float_t         BSdydz;
   Float_t         BSdxdzErr;
   Float_t         BSdydzErr;
   Float_t         BSWidthX;
   Float_t         BSWidthXErr;
   Float_t         BSWidthY;
   Float_t         BSWidthYErr;
   Int_t           Dindex[MAX_XB];   //[Dsize]
   Int_t           Dtype[MAX_XB];   //[Dsize]
   Float_t         Dmass[MAX_XB];   //[Dsize]
   Float_t         D_unfitted_mass[MAX_XB];   //[Dsize]
   Float_t         Dpt[MAX_XB];   //[Dsize]
   Float_t         D_unfitted_pt[MAX_XB];   //[Dsize]
   Float_t         Deta[MAX_XB];   //[Dsize]
   Float_t         Dphi[MAX_XB];   //[Dsize]
   Float_t         Dy[MAX_XB];   //[Dsize]
   Float_t         DvtxX[MAX_XB];   //[Dsize]
   Float_t         DvtxY[MAX_XB];   //[Dsize]
   Float_t         DvtxZ[MAX_XB];   //[Dsize]
   Float_t         DvtxXErr[MAX_XB];   //[Dsize]
   Float_t         DvtxYErr[MAX_XB];   //[Dsize]
   Float_t         DvtxZErr[MAX_XB];   //[Dsize]
   Float_t         Dd0[MAX_XB];   //[Dsize]
   Float_t         Dd0Err[MAX_XB];   //[Dsize]
   Float_t         Ddxyz[MAX_XB];   //[Dsize]
   Float_t         DdxyzErr[MAX_XB];   //[Dsize]
   Float_t         Dchi2ndf[MAX_XB];   //[Dsize]
   Float_t         Dchi2cl[MAX_XB];   //[Dsize]
   Float_t         Ddtheta[MAX_XB];   //[Dsize]
   Float_t         Dlxy[MAX_XB];   //[Dsize]
   Float_t         Dalpha[MAX_XB];   //[Dsize]
   Float_t         Dalpha_BS_2D[MAX_XB];   //[Dsize]
   Float_t         DsvpvDistance[MAX_XB];   //[Dsize]
   Float_t         DsvpvDisErr[MAX_XB];   //[Dsize]
   Float_t         DsvpvDistance_2D[MAX_XB];   //[Dsize]
   Float_t         DsvpvDisErr_2D[MAX_XB];   //[Dsize]
   Float_t         Ddca[MAX_XB];   //[Dsize]
   Float_t         DlxyBS[MAX_XB];   //[Dsize]
   Float_t         DlxyBSErr[MAX_XB];   //[Dsize]
   Float_t         DMaxDoca[MAX_XB];   //[Dsize]
   Float_t         DMaxTkPt[MAX_XB];   //[Dsize]
   Float_t         DMinTkPt[MAX_XB];   //[Dsize]
   Float_t         Dtrk1Pt[MAX_XB];   //[Dsize]
   Float_t         Dtrk2Pt[MAX_XB];   //[Dsize]
   Float_t         Dtrk1PtErr[MAX_XB];   //[Dsize]
   Float_t         Dtrk2PtErr[MAX_XB];   //[Dsize]
   Float_t         Dtrk1Eta[MAX_XB];   //[Dsize]
   Float_t         Dtrk2Eta[MAX_XB];   //[Dsize]
   Float_t         Dtrk1Phi[MAX_XB];   //[Dsize]
   Float_t         Dtrk2Phi[MAX_XB];   //[Dsize]
   Float_t         Dtrk1P[MAX_XB];   //[Dsize]
   Float_t         Dtrk2P[MAX_XB];   //[Dsize]
   Float_t         Dtrk1Dz[MAX_XB];   //[Dsize]
   Float_t         Dtrk2Dz[MAX_XB];   //[Dsize]
   Float_t         Dtrk1Dxy[MAX_XB];   //[Dsize]
   Float_t         Dtrk2Dxy[MAX_XB];   //[Dsize]
   Float_t         Dtrk1MassHypo[MAX_XB];   //[Dsize]
   Float_t         Dtrk2MassHypo[MAX_XB];   //[Dsize]
   Int_t           Dtrk1originalAlgo[MAX_XB];   //[Dsize]
   Int_t           Dtrk2originalAlgo[MAX_XB];   //[Dsize]
   Bool_t          Dtrk1highPurity[MAX_XB];   //[Dsize]
   Bool_t          Dtrk2highPurity[MAX_XB];   //[Dsize]
   Float_t         Dtrk1dedx[MAX_XB];   //[Dsize]
   Float_t         Dtrk2dedx[MAX_XB];   //[Dsize]
   Float_t         Dtrk1thetastar[MAX_XB];   //[Dsize]
   Float_t         Dtrk2thetastar[MAX_XB];   //[Dsize]
   Float_t         Dtrk1thetastar_uf[MAX_XB];   //[Dsize]
   Float_t         Dtrk2thetastar_uf[MAX_XB];   //[Dsize]
   Float_t         Dtrk1PixelHit[MAX_XB];   //[Dsize]
   Float_t         Dtrk2PixelHit[MAX_XB];   //[Dsize]
   Float_t         Dtrk1StripHit[MAX_XB];   //[Dsize]
   Float_t         Dtrk2StripHit[MAX_XB];   //[Dsize]
   Float_t         Dtrk1nStripLayer[MAX_XB];   //[Dsize]
   Float_t         Dtrk2nStripLayer[MAX_XB];   //[Dsize]
   Float_t         Dtrk1nPixelLayer[MAX_XB];   //[Dsize]
   Float_t         Dtrk2nPixelLayer[MAX_XB];   //[Dsize]
   Float_t         Dtrk1Chi2ndf[MAX_XB];   //[Dsize]
   Float_t         Dtrk2Chi2ndf[MAX_XB];   //[Dsize]
   Int_t           Dtrk1Algo[MAX_XB];   //[Dsize]
   Int_t           Dtrk2Algo[MAX_XB];   //[Dsize]
   Float_t         Dtrk3Pt[MAX_XB];   //[Dsize]
   Float_t         Dtrk4Pt[MAX_XB];   //[Dsize]
   Float_t         Dtrk3PtErr[MAX_XB];   //[Dsize]
   Float_t         Dtrk4PtErr[MAX_XB];   //[Dsize]
   Float_t         Dtrk3Eta[MAX_XB];   //[Dsize]
   Float_t         Dtrk4Eta[MAX_XB];   //[Dsize]
   Float_t         Dtrk3Phi[MAX_XB];   //[Dsize]
   Float_t         Dtrk4Phi[MAX_XB];   //[Dsize]
   Float_t         Dtrk3P[MAX_XB];   //[Dsize]
   Float_t         Dtrk4P[MAX_XB];   //[Dsize]
   Float_t         Dtrk3MassHypo[MAX_XB];   //[Dsize]
   Float_t         Dtrk4MassHypo[MAX_XB];   //[Dsize]
   Float_t         Dtrk3Dz[MAX_XB];   //[Dsize]
   Float_t         Dtrk4Dz[MAX_XB];   //[Dsize]
   Float_t         Dtrk3Dxy[MAX_XB];   //[Dsize]
   Float_t         Dtrk4Dxy[MAX_XB];   //[Dsize]
   Int_t           Dtrk3originalAlgo[MAX_XB];   //[Dsize]
   Int_t           Dtrk4originalAlgo[MAX_XB];   //[Dsize]
   Bool_t          Dtrk3highPurity[MAX_XB];   //[Dsize]
   Bool_t          Dtrk4highPurity[MAX_XB];   //[Dsize]
   Float_t         Dtrk3dedx[MAX_XB];   //[Dsize]
   Float_t         Dtrk4dedx[MAX_XB];   //[Dsize]
   Float_t         Dtrk3thetastar[MAX_XB];   //[Dsize]
   Float_t         Dtrk4thetastar[MAX_XB];   //[Dsize]
   Float_t         Dtrk3thetastar_uf[MAX_XB];   //[Dsize]
   Float_t         Dtrk4thetastar_uf[MAX_XB];   //[Dsize]
   Float_t         Dtrk3PixelHit[MAX_XB];   //[Dsize]
   Float_t         Dtrk4PixelHit[MAX_XB];   //[Dsize]
   Float_t         Dtrk3StripHit[MAX_XB];   //[Dsize]
   Float_t         Dtrk4StripHit[MAX_XB];   //[Dsize]
   Float_t         Dtrk3nStripLayer[MAX_XB];   //[Dsize]
   Float_t         Dtrk4nStripLayer[MAX_XB];   //[Dsize]
   Float_t         Dtrk3nPixelLayer[MAX_XB];   //[Dsize]
   Float_t         Dtrk4nPixelLayer[MAX_XB];   //[Dsize]
   Float_t         Dtrk3Chi2ndf[MAX_XB];   //[Dsize]
   Float_t         Dtrk4Chi2ndf[MAX_XB];   //[Dsize]
   Int_t           Dtrk3Algo[MAX_XB];   //[Dsize]
   Int_t           Dtrk4Algo[MAX_XB];   //[Dsize]
   Float_t         Dtrk1Y[MAX_XB];   //[Dsize]
   Float_t         Dtrk2Y[MAX_XB];   //[Dsize]
   Float_t         Dtrk1D0[MAX_XB];   //[Dsize]
   Float_t         Dtrk2D0[MAX_XB];   //[Dsize]
   Float_t         Dtrk1D0Err[MAX_XB];   //[Dsize]
   Float_t         Dtrk2D0Err[MAX_XB];   //[Dsize]
   Int_t           Dtrk1Quality[MAX_XB];   //[Dsize]
   Int_t           Dtrk2Quality[MAX_XB];   //[Dsize]
   Int_t           Dtrk1Idx[MAX_XB];   //[Dsize]
   Int_t           Dtrk2Idx[MAX_XB];   //[Dsize]
   Float_t         Dtrk1EtaErr[MAX_XB];   //[Dsize]
   Float_t         Dtrk2EtaErr[MAX_XB];   //[Dsize]
   Float_t         Dtrk1PhiErr[MAX_XB];   //[Dsize]
   Float_t         Dtrk2PhiErr[MAX_XB];   //[Dsize]
   Float_t         Dtrk1MVAVal[MAX_XB];   //[Dsize]
   Float_t         Dtrk2MVAVal[MAX_XB];   //[Dsize]
   Float_t         Dtrk3Y[MAX_XB];   //[Dsize]
   Float_t         Dtrk4Y[MAX_XB];   //[Dsize]
   Float_t         Dtrk3D0[MAX_XB];   //[Dsize]
   Float_t         Dtrk4D0[MAX_XB];   //[Dsize]
   Float_t         Dtrk3D0Err[MAX_XB];   //[Dsize]
   Float_t         Dtrk4D0Err[MAX_XB];   //[Dsize]
   Int_t           Dtrk3Quality[MAX_XB];   //[Dsize]
   Int_t           Dtrk4Quality[MAX_XB];   //[Dsize]
   Int_t           Dtrk3Idx[MAX_XB];   //[Dsize]
   Int_t           Dtrk4Idx[MAX_XB];   //[Dsize]
   Float_t         Dtrk3EtaErr[MAX_XB];   //[Dsize]
   Float_t         Dtrk4EtaErr[MAX_XB];   //[Dsize]
   Float_t         Dtrk3PhiErr[MAX_XB];   //[Dsize]
   Float_t         Dtrk4PhiErr[MAX_XB];   //[Dsize]
   Float_t         Dtrk3MVAVal[MAX_XB];   //[Dsize]
   Float_t         Dtrk4MVAVal[MAX_XB];   //[Dsize]
   Float_t         DtktkResmass[MAX_XB];   //[Dsize]
   Float_t         DtktkRes_unfitted_mass[MAX_XB];   //[Dsize]
   Float_t         DtktkRespt[MAX_XB];   //[Dsize]
   Float_t         DtktkRes_unfitted_pt[MAX_XB];   //[Dsize]
   Float_t         DtktkReseta[MAX_XB];   //[Dsize]
   Float_t         DtktkResphi[MAX_XB];   //[Dsize]
   Float_t         DtktkRes_chi2ndf[MAX_XB];   //[Dsize]
   Float_t         DtktkRes_chi2cl[MAX_XB];   //[Dsize]
   Float_t         DtktkRes_alpha[MAX_XB];   //[Dsize]
   Float_t         DtktkRes_alphaToSV[MAX_XB];   //[Dsize]
   Float_t         DtktkRes_svpvDistance[MAX_XB];   //[Dsize]
   Float_t         DtktkRes_svpvDisErr[MAX_XB];   //[Dsize]
   Float_t         DtktkRes_svpvDistanceToSV[MAX_XB];   //[Dsize]
   Float_t         DtktkRes_svpvDisErrToSV[MAX_XB];   //[Dsize]
   Float_t         DtktkRes_dca[MAX_XB];   //[Dsize]
   Float_t         DtktkRes_dcaToSV[MAX_XB];   //[Dsize]
   Float_t         DtktkRes_lxyBS[MAX_XB];   //[Dsize]
   Float_t         DtktkRes_lxyBSErr[MAX_XB];   //[Dsize]
   Float_t         DtktkRes_angleToTrk1[MAX_XB];   //[Dsize]
   Float_t         DtktkRes_unfitted_angleToTrk1[MAX_XB];   //[Dsize]
   Float_t         DtktkRes_ptAsymToTrk1[MAX_XB];   //[Dsize]
   Float_t         DtktkRes_unfitter_ptAsymToTrk1[MAX_XB];   //[Dsize]
   Float_t         DRestrk1Pt[MAX_XB];   //[Dsize]
   Float_t         DRestrk2Pt[MAX_XB];   //[Dsize]
   Float_t         DRestrk3Pt[MAX_XB];   //[Dsize]
   Float_t         DRestrk4Pt[MAX_XB];   //[Dsize]
   Float_t         DRestrk1PtErr[MAX_XB];   //[Dsize]
   Float_t         DRestrk2PtErr[MAX_XB];   //[Dsize]
   Float_t         DRestrk3PtErr[MAX_XB];   //[Dsize]
   Float_t         DRestrk4PtErr[MAX_XB];   //[Dsize]
   Float_t         DRestrk1Eta[MAX_XB];   //[Dsize]
   Float_t         DRestrk2Eta[MAX_XB];   //[Dsize]
   Float_t         DRestrk3Eta[MAX_XB];   //[Dsize]
   Float_t         DRestrk4Eta[MAX_XB];   //[Dsize]
   Float_t         DRestrk1Phi[MAX_XB];   //[Dsize]
   Float_t         DRestrk2Phi[MAX_XB];   //[Dsize]
   Float_t         DRestrk3Phi[MAX_XB];   //[Dsize]
   Float_t         DRestrk4Phi[MAX_XB];   //[Dsize]
   Float_t         DRestrk1P[MAX_XB];   //[Dsize]
   Float_t         DRestrk2P[MAX_XB];   //[Dsize]
   Float_t         DRestrk3P[MAX_XB];   //[Dsize]
   Float_t         DRestrk4P[MAX_XB];   //[Dsize]
   Float_t         DRestrk1MassHypo[MAX_XB];   //[Dsize]
   Float_t         DRestrk2MassHypo[MAX_XB];   //[Dsize]
   Float_t         DRestrk3MassHypo[MAX_XB];   //[Dsize]
   Float_t         DRestrk4MassHypo[MAX_XB];   //[Dsize]
   Float_t         DRestrk1Dz[MAX_XB];   //[Dsize]
   Float_t         DRestrk2Dz[MAX_XB];   //[Dsize]
   Float_t         DRestrk3Dz[MAX_XB];   //[Dsize]
   Float_t         DRestrk4Dz[MAX_XB];   //[Dsize]
   Float_t         DRestrk1Dxy[MAX_XB];   //[Dsize]
   Float_t         DRestrk2Dxy[MAX_XB];   //[Dsize]
   Float_t         DRestrk3Dxy[MAX_XB];   //[Dsize]
   Float_t         DRestrk4Dxy[MAX_XB];   //[Dsize]
   Int_t           DRestrk1originalAlgo[MAX_XB];   //[Dsize]
   Int_t           DRestrk2originalAlgo[MAX_XB];   //[Dsize]
   Int_t           DRestrk3originalAlgo[MAX_XB];   //[Dsize]
   Int_t           DRestrk4originalAlgo[MAX_XB];   //[Dsize]
   Bool_t          DRestrk1highPurity[MAX_XB];   //[Dsize]
   Bool_t          DRestrk2highPurity[MAX_XB];   //[Dsize]
   Bool_t          DRestrk3highPurity[MAX_XB];   //[Dsize]
   Bool_t          DRestrk4highPurity[MAX_XB];   //[Dsize]
   Float_t         DRestrk1dedx[MAX_XB];   //[Dsize]
   Float_t         DRestrk2dedx[MAX_XB];   //[Dsize]
   Float_t         DRestrk3dedx[MAX_XB];   //[Dsize]
   Float_t         DRestrk4dedx[MAX_XB];   //[Dsize]
   Float_t         DRestrk1thetastar[MAX_XB];   //[Dsize]
   Float_t         DRestrk2thetastar[MAX_XB];   //[Dsize]
   Float_t         DRestrk3thetastar[MAX_XB];   //[Dsize]
   Float_t         DRestrk4thetastar[MAX_XB];   //[Dsize]
   Float_t         DRestrk1thetastar_uf[MAX_XB];   //[Dsize]
   Float_t         DRestrk2thetastar_uf[MAX_XB];   //[Dsize]
   Float_t         DRestrk3thetastar_uf[MAX_XB];   //[Dsize]
   Float_t         DRestrk4thetastar_uf[MAX_XB];   //[Dsize]
   Float_t         DRestrk1Y[MAX_XB];   //[Dsize]
   Float_t         DRestrk2Y[MAX_XB];   //[Dsize]
   Float_t         DRestrk3Y[MAX_XB];   //[Dsize]
   Float_t         DRestrk4Y[MAX_XB];   //[Dsize]
   Float_t         DRestrk1D0[MAX_XB];   //[Dsize]
   Float_t         DRestrk2D0[MAX_XB];   //[Dsize]
   Float_t         DRestrk3D0[MAX_XB];   //[Dsize]
   Float_t         DRestrk4D0[MAX_XB];   //[Dsize]
   Float_t         DRestrk1D0Err[MAX_XB];   //[Dsize]
   Float_t         DRestrk2D0Err[MAX_XB];   //[Dsize]
   Float_t         DRestrk3D0Err[MAX_XB];   //[Dsize]
   Float_t         DRestrk4D0Err[MAX_XB];   //[Dsize]
   Int_t           DRestrk1Quality[MAX_XB];   //[Dsize]
   Int_t           DRestrk2Quality[MAX_XB];   //[Dsize]
   Int_t           DRestrk3Quality[MAX_XB];   //[Dsize]
   Int_t           DRestrk4Quality[MAX_XB];   //[Dsize]
   Int_t           Dgen[MAX_XB];   //[Dsize]
   Int_t           DsGen[MAX_XB];   //[Dsize]
   Int_t           DgenIndex[MAX_XB];   //[Dsize]
   Int_t           DgennDa[MAX_XB];   //[Dsize]
   Float_t         Dgenpt[MAX_XB];   //[Dsize]
   Float_t         Dgeneta[MAX_XB];   //[Dsize]
   Float_t         Dgenphi[MAX_XB];   //[Dsize]
   Float_t         Dgeny[MAX_XB];   //[Dsize]
   Int_t           DgencollisionId[MAX_XB];   //[Dsize]
   Float_t         DgenBAncestorpt[MAX_XB];   //[Dsize]
   Int_t           DgenBAncestorpdgId[MAX_XB];   //[Dsize]
   Float_t         DgenprodvtxX[MAX_XB];   //[Dsize]
   Float_t         DgenprodvtxY[MAX_XB];   //[Dsize]
   Float_t         DgenprodvtxZ[MAX_XB];   //[Dsize]
   Float_t         DgendecayvtxX[MAX_XB];   //[Dsize]
   Float_t         DgendecayvtxY[MAX_XB];   //[Dsize]
   Float_t         DgendecayvtxZ[MAX_XB];   //[Dsize]
   Float_t         DgendecayvtxXtoPv[MAX_XB];   //[Dsize]
   Float_t         DgendecayvtxYtoPv[MAX_XB];   //[Dsize]
   Float_t         DgendecayvtxZtoPv[MAX_XB];   //[Dsize]
   Int_t           DgenfromgenPV[MAX_XB];   //[Dsize]

   // List of branches
   TBranch        *b_RunNo;   //!
   TBranch        *b_EvtNo;   //!
   TBranch        *b_LumiNo;   //!
   TBranch        *b_Dsize;   //!
   TBranch        *b_PVx;   //!
   TBranch        *b_PVy;   //!
   TBranch        *b_PVz;   //!
   TBranch        *b_PVnchi2;   //!
   TBranch        *b_BSx;   //!
   TBranch        *b_BSy;   //!
   TBranch        *b_BSz;   //!
   TBranch        *b_PVxE;   //!
   TBranch        *b_PVyE;   //!
   TBranch        *b_PVzE;   //!
   TBranch        *b_BSxErr;   //!
   TBranch        *b_BSyErr;   //!
   TBranch        *b_BSzErr;   //!
   TBranch        *b_BSdxdz;   //!
   TBranch        *b_BSdydz;   //!
   TBranch        *b_BSdxdzErr;   //!
   TBranch        *b_BSdydzErr;   //!
   TBranch        *b_BSWidthX;   //!
   TBranch        *b_BSWidthXErr;   //!
   TBranch        *b_BSWidthY;   //!
   TBranch        *b_BSWidthYErr;   //!
   TBranch        *b_Dindex;   //!
   TBranch        *b_Dtype;   //!
   TBranch        *b_Dmass;   //!
   TBranch        *b_D_unfitted_mass;   //!
   TBranch        *b_Dpt;   //!
   TBranch        *b_D_unfitted_pt;   //!
   TBranch        *b_Deta;   //!
   TBranch        *b_Dphi;   //!
   TBranch        *b_Dy;   //!
   TBranch        *b_DvtxX;   //!
   TBranch        *b_DvtxY;   //!
   TBranch        *b_DvtxZ;   //!
   TBranch        *b_DvtxXErr;   //!
   TBranch        *b_DvtxYErr;   //!
   TBranch        *b_DvtxZErr;   //!
   TBranch        *b_Dd0;   //!
   TBranch        *b_Dd0Err;   //!
   TBranch        *b_Ddxyz;   //!
   TBranch        *b_DdxyzErr;   //!
   TBranch        *b_Dchi2ndf;   //!
   TBranch        *b_Dchi2cl;   //!
   TBranch        *b_Ddtheta;   //!
   TBranch        *b_Dlxy;   //!
   TBranch        *b_Dalpha;   //!
   TBranch        *b_Dalpha_BS_2D;   //!
   TBranch        *b_DsvpvDistance;   //!
   TBranch        *b_DsvpvDisErr;   //!
   TBranch        *b_DsvpvDistance_2D;   //!
   TBranch        *b_DsvpvDisErr_2D;   //!
   TBranch        *b_Ddca;   //!
   TBranch        *b_DlxyBS;   //!
   TBranch        *b_DlxyBSErr;   //!
   TBranch        *b_DMaxDoca;   //!
   TBranch        *b_DMaxTkPt;   //!
   TBranch        *b_DMinTkPt;   //!
   TBranch        *b_Dtrk1Pt;   //!
   TBranch        *b_Dtrk2Pt;   //!
   TBranch        *b_Dtrk1PtErr;   //!
   TBranch        *b_Dtrk2PtErr;   //!
   TBranch        *b_Dtrk1Eta;   //!
   TBranch        *b_Dtrk2Eta;   //!
   TBranch        *b_Dtrk1Phi;   //!
   TBranch        *b_Dtrk2Phi;   //!
   TBranch        *b_Dtrk1P;   //!
   TBranch        *b_Dtrk2P;   //!
   TBranch        *b_Dtrk1Dz;   //!
   TBranch        *b_Dtrk2Dz;   //!
   TBranch        *b_Dtrk1Dxy;   //!
   TBranch        *b_Dtrk2Dxy;   //!
   TBranch        *b_Dtrk1MassHypo;   //!
   TBranch        *b_Dtrk2MassHypo;   //!
   TBranch        *b_Dtrk1originalAlgo;   //!
   TBranch        *b_Dtrk2originalAlgo;   //!
   TBranch        *b_Dtrk1highPurity;   //!
   TBranch        *b_Dtrk2highPurity;   //!
   TBranch        *b_Dtrk1dedx;   //!
   TBranch        *b_Dtrk2dedx;   //!
   TBranch        *b_Dtrk1thetastar;   //!
   TBranch        *b_Dtrk2thetastar;   //!
   TBranch        *b_Dtrk1thetastar_uf;   //!
   TBranch        *b_Dtrk2thetastar_uf;   //!
   TBranch        *b_Dtrk1PixelHit;   //!
   TBranch        *b_Dtrk2PixelHit;   //!
   TBranch        *b_Dtrk1StripHit;   //!
   TBranch        *b_Dtrk2StripHit;   //!
   TBranch        *b_Dtrk1nStripLayer;   //!
   TBranch        *b_Dtrk2nStripLayer;   //!
   TBranch        *b_Dtrk1nPixelLayer;   //!
   TBranch        *b_Dtrk2nPixelLayer;   //!
   TBranch        *b_Dtrk1Chi2ndf;   //!
   TBranch        *b_Dtrk2Chi2ndf;   //!
   TBranch        *b_Dtrk1Algo;   //!
   TBranch        *b_Dtrk2Algo;   //!
   TBranch        *b_Dtrk3Pt;   //!
   TBranch        *b_Dtrk4Pt;   //!
   TBranch        *b_Dtrk3PtErr;   //!
   TBranch        *b_Dtrk4PtErr;   //!
   TBranch        *b_Dtrk3Eta;   //!
   TBranch        *b_Dtrk4Eta;   //!
   TBranch        *b_Dtrk3Phi;   //!
   TBranch        *b_Dtrk4Phi;   //!
   TBranch        *b_Dtrk3P;   //!
   TBranch        *b_Dtrk4P;   //!
   TBranch        *b_Dtrk3MassHypo;   //!
   TBranch        *b_Dtrk4MassHypo;   //!
   TBranch        *b_Dtrk3Dz;   //!
   TBranch        *b_Dtrk4Dz;   //!
   TBranch        *b_Dtrk3Dxy;   //!
   TBranch        *b_Dtrk4Dxy;   //!
   TBranch        *b_Dtrk3originalAlgo;   //!
   TBranch        *b_Dtrk4originalAlgo;   //!
   TBranch        *b_Dtrk3highPurity;   //!
   TBranch        *b_Dtrk4highPurity;   //!
   TBranch        *b_Dtrk3dedx;   //!
   TBranch        *b_Dtrk4dedx;   //!
   TBranch        *b_Dtrk3thetastar;   //!
   TBranch        *b_Dtrk4thetastar;   //!
   TBranch        *b_Dtrk3thetastar_uf;   //!
   TBranch        *b_Dtrk4thetastar_uf;   //!
   TBranch        *b_Dtrk3PixelHit;   //!
   TBranch        *b_Dtrk4PixelHit;   //!
   TBranch        *b_Dtrk3StripHit;   //!
   TBranch        *b_Dtrk4StripHit;   //!
   TBranch        *b_Dtrk3nStripLayer;   //!
   TBranch        *b_Dtrk4nStripLayer;   //!
   TBranch        *b_Dtrk3nPixelLayer;   //!
   TBranch        *b_Dtrk4nPixelLayer;   //!
   TBranch        *b_Dtrk3Chi2ndf;   //!
   TBranch        *b_Dtrk4Chi2ndf;   //!
   TBranch        *b_Dtrk3Algo;   //!
   TBranch        *b_Dtrk4Algo;   //!
   TBranch        *b_Dtrk1Y;   //!
   TBranch        *b_Dtrk2Y;   //!
   TBranch        *b_Dtrk1D0;   //!
   TBranch        *b_Dtrk2D0;   //!
   TBranch        *b_Dtrk1D0Err;   //!
   TBranch        *b_Dtrk2D0Err;   //!
   TBranch        *b_Dtrk1Quality;   //!
   TBranch        *b_Dtrk2Quality;   //!
   TBranch        *b_Dtrk1Idx;   //!
   TBranch        *b_Dtrk2Idx;   //!
   TBranch        *b_Dtrk1EtaErr;   //!
   TBranch        *b_Dtrk2EtaErr;   //!
   TBranch        *b_Dtrk1PhiErr;   //!
   TBranch        *b_Dtrk2PhiErr;   //!
   TBranch        *b_Dtrk1MVAVal;   //!
   TBranch        *b_Dtrk2MVAVal;   //!
   TBranch        *b_Dtrk3Y;   //!
   TBranch        *b_Dtrk4Y;   //!
   TBranch        *b_Dtrk3D0;   //!
   TBranch        *b_Dtrk4D0;   //!
   TBranch        *b_Dtrk3D0Err;   //!
   TBranch        *b_Dtrk4D0Err;   //!
   TBranch        *b_Dtrk3Quality;   //!
   TBranch        *b_Dtrk4Quality;   //!
   TBranch        *b_Dtrk3Idx;   //!
   TBranch        *b_Dtrk4Idx;   //!
   TBranch        *b_Dtrk3EtaErr;   //!
   TBranch        *b_Dtrk4EtaErr;   //!
   TBranch        *b_Dtrk3PhiErr;   //!
   TBranch        *b_Dtrk4PhiErr;   //!
   TBranch        *b_Dtrk3MVAVal;   //!
   TBranch        *b_Dtrk4MVAVal;   //!
   TBranch        *b_DtktkResmass;   //!
   TBranch        *b_DtktkRes_unfitted_mass;   //!
   TBranch        *b_DtktkRespt;   //!
   TBranch        *b_DtktkRes_unfitted_pt;   //!
   TBranch        *b_DtktkReseta;   //!
   TBranch        *b_DtktkResphi;   //!
   TBranch        *b_DtktkRes_chi2ndf;   //!
   TBranch        *b_DtktkRes_chi2cl;   //!
   TBranch        *b_DtktkRes_alpha;   //!
   TBranch        *b_DtktkRes_alphaToSV;   //!
   TBranch        *b_DtktkRes_svpvDistance;   //!
   TBranch        *b_DtktkRes_svpvDisErr;   //!
   TBranch        *b_DtktkRes_svpvDistanceToSV;   //!
   TBranch        *b_DtktkRes_svpvDisErrToSV;   //!
   TBranch        *b_DtktkRes_dca;   //!
   TBranch        *b_DtktkRes_dcaToSV;   //!
   TBranch        *b_DtktkRes_lxyBS;   //!
   TBranch        *b_DtktkRes_lxyBSErr;   //!
   TBranch        *b_DtktkRes_angleToTrk1;   //!
   TBranch        *b_DtktkRes_unfitted_angleToTrk1;   //!
   TBranch        *b_DtktkRes_ptAsymToTrk1;   //!
   TBranch        *b_DtktkRes_unfitter_ptAsymToTrk1;   //!
   TBranch        *b_DRestrk1Pt;   //!
   TBranch        *b_DRestrk2Pt;   //!
   TBranch        *b_DRestrk3Pt;   //!
   TBranch        *b_DRestrk4Pt;   //!
   TBranch        *b_DRestrk1PtErr;   //!
   TBranch        *b_DRestrk2PtErr;   //!
   TBranch        *b_DRestrk3PtErr;   //!
   TBranch        *b_DRestrk4PtErr;   //!
   TBranch        *b_DRestrk1Eta;   //!
   TBranch        *b_DRestrk2Eta;   //!
   TBranch        *b_DRestrk3Eta;   //!
   TBranch        *b_DRestrk4Eta;   //!
   TBranch        *b_DRestrk1Phi;   //!
   TBranch        *b_DRestrk2Phi;   //!
   TBranch        *b_DRestrk3Phi;   //!
   TBranch        *b_DRestrk4Phi;   //!
   TBranch        *b_DRestrk1P;   //!
   TBranch        *b_DRestrk2P;   //!
   TBranch        *b_DRestrk3P;   //!
   TBranch        *b_DRestrk4P;   //!
   TBranch        *b_DRestrk1MassHypo;   //!
   TBranch        *b_DRestrk2MassHypo;   //!
   TBranch        *b_DRestrk3MassHypo;   //!
   TBranch        *b_DRestrk4MassHypo;   //!
   TBranch        *b_DRestrk1Dz;   //!
   TBranch        *b_DRestrk2Dz;   //!
   TBranch        *b_DRestrk3Dz;   //!
   TBranch        *b_DRestrk4Dz;   //!
   TBranch        *b_DRestrk1Dxy;   //!
   TBranch        *b_DRestrk2Dxy;   //!
   TBranch        *b_DRestrk3Dxy;   //!
   TBranch        *b_DRestrk4Dxy;   //!
   TBranch        *b_DRestrk1originalAlgo;   //!
   TBranch        *b_DRestrk2originalAlgo;   //!
   TBranch        *b_DRestrk3originalAlgo;   //!
   TBranch        *b_DRestrk4originalAlgo;   //!
   TBranch        *b_DRestrk1highPurity;   //!
   TBranch        *b_DRestrk2highPurity;   //!
   TBranch        *b_DRestrk3highPurity;   //!
   TBranch        *b_DRestrk4highPurity;   //!
   TBranch        *b_DRestrk1dedx;   //!
   TBranch        *b_DRestrk2dedx;   //!
   TBranch        *b_DRestrk3dedx;   //!
   TBranch        *b_DRestrk4dedx;   //!
   TBranch        *b_DRestrk1thetastar;   //!
   TBranch        *b_DRestrk2thetastar;   //!
   TBranch        *b_DRestrk3thetastar;   //!
   TBranch        *b_DRestrk4thetastar;   //!
   TBranch        *b_DRestrk1thetastar_uf;   //!
   TBranch        *b_DRestrk2thetastar_uf;   //!
   TBranch        *b_DRestrk3thetastar_uf;   //!
   TBranch        *b_DRestrk4thetastar_uf;   //!
   TBranch        *b_DRestrk1Y;   //!
   TBranch        *b_DRestrk2Y;   //!
   TBranch        *b_DRestrk3Y;   //!
   TBranch        *b_DRestrk4Y;   //!
   TBranch        *b_DRestrk1D0;   //!
   TBranch        *b_DRestrk2D0;   //!
   TBranch        *b_DRestrk3D0;   //!
   TBranch        *b_DRestrk4D0;   //!
   TBranch        *b_DRestrk1D0Err;   //!
   TBranch        *b_DRestrk2D0Err;   //!
   TBranch        *b_DRestrk3D0Err;   //!
   TBranch        *b_DRestrk4D0Err;   //!
   TBranch        *b_DRestrk1Quality;   //!
   TBranch        *b_DRestrk2Quality;   //!
   TBranch        *b_DRestrk3Quality;   //!
   TBranch        *b_DRestrk4Quality;   //!
   TBranch        *b_Dgen;   //!
   TBranch        *b_DsGen;   //!
   TBranch        *b_DgenIndex;   //!
   TBranch        *b_DgennDa;   //!
   TBranch        *b_Dgenpt;   //!
   TBranch        *b_Dgeneta;   //!
   TBranch        *b_Dgenphi;   //!
   TBranch        *b_Dgeny;   //!
   TBranch        *b_DgencollisionId;   //!
   TBranch        *b_DgenBAncestorpt;   //!
   TBranch        *b_DgenBAncestorpdgId;   //!
   TBranch        *b_DgenprodvtxX;   //!
   TBranch        *b_DgenprodvtxY;   //!
   TBranch        *b_DgenprodvtxZ;   //!
   TBranch        *b_DgendecayvtxX;   //!
   TBranch        *b_DgendecayvtxY;   //!
   TBranch        *b_DgendecayvtxZ;   //!
   TBranch        *b_DgendecayvtxXtoPv;   //!
   TBranch        *b_DgendecayvtxYtoPv;   //!
   TBranch        *b_DgendecayvtxZtoPv;   //!
   TBranch        *b_DgenfromgenPV;   //!


void SetRecoBranches(TTree* fChain, Bool_t isData=true, Bool_t detailMode=false)
{

   fChain->SetBranchAddress("RunNo", &RunNo, &b_RunNo);
   fChain->SetBranchAddress("EvtNo", &EvtNo, &b_EvtNo);
   fChain->SetBranchAddress("LumiNo", &LumiNo, &b_LumiNo);
   fChain->SetBranchAddress("Dsize", &Dsize, &b_Dsize);
   fChain->SetBranchAddress("PVx", &PVx, &b_PVx);
   fChain->SetBranchAddress("PVy", &PVy, &b_PVy);
   fChain->SetBranchAddress("PVz", &PVz, &b_PVz);
   fChain->SetBranchAddress("PVnchi2", &PVnchi2, &b_PVnchi2);
   fChain->SetBranchAddress("BSx", &BSx, &b_BSx);
   fChain->SetBranchAddress("BSy", &BSy, &b_BSy);
   fChain->SetBranchAddress("BSz", &BSz, &b_BSz);
   fChain->SetBranchAddress("PVxE", &PVxE, &b_PVxE);
   fChain->SetBranchAddress("PVyE", &PVyE, &b_PVyE);
   fChain->SetBranchAddress("PVzE", &PVzE, &b_PVzE);
   if(detailMode)
   {
   // fChain->SetBranchAddress("PVxE", &PVxE, &b_PVxE);
   // fChain->SetBranchAddress("PVyE", &PVyE, &b_PVyE);
   // fChain->SetBranchAddress("PVzE", &PVzE, &b_PVzE);
   fChain->SetBranchAddress("BSxErr", &BSxErr, &b_BSxErr);
   fChain->SetBranchAddress("BSyErr", &BSyErr, &b_BSyErr);
   fChain->SetBranchAddress("BSzErr", &BSzErr, &b_BSzErr);
   fChain->SetBranchAddress("BSdxdz", &BSdxdz, &b_BSdxdz);
   fChain->SetBranchAddress("BSdydz", &BSdydz, &b_BSdydz);
   fChain->SetBranchAddress("BSdxdzErr", &BSdxdzErr, &b_BSdxdzErr);
   fChain->SetBranchAddress("BSdydzErr", &BSdydzErr, &b_BSdydzErr);
   fChain->SetBranchAddress("BSWidthX", &BSWidthX, &b_BSWidthX);
   fChain->SetBranchAddress("BSWidthXErr", &BSWidthXErr, &b_BSWidthXErr);
   fChain->SetBranchAddress("BSWidthY", &BSWidthY, &b_BSWidthY);
   fChain->SetBranchAddress("BSWidthYErr", &BSWidthYErr, &b_BSWidthYErr);
	 }
   fChain->SetBranchAddress("Dindex", Dindex, &b_Dindex);
   fChain->SetBranchAddress("Dtype", Dtype, &b_Dtype);
   fChain->SetBranchAddress("Dmass", Dmass, &b_Dmass);
   fChain->SetBranchAddress("D_unfitted_mass", D_unfitted_mass, &b_D_unfitted_mass);
   fChain->SetBranchAddress("Dpt", Dpt, &b_Dpt);
   fChain->SetBranchAddress("D_unfitted_pt", D_unfitted_pt, &b_D_unfitted_pt);
   fChain->SetBranchAddress("Deta", Deta, &b_Deta);
   fChain->SetBranchAddress("Dphi", Dphi, &b_Dphi);
   fChain->SetBranchAddress("Dy", Dy, &b_Dy);
   fChain->SetBranchAddress("DvtxX", DvtxX, &b_DvtxX);
   fChain->SetBranchAddress("DvtxY", DvtxY, &b_DvtxY);
   fChain->SetBranchAddress("DvtxZ", DvtxZ, &b_DvtxZ);
   fChain->SetBranchAddress("DvtxXErr", DvtxXErr, &b_DvtxXErr);
   fChain->SetBranchAddress("DvtxYErr", DvtxYErr, &b_DvtxYErr);
   fChain->SetBranchAddress("DvtxZErr", DvtxZErr, &b_DvtxZErr);
   fChain->SetBranchAddress("Dd0", Dd0, &b_Dd0);
   fChain->SetBranchAddress("Dd0Err", Dd0Err, &b_Dd0Err);
   fChain->SetBranchAddress("Ddxyz", Ddxyz, &b_Ddxyz);
   fChain->SetBranchAddress("DdxyzErr", DdxyzErr, &b_DdxyzErr);
   fChain->SetBranchAddress("Dchi2ndf", Dchi2ndf, &b_Dchi2ndf);
   fChain->SetBranchAddress("Dchi2cl", Dchi2cl, &b_Dchi2cl);
   fChain->SetBranchAddress("Ddtheta", Ddtheta, &b_Ddtheta);
   fChain->SetBranchAddress("Dlxy", Dlxy, &b_Dlxy);
   fChain->SetBranchAddress("Dalpha", Dalpha, &b_Dalpha);
   fChain->SetBranchAddress("Dalpha_BS_2D", Dalpha_BS_2D, &b_Dalpha_BS_2D);
   fChain->SetBranchAddress("DsvpvDistance", DsvpvDistance, &b_DsvpvDistance);
   fChain->SetBranchAddress("DsvpvDisErr", DsvpvDisErr, &b_DsvpvDisErr);
   fChain->SetBranchAddress("DsvpvDistance_2D", DsvpvDistance_2D, &b_DsvpvDistance_2D);
   fChain->SetBranchAddress("DsvpvDisErr_2D", DsvpvDisErr_2D, &b_DsvpvDisErr_2D);
   fChain->SetBranchAddress("Ddca", Ddca, &b_Ddca);
   fChain->SetBranchAddress("DlxyBS", DlxyBS, &b_DlxyBS);
   fChain->SetBranchAddress("DlxyBSErr", DlxyBSErr, &b_DlxyBSErr);
   fChain->SetBranchAddress("DMaxDoca", DMaxDoca, &b_DMaxDoca);
   fChain->SetBranchAddress("DMaxTkPt", DMaxTkPt, &b_DMaxTkPt);
   fChain->SetBranchAddress("DMinTkPt", DMinTkPt, &b_DMinTkPt);
   fChain->SetBranchAddress("Dtrk1Pt", Dtrk1Pt, &b_Dtrk1Pt);
   fChain->SetBranchAddress("Dtrk2Pt", Dtrk2Pt, &b_Dtrk2Pt);
   fChain->SetBranchAddress("Dtrk1PtErr", Dtrk1PtErr, &b_Dtrk1PtErr);
   fChain->SetBranchAddress("Dtrk2PtErr", Dtrk2PtErr, &b_Dtrk2PtErr);
   fChain->SetBranchAddress("Dtrk1Eta", Dtrk1Eta, &b_Dtrk1Eta);
   fChain->SetBranchAddress("Dtrk2Eta", Dtrk2Eta, &b_Dtrk2Eta);
   fChain->SetBranchAddress("Dtrk1Phi", Dtrk1Phi, &b_Dtrk1Phi);
   fChain->SetBranchAddress("Dtrk2Phi", Dtrk2Phi, &b_Dtrk2Phi);
   fChain->SetBranchAddress("Dtrk1P", Dtrk1P, &b_Dtrk1P);
   fChain->SetBranchAddress("Dtrk2P", Dtrk2P, &b_Dtrk2P);
   fChain->SetBranchAddress("Dtrk1Dz", Dtrk1Dz, &b_Dtrk1Dz);
   fChain->SetBranchAddress("Dtrk2Dz", Dtrk2Dz, &b_Dtrk2Dz);
   fChain->SetBranchAddress("Dtrk1Dxy", Dtrk1Dxy, &b_Dtrk1Dxy);
   fChain->SetBranchAddress("Dtrk2Dxy", Dtrk2Dxy, &b_Dtrk2Dxy);
   fChain->SetBranchAddress("Dtrk1MassHypo", Dtrk1MassHypo, &b_Dtrk1MassHypo);
   fChain->SetBranchAddress("Dtrk2MassHypo", Dtrk2MassHypo, &b_Dtrk2MassHypo);
   fChain->SetBranchAddress("Dtrk1originalAlgo", Dtrk1originalAlgo, &b_Dtrk1originalAlgo);
   fChain->SetBranchAddress("Dtrk2originalAlgo", Dtrk2originalAlgo, &b_Dtrk2originalAlgo);
   fChain->SetBranchAddress("Dtrk1highPurity", Dtrk1highPurity, &b_Dtrk1highPurity);
   fChain->SetBranchAddress("Dtrk2highPurity", Dtrk2highPurity, &b_Dtrk2highPurity);
   fChain->SetBranchAddress("Dtrk1dedx", Dtrk1dedx, &b_Dtrk1dedx);
   fChain->SetBranchAddress("Dtrk2dedx", Dtrk2dedx, &b_Dtrk2dedx);
   fChain->SetBranchAddress("Dtrk1thetastar", Dtrk1thetastar, &b_Dtrk1thetastar);
   fChain->SetBranchAddress("Dtrk2thetastar", Dtrk2thetastar, &b_Dtrk2thetastar);
   fChain->SetBranchAddress("Dtrk1thetastar_uf", Dtrk1thetastar_uf, &b_Dtrk1thetastar_uf);
   fChain->SetBranchAddress("Dtrk2thetastar_uf", Dtrk2thetastar_uf, &b_Dtrk2thetastar_uf);
   fChain->SetBranchAddress("Dtrk1PixelHit", Dtrk1PixelHit, &b_Dtrk1PixelHit);
   fChain->SetBranchAddress("Dtrk2PixelHit", Dtrk2PixelHit, &b_Dtrk2PixelHit);
   fChain->SetBranchAddress("Dtrk1StripHit", Dtrk1StripHit, &b_Dtrk1StripHit);
   fChain->SetBranchAddress("Dtrk2StripHit", Dtrk2StripHit, &b_Dtrk2StripHit);
   fChain->SetBranchAddress("Dtrk1nStripLayer", Dtrk1nStripLayer, &b_Dtrk1nStripLayer);
   fChain->SetBranchAddress("Dtrk2nStripLayer", Dtrk2nStripLayer, &b_Dtrk2nStripLayer);
   fChain->SetBranchAddress("Dtrk1nPixelLayer", Dtrk1nPixelLayer, &b_Dtrk1nPixelLayer);
   fChain->SetBranchAddress("Dtrk2nPixelLayer", Dtrk2nPixelLayer, &b_Dtrk2nPixelLayer);
   fChain->SetBranchAddress("Dtrk1Chi2ndf", Dtrk1Chi2ndf, &b_Dtrk1Chi2ndf);
   fChain->SetBranchAddress("Dtrk2Chi2ndf", Dtrk2Chi2ndf, &b_Dtrk2Chi2ndf);
   fChain->SetBranchAddress("Dtrk1Algo", Dtrk1Algo, &b_Dtrk1Algo);
   fChain->SetBranchAddress("Dtrk2Algo", Dtrk2Algo, &b_Dtrk2Algo);
   fChain->SetBranchAddress("Dtrk3Pt", Dtrk3Pt, &b_Dtrk3Pt);
   fChain->SetBranchAddress("Dtrk4Pt", Dtrk4Pt, &b_Dtrk4Pt);
   fChain->SetBranchAddress("Dtrk3PtErr", Dtrk3PtErr, &b_Dtrk3PtErr);
   fChain->SetBranchAddress("Dtrk4PtErr", Dtrk4PtErr, &b_Dtrk4PtErr);
   fChain->SetBranchAddress("Dtrk3Eta", Dtrk3Eta, &b_Dtrk3Eta);
   fChain->SetBranchAddress("Dtrk4Eta", Dtrk4Eta, &b_Dtrk4Eta);
   fChain->SetBranchAddress("Dtrk3Phi", Dtrk3Phi, &b_Dtrk3Phi);
   fChain->SetBranchAddress("Dtrk4Phi", Dtrk4Phi, &b_Dtrk4Phi);
   fChain->SetBranchAddress("Dtrk3P", Dtrk3P, &b_Dtrk3P);
   fChain->SetBranchAddress("Dtrk4P", Dtrk4P, &b_Dtrk4P);
   fChain->SetBranchAddress("Dtrk3MassHypo", Dtrk3MassHypo, &b_Dtrk3MassHypo);
   fChain->SetBranchAddress("Dtrk4MassHypo", Dtrk4MassHypo, &b_Dtrk4MassHypo);
   fChain->SetBranchAddress("Dtrk3Dz", Dtrk3Dz, &b_Dtrk3Dz);
   fChain->SetBranchAddress("Dtrk4Dz", Dtrk4Dz, &b_Dtrk4Dz);
   fChain->SetBranchAddress("Dtrk3Dxy", Dtrk3Dxy, &b_Dtrk3Dxy);
   fChain->SetBranchAddress("Dtrk4Dxy", Dtrk4Dxy, &b_Dtrk4Dxy);
   fChain->SetBranchAddress("Dtrk3originalAlgo", Dtrk3originalAlgo, &b_Dtrk3originalAlgo);
   fChain->SetBranchAddress("Dtrk4originalAlgo", Dtrk4originalAlgo, &b_Dtrk4originalAlgo);
   fChain->SetBranchAddress("Dtrk3highPurity", Dtrk3highPurity, &b_Dtrk3highPurity);
   fChain->SetBranchAddress("Dtrk4highPurity", Dtrk4highPurity, &b_Dtrk4highPurity);
   fChain->SetBranchAddress("Dtrk3dedx", Dtrk3dedx, &b_Dtrk3dedx);
   fChain->SetBranchAddress("Dtrk4dedx", Dtrk4dedx, &b_Dtrk4dedx);
   fChain->SetBranchAddress("Dtrk3thetastar", Dtrk3thetastar, &b_Dtrk3thetastar);
   fChain->SetBranchAddress("Dtrk4thetastar", Dtrk4thetastar, &b_Dtrk4thetastar);
   fChain->SetBranchAddress("Dtrk3thetastar_uf", Dtrk3thetastar_uf, &b_Dtrk3thetastar_uf);
   fChain->SetBranchAddress("Dtrk4thetastar_uf", Dtrk4thetastar_uf, &b_Dtrk4thetastar_uf);
   fChain->SetBranchAddress("Dtrk3PixelHit", Dtrk3PixelHit, &b_Dtrk3PixelHit);
   fChain->SetBranchAddress("Dtrk4PixelHit", Dtrk4PixelHit, &b_Dtrk4PixelHit);
   fChain->SetBranchAddress("Dtrk3StripHit", Dtrk3StripHit, &b_Dtrk3StripHit);
   fChain->SetBranchAddress("Dtrk4StripHit", Dtrk4StripHit, &b_Dtrk4StripHit);
   fChain->SetBranchAddress("Dtrk3nStripLayer", Dtrk3nStripLayer, &b_Dtrk3nStripLayer);
   fChain->SetBranchAddress("Dtrk4nStripLayer", Dtrk4nStripLayer, &b_Dtrk4nStripLayer);
   fChain->SetBranchAddress("Dtrk3nPixelLayer", Dtrk3nPixelLayer, &b_Dtrk3nPixelLayer);
   fChain->SetBranchAddress("Dtrk4nPixelLayer", Dtrk4nPixelLayer, &b_Dtrk4nPixelLayer);
   fChain->SetBranchAddress("Dtrk3Chi2ndf", Dtrk3Chi2ndf, &b_Dtrk3Chi2ndf);
   fChain->SetBranchAddress("Dtrk4Chi2ndf", Dtrk4Chi2ndf, &b_Dtrk4Chi2ndf);
   fChain->SetBranchAddress("Dtrk3Algo", Dtrk3Algo, &b_Dtrk3Algo);
   fChain->SetBranchAddress("Dtrk4Algo", Dtrk4Algo, &b_Dtrk4Algo);
    if(detailMode)
	 {
   fChain->SetBranchAddress("Dtrk1Y", Dtrk1Y, &b_Dtrk1Y);
   fChain->SetBranchAddress("Dtrk2Y", Dtrk2Y, &b_Dtrk2Y);
   fChain->SetBranchAddress("Dtrk1D0", Dtrk1D0, &b_Dtrk1D0);
   fChain->SetBranchAddress("Dtrk2D0", Dtrk2D0, &b_Dtrk2D0);
   fChain->SetBranchAddress("Dtrk1D0Err", Dtrk1D0Err, &b_Dtrk1D0Err);
   fChain->SetBranchAddress("Dtrk2D0Err", Dtrk2D0Err, &b_Dtrk2D0Err);
   fChain->SetBranchAddress("Dtrk1Quality", Dtrk1Quality, &b_Dtrk1Quality);
   fChain->SetBranchAddress("Dtrk2Quality", Dtrk2Quality, &b_Dtrk2Quality);
   fChain->SetBranchAddress("Dtrk1Idx", Dtrk1Idx, &b_Dtrk1Idx);
   fChain->SetBranchAddress("Dtrk2Idx", Dtrk2Idx, &b_Dtrk2Idx);
   fChain->SetBranchAddress("Dtrk1EtaErr", Dtrk1EtaErr, &b_Dtrk1EtaErr);
   fChain->SetBranchAddress("Dtrk2EtaErr", Dtrk2EtaErr, &b_Dtrk2EtaErr);
   fChain->SetBranchAddress("Dtrk1PhiErr", Dtrk1PhiErr, &b_Dtrk1PhiErr);
   fChain->SetBranchAddress("Dtrk2PhiErr", Dtrk2PhiErr, &b_Dtrk2PhiErr);
   fChain->SetBranchAddress("Dtrk1MVAVal", Dtrk1MVAVal, &b_Dtrk1MVAVal);
   fChain->SetBranchAddress("Dtrk2MVAVal", Dtrk2MVAVal, &b_Dtrk2MVAVal);
   fChain->SetBranchAddress("Dtrk3Y", Dtrk3Y, &b_Dtrk3Y);
   fChain->SetBranchAddress("Dtrk4Y", Dtrk4Y, &b_Dtrk4Y);
   fChain->SetBranchAddress("Dtrk3D0", Dtrk3D0, &b_Dtrk3D0);
   fChain->SetBranchAddress("Dtrk4D0", Dtrk4D0, &b_Dtrk4D0);
   fChain->SetBranchAddress("Dtrk3D0Err", Dtrk3D0Err, &b_Dtrk3D0Err);
   fChain->SetBranchAddress("Dtrk4D0Err", Dtrk4D0Err, &b_Dtrk4D0Err);
   fChain->SetBranchAddress("Dtrk3Quality", Dtrk3Quality, &b_Dtrk3Quality);
   fChain->SetBranchAddress("Dtrk4Quality", Dtrk4Quality, &b_Dtrk4Quality);
   fChain->SetBranchAddress("Dtrk3Idx", Dtrk3Idx, &b_Dtrk3Idx);
   fChain->SetBranchAddress("Dtrk4Idx", Dtrk4Idx, &b_Dtrk4Idx);
   fChain->SetBranchAddress("Dtrk3EtaErr", Dtrk3EtaErr, &b_Dtrk3EtaErr);
   fChain->SetBranchAddress("Dtrk4EtaErr", Dtrk4EtaErr, &b_Dtrk4EtaErr);
   fChain->SetBranchAddress("Dtrk3PhiErr", Dtrk3PhiErr, &b_Dtrk3PhiErr);
   fChain->SetBranchAddress("Dtrk4PhiErr", Dtrk4PhiErr, &b_Dtrk4PhiErr);
   fChain->SetBranchAddress("Dtrk3MVAVal", Dtrk3MVAVal, &b_Dtrk3MVAVal);
   fChain->SetBranchAddress("Dtrk4MVAVal", Dtrk4MVAVal, &b_Dtrk4MVAVal);
	 }
   fChain->SetBranchAddress("DtktkResmass", DtktkResmass, &b_DtktkResmass);
   fChain->SetBranchAddress("DtktkRes_unfitted_mass", DtktkRes_unfitted_mass, &b_DtktkRes_unfitted_mass);
   fChain->SetBranchAddress("DtktkRespt", DtktkRespt, &b_DtktkRespt);
   fChain->SetBranchAddress("DtktkRes_unfitted_pt", DtktkRes_unfitted_pt, &b_DtktkRes_unfitted_pt);
   fChain->SetBranchAddress("DtktkReseta", DtktkReseta, &b_DtktkReseta);
   fChain->SetBranchAddress("DtktkResphi", DtktkResphi, &b_DtktkResphi);
   fChain->SetBranchAddress("DtktkRes_chi2ndf", DtktkRes_chi2ndf, &b_DtktkRes_chi2ndf);
   fChain->SetBranchAddress("DtktkRes_chi2cl", DtktkRes_chi2cl, &b_DtktkRes_chi2cl);
   fChain->SetBranchAddress("DtktkRes_alpha", DtktkRes_alpha, &b_DtktkRes_alpha);
   fChain->SetBranchAddress("DtktkRes_alphaToSV", DtktkRes_alphaToSV, &b_DtktkRes_alphaToSV);
   fChain->SetBranchAddress("DtktkRes_svpvDistance", DtktkRes_svpvDistance, &b_DtktkRes_svpvDistance);
   fChain->SetBranchAddress("DtktkRes_svpvDisErr", DtktkRes_svpvDisErr, &b_DtktkRes_svpvDisErr);
   fChain->SetBranchAddress("DtktkRes_svpvDistanceToSV", DtktkRes_svpvDistanceToSV, &b_DtktkRes_svpvDistanceToSV);
   fChain->SetBranchAddress("DtktkRes_svpvDisErrToSV", DtktkRes_svpvDisErrToSV, &b_DtktkRes_svpvDisErrToSV);
   fChain->SetBranchAddress("DtktkRes_dca", DtktkRes_dca, &b_DtktkRes_dca);
   fChain->SetBranchAddress("DtktkRes_dcaToSV", DtktkRes_dcaToSV, &b_DtktkRes_dcaToSV);
   fChain->SetBranchAddress("DtktkRes_lxyBS", DtktkRes_lxyBS, &b_DtktkRes_lxyBS);
   fChain->SetBranchAddress("DtktkRes_lxyBSErr", DtktkRes_lxyBSErr, &b_DtktkRes_lxyBSErr);
   fChain->SetBranchAddress("DtktkRes_angleToTrk1", DtktkRes_angleToTrk1, &b_DtktkRes_angleToTrk1);
   fChain->SetBranchAddress("DtktkRes_unfitted_angleToTrk1", DtktkRes_unfitted_angleToTrk1, &b_DtktkRes_unfitted_angleToTrk1);
   fChain->SetBranchAddress("DtktkRes_ptAsymToTrk1", DtktkRes_ptAsymToTrk1, &b_DtktkRes_ptAsymToTrk1);
   fChain->SetBranchAddress("DtktkRes_unfitter_ptAsymToTrk1", DtktkRes_unfitter_ptAsymToTrk1, &b_DtktkRes_unfitter_ptAsymToTrk1);
   fChain->SetBranchAddress("DRestrk1Pt", DRestrk1Pt, &b_DRestrk1Pt);
   fChain->SetBranchAddress("DRestrk2Pt", DRestrk2Pt, &b_DRestrk2Pt);
   fChain->SetBranchAddress("DRestrk3Pt", DRestrk3Pt, &b_DRestrk3Pt);
   fChain->SetBranchAddress("DRestrk4Pt", DRestrk4Pt, &b_DRestrk4Pt);
   fChain->SetBranchAddress("DRestrk1PtErr", DRestrk1PtErr, &b_DRestrk1PtErr);
   fChain->SetBranchAddress("DRestrk2PtErr", DRestrk2PtErr, &b_DRestrk2PtErr);
   fChain->SetBranchAddress("DRestrk3PtErr", DRestrk3PtErr, &b_DRestrk3PtErr);
   fChain->SetBranchAddress("DRestrk4PtErr", DRestrk4PtErr, &b_DRestrk4PtErr);
   fChain->SetBranchAddress("DRestrk1Eta", DRestrk1Eta, &b_DRestrk1Eta);
   fChain->SetBranchAddress("DRestrk2Eta", DRestrk2Eta, &b_DRestrk2Eta);
   fChain->SetBranchAddress("DRestrk3Eta", DRestrk3Eta, &b_DRestrk3Eta);
   fChain->SetBranchAddress("DRestrk4Eta", DRestrk4Eta, &b_DRestrk4Eta);
   fChain->SetBranchAddress("DRestrk1Phi", DRestrk1Phi, &b_DRestrk1Phi);
   fChain->SetBranchAddress("DRestrk2Phi", DRestrk2Phi, &b_DRestrk2Phi);
   fChain->SetBranchAddress("DRestrk3Phi", DRestrk3Phi, &b_DRestrk3Phi);
   fChain->SetBranchAddress("DRestrk4Phi", DRestrk4Phi, &b_DRestrk4Phi);
   fChain->SetBranchAddress("DRestrk1P", DRestrk1P, &b_DRestrk1P);
   fChain->SetBranchAddress("DRestrk2P", DRestrk2P, &b_DRestrk2P);
   fChain->SetBranchAddress("DRestrk3P", DRestrk3P, &b_DRestrk3P);
   fChain->SetBranchAddress("DRestrk4P", DRestrk4P, &b_DRestrk4P);
   fChain->SetBranchAddress("DRestrk1MassHypo", DRestrk1MassHypo, &b_DRestrk1MassHypo);
   fChain->SetBranchAddress("DRestrk2MassHypo", DRestrk2MassHypo, &b_DRestrk2MassHypo);
   fChain->SetBranchAddress("DRestrk3MassHypo", DRestrk3MassHypo, &b_DRestrk3MassHypo);
   fChain->SetBranchAddress("DRestrk4MassHypo", DRestrk4MassHypo, &b_DRestrk4MassHypo);
   fChain->SetBranchAddress("DRestrk1Dz", DRestrk1Dz, &b_DRestrk1Dz);
   fChain->SetBranchAddress("DRestrk2Dz", DRestrk2Dz, &b_DRestrk2Dz);
   fChain->SetBranchAddress("DRestrk3Dz", DRestrk3Dz, &b_DRestrk3Dz);
   fChain->SetBranchAddress("DRestrk4Dz", DRestrk4Dz, &b_DRestrk4Dz);
   fChain->SetBranchAddress("DRestrk1Dxy", DRestrk1Dxy, &b_DRestrk1Dxy);
   fChain->SetBranchAddress("DRestrk2Dxy", DRestrk2Dxy, &b_DRestrk2Dxy);
   fChain->SetBranchAddress("DRestrk3Dxy", DRestrk3Dxy, &b_DRestrk3Dxy);
   fChain->SetBranchAddress("DRestrk4Dxy", DRestrk4Dxy, &b_DRestrk4Dxy);
   fChain->SetBranchAddress("DRestrk1originalAlgo", DRestrk1originalAlgo, &b_DRestrk1originalAlgo);
   fChain->SetBranchAddress("DRestrk2originalAlgo", DRestrk2originalAlgo, &b_DRestrk2originalAlgo);
   fChain->SetBranchAddress("DRestrk3originalAlgo", DRestrk3originalAlgo, &b_DRestrk3originalAlgo);
   fChain->SetBranchAddress("DRestrk4originalAlgo", DRestrk4originalAlgo, &b_DRestrk4originalAlgo);
   fChain->SetBranchAddress("DRestrk1highPurity", DRestrk1highPurity, &b_DRestrk1highPurity);
   fChain->SetBranchAddress("DRestrk2highPurity", DRestrk2highPurity, &b_DRestrk2highPurity);
   fChain->SetBranchAddress("DRestrk3highPurity", DRestrk3highPurity, &b_DRestrk3highPurity);
   fChain->SetBranchAddress("DRestrk4highPurity", DRestrk4highPurity, &b_DRestrk4highPurity);
   fChain->SetBranchAddress("DRestrk1dedx", DRestrk1dedx, &b_DRestrk1dedx);
   fChain->SetBranchAddress("DRestrk2dedx", DRestrk2dedx, &b_DRestrk2dedx);
   fChain->SetBranchAddress("DRestrk3dedx", DRestrk3dedx, &b_DRestrk3dedx);
   fChain->SetBranchAddress("DRestrk4dedx", DRestrk4dedx, &b_DRestrk4dedx);
   fChain->SetBranchAddress("DRestrk1thetastar", DRestrk1thetastar, &b_DRestrk1thetastar);
   fChain->SetBranchAddress("DRestrk2thetastar", DRestrk2thetastar, &b_DRestrk2thetastar);
   fChain->SetBranchAddress("DRestrk3thetastar", DRestrk3thetastar, &b_DRestrk3thetastar);
   fChain->SetBranchAddress("DRestrk4thetastar", DRestrk4thetastar, &b_DRestrk4thetastar);
   fChain->SetBranchAddress("DRestrk1thetastar_uf", DRestrk1thetastar_uf, &b_DRestrk1thetastar_uf);
   fChain->SetBranchAddress("DRestrk2thetastar_uf", DRestrk2thetastar_uf, &b_DRestrk2thetastar_uf);
   fChain->SetBranchAddress("DRestrk3thetastar_uf", DRestrk3thetastar_uf, &b_DRestrk3thetastar_uf);
   fChain->SetBranchAddress("DRestrk4thetastar_uf", DRestrk4thetastar_uf, &b_DRestrk4thetastar_uf);
   if(detailMode)
	 {	 
   fChain->SetBranchAddress("DRestrk1Y", DRestrk1Y, &b_DRestrk1Y);
   fChain->SetBranchAddress("DRestrk2Y", DRestrk2Y, &b_DRestrk2Y);
   fChain->SetBranchAddress("DRestrk3Y", DRestrk3Y, &b_DRestrk3Y);
   fChain->SetBranchAddress("DRestrk4Y", DRestrk4Y, &b_DRestrk4Y);
   fChain->SetBranchAddress("DRestrk1D0", DRestrk1D0, &b_DRestrk1D0);
   fChain->SetBranchAddress("DRestrk2D0", DRestrk2D0, &b_DRestrk2D0);
   fChain->SetBranchAddress("DRestrk3D0", DRestrk3D0, &b_DRestrk3D0);
   fChain->SetBranchAddress("DRestrk4D0", DRestrk4D0, &b_DRestrk4D0);
   fChain->SetBranchAddress("DRestrk1D0Err", DRestrk1D0Err, &b_DRestrk1D0Err);
   fChain->SetBranchAddress("DRestrk2D0Err", DRestrk2D0Err, &b_DRestrk2D0Err);
   fChain->SetBranchAddress("DRestrk3D0Err", DRestrk3D0Err, &b_DRestrk3D0Err);
   fChain->SetBranchAddress("DRestrk4D0Err", DRestrk4D0Err, &b_DRestrk4D0Err);
   fChain->SetBranchAddress("DRestrk1Quality", DRestrk1Quality, &b_DRestrk1Quality);
   fChain->SetBranchAddress("DRestrk2Quality", DRestrk2Quality, &b_DRestrk2Quality);
   fChain->SetBranchAddress("DRestrk3Quality", DRestrk3Quality, &b_DRestrk3Quality);
   fChain->SetBranchAddress("DRestrk4Quality", DRestrk4Quality, &b_DRestrk4Quality);
	 }
	 if(!isData){
   fChain->SetBranchAddress("Dgen", Dgen, &b_Dgen);
   fChain->SetBranchAddress("DsGen", DsGen, &b_DsGen);
   fChain->SetBranchAddress("DgenIndex", DgenIndex, &b_DgenIndex);
   fChain->SetBranchAddress("DgennDa", DgennDa, &b_DgennDa);
   fChain->SetBranchAddress("Dgenpt", Dgenpt, &b_Dgenpt);
   fChain->SetBranchAddress("Dgeneta", Dgeneta, &b_Dgeneta);
   fChain->SetBranchAddress("Dgenphi", Dgenphi, &b_Dgenphi);
   fChain->SetBranchAddress("Dgeny", Dgeny, &b_Dgeny);
   fChain->SetBranchAddress("DgencollisionId", DgencollisionId, &b_DgencollisionId);
   fChain->SetBranchAddress("DgenBAncestorpt", DgenBAncestorpt, &b_DgenBAncestorpt);
   fChain->SetBranchAddress("DgenBAncestorpdgId", DgenBAncestorpdgId, &b_DgenBAncestorpdgId);
   fChain->SetBranchAddress("DgenprodvtxX", DgenprodvtxX, &b_DgenprodvtxX);
   fChain->SetBranchAddress("DgenprodvtxY", DgenprodvtxY, &b_DgenprodvtxY);
   fChain->SetBranchAddress("DgenprodvtxZ", DgenprodvtxZ, &b_DgenprodvtxZ);
   fChain->SetBranchAddress("DgendecayvtxX", DgendecayvtxX, &b_DgendecayvtxX);
   fChain->SetBranchAddress("DgendecayvtxY", DgendecayvtxY, &b_DgendecayvtxY);
   fChain->SetBranchAddress("DgendecayvtxZ", DgendecayvtxZ, &b_DgendecayvtxZ);
   fChain->SetBranchAddress("DgendecayvtxXtoPv", DgendecayvtxXtoPv, &b_DgendecayvtxXtoPv);
   fChain->SetBranchAddress("DgendecayvtxYtoPv", DgendecayvtxYtoPv, &b_DgendecayvtxYtoPv);
   fChain->SetBranchAddress("DgendecayvtxZtoPv", DgendecayvtxZtoPv, &b_DgendecayvtxZtoPv);
   fChain->SetBranchAddress("DgenfromgenPV", DgenfromgenPV, &b_DgenfromgenPV);
	 } 

// remove later
//  if(hasFonllDptWeight){
//  nt->SetBranchAddress("FonllDptWeight",FonllDptWeight);
//  }

}


/*
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
*/

// Hlt Tree

/*
Int_t     HLT_HIL1MinimumBiasHF1AND_v1;
Int_t     HLT_HIL1MinimumBiasHF2AND_v1;
Int_t     HLT_HIL1MinimumBiasHF2AND_part1_v1;
Int_t     HLT_HIL1MinimumBiasHF2AND_part2_v1;
Int_t     HLT_HIL1MinimumBiasHF2AND_part3_v1;
Int_t     HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part1_v1;
Int_t     HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part2_v1;
Int_t     HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part3_v1;
Int_t     HLT_HIL1Centralityext50100MinimumumBiasHF2AND_v1;
*/


// for PbPb
   // Int_t           HLT_HIL1MinimumBiasHF1OR_v1;
   // Int_t           HLT_HIL1MinimumBiasHF2OR_v1;
   // Int_t           HLT_HIL1MinimumBiasHF1AND_v1;
   Int_t           HLT_HIL1MinimumBiasHF2AND_part1_v1;
   Int_t           HLT_HIL1MinimumBiasHF2AND_part2_v1;
   Int_t           HLT_HIL1MinimumBiasHF2AND_part3_v1; // Golden json only use 123, track only has 1-11
   Int_t           HLT_HIL1MinimumBiasHF2AND_part4_v1;
   Int_t           HLT_HIL1MinimumBiasHF2AND_part5_v1;
   Int_t           HLT_HIL1MinimumBiasHF2AND_part6_v1;
   Int_t           HLT_HIL1MinimumBiasHF2AND_part7_v1;
   Int_t           HLT_HIL1MinimumBiasHF2AND_part8_v1;
   Int_t           HLT_HIL1MinimumBiasHF2AND_part9_v1;
   Int_t           HLT_HIL1MinimumBiasHF2AND_part10_v1;
   Int_t           HLT_HIL1MinimumBiasHF2AND_part11_v1;

   // Int_t           HLT_HIL1Centralityext30100MinimumumBiasHF1AND_v1;  // for HIMB567, PbPb 30-100 trigger
   Int_t           HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part1_v1;
   Int_t           HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part2_v1;
   Int_t           HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part3_v1;

	// for pp
	 Int_t           HLT_L1MinimumBiasHF1OR_part0_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part1_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part2_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part3_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part4_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part5_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part6_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part7_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part8_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part9_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part10_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part11_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part12_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part13_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part14_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part15_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part16_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part17_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part18_v1;
	 Int_t           HLT_L1MinimumBiasHF1OR_part19_v1;

void SetHltBranches(TTree* nt, Bool_t isPbPb)
{
	if( isPbPb){
		nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part1_v1",&HLT_HIL1MinimumBiasHF2AND_part1_v1);
		nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part2_v1",&HLT_HIL1MinimumBiasHF2AND_part2_v1);
		nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part3_v1",&HLT_HIL1MinimumBiasHF2AND_part3_v1);
		nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part4_v1",&HLT_HIL1MinimumBiasHF2AND_part4_v1);
		nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part5_v1",&HLT_HIL1MinimumBiasHF2AND_part5_v1);
		nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part6_v1",&HLT_HIL1MinimumBiasHF2AND_part6_v1);
		nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part7_v1",&HLT_HIL1MinimumBiasHF2AND_part7_v1);
		nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part8_v1",&HLT_HIL1MinimumBiasHF2AND_part8_v1);
		nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part9_v1",&HLT_HIL1MinimumBiasHF2AND_part9_v1);
		nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part10_v1",&HLT_HIL1MinimumBiasHF2AND_part10_v1);
		nt->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part11_v1",&HLT_HIL1MinimumBiasHF2AND_part11_v1);

		nt->SetBranchAddress("HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part1_v1",&HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part1_v1);
		nt->SetBranchAddress("HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part2_v1",&HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part2_v1);
		nt->SetBranchAddress("HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part3_v1",&HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part3_v1);
	}
	else{
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part0_v1",&HLT_L1MinimumBiasHF1OR_part0_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part1_v1",&HLT_L1MinimumBiasHF1OR_part1_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part2_v1",&HLT_L1MinimumBiasHF1OR_part2_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part3_v1",&HLT_L1MinimumBiasHF1OR_part3_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part4_v1",&HLT_L1MinimumBiasHF1OR_part4_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part5_v1",&HLT_L1MinimumBiasHF1OR_part5_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part6_v1",&HLT_L1MinimumBiasHF1OR_part6_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part7_v1",&HLT_L1MinimumBiasHF1OR_part7_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part8_v1",&HLT_L1MinimumBiasHF1OR_part8_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part9_v1",&HLT_L1MinimumBiasHF1OR_part9_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part10_v1",&HLT_L1MinimumBiasHF1OR_part10_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part11_v1",&HLT_L1MinimumBiasHF1OR_part11_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part12_v1",&HLT_L1MinimumBiasHF1OR_part12_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part13_v1",&HLT_L1MinimumBiasHF1OR_part13_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part14_v1",&HLT_L1MinimumBiasHF1OR_part14_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part15_v1",&HLT_L1MinimumBiasHF1OR_part15_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part16_v1",&HLT_L1MinimumBiasHF1OR_part16_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part17_v1",&HLT_L1MinimumBiasHF1OR_part17_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part18_v1",&HLT_L1MinimumBiasHF1OR_part18_v1);	
		nt->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part19_v1",&HLT_L1MinimumBiasHF1OR_part19_v1);	
	}
}

/*
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
*/



Int_t     hiBin;
// Int_t     hiNevtPlane;
// Float_t   hiEvtPlanes[kMaxEvtPlanes];
// Float_t   hiEvtPlanesqx[kMaxEvtPlanes];
// Float_t   hiEvtPlanesqy[kMaxEvtPlanes];
// Float_t   pthatweight;
Float_t   pthat;
Float_t   weight;
Float_t   vz;
Float_t   vx;
Float_t   vy;
Float_t   Npart; // only PbPb
Float_t   Ncoll; // only PbPb
// Int_t     DPtSample;

void SetHiBranches(TTree* nt, Bool_t isData=true, Int_t isPbPb=0)
{
  if(!isData){
		nt->SetBranchAddress("pthat",&pthat);
		nt->SetBranchAddress("weight",&weight);
		if(isPbPb){
			nt->SetBranchAddress("Npart",&Npart);
			nt->SetBranchAddress("Ncoll",&Ncoll);
		}
	}
		nt->SetBranchAddress("hiBin",&hiBin);
		nt->SetBranchAddress("vx",&vx);
		nt->SetBranchAddress("vy",&vy);
		nt->SetBranchAddress("vz",&vz);
		// nt->SetBranchAddress("hiNevtPlane",&hiNevtPlane);
		// nt->SetBranchAddress("hiEvtPlanes",hiEvtPlanes);
		// nt->SetBranchAddress("hiEvtPlanesqx",hiEvtPlanesqx);
		// nt->SetBranchAddress("hiEvtPlanesqy",hiEvtPlanesqy);
}

/*
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
*/


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


/*
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
*/


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




void SetGenBranches(TTree* nt, Bool_t SelectBranchOnly=false, Bool_t DsGenTree=true)
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

	nt->SetBranchAddress("GRestk1pt",    GRestk1pt		);
  nt->SetBranchAddress("GRestk1eta",   GRestk1eta);
  nt->SetBranchAddress("GRestk1y",     GRestk1y);
  nt->SetBranchAddress("GRestk1phi",   GRestk1phi);
  nt->SetBranchAddress("GRestk2pt",    GRestk2pt);
  nt->SetBranchAddress("GRestk2eta",   GRestk2eta);
  nt->SetBranchAddress("GRestk2y",     GRestk2y);
  nt->SetBranchAddress("GRestk2phi",   GRestk2phi);

	if(DsGenTree){
  nt->SetBranchAddress("GSignalType", GSignalType);
	}

	if(SelectBranchOnly)
	{
		nt->SetBranchStatus("*",0);
    nt->SetBranchStatus("GPVx", 1);
    nt->SetBranchStatus("GPVy", 1);
    nt->SetBranchStatus("GPVz", 1);
    nt->SetBranchStatus("Gsize", 1);
    nt->SetBranchStatus("Gy", 1);
    nt->SetBranchStatus("Geta", 1);
    nt->SetBranchStatus("Gphi", 1);
    nt->SetBranchStatus("Gpt", 1);
    nt->SetBranchStatus("GpdgId",1 );
    nt->SetBranchStatus("GcollisionId",1 );
    nt->SetBranchStatus("GisSignal",1 );
    // nt->SetBranchStatus("GSignalType",1 );
    nt->SetBranchStatus("GBAncestorpt",1 );
    nt->SetBranchStatus("GBAncestorpdgId",1 );
    nt->SetBranchStatus("Gtk1pt",1 );
    nt->SetBranchStatus("Gtk1eta",1 );
    // nt->SetBranchStatus("Gtk1y", 1);
    // nt->SetBranchStatus("Gtk1phi", 1);
    nt->SetBranchStatus("Gtk2pt", 1);
    nt->SetBranchStatus("Gtk2eta", 1);
    // nt->SetBranchStatus("Gtk2y", 1);
    // nt->SetBranchStatus("Gtk2phi", 1);

	  nt->SetBranchStatus("GRestk1pt",1   );
    nt->SetBranchStatus("GRestk1eta",1  );
    // nt->SetBranchStatus("GRestk1y",1 );
    // nt->SetBranchStatus("GRestk1phi",1  );
    nt->SetBranchStatus("GRestk2pt",1  );
    nt->SetBranchStatus("GRestk2eta",1  );
    // nt->SetBranchStatus("GRestk2y",1  );
    // nt->SetBranchStatus("GRestk2phi",1  );
		if(DsGenTree){
    nt->SetBranchStatus("GSignalType",1 );
		}
	
	}

/*
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
*/

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
