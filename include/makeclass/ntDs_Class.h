//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb 26 17:53:40 2018 by ROOT version 6.02/13
// from TTree ntDPhikkpi/
// found on file: Dntuple_finder_PbPb_1111.root_1104.root
//////////////////////////////////////////////////////////

#ifndef ntDs_Class_h
#define ntDs_Class_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ntDs_Class {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
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
   Int_t           Dindex[280];   //[Dsize]
   Int_t           Dtype[280];   //[Dsize]
   Float_t         Dmass[280];   //[Dsize]
   Float_t         D_unfitted_mass[280];   //[Dsize]
   Float_t         Dpt[280];   //[Dsize]
   Float_t         D_unfitted_pt[280];   //[Dsize]
   Float_t         Deta[280];   //[Dsize]
   Float_t         Dphi[280];   //[Dsize]
   Float_t         Dy[280];   //[Dsize]
   Float_t         DvtxX[280];   //[Dsize]
   Float_t         DvtxY[280];   //[Dsize]
   Float_t         DvtxZ[280];   //[Dsize]
   Float_t         Dd0[280];   //[Dsize]
   Float_t         Dd0Err[280];   //[Dsize]
   Float_t         Ddxyz[280];   //[Dsize]
   Float_t         DdxyzErr[280];   //[Dsize]
   Float_t         Dchi2ndf[280];   //[Dsize]
   Float_t         Dchi2cl[280];   //[Dsize]
   Float_t         Ddtheta[280];   //[Dsize]
   Float_t         Dlxy[280];   //[Dsize]
   Float_t         Dalpha[280];   //[Dsize]
   Float_t         DsvpvDistance[280];   //[Dsize]
   Float_t         DsvpvDisErr[280];   //[Dsize]
   Float_t         DsvpvDistance_2D[280];   //[Dsize]
   Float_t         DsvpvDisErr_2D[280];   //[Dsize]
   Float_t         Ddca[280];   //[Dsize]
   Float_t         DlxyBS[280];   //[Dsize]
   Float_t         DlxyBSErr[280];   //[Dsize]
   Float_t         DMaxDoca[280];   //[Dsize]
   Float_t         DMaxTkPt[280];   //[Dsize]
   Float_t         DMinTkPt[280];   //[Dsize]
   Float_t         Dtrk1Pt[280];   //[Dsize]
   Float_t         Dtrk2Pt[280];   //[Dsize]
   Float_t         Dtrk1PtErr[280];   //[Dsize]
   Float_t         Dtrk2PtErr[280];   //[Dsize]
   Float_t         Dtrk1Eta[280];   //[Dsize]
   Float_t         Dtrk2Eta[280];   //[Dsize]
   Float_t         Dtrk1Phi[280];   //[Dsize]
   Float_t         Dtrk2Phi[280];   //[Dsize]
   Float_t         Dtrk1P[280];   //[Dsize]
   Float_t         Dtrk2P[280];   //[Dsize]
   Float_t         Dtrk1Dz[280];   //[Dsize]
   Float_t         Dtrk2Dz[280];   //[Dsize]
   Float_t         Dtrk1Dxy[280];   //[Dsize]
   Float_t         Dtrk2Dxy[280];   //[Dsize]
   Float_t         Dtrk1MassHypo[280];   //[Dsize]
   Float_t         Dtrk2MassHypo[280];   //[Dsize]
   Int_t           Dtrk1originalAlgo[280];   //[Dsize]
   Int_t           Dtrk2originalAlgo[280];   //[Dsize]
   Bool_t          Dtrk1highPurity[280];   //[Dsize]
   Bool_t          Dtrk2highPurity[280];   //[Dsize]
   Float_t         Dtrk1dedx[280];   //[Dsize]
   Float_t         Dtrk2dedx[280];   //[Dsize]
   Float_t         Dtrk1thetastar[280];   //[Dsize]
   Float_t         Dtrk2thetastar[280];   //[Dsize]
   Float_t         Dtrk1thetastar_uf[280];   //[Dsize]
   Float_t         Dtrk2thetastar_uf[280];   //[Dsize]
   Float_t         Dtrk1PixelHit[280];   //[Dsize]
   Float_t         Dtrk2PixelHit[280];   //[Dsize]
   Float_t         Dtrk1StripHit[280];   //[Dsize]
   Float_t         Dtrk2StripHit[280];   //[Dsize]
   Float_t         Dtrk1nStripLayer[280];   //[Dsize]
   Float_t         Dtrk2nStripLayer[280];   //[Dsize]
   Float_t         Dtrk1nPixelLayer[280];   //[Dsize]
   Float_t         Dtrk2nPixelLayer[280];   //[Dsize]
   Float_t         Dtrk1Chi2ndf[280];   //[Dsize]
   Float_t         Dtrk2Chi2ndf[280];   //[Dsize]
   Int_t           Dtrk1Algo[280];   //[Dsize]
   Int_t           Dtrk2Algo[280];   //[Dsize]
   Float_t         Dtrk3Pt[280];   //[Dsize]
   Float_t         Dtrk4Pt[280];   //[Dsize]
   Float_t         Dtrk3PtErr[280];   //[Dsize]
   Float_t         Dtrk4PtErr[280];   //[Dsize]
   Float_t         Dtrk3Eta[280];   //[Dsize]
   Float_t         Dtrk4Eta[280];   //[Dsize]
   Float_t         Dtrk3Phi[280];   //[Dsize]
   Float_t         Dtrk4Phi[280];   //[Dsize]
   Float_t         Dtrk3P[280];   //[Dsize]
   Float_t         Dtrk4P[280];   //[Dsize]
   Float_t         Dtrk3MassHypo[280];   //[Dsize]
   Float_t         Dtrk4MassHypo[280];   //[Dsize]
   Float_t         Dtrk3Dz[280];   //[Dsize]
   Float_t         Dtrk4Dz[280];   //[Dsize]
   Float_t         Dtrk3Dxy[280];   //[Dsize]
   Float_t         Dtrk4Dxy[280];   //[Dsize]
   Int_t           Dtrk3originalAlgo[280];   //[Dsize]
   Int_t           Dtrk4originalAlgo[280];   //[Dsize]
   Bool_t          Dtrk3highPurity[280];   //[Dsize]
   Bool_t          Dtrk4highPurity[280];   //[Dsize]
   Float_t         Dtrk3dedx[280];   //[Dsize]
   Float_t         Dtrk4dedx[280];   //[Dsize]
   Float_t         Dtrk3thetastar[280];   //[Dsize]
   Float_t         Dtrk4thetastar[280];   //[Dsize]
   Float_t         Dtrk3thetastar_uf[280];   //[Dsize]
   Float_t         Dtrk4thetastar_uf[280];   //[Dsize]
   Float_t         Dtrk3PixelHit[280];   //[Dsize]
   Float_t         Dtrk4PixelHit[280];   //[Dsize]
   Float_t         Dtrk3StripHit[280];   //[Dsize]
   Float_t         Dtrk4StripHit[280];   //[Dsize]
   Float_t         Dtrk3nStripLayer[280];   //[Dsize]
   Float_t         Dtrk4nStripLayer[280];   //[Dsize]
   Float_t         Dtrk3nPixelLayer[280];   //[Dsize]
   Float_t         Dtrk4nPixelLayer[280];   //[Dsize]
   Float_t         Dtrk3Chi2ndf[280];   //[Dsize]
   Float_t         Dtrk4Chi2ndf[280];   //[Dsize]
   Int_t           Dtrk3Algo[280];   //[Dsize]
   Int_t           Dtrk4Algo[280];   //[Dsize]
   Float_t         Dtrk1Y[280];   //[Dsize]
   Float_t         Dtrk2Y[280];   //[Dsize]
   Float_t         Dtrk1D0[280];   //[Dsize]
   Float_t         Dtrk2D0[280];   //[Dsize]
   Float_t         Dtrk1D0Err[280];   //[Dsize]
   Float_t         Dtrk2D0Err[280];   //[Dsize]
   Int_t           Dtrk1Quality[280];   //[Dsize]
   Int_t           Dtrk2Quality[280];   //[Dsize]
   Int_t           Dtrk1Idx[280];   //[Dsize]
   Int_t           Dtrk2Idx[280];   //[Dsize]
   Float_t         Dtrk1EtaErr[280];   //[Dsize]
   Float_t         Dtrk2EtaErr[280];   //[Dsize]
   Float_t         Dtrk1PhiErr[280];   //[Dsize]
   Float_t         Dtrk2PhiErr[280];   //[Dsize]
   Float_t         Dtrk1MVAVal[280];   //[Dsize]
   Float_t         Dtrk2MVAVal[280];   //[Dsize]
   Float_t         Dtrk3Y[280];   //[Dsize]
   Float_t         Dtrk4Y[280];   //[Dsize]
   Float_t         Dtrk3D0[280];   //[Dsize]
   Float_t         Dtrk4D0[280];   //[Dsize]
   Float_t         Dtrk3D0Err[280];   //[Dsize]
   Float_t         Dtrk4D0Err[280];   //[Dsize]
   Int_t           Dtrk3Quality[280];   //[Dsize]
   Int_t           Dtrk4Quality[280];   //[Dsize]
   Int_t           Dtrk3Idx[280];   //[Dsize]
   Int_t           Dtrk4Idx[280];   //[Dsize]
   Float_t         Dtrk3EtaErr[280];   //[Dsize]
   Float_t         Dtrk4EtaErr[280];   //[Dsize]
   Float_t         Dtrk3PhiErr[280];   //[Dsize]
   Float_t         Dtrk4PhiErr[280];   //[Dsize]
   Float_t         Dtrk3MVAVal[280];   //[Dsize]
   Float_t         Dtrk4MVAVal[280];   //[Dsize]
   Float_t         DtktkResmass[280];   //[Dsize]
   Float_t         DtktkRes_unfitted_mass[280];   //[Dsize]
   Float_t         DtktkRespt[280];   //[Dsize]
   Float_t         DtktkRes_unfitted_pt[280];   //[Dsize]
   Float_t         DtktkReseta[280];   //[Dsize]
   Float_t         DtktkResphi[280];   //[Dsize]
   Float_t         DtktkRes_chi2ndf[280];   //[Dsize]
   Float_t         DtktkRes_chi2cl[280];   //[Dsize]
   Float_t         DtktkRes_alpha[280];   //[Dsize]
   Float_t         DtktkRes_alphaToSV[280];   //[Dsize]
   Float_t         DtktkRes_svpvDistance[280];   //[Dsize]
   Float_t         DtktkRes_svpvDisErr[280];   //[Dsize]
   Float_t         DtktkRes_svpvDistanceToSV[280];   //[Dsize]
   Float_t         DtktkRes_svpvDisErrToSV[280];   //[Dsize]
   Float_t         DtktkRes_dca[280];   //[Dsize]
   Float_t         DtktkRes_dcaToSV[280];   //[Dsize]
   Float_t         DtktkRes_lxyBS[280];   //[Dsize]
   Float_t         DtktkRes_lxyBSErr[280];   //[Dsize]
   Float_t         DtktkRes_angleToTrk1[280];   //[Dsize]
   Float_t         DtktkRes_unfitted_angleToTrk1[280];   //[Dsize]
   Float_t         DtktkRes_ptAsymToTrk1[280];   //[Dsize]
   Float_t         DtktkRes_unfitter_ptAsymToTrk1[280];   //[Dsize]
   Float_t         DRestrk1Pt[280];   //[Dsize]
   Float_t         DRestrk2Pt[280];   //[Dsize]
   Float_t         DRestrk3Pt[280];   //[Dsize]
   Float_t         DRestrk4Pt[280];   //[Dsize]
   Float_t         DRestrk1PtErr[280];   //[Dsize]
   Float_t         DRestrk2PtErr[280];   //[Dsize]
   Float_t         DRestrk3PtErr[280];   //[Dsize]
   Float_t         DRestrk4PtErr[280];   //[Dsize]
   Float_t         DRestrk1Eta[280];   //[Dsize]
   Float_t         DRestrk2Eta[280];   //[Dsize]
   Float_t         DRestrk3Eta[280];   //[Dsize]
   Float_t         DRestrk4Eta[280];   //[Dsize]
   Float_t         DRestrk1Phi[280];   //[Dsize]
   Float_t         DRestrk2Phi[280];   //[Dsize]
   Float_t         DRestrk3Phi[280];   //[Dsize]
   Float_t         DRestrk4Phi[280];   //[Dsize]
   Float_t         DRestrk1P[280];   //[Dsize]
   Float_t         DRestrk2P[280];   //[Dsize]
   Float_t         DRestrk3P[280];   //[Dsize]
   Float_t         DRestrk4P[280];   //[Dsize]
   Float_t         DRestrk1MassHypo[280];   //[Dsize]
   Float_t         DRestrk2MassHypo[280];   //[Dsize]
   Float_t         DRestrk3MassHypo[280];   //[Dsize]
   Float_t         DRestrk4MassHypo[280];   //[Dsize]
   Float_t         DRestrk1Dz[280];   //[Dsize]
   Float_t         DRestrk2Dz[280];   //[Dsize]
   Float_t         DRestrk3Dz[280];   //[Dsize]
   Float_t         DRestrk4Dz[280];   //[Dsize]
   Float_t         DRestrk1Dxy[280];   //[Dsize]
   Float_t         DRestrk2Dxy[280];   //[Dsize]
   Float_t         DRestrk3Dxy[280];   //[Dsize]
   Float_t         DRestrk4Dxy[280];   //[Dsize]
   Int_t           DRestrk1originalAlgo[280];   //[Dsize]
   Int_t           DRestrk2originalAlgo[280];   //[Dsize]
   Int_t           DRestrk3originalAlgo[280];   //[Dsize]
   Int_t           DRestrk4originalAlgo[280];   //[Dsize]
   Bool_t          DRestrk1highPurity[280];   //[Dsize]
   Bool_t          DRestrk2highPurity[280];   //[Dsize]
   Bool_t          DRestrk3highPurity[280];   //[Dsize]
   Bool_t          DRestrk4highPurity[280];   //[Dsize]
   Float_t         DRestrk1dedx[280];   //[Dsize]
   Float_t         DRestrk2dedx[280];   //[Dsize]
   Float_t         DRestrk3dedx[280];   //[Dsize]
   Float_t         DRestrk4dedx[280];   //[Dsize]
   Float_t         DRestrk1thetastar[280];   //[Dsize]
   Float_t         DRestrk2thetastar[280];   //[Dsize]
   Float_t         DRestrk3thetastar[280];   //[Dsize]
   Float_t         DRestrk4thetastar[280];   //[Dsize]
   Float_t         DRestrk1thetastar_uf[280];   //[Dsize]
   Float_t         DRestrk2thetastar_uf[280];   //[Dsize]
   Float_t         DRestrk3thetastar_uf[280];   //[Dsize]
   Float_t         DRestrk4thetastar_uf[280];   //[Dsize]
   Float_t         DRestrk1Y[280];   //[Dsize]
   Float_t         DRestrk2Y[280];   //[Dsize]
   Float_t         DRestrk3Y[280];   //[Dsize]
   Float_t         DRestrk4Y[280];   //[Dsize]
   Float_t         DRestrk1D0[280];   //[Dsize]
   Float_t         DRestrk2D0[280];   //[Dsize]
   Float_t         DRestrk3D0[280];   //[Dsize]
   Float_t         DRestrk4D0[280];   //[Dsize]
   Float_t         DRestrk1D0Err[280];   //[Dsize]
   Float_t         DRestrk2D0Err[280];   //[Dsize]
   Float_t         DRestrk3D0Err[280];   //[Dsize]
   Float_t         DRestrk4D0Err[280];   //[Dsize]
   Int_t           DRestrk1Quality[280];   //[Dsize]
   Int_t           DRestrk2Quality[280];   //[Dsize]
   Int_t           DRestrk3Quality[280];   //[Dsize]
   Int_t           DRestrk4Quality[280];   //[Dsize]
   Float_t         Dgen[280];   //[Dsize]
   Int_t           DgenIndex[280];   //[Dsize]
   Int_t           DgennDa[280];   //[Dsize]
   Float_t         Dgenpt[280];   //[Dsize]
   Float_t         Dgeneta[280];   //[Dsize]
   Float_t         Dgenphi[280];   //[Dsize]
   Float_t         Dgeny[280];   //[Dsize]
   Int_t           DgencollisionId[280];   //[Dsize]
   Float_t         DgenBAncestorpt[280];   //[Dsize]
   Int_t           DgenBAncestorpdgId[280];   //[Dsize]
   Float_t         DgenprodvtxX[280];   //[Dsize]
   Float_t         DgenprodvtxY[280];   //[Dsize]
   Float_t         DgenprodvtxZ[280];   //[Dsize]
   Float_t         DgendecayvtxX[280];   //[Dsize]
   Float_t         DgendecayvtxY[280];   //[Dsize]
   Float_t         DgendecayvtxZ[280];   //[Dsize]
   Int_t           DgenfromgenPV[280];   //[Dsize]

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
   TBranch        *b_Dd0;   //!
   TBranch        *b_Dd0Err;   //!
   TBranch        *b_Ddxyz;   //!
   TBranch        *b_DdxyzErr;   //!
   TBranch        *b_Dchi2ndf;   //!
   TBranch        *b_Dchi2cl;   //!
   TBranch        *b_Ddtheta;   //!
   TBranch        *b_Dlxy;   //!
   TBranch        *b_Dalpha;   //!
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
   TBranch        *b_DgenfromgenPV;   //!

   ntDs_Class(TTree *tree=0);
   virtual ~ntDs_Class();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ntDs_Class_cxx
ntDs_Class::ntDs_Class(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Dntuple_finder_PbPb_1111.root_1104.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Dntuple_finder_PbPb_1111.root_1104.root");
      }
      f->GetObject("ntDPhikkpi",tree);

   }
   Init(tree);
}

ntDs_Class::~ntDs_Class()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ntDs_Class::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ntDs_Class::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ntDs_Class::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

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
   fChain->SetBranchAddress("Dd0", Dd0, &b_Dd0);
   fChain->SetBranchAddress("Dd0Err", Dd0Err, &b_Dd0Err);
   fChain->SetBranchAddress("Ddxyz", Ddxyz, &b_Ddxyz);
   fChain->SetBranchAddress("DdxyzErr", DdxyzErr, &b_DdxyzErr);
   fChain->SetBranchAddress("Dchi2ndf", Dchi2ndf, &b_Dchi2ndf);
   fChain->SetBranchAddress("Dchi2cl", Dchi2cl, &b_Dchi2cl);
   fChain->SetBranchAddress("Ddtheta", Ddtheta, &b_Ddtheta);
   fChain->SetBranchAddress("Dlxy", Dlxy, &b_Dlxy);
   fChain->SetBranchAddress("Dalpha", Dalpha, &b_Dalpha);
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
   fChain->SetBranchAddress("Dgen", Dgen, &b_Dgen);
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
   fChain->SetBranchAddress("DgenfromgenPV", DgenfromgenPV, &b_DgenfromgenPV);
   Notify();
}

Bool_t ntDs_Class::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ntDs_Class::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ntDs_Class::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ntDs_Class_cxx
