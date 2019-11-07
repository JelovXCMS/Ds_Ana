#!/bin/bash

DOSAVEHISTPP_ALPHA=0
DOSAVEHISTPbPb_ALPHA=0
DOSAVEHISTPPMB_ALPHA=0
DOSAVEHISTPbPbMB_ALPHA=0
DOSAVEHISTPPMBLowPt_ALPHA=0
DOSAVEHISTPbPbMBLowPt_ALPHA=0

DOSAVEHISTPP_DECAYLENGTH=0
DOSAVEHISTPbPb_DECAYLENGTH=0
DOSAVEHISTPPMB_DECAYLENGTH=0
DOSAVEHISTPbPbMB_DECAYLENGTH=0
DOSAVEHISTPPMBLowPt_DECAYLENGTH=0
DOSAVEHISTPbPbMBLowPt_DECAYLENGTH=0

DOSAVEHISTPP_CHI2=0
DOSAVEHISTPbPb_CHI2=0
DOSAVEHISTPPMB_CHI2=0
DOSAVEHISTPbPbMB_CHI2=0
DOSAVEHISTPPMBLowPt_CHI2=0
DOSAVEHISTPbPbMBLowPt_CHI2=0

DOSAVEHISTPP_TRKDCA=0
DOSAVEHISTPbPb_TRKDCA=0
DOSAVEHISTPPMB_TRKDCA=0
DOSAVEHISTPbPbMB_TRKDCA=0
DOSAVEHISTPPMBLowPt_TRKDCA=0
DOSAVEHISTPbPbMBLowPt_TRKDCA=0

#
DOFITPP_ALPHA=1
DOFITPbPb_ALPHA=1
DOFITPPMB_ALPHA=1
DOFITPbPbMB_ALPHA=1
DOFITPPMBLowPt_ALPHA=1
DOFITPbPbMBLowPt_ALPHA=1

DOFITPP_DECAYLENGTH=1
DOFITPbPb_DECAYLENGTH=1
DOFITPPMB_DECAYLENGTH=1
DOFITPbPbMB_DECAYLENGTH=1
DOFITPPMBLowPt_DECAYLENGTH=1
DOFITPbPbMBLowPt_DECAYLENGTH=1

DOFITPP_CHI2=1
DOFITPbPb_CHI2=1
DOFITPPMB_CHI2=1
DOFITPbPbMB_CHI2=1
DOFITPPMBLowPt_CHI2=1
DOFITPbPbMBLowPt_CHI2=1

DOFITPP_TRKDCA=1
DOFITPbPb_TRKDCA=1
DOFITPPMB_TRKDCA=1
DOFITPbPbMB_TRKDCA=1
DOFITPPMBLowPt_TRKDCA=1
DOFITPbPbMBLowPt_TRKDCA=1

#
DORATIO_ALPHA=0
DORATIO_DECAYLENGTH=0
DORATIO_CHI2=1
DORATIO_TRKDCA=0

###

LABEL_ALPHA="Alpha"
LABEL_DECAYLENGTH="Decaylength"
LABEL_CHI2="Chi2"
LABEL_TRKDCA="TrkDca"
VAR_ALPHA="Dalpha"
VAR_DECAYLENGTH="(DsvpvDistance/DsvpvDisErr)"
VAR_CHI2="Dchi2cl"
VAR_TRKDCA="min(sqrt(pow(Dtrk1Dxy,2)+pow(Dtrk1Dsz,2))/Dtrk1DxyErr,sqrt(pow(Dtrk2Dxy,2)+pow(Dtrk2Dsz,2))/Dtrk2DxyErr)"

VARBINS=5
VARMIN_ALPHA=0.16
VARMAX_ALPHA=0.2
VARMIN_DECAYLENGTH=4.
VARMAX_DECAYLENGTH=12.
VARMIN_CHI2=0.05
VARMAX_CHI2=0.45
VARMIN_TRKDCA=0
VARMAX_TRKDCA=6
VARSIGN_ALPHA=0
VARSIGN_DECAYLENGTH=1
VARSIGN_CHI2=1
VARSIGN_TRKDCA=1

###

INPUTDATAPP="/mnt/hadoop/store/user/hqiu/DfinderProduction_data/skimMerge/Dntuple_PP_DTriggers_skimEvents.root"
INPUTMCPP_P="/home/qiu110/myHadoop/DfinderProduction_MCMatchOnly/MCWeightMerge/DNtuple_PromptPP_selectMergeWeight.root"
INPUTMCPP_NP="/home/qiu110/myHadoop/DfinderProduction_MCMatchOnly/MCWeightMerge/DNtuple_PromptPP_selectMergeWeight.root"
WEIGHTPP_P="DPtSampleWeight*MeasuredDPtWeight"
WEIGHTPP_NP="DPtSampleWeight*FonllPtWeight*(0.55369+1.22995*log(DgenBAncestorpt)+-0.49766*log(DgenBAncestorpt)*log(DgenBAncestorpt)+0.05655*log(DgenBAncestorpt)*log(DgenBAncestorpt)*log(DgenBAncestorpt))"
FILTERCUTPP="abs(PVz)<15&&pBeamScrapingFilter&&pPAprimaryVertexFilter&&Dy>-1.&&Dy<1.&&Dtrk1highPurity&&Dtrk2highPurity&&Dtrk1Pt>2.0&&Dtrk2Pt>2.0&&(DlxyBS/DlxyBSErr)>2.5&&abs(Dtrk1Eta)<1.5&&abs(Dtrk2Eta)<1.5&&Dtrk1PtErr/Dtrk1Pt<0.3&&Dtrk2PtErr/Dtrk2Pt<0.3"
CENTRALCUTPP_ALPHA="Dalpha<0.2"
CENTRALCUTPP_DECAYLENGTH="(DsvpvDistance/DsvpvDisErr)>12"
CENTRALCUTPP_CHI2="Dchi2cl>0.05"
CENTRALCUTPP_TRKDCA="min(sqrt(pow(Dtrk1Dxy,2)+pow(Dtrk1Dsz,2))/Dtrk1DxyErr,sqrt(pow(Dtrk2Dxy,2)+pow(Dtrk2Dsz,2))/Dtrk2DxyErr)>4.5"
TRGPP="((HLT_DmesonPPTrackingGlobal_Dpt15_v1&&Dpt>20&&Dpt<40)||(HLT_DmesonPPTrackingGlobal_Dpt30_v1&&Dpt>40&&Dpt<60)||(HLT_DmesonPPTrackingGlobal_Dpt50_v1&&Dpt>60))"
LABELPP="PP"
PTMINPP=20
PTMAXPP=100

INPUTDATAPbPb="/mnt/hadoop/store/user/hqiu/DfinderProduction_data/skimMerge/Dntuple_HardProbes.root"
INPUTMCPbPb_P="/home/qiu110/myHadoop/DfinderProduction_MCMatchOnly/MCWeightMerge/DNtuple_PromptPbPb_selectMergeWeight.root"
INPUTMCPbPb_NP="/home/qiu110/myHadoop/DfinderProduction_MCMatchOnly/MCWeightMerge/DNtuple_NonPromptPbPb_selectMergeWeight.root"
WEIGHTPbPb_P="DPtSampleWeight*MeasuredDPtWeight*centralityNCollWeight*vertexZWeight"
WEIGHTPbPb_NP="DPtSampleWeight*FonllPtWeight*centralityNCollWeight*vertexZWeight*(1.43138+0.37192*log(DgenBAncestorpt)+-0.71770*log(DgenBAncestorpt)*log(DgenBAncestorpt)+0.23269*log(DgenBAncestorpt)*log(DgenBAncestorpt)*log(DgenBAncestorpt)+-0.02181*log(DgenBAncestorpt)*log(DgenBAncestorpt)*log(DgenBAncestorpt)*log(DgenBAncestorpt))"
FILTERCUTPbPb="pclusterCompatibilityFilter&&pprimaryVertexFilter&&phfCoincFilter3&&abs(PVz)<15&&PVzE<4&&Dy>-1.&&Dy<1.&&Dtrk1highPurity&&Dtrk2highPurity&&Dtrk1Pt>8.5&&Dtrk2Pt>8.5&&abs(Dtrk1Eta)<1.5&&abs(Dtrk2Eta)<1.5&&Dtrk1PtErr/Dtrk1Pt<0.3&&Dtrk2PtErr/Dtrk2Pt<0.3&&(DlxyBS/DlxyBSErr)>2.5"
CENTRALCUTPbPb_ALPHA="Dalpha<0.2"
CENTRALCUTPbPb_DECAYLENGTH="(DsvpvDistance/DsvpvDisErr)>8."
CENTRALCUTPbPb_CHI2="Dchi2cl>0.05"
CENTRALCUTPbPb_TRKDCA="min(sqrt(pow(Dtrk1Dxy,2)+pow(Dtrk1Dsz,2))/Dtrk1DxyErr,sqrt(pow(Dtrk2Dxy,2)+pow(Dtrk2Dsz,2))/Dtrk2DxyErr)>3"
TRGPbPb="run>262922&&((HLT_HIDmesonHITrackingGlobal_Dpt20_v2&&Dpt>20&&Dpt<40)||(HLT_HIDmesonHITrackingGlobal_Dpt40_v1&&Dpt>40&&Dpt<60)||(HLT_HIDmesonHITrackingGlobal_Dpt60_v1&&Dpt>60))"
LABELPbPb="PbPb"
PTMINPbPb=20
PTMAXPbPb=100

#

INPUTDATAPPMB="/mnt/hadoop/store/user/hqiu/DfinderProduction_data/skimMerge/Dntuple_PP_MinimumBias_skimEvents.root"
# INPUTMCPPMB_P="/home/qiu110/myHadoop/DfinderProduction_MCMatchOnly/MCWeightMerge/DNtuple_PromptPP_selectMergeWeight.root"
INPUTMCPPMB_P="/mnt/hadoop/store/user/hqiu/DfinderProduction_MCMatchOnly/MCWeightMerge/DNtuple_PromptPP_selectMergeWeight.root"
# INPUTMCPPMB_NP="/home/qiu110/myHadoop/DfinderProduction_MCMatchOnly/MCWeightMerge/DNtuple_PromptPP_selectMergeWeight.root"
INPUTMCPPMB_NP="/mnt/hadoop/store/user/hqiu/DfinderProduction_MCMatchOnly/MCWeightMerge/DNtuple_PromptPP_selectMergeWeight.root"
WEIGHTPPMB_P="DPtSampleWeight*MeasuredDPtWeight"
WEIGHTPPMB_NP="DPtSampleWeight*FonllPtWeight*(0.55369+1.22995*log(DgenBAncestorpt)+-0.49766*log(DgenBAncestorpt)*log(DgenBAncestorpt)+0.05655*log(DgenBAncestorpt)*log(DgenBAncestorpt)*log(DgenBAncestorpt))"
FILTERCUTPPMB="abs(PVz)<15&&pBeamScrapingFilter&&pPAprimaryVertexFilter&&Dy>-1.&&Dy<1.&&Dtrk1highPurity&&Dtrk2highPurity&&Dtrk1Pt>1.0&&Dtrk2Pt>1.0&&Dtrk1PtErr/Dtrk1Pt<0.3&&Dtrk2PtErr/Dtrk2Pt<0.3&&abs(Dtrk1Eta)<1.5&&abs(Dtrk2Eta)<1.5&&(DlxyBS/DlxyBSErr)>2.5"
CENTRALCUTPPMB_ALPHA="Dalpha<0.2"
CENTRALCUTPPMB_DECAYLENGTH="(DsvpvDistance/DsvpvDisErr)>8"
CENTRALCUTPPMB_CHI2="Dchi2cl>0.05" 
CENTRALCUTPPMB_TRKDCA="min(sqrt(pow(Dtrk1Dxy,2)+pow(Dtrk1Dsz,2))/Dtrk1DxyErr,sqrt(pow(Dtrk2Dxy,2)+pow(Dtrk2Dsz,2))/Dtrk2DxyErr)>4.5"
TRGPPMB="(HLT_L1MinimumBiasHF1OR_part1_v1||HLT_L1MinimumBiasHF1OR_part2_v1||HLT_L1MinimumBiasHF1OR_part3_v1||HLT_L1MinimumBiasHF1OR_part4_v1||HLT_L1MinimumBiasHF1OR_part5_v1||HLT_L1MinimumBiasHF1OR_part6_v1||HLT_L1MinimumBiasHF1OR_part7_v1||HLT_L1MinimumBiasHF1OR_part8_v1||HLT_L1MinimumBiasHF1OR_part9_v1||HLT_L1MinimumBiasHF1OR_part10_v1||HLT_L1MinimumBiasHF1OR_part11_v1||HLT_L1MinimumBiasHF1OR_part12_v1||HLT_L1MinimumBiasHF1OR_part13_v1||HLT_L1MinimumBiasHF1OR_part14_v1||HLT_L1MinimumBiasHF1OR_part15_v1||HLT_L1MinimumBiasHF1OR_part16_v1||HLT_L1MinimumBiasHF1OR_part17_v1||HLT_L1MinimumBiasHF1OR_part18_v1||HLT_L1MinimumBiasHF1OR_part19_v1)"
LABELPPMB="PPMB"
PTMINPPMB=10
PTMAXPPMB=20

LABELPPMBLowPt="PPMBLowPt"
PTMINPPMBLowPt=2
PTMAXPPMBLowPt=6

INPUTDATAPbPbMB="/mnt/hadoop/store/user/hqiu/DfinderProduction_data/skimMerge/Dntuple_HIMinimumBias*.root"
INPUTMCPbPbMB_P="/mnt/hadoop/store/user/hqiu/DfinderProduction_MCMatchOnly/MCWeightMerge/DNtuple_PromptPbPb_selectMergeWeight.root"
# INPUTMCPbPbMB_P="/home/qiu110/myHadoop/DfinderProduction_MCMatchOnly/MCWeightMerge/DNtuple_PromptPbPb_selectMergeWeight.root"
INPUTMCPbPbMB_NP="/mnt/hadoop/store/user/hqiu/DfinderProduction_MCMatchOnly/MCWeightMerge/DNtuple_NonPromptPbPb_selectMergeWeight.root"
# INPUTMCPbPbMB_NP="/home/qiu110/myHadoop/DfinderProduction_MCMatchOnly/MCWeightMerge/DNtuple_NonPromptPbPb_selectMergeWeight.root"
WEIGHTPbPbMB_P="DPtSampleWeight*MeasuredDPtWeight*centralityNCollWeight*vertexZWeight"
WEIGHTPbPbMB_NP="DPtSampleWeight*FonllPtWeight*centralityNCollWeight*vertexZWeight*(1.43138+0.37192*log(DgenBAncestorpt)+-0.71770*log(DgenBAncestorpt)*log(DgenBAncestorpt)+0.23269*log(DgenBAncestorpt)*log(DgenBAncestorpt)*log(DgenBAncestorpt)+-0.02181*log(DgenBAncestorpt)*log(DgenBAncestorpt)*log(DgenBAncestorpt)*log(DgenBAncestorpt))"
FILTERCUTPbPbMB="pclusterCompatibilityFilter&&pprimaryVertexFilter&&abs(PVz)<15&&PVzE<4&&phfCoincFilter3&&Dy>-1.&&Dy<1.&&Dtrk1highPurity&&Dtrk2highPurity&&Dtrk1Pt>1.0&&Dtrk2Pt>1.0&&Dtrk1PtErr/Dtrk1Pt<0.3&&Dtrk2PtErr/Dtrk2Pt<0.3&&abs(Dtrk1Eta)<1.5&&abs(Dtrk2Eta)<1.5&&(DlxyBS/DlxyBSErr)>1.5"
CENTRALCUTPbPbMB_ALPHA="Dalpha<0.2"
CENTRALCUTPbPbMB_DECAYLENGTH="(DsvpvDistance/DsvpvDisErr)>6."
CENTRALCUTPbPbMB_CHI2="Dchi2cl>0.05"
CENTRALCUTPbPbMB_TRKDCA="min(sqrt(pow(Dtrk1Dxy,2)+pow(Dtrk1Dsz,2))/Dtrk1DxyErr,sqrt(pow(Dtrk2Dxy,2)+pow(Dtrk2Dsz,2))/Dtrk2DxyErr)>4.5"
TRGPbPbMB="(HLT_HIL1MinimumBiasHF2AND_part1_v1||HLT_HIL1MinimumBiasHF2AND_part2_v1||HLT_HIL1MinimumBiasHF2AND_part3_v1||HLT_HIL1MinimumBiasHF2AND_part4_v1||HLT_HIL1MinimumBiasHF2AND_part5_v1||HLT_HIL1MinimumBiasHF2AND_part6_v1||HLT_HIL1MinimumBiasHF2AND_part7_v1||HLT_HIL1MinimumBiasHF2AND_part8_v1||HLT_HIL1MinimumBiasHF2AND_part9_v1||HLT_HIL1MinimumBiasHF2AND_part10_v1||HLT_HIL1MinimumBiasHF2AND_part11_v1)"
LABELPbPbMB="PbPbMB"
PTMINPbPbMB=10
PTMAXPbPbMB=20

CENTRALCUTPbPbMBLowPt_ALPHA="Dalpha<0.2"
CENTRALCUTPbPbMBLowPt_DECAYLENGTH="(DsvpvDistance/DsvpvDisErr)>8."
CENTRALCUTPbPbMBLowPt_CHI2="Dchi2cl>0.09"
CENTRALCUTPbPbMBLowPt_TRKDCA="min(sqrt(pow(Dtrk1Dxy,2)+pow(Dtrk1Dsz,2))/Dtrk1DxyErr,sqrt(pow(Dtrk2Dxy,2)+pow(Dtrk2Dsz,2))/Dtrk2DxyErr)>6"
LABELPbPbMBLowPt="PbPbMBLowPt"
PTMINPbPbMBLowPt=2
PTMAXPbPbMBLowPt=6

###

if [ $DOSAVEHISTPP_ALPHA -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe 
./projectVariableMD0Dca.exe "$INPUTDATAPP" "$INPUTMCPP_P" "$INPUTMCPP_NP" "$TRGPP" "${FILTERCUTPP}&&${CENTRALCUTPP_DECAYLENGTH}&&${CENTRALCUTPP_CHI2}&&${CENTRALCUTPP_TRKDCA}" "$WEIGHTPP_P" "$WEIGHTPP_NP" "$LABELPP" "$LABEL_ALPHA" "$VAR_ALPHA" "$VARBINS" "$VARMIN_ALPHA" "$VARMAX_ALPHA" "$VARSIGN_ALPHA" "$PTMINPP" "$PTMAXPP"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPbPb_ALPHA -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe 
./projectVariableMD0Dca.exe "$INPUTDATAPbPb" "$INPUTMCPbPb_P" "$INPUTMCPbPb_NP" "$TRGPbPb" "${FILTERCUTPbPb}&&${CENTRALCUTPbPb_DECAYLENGTH}&&${CENTRALCUTPbPb_CHI2}&&${CENTRALCUTPbPb_TRKDCA}" "$WEIGHTPbPb_P" "$WEIGHTPbPb_NP" "$LABELPbPb" "$LABEL_ALPHA" "$VAR_ALPHA" "$VARBINS" "$VARMIN_ALPHA" "$VARMAX_ALPHA" "$VARSIGN_ALPHA" "$PTMINPbPb" "$PTMAXPbPb"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPPMB_ALPHA -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe 
./projectVariableMD0Dca.exe "$INPUTDATAPPMB" "$INPUTMCPPMB_P" "$INPUTMCPPMB_NP" "$TRGPPMB" "${FILTERCUTPPMB}&&${CENTRALCUTPPMB_DECAYLENGTH}&&${CENTRALCUTPPMB_CHI2}&&${CENTRALCUTPPMB_TRKDCA}" "$WEIGHTPPMB_P" "$WEIGHTPPMB_NP" "$LABELPPMB" "$LABEL_ALPHA" "$VAR_ALPHA" "$VARBINS" "$VARMIN_ALPHA" "$VARMAX_ALPHA" "$VARSIGN_ALPHA" "$PTMINPPMB" "$PTMAXPPMB"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPbPbMB_ALPHA -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe
./projectVariableMD0Dca.exe "$INPUTDATAPbPbMB" "$INPUTMCPbPbMB_P" "$INPUTMCPbPbMB_NP" "$TRGPbPbMB" "${FILTERCUTPbPbMB}&&${CENTRALCUTPbPbMB_DECAYLENGTH}&&${CENTRALCUTPbPbMB_CHI2}&&${CENTRALCUTPbPbMB_TRKDCA}" "$WEIGHTPbPbMB_P" "$WEIGHTPbPbMB_NP" "$LABELPbPbMB" "$LABEL_ALPHA" "$VAR_ALPHA" "$VARBINS" "$VARMIN_ALPHA" "$VARMAX_ALPHA" "$VARSIGN_ALPHA" "$PTMINPbPbMB" "$PTMAXPbPbMB"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPPMBLowPt_ALPHA -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe
./projectVariableMD0Dca.exe "$INPUTDATAPPMB" "$INPUTMCPPMB_P" "$INPUTMCPPMB_NP" "$TRGPPMB" "${FILTERCUTPPMB}&&${CENTRALCUTPPMB_DECAYLENGTH}&&${CENTRALCUTPPMB_CHI2}&&${CENTRALCUTPPMB_TRKDCA}" "$WEIGHTPPMB_P" "$WEIGHTPPMB_NP" "$LABELPPMBLowPt" "$LABEL_ALPHA" "$VAR_ALPHA" "$VARBINS" "$VARMIN_ALPHA" "$VARMAX_ALPHA" "$VARSIGN_ALPHA" "$PTMINPPMBLowPt" "$PTMAXPPMBLowPt"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPbPbMBLowPt_ALPHA -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe
./projectVariableMD0Dca.exe "$INPUTDATAPbPbMB" "$INPUTMCPbPbMB_P" "$INPUTMCPbPbMB_NP" "$TRGPbPbMB" "${FILTERCUTPbPbMB}&&${CENTRALCUTPbPbMBLowPt_DECAYLENGTH}&&${CENTRALCUTPbPbMBLowPt_CHI2}&&${CENTRALCUTPbPbMBLowPt_TRKDCA}" "$WEIGHTPbPbMB_P" "$WEIGHTPbPbMB_NP" "$LABELPbPbMBLowPt" "$LABEL_ALPHA" "$VAR_ALPHA" "$VARBINS" "$VARMIN_ALPHA" "$VARMAX_ALPHA" "$VARSIGN_ALPHA" "$PTMINPbPbMBLowPt" "$PTMAXPbPbMBLowPt"
rm projectVariableMD0Dca.exe
fi

#

if [ $DOSAVEHISTPP_DECAYLENGTH -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe 
./projectVariableMD0Dca.exe "$INPUTDATAPP" "$INPUTMCPP_P" "$INPUTMCPP_NP" "$TRGPP" "${FILTERCUTPP}&&${CENTRALCUTPP_ALPHA}&&${CENTRALCUTPP_CHI2}&&${CENTRALCUTPP_TRKDCA}" "$WEIGHTPP_P" "$WEIGHTPP_NP" "$LABELPP" "$LABEL_DECAYLENGTH" "$VAR_DECAYLENGTH" "$VARBINS" "$VARMIN_DECAYLENGTH" "$VARMAX_DECAYLENGTH" "$VARSIGN_DECAYLENGTH" "$PTMINPP" "$PTMAXPP"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPbPb_DECAYLENGTH -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe 
./projectVariableMD0Dca.exe "$INPUTDATAPbPb" "$INPUTMCPbPb_P" "$INPUTMCPbPb_NP" "$TRGPbPb" "${FILTERCUTPbPb}&&${CENTRALCUTPbPb_ALPHA}&&${CENTRALCUTPbPb_CHI2}&&${CENTRALCUTPbPb_TRKDCA}" "$WEIGHTPbPb_P" "$WEIGHTPbPb_NP" "$LABELPbPb" "$LABEL_DECAYLENGTH" "$VAR_DECAYLENGTH" "$VARBINS" "$VARMIN_DECAYLENGTH" "$VARMAX_DECAYLENGTH" "$VARSIGN_DECAYLENGTH" "$PTMINPbPb" "$PTMAXPbPb"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPPMB_DECAYLENGTH -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe 
./projectVariableMD0Dca.exe "$INPUTDATAPPMB" "$INPUTMCPPMB_P" "$INPUTMCPPMB_NP" "$TRGPPMB" "${FILTERCUTPPMB}&&${CENTRALCUTPPMB_ALPHA}&&${CENTRALCUTPPMB_CHI2}&&${CENTRALCUTPPMB_TRKDCA}" "$WEIGHTPPMB_P" "$WEIGHTPPMB_NP" "$LABELPPMB" "$LABEL_DECAYLENGTH" "$VAR_DECAYLENGTH" "$VARBINS" "$VARMIN_DECAYLENGTH" "$VARMAX_DECAYLENGTH" "$VARSIGN_DECAYLENGTH" "$PTMINPPMB" "$PTMAXPPMB"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPbPbMB_DECAYLENGTH -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe 
./projectVariableMD0Dca.exe "$INPUTDATAPbPbMB" "$INPUTMCPbPbMB_P" "$INPUTMCPbPbMB_NP" "$TRGPbPbMB" "${FILTERCUTPbPbMB}&&${CENTRALCUTPbPbMB_ALPHA}&&${CENTRALCUTPbPbMB_CHI2}&&${CENTRALCUTPbPbMB_TRKDCA}" "$WEIGHTPbPbMB_P" "$WEIGHTPbPbMB_NP" "$LABELPbPbMB" "$LABEL_DECAYLENGTH" "$VAR_DECAYLENGTH" "$VARBINS" "$VARMIN_DECAYLENGTH" "$VARMAX_DECAYLENGTH" "$VARSIGN_DECAYLENGTH" "$PTMINPbPbMB" "$PTMAXPbPbMB"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPPMBLowPt_DECAYLENGTH -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe
./projectVariableMD0Dca.exe "$INPUTDATAPPMB" "$INPUTMCPPMB_P" "$INPUTMCPPMB_NP" "$TRGPPMB" "${FILTERCUTPPMB}&&${CENTRALCUTPPMB_ALPHA}&&${CENTRALCUTPPMB_CHI2}&&${CENTRALCUTPPMB_TRKDCA}" "$WEIGHTPPMB_P" "$WEIGHTPPMB_NP" "$LABELPPMBLowPt" "$LABEL_DECAYLENGTH" "$VAR_DECAYLENGTH" "$VARBINS" "$VARMIN_DECAYLENGTH" "$VARMAX_DECAYLENGTH" "$VARSIGN_DECAYLENGTH" "$PTMINPPMBLowPt" "$PTMAXPPMBLowPt"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPbPbMBLowPt_DECAYLENGTH -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe
./projectVariableMD0Dca.exe "$INPUTDATAPbPbMB" "$INPUTMCPbPbMB_P" "$INPUTMCPbPbMB_NP" "$TRGPbPbMB" "${FILTERCUTPbPbMB}&&${CENTRALCUTPbPbMBLowPt_ALPHA}&&${CENTRALCUTPbPbMBLowPt_CHI2}&&${CENTRALCUTPbPbMBLowPt_TRKDCA}" "$WEIGHTPbPbMB_P" "$WEIGHTPbPbMB_NP" "$LABELPbPbMBLowPt" "$LABEL_DECAYLENGTH" "$VAR_DECAYLENGTH" "$VARBINS" "$VARMIN_DECAYLENGTH" "$VARMAX_DECAYLENGTH" "$VARSIGN_DECAYLENGTH" "$PTMINPbPbMBLowPt" "$PTMAXPbPbMBLowPt"
rm projectVariableMD0Dca.exe
fi

#

if [ $DOSAVEHISTPP_CHI2 -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe 
./projectVariableMD0Dca.exe "$INPUTDATAPP" "$INPUTMCPP_P" "$INPUTMCPP_NP" "$TRGPP" "${FILTERCUTPP}&&${CENTRALCUTPP_ALPHA}&&${CENTRALCUTPP_DECAYLENGTH}&&${CENTRALCUTPP_TRKDCA}" "$WEIGHTPP_P" "$WEIGHTPP_NP" "$LABELPP" "$LABEL_CHI2" "$VAR_CHI2" "$VARBINS" "$VARMIN_CHI2" "$VARMAX_CHI2" "$VARSIGN_CHI2" "$PTMINPP" "$PTMAXPP"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPbPb_CHI2 -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe 
./projectVariableMD0Dca.exe "$INPUTDATAPbPb" "$INPUTMCPbPb_P" "$INPUTMCPbPb_NP" "$TRGPbPb" "${FILTERCUTPbPb}&&${CENTRALCUTPbPb_ALPHA}&&${CENTRALCUTPbPb_DECAYLENGTH}&&${CENTRALCUTPbPb_TRKDCA}" "$WEIGHTPbPb_P" "$WEIGHTPbPb_NP" "$LABELPbPb" "$LABEL_CHI2" "$VAR_CHI2" "$VARBINS" "$VARMIN_CHI2" "$VARMAX_CHI2" "$VARSIGN_CHI2" "$PTMINPbPb" "$PTMAXPbPb"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPPMB_CHI2 -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe 
./projectVariableMD0Dca.exe "$INPUTDATAPPMB" "$INPUTMCPPMB_P" "$INPUTMCPPMB_NP" "$TRGPPMB" "${FILTERCUTPPMB}&&${CENTRALCUTPPMB_ALPHA}&&${CENTRALCUTPPMB_DECAYLENGTH}&&${CENTRALCUTPPMB_TRKDCA}" "$WEIGHTPPMB_P" "$WEIGHTPPMB_NP" "$LABELPPMB" "$LABEL_CHI2" "$VAR_CHI2" "$VARBINS" "$VARMIN_CHI2" "$VARMAX_CHI2" "$VARSIGN_CHI2" "$PTMINPPMB" "$PTMAXPPMB"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPbPbMB_CHI2 -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe 
./projectVariableMD0Dca.exe "$INPUTDATAPbPbMB" "$INPUTMCPbPbMB_P" "$INPUTMCPbPbMB_NP" "$TRGPbPbMB" "${FILTERCUTPbPbMB}&&${CENTRALCUTPbPbMB_ALPHA}&&${CENTRALCUTPbPbMB_DECAYLENGTH}&&${CENTRALCUTPbPbMB_TRKDCA}" "$WEIGHTPbPbMB_P" "$WEIGHTPbPbMB_NP" "$LABELPbPbMB" "$LABEL_CHI2" "$VAR_CHI2" "$VARBINS" "$VARMIN_CHI2" "$VARMAX_CHI2" "$VARSIGN_CHI2" "$PTMINPbPbMB" "$PTMAXPbPbMB"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPPMBLowPt_CHI2 -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe
./projectVariableMD0Dca.exe "$INPUTDATAPPMB" "$INPUTMCPPMB_P" "$INPUTMCPPMB_NP" "$TRGPPMB" "${FILTERCUTPPMB}&&${CENTRALCUTPPMB_ALPHA}&&${CENTRALCUTPPMB_DECAYLENGTH}&&${CENTRALCUTPPMB_TRKDCA}" "$WEIGHTPPMB_P" "$WEIGHTPPMB_NP" "$LABELPPMBLowPt" "$LABEL_CHI2" "$VAR_CHI2" "$VARBINS" "$VARMIN_CHI2" "$VARMAX_CHI2" "$VARSIGN_CHI2" "$PTMINPPMBLowPt" "$PTMAXPPMBLowPt"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPbPbMBLowPt_CHI2 -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe
./projectVariableMD0Dca.exe "$INPUTDATAPbPbMB" "$INPUTMCPbPbMB_P" "$INPUTMCPbPbMB_NP" "$TRGPbPbMB" "${FILTERCUTPbPbMB}&&${CENTRALCUTPbPbMBLowPt_ALPHA}&&${CENTRALCUTPbPbMBLowPt_DECAYLENGTH}&&${CENTRALCUTPbPbMBLowPt_TRKDCA}" "$WEIGHTPbPbMB_P" "$WEIGHTPbPbMB_NP" "$LABELPbPbMBLowPt" "$LABEL_CHI2" "$VAR_CHI2" "$VARBINS" "$VARMIN_CHI2" "$VARMAX_CHI2" "$VARSIGN_CHI2" "$PTMINPbPbMBLowPt" "$PTMAXPbPbMBLowPt"
rm projectVariableMD0Dca.exe
fi

#

if [ $DOSAVEHISTPP_TRKDCA -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe
./projectVariableMD0Dca.exe "$INPUTDATAPP" "$INPUTMCPP_P" "$INPUTMCPP_NP" "$TRGPP" "${FILTERCUTPP}&&${CENTRALCUTPP_ALPHA}&&${CENTRALCUTPP_DECAYLENGTH}&&${CENTRALCUTPP_CHI2}" "$WEIGHTPP_P" "$WEIGHTPP_NP" "$LABELPP" "$LABEL_TRKDCA" "$VAR_TRKDCA" "$VARBINS" "$VARMIN_TRKDCA" "$VARMAX_TRKDCA" "$VARSIGN_TRKDCA" "$PTMINPP" "$PTMAXPP"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPbPb_TRKDCA -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe
./projectVariableMD0Dca.exe "$INPUTDATAPbPb" "$INPUTMCPbPb_P" "$INPUTMCPbPb_NP" "$TRGPbPb" "${FILTERCUTPbPb}&&${CENTRALCUTPbPb_ALPHA}&&${CENTRALCUTPbPb_DECAYLENGTH}&&${CENTRALCUTPbPb_CHI2}" "$WEIGHTPbPb_P" "$WEIGHTPbPb_NP" "$LABELPbPb" "$LABEL_TRKDCA" "$VAR_TRKDCA" "$VARBINS" "$VARMIN_TRKDCA" "$VARMAX_TRKDCA" "$VARSIGN_TRKDCA" "$PTMINPbPb" "$PTMAXPbPb"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPPMB_TRKDCA -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe
./projectVariableMD0Dca.exe "$INPUTDATAPPMB" "$INPUTMCPPMB_P" "$INPUTMCPPMB_NP" "$TRGPPMB" "${FILTERCUTPPMB}&&${CENTRALCUTPPMB_ALPHA}&&${CENTRALCUTPPMB_DECAYLENGTH}&&${CENTRALCUTPPMB_CHI2}" "$WEIGHTPPMB_P" "$WEIGHTPPMB_NP" "$LABELPPMB" "$LABEL_TRKDCA" "$VAR_TRKDCA" "$VARBINS" "$VARMIN_TRKDCA" "$VARMAX_TRKDCA" "$VARSIGN_TRKDCA" "$PTMINPPMB" "$PTMAXPPMB"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPbPbMB_TRKDCA -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe
./projectVariableMD0Dca.exe "$INPUTDATAPbPbMB" "$INPUTMCPbPbMB_P" "$INPUTMCPbPbMB_NP" "$TRGPbPbMB" "${FILTERCUTPbPbMB}&&${CENTRALCUTPbPbMB_ALPHA}&&${CENTRALCUTPbPbMB_DECAYLENGTH}&&${CENTRALCUTPbPbMB_CHI2}" "$WEIGHTPbPbMB_P" "$WEIGHTPbPbMB_NP" "$LABELPbPbMB" "$LABEL_TRKDCA" "$VAR_TRKDCA" "$VARBINS" "$VARMIN_TRKDCA" "$VARMAX_TRKDCA" "$VARSIGN_TRKDCA" "$PTMINPbPbMB" "$PTMAXPbPbMB"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPPMBLowPt_TRKDCA -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe
./projectVariableMD0Dca.exe "$INPUTDATAPPMB" "$INPUTMCPPMB_P" "$INPUTMCPPMB_NP" "$TRGPPMB" "${FILTERCUTPPMB}&&${CENTRALCUTPPMB_ALPHA}&&${CENTRALCUTPPMB_DECAYLENGTH}&&${CENTRALCUTPPMB_CHI2}" "$WEIGHTPPMB_P" "$WEIGHTPPMB_NP" "$LABELPPMBLowPt" "$LABEL_TRKDCA" "$VAR_TRKDCA" "$VARBINS" "$VARMIN_TRKDCA" "$VARMAX_TRKDCA" "$VARSIGN_TRKDCA" "$PTMINPPMBLowPt" "$PTMAXPPMBLowPt"
rm projectVariableMD0Dca.exe
fi

if [ $DOSAVEHISTPbPbMBLowPt_TRKDCA -eq 1 ]; then
g++ projectVariableMD0Dca.C $(root-config --cflags --libs) -g -o projectVariableMD0Dca.exe
./projectVariableMD0Dca.exe "$INPUTDATAPbPbMB" "$INPUTMCPbPbMB_P" "$INPUTMCPbPbMB_NP" "$TRGPbPbMB" "${FILTERCUTPbPbMB}&&${CENTRALCUTPbPbMBLowPt_ALPHA}&&${CENTRALCUTPbPbMBLowPt_DECAYLENGTH}&&${CENTRALCUTPbPbMBLowPt_CHI2}" "$WEIGHTPbPbMB_P" "$WEIGHTPbPbMB_NP" "$LABELPbPbMBLowPt" "$LABEL_TRKDCA" "$VAR_TRKDCA" "$VARBINS" "$VARMIN_TRKDCA" "$VARMAX_TRKDCA" "$VARSIGN_TRKDCA" "$PTMINPbPbMBLowPt" "$PTMAXPbPbMBLowPt"
rm projectVariableMD0Dca.exe
fi

##

## 
## if [ $DOFITPP_ALPHA -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPP'", "'$LABEL_ALPHA'", '$VARBINS', '$VARMIN_ALPHA', '$VARMAX_ALPHA', '$VARSIGN_ALPHA')'
## fi
## 
## if [ $DOFITPbPb_ALPHA -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPbPb'", "'$LABEL_ALPHA'", '$VARBINS', '$VARMIN_ALPHA', '$VARMAX_ALPHA', '$VARSIGN_ALPHA')'
## fi
## 
## if [ $DOFITPPMB_ALPHA -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPPMB'", "'$LABEL_ALPHA'", '$VARBINS', '$VARMIN_ALPHA', '$VARMAX_ALPHA', '$VARSIGN_ALPHA')'
## fi
## 
## if [ $DOFITPbPbMB_ALPHA -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPbPbMB'", "'$LABEL_ALPHA'", '$VARBINS', '$VARMIN_ALPHA', '$VARMAX_ALPHA', '$VARSIGN_ALPHA')'
## fi
## 
## if [ $DOFITPPMBLowPt_ALPHA -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPPMBLowPt'", "'$LABEL_ALPHA'", '$VARBINS', '$VARMIN_ALPHA', '$VARMAX_ALPHA', '$VARSIGN_ALPHA')'
## fi
## 
## if [ $DOFITPbPbMBLowPt_ALPHA -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPbPbMBLowPt'", "'$LABEL_ALPHA'", '$VARBINS', '$VARMIN_ALPHA', '$VARMAX_ALPHA', '$VARSIGN_ALPHA')'
## fi
## 
## ###
## 
## if [ $DOFITPP_DECAYLENGTH -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPP'", "'$LABEL_DECAYLENGTH'", '$VARBINS', '$VARMIN_DECAYLENGTH', '$VARMAX_DECAYLENGTH', '$VARSIGN_DECAYLENGTH')'
## fi
## 
## if [ $DOFITPbPb_DECAYLENGTH -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPbPb'", "'$LABEL_DECAYLENGTH'", '$VARBINS', '$VARMIN_DECAYLENGTH', '$VARMAX_DECAYLENGTH', '$VARSIGN_DECAYLENGTH')'
## fi
## 
## if [ $DOFITPPMB_DECAYLENGTH -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPPMB'", "'$LABEL_DECAYLENGTH'", '$VARBINS', '$VARMIN_DECAYLENGTH', '$VARMAX_DECAYLENGTH', '$VARSIGN_DECAYLENGTH')'
## fi
## 
## if [ $DOFITPbPbMB_DECAYLENGTH -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPbPbMB'", "'$LABEL_DECAYLENGTH'", '$VARBINS', '$VARMIN_DECAYLENGTH', '$VARMAX_DECAYLENGTH', '$VARSIGN_DECAYLENGTH')'
## fi
## 
## if [ $DOFITPPMBLowPt_DECAYLENGTH -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPPMBLowPt'", "'$LABEL_DECAYLENGTH'", '$VARBINS', '$VARMIN_DECAYLENGTH', '$VARMAX_DECAYLENGTH', '$VARSIGN_DECAYLENGTH')'
## fi
## 
## if [ $DOFITPbPbMBLowPt_DECAYLENGTH -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPbPbMBLowPt'", "'$LABEL_DECAYLENGTH'", '$VARBINS', '$VARMIN_DECAYLENGTH', '$VARMAX_DECAYLENGTH', '$VARSIGN_DECAYLENGTH')'
## fi
## 
## ###
## 
## if [ $DOFITPP_CHI2 -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPP'", "'$LABEL_CHI2'", '$VARBINS', '$VARMIN_CHI2', '$VARMAX_CHI2', '$VARSIGN_CHI2')'
## fi
## 
## if [ $DOFITPbPb_CHI2 -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPbPb'", "'$LABEL_CHI2'", '$VARBINS', '$VARMIN_CHI2', '$VARMAX_CHI2', '$VARSIGN_CHI2')'
## fi
## 
if [ $DOFITPPMB_CHI2 -eq 1 ]; then
root -l -b -q 'calculateYieldRatio.C("'$LABELPPMB'", "'$LABEL_CHI2'", '$VARBINS', '$VARMIN_CHI2', '$VARMAX_CHI2', '$VARSIGN_CHI2')'
fi
## 
## if [ $DOFITPbPbMB_CHI2 -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPbPbMB'", "'$LABEL_CHI2'", '$VARBINS', '$VARMIN_CHI2', '$VARMAX_CHI2', '$VARSIGN_CHI2')'
## fi
## 
## if [ $DOFITPPMBLowPt_CHI2 -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPPMBLowPt'", "'$LABEL_CHI2'", '$VARBINS', '$VARMIN_CHI2', '$VARMAX_CHI2', '$VARSIGN_CHI2')'
## fi
## 
## if [ $DOFITPbPbMBLowPt_CHI2 -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPbPbMBLowPt'", "'$LABEL_CHI2'", '$VARBINS', '$VARMIN_CHI2', '$VARMAX_CHI2', '$VARSIGN_CHI2')'
## fi
## 
## ###
## 
## if [ $DOFITPP_TRKDCA -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPP'", "'$LABEL_TRKDCA'", '$VARBINS', '$VARMIN_TRKDCA', '$VARMAX_TRKDCA', '$VARSIGN_TRKDCA')'
## fi
## 
## if [ $DOFITPbPb_TRKDCA -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPbPb'", "'$LABEL_TRKDCA'", '$VARBINS', '$VARMIN_TRKDCA', '$VARMAX_TRKDCA', '$VARSIGN_TRKDCA')'
## fi
## 
## if [ $DOFITPPMB_TRKDCA -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPPMB'", "'$LABEL_TRKDCA'", '$VARBINS', '$VARMIN_TRKDCA', '$VARMAX_TRKDCA', '$VARSIGN_TRKDCA')'
## fi
## 
## if [ $DOFITPbPbMB_TRKDCA -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPbPbMB'", "'$LABEL_TRKDCA'", '$VARBINS', '$VARMIN_TRKDCA', '$VARMAX_TRKDCA', '$VARSIGN_TRKDCA')'
## fi
## 
## if [ $DOFITPPMBLowPt_TRKDCA -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPPMBLowPt'", "'$LABEL_TRKDCA'", '$VARBINS', '$VARMIN_TRKDCA', '$VARMAX_TRKDCA', '$VARSIGN_TRKDCA')'
## fi
## 
## if [ $DOFITPbPbMBLowPt_TRKDCA -eq 1 ]; then
## root -l -b -q 'calculateYieldRatio.C("'$LABELPbPbMBLowPt'", "'$LABEL_TRKDCA'", '$VARBINS', '$VARMIN_TRKDCA', '$VARMAX_TRKDCA', '$VARSIGN_TRKDCA')'
## fi
## 
## ##
## 
## if [ $DORATIO_ALPHA -eq 1 ]; then
## g++ plotRatioDoubleratio.C $(root-config --cflags --libs) -g -o plotRatioDoubleratio.exe 
## ./plotRatioDoubleratio.exe "$LABEL_ALPHA" "VAR_ALPHA"
## rm plotRatioDoubleratio.exe
## fi
## 
## if [ $DORATIO_DECAYLENGTH -eq 1 ]; then
## g++ plotRatioDoubleratio.C $(root-config --cflags --libs) -g -o plotRatioDoubleratio.exe 
## ./plotRatioDoubleratio.exe "$LABEL_DECAYLENGTH" "VAR_DECAYLENGTH"
## rm plotRatioDoubleratio.exe
## fi
## 
## if [ $DORATIO_CHI2 -eq 1 ]; then
## g++ plotRatioDoubleratio.C $(root-config --cflags --libs) -g -o plotRatioDoubleratio.exe 
## ./plotRatioDoubleratio.exe "$LABEL_CHI2" "VAR__CHI2"
## rm plotRatioDoubleratio.exe
## fi
