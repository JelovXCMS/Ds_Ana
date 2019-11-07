#!/bin/sh


#int th2_correlation(TString Svar="Dtrk1PtErr/Dtrk1Pt", int nbin=100, double binlow =0, double binhi=0.1

root -b -q 'th2_correlation.C("Dtrk1PtErr/Dtrk1Pt",100, 0, 0.1)' >> out.txt &
root -b -q 'th2_correlation.C("Pvnchi2",100, 0, 2)' >> out.txt &
