#ifndef PLOTFUN_H
#define PLOTFUN_H

#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TString.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TColor.h"
#include "TLine.h"
#include "TTree.h"
#include "TTimeStamp.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "config.h"

using namespace std;


const int kRainBow = 55; //just in case

int darkred = TColor::GetColorDark(2);
int darkgreen = TColor::GetColorDark(3);
int darkblue = TColor::GetColorDark(4);
int lightblue = TColor::GetColorBright(4);
int darkviolet = TColor::GetColorDark(kMagenta);

bool buildFromVector = false;
vector<float> buildvector;
int buildnbins = 20;
float buildxmin = 0;
float buildxmax = 1;
int buildnbinsy = 20;
float buildymin = 0;
float buildymax = 1;
int buildndiv = -1;
vector<TH1 *> allhists;







