#ifndef _CMS_LUMI_H_
#define _CMS_LUMI_H_


#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TASImage.h"

//
// Global variables
//

TString cmsText     = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold

// bool writeExtraText = false;
Int_t writeExtraText = 1;

TString extraText   = "Preliminary";

float extraTextFont = 52;  // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
float lumiTextSize     = 0.6;
float lumiTextOffset   = 0.2;
float cmsTextSize      = 0.75;
float cmsTextOffset    = 0.1;  // only used in outOfFrame version

float relPosX    = 0.045;
float relPosY    = 0.035;
float relExtraDY = 1.2;

// ratio of "CMS" and extra text size
float extraOverCmsTextSize  = 0.76;

TString lumi_13TeV = "20.1 fb^{-1}";
TString lumi_8TeV  = "19.7 fb^{-1}";
TString lumi_7TeV  = "5.1 fb^{-1}";
TString lumi_sqrtS = "";

TString lumi_5TeVpp2018  = "650 pb^{-1}";
TString lumi_5TeVPbPb2018  = "0.2 nb^{-1}";

TString lumi_5TeVpp  = "27.4 pb^{-1}";
TString lumi_5TeVPbPb  = "530 #mub^{-1}";

bool drawLogo      = false;

void CMS_lumi( TPad* pad, int iPeriod=3, int iPosX=10 );


#endif
