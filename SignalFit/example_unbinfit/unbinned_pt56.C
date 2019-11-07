#include <TCanvas.h>

#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLorentzVector.h>

#include <TMath.h>
#include <TNamed.h>
#include <TNtuple.h>

#include <TStyle.h>
#include <TTree.h>

#include <RooAddPdf.h>
#include <RooDataSet.h>
#include <RooFormulaVar.h>
#include <RooGaussian.h>
#include <RooPlot.h>

#include <RooChebychev.h>

#include <RooRealVar.h>

using namespace RooFit;


void unbinned_pt56(){

TFile *file      = new TFile("mass_tuple_whole.root");
TFile *file_pt56 = new TFile("5_6_pp_officialTMVAcuts_2gaus_pol3_weighted_ptcuts_eventselection_changefit_L_changefit_function.root");

TH1F *h1     = (TH1F*) file_pt56->Get("h1");
TF1 *f_pt56  = (TF1*) file_pt56->Get("f3");

// from 5-6 TMVA file
double w1 = f_pt56->GetParameter(2);
double w2 = f_pt56->GetParameter(4);
double r1 = f_pt56->GetParameter(3);
double mean_value = f_pt56->GetParameter(1);
double float_width_parameter = f_pt56->GetParameter(5);

 
double numbers = f_pt56->GetParameter(0);// signal
double numberb = f_pt56->Integral(2.2,2.4)/h1->GetBinWidth(0)-numbers;//bkgd
TH1F *h_output = new TH1F ("h_output","output",2,0,2);

 cout<<"Initial parameter values, from the TMVAcuts*.root file: \n";
 cout<<"w1 = f_pt56->GetParameter(2): "<<w1<<endl;
 cout <<"w2 = f_pt56->GetParameter(4): "<< w2<<endl;
 cout<<"r1 = f_pt56->GetParameter(3): "<<r1<<endl;
 cout<<"mean_value = f_pt56->GetParameter(1): "<<mean_value<<endl;
 cout<<"float_width_parameter = f_pt56->GetParameter(5): "<<float_width_parameter<<endl;
 
cout<<"numbers= f_pt56->GetParameter(0): "<<numbers<<endl;
cout<<"numberb = f_pt56->Integral(2.2,2.4)/h1->GetBinWidth(0)-numbers:  "<<numberb<<endl;


TNtuple *nt4 = (TNtuple*)file->Get("nt4");
RooRealVar mass("mass","mass",2.2,2.4);

RooDataSet ds("ds","ds",nt4,mass);
ds.Print();

RooPlot *massframe = new RooPlot("mass","mass",mass,2.2,2.4,25);
ds.plotOn(massframe,MarkerStyle(20),MarkerColor(2),MarkerSize(1.5),LineColor(2));
TCanvas *c1 = new TCanvas("c1","c1",600,600);
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);
gStyle->SetTitleX(0.8f);
gStyle->SetTitleY(0.95);
gStyle->SetTitleW(0.2f);
gStyle->SetTitleFontSize(0.06);

massframe->Draw();
cout<<"begin fitting"<<endl;

RooRealVar mean("mean","mean",mean_value,2.2,2.4);
RooRealVar width1("width1","width1",w1,w1,w1);
RooRealVar width2("width2","width2",w2,w2,w2);
RooRealVar ratio("ratio","ratio",r1,r1,r1);
RooRealVar p1("p1","p1",0.,-1.0,1.0);
RooRealVar p2("p2","p2",0.,-1.0,1.0);
RooRealVar p3("p3","p3",0.,-1.0,1.0);

 
RooRealVar float_width("float_width","float_width",float_width_parameter,-0.9,0.9);
 
RooFormulaVar scale_width1("scale width1","scaled width1","width1*(1+float_width)",RooArgSet(width1,float_width));
RooFormulaVar scale_width2("scale width2","scaled width2","width2*(1+float_width)",RooArgSet(width2,float_width));
RooGaussian sigma1("sigma1","gauss(mass,mean,scale_width1)",mass,mean,scale_width1);
RooGaussian sigma2("sigma2","gauss(mass,mean,scale_width2)",mass,mean,scale_width2);
RooAddPdf signal("signal","signal",RooArgList(sigma1,sigma2),ratio);

 
RooChebychev back("back","second polynomial",mass,RooArgList(p1,p2,p3));
RooRealVar NumSig("NumSig","Number of Signal",numbers,-1e4,1e4);
RooRealVar NumBkg("NumBkg","Number of Background",numberb,-1e6,1e6);

RooAddPdf model("model","model",RooArgList(signal,back),RooArgList(NumSig,NumBkg));

model.fitTo(ds,Extended(kTRUE),Range(2.2,2.4)); 


model.plotOn(massframe,Range("full"),LineColor(2));
model.plotOn(massframe,Components(back),LineColor(3));

massframe->SetXTitle("m_{pk#pi} (GeV/c^{2})");
massframe->SetTitleOffset(1.5,"Y");
massframe->SetTitleOffset(1.0,"X");
massframe->SetTitleFont(42,"X");
massframe->SetTitleFont(42,"Y");
massframe->SetLabelFont(42,"Y");
massframe->SetLabelFont(42,"X");
massframe->SetLabelSize(0.033,"X");
massframe->SetLabelSize(0.032,"Y");
massframe->SetTitleSize(0.034,"X");
massframe->SetTitleSize(0.033,"Y");
massframe->SetMinimum(18000);
//massframe->SetMinimum(18000);

massframe->Draw();


TLatex* tex;
tex = new TLatex(0.52,0.85,"pp: 5GeV < P_{T} < 6GeV");
tex->SetNDC();
tex->SetTextFont(42);
tex->SetTextSize(0.033);
tex->SetLineWidth(2);
tex->Draw();

int a = round(NumSig.getVal());
int b = round(NumSig.getError());
tex = new TLatex(0.52,0.8,Form("yield: %i #pm %i ",a,b));
tex->SetNDC();
tex->SetTextFont(42);
tex->SetTextSize(0.033);
tex->SetLineWidth(2);
tex->Draw();

tex = new TLatex(0.52,0.75,"|y| < 1.0");
tex->SetNDC();
tex->SetTextFont(42);
tex->SetTextSize(0.033);
tex->SetLineWidth(2);
tex->Draw();

TLatex Tl; 
Tl.SetNDC();
Tl.SetTextAlign(12);
Tl.SetTextSize(0.045);
Tl.SetTextFont(42);
Tl.DrawLatex(0.1,0.93, "#font[61]{CMS }");
Tl.DrawLatex(0.22,0.93, "#scale[0.8]{Preliminary}");
Tl.DrawLatex(0.59,0.93, "#scale[0.8]{38.2 nb^{-1} (5.02 TeV)}");//pp
//Tl.DrawLatex(0.57,0.93, "#scale[0.8]{101.573 #mub^{-1} (5.02 TeV)}");//PbPb cen30_100
h_output->SetBinContent(1,a);
h_output->SetBinError(1,b);
c1->Range(2.2,9000,2.4,12500);

c1->SaveAs("test_pp_pt56.gif");

}
