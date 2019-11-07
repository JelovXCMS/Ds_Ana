#include <TMath.h>
#include <TF1.h>
#include <TNamed.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLorentzVector.h>
#include <RooDataSet.h>
#include <RooFormulaVar.h>
#include <RooRealVar.h>
using namespace RooFit;
void unbinned_fit_test_modified_PbPb_trytomanytimes_PbPb_realPbPb_test_massrange_willdelet_30_100_Chebychev_pol2(){
//pp 10-20
//TFile * f = new TFile("/scratch/halstead/x/xiao147/pp2015_lambdaC_pt8_makeup_halstead/MB1_20_pp_2015_adddaughterpt/pp_eventselection_trackquality_trigger/TMVA_histogram_version2/mass_tuple_whole.root");
//pp 10-20 pol2
//TFile *f_fun = new TFile("/home/xiao147/private/newchannel_lambda_CtoproduceDntuple/CMSSW_7_5_8_patch3/src/Bfinder/Bfinder/Dntuple/TMVA_official_MC/signal_exaction/ROOT_signalextraction_witheventselection/signalectraction/forAN_plot/change_fit_function/pp/pol2/10_20_pp_officialTMVAcuts_2gaus_pol2_weighted_ptcuts_eventselection_withLoption_change_fit_funcition.root");
//TFile *f_fun =  new TFile("/home/xiao147/private/newchannel_lambda_CtoproduceDntuple/CMSSW_7_5_8_patch3/src/Bfinder/Bfinder/Dntuple/TMVA_official_MC/signal_exaction/ROOT_signalextraction_witheventselection/signalectraction/forAN_plot/change_fit_function/pp/pol3/10_20_pp_officialTMVAcuts_2gaus_pol3_weighted_ptcuts_eventselection_changefit_L_changefit_function.root");
//8-10
//TFile *f_fun = new TFile("/home/xiao147/private/newchannel_lambda_CtoproduceDntuple/CMSSW_7_5_8_patch3/src/Bfinder/Bfinder/Dntuple/TMVA_official_MC/signal_exaction/ROOT_signalextraction_witheventselection/signalectraction/forAN_plot/change_fit_function/pp/pol3/8_10_pp_officialTMVAcuts_2gaus_pol3_weighted_ptcuts_eventselection_changefit_L_changefit_function.root");

//6-8
//TFile *f_fun = new TFile("/home/xiao147/private/newchannel_lambda_CtoproduceDntuple/CMSSW_7_5_8_patch3/src/Bfinder/Bfinder/Dntuple/TMVA_official_MC/signal_exaction/ROOT_signalextraction_witheventselection/signalectraction/forAN_plot/change_fit_function/pp/pol3/6_8_pp_officialTMVAcuts_2gaus_pol3_weighted_ptcuts_eventselection_changefit_L_changefit_function.root");
//5-6
//TFile *f_fun = new TFile("/home/xiao147/private/newchannel_lambda_CtoproduceDntuple/CMSSW_7_5_8_patch3/src/Bfinder/Bfinder/Dntuple/TMVA_official_MC/signal_exaction/ROOT_signalextraction_witheventselection/signalectraction/forAN_plot/change_fit_function/pp/pol3/5_6_pp_officialTMVAcuts_2gaus_pol3_weighted_ptcuts_eventselection_changefit_L_changefit_function.root");

////PbPb cen30-100%

TFile * f = new TFile("/scratch/halstead/x/xiao147/whole_2015_PbPb_goldernjson_trackingonly_lambdaC_pkpi/whole_addbranch/firetrigger/newest_addeventandtrigger_trackquality_branch/After_apply_ARC_comments_massplots_0718/TMVA_cen0_100_traincuts_alltrackptratio_0828/Kevin_cuts_fromthisTMVA_alltrackptratio/cen30_100/cen30_100_Whole_Whole.root");
//TFile *f_fun = new TFile("/home/xiao147/private/newchannel_lambda_CtoproduceDntuple/CMSSW_7_5_8_patch3/src/Bfinder/Bfinder/Dntuple/TMVA_official_MC/signal_exaction/ROOT_signalextraction_witheventselection/signalectraction/forAN_plot/change_fit_function/PbPb30_100/cen30_100_10_20_PbPb_official_TMVA_2gaus_pol3_1mean_PbPbMC_2pt_fitoption_IL_0310_change_fitfunction.root");
TFile *f_fun = new TFile("/home/xiao147/private/newchannel_lambda_CtoproduceDntuple/CMSSW_7_5_8_patch3/src/Bfinder/Bfinder/Dntuple/TMVA_official_MC/ARC_comments/signal/ROOT_files_AfterARC/should_be_final_apply_ARC_comments/2gaus_1mean/PbPb_uniformcuts_0828_alltrackptraito_TMVA/Kevin_PbPb_cuts_basedon0828_TMVA_cuts/PbPb_10_20_2gaus_1mean_cen30_100_uniformcuts_TMVA_basedon_alltrackptratio_0828_sameaspp_daughtertrackptratio.root");

/*
TFile * f = new TFile("/scratch/halstead/x/xiao147/2015triggerdata_567_resubmit_pt8/whole_addbranch/MB5_7_addtrigger_eventselection_trackquality_branch/TMVA_histogram/8_10_cen30_100_whole.root");
TFile *f_fun =  new TFile("/home/xiao147/private/newchannel_lambda_CtoproduceDntuple/CMSSW_7_5_8_patch3/src/Bfinder/Bfinder/Dntuple/TMVA_official_MC/signal_exaction/ROOT_signalextraction_witheventselection/signalectraction/forAN_plot/change_fit_function/PbPb30_100/cen30_100_8_10_PbPb_official_TMVA_2gaus_pol3_1mean_PbPbMC_2pt_fitoption_IL_0310_change_fitfunction.root");
*/
/////PbPb cen0_100%
/*
TFile *f = new TFile("/scratch/halstead/x/xiao147/whole_2015_PbPb_goldernjson_trackingonly_lambdaC_pkpi/whole_addbranch/firetrigger/newest_addeventandtrigger_trackquality_branch/cen0_100_histograms/cen0_100_TMVA_whole.root");
TFile *f_fun = new TFile("/home/xiao147/private/newchannel_lambda_CtoproduceDntuple/CMSSW_7_5_8_patch3/src/Bfinder/Bfinder/Dntuple/TMVA_official_MC/signal_exaction/ROOT_signalextraction_witheventselection/signalectraction/forAN_plot/change_fit_function/PbPb_0_100/cen0_100_10_20_PbPb_official_TMVA_2gaus_pol3_1mean_PbPbMC_2pt_fitoption_IL_0310_change_fitfunction.root");
*/
///PbPb cen0_30%
/*
TFile *f = new TFile("/scratch/halstead/x/xiao147/whole_2015_PbPb_goldernjson_trackingonly_lambdaC_pkpi/whole_addbranch/firetrigger/newest_addeventandtrigger_trackquality_branch/cen0_100_histograms/cen0_30_TMVA_whole.root");
TFile *f_fun =  new TFile("/home/xiao147/private/newchannel_lambda_CtoproduceDntuple/CMSSW_7_5_8_patch3/src/Bfinder/Bfinder/Dntuple/TMVA_official_MC/signal_exaction/ROOT_signalextraction_witheventselection/signalectraction/forAN_plot/change_fit_function/PbPb_0_30/cen0_30_10_20_PbPb_official_TMVA_2gaus_pol3_1mean_PbPbMC_2pt_fitoption_IL_0310_change_fitfunction.root");
*/

TF1 *f3 = (TF1*) f_fun->Get("f3");
//TH1F *h1 = (TH1F*) f_fun->Get("h1");
double w1 = f3->GetParameter(2);
double w2 = f3->GetParameter(4);
double r1 = f3->GetParameter(3);
double mean_value = f3->GetParameter(1);
double pol1_value = f3->GetParameter(7);
double pol2_value = f3->GetParameter(8);
double pol3_value = f3->GetParameter(9);
double float_width_parameter = f3->GetParameter(5);
//double numbers = f3->GetParameter(0);
//double numberb = f3->Integral(2.2,2.4)/h1->GetBinWidth(0)-f3->GetParameter(0);
//double numbers = 3800;
//double numberb = 5.780790e+05;
TH1F *h_output = new TH1F ("h_output","output",2,0,2);
//cout<<"numbers"<<numbers<<endl;
//cout<<"numberb"<<numberb<<endl;
/*
ifstream file_stream("/home/xiao147/private/lambdaC_mix_toDntuple/CMSSW_7_5_8_patch3/src/Bfinder/Bfinder/Dntuple/5_6_7_files_tracker/MB7_resubmit.txt");

char* filename = new char[10000];
int nfile_tot=20;
int ifile = 0;
TChain *ch1 = new TChain("nt");
for(int ii = 0; ii<nfile_tot; ii++) {

	file_stream >> filename;

	ifile++;

	cout << "file  = " << ifile<<" name = "<<filename <<endl;

	TFile fin(filename);
	ch1->Add(filename);
	fin.Close();
}
*/



TNtuple *nt = (TNtuple*)f->Get("nt1");
TH1F *hb = new TH1F("hb","hb",42,2.108,2.444);
nt->Draw("mass>>hb");
double numbers = 2100;
double numberb = hb->Integral();
RooRealVar mass("mass","mass",2.1,2.444);
//RooRealVar mass("mass","mass",2.16,2.42);
RooDataSet ds("ds","ds",nt,mass);
ds.Print();
//RooPlot* massframe = mass.frame(Title("mass"));
/*
RooPlot *massframe = new RooPlot("mass","mass",mass,2.2,2.4,50);
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
*/
cout<<"begin fitting"<<endl;
RooRealVar mean("mean","mean",mean_value,2.26,2.31);
RooRealVar width1("width1","width1",w1,w1,w1);
RooRealVar width2("width2","width2",w2,w2,w2);
RooRealVar ratio("ratio","ratio",r1,r1,r1);
RooRealVar p1("p1","p1",0,-1,1);
RooRealVar p2("p2","p2",0,-1,1);
RooRealVar p3("p3","p3",0,-1,1);
RooRealVar float_width("float_width","float_width",float_width_parameter,-0.5,0.5);
RooFormulaVar scale_width1("scale width1","scaled width1","width1*(1+float_width)",RooArgSet(width1,float_width));
RooFormulaVar scale_width2("scale width2","scaled width2","width2*(1+float_width)",RooArgSet(width2,float_width));
RooGaussian sigma1("sigma1","gauss(mass,mean,scale_width1)",mass,mean,scale_width1);
RooGaussian sigma2("sigma2","gauss(mass,mean,scale_width2)",mass,mean,scale_width2);
RooAddPdf signal("signal","signal",RooArgList(sigma1,sigma2),ratio);
RooChebychev back("back","background",mass,RooArgList(p1,p2,p3));
RooRealVar NumSig("NumSig","Number of Signal",numbers,-1e4,1e4);
RooRealVar NumBkg("NumBkg","Number of Background",numberb,-1e6,1e6);
RooAddPdf model("model","model",RooArgList(signal,back),RooArgList(NumSig,NumBkg));

model.fitTo(ds,Extended(kTRUE),Range(2.108,2.444),Minos(kTRUE)); 

RooPlot *massframe = new RooPlot("mass","mass",mass,2.108,2.444,42);
ds.plotOn(massframe,Name("ds"),MarkerStyle(20),MarkerColor(1),MarkerSize(1),LineColor(1));
TCanvas *c1 = new TCanvas("c1","c1",600,600);
c1->SetLeftMargin(0.12);
c1->SetRightMargin(0.05);
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);
gStyle->SetTitleX(0.8f);
gStyle->SetTitleY(0.95);
gStyle->SetTitleW(0.2f);
gStyle->SetTitleFontSize(0.06);
//try the TStyle setting for canvas
gStyle->SetCanvasBorderMode(0);
gStyle->SetCanvasColor(kWhite);
gStyle->SetCanvasDefX(0);//position on screen
gStyle->SetCanvasDefY(0);
///TStyle setting for the Margins:
gStyle->SetPadTopMargin(0.05);
gStyle->SetPadBottomMargin(0.13);
gStyle->SetPadLeftMargin(0.16);
gStyle->SetPadRightMargin(0.02);
//TStyle setting about the Postscript options
gStyle->SetPaperSize(20.,20.);
//TStyle setting about the frame
gStyle->SetFrameBorderMode(0);
gStyle->SetFrameBorderSize(1);
gStyle->SetFrameFillColor(0);
gStyle->SetFrameFillStyle(0);
gStyle->SetFrameLineColor(1);
gStyle->SetFrameLineStyle(1);
gStyle->SetFrameLineWidth(1);


massframe->Draw();

model.plotOn(massframe,Name("model"),Range(2.108,2.444),LineColor(2));
model.plotOn(massframe,Name("back"),Components(back),LineColor(9));


massframe->GetXaxis()->CenterTitle();
massframe->GetYaxis()->CenterTitle();
massframe->SetXTitle("m_{pk#pi} (GeV/c^{2})");
massframe->SetTitleOffset(1.6,"Y");
massframe->SetTitleOffset(0.9,"X");
massframe->SetTitleFont(42,"X");
massframe->SetTitleFont(42,"Y");
massframe->SetLabelFont(42,"Y");
massframe->SetLabelFont(42,"X");
massframe->SetLabelSize(0.033,"X");
massframe->SetLabelSize(0.033,"Y");
massframe->SetTitleSize(0.045,"X");
massframe->SetTitleSize(0.04,"Y");
massframe->SetAxisRange(4200,7300,"Y");
massframe->GetXaxis()->SetRangeUser(2.108,2.444);
massframe->Draw();


TLatex* tex;
tex = new TLatex(0.16,0.85,"10 < p_{T} < 20 GeV/c");
tex->SetNDC();
tex->SetTextFont(42);
tex->SetTextSize(0.04);
tex->SetLineWidth(2);
tex->Draw();

tex =  new TLatex(0.75,0.85,"#font[61]{#Lambda_{C}^{+} + #Lambda_{C}^{-}}");
tex->SetNDC();
tex->SetTextFont(42);
tex->SetTextSize(0.055);
tex->SetLineWidth(2);
tex->Draw();

int a = round(NumSig.getVal());
int b = round(NumSig.getAsymErrorLo());
int c = round(NumSig.getAsymErrorHi());
int d = (c-b)/2;
cout<<"low: "<<b<<endl;
cout<<"high: "<<c<<endl;
cout<<"err: "<<d<<endl;
tex = new TLatex(0.16,0.7,Form("Raw yield: %i^{+%i}_{%i} ",a,c,b));
tex->SetNDC();
tex->SetTextFont(42);
tex->SetTextSize(0.04);
tex->SetLineWidth(2);
tex->Draw();

tex = new TLatex(0.16,0.75,"Cent. 30-100%");
tex->SetNDC();
tex->SetTextFont(42);
tex->SetTextSize(0.04);
tex->SetLineWidth(2);
tex->Draw();


tex = new TLatex(0.16,0.8,"|y| < 1.0");
tex->SetNDC();
tex->SetTextFont(42);
tex->SetTextSize(0.04);
tex->SetLineWidth(2);
tex->Draw();
/*
TLatex Tl; 
Tl.SetNDC();
Tl.SetTextAlign(12);
Tl.SetTextSize(0.045);
Tl.SetTextFont(42);
Tl.DrawLatex(0.1,0.93, "#font[61]{CMS }");
Tl.DrawLatex(0.22,0.93, "#scale[0.8]{Preliminary}");
//Tl.DrawLatex(0.59,0.93, "#scale[0.8]{0.0382 pb^{-1} (5.02 TeV)}");//pp
//Tl.DrawLatex(0.59,0.93, "#scale[0.8]{43.9 #mub^{-1} (5.02 TeV)}");//PbPb cen0_30
Tl.DrawLatex(0.57,0.93, "#scale[0.8]{101.573 #mub^{-1} (5.02 TeV)}");//PbPb cen30_100
*/
TLatex Tl;
Tl.SetNDC();
Tl.SetTextAlign(12);
Tl.SetTextSize(0.05);
Tl.SetTextFont(42);
Tl.DrawLatex(0.12,0.93, "#font[61]{CMS }");
//Tl.DrawLatex(0.57,0.93, "#scale[0.8]{101.6 #mub^{-1} (5.02 TeV PbPb)}");//PbPb cen30_100
Tl.DrawLatex(0.56,0.93, "#scale[0.8]{102 #mub^{-1} (5.02 TeV PbPb)}");//PbPb cen30_100

TLatex Tl2;
Tl.SetNDC();
Tl.SetTextAlign(12);
Tl.SetTextSize(0.05*0.75);
Tl.SetTextFont(42);
Tl.DrawLatex(0.23,0.93, "#font[52]{Preliminary}");

leg = new TLegend(0.16,0.16,0.3,0.3);
leg->SetTextSize(0.04);
leg->AddEntry("ds","Data","LP");
leg->AddEntry("model","Signal+Background","L");
leg->AddEntry("back","Background","L");
leg->SetBorderSize(0);
leg->Draw();






h_output->SetBinContent(1,a);
h_output->SetBinError(1,c);
h_output->SetBinContent(2,a);
h_output->SetBinError(2,b);





//c1->SaveAs("/home/xiao147/private/newchannel_lambda_CtoproduceDntuple/CMSSW_7_5_8_patch3/src/Bfinder/Bfinder/Dntuple/TMVA_official_MC/ARC_comments/signal/ROOT_files_AfterARC/should_be_final_apply_ARC_comments/PbPb_uniformcuts_invariantmass_fit_basedonMC_expandmassrange/2gaus_1mean_Chebychevpol3_expandmassrange_avoidedge/PbPb_10_20_MC_2gaus_1mean_Chebychevpol3_Minos_daughterptraio_unioformcuts_cen30_100_expandmassrange.gif");
c1->SaveAs("/home/xiao147/private/newchannel_lambda_CtoproduceDntuple/CMSSW_7_5_8_patch3/src/Bfinder/Bfinder/Dntuple/TMVA_official_MC/ARC_comments/signal/ROOT_files_AfterARC/should_be_final_apply_ARC_comments/0828_PbPb_uniformcuts_Kevincuts_basedon0828_TMVA_expandmassrange_withalltrackpt_ratio/Chebychevpol3/PbPb_10_20_MC_2gaus_1mean_Chebychevpol3_Minos_daughterptraio_unioformcuts_cen30_100_expandmassrange_Kevincut_basedon0828TMVA.pdf");

TFile *result = new TFile("/home/xiao147/private/newchannel_lambda_CtoproduceDntuple/CMSSW_7_5_8_patch3/src/Bfinder/Bfinder/Dntuple/TMVA_official_MC/ARC_comments/signal/ROOT_files_AfterARC/should_be_final_apply_ARC_comments/0828_PbPb_uniformcuts_Kevincuts_basedon0828_TMVA_expandmassrange_withalltrackpt_ratio/Chebychevpol3/PbPb_10_20_MC_2gaus_1mean_Chebychevpol3_Minos_daughterptratio_uniformcuts_cen30_100_expandmassrange_basedon0828TMVA.root","RECREATE");
h_output->Write();
result->Close();

}
