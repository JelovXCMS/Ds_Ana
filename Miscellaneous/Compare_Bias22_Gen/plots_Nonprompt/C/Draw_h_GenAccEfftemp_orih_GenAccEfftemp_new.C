void Draw_h_GenAccEfftemp_orih_GenAccEfftemp_new()
{
//=========Macro generated from canvas: 15/15
//=========  (Wed Apr 11 18:26:30 2018) by ROOT version6.02/13
   TCanvas *15 = new TCanvas("15", "15",0,0,600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   15->Range(5.571429,0.7905301,8.428571,0.9484997);
   15->SetFillColor(0);
   15->SetBorderMode(0);
   15->SetBorderSize(0);
   15->SetTickx(1);
   15->SetTicky(1);
   15->SetLeftMargin(0.15);
   15->SetRightMargin(0.15);
   15->SetTopMargin(0.15);
   15->SetBottomMargin(0.15);
   15->SetFrameLineColor(0);
   15->SetFrameBorderMode(0);
   15->SetFrameLineColor(0);
   15->SetFrameBorderMode(0);
   
   TH1D *h_GenAccEfftemp_ori85 = new TH1D("h_GenAccEfftemp_ori85","h_GenAcctemp_ori",10,6,8);
   h_GenAccEfftemp_ori85->SetBinContent(1,0.8336222);
   h_GenAccEfftemp_ori85->SetBinContent(2,0.8311321);
   h_GenAccEfftemp_ori85->SetBinContent(3,0.8333333);
   h_GenAccEfftemp_ori85->SetBinContent(4,0.8390663);
   h_GenAccEfftemp_ori85->SetBinContent(5,0.8558201);
   h_GenAccEfftemp_ori85->SetBinContent(6,0.8625954);
   h_GenAccEfftemp_ori85->SetBinContent(7,0.8458961);
   h_GenAccEfftemp_ori85->SetBinContent(8,0.9063745);
   h_GenAccEfftemp_ori85->SetBinContent(9,0.8734694);
   h_GenAccEfftemp_ori85->SetBinContent(10,0.8958838);
   h_GenAccEfftemp_ori85->SetBinError(1,0.010963);
   h_GenAccEfftemp_ori85->SetBinError(2,0.01150683);
   h_GenAccEfftemp_ori85->SetBinError(3,0.01254872);
   h_GenAccEfftemp_ori85->SetBinError(4,0.01287981);
   h_GenAccEfftemp_ori85->SetBinError(5,0.01277564);
   h_GenAccEfftemp_ori85->SetBinError(6,0.0134519);
   h_GenAccEfftemp_ori85->SetBinError(7,0.01477672);
   h_GenAccEfftemp_ori85->SetBinError(8,0.01300167);
   h_GenAccEfftemp_ori85->SetBinError(9,0.0150184);
   h_GenAccEfftemp_ori85->SetBinError(10,0.01502831);
   h_GenAccEfftemp_ori85->SetMinimum(0.8142256);
   h_GenAccEfftemp_ori85->SetMaximum(0.9248043);
   h_GenAccEfftemp_ori85->SetEntries(7323);
   h_GenAccEfftemp_ori85->SetStats(0);
   h_GenAccEfftemp_ori85->SetFillColor(1);
   h_GenAccEfftemp_ori85->SetFillStyle(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000082");
   h_GenAccEfftemp_ori85->SetLineColor(ci);
   h_GenAccEfftemp_ori85->SetLineStyle(0);

   ci = TColor::GetColor("#000082");
   h_GenAccEfftemp_ori85->SetMarkerColor(ci);
   h_GenAccEfftemp_ori85->SetMarkerStyle(20);
   h_GenAccEfftemp_ori85->GetXaxis()->SetLabelFont(43);
   h_GenAccEfftemp_ori85->GetXaxis()->SetLabelOffset(0.01);
   h_GenAccEfftemp_ori85->GetXaxis()->SetLabelSize(20);
   h_GenAccEfftemp_ori85->GetXaxis()->SetTitleSize(25);
   h_GenAccEfftemp_ori85->GetXaxis()->SetTitleOffset(1.5);
   h_GenAccEfftemp_ori85->GetXaxis()->SetTitleFont(43);
   h_GenAccEfftemp_ori85->GetYaxis()->SetLabelFont(43);
   h_GenAccEfftemp_ori85->GetYaxis()->SetLabelOffset(0.01);
   h_GenAccEfftemp_ori85->GetYaxis()->SetLabelSize(20);
   h_GenAccEfftemp_ori85->GetYaxis()->SetTitleSize(25);
   h_GenAccEfftemp_ori85->GetYaxis()->SetTitleOffset(1.7);
   h_GenAccEfftemp_ori85->GetYaxis()->SetTitleFont(43);
   h_GenAccEfftemp_ori85->GetZaxis()->SetLabelFont(43);
   h_GenAccEfftemp_ori85->GetZaxis()->SetLabelSize(20);
   h_GenAccEfftemp_ori85->GetZaxis()->SetTitleSize(25);
   h_GenAccEfftemp_ori85->GetZaxis()->SetTitleFont(43);
   h_GenAccEfftemp_ori85->Draw("E1");
   
   TH1D *h_GenAccEfftemp_new86 = new TH1D("h_GenAccEfftemp_new86","h_GenAcctemp_new",10,6,8);
   h_GenAccEfftemp_new86->SetBinContent(1,0.8240049);
   h_GenAccEfftemp_new86->SetBinContent(2,0.8142256);
   h_GenAccEfftemp_new86->SetBinContent(3,0.82214);
   h_GenAccEfftemp_new86->SetBinContent(4,0.8279969);
   h_GenAccEfftemp_new86->SetBinContent(5,0.8679869);
   h_GenAccEfftemp_new86->SetBinContent(6,0.8682656);
   h_GenAccEfftemp_new86->SetBinContent(7,0.8648249);
   h_GenAccEfftemp_new86->SetBinContent(8,0.8765661);
   h_GenAccEfftemp_new86->SetBinContent(9,0.9014437);
   h_GenAccEfftemp_new86->SetBinContent(10,0.8741055);
   h_GenAccEfftemp_new86->SetBinError(1,0.01153667);
   h_GenAccEfftemp_new86->SetBinError(2,0.01234196);
   h_GenAccEfftemp_new86->SetBinError(3,0.01258077);
   h_GenAccEfftemp_new86->SetBinError(4,0.01261187);
   h_GenAccEfftemp_new86->SetBinError(5,0.0121972);
   h_GenAccEfftemp_new86->SetBinError(6,0.01191956);
   h_GenAccEfftemp_new86->SetBinError(7,0.01443378);
   h_GenAccEfftemp_new86->SetBinError(8,0.01264159);
   h_GenAccEfftemp_new86->SetBinError(9,0.01213916);
   h_GenAccEfftemp_new86->SetBinError(10,0.01382916);
   h_GenAccEfftemp_new86->SetEntries(11771);
   h_GenAccEfftemp_new86->SetStats(0);
   h_GenAccEfftemp_new86->SetFillColor(1);
   h_GenAccEfftemp_new86->SetFillStyle(0);

   ci = TColor::GetColor("#87ec76");
   h_GenAccEfftemp_new86->SetLineColor(ci);
   h_GenAccEfftemp_new86->SetLineStyle(0);

   ci = TColor::GetColor("#87ec76");
   h_GenAccEfftemp_new86->SetMarkerColor(ci);
   h_GenAccEfftemp_new86->SetMarkerStyle(20);
   h_GenAccEfftemp_new86->GetXaxis()->SetLabelFont(43);
   h_GenAccEfftemp_new86->GetXaxis()->SetLabelOffset(0.01);
   h_GenAccEfftemp_new86->GetXaxis()->SetLabelSize(20);
   h_GenAccEfftemp_new86->GetXaxis()->SetTitleSize(25);
   h_GenAccEfftemp_new86->GetXaxis()->SetTitleOffset(1.5);
   h_GenAccEfftemp_new86->GetXaxis()->SetTitleFont(43);
   h_GenAccEfftemp_new86->GetYaxis()->SetLabelFont(43);
   h_GenAccEfftemp_new86->GetYaxis()->SetLabelOffset(0.01);
   h_GenAccEfftemp_new86->GetYaxis()->SetLabelSize(20);
   h_GenAccEfftemp_new86->GetYaxis()->SetTitleSize(25);
   h_GenAccEfftemp_new86->GetYaxis()->SetTitleOffset(1.7);
   h_GenAccEfftemp_new86->GetYaxis()->SetTitleFont(43);
   h_GenAccEfftemp_new86->GetZaxis()->SetLabelFont(43);
   h_GenAccEfftemp_new86->GetZaxis()->SetLabelSize(20);
   h_GenAccEfftemp_new86->GetZaxis()->SetTitleSize(25);
   h_GenAccEfftemp_new86->GetZaxis()->SetTitleFont(43);
   h_GenAccEfftemp_new86->Draw("E1same");
   TLatex *   tex = new TLatex(0.55,0.79,"");
tex->SetNDC();
   tex->SetTextFont(43);
   tex->SetTextSize(20);
   tex->SetLineWidth(2);
   tex->Draw();
   
   TLegend *leg = new TLegend(0.48,0.6,0.85,0.8,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(43);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(4000);
   TLegendEntry *entry=leg->AddEntry("h_GenAccEfftemp_ori","h_GenAcctemp_ori","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#000082");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("h_GenAccEfftemp_new","h_GenAcctemp_new","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#87ec76");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   15->Modified();
   15->cd();
   15->SetSelected(15);
}
