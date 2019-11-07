void Draw_h_RecoNormEfftemp_orih_RecoNormEfftemp_new()
{
//=========Macro generated from canvas: 16/16
//=========  (Wed Apr 11 18:26:30 2018) by ROOT version6.02/13
   TCanvas *16 = new TCanvas("16", "16",0,0,600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   16->Range(5.571429,0.2527679,8.428571,0.527671);
   16->SetFillColor(0);
   16->SetBorderMode(0);
   16->SetBorderSize(0);
   16->SetTickx(1);
   16->SetTicky(1);
   16->SetLeftMargin(0.15);
   16->SetRightMargin(0.15);
   16->SetTopMargin(0.15);
   16->SetBottomMargin(0.15);
   16->SetFrameLineColor(0);
   16->SetFrameBorderMode(0);
   16->SetFrameLineColor(0);
   16->SetFrameBorderMode(0);
   
   TH1D *h_RecoNormEfftemp_ori91 = new TH1D("h_RecoNormEfftemp_ori91","h_RecoNormtemp_ori",10,6,8);
   h_RecoNormEfftemp_ori91->SetBinContent(1,0.3232236);
   h_RecoNormEfftemp_ori91->SetBinContent(2,0.3433962);
   h_RecoNormEfftemp_ori91->SetBinContent(3,0.3344671);
   h_RecoNormEfftemp_ori91->SetBinContent(4,0.3427518);
   h_RecoNormEfftemp_ori91->SetBinContent(5,0.4126984);
   h_RecoNormEfftemp_ori91->SetBinContent(6,0.4045802);
   h_RecoNormEfftemp_ori91->SetBinContent(7,0.3819095);
   h_RecoNormEfftemp_ori91->SetBinContent(8,0.440239);
   h_RecoNormEfftemp_ori91->SetBinContent(9,0.4367347);
   h_RecoNormEfftemp_ori91->SetBinContent(10,0.4455206);
   h_RecoNormEfftemp_ori91->SetBinError(1,0.01376801);
   h_RecoNormEfftemp_ori91->SetBinError(2,0.01458467);
   h_RecoNormEfftemp_ori91->SetBinError(3,0.01588646);
   h_RecoNormEfftemp_ori91->SetBinError(4,0.01663575);
   h_RecoNormEfftemp_ori91->SetBinError(5,0.01790548);
   h_RecoNormEfftemp_ori91->SetBinError(6,0.01917756);
   h_RecoNormEfftemp_ori91->SetBinError(7,0.0198847);
   h_RecoNormEfftemp_ori91->SetBinError(8,0.02215612);
   h_RecoNormEfftemp_ori91->SetBinError(9,0.02240615);
   h_RecoNormEfftemp_ori91->SetBinError(10,0.02445691);
   h_RecoNormEfftemp_ori91->SetMinimum(0.2940033);
   h_RecoNormEfftemp_ori91->SetMaximum(0.4864356);
   h_RecoNormEfftemp_ori91->SetEntries(7323);
   h_RecoNormEfftemp_ori91->SetStats(0);
   h_RecoNormEfftemp_ori91->SetFillColor(1);
   h_RecoNormEfftemp_ori91->SetFillStyle(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000082");
   h_RecoNormEfftemp_ori91->SetLineColor(ci);
   h_RecoNormEfftemp_ori91->SetLineStyle(0);

   ci = TColor::GetColor("#000082");
   h_RecoNormEfftemp_ori91->SetMarkerColor(ci);
   h_RecoNormEfftemp_ori91->SetMarkerStyle(20);
   h_RecoNormEfftemp_ori91->GetXaxis()->SetLabelFont(43);
   h_RecoNormEfftemp_ori91->GetXaxis()->SetLabelOffset(0.01);
   h_RecoNormEfftemp_ori91->GetXaxis()->SetLabelSize(20);
   h_RecoNormEfftemp_ori91->GetXaxis()->SetTitleSize(25);
   h_RecoNormEfftemp_ori91->GetXaxis()->SetTitleOffset(1.5);
   h_RecoNormEfftemp_ori91->GetXaxis()->SetTitleFont(43);
   h_RecoNormEfftemp_ori91->GetYaxis()->SetLabelFont(43);
   h_RecoNormEfftemp_ori91->GetYaxis()->SetLabelOffset(0.01);
   h_RecoNormEfftemp_ori91->GetYaxis()->SetLabelSize(20);
   h_RecoNormEfftemp_ori91->GetYaxis()->SetTitleSize(25);
   h_RecoNormEfftemp_ori91->GetYaxis()->SetTitleOffset(1.7);
   h_RecoNormEfftemp_ori91->GetYaxis()->SetTitleFont(43);
   h_RecoNormEfftemp_ori91->GetZaxis()->SetLabelFont(43);
   h_RecoNormEfftemp_ori91->GetZaxis()->SetLabelSize(20);
   h_RecoNormEfftemp_ori91->GetZaxis()->SetTitleSize(25);
   h_RecoNormEfftemp_ori91->GetZaxis()->SetTitleFont(43);
   h_RecoNormEfftemp_ori91->Draw("E1");
   
   TH1D *h_RecoNormEfftemp_new92 = new TH1D("h_RecoNormEfftemp_new92","h_RecoNormtemp_new",10,6,8);
   h_RecoNormEfftemp_new92->SetBinContent(1,0.2940033);
   h_RecoNormEfftemp_new92->SetBinContent(2,0.3336796);
   h_RecoNormEfftemp_new92->SetBinContent(3,0.3168845);
   h_RecoNormEfftemp_new92->SetBinContent(4,0.3832166);
   h_RecoNormEfftemp_new92->SetBinContent(5,0.4172689);
   h_RecoNormEfftemp_new92->SetBinContent(6,0.4136793);
   h_RecoNormEfftemp_new92->SetBinContent(7,0.4409714);
   h_RecoNormEfftemp_new92->SetBinContent(8,0.4529542);
   h_RecoNormEfftemp_new92->SetBinContent(9,0.4197908);
   h_RecoNormEfftemp_new92->SetBinContent(10,0.4543635);
   h_RecoNormEfftemp_new92->SetBinError(1,0.01350804);
   h_RecoNormEfftemp_new92->SetBinError(2,0.0148923);
   h_RecoNormEfftemp_new92->SetBinError(3,0.01503586);
   h_RecoNormEfftemp_new92->SetBinError(4,0.01733445);
   h_RecoNormEfftemp_new92->SetBinError(5,0.01787033);
   h_RecoNormEfftemp_new92->SetBinError(6,0.01804145);
   h_RecoNormEfftemp_new92->SetBinError(7,0.01953653);
   h_RecoNormEfftemp_new92->SetBinError(8,0.01967643);
   h_RecoNormEfftemp_new92->SetBinError(9,0.02038585);
   h_RecoNormEfftemp_new92->SetBinError(10,0.02100875);
   h_RecoNormEfftemp_new92->SetEntries(11771);
   h_RecoNormEfftemp_new92->SetStats(0);
   h_RecoNormEfftemp_new92->SetFillColor(1);
   h_RecoNormEfftemp_new92->SetFillStyle(0);

   ci = TColor::GetColor("#87ec76");
   h_RecoNormEfftemp_new92->SetLineColor(ci);
   h_RecoNormEfftemp_new92->SetLineStyle(0);

   ci = TColor::GetColor("#87ec76");
   h_RecoNormEfftemp_new92->SetMarkerColor(ci);
   h_RecoNormEfftemp_new92->SetMarkerStyle(20);
   h_RecoNormEfftemp_new92->GetXaxis()->SetLabelFont(43);
   h_RecoNormEfftemp_new92->GetXaxis()->SetLabelOffset(0.01);
   h_RecoNormEfftemp_new92->GetXaxis()->SetLabelSize(20);
   h_RecoNormEfftemp_new92->GetXaxis()->SetTitleSize(25);
   h_RecoNormEfftemp_new92->GetXaxis()->SetTitleOffset(1.5);
   h_RecoNormEfftemp_new92->GetXaxis()->SetTitleFont(43);
   h_RecoNormEfftemp_new92->GetYaxis()->SetLabelFont(43);
   h_RecoNormEfftemp_new92->GetYaxis()->SetLabelOffset(0.01);
   h_RecoNormEfftemp_new92->GetYaxis()->SetLabelSize(20);
   h_RecoNormEfftemp_new92->GetYaxis()->SetTitleSize(25);
   h_RecoNormEfftemp_new92->GetYaxis()->SetTitleOffset(1.7);
   h_RecoNormEfftemp_new92->GetYaxis()->SetTitleFont(43);
   h_RecoNormEfftemp_new92->GetZaxis()->SetLabelFont(43);
   h_RecoNormEfftemp_new92->GetZaxis()->SetLabelSize(20);
   h_RecoNormEfftemp_new92->GetZaxis()->SetTitleSize(25);
   h_RecoNormEfftemp_new92->GetZaxis()->SetTitleFont(43);
   h_RecoNormEfftemp_new92->Draw("E1same");
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
   TLegendEntry *entry=leg->AddEntry("h_RecoNormEfftemp_ori","h_RecoNormtemp_ori","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#000082");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("h_RecoNormEfftemp_new","h_RecoNormtemp_new","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#87ec76");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   16->Modified();
   16->cd();
   16->SetSelected(16);
}
