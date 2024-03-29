void Compare_h_ori_3_h_new_3()
{
//=========Macro generated from canvas: Compare_h_ori_3_h_new_3/Gtk1pt_oriGtk1pt_new
//=========  (Wed Apr 11 18:26:27 2018) by ROOT version6.02/13
   TCanvas *Compare_h_ori_3_h_new_3 = new TCanvas("Compare_h_ori_3_h_new_3", "Gtk1pt_oriGtk1pt_new",0,0,600,700);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   Compare_h_ori_3_h_new_3->Range(0,0,1,1);
   Compare_h_ori_3_h_new_3->SetFillColor(0);
   Compare_h_ori_3_h_new_3->SetBorderMode(0);
   Compare_h_ori_3_h_new_3->SetBorderSize(0);
   Compare_h_ori_3_h_new_3->SetTickx(1);
   Compare_h_ori_3_h_new_3->SetTicky(1);
   Compare_h_ori_3_h_new_3->SetLeftMargin(0.15);
   Compare_h_ori_3_h_new_3->SetRightMargin(0.15);
   Compare_h_ori_3_h_new_3->SetTopMargin(0.15);
   Compare_h_ori_3_h_new_3->SetBottomMargin(0.15);
   Compare_h_ori_3_h_new_3->SetFrameLineColor(0);
   Compare_h_ori_3_h_new_3->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: pad1
   TPad *pad1 = new TPad("pad1", "pad1",0,0.35,1,1);
   pad1->Draw();
   pad1->cd();
   pad1->Range(-2.142857,-0.001016655,12.14286,0.1575738);
   pad1->SetFillColor(0);
   pad1->SetBorderMode(0);
   pad1->SetBorderSize(0);
   pad1->SetTickx(1);
   pad1->SetTicky(1);
   pad1->SetLeftMargin(0.15);
   pad1->SetRightMargin(0.15);
   pad1->SetTopMargin(0.15);
   pad1->SetBottomMargin(0.02);
   pad1->SetFrameLineColor(0);
   pad1->SetFrameBorderMode(0);
   pad1->SetFrameLineColor(0);
   pad1->SetFrameBorderMode(0);
   
   TH1F *h_new_315 = new TH1F("h_new_315","Gtk1pt_new",20,0,10);
   h_new_315->SetBinContent(0,3.139087);
   h_new_315->SetBinContent(1,0.02125203);
   h_new_315->SetBinContent(2,0.03103538);
   h_new_315->SetBinContent(3,0.02613379);
   h_new_315->SetBinContent(4,0.05363607);
   h_new_315->SetBinContent(5,0.07339138);
   h_new_315->SetBinContent(6,0.07951222);
   h_new_315->SetBinContent(7,0.08576944);
   h_new_315->SetBinContent(8,0.08723883);
   h_new_315->SetBinContent(9,0.08368294);
   h_new_315->SetBinContent(10,0.08071031);
   h_new_315->SetBinContent(11,0.08254959);
   h_new_315->SetBinContent(12,0.08088773);
   h_new_315->SetBinContent(13,0.05895385);
   h_new_315->SetBinContent(14,0.04598147);
   h_new_315->SetBinContent(15,0.03270571);
   h_new_315->SetBinContent(16,0.0247635);
   h_new_315->SetBinContent(17,0.01772635);
   h_new_315->SetBinContent(18,0.01506399);
   h_new_315->SetBinContent(19,0.01013091);
   h_new_315->SetBinContent(20,0.008874485);
   h_new_315->SetBinContent(21,0.04084381);
   h_new_315->SetBinError(0,0.01172593);
   h_new_315->SetBinError(1,0.001067441);
   h_new_315->SetBinError(2,0.001280365);
   h_new_315->SetBinError(3,0.001115208);
   h_new_315->SetBinError(4,0.00168736);
   h_new_315->SetBinError(5,0.001870935);
   h_new_315->SetBinError(6,0.001926794);
   h_new_315->SetBinError(7,0.002089507);
   h_new_315->SetBinError(8,0.002006772);
   h_new_315->SetBinError(9,0.001970186);
   h_new_315->SetBinError(10,0.001895052);
   h_new_315->SetBinError(11,0.001950576);
   h_new_315->SetBinError(12,0.001904106);
   h_new_315->SetBinError(13,0.001539585);
   h_new_315->SetBinError(14,0.001317273);
   h_new_315->SetBinError(15,0.001086232);
   h_new_315->SetBinError(16,0.0008847138);
   h_new_315->SetBinError(17,0.0007386104);
   h_new_315->SetBinError(18,0.0006370623);
   h_new_315->SetBinError(19,0.00049268);
   h_new_315->SetBinError(20,0.0004512331);
   h_new_315->SetBinError(21,0.0008141919);
   h_new_315->SetMaximum(0.1337852);
   h_new_315->SetEntries(174188);
   h_new_315->SetStats(0);
   h_new_315->SetFillColor(1);
   h_new_315->SetFillStyle(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#87ec76");
   h_new_315->SetLineColor(ci);
   h_new_315->SetLineStyle(0);

   ci = TColor::GetColor("#87ec76");
   h_new_315->SetMarkerColor(ci);
   h_new_315->SetMarkerStyle(20);
   h_new_315->GetXaxis()->SetLabelFont(43);
   h_new_315->GetXaxis()->SetLabelOffset(0.01);
   h_new_315->GetXaxis()->SetLabelSize(0);
   h_new_315->GetXaxis()->SetTitleSize(25);
   h_new_315->GetXaxis()->SetTitleOffset(1.5);
   h_new_315->GetXaxis()->SetTitleFont(43);
   h_new_315->GetYaxis()->SetLabelFont(43);
   h_new_315->GetYaxis()->SetLabelOffset(0.01);
   h_new_315->GetYaxis()->SetLabelSize(20);
   h_new_315->GetYaxis()->SetTitleSize(25);
   h_new_315->GetYaxis()->SetTitleOffset(1.7);
   h_new_315->GetYaxis()->SetTitleFont(43);
   h_new_315->GetZaxis()->SetLabelFont(43);
   h_new_315->GetZaxis()->SetLabelSize(20);
   h_new_315->GetZaxis()->SetTitleSize(25);
   h_new_315->GetZaxis()->SetTitleFont(43);
   h_new_315->Draw("");
   
   TH1F *h_new_316 = new TH1F("h_new_316","Gtk1pt_new",20,0,10);
   h_new_316->SetBinContent(0,3.139087);
   h_new_316->SetBinContent(1,0.02125203);
   h_new_316->SetBinContent(2,0.03103538);
   h_new_316->SetBinContent(3,0.02613379);
   h_new_316->SetBinContent(4,0.05363607);
   h_new_316->SetBinContent(5,0.07339138);
   h_new_316->SetBinContent(6,0.07951222);
   h_new_316->SetBinContent(7,0.08576944);
   h_new_316->SetBinContent(8,0.08723883);
   h_new_316->SetBinContent(9,0.08368294);
   h_new_316->SetBinContent(10,0.08071031);
   h_new_316->SetBinContent(11,0.08254959);
   h_new_316->SetBinContent(12,0.08088773);
   h_new_316->SetBinContent(13,0.05895385);
   h_new_316->SetBinContent(14,0.04598147);
   h_new_316->SetBinContent(15,0.03270571);
   h_new_316->SetBinContent(16,0.0247635);
   h_new_316->SetBinContent(17,0.01772635);
   h_new_316->SetBinContent(18,0.01506399);
   h_new_316->SetBinContent(19,0.01013091);
   h_new_316->SetBinContent(20,0.008874485);
   h_new_316->SetBinContent(21,0.04084381);
   h_new_316->SetBinError(0,0.01172593);
   h_new_316->SetBinError(1,0.001067441);
   h_new_316->SetBinError(2,0.001280365);
   h_new_316->SetBinError(3,0.001115208);
   h_new_316->SetBinError(4,0.00168736);
   h_new_316->SetBinError(5,0.001870935);
   h_new_316->SetBinError(6,0.001926794);
   h_new_316->SetBinError(7,0.002089507);
   h_new_316->SetBinError(8,0.002006772);
   h_new_316->SetBinError(9,0.001970186);
   h_new_316->SetBinError(10,0.001895052);
   h_new_316->SetBinError(11,0.001950576);
   h_new_316->SetBinError(12,0.001904106);
   h_new_316->SetBinError(13,0.001539585);
   h_new_316->SetBinError(14,0.001317273);
   h_new_316->SetBinError(15,0.001086232);
   h_new_316->SetBinError(16,0.0008847138);
   h_new_316->SetBinError(17,0.0007386104);
   h_new_316->SetBinError(18,0.0006370623);
   h_new_316->SetBinError(19,0.00049268);
   h_new_316->SetBinError(20,0.0004512331);
   h_new_316->SetBinError(21,0.0008141919);
   h_new_316->SetMaximum(0.1337852);
   h_new_316->SetEntries(174188);
   h_new_316->SetStats(0);
   h_new_316->SetFillColor(1);
   h_new_316->SetFillStyle(0);

   ci = TColor::GetColor("#87ec76");
   h_new_316->SetLineColor(ci);
   h_new_316->SetLineStyle(0);

   ci = TColor::GetColor("#87ec76");
   h_new_316->SetMarkerColor(ci);
   h_new_316->SetMarkerStyle(20);
   h_new_316->GetXaxis()->SetLabelFont(43);
   h_new_316->GetXaxis()->SetLabelOffset(0.01);
   h_new_316->GetXaxis()->SetLabelSize(0);
   h_new_316->GetXaxis()->SetTitleSize(25);
   h_new_316->GetXaxis()->SetTitleOffset(1.5);
   h_new_316->GetXaxis()->SetTitleFont(43);
   h_new_316->GetYaxis()->SetLabelFont(43);
   h_new_316->GetYaxis()->SetLabelOffset(0.01);
   h_new_316->GetYaxis()->SetLabelSize(20);
   h_new_316->GetYaxis()->SetTitleSize(25);
   h_new_316->GetYaxis()->SetTitleOffset(1.7);
   h_new_316->GetYaxis()->SetTitleFont(43);
   h_new_316->GetZaxis()->SetLabelFont(43);
   h_new_316->GetZaxis()->SetLabelSize(20);
   h_new_316->GetZaxis()->SetTitleSize(25);
   h_new_316->GetZaxis()->SetTitleFont(43);
   h_new_316->Draw("same");
   
   TH1F *h_ori_317 = new TH1F("h_ori_317","Gtk1pt_ori",20,0,10);
   h_ori_317->SetBinContent(0,3.145416);
   h_ori_317->SetBinContent(1,0.0201113);
   h_ori_317->SetBinContent(2,0.03163136);
   h_ori_317->SetBinContent(3,0.02489505);
   h_ori_317->SetBinContent(4,0.05164503);
   h_ori_317->SetBinContent(5,0.07488041);
   h_ori_317->SetBinContent(6,0.07780924);
   h_ori_317->SetBinContent(7,0.08342282);
   h_ori_317->SetBinContent(8,0.08337401);
   h_ori_317->SetBinContent(9,0.08327638);
   h_ori_317->SetBinContent(10,0.08464317);
   h_ori_317->SetBinContent(11,0.08669335);
   h_ori_317->SetBinContent(12,0.07824856);
   h_ori_317->SetBinContent(13,0.06067558);
   h_ori_317->SetBinContent(14,0.04549449);
   h_ori_317->SetBinContent(15,0.03314459);
   h_ori_317->SetBinContent(16,0.02513912);
   h_ori_317->SetBinContent(17,0.01937909);
   h_ori_317->SetBinContent(18,0.01537635);
   h_ori_317->SetBinContent(19,0.01117837);
   h_ori_317->SetBinContent(20,0.008981744);
   h_ori_317->SetBinContent(21,0.03807478);
   h_ori_317->SetBinError(0,0.01239112);
   h_ori_317->SetBinError(1,0.0009908124);
   h_ori_317->SetBinError(2,0.001242597);
   h_ori_317->SetBinError(3,0.001102371);
   h_ori_317->SetBinError(4,0.001587763);
   h_ori_317->SetBinError(5,0.001911857);
   h_ori_317->SetBinError(6,0.001948888);
   h_ori_317->SetBinError(7,0.002017966);
   h_ori_317->SetBinError(8,0.002017376);
   h_ori_317->SetBinError(9,0.002016194);
   h_ori_317->SetBinError(10,0.002032672);
   h_ori_317->SetBinError(11,0.002057142);
   h_ori_317->SetBinError(12,0.001954383);
   h_ori_317->SetBinError(13,0.00172099);
   h_ori_317->SetBinError(14,0.001490221);
   h_ori_317->SetBinError(15,0.001271972);
   h_ori_317->SetBinError(16,0.001107762);
   h_ori_317->SetBinError(17,0.0009726086);
   h_ori_317->SetBinError(18,0.0008663594);
   h_ori_317->SetBinError(19,0.0007386872);
   h_ori_317->SetBinError(20,0.0006621429);
   h_ori_317->SetBinError(21,0.001363296);
   h_ori_317->SetMinimum(0.008874485);
   h_ori_317->SetMaximum(0.1337852);
   h_ori_317->SetEntries(85703);
   h_ori_317->SetStats(0);
   h_ori_317->SetFillColor(1);
   h_ori_317->SetFillStyle(0);

   ci = TColor::GetColor("#000082");
   h_ori_317->SetLineColor(ci);
   h_ori_317->SetLineStyle(0);

   ci = TColor::GetColor("#000082");
   h_ori_317->SetMarkerColor(ci);
   h_ori_317->SetMarkerStyle(20);
   h_ori_317->GetXaxis()->SetLabelFont(43);
   h_ori_317->GetXaxis()->SetLabelOffset(0.01);
   h_ori_317->GetXaxis()->SetLabelSize(20);
   h_ori_317->GetXaxis()->SetTitleSize(25);
   h_ori_317->GetXaxis()->SetTitleOffset(1.5);
   h_ori_317->GetXaxis()->SetTitleFont(43);
   h_ori_317->GetYaxis()->SetLabelFont(43);
   h_ori_317->GetYaxis()->SetLabelOffset(0.01);
   h_ori_317->GetYaxis()->SetLabelSize(20);
   h_ori_317->GetYaxis()->SetTitleSize(25);
   h_ori_317->GetYaxis()->SetTitleOffset(1.7);
   h_ori_317->GetYaxis()->SetTitleFont(43);
   h_ori_317->GetZaxis()->SetLabelFont(43);
   h_ori_317->GetZaxis()->SetLabelSize(20);
   h_ori_317->GetZaxis()->SetTitleSize(25);
   h_ori_317->GetZaxis()->SetTitleFont(43);
   h_ori_317->Draw("same");
   
   TLegend *leg = new TLegend(0.57,0.68,0.84,0.88,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(43);
   leg->SetTextSize(20);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(4000);
   TLegendEntry *entry=leg->AddEntry("NULL","","h");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(43);
   entry=leg->AddEntry("h_ori_3","Gtk1pt_ori","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#000082");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(43);
   entry=leg->AddEntry("h_new_3","Gtk1pt_new","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#87ec76");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(43);
   leg->Draw();
   TLatex *   tex = new TLatex(0.2,0.77,"");
tex->SetNDC();
   tex->SetTextFont(43);
   tex->SetTextSize(20);
   tex->SetLineWidth(2);
   tex->Draw();
   pad1->Modified();
   Compare_h_ori_3_h_new_3->cd();
  
// ------------>Primitives in pad: pad2
   TPad *pad2 = new TPad("pad2", "pad2",0,0,1,0.35);
   pad2->Draw();
   pad2->cd();
   pad2->Range(-2.142857,-1.228571,12.14286,2.2);
   pad2->SetFillColor(0);
   pad2->SetBorderMode(0);
   pad2->SetBorderSize(0);
   pad2->SetTickx(1);
   pad2->SetTicky(1);
   pad2->SetLeftMargin(0.15);
   pad2->SetRightMargin(0.15);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.3);
   pad2->SetFrameLineColor(0);
   pad2->SetFrameBorderMode(0);
   pad2->SetFrameLineColor(0);
   pad2->SetFrameBorderMode(0);
   
   TH1F *h_ori_318 = new TH1F("h_ori_318","Gtk1pt_ori",20,0,10);
   h_ori_318->SetBinContent(0,1.002016);
   h_ori_318->SetBinContent(1,0.9463233);
   h_ori_318->SetBinContent(2,1.019203);
   h_ori_318->SetBinContent(3,0.9525999);
   h_ori_318->SetBinContent(4,0.9628786);
   h_ori_318->SetBinContent(5,1.020289);
   h_ori_318->SetBinContent(6,0.9785821);
   h_ori_318->SetBinContent(7,0.9726405);
   h_ori_318->SetBinContent(8,0.9556984);
   h_ori_318->SetBinContent(9,0.9951417);
   h_ori_318->SetBinContent(10,1.048728);
   h_ori_318->SetBinContent(11,1.050197);
   h_ori_318->SetBinContent(12,0.9673724);
   h_ori_318->SetBinContent(13,1.029205);
   h_ori_318->SetBinContent(14,0.9894092);
   h_ori_318->SetBinContent(15,1.013419);
   h_ori_318->SetBinContent(16,1.015168);
   h_ori_318->SetBinContent(17,1.093236);
   h_ori_318->SetBinContent(18,1.020736);
   h_ori_318->SetBinContent(19,1.103392);
   h_ori_318->SetBinContent(20,1.012086);
   h_ori_318->SetBinContent(21,0.9322046);
   h_ori_318->SetBinError(0,0.005439821);
   h_ori_318->SetBinError(1,0.06657981);
   h_ori_318->SetBinError(2,0.05806048);
   h_ori_318->SetBinError(3,0.05858119);
   h_ori_318->SetBinError(4,0.04235434);
   h_ori_318->SetBinError(5,0.03681196);
   h_ori_318->SetBinError(6,0.03410433);
   h_ori_318->SetBinError(7,0.03339204);
   h_ori_318->SetBinError(8,0.03190698);
   h_ori_318->SetBinError(9,0.03360665);
   h_ori_318->SetBinError(10,0.03522222);
   h_ori_318->SetBinError(11,0.03516827);
   h_ori_318->SetBinError(12,0.03320169);
   h_ori_318->SetBinError(13,0.0396812);
   h_ori_318->SetBinError(14,0.04305538);
   h_ori_318->SetBinError(15,0.0514335);
   h_ori_318->SetBinError(16,0.05758905);
   h_ori_318->SetBinError(17,0.07131269);
   h_ori_318->SetBinError(18,0.07190999);
   h_ori_318->SetBinError(19,0.09053072);
   h_ori_318->SetBinError(20,0.09063744);
   h_ori_318->SetBinError(21,0.0382025);
   h_ori_318->SetMinimum(-0.2);
   h_ori_318->SetMaximum(2.2);
   h_ori_318->SetEntries(6944.209);
   h_ori_318->SetStats(0);
   h_ori_318->SetFillColor(1);
   h_ori_318->SetFillStyle(0);

   ci = TColor::GetColor("#333333");
   h_ori_318->SetLineColor(ci);

   ci = TColor::GetColor("#333333");
   h_ori_318->SetMarkerColor(ci);
   h_ori_318->SetMarkerStyle(21);
   h_ori_318->GetXaxis()->SetTitle("Gtk1pt");
   h_ori_318->GetXaxis()->CenterTitle(true);
   h_ori_318->GetXaxis()->SetLabelFont(43);
   h_ori_318->GetXaxis()->SetLabelOffset(0.01);
   h_ori_318->GetXaxis()->SetLabelSize(20);
   h_ori_318->GetXaxis()->SetTitleSize(25);
   h_ori_318->GetXaxis()->SetTitleOffset(3.5);
   h_ori_318->GetXaxis()->SetTitleFont(43);
   h_ori_318->GetYaxis()->SetTitle("ratio");
   h_ori_318->GetYaxis()->CenterTitle(true);
   h_ori_318->GetYaxis()->SetLabelFont(43);
   h_ori_318->GetYaxis()->SetLabelOffset(0.01);
   h_ori_318->GetYaxis()->SetLabelSize(20);
   h_ori_318->GetYaxis()->SetTitleSize(25);
   h_ori_318->GetYaxis()->SetTitleOffset(1.7);
   h_ori_318->GetYaxis()->SetTitleFont(43);
   h_ori_318->GetZaxis()->SetLabelFont(43);
   h_ori_318->GetZaxis()->SetLabelSize(20);
   h_ori_318->GetZaxis()->SetTitleSize(25);
   h_ori_318->GetZaxis()->SetTitleFont(43);
   h_ori_318->Draw("ep");
   TLine *line = new TLine(0,1,10,1);
   line->SetLineStyle(2);
   line->Draw();
   pad2->Modified();
   Compare_h_ori_3_h_new_3->cd();
   Compare_h_ori_3_h_new_3->Modified();
   Compare_h_ori_3_h_new_3->cd();
   Compare_h_ori_3_h_new_3->SetSelected(Compare_h_ori_3_h_new_3);
}
