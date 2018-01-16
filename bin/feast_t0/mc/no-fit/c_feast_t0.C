void c_feast_t0()
{
//=========Macro generated from canvas: c_feast_t0/c_feast_t0
//=========  (Tue Jan 16 17:30:12 2018) by ROOT version6.08/06
   TCanvas *c_feast_t0 = new TCanvas("c_feast_t0", "c_feast_t0",0,0,800,600);
   c_feast_t0->SetHighLightColor(2);
   c_feast_t0->Range(0,0,1,1);
   c_feast_t0->SetFillColor(0);
   c_feast_t0->SetBorderMode(0);
   c_feast_t0->SetBorderSize(2);
   c_feast_t0->SetFrameBorderMode(0);
   
   TH1F *h_feast_t0__2 = new TH1F("h_feast_t0__2","h_feast_t0",50,-1,7);
   h_feast_t0__2->SetBinContent(7,1618);
   h_feast_t0__2->SetBinContent(8,2263);
   h_feast_t0__2->SetBinContent(9,1542);
   h_feast_t0__2->SetBinContent(10,1269);
   h_feast_t0__2->SetBinContent(11,1114);
   h_feast_t0__2->SetBinContent(12,1006);
   h_feast_t0__2->SetBinContent(13,874);
   h_feast_t0__2->SetBinContent(14,823);
   h_feast_t0__2->SetBinContent(15,766);
   h_feast_t0__2->SetBinContent(16,749);
   h_feast_t0__2->SetBinContent(17,651);
   h_feast_t0__2->SetBinContent(18,623);
   h_feast_t0__2->SetBinContent(19,650);
   h_feast_t0__2->SetBinContent(20,587);
   h_feast_t0__2->SetBinContent(21,537);
   h_feast_t0__2->SetBinContent(22,576);
   h_feast_t0__2->SetBinContent(23,484);
   h_feast_t0__2->SetBinContent(24,500);
   h_feast_t0__2->SetBinContent(25,518);
   h_feast_t0__2->SetBinContent(26,472);
   h_feast_t0__2->SetBinContent(27,503);
   h_feast_t0__2->SetBinContent(28,443);
   h_feast_t0__2->SetBinContent(29,443);
   h_feast_t0__2->SetBinContent(30,374);
   h_feast_t0__2->SetBinContent(31,268);
   h_feast_t0__2->SetBinContent(32,129);
   h_feast_t0__2->SetBinContent(33,78);
   h_feast_t0__2->SetBinContent(34,40);
   h_feast_t0__2->SetBinContent(35,16);
   h_feast_t0__2->SetBinContent(36,7);
   h_feast_t0__2->SetBinContent(37,8);
   h_feast_t0__2->SetBinContent(38,11);
   h_feast_t0__2->SetBinContent(40,9);
   h_feast_t0__2->SetBinContent(41,4);
   h_feast_t0__2->SetBinContent(42,4);
   h_feast_t0__2->SetBinContent(43,5);
   h_feast_t0__2->SetBinContent(44,2);
   h_feast_t0__2->SetBinContent(45,7);
   h_feast_t0__2->SetBinContent(46,4);
   h_feast_t0__2->SetBinContent(47,5);
   h_feast_t0__2->SetBinContent(48,6);
   h_feast_t0__2->SetBinContent(49,2);
   h_feast_t0__2->SetBinContent(50,3);
   h_feast_t0__2->SetBinContent(51,30);
   h_feast_t0__2->SetEntries(20023);
   h_feast_t0__2->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   h_feast_t0__2->SetLineColor(ci);
   h_feast_t0__2->GetXaxis()->SetLabelFont(42);
   h_feast_t0__2->GetXaxis()->SetLabelSize(0.035);
   h_feast_t0__2->GetXaxis()->SetTitleSize(0.035);
   h_feast_t0__2->GetXaxis()->SetTitleFont(42);
   h_feast_t0__2->GetYaxis()->SetLabelFont(42);
   h_feast_t0__2->GetYaxis()->SetLabelSize(0.035);
   h_feast_t0__2->GetYaxis()->SetTitleSize(0.035);
   h_feast_t0__2->GetYaxis()->SetTitleFont(42);
   h_feast_t0__2->GetZaxis()->SetLabelFont(42);
   h_feast_t0__2->GetZaxis()->SetLabelSize(0.035);
   h_feast_t0__2->GetZaxis()->SetTitleSize(0.035);
   h_feast_t0__2->GetZaxis()->SetTitleFont(42);
   h_feast_t0__2->Draw("");
   c_feast_t0->Modified();
   c_feast_t0->cd();
   c_feast_t0->SetSelected(c_feast_t0);
}
