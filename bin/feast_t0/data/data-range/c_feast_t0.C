void c_feast_t0()
{
//=========Macro generated from canvas: c_feast_t0/c_feast_t0
//=========  (Tue Jan 16 17:32:55 2018) by ROOT version6.08/06
   TCanvas *c_feast_t0 = new TCanvas("c_feast_t0", "c_feast_t0",0,0,800,600);
   c_feast_t0->SetHighLightColor(2);
   c_feast_t0->Range(0,0,1,1);
   c_feast_t0->SetFillColor(0);
   c_feast_t0->SetBorderMode(0);
   c_feast_t0->SetBorderSize(2);
   c_feast_t0->SetFrameBorderMode(0);
   
   TH1F *h_feast_t0__2 = new TH1F("h_feast_t0__2","h_feast_t0",50,4.76,4.86);
   h_feast_t0__2->SetBinContent(9,2);
   h_feast_t0__2->SetBinContent(10,3);
   h_feast_t0__2->SetBinContent(11,9);
   h_feast_t0__2->SetBinContent(12,42);
   h_feast_t0__2->SetBinContent(13,127);
   h_feast_t0__2->SetBinContent(14,169);
   h_feast_t0__2->SetBinContent(15,182);
   h_feast_t0__2->SetBinContent(16,200);
   h_feast_t0__2->SetBinContent(17,208);
   h_feast_t0__2->SetBinContent(18,229);
   h_feast_t0__2->SetBinContent(19,199);
   h_feast_t0__2->SetBinContent(20,203);
   h_feast_t0__2->SetBinContent(21,221);
   h_feast_t0__2->SetBinContent(22,195);
   h_feast_t0__2->SetBinContent(23,180);
   h_feast_t0__2->SetBinContent(24,168);
   h_feast_t0__2->SetBinContent(25,219);
   h_feast_t0__2->SetBinContent(26,194);
   h_feast_t0__2->SetBinContent(27,212);
   h_feast_t0__2->SetBinContent(28,203);
   h_feast_t0__2->SetBinContent(29,183);
   h_feast_t0__2->SetBinContent(30,205);
   h_feast_t0__2->SetBinContent(31,210);
   h_feast_t0__2->SetBinContent(32,185);
   h_feast_t0__2->SetBinContent(33,227);
   h_feast_t0__2->SetBinContent(34,211);
   h_feast_t0__2->SetBinContent(35,226);
   h_feast_t0__2->SetBinContent(36,155);
   h_feast_t0__2->SetBinContent(37,77);
   h_feast_t0__2->SetBinContent(38,20);
   h_feast_t0__2->SetBinContent(39,2);
   h_feast_t0__2->SetEntries(4866);
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
