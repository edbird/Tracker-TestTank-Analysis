void c_feast_t0()
{
//=========Macro generated from canvas: c_feast_t0/c_feast_t0
//=========  (Tue Jan 16 17:31:08 2018) by ROOT version6.08/06
   TCanvas *c_feast_t0 = new TCanvas("c_feast_t0", "c_feast_t0",0,0,800,600);
   c_feast_t0->SetHighLightColor(2);
   c_feast_t0->Range(0,0,1,1);
   c_feast_t0->SetFillColor(0);
   c_feast_t0->SetBorderMode(0);
   c_feast_t0->SetBorderSize(2);
   c_feast_t0->SetFrameBorderMode(0);
   
   TH1F *h_feast_t0__2 = new TH1F("h_feast_t0__2","h_feast_t0",50,-1,7);
   h_feast_t0__2->SetBinContent(37,4866);
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
