void quick_root_plot()
{
   static Int_t c_AllLine     = TColor::GetColor("#0000ee");
   static Int_t c_AllFill     = TColor::GetColor("#7d99d1");
   static Int_t c_SignalLine = TColor::GetColor("#ff0000");
   static Int_t c_SignalFill = TColor::GetColor("#ff0000");

   gStyle->SetOptStat(0);

   TCanvas *c1 = new TCanvas("c1","xmax",900,600);
   TNtuple *sig = new TNtuple("sig","signal","a:b:c:d:x");
   TNtuple *all = new TNtuple("all","all","a:b:c:d:x");
   sig->ReadFile("root_mva/plots/gkm_simple_signal.txt");
   all->ReadFile("root_mva/plots/gkm_simple_all.txt");
   sig->SetLineColor(c_SignalLine);
   sig->SetLineWidth(2);
   sig->SetFillColor(c_SignalFill);
   sig->SetFillStyle(3554);
   all->SetLineColor(c_AllLine);
   all->SetLineWidth(2);
   all->SetFillColor(c_AllFill);
   all->SetFillStyle(1001);

   // First observable plot
   all->Draw("a");
   sig->Draw("a","","SAME");

   // Draw legend
   TLegend *legend= new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-.10, gPad->GetLeftMargin()+.20, 1-gPad->GetTopMargin());
   legend->SetFillStyle(1);
   legend->AddEntry(all,"Before MVA cut","f");
   legend->AddEntry(sig,"After MVA cut","f");
   legend->SetBorderSize(1);
   legend->SetMargin(0.3);
   legend->Draw("same");

   // Second observable plot
   TCanvas *c2 = new TCanvas("c2","shfoot",900,600);
   all->Draw("b");
   sig->Draw("b","","SAME");
   legend->Draw("same");

   // Third observable plot
   TCanvas *c3 = new TCanvas("c3","shwsize",900,600);
   all->Draw("c");
   sig->Draw("c","","SAME");
   legend->Draw("same");

   // Fourth observable plot
   TCanvas *c4 = new TCanvas("c4","risetime",900,600);
   all->Draw("d");
   sig->Draw("d","","SAME");
   legend->Draw("same");

   // MVA observable plot
   TCanvas *c5 = new TCanvas("c5","mva",900,600);
   all->Draw("x");
   sig->Draw("x","","SAME");
   legend->Draw("same");
}
