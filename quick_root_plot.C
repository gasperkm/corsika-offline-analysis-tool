void quick_root_plot()
{
   static Int_t c_AllLine     = TColor::GetColor("#0000ee");
   static Int_t c_AllFill     = TColor::GetColor("#7d99d1");
   static Int_t c_SignalLine = TColor::GetColor("#ff0000");
   static Int_t c_SignalFill = TColor::GetColor("#ff0000");
   static Int_t c_Back      = TColor::GetColor("#ffff66");

   gStyle->SetOptStat(0);

   TCanvas *c1 = new TCanvas("c1","",2000,1000);
   c1->SetRightMargin(0.03);
   c1->SetLeftMargin(0.08);
   c1->SetTopMargin(0.03);   
   c1->SetGrid();   

   TH1F *h1 = new TH1F("h1","h1",120,0,600);
   TH1F *h2 = new TH1F("h1","h1",120,0,600);
   TFile *file;
   TTree *rtree;
   float risetime;
   int entries;
   string filenames[2];
   filenames[0] = "tmva_newrisetime_proton_vs_iron_vs_data_10e19eV.root";
   filenames[1] = "tmva_proton_vs_iron_vs_data_10e19eV.root";
   
   for(int i = 0; i < 2; i++)
   {
      cout << "# New file -----------------------------" << endl;
      file = TFile::Open(filenames[i].c_str(),"READ");
      file->GetObject("TreeS3", rtree);

      rtree->SetBranchAddress("risetime", &risetime);
      entries = rtree->GetEntries();

      for(int j = 0; j < entries; j++)
      {
         rtree->GetEntry(j);
	 cout << j << ": " << risetime << endl;
	 if(i == 0)
	    h1->Fill(risetime);
	 else if(i == 1)
	    h2->Fill(risetime);
      }
   }

   h1->SetLineColor(c_SignalLine);
   h1->SetLineWidth(2);
   h1->SetFillColor(c_SignalFill);
   h1->SetFillStyle(3554);
   h2->SetLineColor(c_AllLine);
   h2->SetLineWidth(2);
   h2->SetFillColor(c_AllFill);
   h2->SetFillStyle(1001);

   h2->SetTitle(";Risetime of SD tanks (ns);Number of events");
   h2->GetXaxis()->SetTitleOffset(1.3);
   h2->GetYaxis()->SetRange(0.,80.);
   h2->GetYaxis()->SetRangeUser(0.,80.);
   h2->GetYaxis()->SetTitleOffset(1.1);

   h2->Draw();
   h1->Draw("same");

   // Draw legend
   TLegend *legend= new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-.10, gPad->GetLeftMargin()+.20, 1-gPad->GetTopMargin());
   legend->SetFillStyle(1001);
   legend->SetFillColor(c_Back);
   legend->AddEntry(h2,"Offline risetime","f");
   legend->AddEntry(h1,"Custom risetime","f");
   legend->SetBorderSize(1);
   legend->SetMargin(0.3);
   legend->Draw("same");

   c1->SaveAs("results/risetime_recalculation.pdf");

/*   static Int_t c_AllLine     = TColor::GetColor("#0000ee");
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
   legend->Draw("same");*/
}
