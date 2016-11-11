void quick_MVA()
{
   static Int_t c_AllLine     = TColor::GetColor("#0000ee");
   static Int_t c_AllFill     = TColor::GetColor("#7d99d1");
   static Int_t c_SignalLine  = TColor::GetColor("#ff0000");
   static Int_t c_SignalFill  = TColor::GetColor("#ff0000");
   static Int_t c_Back      = TColor::GetColor("#ffff66");

   gStyle->SetOptStat(0);

   TFile *ofile = TFile::Open("./tmva_output_quick.root","RECREATE");

   TMVA::Factory *factory = new TMVA::Factory("TMVAClassification",ofile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

   factory->AddVariable("xmax", 'F');
   factory->AddVariable("shfoot", 'F');
   factory->AddVariable("shwsize", 'F');
   factory->AddVariable("risetime", 'F');

   // TreeS1 = proton, TreeS2 = iron, TreeS3 = data
//   TFile *ifile = TFile::Open("tmva_newrisetime_proton_vs_iron_vs_data_10e19eV.root");
   TFile *ifile = TFile::Open("tmva_newrisetime_proton_vs_iron_vs_data_10deg_10e19eV.root");
   TTree *signal = (TTree*)ifile->Get("TreeS1");
   TTree *background = (TTree*)ifile->Get("TreeS2");

   factory->AddSignalTree(signal, 1.0);
   factory->AddBackgroundTree(background, 1.0);

   factory->PrepareTrainingAndTestTree("", "", "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V");

   factory->BookMethod(TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator");

   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();

   ifile->Close();
   delete factory;
   ofile->Close();

   system("./tmvagui tmva_output_quick.root");

   TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");
   float obsvars[4];
   reader->AddVariable("xmax", &obsvars[0]);
   reader->AddVariable("shfoot", &obsvars[1]);
   reader->AddVariable("shwsize", &obsvars[2]);
   reader->AddVariable("risetime", &obsvars[3]);
   reader->BookMVA("MLPBNN method", "./weights/TMVAClassification_MLPBNN.weights.xml");

//   ifile = TFile::Open("tmva_newrisetime_proton_vs_iron_vs_data_10e19eV.root");
   ifile = TFile::Open("tmva_newrisetime_proton_vs_iron_vs_data_10deg_10e19eV.root");

   string evalTree;
   cout << "Select the tree to evaluate: ";
   cin >> evalTree;

   TTree *dataapp = (TTree*)ifile->Get(evalTree.c_str());
   dataapp->SetBranchAddress("xmax", &obsvars[0]);
   dataapp->SetBranchAddress("shfoot", &obsvars[1]);
   dataapp->SetBranchAddress("shwsize", &obsvars[2]);
   dataapp->SetBranchAddress("risetime", &obsvars[3]);

   double cutmva;
   cout << "Select the cut to be performed on the MVA variable: ";
   cin >> cutmva;

   TCanvas *c1 = new TCanvas("c1","",1200,900);
   c1->SetGrid();
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);
   TH1F *hsig[4];
   hsig[0] = new TH1F("hsig0","xmax",100,400.,1200.);
   hsig[1] = new TH1F("hsig1","shfoot",100,200.,900.);
   hsig[2] = new TH1F("hsig2","shwsize",100,0.,160.);
   hsig[3] = new TH1F("hsig3","risetime",100,0.,600.);
   TH1F *hback[4];
   hback[0] = new TH1F("hback0","xmax",100,400.,1200.);
   hback[1] = new TH1F("hback1","shfoot",100,200.,900.);
   hback[2] = new TH1F("hback2","shwsize",100,0.,160.);
   hback[3] = new TH1F("hback3","risetime",100,0.,600.);

   float max[] = {0.0,0.0,0.0,0.0};
   int sigcount[] = {0,0,0,0};
   int backcount[] = {0,0,0,0};

   for(int ievt=0; ievt < dataapp->GetEntries(); ievt++)
   {
      dataapp->GetEntry(ievt);

      for(int i = 0; i < 4; i++)
      {
         if(reader->EvaluateMVA("MLPBNN method") >= cutmva)
	 {
//	    cout << ievt << ": Signal value =     " << obsvars[i] << endl;
	    hsig[i]->Fill(obsvars[i]);
	    sigcount[i]++;
	 }
	 else
	 {
	    cout << ievt << ": Background value = " << obsvars[i] << endl;
	    hback[i]->Fill(obsvars[i]);
	    backcount[i]++;
	 }
/*//         cout << obsvars[i] << "\t";
         fprintf(fpall, "%lf\t", obsvars[i]);
         fprintf(fpsigstart, "%lf\t", obsvars[i]);

         if(reader->EvaluateMVA("MLPBNN") >= cutmva)
            fprintf(fpsig, "%lf\t", obsvars[i]);
         else
            fprintf(fpback, "%lf\t", obsvars[i]);*/
      }

/*//    cout << reader->EvaluateMVA("MLPBNN") << "\n";
      fprintf(fpall, "%lf\n", reader->EvaluateMVA("MLPBNN"));
      fprintf(fpsigstart, "%lf\n", reader->EvaluateMVA("MLPBNN"));

      if(reader->EvaluateMVA("MLPBNN") >= cutmva)
         fprintf(fpsig, "%lf\n", reader->EvaluateMVA("MLPBNN"));
      else
         fprintf(fpback, "%lf\n", reader->EvaluateMVA("MLPBNN"));

      if(reader->EvaluateMVA("MLPBNN") >= cutmva)
         sigval[0]++;
      else
         backval[0]++;*/
   } 

   cout << "Signal vs. background:" << endl;
   cout << " - xmax = " << sigcount[0] << " vs. " << backcount[0] << endl;
   cout << " - shfoot = " << sigcount[1] << " vs. " << backcount[1] << endl;
   cout << " - shwsize = " << sigcount[2] << " vs. " << backcount[2] << endl;
   cout << " - risetime = " << sigcount[3] << " vs. " << backcount[3] << endl;

   string outname;

   for(int i = 0; i < 4; i++)
   {
      if(hsig[i]->GetMaximum() > max[i]) max[i] = hsig[i]->GetMaximum();
      if(hback[i]->GetMaximum() > max[i]) max[i] = hback[i]->GetMaximum();

      cout << "Maximum = " << max[i] << endl;

      hsig[i]->SetLineColor(c_SignalLine);
      hsig[i]->SetLineWidth(2);
      hsig[i]->SetFillColor(c_SignalFill);
      hsig[i]->SetFillStyle(3554);
      hback[i]->SetLineColor(c_AllLine);
      hback[i]->SetLineWidth(2);
      hback[i]->SetFillColor(c_AllFill);
      hback[i]->SetFillStyle(1001);

      if(i == 0)
         hback[i]->SetTitle(";X_{max} (g/cm^{2});Number of events");
      else if(i == 1)
         hback[i]->SetTitle(";Shower foot (g/cm^{2});Number of events");
      else if(i == 2)
         hback[i]->SetTitle(";SD signal at 1000m from axis (VEM);Number of events");
      else if(i == 3)
         hback[i]->SetTitle(";SD tank risetime (ns);Number of events");

      hback[i]->Draw();
      hsig[i]->Draw("same");

      hback[i]->GetYaxis()->SetRangeUser(0.,max[i]*1.2);
      hback[i]->SetMaximum(max[i]*1.2);

      // Draw legend
      TLegend *legend= new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-.10, gPad->GetLeftMargin()+.20, 1-gPad->GetTopMargin());
      legend->SetFillStyle(1001);
      legend->SetFillColor(c_Back);
      legend->AddEntry(hback[i],"MVA cut iron events","f");
      legend->AddEntry(hsig[i],"MVA cut proton events","f");
      legend->SetBorderSize(1);
      legend->SetMargin(0.3);
      legend->Draw("same");

      if(i == 0)
         c1->SaveAs("results/tmva_quick_xmax.pdf");
      else if(i == 1)
         c1->SaveAs("results/tmva_quick_shfoot.pdf");
      else if(i == 2)
         c1->SaveAs("results/tmva_quick_shwsize.pdf");
      else if(i == 3)
         c1->SaveAs("results/tmva_quick_risetime.pdf");
   }

//   ifile->Close();
}
