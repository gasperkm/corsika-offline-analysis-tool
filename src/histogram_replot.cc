// Compiling: g++ -o histogram_replot.cc histogram_replot.cc `root-config --cflags --glibs`

#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TApplication.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMarker.h"
#include "TBox.h"
#include "TTree.h"
#include "TFile.h"
#include "TColor.h"
#include "TStyle.h"

#define debug false

using namespace std;

string remove_ext(string inname)
{
   string stemp;
   for(int i = 0; i < inname.length(); i++)
   {
      if( (inname[i] == '.') && (i > (int)(inname.length()-6)) )
      {
         stemp = inname.erase(i,string::npos);
         break;
      }
   }

   if(debug)
     cout << "Outfile (remove_ext): " << stemp << endl;
   return stemp;
}

string replace_ext(string inname, string newext)
{
   string stemp = remove_ext(inname);
   stemp = stemp + newext;
  
   if(debug)
      cout << "Outfile (replace_ext): " << stemp << endl;
   return stemp;
}

int main(int argc, char **argv)
{
   string stemp;
   int rebinit = 1;
   double range[2] = {0,0};
   int nrbins;
   double val;
   ifstream infile;
   int crpdf = 0;
   int argtype = 0;
   vector<string> innames;
   int nrfiles = 0;
   double meanval[4];
   double sigval[4];
   char obsdir[2];
   double maxrange = 0;

   // Check input filename
   if(argc > 1)
   {
      // Normal input with multiple files
      if(strcmp("-f",argv[1]) == 0)
      {
	 argtype = 0;
         for(int i = 2; i < argc; i++)
	 {
            stemp = string(argv[i]);
	    innames.push_back(stemp);
	    if(debug)
	       cout << "Filename: " << innames[nrfiles] << endl;
	    nrfiles++;
	 }
      }
      // Input for statistical analysis
      else if(strcmp("-s",argv[1]) == 0)
      {
         // Only accept two input files
	 if( argc != 4 )
	 {
            cout << "Statistical analysis only possible on two input files." << endl;
            return 1;
	 }
	 argtype = 1;
         for(int i = 2; i < argc; i++)
	 {
            stemp = string(argv[i]);
	    innames.push_back(stemp);
	    if(debug)
	       cout << "Filename: " << innames[nrfiles] << endl;
	    nrfiles++;
	 }
      }
      else
      {
         cout << "Wrong or no options selected." << endl;
         return 1;
      }
   }
   else
   {
      cout << "No arguments were supplied." << endl;
      return 1;
   }

   TApplication theApp("App", &argc, argv);

   TCanvas *c1 = new TCanvas("c1","c1",1100,700);
   c1->SetGrid();

   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   TH1F *histf[nrfiles];
   TH1F *histfcdf[2];

   for(int i = 0; i < nrfiles; i++)
   {
      // Open the input file
      infile.open((innames[i]).c_str(), ifstream::in);

      // Count the number of data points in the file
      int nrlines = count(istreambuf_iterator<char>(infile),istreambuf_iterator<char>(), '\n');
      cout << "Number of all data points = " << nrlines << endl;

      infile.seekg(0);

      if(i == 0)
      {
         // Settings for bin number and sizes
         cout << "Manual bin set (yes = 1, no = 0): ";
         cin >> rebinit;
	 while( (rebinit != 0) && (rebinit != 1) )
	 {
            cout << "Value different than 1 or 0 inserted. Please reenter manual bin set (yes = 1, no = 0): ";
            cin >> rebinit;
	 }
         if(rebinit == 1)
         {
            cout << "Set the histogram range (separate with space): ";
            cin >> range[0] >> range[1];
	    while(range[0] >= range[1])
	    {
	       range[0] = 0; range[1] = 0;
               cout << "Minimum range should be larger than maximum range." << endl << "Set the histogram range (separate with space): ";
               cin >> range[0] >> range[1];
	    }
            cout << "Histogram range is (" << range[0] << "," << range[1] << ")" << endl;
            cout << "Set the number of bins (choose -1 for sqrt(N), -2 for log2(n)+1): ";
            cin >> nrbins;
            if(nrbins == -1)
               nrbins = sqrt(nrlines);
            else if(nrbins == -2)
               nrbins = log2(nrlines) + 1;
            cout << "Number of histogram bins is " << nrbins << endl;
         }
         else if(rebinit != 1)
         {
            nrbins = 50;
            range[0] = 1;
            range[1] = 1;
         }

         // Settings for pdf creation
         cout << "Create a normal histogram (0), PDF (1), or CDF (2): ";
         cin >> crpdf;
	 while( (crpdf != 0) && (crpdf != 1) && (crpdf != 2) )
	 {
            cout << "Value different than 2, 1 or 0 inserted. Please reenter the plot setting (normal = 0, PDF = 1, CDF = 2): ";
            cin >> crpdf;
	 }
      }

      cout << endl << "Opening file: " << innames[i] << endl;

      // Histogram variable
      histf[i] = new TH1F("th1","",nrbins,range[0],range[1]);
      if(rebinit != 1)
         histf[i]->SetBit(TH1::kCanRebin);

      // Read the data file and fill the histogram
      while(infile.good())
      {
         if(infile.eof()) break;

         infile >> val;
         if(debug) cout << val << endl;

         infile.ignore(1,' ');

         histf[i]->Fill(val);
      }

      // Rescale the histogram to produce a Probability Density Function
      if(crpdf == 1)
         histf[i]->Scale(1./nrlines);

      // Find the maximal value to determine the range for the double plot
      if(maxrange < histf[i]->GetMaximum())
      {
         if(debug)
            cout << "Maximum: " << maxrange << ", Current value: " << histf[i]->GetMaximum() << endl;
         maxrange = histf[i]->GetMaximum();
      }

      // Plot the histogram and save it as a root macro
      histf[i]->Draw();
      stemp = replace_ext(innames[i],".C");
      histf[i]->SaveAs(stemp.c_str());

      // Additional statistical analysis (only possible for two files)
      if(argtype == 1)
      {
	 meanval[i] = histf[i]->GetMean();
	 meanval[i+2] = histf[i]->GetMeanError();
	 sigval[i] = histf[i]->GetStdDev();
	 sigval[i+2] = histf[i]->GetStdDevError();

	 cout << "The mean value is: " << meanval[i] << " +- " << meanval[i+2] << endl;
	 cout << "The RMS is: " << sigval[i] << " +- " << sigval[i+2] << endl;

         // Set observation direction for the CDF: PDF that peaks earlier, will have CDF value 1 on the left (left observation direction), PDF that peaks later, will have CDF value 2 on the right (right observation direction)
	 if( i == nrfiles-1 )
	 {
	    if( meanval[0] > meanval[1] )
	    {
	       obsdir[0] = 'r';
	       obsdir[1] = 'l';

	       histfcdf[0] = (TH1F*)histf[0]->GetCumulative(kTRUE);
	       histfcdf[1] = (TH1F*)histf[1]->GetCumulative(kFALSE);
	    }
	    else if( meanval[1] > meanval[0] )
	    {
	       obsdir[0] = 'l';
	       obsdir[1] = 'r';

	       histfcdf[0] = (TH1F*)histf[0]->GetCumulative(kFALSE);
	       histfcdf[1] = (TH1F*)histf[1]->GetCumulative(kTRUE);
	    }

//            if(debug)
	       cout << "Observation directions: " << obsdir[0] << obsdir[1] << endl;

            if(crpdf == 2)
	    {
	       maxrange = 0;
               maxrange = histfcdf[0]->GetMaximum();
               histfcdf[0]->Scale(1./maxrange);
               maxrange = histfcdf[1]->GetMaximum();
               histfcdf[1]->Scale(1./maxrange);
            }
	 }
      }

      infile.close();
   }

   Int_t ci = 1756;
   TColor *color;

   // To finish, plot all created histograms on a single plot
   for(int i = 0; i < nrfiles; i++)
   {
      // Create color effect from red to blue
      if(nrfiles == 1)
         color = new TColor(ci, 1., 0., 0.);
      else
         color = new TColor(ci, 1.-((double)i/(double)(nrfiles-1.)), 0, (double)i/(double)(nrfiles-1.));

      if(crpdf < 2)
      {
         histf[i]->SetLineColor(ci);
         histf[0]->GetYaxis()->SetRange(0,maxrange+(0.1)*maxrange);
         histf[0]->GetYaxis()->SetRangeUser(0,maxrange+(0.1)*maxrange);
         histf[0]->GetYaxis()->SetLimits(0,maxrange+(0.1)*maxrange);

         if(i == 0)
            histf[i]->Draw();
         else
            histf[i]->Draw("same");
      }
      else if(crpdf == 2)
      {
         histfcdf[i]->SetLineColor(ci);
         histfcdf[0]->GetYaxis()->SetRange(0,1);
         histfcdf[0]->GetYaxis()->SetRangeUser(0,1.05);
         histfcdf[0]->GetYaxis()->SetLimits(0,1.05);

         if(i == 0)
            histfcdf[i]->Draw();
         else
            histfcdf[i]->Draw("same");
      }

      ci++;
   }

   // Run and exit the program
   theApp.Run();
   cout << "Exiting program..." << endl;

   delete[] *histf;
   delete color;

   return 0;
}
