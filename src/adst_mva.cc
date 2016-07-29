#include "adst_mva.h"
#include "workstation.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <stdlib.h>

using namespace std;

// Class that holds the root file structure -------------------------------------------------------------
Observables::Observables()
{
   xmax = -1;
   x0 = -1;
   x1 = -1;
   lambda = -1;
   fdenergy = -1;
   shfoot = -1;
}

Observables::~Observables()
{
}
// ------------------------------------------------------------------------------------------------------

// Analysis tools class constructor and destructor ------------------------------------------------------
//    Define the TTree names, set up and delete the canvas
AdstMva::AdstMva()
{
   gSystem->Load("libRecEventKG.so");
   gSystem->Load("libAnalysisKG.so");

/*   c1 = new TCanvas("c1AA","c1AA",1200,800);
   c1->SetGrid();
*/
   fRecEvent = new RecEvent();
   fDetGeo = new DetectorGeometry();

   outname = "output.root";

   shfootlimit = 0.1;

/*   xmaxlimit = 0;
   slantset = -1;
   shfootlimit = 0;
   distlimit = 0;
   pointcnt = 0;
   btemp = SUBPLOT;
   readevent = 0;

   for(int i = 0; i < ALLEYES; i++)
      eyevalid[i] = 0;*/
}

AdstMva::~AdstMva()
{
//   delete c1;
}
// ------------------------------------------------------------------------------------------------------

void AdstMva::RewriteObservables(int innr, Observables sig, Observables back)
{
   string stemp, stemp2;
   int itemp;
   
   cout << "# Entering function AdstMva::RewriteObservables()..." << endl;

   cout << "# New input file (" << inname[innr] << ") ---------------------------------" << endl;

   // Prepare signal and background trees
   stemp = "TreeS" + IntToStr(innr+1);
   stemp2 = "Signal tree from file " + inname[innr] + ".";
   
   sig_tree = new TTree(stemp.c_str(), stemp2.c_str());
   sig_tree->Branch("xmax", &(sig.xmax), "xmax/D");
   sig_tree->Branch("x0", &(sig.x0), "x0/D");
   sig_tree->Branch("x1", &(sig.x1), "x1/D");
   sig_tree->Branch("lambda", &(sig.lambda), "lambda/D");
   sig_tree->Branch("fdenergy", &(sig.fdenergy), "fdenergy/D");
   sig_tree->Branch("shfoot", &(sig.shfoot), "shfoot/D");

   back_tree->Branch("xmax", &(back.xmax), "xmax/D");
   back_tree->Branch("x0", &(back.x0), "x0/D");
   back_tree->Branch("x1", &(back.x1), "x1/D");
   back_tree->Branch("lambda", &(back.lambda), "lambda/D");
   back_tree->Branch("fdenergy", &(back.fdenergy), "fdenergy/D");
   back_tree->Branch("shfoot", &(back.shfoot), "shfoot/D");

   // Open and prepare the ADST files that we wish to read
   fFile = new RecEventFile(inname[innr].c_str(), RecEventFile::eRead);
   fFile->SetBuffers(&fRecEvent);
   fFile->ReadDetectorGeometry(*fDetGeo);

   cout << "Number of events: " << fFile->GetNEvents() << endl;

   // Go over all events in the ADST file and write them out to the output file
   for(int j = 0; j < fFile->GetNEvents(); j++)
   {
      cout << "# New event (" << j+1 << ") ---------------------------------" << endl;
      fFile->ReadEvent(j);

      if(acteyes.size() != 0)
         acteyes.erase(acteyes.begin(),acteyes.end());

      cout << "Number of eyes: " << fRecEvent->GetNEyes() << endl;
      if(fRecEvent->GetNEyes() == 0)
         cout << "Error! No reconstructed eyes for this event." << endl;
      else
      {
         vector<FDEvent> fdevt = fRecEvent->GetFDEvents();
         cout << "Size: " << fdevt.size() << endl;

         for(int i = 0; i < fdevt.size(); i++)
            acteyes.push_back(fdevt[i].GetFdRecShower());

         cout << "Number of active eyes: " << acteyes.size() << endl;

         itemp = GetEyeLongestTrack();
	 if( (itemp == -1) || (acteyes[itemp].GetEnergy() == 0) )
	    cout << "Error! The selected eye has no valid reconstructions." << endl;
	 else
	 {
            sig.xmax = acteyes[itemp].GetXmax();
            sig.x0 = acteyes[itemp].GetX0();
            sig.x1 = acteyes[itemp].GetX1();
            sig.lambda = acteyes[itemp].GetLambda();
            sig.fdenergy = acteyes[itemp].GetEnergy();

	    back.xmax = acteyes[itemp].GetXmax();
            back.x0 = acteyes[itemp].GetX0();
            back.x1 = acteyes[itemp].GetX1();
            back.lambda = acteyes[itemp].GetLambda();
            back.fdenergy = acteyes[itemp].GetEnergy();

	    if(GetShowerFoot(itemp, fdevt) == 0)
	    {
	       sig.shfoot = shfoot;
	       back.shfoot = shfoot;
	    }

            cout << "Values to save: " << endl
	         << "\t- Xmax = " << sig.xmax << endl
	         << "\t- X0 = " << sig.x0 << endl
	         << "\t- X1 = " << sig.x1 << endl
	         << "\t- Lambda = " << sig.lambda << endl
		 << "\t- FD Energy = " << sig.fdenergy << endl
		 << "\t- Shower foot = " << sig.shfoot << endl;

            sig_tree->Fill();
	    back_tree->Fill();
	 }
      }
   }

   sig_tree->Write();
}

// Check for the eye with the longest FD track
int AdstMva::GetEyeLongestTrack()
{
   cout << "# Entering function AdstMva::GetEyeLongestTrack()..." << endl;
   double longtrack = 0.0;
   int itemp;

   itemp = -1;

   for(int i = 0; i < acteyes.size(); i++)
   {
      cout << "Current track length: " << acteyes[i].GetTrackLength() << endl;
      if( acteyes[i].GetTrackLength() > longtrack )
      {
         itemp = i;
	 longtrack = acteyes[i].GetTrackLength();
      }
   }
   
   cout << "The longest track is " << longtrack << endl;

   return itemp;
}

int AdstMva::GetShowerFoot(int longestEye, vector<FDEvent> fdevt)
{
   double *x, *xerr;
   int *itemp;
   double *dtemp1, *dtemp2, *dtemp3;
//   double *outerr;

   x = new double;
   xerr = new double;
   itemp = new int;
   dtemp1 = new double;
   dtemp2 = new double;
   dtemp3 = new double;
//   outerr = new double[2];

   *x = 0;
   *xerr = 0;

   vector<double> slantDepth;
   vector<double> profiledEdX;
   vector<double> profiledEdXerr;

   slantDepth = acteyes[longestEye].GetDepth();
   profiledEdX = acteyes[longestEye].GetEnergyDeposit();
   profiledEdXerr = acteyes[longestEye].GetEnergyDepositError();

   vector<double> xfoot;
   vector<double> yfoot;
   vector<double> yerrfoot;

   if( xfoot.size() != 0 )
      xfoot.erase(xfoot.begin(),xfoot.end());
   if( yfoot.size() != 0 )
      yfoot.erase(yfoot.begin(),yfoot.end());
   if( yerrfoot.size() != 0 )
      yerrfoot.erase(yerrfoot.begin(),yerrfoot.end());

   for(int i = 0; i < slantDepth.size(); i++)
   {
      *x += profiledEdX[i];
      *xerr += profiledEdXerr[i];
   
      xfoot.push_back(slantDepth[i]);
      yfoot.push_back(*x);
      yerrfoot.push_back(*xerr);

//      cout << xfoot[i] << "\t" << yfoot[i] << "\t" << yerrfoot[i] << endl;
   }
   cout << endl;
   
   *itemp = 0;

   for(int i = 0; i < yfoot.size(); i++)
   {
      if( (yfoot[i] >= shfootlimit*(yfoot[yfoot.size()-1])) && (*itemp == 0) )
      {
         *itemp = 1;

	 // Find the x value of point with y value that is a fraction of the maximum, that lies on a line between two points
	 // y = k*x + a
	 //    k = (y2 - y1)/(x2 - x1)
	 //    a = y2 - (y2 - y1)/(x2 - x1)*x2
	 // x = ((x2 - x1)/(y2 - y1))*(y - y2) + x2
//	 cout << "Izr.1a: P1 = (" << xfoot[i-1] << "," << yfoot[i-1] << "), P2 = (" << xfoot[i] << "," << yfoot[i] << ")" << endl;

	 *dtemp1 = (xfoot[i] - xfoot[i-1])/(yfoot[i] - yfoot[i-1]); // 1/k = (x2 - x1)/(y2 - y1)
	 *dtemp2 = shfootlimit*(yfoot[yfoot.size()-1]) - yfoot[i]; // y - y2
         *x = (*dtemp1)*(*dtemp2) + xfoot[i]; // x = (1/k)*(y - y2) + x2

//	 cout << "Izr.1b: 1/k = " << *dtemp1 << ", max = " << yfoot[yfoot.size()-1] << ", y - y2 = " << *dtemp2 << ", x = " << *x << endl;
//	 cout << "Izr.2a: P1err = (" << xfoot[i-1] << "," << yerrfoot[i-1] << "), P2err = (" << xfoot[i] << "," << yerrfoot[i] << ")" << endl;

         *dtemp1 = (xfoot[i] - xfoot[i-1])/(yfoot[i]+yerrfoot[i] - (yfoot[i-1]+yerrfoot[i-1])); // 1/kerr = (x2 - x1)/(y2err - y1err)
	 *dtemp2 = (yfoot[i]+yerrfoot[i]) - (1/(*dtemp1))*(xfoot[i]); // aerr = y2err - (y2err - y1err)/(x2 - x1)*x2
	 *dtemp3 = (1/(*dtemp1))*(*x) + (*dtemp2); // yerr = kerr*x + aerr
	 *xerr = (*dtemp3) - (((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(*x) + (yfoot[i] - ((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(xfoot[i]))); // Dy = yerr - y

//	 cout << "Izr.2b: 1/kerr = " << *dtemp1 << ", a = " << *dtemp2 << ", yerr = " << *dtemp3 << ", Dy = " << *xerr << endl;

	 *dtemp1 = ((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(*x) + (yfoot[i] - ((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(xfoot[i])); // y = k*x + a

         for(int j = i; ; j++)
	 {
	    if( yfoot[j] >= (*dtemp1)+(*xerr) )
	    {
               *dtemp2 = xfoot[j];
//             cout << "Upper limit is = " << xfoot[j] << endl;
     	       break;
            }
         }
         for(int j = i; ; j--)
         {
            if(j == 0)
            {
	       *dtemp3 = xfoot[j];
//	       cout << "Lower limit is = " << xfoot[j] << endl;
     	       break;
            }

	    if( yfoot[j] <= (*dtemp1)-(*xerr) )
	    {
	       *dtemp3 = xfoot[j];
//	       cout << "Lower limit is = " << xfoot[j] << endl;
	       break;
	    }
	 }
//	 cout << "Izr.2c: x = " << *x << ", y = " << *dtemp1 << ", x+Dx = " << *dtemp2 << ", x-Dx = " << *dtemp3 << ", Dx+ = " << (*dtemp2)-(*x) << ", Dx- = " << (*x)-(*dtemp3) << endl;

         shfoot = *x;

//         *out += *x;//xfoot[i-1];
//         outerr[0] += (*x)-(*dtemp3);
//         outerr[1] += (*dtemp2)-(*x);
      }
   }

   delete x;
   delete xerr;
   delete itemp;
   delete dtemp1;
   delete dtemp2;
   delete dtemp3;
//   delete[] outerr;

   return 0;
}
