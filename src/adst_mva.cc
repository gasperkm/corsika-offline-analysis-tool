#include "adst_mva.h"
#include "workstation.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <stdlib.h>

#define debug false

using namespace std;

// Class that holds the root file structure -------------------------------------------------------------
Observables::Observables()
{
   xmax = -1;
//   xmaxmu = -1;
   x0 = -1;
   lambda = -1;
   shfoot = -1;
   fdenergy = -1;
   nrmu = -1;
   shwsize = -1;
   ldfbeta = -1;
   curvature = -1;
   risetime = -1;
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
   genshw = new GenShower();
   sdrecshw = new SdRecShower();
#ifdef OFFLINENEW
   unishw = new UnivRecShower();
#endif

   outname = "tmva_output.root";

   shfootlimit = 0.1;

   goodrec = true;
   graphical = false;

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

void AdstMva::RewriteObservables(int innr, Observables sig, Observables back, TTree *back_tree)
{
   string stemp, stemp2;
   int itemp;
   bool singlerun = false;
   double dtemp;
   
   cout << "# Entering function AdstMva::RewriteObservables()..." << endl;

   if(inname.size() == 1)
   {
      innr = 0;
      singlerun = true;
   }

   cout << "# New input file (" << inname[innr] << ") ---------------------------------" << endl;

   // Prepare signal and background trees (all = complete set of events, back = only events that are not signal)
#ifdef OFFLINEOLD
   stemp = "TreeOldS" + IntToStr(innr+1);
#elif defined OFFLINENEW
   stemp = "TreeNewS" + IntToStr(innr+1);
#endif
   stemp2 = "Signal tree from file " + inname[innr] + ".";
   
   sig_tree = new TTree(stemp.c_str(), stemp2.c_str());
   sig_tree->Branch("xmax", &(sig.xmax), "xmax/F");
//   sig_tree->Branch("xmaxmu", &(sig.xmaxmu), "xmaxmu/F");
   sig_tree->Branch("x0", &(sig.x0), "x0/F");
   sig_tree->Branch("lambda", &(sig.lambda), "lambda/F");
   sig_tree->Branch("fdenergy", &(sig.fdenergy), "fdenergy/F");
   sig_tree->Branch("shfoot", &(sig.shfoot), "shfoot/F");
   sig_tree->Branch("shwsize", &(sig.shwsize), "shwsize/F");
   sig_tree->Branch("ldfbeta", &(sig.ldfbeta), "ldfbeta/F");
   sig_tree->Branch("curvature", &(sig.curvature), "curvature/F");
   sig_tree->Branch("nrmu", &(sig.nrmu), "nrmu/F");
   sig_tree->Branch("risetime", &(sig.risetime), "risetime/F");

   all_tree->Branch("xmax", &(back.xmax), "xmax/F");
//   all_tree->Branch("xmaxmu", &(back.xmaxmu), "xmaxmu/F");
   all_tree->Branch("x0", &(back.x0), "x0/F");
   all_tree->Branch("lambda", &(back.lambda), "lambda/F");
   all_tree->Branch("fdenergy", &(back.fdenergy), "fdenergy/F");
   all_tree->Branch("shfoot", &(back.shfoot), "shfoot/F");
   all_tree->Branch("shwsize", &(back.shwsize), "shwsize/F");
   all_tree->Branch("ldfbeta", &(back.ldfbeta), "ldfbeta/F");
   all_tree->Branch("curvature", &(back.curvature), "curvature/F");
   all_tree->Branch("nrmu", &(back.nrmu), "nrmu/F");
   all_tree->Branch("risetime", &(back.risetime), "risetime/F");

   for(int i = 0; i < inname.size(); i++)
   {
      back_tree[i].Branch("xmax", &(back.xmax), "xmax/F");
//      back_tree[i]->Branch("xmaxmu", &(back.xmaxmu), "xmaxmu/F");
      back_tree[i].Branch("x0", &(back.x0), "x0/F");
      back_tree[i].Branch("lambda", &(back.lambda), "lambda/F");
      back_tree[i].Branch("fdenergy", &(back.fdenergy), "fdenergy/F");
      back_tree[i].Branch("shfoot", &(back.shfoot), "shfoot/F");
      back_tree[i].Branch("shwsize", &(back.shwsize), "shwsize/F");
      back_tree[i].Branch("ldfbeta", &(back.ldfbeta), "ldfbeta/F");
      back_tree[i].Branch("curvature", &(back.curvature), "curvature/F");
      back_tree[i].Branch("nrmu", &(back.nrmu), "nrmu/F");
      back_tree[i].Branch("risetime", &(back.risetime), "risetime/F");
   }

   // Open and prepare the ADST files that we wish to read
   fFile = new RecEventFile(inname[innr].c_str(), RecEventFile::eRead);
   fFile->SetBuffers(&fRecEvent);
   fFile->ReadDetectorGeometry(*fDetGeo);

   cout << "Number of events: " << fFile->GetNEvents() << endl;

   int recfail[7];
   recfail[0] = 0;
   recfail[1] = 0;
   recfail[2] = 0;
   recfail[3] = 0;
   recfail[4] = 0;
   recfail[5] = 0;
   recfail[6] = 0;

   // Go over all events in the ADST file and write them out to the output file
   for(int j = 0; j < fFile->GetNEvents(); j++)
   {
      // goodrec variable determines if FD or SD reconstructions failed - not using results from those
      goodrec = true;

      cout << "# New event (" << j+1 << ") ---------------------------------" << endl;
      fFile->ReadEvent(j);

      if(acteyes.size() != 0)
         acteyes.erase(acteyes.begin(),acteyes.end());

      // Go over the FD eye events
      cout << "Number of eyes: " << fRecEvent->GetNEyes() << endl;
      if(fRecEvent->GetNEyes() == 0)
      {
         cout << "Error! No reconstructed eyes for this event." << endl;
	 if(goodrec) recfail[0]++;
	 goodrec = false;
      }
      else
      {
         vector<FDEvent> fdevt = fRecEvent->GetFDEvents();
         cout << "Size: " << fdevt.size() << endl;

         for(int i = 0; i < fdevt.size(); i++)
            acteyes.push_back(fdevt[i].GetFdRecShower());

         cout << "Number of active eyes: " << acteyes.size() << endl;

         for(int i = 0; i < fdevt.size(); i++)
	 {
	    if(!fdevt[i].IsHybridEvent())
	    {
	       cout << "Error! Eye " << i << " not a hybrid event." << endl;
	       goodrec = false;
	    }
	    else
	    {
	       cout << "Eye " << i << " is a hybrid event." << endl;
	       goodrec = true;
	       break;
	    }
	 }

	 if(!goodrec) recfail[1]++;

         itemp = GetEyeLongestTrack();
	 if( (itemp == -1) || (acteyes[itemp].GetEnergy() == 0) )
	 {
	    cout << "Error! The selected eye has no valid reconstructions." << endl;
	    if(goodrec) recfail[2]++;
	    goodrec = false;
	 }
	 else
	 {
            sig.xmax = acteyes[itemp].GetXmax();
            sig.x0 = acteyes[itemp].GetX0();
            sig.lambda = acteyes[itemp].GetLambda();
            sig.fdenergy = acteyes[itemp].GetEnergy();

	    back.xmax = acteyes[itemp].GetXmax();
            back.x0 = acteyes[itemp].GetX0();
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
	         << "\t- Lambda = " << sig.lambda << endl
		 << "\t- FD Energy = " << sig.fdenergy << endl
		 << "\t- Shower foot = " << sig.shfoot << endl;
	 }
      }

      // Go over the SD tank events
      *sdrecshw = fRecEvent->GetSDEvent().GetSdRecShower();

      if(!(fRecEvent->GetSDEvent().HasTriggeredStations()))
      {
         cout << "Error! No triggered stations in SD reconstruction." << endl;
	 if(goodrec) recfail[3]++;
         goodrec = false;
      }

      if(!(fRecEvent->GetSDEvent().HasStations()))
      {
         cout << "Error! No stations in SD reconstruction." << endl;
	 if(goodrec) recfail[4]++;
         goodrec = false;
      }

      if(!(fRecEvent->GetSDEvent().HasVEMTraces()))
      {
         cout << "Error! No VEM traces in SD tanks." << endl;
	 if(goodrec) recfail[5]++;
         goodrec = false;
      }

      sig.shwsize = sdrecshw->GetShowerSize();
      sig.ldfbeta = sdrecshw->GetBeta();
      sig.curvature = sdrecshw->GetCurvature();
      sig.risetime = sdrecshw->GetRiseTimeResults().GetRiseTime1000();
      back.shwsize = sdrecshw->GetShowerSize();
      back.ldfbeta = sdrecshw->GetBeta();
      back.curvature = sdrecshw->GetCurvature();
      back.risetime = sdrecshw->GetRiseTimeResults().GetRiseTime1000();
      
      cout << "\t- Shower size (replacement for S1000) = " << sig.shwsize << endl
           << "\t- LDF Beta = " << sig.ldfbeta << endl
           << "\t- Curvature R = " << sig.curvature << endl
           << "\t- Risetime at 1000m = " << sig.risetime << endl;

      if(GetRisetime(j, sig.risetime, &dtemp) == -1)
      {
         cout << "Error: Getting risetime from VEM signal not possible." << endl;
	 if(goodrec) recfail[6]++;
	 goodrec = false;
      }
      else
      {
         sig.risetime = dtemp;
         back.risetime = dtemp;
	 cout << "\t- Risetime (recalculated) = " << sig.risetime << endl;
      }

      if(sig.risetime < 0) // TODO: For some reason, some of the Risetime results are -1 -> check why and maybe try to calculate them from actual VEM traces
      {
         cout << "Error! Risetime not calculated. " << sdrecshw->GetRiseTimeResults().GetRiseTime1000() << endl;
	 if(goodrec) recfail[6]++;
	 goodrec = false;
      }

      // Go over the simulated events (Muon number at ground level) - only if we have simulations!
      *genshw = fRecEvent->GetGenShower();
      sig.nrmu = genshw->GetMuonNumber();
      back.nrmu = genshw->GetMuonNumber();
      cout << "\t- Nr. of muons = " << sig.nrmu << endl;

#ifdef OFFLINENEW
      // Go over simulation reconstruction (Muon number at ground level) - only if we have PAO data!
//      *unishw = fRecEvent->GetSDEvent().GetUnivRecShower();
//      sig.xmaxmu = genshw->GetXmaxMu();
//      back.xmaxmu = genshw->GetXmaxMu();
//      cout << "\t- Muon Xmax = " << sig.xmaxmu << endl;

      // In case we have PAO data, calculate the number of muons at ground level
      cout << "### Checking for muons from PAO data ###" << endl;
      vector<SdRecStation> stationVector = fRecEvent->GetSDEvent().GetStationVector();
      for(int i = 0; i < stationVector.size(); i++)
      {
         // Only use stations that are valid candidates
         if( (stationVector[i].IsCandidate()) )
         {
	    cout << "station ID = " << stationVector[i].GetId() << ", muon component = " << stationVector[i].GetMuonComponent() << ", muon signal = " << stationVector[i].GetMuonSignal() << ", electron signal = " << stationVector[i].GetElectronSignal() << ", photon signal = " << stationVector[i].GetPhotonSignal() << endl;
	 }
      }
#endif

      if(goodrec)
      {
         sig_tree->Fill();
         all_tree->Fill();

         for(int i = 0; i < inname.size(); i++)
         {
            if(i != innr)
               back_tree[i].Fill();
         }
      }
   }

   cout << recfail[0]+recfail[1]+recfail[2]+recfail[3]+recfail[4]+recfail[5]+recfail[6] << " events have been removed:" << endl;
   cout << " - No reconstructed FD eyes:      " << recfail[0] << " events" << endl;
   cout << " - Not a hybrid event:            " << recfail[1] << " events" << endl;
   cout << " - No valid FD reconstructions:   " << recfail[2] << " events" << endl;
   cout << " - No triggered SD stations:      " << recfail[3] << " events" << endl;
   cout << " - No reconstructed SD stations:  " << recfail[4] << " events" << endl;
   cout << " - No reconstructed VEM traces:   " << recfail[5] << " events" << endl;
   cout << " - No reconstructed risetime:     " << recfail[6] << " events" << endl;

   sig_tree->Write();
}

void AdstMva::PrepareOtherTrees(int nr, int sig)
{
   if(sig == 1)
   {
      Observables othersig;
      TTree *other_sig_tree;
      other_sig_tree = new TTree[nr];

      for(int i = 0; i < nr; i++)
      {
#ifdef OFFLINEOLD
         other_sig_tree[i].SetNameTitle(("TreeNewS" + IntToStr(i+1)).c_str(), "Signal tree from file new file");
#elif defined OFFLINENEW
         other_sig_tree[i].SetNameTitle(("TreeOldS" + IntToStr(i+1)).c_str(), "Signal tree from old file");
#endif

         other_sig_tree[i].Branch("xmax", &(othersig.xmax), "xmax/F");
         other_sig_tree[i].Branch("x0", &(othersig.x0), "x0/F");
         other_sig_tree[i].Branch("lambda", &(othersig.lambda), "lambda/F");
         other_sig_tree[i].Branch("fdenergy", &(othersig.fdenergy), "fdenergy/F");
         other_sig_tree[i].Branch("shfoot", &(othersig.shfoot), "shfoot/F");
         other_sig_tree[i].Branch("shwsize", &(othersig.shwsize), "shwsize/F");
         other_sig_tree[i].Branch("ldfbeta", &(othersig.ldfbeta), "ldfbeta/F");
         other_sig_tree[i].Branch("curvature", &(othersig.curvature), "curvature/F");
         other_sig_tree[i].Branch("nrmu", &(othersig.nrmu), "nrmu/F");
         other_sig_tree[i].Branch("risetime", &(othersig.risetime), "risetime/F");

//         other_sig_tree[i].Fill();

         other_sig_tree[i].Write();
      }
   }
   else if(sig == 0)
   {
      Observables otherback;
      TTree *other_back_tree;
      other_back_tree = new TTree[nr];

      for(int i = 0; i < nr; i++)
      {
#ifdef OFFLINEOLD
         other_back_tree[i].SetNameTitle(("TreeNewB" + IntToStr(i+1)).c_str(), "Background tree without events from new file");
#elif defined OFFLINENEW
         other_back_tree[i].SetNameTitle(("TreeOldB" + IntToStr(i+1)).c_str(), "Background tree without events from old file");
#endif

         other_back_tree[i].Branch("xmax", &(otherback.xmax), "xmax/F");
         other_back_tree[i].Branch("x0", &(otherback.x0), "x0/F");
         other_back_tree[i].Branch("lambda", &(otherback.lambda), "lambda/F");
         other_back_tree[i].Branch("fdenergy", &(otherback.fdenergy), "fdenergy/F");
         other_back_tree[i].Branch("shfoot", &(otherback.shfoot), "shfoot/F");
         other_back_tree[i].Branch("shwsize", &(otherback.shwsize), "shwsize/F");
         other_back_tree[i].Branch("ldfbeta", &(otherback.ldfbeta), "ldfbeta/F");
         other_back_tree[i].Branch("curvature", &(otherback.curvature), "curvature/F");
         other_back_tree[i].Branch("nrmu", &(otherback.nrmu), "nrmu/F");
         other_back_tree[i].Branch("risetime", &(otherback.risetime), "risetime/F");

//         other_back_tree[i].Fill();

         other_back_tree[i].Write();
      }
   }
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

   x = new double;
   xerr = new double;
   itemp = new int;
   dtemp1 = new double;
   dtemp2 = new double;
   dtemp3 = new double;

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
      }
   }

   delete x;
   delete xerr;
   delete itemp;
   delete dtemp1;
   delete dtemp2;
   delete dtemp3;

   return 0;
}

int AdstMva::GetRisetime(int event, double inRisetime, double *outRisetime)
{
   if(!(fRecEvent->GetSDEvent().HasVEMTraces()))
      return -1;

   cout << "# Entering function AdstMva::GetRisetime()..." << endl;
   double *x;
   double *y;
   double *maxval;
   int *nrpoints;
   double *dtemp;
   int *itemp;
   int start_bin, stop_bin;
   string stemp;

   double byrange[2];
   double bzrange[2];
   byrange[0] = 1.e+40;
   byrange[1] = -1.e+40;
   bzrange[0] = 1.e+40;
   bzrange[1] = -1.e+40;

   x = new double;
   y = new double;
   maxval = new double;
   nrpoints = new int;
   dtemp = new double[2];
   itemp = new int[2];

   vector<SdRecStation> stationVector = fRecEvent->GetSDEvent().GetStationVector();
   cout << "Number of triggered stations: " << stationVector.size() << ", Zenith angle = " << fRecEvent->GetSDEvent().GetSdRecShower().GetZenith() << ", Cos Zenith angle = " << fRecEvent->GetSDEvent().GetSdRecShower().GetCosZenith() << ", Energy = " << fRecEvent->GetSDEvent().GetSdRecShower().GetEnergy() << endl;

   vector<float> time;
   vector<float> vemtrace;

   vector<float> yvalue;

   vector<double> riseVect;
   vector<double> riseVectErr;
   vector<double> distVect;

   // Weighting function for determining the error bars (Risetime1000LL.xml)
   TFormula *fRTWeights = new TFormula("RiseTimeWeights", "(80.0+(5.071e-7+6.48e-4*y-3.051e-4*y*y)*x*x)/z-16.46*y+36.16");
   double limitTankDistance[2];
   limitTankDistance[0] = 0.;
   limitTankDistance[1] = 1800.;

   dtemp[0] = 0;
   itemp[0] = 0;
   *outRisetime;

   // Check all stations
   for(int i = 0; i < stationVector.size(); i++)
   {
      // Only use stations that are valid candidates
      if( (stationVector[i].IsCandidate()) /*&& (stationVector[i].GetSPDistance() < 1500.)*/ )
      {
	 start_bin = stationVector[i].GetSignalStartSlot() - 4;
	 stop_bin = stationVector[i].GetSignalEndSlot();

	 if( (start_bin >= stop_bin) || (start_bin < 0) || (start_bin > 5000) ) start_bin = 0;

//         cout << "Tank " << i << " is a candidate (" << start_bin << "," << stop_bin << ")." << endl;

	 dtemp[1] = 0;
	 itemp[1] = 0;

         // Check all PMTs
         for(int j = 1; j <= 3; j++)
	 {
            if(time.size() != 0)
	       time.erase(time.begin(), time.end());
            if(yvalue.size() != 0)
	       yvalue.erase(yvalue.begin(), yvalue.end());
	    *y = 0;
	    *maxval = -1.e40;

            vemtrace = stationVector[i].GetVEMTrace(j);
	    *nrpoints = vemtrace.size();
//	    cout << "PMT " << j << ": Number of points in the VEM trace: " << *nrpoints << " --------------------------------------------------------" << endl;

	    // Continue if there is a VEM trace
	    if( *nrpoints > 0 )
	    {
               itemp[0]++;
	       itemp[1]++;

	       // Prepare the time vector (each point is multiplied by 25 to get nanoseconds)
	       for(int k = 0; k < *nrpoints; k++)
	       {
                  if( (k >= start_bin) && (k <= stop_bin) )
		  {
	             time.push_back((float)k*25.);

		     *y += vemtrace[k];
		     if(*y > *maxval) *maxval = *y;
		  
		     yvalue.push_back(*y);
		  }
	       }

//	       cout << "Number of points in the signal slot: " << yvalue.size() << endl;

               for(int k = 0; k < yvalue.size(); k++)
	       {
//	          cout << time[k]/25. << "\t" << yvalue[k]/(*maxval) << endl;

		  if(yvalue[k]/(*maxval) <= 0.10)
		  {
		     byrange[0] = yvalue[k]/(*maxval);
		     byrange[1] = yvalue[k+1]/(*maxval);

		     *y = 0.1;
                     // Find the x value of point with y value = *y = 0.1, that lies on a line between two points
                     // y = k*x + a
                     //    k = (y2 - y1)/(x2 - x1)
                     //    a = y2 - (y2 - y1)/(x2 - x1)*x2
                     // x = ((x2 - x1)/(y2 - y1))*(y - y2) + x2
                     *x = ((time[k+1] - time[k])*((*y) - byrange[1]))/(byrange[1] - byrange[0]) + time[k+1];

                     byrange[0] = *x;
                     byrange[1] = *y;
		  }

		  if(yvalue[k]/(*maxval) <= 0.50)
		  {
                     bzrange[0] = yvalue[k]/(*maxval);
                     bzrange[1] = yvalue[k+1]/(*maxval);

                     *y = 0.5;
                     // Find the x value of point with y value = *y = 0.5, that lies on a line between two points
                     // y = k*x + a
                     //    k = (y2 - y1)/(x2 - x1)
                     //    a = y2 - (y2 - y1)/(x2 - x1)*x2
                     // x = ((x2 - x1)/(y2 - y1))*(y - y2) + x2
                     *x = ((time[k+1] - time[k])*((*y) - bzrange[1]))/(bzrange[1] - bzrange[0]) + time[k+1];

                     bzrange[0] = *x;
                     bzrange[1] = *y;
		  }
	       }

//	       cout << "Reconstructed risetime = " << inRisetime << ", calculated risetime (" << byrange[0]/25. << "," << bzrange[0]/25. << ") = " << bzrange[0] - byrange[0] << endl;
	       dtemp[0] += bzrange[0] - byrange[0];
	       dtemp[1] += bzrange[0] - byrange[0];
	    }
	 }

         dtemp[1] = dtemp[1]/itemp[1];
	 cout << "Tank " << stationVector[i].GetId() << ", " << stationVector[i].GetSPDistance() << " m: Calculated average risetime (for " << itemp[1] << " PMTs in the tank) = " << dtemp[1] << endl;

         // Asymmetry correction
         double eventThetaRec = fRecEvent->GetSDEvent().GetSdRecShower().GetZenith();
         double secZenith = 1/cos(eventThetaRec);
         const double alpha = 96.73 + secZenith*(-282.40 + secZenith*(241.80 - 62.61*secZenith));
         const double gamma = -0.0009572 + secZenith*(0.002068 + secZenith*(-0.001362 + 0.0002861*secZenith));
         const double g = alpha + gamma * stationVector[i].GetSPDistance()*stationVector[i].GetSPDistance();
         const double zeta = stationVector[i].GetAzimuthSP();

         *outRisetime = dtemp[1] - g*cos(zeta);
	 cout << "Asymmetry corrected risetime = " << *outRisetime << endl;
	 riseVect.push_back(*outRisetime);
	 distVect.push_back(stationVector[i].GetSPDistance());
	 riseVectErr.push_back(fRTWeights->Eval(stationVector[i].GetSPDistance(), secZenith, stationVector[i].GetTotalSignal()));
      }
   }

   int z = 0;
   for(int i = 0; i < riseVect.size(); i++)
      if( (distVect[i] >= limitTankDistance[0]) && (distVect[i] <= limitTankDistance[1]) )
         z++;

   TCanvas *c1;
   if(event < 20)
      c1 = new TCanvas("c1","c1",1200,900);

   TGraphErrors riseGraph(riseVect.size(), &distVect.front(), &riseVect.front(), 0, &riseVectErr.front());
   TF1 risetimeFit("RisetimeFit", "40+[0]*x+[1]*x*x", limitTankDistance[0], limitTankDistance[1]);
   risetimeFit.SetParLimits(0, 0, 10000);
   risetimeFit.SetParLimits(1, 0, 10000);
   RiseTimeFunction(fRecEvent->GetSDEvent().GetSdRecShower().GetZenith(), fRecEvent->GetSDEvent().GetSdRecShower().GetEnergy(), &risetimeFit);
   riseGraph.Fit(&risetimeFit, "Q", "", limitTankDistance[0], limitTankDistance[1]);
   *outRisetime = risetimeFit.Eval(1000.);

   cout << "Reconstructed risetime = " << inRisetime << ", Calculated average risetime (for " << itemp[0] << " PMTs in all tanks, " << z << " fitting points) = " << *outRisetime << endl;

   if(event < 20)
   {
      c1->SetGrid();
      c1->cd();
      riseGraph.SetMarkerSize(1);
      riseGraph.SetMarkerStyle(20);
      riseGraph.Draw("AP");
      risetimeFit.Draw("same");

      TLine *line = new TLine(limitTankDistance[0], inRisetime, 1000., inRisetime);
      line->SetLineWidth(2);
      line->SetLineStyle(7);
      line->SetLineColor(kBlue+2);
      line->Draw("same");

      line = new TLine(limitTankDistance[0], *outRisetime, 1000., *outRisetime);
      line->SetLineWidth(2);
      line->SetLineStyle(7);
      line->SetLineColor(kOrange+2);
      line->Draw("same");

      riseGraph.GetXaxis()->SetRange(0., 2300.);
      riseGraph.GetXaxis()->SetRangeUser(0., 2300.);
      riseGraph.GetXaxis()->SetLimits(0., 2300.);
      riseGraph.GetYaxis()->SetRange(0., 500.);
      riseGraph.GetYaxis()->SetRangeUser(0., 500.);
      riseGraph.GetYaxis()->SetLimits(0., 500.);
      riseGraph.SetTitle(";Distance from shower axis (m);Risetime of SD tanks (ns)");
      riseGraph.GetXaxis()->SetTitleOffset(1.2);
      riseGraph.GetYaxis()->SetTitleOffset(1.2);

      gPad->Update();

      if(event < 10)
         stemp = string(BASEDIR) + "/results/risetime_graph_0" + IntToStr(event) + ".pdf";
      else
         stemp = string(BASEDIR) + "/results/risetime_graph_" + IntToStr(event) + ".pdf";
      c1->SaveAs(stemp.c_str());
      delete c1;
   }

   delete x;
   delete y;
   delete nrpoints;
   delete maxval;
   delete[] dtemp;
   delete[] itemp;

   return 0;
}

void AdstMva::RiseTimeFunction(double zenith, double energy, TF1 *risetimeFit)
{
  //Karen Mora paremeterisation:
  const double secTheta= 1./cos(zenith); 
  const double a_par1= -0.141152; 
  const double a_par2=0.0141074;   
  const double a_par3=1.25107;  
  const double a_par4=-0.405333;   
  
  const double b_par1=0.000904323;   
  const double b_par2=6.4291e-06; 
  const double b_par3=-1.09992; 
  const double b_par4=0.30987;  
  
  double alpha=(a_par1+a_par2*log10(energy))*exp(-0.5*pow((secTheta-a_par3)/a_par4, 2));
  double beta=(b_par1+b_par2*log10(energy))*(1+b_par3* secTheta+b_par4*pow(secTheta, 2));
  
  risetimeFit->SetParameter(0,alpha);
  risetimeFit->SetParameter(1,beta);
}

/*void AdstMva::CreateMVAPlots(double cut, vector<string> obs)
{
   // Prepare colors for signal, background and MVA cut line
   static Int_t c_AllLine     = TColor::GetColor("#0000ee");
   static Int_t c_AllFill     = TColor::GetColor("#7d99d1");
   static Int_t c_SignalLine  = TColor::GetColor("#ff0000");
   static Int_t c_SignalFill  = TColor::GetColor("#ff0000");
   static Int_t c_MvaCut      = TColor::GetColor("#ffff66");

   // All additional things we need for plotting
   TLegend *legend;
   TNtuple *sigtuple;
   TNtuple *backtuple;
   TLine *line;
   TH1F *basehist;

   // Variables for setting range and fill style of the legend
   double xrangeset = 1.1;
   double yrangeset = 1.15;
   int legendFill = 1001;
   double xhistmax[2];
   double yhistmax[2];

   string obslist = "";

   gStyle->SetOptStat(0);

   TCanvas *c1 = new TCanvas("c1","",1200,900);
   c1->SetGrid();
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);

   // Prepare a semicolon separated list of observables to be used for plotting
   if(obs.size() <= 0)
   {
      cout << "AdstMva::CreateMVAPlots(): Error! Incorrect number of observables." << endl;
      return;
   }

   for(int i = 0; i < obs.size(); i++)
      obslist = obslist + obs[i] + ":";
   obslist = obslist + "MVA";

   // Plotting Signal events before MVA and Signal events + false background events after MVA
   sigtuple = new TNtuple("sig","signal",obslist.c_str());
   backtuple = new TNtuple("back","back",obslist.c_str());

   backtuple->ReadFile((string(BASEDIR) + "/root_mva/plots/gkm_simple_signal_start.txt").c_str());	// signal before MVA cut
   sigtuple->ReadFile((string(BASEDIR) + "/root_mva/plots/gkm_simple_signal.txt").c_str());		// signal + wrong back after MVA cut
   sigtuple->SetLineColor(c_SignalLine);
   sigtuple->SetLineWidth(2);
   sigtuple->SetFillColor(c_SignalFill);
   sigtuple->SetFillStyle(3554);
   backtuple->SetLineColor(c_AllLine);
   backtuple->SetLineWidth(2);
   backtuple->SetFillColor(c_AllFill);
   backtuple->SetFillStyle(1001);

   for(int i = 0; i <= obs.size(); i++)
   {
      basehist = new TH1F("t1","",100,-50.,50.);
      if(i == obs.size())
      {
         // Setup X range
         SetPlotRange(sigtuple, backtuple, basehist, xhistmax, xrangeset, "MVA", 'x');
         backtuple->Draw("MVA","","SAME");
         sigtuple->Draw("MVA","","SAME");
      }
      else
      {
         // Setup X range
         SetPlotRange(sigtuple, backtuple, basehist, xhistmax, xrangeset, obs[i], 'x');
         backtuple->Draw(obs[i].c_str(),"","SAME");
         sigtuple->Draw(obs[i].c_str(),"","SAME");
      }

      legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-.08, gPad->GetLeftMargin()+.30, 1-gPad->GetTopMargin());
      legend->SetFillStyle(legendFill);
      legend->SetFillColor(c_MvaCut);
      legend->AddEntry(backtuple,"Signal before MVA cut","f");
      legend->AddEntry(sigtuple,"Signal + false signal after MVA cut","f");
      legend->SetBorderSize(1);
      legend->SetMargin(0.3);
      if(i < obs.size())
         legend->Draw("same");

      gPad->Update();

      if(i < obs.size())
      {
         // Setup Y range
         SetPlotRange(sigtuple, backtuple, basehist, yhistmax, yrangeset, obs[i], 'y');
         SetupAxis(basehist, obs[i]);

         basehist->GetXaxis()->SetRange(xhistmax[0], xhistmax[1]);
         basehist->GetXaxis()->SetRangeUser(xhistmax[0], xhistmax[1]);
         basehist->GetXaxis()->SetLimits(xhistmax[0], xhistmax[1]);

         c1->SaveAs((string(BASEDIR) + "/results/mva_analysis_signal_" + obs[i] + ".pdf").c_str());
//         c1->SaveAs((string(BASEDIR) + "/results/mva_analysis_signal_" + obs[i] + ".C").c_str());
      }
      else
      {
         // Setup Y range
         SetPlotRange(sigtuple, backtuple, basehist, yhistmax, yrangeset, "MVA", 'y');

         line = new TLine(cut, yhistmax[0], cut, yhistmax[1]);
         line->SetLineWidth(2);
         line->SetLineStyle(7);
         line->SetLineColor(kOrange+2);
         line->Draw("same");

         legend->Draw("same");

         basehist->GetXaxis()->SetRange(xhistmax[0], xhistmax[1]);
         basehist->GetXaxis()->SetRangeUser(xhistmax[0], xhistmax[1]);
         basehist->GetXaxis()->SetLimits(xhistmax[0], xhistmax[1]);

	 SetupAxis(basehist, "MVA");

         gPad->Update();

         c1->SaveAs((string(BASEDIR) + "/results/mva_analysis_signal_MVA.pdf").c_str());
//         c1->SaveAs((string(BASEDIR) + "/results/mva_analysis_signal_MVA.C").c_str());
      }

      delete basehist;
   }

   delete sigtuple;
   delete backtuple;

   // Plotting All events before MVA and Signal events + false background events after MVA
   sigtuple = new TNtuple("sig","signal",obslist.c_str());
   backtuple = new TNtuple("all","all",obslist.c_str());

   backtuple->ReadFile((string(BASEDIR) + "/root_mva/plots/gkm_simple_all.txt").c_str());		// all events before MVA cut
   sigtuple->ReadFile((string(BASEDIR) + "/root_mva/plots/gkm_simple_signal.txt").c_str());		// signal + wrong back after MVA cut
   sigtuple->SetLineColor(c_SignalLine);
   sigtuple->SetLineWidth(2);
   sigtuple->SetFillColor(c_SignalFill);
   sigtuple->SetFillStyle(3554);
   backtuple->SetLineColor(c_AllLine);
   backtuple->SetLineWidth(2);
   backtuple->SetFillColor(c_AllFill);
   backtuple->SetFillStyle(1001);

   for(int i = 0; i <= obs.size(); i++)
   {
      basehist = new TH1F("t1","",100,-50.,50.);
      if(i == obs.size())
      {
         // Setup X range
         SetPlotRange(sigtuple, backtuple, basehist, xhistmax, xrangeset, "MVA", 'x');
         backtuple->Draw("MVA","","SAME");
         sigtuple->Draw("MVA","","SAME");
      }
      else
      {
         // Setup X range
         SetPlotRange(sigtuple, backtuple, basehist, xhistmax, xrangeset, obs[i], 'x');
         backtuple->Draw(obs[i].c_str(),"","SAME");
         sigtuple->Draw(obs[i].c_str(),"","SAME");
      }

      legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-.08, gPad->GetLeftMargin()+.30, 1-gPad->GetTopMargin());
      legend->SetFillStyle(legendFill);
      legend->SetFillColor(c_MvaCut);
      legend->AddEntry(backtuple,"All events before MVA cut","f");
      legend->AddEntry(sigtuple,"Signal + false signal after MVA cut","f");
      legend->SetBorderSize(1);
      legend->SetMargin(0.3);
      if(i < obs.size())
         legend->Draw("same");

      gPad->Update();

      if(i < obs.size())
      {
         // Setup Y range
         SetPlotRange(sigtuple, backtuple, basehist, yhistmax, yrangeset, obs[i], 'y');
         SetupAxis(basehist, obs[i]);

         basehist->GetXaxis()->SetRange(xhistmax[0], xhistmax[1]);
         basehist->GetXaxis()->SetRangeUser(xhistmax[0], xhistmax[1]);
         basehist->GetXaxis()->SetLimits(xhistmax[0], xhistmax[1]);

         c1->SaveAs((string(BASEDIR) + "/results/mva_analysis_all_" + obs[i] + ".pdf").c_str());
//         c1->SaveAs((string(BASEDIR) + "/results/mva_analysis_all_" + obs[i] + ".C").c_str());
      }
      else
      {
         // Setup Y range
         SetPlotRange(sigtuple, backtuple, basehist, yhistmax, yrangeset, "MVA", 'y');

         line = new TLine(cut, yhistmax[0], cut, yhistmax[1]);
         line->SetLineWidth(2);
         line->SetLineStyle(7);
         line->SetLineColor(kOrange+2);
         line->Draw("same");

         legend->Draw("same");

         basehist->GetXaxis()->SetRange(xhistmax[0], xhistmax[1]);
         basehist->GetXaxis()->SetRangeUser(xhistmax[0], xhistmax[1]);
         basehist->GetXaxis()->SetLimits(xhistmax[0], xhistmax[1]);

	 SetupAxis(basehist, "MVA");

         gPad->Update();

         c1->SaveAs((string(BASEDIR) + "/results/mva_analysis_all_MVA.pdf").c_str());
//         c1->SaveAs((string(BASEDIR) + "/results/mva_analysis_all_MVA.C").c_str());
      }

      delete basehist;
   }

   delete sigtuple;
   delete backtuple;

   // Plotting only starting signal and background before MVA
   sigtuple = new TNtuple("sig","signal",obslist.c_str());
   backtuple = new TNtuple("back","back",obslist.c_str());

   backtuple->ReadFile((string(BASEDIR) + "/root_mva/plots/gkm_simple_back_start.txt").c_str());		// back before MVA cut
   sigtuple->ReadFile((string(BASEDIR) + "/root_mva/plots/gkm_simple_signal_start.txt").c_str());		// signal before MVA cut
   sigtuple->SetLineColor(c_SignalLine);
   sigtuple->SetLineWidth(2);
   sigtuple->SetFillColor(c_SignalFill);
   sigtuple->SetFillStyle(3554);
   backtuple->SetLineColor(c_AllLine);
   backtuple->SetLineWidth(2);
   backtuple->SetFillColor(c_AllFill);
   backtuple->SetFillStyle(1001);

   for(int i = 0; i <= obs.size(); i++)
   {
      basehist = new TH1F("t1","",100,-50.,50.);
      if(i == obs.size())
      {
         // Setup X range
         SetPlotRange(sigtuple, backtuple, basehist, xhistmax, xrangeset, "MVA", 'x');
         backtuple->Draw("MVA","","SAME");
         sigtuple->Draw("MVA","","SAME");
      }
      else
      {
         // Setup X range
         SetPlotRange(sigtuple, backtuple, basehist, xhistmax, xrangeset, obs[i], 'x');
         backtuple->Draw(obs[i].c_str(),"","SAME");
         sigtuple->Draw(obs[i].c_str(),"","SAME");
      }

      legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-.08, gPad->GetLeftMargin()+.30, 1-gPad->GetTopMargin());
      legend->SetFillStyle(legendFill);
      legend->SetFillColor(c_MvaCut);
      legend->AddEntry(backtuple,"Background before MVA cut","f");
      legend->AddEntry(sigtuple,"Signal before MVA cut","f");
      legend->SetBorderSize(1);
      legend->SetMargin(0.3);
      if(i < obs.size())
         legend->Draw("same");

      gPad->Update();

      if(i < obs.size())
      {
         // Setup Y range
         SetPlotRange(sigtuple, backtuple, basehist, yhistmax, yrangeset, obs[i], 'y');
         SetupAxis(basehist, obs[i]);

         basehist->GetXaxis()->SetRange(xhistmax[0], xhistmax[1]);
         basehist->GetXaxis()->SetRangeUser(xhistmax[0], xhistmax[1]);
         basehist->GetXaxis()->SetLimits(xhistmax[0], xhistmax[1]);

         c1->SaveAs((string(BASEDIR) + "/results/mva_analysis_back_" + obs[i] + ".pdf").c_str());
//         c1->SaveAs((string(BASEDIR) + "/results/mva_analysis_back_" + obs[i] + ".C").c_str());
      }
      else
      {
         // Setup Y range
         SetPlotRange(sigtuple, backtuple, basehist, yhistmax, yrangeset, "MVA", 'y');

         line = new TLine(cut, yhistmax[0], cut, yhistmax[1]);
         line->SetLineWidth(2);
         line->SetLineStyle(7);
         line->SetLineColor(kOrange+2);
         line->Draw("same");

         legend->Draw("same");

         basehist->GetXaxis()->SetRange(xhistmax[0], xhistmax[1]);
         basehist->GetXaxis()->SetRangeUser(xhistmax[0], xhistmax[1]);
         basehist->GetXaxis()->SetLimits(xhistmax[0], xhistmax[1]);

	 SetupAxis(basehist, "MVA");

         gPad->Update();

         c1->SaveAs((string(BASEDIR) + "/results/mva_analysis_back_MVA.pdf").c_str());
//         c1->SaveAs((string(BASEDIR) + "/results/mva_analysis_back_MVA.C").c_str());
      }

      delete basehist;
   }

   delete sigtuple;
   delete backtuple;
}*/

void AdstMva::CreateMVAPlots(double cut, vector<string> obs, TTree *app, TMVA::Reader *reader, string mvamethod, float *obsvars, string signalName)
{
   if(obs.size() <= 0)
   {
      cout << "AdstMva::CreateMVAPlots(): Error! Incorrect number of observables." << endl;
      return;
   }

   // Prepare colors for signal, background and MVA cut line
   static Int_t c_SignalLine     = TColor::GetColor("#0000ee");
   static Int_t c_SignalFill     = TColor::GetColor("#7d99d1");
   static Int_t c_AllLine 	 = TColor::GetColor("#ff0000");
   static Int_t c_AllFill 	 = TColor::GetColor("#ff0000");
   static Int_t c_MvaCut     	 = TColor::GetColor("#ffff66");

   // All additional things we need for plotting
   TLegend *legend;
   TLine *line;
   TH1F *basesig[obs.size()+1];
   TH1F *baseback[obs.size()+1];

   string stemp;

   int legendFill = 1001;
   float yhistlimit[2];

   float *max;
   int *sigcount, *backcount;
   max = new float[obs.size()+1];
   sigcount = new int[obs.size()+1];
   backcount = new int[obs.size()+1];

//   string obslist = "";

   gStyle->SetOptStat(0);

   TCanvas *c1 = new TCanvas("c1","",1200,900);
   c1->SetGrid();
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);

   for(int i = 0; i <= obs.size(); i++)
   {
      if(i < obs.size())
      {
         SetupBinning(obs[i], yhistlimit);
         stemp = "basesig" + IntToStr(i);
         basesig[i] = new TH1F(stemp.c_str(), obs[i].c_str(), 100, yhistlimit[0], yhistlimit[1]);
         basesig[i]->SetBit(TH1::kCanRebin);
         stemp = "baseback" + IntToStr(i);
         baseback[i] = new TH1F(stemp.c_str(), obs[i].c_str(), 100, yhistlimit[0], yhistlimit[1]);
         baseback[i]->SetBit(TH1::kCanRebin);
      }
      else
      {
         SetupBinning("MVA", yhistlimit);
         stemp = "basesig" + IntToStr(i);
         basesig[i] = new TH1F(stemp.c_str(), "mva", 100, yhistlimit[0], yhistlimit[1]);
         basesig[i]->SetBit(TH1::kCanRebin);
         stemp = "baseback" + IntToStr(i);
         baseback[i] = new TH1F(stemp.c_str(), "mva", 100, yhistlimit[0], yhistlimit[1]);
         baseback[i]->SetBit(TH1::kCanRebin);
      }

      max[i] = 0.0;
      sigcount[i] = 0;
      backcount[i] = 0;
   }

   for(int ievt = 0; ievt < app->GetEntries(); ievt++)
   {
      app->GetEntry(ievt);

      for(int i = 0; i <= obs.size(); i++)
      {
         if(i < obs.size())
	 {
            if(reader->EvaluateMVA(mvamethod) >= cut)
	    {
	       basesig[i]->Fill(obsvars[i]);
	       sigcount[i]++;
	    }
            else
	    {
//	       cout << ievt << ", " << obs[i] << ": Background value = " << obsvars[i] << endl;
	       baseback[i]->Fill(obsvars[i]);
	       backcount[i]++;
	    }
	 }
	 else
	 {
            if(reader->EvaluateMVA(mvamethod) >= cut)
	    {
	       basesig[i]->Fill(reader->EvaluateMVA(mvamethod));
	       sigcount[i]++;
	    }
            else
	    {
//	       cout << ievt << ", MVA: Background value = " << reader->EvaluateMVA(mvamethod) << endl;
	       baseback[i]->Fill(reader->EvaluateMVA(mvamethod));
	       backcount[i]++;
	    }
	 }
      }
   }

   cout << "Signal vs. background:" << endl;
   for(int i = 0; i <= obs.size(); i++)
   {
      if(i < obs.size())
         cout << " - " << obs[i] << " = " << sigcount[i] << " vs. " << backcount[i] << endl;
      else
         cout << " - MVA = " << sigcount[i] << " vs. " << backcount[i] << endl;
   }
   cout << endl;

   string outname;

   for(int i = 0; i <= obs.size(); i++)
   {
      if(basesig[i]->GetMaximum() > max[i]) max[i] = basesig[i]->GetMaximum();
      if(baseback[i]->GetMaximum() > max[i]) max[i] = baseback[i]->GetMaximum();

      cout << "Maximum = " << max[i] << endl;

      basesig[i]->SetLineColor(c_SignalLine);
      basesig[i]->SetLineWidth(2);
      basesig[i]->SetFillColor(c_SignalFill);
      basesig[i]->SetFillStyle(1001);
      baseback[i]->SetLineColor(c_AllLine);
      baseback[i]->SetLineWidth(2);
      baseback[i]->SetFillColor(c_AllFill);
      baseback[i]->SetFillStyle(3554);

      basesig[i]->Draw();
      baseback[i]->Draw("same");

      basesig[i]->GetYaxis()->SetRangeUser(0.,max[i]*1.2);
      basesig[i]->SetMaximum(max[i]*1.2);

      if(i < obs.size())
         SetupAxis(basesig[i], obs[i]);
      else
      {
         SetupAxis(basesig[i], "MVA");
         line = new TLine(cut, 0., cut, basesig[i]->GetMaximum());
         line->SetLineWidth(2);
         line->SetLineStyle(7);
         line->SetLineColor(kOrange+2);
         line->Draw("same");
      }

      // Draw legend
      TLegend *legend= new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-.10, gPad->GetLeftMargin()+.25, 1-gPad->GetTopMargin());
      legend->SetFillStyle(legendFill);
      legend->SetFillColor(c_MvaCut);
      stemp = "MVA cut iron events (" + IntToStr(backcount[i]) + ")";
      legend->AddEntry(baseback[i],stemp.c_str(),"f");
      stemp = "MVA cut proton events (" + IntToStr(sigcount[i]) + ")";
      legend->AddEntry(basesig[i],stemp.c_str(),"f");
      legend->SetBorderSize(1);
      legend->SetMargin(0.3);
      legend->Draw("same");

      if(i < obs.size())
         stemp = "results/mva_analysis_" + signalName + "_" + obs[i] + ".pdf";
      else
         stemp = "results/mva_analysis_" + signalName + "_MVA.pdf";
      c1->SaveAs(stemp.c_str());
   }

   delete[] max;
   delete[] sigcount;
   delete[] backcount;
   for(int i = 0; i <= obs.size(); i++)
   {
      delete basesig[i];
      delete baseback[i];
   }
   delete c1;
}

void AdstMva::SetPlotRange(TNtuple *sig, TNtuple *back, TH1F *base, double *val, double rangeboost, string obs, char xory)
{
   double *dtemp;

   if(xory == 'x')
   {
      dtemp = new double[2];

      if(back->GetMinimum(obs.c_str()) < sig->GetMinimum(obs.c_str()))
         val[0] = back->GetMinimum(obs.c_str());
      else
         val[0] = sig->GetMinimum(obs.c_str());

      if(back->GetMaximum(obs.c_str()) > sig->GetMaximum(obs.c_str()))
         val[1] = back->GetMaximum(obs.c_str());
      else
         val[1] = sig->GetMaximum(obs.c_str());

      if(debug)
         cout << obs << ": start X min = " << val[0] << ", start X max = " << val[1] << endl;
      base->SetBit(TH1::kCanRebin);
      for(int i = 0; i <= base->GetNbinsX(); i++)
         base->SetBinContent(i, 0);
      base->SetLineColor(1);
      base->SetLineWidth(1);
      base->Draw();

      gPad->Update();

      // Set a small empty space on each side of the X axis
      if(val[0] < 0)
         dtemp[0] = val[0]*rangeboost;
      else if(val[0] == 0)
         dtemp[0] = val[1] - val[1]*rangeboost;
      else
         dtemp[0] = val[0]*(2.-rangeboost);

      if(val[1] < 0)
         dtemp[1] = val[1]*(2.-rangeboost);
      else if(val[1] == 0)
         dtemp[1] = val[0] - val[0]*rangeboost;
      else
         dtemp[1] = val[1]*rangeboost;

      if(debug)
         cout << obs << ": boost X min = " << dtemp[0] << " (diff = " << TMath::Abs(dtemp[0] - val[0]) << "), boost X max = " << dtemp[1] << " (diff = " << TMath::Abs(dtemp[1] - val[1]) << ")" << endl;

      val[0] = dtemp[0];
      val[1] = dtemp[1];

      if(debug)
         cout << obs << ": after X min = " << val[0] << ", after X max = " << val[1] << endl;

      // Set ranges for the histogram
      base->GetXaxis()->SetRange(val[0], val[1]);
      base->GetXaxis()->SetRangeUser(val[0], val[1]);
      base->GetXaxis()->SetLimits(val[0], val[1]);

      gPad->Update();

      delete[] dtemp;
   }
   else if(xory == 'y')
   {
      TH1F *hist1 = (TH1F*)back->GetHistogram();
      TH1F *hist2 = (TH1F*)sig->GetHistogram();

      // Get Y axis minimum
      val[0] = 0;
      // Get Y axis maximum
      if(hist1->GetMaximum() < hist2->GetMaximum())
         val[1] = hist2->GetMaximum();
      else
         val[1] = hist1->GetMaximum();

      if(debug)
         cout << obs << ": before Y min = " << val[0] << ", before Y max = " << val[1] << endl;

      // Set a small empty space on the top of the Y axis
      val[1] = val[1]*rangeboost;

      if(debug)
         cout << obs << ": after Y min = " << val[0] << ", after Y max = " << val[1] << endl;

      // Set ranges for the histogram
      base->GetYaxis()->SetRange(val[0], val[1]);
      base->GetYaxis()->SetRangeUser(val[0], val[1]);
      base->GetYaxis()->SetLimits(val[0], val[1]);

      gPad->Update();
   }
}

void AdstMva::SetupAxis(TH1F *hist, string obs)
{
   string obsdesc;
   // Determine which observable we have
   if(obs == "xmax")
      obsdesc = "X_{max} (g/cm^{2})";
   else if(obs == "x0")
      obsdesc = "GH parameter X_{0} (g/cm^{2})";
   else if(obs == "lambda")
      obsdesc = "GH parameter lambda (g/cm^{2})";
   else if(obs == "fdenergy")
      obsdesc = "FD reconstructed energy (eV)";
   else if(obs == "shfoot")
      obsdesc = "Shower foot (g/cm^{2})";
   else if(obs == "shwsize")
      obsdesc = "SD signal at 1000m from axis (VEM)";
   else if(obs == "ldfbeta")
      obsdesc = "LDF beta parameter";
   else if(obs == "curvature")
      obsdesc = "Curvature of shower plane (m^{-1})";
   else if(obs == "nrmu")
      obsdesc = "Number of muons at ground level";
   else if(obs == "risetime")
      obsdesc = "SD tank risetime (ns)";
   else if(obs == "MVA")
      obsdesc = "MVA observable";
   else
      return;

   hist->GetYaxis()->SetTitle("Number of events");
   hist->GetXaxis()->SetTitle(obsdesc.c_str());

   hist->GetXaxis()->SetTitleOffset(1.2);
   hist->GetXaxis()->CenterTitle(kTRUE);
   hist->GetXaxis()->SetLabelSize(0.028);
   hist->GetXaxis()->SetLabelOffset(0.015);
   hist->GetYaxis()->SetTitleOffset(1.3);
   hist->GetYaxis()->CenterTitle(kTRUE);
   hist->GetYaxis()->SetLabelSize(0.028);
   hist->GetYaxis()->SetLabelOffset(0.015);

   hist->SetTitle("");
}

void AdstMva::SetupBinning(std::string obs, float *limit)
{
   // Determine which observable we have
   if(obs == "xmax")
   {
      limit[0] = 400.;
      limit[1] = 1200.;
   }
   else if(obs == "x0")
   {
      limit[0] = -500.;
      limit[1] = 100.;
   }
   else if(obs == "lambda")
   {
      limit[0] = 40.;
      limit[1] = 80.;
   }
   else if(obs == "fdenergy")
   {
      limit[0] = 1.e+17;
      limit[1] = 1.e+20;
   }
   else if(obs == "shfoot")
   {
      limit[0] = 300.;
      limit[1] = 900.;
   }
   else if(obs == "shwsize")
   {
      limit[0] = 0.;
      limit[1] = 70.;
   }
   else if(obs == "ldfbeta")
   {
      limit[0] = -2.45;
      limit[1] = -2.25;
   }
   else if(obs == "curvature")
   {
      limit[0] = 1.e-4;
      limit[1] = 0.25e-3;
   }
   else if(obs == "nrmu")
   {
      limit[0] = 0.;
      limit[1] = 100.;
   }
   else if(obs == "risetime")
   {
      limit[0] = 0.;
      limit[1] = 500.;
   }
   else if(obs == "MVA")
   {
      limit[0] = -0.5;
      limit[1] = 1.5;
   }
   else
      return;
}
