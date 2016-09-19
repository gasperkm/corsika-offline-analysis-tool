#include "adstanalyse.h"
#include "workstation.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <stdlib.h>

using namespace std;

// Analysis tools class constructor and destructor ------------------------------------------------------
//    Define the TTree names, set up and delete the canvas
AdstAnalyseTool::AdstAnalyseTool()
{
   gSystem->Load("libRecEventKG.so");

   c1 = new TCanvas("c1AA","c1AA",1200,800);
   c1->SetGrid();

   fRecEvent = new RecEvent();
   fDetGeo = new DetectorGeometry();

   xmaxlimit = 0;
   slantset = -1;
   shfootlimit = 0;
   distlimit = 0;
   pointcnt = 0;
   btemp = SUBPLOT;
   readevent = 0;

   for(int i = 0; i < ALLEYES; i++)
      eyevalid[i] = 0;
}

AdstAnalyseTool::~AdstAnalyseTool()
{
//   outdata.close();
   delete c1;
}
// ------------------------------------------------------------------------------------------------------

// User given control to define type of analysis --------------------------------------------------------
//    Before using the Hello option, select if we want to use the program on a specific event or on all events
void AdstAnalyseTool::ReadOption()
{
   if(fFile->GetNEvents() == 1)
   {
      readopt = 1;
      readevent = 0;
      return;
   }

   cout << "# Entering function AdstAnalyseTool::ReadOption()..." << endl;

   cout << " The opened file has " << fFile->GetNEvents() << " events. Read all of them (0) or select a specific event (1)? ";
   cin >> readopt;

   while( (readopt != 0) && (readopt != 1) )
   {
      cout << " Invalid choice of read options. Please select again (0|1): ";
      cin >> readopt;
   }

   if(readopt == 1)
   {
      cout << " Which event do you want to read (from 0 to " << fFile->GetNEvents()-1 << ")? ";
      cin >> readevent;

      while( readevent >= fFile->GetNEvents() )
      {
         cout << "Invalid choice of event. Please select again (from 0 to " << fFile->GetNEvents()-1 << "): ";
	 cin >> readevent;
      }
   }
   else if(readopt == 0)
      readevent = -1;
}

//    Gives quick information on the program and gives possibilities for the analysis
void AdstAnalyseTool::Hello()
{
   cout << "# Entering function AdstAnalyseTool::Hello()..." << endl;
   string stemp;
   analysistype = 0;
   int itemp[2];

   // Hello display
   if(filesel == 0)
   {
      cout << endl;
      cout << "------------------------------------------------------------------------------------" << endl;
      cout << "--- Analysis tool for Offline (single ROOT analysis file mode)----------------------" << endl;
      cout << "------------------------------------------------------------------------------------" << endl;
      cout << " This is the analysis tool for data saved to a ROOT analysis file. It handles       " << endl
           << " plotting of all saved observables versus any other observable and plotting of tank " << endl
	   << " VEM signals, eye longitudinal profiles and SD reconstructed LDF functions.         " << endl << endl;
      cout << " The following analysis is possible:                                                " << endl
           << "    1 = Plot observables vs. other observables.                                     " << endl
	   << "    2 = Plot profiles (VEM signal, LDF profile, longitudinal profile).              " << endl
	   << "    3 = More to come                                                                " << endl;
      cout << "------------------------------------------------------------------------------------" << endl;
      cout << " Please select the analysis to perform: ";
      cin >> analysistype;
      cout << "------------------------------------------------------------------------------------" << endl;
   }
   
   if(filesel == 1)
   {
      cout << endl;
      cout << "------------------------------------------------------------------------------------" << endl;
      cout << "--- Analysis tool for Offline (multiple ROOT analysis file mode)--------------------" << endl;
      cout << "------------------------------------------------------------------------------------" << endl;
      cout << " This is the analysis tool for data saved to multiple ROOT analysis files. It       " << endl
           << " handles plotting of all saved observables versus any other observable and plotting " << endl
	   << " of tank VEM signals, eye longitudinal profiles and SD reconstructed LDF functions. " << endl
	   << " The program produces an average of observables over all tanks and eyes, enabling   " << endl
	   << " comparison of observables from multiple shower simulations.                        " << endl;
      cout << "------------------------------------------------------------------------------------" << endl;
      analysistype = 1;
   }

   // Plot observables
   if(analysistype == 1)
   {
      if(filesel == 0)
      {
         cout << " To plot observables, use the following notation: [obs1]:[obs2]                     " << endl; 
         cout << " Using only one will produce a histogram, using two will produce a scatter plot.    " << endl
              << " [obs1] will in this case be the quantity on the y-axis and [obs2] the quantity on  " << endl
              << " the x-axis.                                                                        " << endl;
         cout << " The possible observables are:                                                      " << endl
              << "    stationid      = Unique ID of the SD tank.                                      " << endl
              << "    totalvem       = Total VEM signal from the SD tank.                             " << endl
              << "    starttime.nsec = Nanosecond time when SD tank was triggered.                    " << endl
              << "    distance       = Distance from the shower core (SD reconstruction).             " << endl
              << "    muoncount      = Number of simulated muons in an SD tank.                       " << endl << endl
              << "    eyeid          = Unique ID of the FD eye.                                       " << endl
              << "    energy         = Reconstructed total shower energy from the FD eye.             " << endl
              << "    emenergy       = Reconstructed EM shower energy from the FD eye.                " << endl
              << "    xmax           = Reconstructed Xmax from the FD eye.                            " << endl;
      }
      else if(filesel == 1)
      {
         cout << " To plot observables, use the following notation: [obs1]:[obs2]:[obs3]              " << endl; 
         cout << " Using only one will produce a histogram, using two will produce a scatter plot,    " << endl
              << " using three will produce a color 2D histogram. [obs1] will in this case be the     " << endl
	      << " quantity on the x-axis, [obs2] the quantity on the y-axis and [obs3] the quantity  " << endl
              << " on the the color axis (z-axis).                                                    " << endl;
         cout << " The possible observables are:                                                      " << endl
              << "    actstation     = Number of activated SD tanks for selected reconstruction.       (SD)" << endl
              << "    totalvemmax    = Total VEM signal from SD tank with maximal VEM signal.          (SD)" << endl
              << "    totalvemaver   = Total VEM signal averaged over all active SD tanks.             (SD average)" << endl
              << "    totalvemsum    = Total VEM signal sum from all active SD tanks.                  (SD)" << endl
	      << "    risetime       = Risetime averaged over all active SD tanks.                     (SD average)" << endl
	      << "    risetimeearly  = Risetime averaged over all active SD tanks (start time < mean). (SD)" << endl
	      << "    risetimelate   = Risetime averaged over all active SD tanks (start time > mean). (SD)" << endl
	      << "    sdenergy       = Reconstructed total energy from active SD tanks.                (SD)" << endl
              << "    S1000          = Reconstructed S1000 value from active SD tanks (LDF at 1000m).  (SD)" << endl
              << "    curvature      = Reconstructed curvature value from active SD tanks.             (SD)" << endl
              << "    distaver       = Average distance from the shower core of all active SD tanks.   (SD average)" << endl
              << "    muoncount      = Total number of simulated muons from all active SD tanks.       (SD)" << endl
              << "    muoncountaver  = Number of simulated muons averaged over all active SD tanks.    (SD average)" << endl << endl
              << "    acteyes        = Number of activated FD eyes for selected reconstruction.        (FD)" << endl
              << "    fdenergy       = Reconstructed total energy averaged over active FD eyes.        (FD average)" << endl
              << "    fdcalenergy    = Reconstructed calorimetric energy averaged over active FD eyes. (FD average)" << endl
              << "    xmax           = Reconstructed Xmax averaged over active FD eyes.                (FD average)" << endl
              << "    xmaxqual       = Reconstructed Xmax of FD eye with smallest error on Xmax.       (FD)" << endl
	      << "    shfoot         = Depth when shower reaches N\% of cumulative long. profile.       (FD average)" << endl;
      }
      cout << "------------------------------------------------------------------------------------" << endl;
      cout << " Please enter the plot directive: ";
      cin >> plotdirective;
      cout << "------------------------------------------------------------------------------------" << endl;
   }
   // Plot profiles
   else if(analysistype == 2)
   {
      if(filesel == 0)
      {
         cout << " To plot profiles, use one the following:                                           " << endl 
              << "    ldf           = Plot the LDF profile (SD reconstruction).                       " << endl
              << "    vemsig:[ID]   = Plot the VEM signal from the selected tank with ID = [ID].      " << endl
              << "    longprof:[ID] = Plot the longitudinal profile from an eye with ID = [ID].       " << endl << endl
              << " When using [ID]=all, profiles for all tanks/eyes will be created and saved.        " << endl;
         cout << "------------------------------------------------------------------------------------" << endl;
         stemp = GetActiveTanks();
         cout << " SD tanks available: " << endl << "    " << stemp << endl;
         stemp = GetActiveEyeList();
         cout << " FD eyes available: " << endl << "    " << stemp << endl;
         if( (nracttanks > 0) || (nracteyes > 0) )
         {
            cout << "------------------------------------------------------------------------------------" << endl;
            cout << " Please enter the plot directive: ";
            cin >> plotdirective;
            cout << "------------------------------------------------------------------------------------" << endl;

            if( btemp && (plotdirective.find("longprof") != string::npos) )
            {
               cout << " Enter fraction of particles in longitudinal profile that still count for shfoot: ";
               cin >> shfootlimit;
               cout << "------------------------------------------------------------------------------------" << endl;
            }
         }
         else
            plotdirective = "-1";
      }
      else if(filesel == 1)
      {
         // TODO
         plotdirective = "-1";
      }
   }
   // More to come
   else if(analysistype == 3)
   {
      cout << " More to come..."  << endl; 
      cout << "------------------------------------------------------------------------------------" << endl;
      cout << " Please enter the plot directive: ";
      cin >> plotdirective;
      cout << "------------------------------------------------------------------------------------" << endl;
      plotdirective = "-1";
   }
   // Default plot directive (will stop any further analysis)
   else
   {
      cout << " Analysis type not supported." << endl;
      cout << "------------------------------------------------------------------------------------" << endl;
      plotdirective = "-1";
   }
}

// Initialise X, Y and Z plot ranges for later use, saves/returns input filename ------------------------
void AdstAnalyseTool::InitRanges()
{
   cout << "# Entering function AdstAnalyseTool::InitRanges()..." << endl;
   xrange[0] = 1.e+40;
   xrange[1] = -1.e+40;
   yrange[0] = 1.e+40;
   yrange[1] = -1.e+40;
   zrange[0] = 1.e+40;
   zrange[1] = -1.e+40;
}

void AdstAnalyseTool::SetAnalysisFilename(char *file)
{
   cout << "# Entering function AdstAnalyseTool::SetAnalysisFilename()..." << endl;
   adstanalysisfile = string(file);
}

string AdstAnalyseTool::GetAnalysisFilename()
{
   cout << "# Entering function AdstAnalyseTool::GetAnalysisFilename()..." << endl;
   return adstanalysisfile;
}
// ------------------------------------------------------------------------------------------------------

// Retrieve and save the tanks that are triggered in the reconstruction ---------------------------------
string AdstAnalyseTool::GetActiveTanks()
{
   cout << "# Entering function AdstAnalyseTool::GetActiveTanks()..." << endl;
   string stemp;

   fFile->ReadEvent(readevent);

   SDEvent &sdevt = fRecEvent->GetSDEvent();
   nracttanks = sdevt.GetNumberOfCandidates();

   if( acttanks.size() != 0 )
      acttanks.erase(acttanks.begin(),acttanks.end());

   acttanks = sdevt.GetStationVector();

   if(nracttanks == 0)
      return "    No reconstructed SD tanks.";
   else if(nracttanks == 1)
   {
      if(acttanks[0].GetId() < 90000)
         stemp = IntToStr(acttanks[0].GetId());
      else
         return "    No reconstructed SD tanks.";
   }
   else
   {
      for(int i = 0; i < (int)acttanks.size(); i++)
      {
         if(acttanks[i].GetId() < 90000)
	 {
            if(i == 0)
               stemp = IntToStr(acttanks[i].GetId());
            else
               stemp = stemp + ", " + IntToStr(acttanks[i].GetId());
	 }
      }
   }

   return stemp;
}

// ------------------------------------------------------------------------------------------------------

// Retrieve and save the FD eyes that are triggered in the reconstruction, produce a list ---------------
void AdstAnalyseTool::GetActiveEyes()
{
   cout << "# Entering function AdstAnalyseTool::GetActiveEyes()..." << endl;
   string stemp;
   double tempmax;
   int iMax;

   fFile->ReadEvent(readevent);

   nracteyes = fRecEvent->GetNEyes();
//cout << " Number of eyes = " << nracteyes << endl;
   vector<FDEvent> fdevt = fRecEvent->GetFDEvents();
//cout << " fdevt.size() = " << fdevt.size() << endl;

   if( acteyes.size() != 0 )
      acteyes.erase(acteyes.begin(),acteyes.end());

   for(int i = 0; i < ALLEYES; i++)
      eyevalid[i] = 0;

   for(int i = 0; i < fdevt.size(); i++)
   {
      acteyes.push_back(fdevt[i].GetFdRecShower());
//cout << " Eye " << fdevt[i].GetEyeId() << ": Energy = " << acteyes[i].GetEnergy() << " eV, Track = " << acteyes[i].GetTrackLength() << " g/cm^2" << endl;

      if( (acteyes[i].GetEnergy() == 0) || (acteyes[i].GetTrackLength() == 0) )
         eyevalid[fdevt[i].GetEyeId()-1] = -1;
      else
      {
         eyevalid[fdevt[i].GetEyeId()-1] = 1;
//         cout << "Eye " << i << " has ID = " << fdevt[i].GetEyeId() << ", energy = " << acteyes[i].GetEnergy() << " eV, track length = " << acteyes[i].GetTrackLength() << " g/cm2, track min = " << acteyes[i].GetXTrackMin() << ", track max = " << acteyes[i].GetXTrackMax() << endl;
      }

      if(i == 0)
         tempmax = acteyes[i].GetTrackLength();
      else
      {
	 if(acteyes[i].GetTrackLength() > tempmax)
	    iMax = fdevt[i].GetEyeId();
         tempmax = TMath::Max(tempmax, acteyes[i].GetTrackLength());
      }
   }

   cout << "Maximum track length is = " << tempmax << " g/cm2 (from eye " << iMax << ")" << endl;

//   cout << "Bool set: ";
//   for(int i = 0; i < ALLEYES; i++)
//      cout << eyevalid[i] << ", ";
//   cout << endl;
}

// Gives a list of all active eyes
string AdstAnalyseTool::GetActiveEyeList()
{
   cout << "# Entering function AdstAnalyseTool::GetActiveEyeList()..." << endl;
   int eyecnt = 0;
   string stemp;

   for(int i = 0; i < ALLEYES; i++)
   {
      if(eyevalid[i] == 1)
      {
         if(eyecnt == 0)
            stemp = IntToStr(i+1);
	 else
	    stemp = stemp + ", " + IntToStr(i+1);
         
         eyecnt++;
      }
   }

   if(eyecnt == 0)
      return "    No reconstructed FD eyes.";
   else
      return stemp;
}
// ------------------------------------------------------------------------------------------------------

// Start of the actual analysis program -----------------------------------------------------------------
//    infilenr = currently used input file number
int AdstAnalyseTool::RunDirective(int infilenr)
{
   cout << "# Entering function AdstAnalyseTool::RunDirective()..." << endl;
   string stemp;
   string dirtemp[3];
   int itemp[3] = {-1,-1,-1};
   double dtemp;
   int nrpoints;

   // Stop analysis if selection of analysis was not correct
   if(!plotdirective.find("-1"))
      return -1;

   nrdirs = sizeof(directivetype);

   // Find the : in the plot directive and split the two sides into dirtemp[0] and dirtemp[1] (check both if they are valid options)
   itemp[0] = plotdirective.find(':');
   if(plotdirective.find(':') != string::npos)
   {
      dirtemp[0] = plotdirective.substr(0,itemp[0]);
      itemp[1] = plotdirective.find(':', itemp[0]+1);
      if( plotdirective.find(':', itemp[0]+1) != string::npos )
      {
         dirtemp[1] = plotdirective.substr(itemp[0]+1,itemp[1] - (itemp[0]+1));
	 dirtemp[2] = plotdirective.substr(itemp[1]+1,plotdirective.length() - (itemp[1]+1));
      }
      else
      {
         dirtemp[1] = plotdirective.substr(itemp[0]+1,plotdirective.length() - (itemp[0]+1));
	 dirtemp[2] = "empty";
      }
    
      for(int i = 0; i < nrdirs; i++)
      {
         if(dirtemp[0].compare(directive[i]) == 0)
            itemp[0] = i;
         if(dirtemp[1].compare(directive[i]) == 0)
            itemp[1] = i;
         if(dirtemp[2].compare(directive[i]) == 0)
            itemp[2] = i;
      }
   }
   else
   {
      dirtemp[0] = plotdirective;
      for(int i = 0; i < nrdirs; i++)
         if(dirtemp[0].compare(directive[i]) == 0)
            itemp[0] = i;
   }

   // Check to see if the selected directive is available for single and multiple file mode
   nrdirs = 0;
   // Single file mode
   if(filesel == 0)
   {
      for(int i = 0; i < 3; i++)
      {
         if(itemp[i] != -1)
            nrdirs++;
         if( directiveavail[itemp[i]] == 1 )
	    itemp[i] = -1;
      }

      if((itemp[1] == -1) && (itemp[2] == -1))
      {
         if( directive[itemp[0]] == "xmax" )
         {
            cout << " Enter maximal Xmax error to be used (set to 0, to select all showers): ";
            cin >> xmaxlimit;
            cout << "------------------------------------------------------------------------------------" << endl;
            cout << " Use vertical (0) or slant depth (1): ";
            cin >> slantset;
            cout << "------------------------------------------------------------------------------------" << endl;
         }
      }
      else if((itemp[1] != -1) && (itemp[2] == -1))
      {
         if( (directive[itemp[0]] == "xmax") || (directive[itemp[1]] == "xmax") )
         {
            cout << " Enter maximal Xmax error to be used (set to 0, to select all showers): ";
            cin >> xmaxlimit;
            cout << "------------------------------------------------------------------------------------" << endl;
            cout << " Use vertical (0) or slant depth (1): ";
            cin >> slantset;
            cout << "------------------------------------------------------------------------------------" << endl;
         }
      }
   }
   // Multiple file mode
   else if(filesel == 1)
   {
      for(int i = 0; i < 3; i++)
      {
         if(itemp[i] != -1)
            nrdirs++;
         if( directiveavail[itemp[i]] == 0 )
	    itemp[i] = -1;
      }

      // Only valid for the multiple file mode
      if((itemp[1] == -1) && (itemp[2] == -1) && (infilenr == 0))
      {
         if( (directiveaffil[itemp[0]] == "showeyedata") || (directiveaffil[itemp[0]] == "eyelong") )
         {
            cout << " Enter maximal Xmax error to be used (set to 0, to select all showers): ";
            cin >> xmaxlimit;
            cout << "------------------------------------------------------------------------------------" << endl;
            cout << " Use vertical (0) or slant depth (1): ";
            cin >> slantset;
            cout << "------------------------------------------------------------------------------------" << endl;
         }
         if( directive[itemp[0]] == "shfoot" )
         {
            cout << " Enter fraction of particles in longitudinal profile that still count for shfoot: ";
            cin >> shfootlimit;
            cout << "------------------------------------------------------------------------------------" << endl;
         }
	 if( directiveaffil[itemp[0]] == "tankvem" )
	 {
            cout << " Enter maximal distance of SD tank from shower core (set to 0, to select all tanks): ";
            cin >> distlimit;
            cout << "------------------------------------------------------------------------------------" << endl;
	 }
      }
      else if((itemp[1] != -1) && (itemp[2] == -1) && (infilenr == 0))
      {
         if( (directiveaffil[itemp[0]] == "showeyedata") || (directiveaffil[itemp[0]] == "eyelong") || (directiveaffil[itemp[1]] == "showeyedata") || (directiveaffil[itemp[1]] == "eyelong") )
         {
            cout << " Enter maximal Xmax error to be used (set to 0, to select all showers): ";
            cin >> xmaxlimit;
            cout << "------------------------------------------------------------------------------------" << endl;
            cout << " Use vertical (0) or slant depth (1): ";
            cin >> slantset;
            cout << "------------------------------------------------------------------------------------" << endl;
         }
         if( (directive[itemp[0]] == "shfoot") || (directive[itemp[1]] == "shfoot") )
         {
            cout << " Enter fraction of particles in longitudinal profile that still count for shfoot: ";
            cin >> shfootlimit;
            cout << "------------------------------------------------------------------------------------" << endl;
         }
	 if( (directiveaffil[itemp[0]] == "tankvem") || (directiveaffil[itemp[1]] == "tankvem") )
	 {
            cout << " Enter maximal distance of SD tank from shower core (set to 0, to select all tanks): ";
            cin >> distlimit;
            cout << "------------------------------------------------------------------------------------" << endl;
	 }
      }

      // Run directive for the multiple file mode
      PrepareDirective(dirtemp, itemp, infilenr);

      return -1;
   }
/*
   // Run directive for the single file mode
   // Plotting observables
   if(analysistype == 1)
   {
      // Make 1D histogram
      if(plotdirective.find(':') == string::npos)
      {
         if( itemp[0] == -1 )
	    return -1;
	 else
	 {
            // 1D histogram directly from observable
	    stemp = dirtemp[0] + ">>" + "obs";

	    if(directiveaffil[itemp[0]] == "showtankdata")
               showtankdata->Draw(stemp.c_str());
	    else if(directiveaffil[itemp[0]] == "showeyedata")
               showeyedata->Draw(stemp.c_str());	// TODO: Still have to make a high limit for Xmax

	    stemp = string(BASEDIR) + "/results/" + dirtemp[0] + ".pdf";
	 }

         histf = (TH1F*)gDirectory->Get("obs");
	 dirtemp[0] = ";" + directivedesc[itemp[0]] + ";Number of events";
	 histf->SetTitle(dirtemp[0].c_str());
         histf->Draw();

         c1->SaveAs(stemp.c_str());

	 delete histf;
         
         return 0;
      }
      // Make 2D histogram
      else
      {
         if( (itemp[0] == -1) || (itemp[1] == -1) )
	    return -1;
	 else
	 {
            // Compare SD observables between eachother
	    if( (directiveaffil[itemp[0]] == "showtankdata") && (directiveaffil[itemp[1]] == "showtankdata") )
	    {
	       stemp = plotdirective + ">>" + "obs";
               showtankdata->Draw(stemp.c_str());

	       stemp = string(BASEDIR) + "/results/" + directive[itemp[0]] + "_" + directive[itemp[1]]+ ".pdf";
	    }
            // Compare FD observables between eachother
	    else if( (directiveaffil[itemp[0]] == "showeyedata") && (directiveaffil[itemp[1]] == "showeyedata") )
	    {
	       stemp = plotdirective + ">>" + "obs";
               showeyedata->Draw(stemp.c_str());

	       stemp = string(BASEDIR) + "/results/" + directive[itemp[0]] + "_" + directive[itemp[1]]+ ".pdf";
	    }
	    else
	    {
	       cout << "Can only compare observables between tanks or between eyes, not mixed ones." << endl;
	       return -1;
	    }

            gr = (TGraph*)gDirectory->Get("obs");
            gr->SetMarkerStyle(20);
            gr->SetMarkerSize(0.7);
	    dirtemp[0] = ";" + directivedesc[itemp[1]] + ";" + directivedesc[itemp[0]];
	    gr->SetTitle(dirtemp[0].c_str());
            gr->Draw();
            
            c1->SaveAs(stemp.c_str());
            
            delete gr;
	 }
         
         return 0;
      }
   }
   // Plotting profiles
   else*/ if(analysistype == 2)
   {
      double *x, *xerr, *y, *yerr;
      // We are plotting the LDF profile
      if( directive[itemp[0]] == "ldf" )
      {
         return PlotLDF(dirtemp, itemp);
      }
      // We are plotting a specific tank/eye profile
      else if( (itemp[0] != -1) && (itemp[1] == -1) )
      {
         // Plotting the VEM signal trace
         if(directive[itemp[0]] == "vemsig")
	 {
	    return PlotVEM(dirtemp, itemp);
	 }
         // Plotting the longitudinal profile
         else if(directive[itemp[0]] == "longprof")
	 {
	    return PlotLongProf(dirtemp, itemp);
	 }
	 else
            return -1;
      }
      // We are plotting tank/eye profiles for all tanks/eyes
      else if( (itemp[0] != -1) && (itemp[1] != -1) )
      {
         nrpoints = 0;
         if(directive[itemp[1]] == "all")
	 {
            // Plotting the VEM signal trace
            if(directive[itemp[0]] == "vemsig")
            {
	       for(int i = 0; i < acttanks.size(); i++)
	       {
	          if(acttanks[i].GetId() < 90000)
		  {
	             dirtemp[1] = IntToStr(acttanks[i].GetId());
	             nrpoints += PlotVEM(dirtemp, itemp);
		  }
	       }
	       return nrpoints;
            }
            // Plotting the longitudinal profile
            else if(directive[itemp[0]] == "longprof")
            {
	       for(int i = 0; i < ALLEYES; i++)
	       {
	          if(eyevalid[i] == 1)
		  {
	             dirtemp[1] = IntToStr(i+1);
	             nrpoints += PlotLongProf(dirtemp, itemp);
		  }
           	  else
                     cout << "No longitudinal profile to plot for eye with ID " << i+1 << endl;
	       }
               return nrpoints;
            }
            else
               return -1;
	 }
	 else
	    return -1;
      }
      else
         return -1;
   }
   // More to come
   else if(analysistype == 3)
   {
      cout << "Plotting more to come" << endl;
      return 0;
   }
   else
      return 1;
}

// Profile plotting instructions ------------------------------------------------------------------------
//    PlotLDF: *sepdir = directive types (verbal, i.e. ldf), *seldir = directive types (int)
//    PlotVEM: *sepdir = directive types (verbal, i.e. vemsig:5423), *seldir = directive types (int)
//    PlotLongProf: *sepdir = directive types (verbal, i.e. longprof:3), *seldir = directive types (int)
int AdstAnalyseTool::PlotLDF(string *sepdir, int *seldir)
{
   cout << "# Entering function AdstAnalyseTool::PlotLDF()..." << endl;
   double *x, *xerr, *y, *yerr;
   string *stemp;

   x = new double;
   xerr = new double;
   y = new double;
   yerr = new double;

   stemp = new string;

   SdRecShower &sdrecshw = fRecEvent->GetSDEvent().GetSdRecShower();
   LDF shwldf = sdrecshw.GetLDF();

   cout << "LDF @ 1000m = " << shwldf.Evaluate(1000, eNKG) << endl;

   grErr = new TGraphErrors(nracttanks);

   pointcnt = 0;

   for(int i = 0; i < acttanks.size(); i++)
   {
      if(acttanks[i].GetId() < 90000)
      {
         *x = acttanks[i].GetSPDistance();
         *xerr = acttanks[i].GetSPDistanceError();
         *y = acttanks[i].GetTotalSignal();
         *yerr = acttanks[i].GetTotalSignalError();

         grErr->SetPoint(pointcnt,*x,*y);
         grErr->SetPointError(pointcnt,*xerr,*yerr);

	 pointcnt++;
      }
   }

   TF1 ldffunc = (TF1)shwldf.GetFunction(eNKG);
   ldffunc.Draw();

   c1->SetLogy(1);
   *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_profile.pdf";
   grErr->SetMarkerStyle(20);
   grErr->SetMarkerSize(0.7);
   sepdir[0] = directivedesc[seldir[0]] + ";Distance from axis (m);Signal (VEM)";
   grErr->SetTitle(sepdir[0].c_str());
   grErr->Draw("P;SAME");
      
   c1->SaveAs((*stemp).c_str());

   *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_profile.C";
   c1->SaveAs((*stemp).c_str());

   delete grErr;
   delete x;
   delete xerr;
   delete y;
   delete yerr;
   delete stemp;

   c1->SetLogy(0);

   return 0;
}

int AdstAnalyseTool::PlotVEM(string *sepdir, int *seldir)
{
   cout << "# Entering function AdstAnalyseTool::PlotVEM()..." << endl;
   string *stemp;
   double *range;

   stemp = new string;
   range = new double[3];

   c1 = new TCanvas("c1","c1",1200,800);
   c1->SetGrid();

   vector<Traces> vemtraces;
   vector<float> trace1;
   vector<float> trace2;
   vector<float> trace3;

   for(int i = 0; i < acttanks.size(); i++)
   {
      if(acttanks[i].GetId() == stoi(sepdir[1]))
      {
	 vemtraces = acttanks[i].GetPMTTraces();
	 trace1 = vemtraces[0].GetVEMComponent();
	 trace2 = vemtraces[1].GetVEMComponent();
	 trace3 = vemtraces[2].GetVEMComponent();

         range[0] = acttanks[i].GetSignalStartSlot();
         range[1] = acttanks[i].GetSignalEndSlot();
      }

      printf("\nSecond = %.7le\n", acttanks[i].GetTimeSecond());
      printf("Nanosecond = %.7le\n", acttanks[i].GetTimeNSecond());
      cout << "Signal start slot = " << acttanks[i].GetSignalStartSlot() << endl;
      cout << "Signal end slot = " << acttanks[i].GetSignalEndSlot() << endl;
   }

   gr = new TGraph(trace1.size());

   for(int i = 0; i < trace1.size(); i++)
   {
      gr->SetPoint(i,i,(trace1[i]+trace2[i]+trace3[i]));
   }

   // Prepare export filename and marker style
   *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_" + sepdir[1] + "_profile.pdf";
   gr->SetMarkerStyle(20);
   gr->SetMarkerSize(0.7);
   // Set the range of the plot and title
   range[2] = (range[1] - range[0])*0.2;
   gr->GetXaxis()->SetRangeUser( range[0]-range[2], range[1]+range[2] );
   sepdir[0] = directivedesc[seldir[0]] + " (ID = " + sepdir[1] + ")" + ";SD tank time (25 ns);Signal (VEM)";
   gr->SetTitle(sepdir[0].c_str());
   gr->Draw("APL");

   c1->SaveAs((*stemp).c_str());

   *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_" + sepdir[1] + "_profile.C";
   c1->SaveAs((*stemp).c_str());
   
   delete gr;
   delete stemp;
   delete range;
//   delete c1;

   return 0;
}

int AdstAnalyseTool::PlotLongProf(string *sepdir, int *seldir)
{
   cout << "# Entering function AdstAnalyseTool::PlotLongProf()..." << endl;
   double *x, *xerr, *y, *yerr;
   string *stemp, *stemp2;
   double *cum, *cumerr;
   int *nrpoints;

   x = new double;
   xerr = new double;
   y = new double;
   yerr = new double;
   cum = new double;
   cumerr = new double;
   nrpoints = new int;

   stemp = new string;
   stemp2 = new string;

   c1 = new TCanvas("c1","c1",1200,800);
   c1->SetGrid();
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);

   vector<FDEvent> fdevt = fRecEvent->GetFDEvents();

   vector<double> slantDepth;
   vector<double> profiledEdX;
   vector<double> profiledEdXerr;

   for(int i = 0; i < fdevt.size(); i++)
   {
      if(fdevt[i].GetEyeId() == stoi(sepdir[1]))
      {
         slantDepth = acteyes[i].GetDepth();
         profiledEdX = acteyes[i].GetEnergyDeposit();
         profiledEdXerr = acteyes[i].GetEnergyDepositError();
      }
   }

   grErr = new TGraphErrors(slantDepth.size());

   if(btemp)	// only valid if we want separate cumulative longitudinal distribution plots
   {
      grErrtemp = new TGraphErrors(slantDepth.size());
      *cum = 0;
      if( xfoot.size() != 0 )
         xfoot.erase(xfoot.begin(),xfoot.end());
      if( yfoot.size() != 0 )
         yfoot.erase(yfoot.begin(),yfoot.end());
      if( yerrfoot.size() != 0 )
         yerrfoot.erase(yerrfoot.begin(),yerrfoot.end());
   }

   for(int i = 0; i < slantDepth.size(); i++)
   {
      grErr->SetPoint(i,slantDepth[i],profiledEdX[i]);
      grErr->SetPointError(i,1.e-7,profiledEdXerr[i]);

      if(btemp)
      {
         *cum += profiledEdX[i];
         *cumerr += profiledEdXerr[i];

	 xfoot.push_back(slantDepth[i]);
	 yfoot.push_back(*cum);
	 yerrfoot.push_back(*cumerr);

	 grErrtemp->SetPoint(i,slantDepth[i],*cum);
         grErrtemp->SetPointError(i,1.e-7,*cumerr);

//	 cout << "Point " << i << ": " << xfoot[i] << "\t" << yfoot[i] << " (" << yerrfoot[i] << ")" << endl;
      }
   }

   // Setting output names for plots
   if(stoi(sepdir[1]) == 1)
   {
      *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_LL_profile.pdf";
      *stemp2 = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_LL_profile.C";
   }
   else if(stoi(sepdir[1]) == 2)
   {
      *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_LM_profile.pdf";
      *stemp2 = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_LM_profile.C";
   }
   else if(stoi(sepdir[1]) == 3)
   {
      *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_LA_profile.pdf";
      *stemp2 = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_LA_profile.C";
   }
   else if(stoi(sepdir[1]) == 4)
   {
      *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_CO_profile.pdf";
      *stemp2 = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_CO_profile.C";
   }
   else
   {
      *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_" + sepdir[1] + "_profile.pdf";
      *stemp2 = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_" + sepdir[1] + "_profile.C";
   }

   grErr->SetMarkerStyle(20);
   grErr->SetMarkerSize(0.7);
   sepdir[0] = directivedesc[seldir[0]] + ";Slant depth (g/cm^{2});dE/dX (PeV/(g/cm^{2}))";
   grErr->SetTitle(sepdir[0].c_str());
   grErr->GetXaxis()->SetTitleOffset(1.2);
   grErr->GetXaxis()->CenterTitle(kTRUE);
   grErr->GetXaxis()->SetLabelSize(0.028);
   grErr->GetXaxis()->SetLabelOffset(0.015);
   grErr->GetYaxis()->SetTitleOffset(1.3);
   grErr->GetYaxis()->CenterTitle(kTRUE);
   grErr->GetYaxis()->SetLabelSize(0.028);
   grErr->GetYaxis()->SetLabelOffset(0.015);
   grErr->SetTitle("");
   grErr->Draw("AP");

   // Plot also the generated profile on top of it
   GenShower genevt = fRecEvent->GetGenShower();
   slantDepth = genevt.GetDepth();
   profiledEdX = genevt.GetEnergyDeposit();
   gr = new TGraph(slantDepth.size());

   for(int i = 0; i < slantDepth.size(); i++)
      gr->SetPoint(i,slantDepth[i],profiledEdX[i]);
   gr->SetLineColor(kRed);
   gr->SetLineWidth(2);
   gr->Draw("LP;SAME");
      
   c1->SaveAs((*stemp).c_str());
   c1->SaveAs((*stemp2).c_str());
//-------------------------------------------
   // Plotting instructions for the cumulative distribution
   if(btemp)
   {
      *nrpoints = 0;

      for(int i = 0; i < yfoot.size(); i++)
      {
         if( (yfoot[i] >= shfootlimit*(yfoot[yfoot.size()-1])) && (*nrpoints == 0) )
         {
	    *nrpoints = 1;

	    // Find the x value of point with y value that is a fraction of the maximum, that lies on a line between two points
	    // y = k*x + a
	    //    k = (y2 - y1)/(x2 - x1)
	    //    a = y2 - (y2 - y1)/(x2 - x1)*x2
	    // x = ((x2 - x1)/(y2 - y1))*(y - y2) + x2
	    *cum = (xfoot[i] - xfoot[i-1])/(yfoot[i] - yfoot[i-1]); // 1/k = (x2 - x1)/(y2 - y1)
	    *cumerr = shfootlimit*(yfoot[yfoot.size()-1]) - yfoot[i]; // y - y2
            *x = (*cum)*(*cumerr) + xfoot[i]; // x = (1/k)*(y - y2) + x2

//	    cout << *cum << "\t" << *cumerr << "\t" << *x << endl;

            *cum = (xfoot[i] - xfoot[i-1])/(yfoot[i]+yerrfoot[i] - (yfoot[i-1]+yerrfoot[i-1])); // 1/kerr = (x2 - x1)/(y2err - y1err)
	    *cumerr = (yfoot[i]+yerrfoot[i]) - (1/(*cum))*(xfoot[i]); // aerr = y2err - (y2err - y1err)/(x2 - x1)*x2
	    *cum = (1/(*cum))*(*x) + (*cumerr); // yerr = kerr*x + aerr
	    *yerr = (*cum) - (((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(*x) + (yfoot[i] - ((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(xfoot[i]))); // Dy = yerr - y

	    *y = ((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(*x) + (yfoot[i] - ((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(xfoot[i])); // y = k*x + a

            for(int j = i; ; j++)
	    {
	       if( yfoot[j] >= (*y)+(*yerr) )
	       {
//	          cout << "Upper limit is = " << xfoot[j] << endl;
		  *cumerr = xfoot[j];
		  break;
	       }
	    }
            for(int j = i; ; j--)
	    {
               if(j == 0)
               {
	          *xerr = xfoot[j];
//	          cout << "Lower limit is = " << xfoot[j] << endl;
     	          break;
               }

	       if( yfoot[j] <= (*y)-(*yerr) )
	       {
//	          cout << "Lower limit is = " << xfoot[j] << endl;
		  *xerr = xfoot[j];
		  break;
	       }
	    }

	    break;
         }
      }

      grErrtemp->SetTitle(";Slant depth (g/cm^{2});Cumulative dE/dX (PeV/(g/cm^{2}))");
      grErrtemp->GetXaxis()->SetTitleOffset(1.2);
      grErrtemp->GetXaxis()->CenterTitle(kTRUE);
      grErrtemp->GetXaxis()->SetLabelSize(0.028);
      grErrtemp->GetXaxis()->SetLabelOffset(0.015);
      grErrtemp->GetYaxis()->SetTitleOffset(1.4);
      grErrtemp->GetYaxis()->CenterTitle(kTRUE);
      grErrtemp->GetYaxis()->SetLabelSize(0.028);
      grErrtemp->GetYaxis()->SetLabelOffset(0.015);
      grErrtemp->Draw("AL");
      grErrtemp->SetLineColor(kGray+2);

      cout << "Marker at = " << *x << ", " << *xerr << " | " << *y << ", " << *yerr << endl;
      cout << "Box at = " << *xerr << ", " << (*y)+(*yerr) << " | " << *cumerr << ", " << (*y)-(*yerr) << endl;

      // Draw point and a red box do represent the shower foot limiting value and its error
/*      TBox *b1 = new TBox(*xerr, (*y)+(*yerr), *cumerr, (*y)-(*yerr));
      b1->SetLineColor(kRed);
      b1->SetFillColorAlpha(kRed, 0.3);
      b1->Draw();*/

      TGraphAsymmErrors *m1 = new TGraphAsymmErrors(1);
      m1->SetPoint(0,*x,*y);
      m1->SetPointError(0, (*x)-(*xerr), (*cumerr)-(*x), (*yerr), (*yerr));
      m1->SetMarkerColor(kRed);
      m1->SetMarkerSize(0.8);
      m1->SetLineColor(kRed);
      m1->Draw("SAME");

      if(stoi(sepdir[1]) == 1)
         *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_LL_cumulative_profile.pdf";
      else if(stoi(sepdir[1]) == 2)
         *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_LM_cumulative_profile.pdf";
      else if(stoi(sepdir[1]) == 3)
         *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_LA_cumulative_profile.pdf";
      else if(stoi(sepdir[1]) == 4)
         *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_CO_cumulative_profile.pdf";
      else
         *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_" + sepdir[1] + "_cumulative_profile.pdf";
      c1->SaveAs((*stemp).c_str());

      delete grErrtemp;
      delete m1;
   }

   delete grErr;
   delete x;
   delete xerr;
   delete y;
   delete yerr;
   delete cum;
   delete cumerr;
   delete nrpoints;
   delete stemp;
   delete stemp2;
//   delete c1;*/

   return 0;
}

// Run directive for the multiple file mode -------------------------------------------------------------
int AdstAnalyseTool::PrepareDirective(string *sepdir, int *seldir, int infilenr)
{
   cout << "# Entering function AdstAnalyseTool::PrepareDirective()..." << endl;

   double *x, *xerr, *y, *yerr;
   double *derr;
   string *stemp;
   int *itemp;
   double *dtemp;
   int *nrpoints;
   double *cordist;

   x = new double;
   xerr = new double;
   y = new double;
   yerr = new double;

   nrpoints = new int;
   stemp = new string;
   itemp = new int;
   dtemp = new double;
   derr = new double[2];
   cordist = new double;

   if(argtype == 0)
      fFile->ReadEvent(infilenr);

   cout << "Eye status: ";
   for(int i = 0; i < ALLEYES; i++)
      cout << eyevalid[i] << ", ";
   cout << endl;

   // Open up a file to print out results in an ASCII file
   if(infilenr == 0)
   {
      if( (seldir[0] != -1) && (seldir[1] == -1) && (seldir[2] == -1) && (nrdirs == 1) )
      {
         *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_multi" +  + ".dat";
         outdata.open((*stemp).c_str(), ofstream::trunc | ofstream::out);
      }
      else if( (seldir[0] != -1) && (seldir[1] != -1) && (seldir[2] == -1) && (nrdirs == 2) )
      {
         *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_vs_" + directive[seldir[1]] + "_multi" +  + ".dat";
         outdata.open((*stemp).c_str(), ofstream::trunc | ofstream::out);
      }
   }

   // Integer values of the selected directives
   cout << "Seldir: " << seldir[0] << ", " << seldir[1] << ", " << seldir[2] << ". (" << nrdirs << ")" << endl;

   if( xmaxerr.size() != 0 )
      xmaxerr.erase(xmaxerr.begin(), xmaxerr.end());

   // 1D histogram
   if( (seldir[0] != -1) && (seldir[1] == -1) && (seldir[2] == -1) && (nrdirs == 1) )
   {
      // Getting separate tank data (totalvem, distance, muoncount)
      if( directiveaffil[seldir[0]] == "showtankdata" )
      {
	 *y = 0;
	 *itemp = 0;

         for(int i = 0; i < acttanks.size(); i++)
         {
	    if(acttanks[i].GetId() < 90000)
	    {
               cout << "Active tank: " << acttanks[i].GetId() << endl;

	       GetData(sepdir, seldir, "Xaxis", y, yerr, itemp, i);
	    }
         }

         if(acttanks.size() > 0)
	 {
            *dtemp = *y/(double)*itemp;

	    // Write out to data file for additional analysis (histogram bin widths)
	    outdata << *dtemp << endl;

            // Set ranges
            if( *dtemp < xrange[0] )
	       xrange[0] = *dtemp;
            if( *dtemp > xrange[1] )
	       xrange[1] = *dtemp;

            cout << "The value to plot " << *dtemp << endl;
            histf->Fill(*dtemp);
	 }

	 if( infilenr+1 == evtcount )
	 {
            histf->GetXaxis()->SetRangeUser(xrange[0]-0.1*(xrange[1]-xrange[0]),xrange[1]+0.1*(xrange[1]-xrange[0]));
	    *stemp = ";" + directivedesc[seldir[0]] + ";Number of events";
  	    histf->SetTitle((*stemp).c_str());
	    histf->Draw("");

	    *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_multi" +  + ".pdf";
            c1->SaveAs((*stemp).c_str());

	    *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_multi" + ".C";
            c1->SaveAs((*stemp).c_str());
	 }
      }
      // Getting separate eye data (energy, emenergy, xmax)
      if( directiveaffil[seldir[0]] == "showeyedata" ) // TODO: In addition, maybe also incorporate goodness of fit in determining a quality xmax
      {
	 *y = 0;
         if(directive[seldir[0]] == "xmaxqual")
	    *yerr = 1.e+40;
	 else
	    *yerr = 0;
	 *itemp = 0;

         for(int i = 0; i < ALLEYES; i++)
         {
	    if(eyevalid[i] != 0)
	    {
	       if(eyevalid[i] == 1)
	       {
                  cout << "Active eye: " << i+1 << endl;

	          GetData(sepdir, seldir, "Xaxis", y, yerr, itemp, i);
	       }
	    }
	    else
	       cout << "No FD data to plot for eye ID " << i+1 << endl;
         }

         if( acteyes.size() > 0 )
	 {
            if(directive[seldir[0]] == "xmaxqual")
	       *dtemp = *y;
	    else
               *dtemp = *y/(double)*itemp;

	    // Write out to data file for additional analysis (histogram bin widths)
	    outdata << *dtemp << endl;

            // Set ranges
            if( *dtemp < xrange[0] )
	       xrange[0] = *dtemp;
            if( *dtemp > xrange[1] )
	       xrange[1] = *dtemp;

            cout << "The value to plot " << *dtemp << endl;
            histf->Fill(*dtemp);
	 }

	 if( infilenr+1 == evtcount )
	 {
            histf->GetXaxis()->SetRangeUser(xrange[0]-0.1*(xrange[1]-xrange[0]),xrange[1]+0.1*(xrange[1]-xrange[0]));
	    *stemp = ";" + directivedesc[seldir[0]] + ";Number of events";
  	    histf->SetTitle((*stemp).c_str());
	    histf->Draw("");

	    *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_multi" +  + ".pdf";
            c1->SaveAs((*stemp).c_str());

	    *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_multi" + ".C";
            c1->SaveAs((*stemp).c_str());
	 }
      }
      // Getting shower foot information (shfoot)
      if( directiveaffil[seldir[0]] == "eyelong" )
      {
	 *y = 0;
	 derr[0] = 0; derr[1] = 0;
	 *itemp = 0;

         for(int i = 0; i < ALLEYES; i++)
         {
	    if(eyevalid[i] != 0)
	    {
	       if(eyevalid[i] == 1)
	       {
                  cout << "Active eye: " << i+1 << endl;

	          GetData(sepdir, seldir, "Xaxis", y, derr, itemp, i);
	       }
	    }
	    else
	       cout << "No FD data to plot for eye ID " << i+1 << endl;
         }

         if( acteyes.size() > 0 )
	 {
            *dtemp = *y/(double)*itemp;

	    // Write out to data file for additional analysis (histogram bin widths)
	    outdata << *dtemp << endl;

            // Set ranges
            if( *dtemp < xrange[0] )
	       xrange[0] = *dtemp;
            if( *dtemp > xrange[1] )
	       xrange[1] = *dtemp;

            cout << "The value to plot " << *dtemp << endl;
            histf->Fill(*dtemp);
	 }

	 if( infilenr+1 == evtcount )
	 {
            histf->GetXaxis()->SetRangeUser(xrange[0]-0.1*(xrange[1]-xrange[0]),xrange[1]+0.1*(xrange[1]-xrange[0]));
            *stemp = ";Depth of " + IntToStr((int)(shfootlimit*100)) + directivedesc[seldir[0]] + ";Number of events";
  	    histf->SetTitle((*stemp).c_str());
	    histf->Draw("");

	    *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_multi" +  + ".pdf";
            c1->SaveAs((*stemp).c_str());

	    *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_multi" + ".C";
            c1->SaveAs((*stemp).c_str());
	 }
      }
      // Getting SD reconstructed data (actstation, energy, S1000, curvature)
      if( directiveaffil[seldir[0]] == "showsrecdata" )
      {
	 *y = 0;
	 *itemp = 0;

	 GetData(sepdir, seldir, "Xaxis", y, yerr, itemp, 0);

         // Set ranges
         if( *y < xrange[0] )
	    xrange[0] = *y;
         if( *y > xrange[1] )
	    xrange[1] = *y;

	 // Write out to data file for additional analysis (histogram bin widths)
	 outdata << *y << endl;

         cout << "The value to plot " << *y << endl;
         histf->Fill(*y);

	 if( infilenr+1 == evtcount )
	 {
            histf->GetXaxis()->SetRangeUser(xrange[0]-0.1*(xrange[1]-xrange[0]),xrange[1]+0.1*(xrange[1]-xrange[0]));
	    *stemp = ";" + directivedesc[seldir[0]] + ";Number of events";
  	    histf->SetTitle((*stemp).c_str());
	    histf->Draw("");

	    *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_multi" +  + ".pdf";
            c1->SaveAs((*stemp).c_str());

	    *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_multi" + ".C";
            c1->SaveAs((*stemp).c_str());
	 }
      }
      // Getting FD reconstructed data (acteyes)
      if( directiveaffil[seldir[0]] == "showfrecdata" )
      {
	 *y = 0;
	 *itemp = 0;

	 GetData(sepdir, seldir, "Xaxis", y, yerr, itemp, 0);

         // Set ranges
         if( *y < xrange[0] )
	    xrange[0] = *y;
         if( *y > xrange[1] )
	    xrange[1] = *y;

	 // Write out to data file for additional analysis (histogram bin widths)
	 outdata << *y << endl;

         cout << "The value to plot " << *y << endl;
         histf->Fill(*y);

	 if( infilenr+1 == evtcount )
	 {
            histf->GetXaxis()->SetRangeUser(xrange[0]-0.1*(xrange[1]-xrange[0]),xrange[1]+0.1*(xrange[1]-xrange[0]));
	    *stemp = ";" + directivedesc[seldir[0]] + ";Number of events";
  	    histf->SetTitle((*stemp).c_str());
	    histf->Draw("");

	    *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_multi" +  + ".pdf";
            c1->SaveAs((*stemp).c_str());

	    *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_multi" + ".C";
            c1->SaveAs((*stemp).c_str());
	 }
      }
      // Getting separate tank calculation (risetime, risetimeearly, risetimelate)
      if( directiveaffil[seldir[0]] == "tankvem" )
      {
	 *y = 10.e+31;
	 *itemp = 0;

         if( sdtimesec.size() != 0 )
            sdtimesec.erase(sdtimesec.begin(),sdtimesec.end());
         if( sdtimensec.size() != 0 )
            sdtimensec.erase(sdtimensec.begin(),sdtimensec.end());

         for(int i = 0; i < acttanks.size(); i++)
	 {
	    if(acttanks[i].GetId() < 90000)
	    {
	       sdtimesec.push_back(acttanks[i].GetTimeSecond());
	       sdtimensec.push_back(acttanks[i].GetTimeNSecond());

//               printf("Second: %lf, NSecond: %lf, Distance: %lf\n", acttanks[i].GetTimeSecond(), acttanks[i].GetTimeNSecond(), acttanks[i].GetSPDistance());

               if( acttanks[i].GetSPDistance() < *y )
	       {
	          *y = acttanks[i].GetSPDistance();
	          sdhotnsec = sdtimensec[i];
               }
	    }
	 }

	 printf("Hottest tank time: %lf\n", sdhotnsec);

	 *y = 0;

         for(int i = 0; i < acttanks.size(); i++)
         {
	    if(acttanks[i].GetId() < 90000)
	    {
               cout << "Active tank: " << acttanks[i].GetId() << endl;

               if( (distlimit != 0 ) && (acttanks[i].GetSPDistance() <= distlimit) )
                  GetData(sepdir, seldir, "Xaxis", y, yerr, itemp, i);
	       else if( distlimit == 0 )
                  GetData(sepdir, seldir, "Xaxis", y, yerr, itemp, i);
	    }
         }

         if(acttanks.size() > 0)
	 {
            *dtemp = *y/(double)*itemp;

	    // Write out to data file for additional analysis (histogram bin widths)
	    outdata << *dtemp << endl;

            // Set ranges
            if( *dtemp < xrange[0] )
	       xrange[0] = *dtemp;
            if( *dtemp > xrange[1] )
	       xrange[1] = *dtemp;

            cout << "The value to plot " << *dtemp << endl;
            histf->Fill(*dtemp);
	 }

	 if( infilenr+1 == evtcount )
	 {
            histf->GetXaxis()->SetRangeUser(xrange[0]-0.1*(xrange[1]-xrange[0]),xrange[1]+0.1*(xrange[1]-xrange[0]));
	    *stemp = ";" + directivedesc[seldir[0]] + ";Number of events";
  	    histf->SetTitle((*stemp).c_str());
	    histf->Draw("");

	    *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_multi" +  + ".pdf";
            c1->SaveAs((*stemp).c_str());

	    *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_multi" + ".C";
            c1->SaveAs((*stemp).c_str());
	 }
      }
   }
   // Scatter plot
   else if( (seldir[0] != -1) && (seldir[1] != -1) && (seldir[2] == -1) && (nrdirs == 2) )
   {
      // Go through both directives (for X and Y axis)
/*------------------------------------------------------------------------------------------*/
      for(int k = 0; k < nrdirs; k++)
      {
         // Getting separate tank data (totalvem, distance, muoncount)
         if( directiveaffil[seldir[k]] == "showtankdata" ) // probably still not working correctly
         {
	    if(k == 0) *x = 0;
            else *y = 0;
            *itemp = 0;

            for(int i = 0; i < acttanks.size(); i++)
            {
               if(acttanks[i].GetId() < 90000)
               {
                  cout << "Active tank: " << acttanks[i].GetId() << endl;

                  if(k == 0)
                     GetData(sepdir, seldir, "Xaxis", x, xerr, itemp, i);
	          else
                     GetData(sepdir, seldir, "Yaxis", y, yerr, itemp, i);
	       }
            }

            if(acttanks.size() > 0)
            {
	       if(k == 0)
	       {
                  *x = *x/(double)*itemp;
                  *xerr = *xerr/(double)*itemp;

                  // Write out to data file for additional analysis (histogram bin widths)
                  outdata << *x << "\t" << *xerr << "\t" << *xerr << "\t";
	       }
	       else
	       {
                  *y = *y/(double)*itemp;
                  *yerr = *yerr/(double)*itemp;

                  // Write out to data file for additional analysis (histogram bin widths)
                  outdata << *y << "\t" << *yerr << "\t" << *yerr << endl;
	       }
            }
         }
         // Getting separate eye data (energy, emenergy, xmax)
         if( directiveaffil[seldir[k]] == "showeyedata" )
         {
	    if(k == 0) *x = 0;
            else *y = 0;
	    if(k == 0) *xerr = 0;
            else *yerr = 0;
            *itemp = 0;

            for(int i = 0; i < ALLEYES; i++)
            {
               if(eyevalid[i] != 0)
               {
                  if(eyevalid[i] == 1)
                  {
                     cout << "Active eye: " << i+1 << endl;

                     if(k == 0)
                        GetData(sepdir, seldir, "Xaxis", x, xerr, itemp, i);
	             else
                        GetData(sepdir, seldir, "Yaxis", y, yerr, itemp, i);
                  }
               }
               else
                  cout << "No FD data to plot for eye ID " << i+1 << endl;
            }

            if( acteyes.size() > 0 )
            {
	       if(k == 0)
	       {
/*	          if(*itemp == 0)
		  {
		     *x = 0;
                     *xerr = *xerr/(double)*itemp;
		  }
		  else
		  {
*/                     if(directive[seldir[k]] == "xmaxqual")
                     {
                        *x = *x;
                        *xerr = *xerr;
	             }
		     else
		     {
                        *x = *x/(double)*itemp;
                        *xerr = *xerr/(double)*itemp;
	             }

                  // Write out to data file for additional analysis (histogram bin widths)
                  outdata << *x << "\t" << *xerr << "\t" << *xerr << "\t";
//	          }
	       }
	       else
	       {
/*	          if(*itemp == 0)
		  {
		     *y = 0;
                     *yerr = *yerr/(double)*itemp;
		  }
		  else
		  {
*/                     if(directive[seldir[k]] == "xmaxqual")
                     {
                        *y = *y;
                        *yerr = *yerr;
		     }
		     else
		     {
                        *y = *y/(double)*itemp;
                        *yerr = *yerr/(double)*itemp;
	             }

                  // Write out to data file for additional analysis (histogram bin widths)
                  outdata << *y << "\t" << *yerr << "\t" << *yerr << endl;
//	          }
	       }
            }
         }
         // Getting shower foot information (shfoot)
         if( directiveaffil[seldir[k]] == "eyelong" )
         {
	    if(k == 0) *x = 0;
            else *y = 0;
	    if(k == 0) {derr[0] = 0; derr[1] = 0;}
            else {derr[0] = 0; derr[1] = 0;}
            *itemp = 0;

            for(int i = 0; i < ALLEYES; i++)
            {
               if(eyevalid[i] != 0)
               {
                  if(eyevalid[i] == 1)
                  {
                     cout << "Active eye: " << i+1 << endl;

                     if( (k == 0) && (directive[seldir[k]] == "shfoot") )
                        GetData(sepdir, seldir, "Xaxis", x, derr, itemp, i);
                     else if( (k == 0) && (directive[seldir[k]] != "shfoot") )
                        GetData(sepdir, seldir, "Xaxis", x, xerr, itemp, i);
                     else if( (k == 1) && (directive[seldir[k]] == "shfoot") )
                        GetData(sepdir, seldir, "Yaxis", y, derr, itemp, i);
                     else if( (k == 1) && (directive[seldir[k]] != "shfoot") )
                        GetData(sepdir, seldir, "Yaxis", y, yerr, itemp, i);
                  }
               }
               else
                  cout << "No FD data to plot for eye ID " << i+1 << endl;
            }

            if( acteyes.size() > 0 )
            {
	       if(k == 0)
	       {
/*	          if(*itemp == 0)
		     *x = 0;
		  else
*/                     *x = *x/(double)*itemp;

                  // Write out to data file for additional analysis (histogram bin widths)
                  outdata << *x << "\t";

		  if(directive[seldir[k]] == "shfoot")
		  {
 		     derr[0] = derr[0]/(double)*itemp;
 		     derr[1] = derr[1]/(double)*itemp;

                     outdata << derr[0] << "\t" << derr[1] << "\t";
		  }
		  else
		  {
                     *xerr = *xerr/(double)*itemp;
                     outdata << *xerr << "\t" << *xerr << "\t";
		  }
	       }
	       else
	       {
/*	          if(*itemp == 0)
		     *y = 0;
		  else
*/                     *y = *y/(double)*itemp;

                  // Write out to data file for additional analysis (histogram bin widths)
                  outdata << *y << "\t";

		  if(directive[seldir[k]] == "shfoot")
		  {
 		     derr[0] = derr[0]/(double)*itemp;
 		     derr[1] = derr[1]/(double)*itemp;

                     outdata << derr[0] << "\t" << derr[1] << endl;
		  }
		  else
		  {
                     *yerr = *yerr/(double)*itemp;
                     outdata << *yerr << "\t" << *yerr << endl;
		  }
	       }
            }
         }
         // Getting SD reconstructed data (actstation, energy, S1000, curvature)
         if( directiveaffil[seldir[k]] == "showsrecdata" )
         {
	    if(k == 0) *x = 0;
	    else *y = 0;
	    if(k == 0) *xerr = 0;
            else *yerr = 0;
	    *itemp = 0;

            if(k == 0)
	    {
	       GetData(sepdir, seldir, "Xaxis", x, xerr, itemp, 0);

               outdata << *x << "\t" << *xerr << "\t" << *xerr << "\t";
	    }
            else
	    {
	       GetData(sepdir, seldir, "Yaxis", y, yerr, itemp, 0);

               outdata << *y << "\t" << *yerr << "\t" << *yerr << endl;
	    }
	 }
         // Getting FD reconstructed data (acteyes)
         if( directiveaffil[seldir[k]] == "showfrecdata" )
         {
	    if(k == 0) *x = 0;
	    else *y = 0;
	    if(k == 0) *xerr = 0;
            else *yerr = 0;
	    *itemp = 0;

            if(k == 0)
	    {
	       GetData(sepdir, seldir, "Xaxis", x, xerr, itemp, 0);

               outdata << *x << "\t" << *xerr << "\t" << *xerr << "\t";
	    }
            else
	    {
	       GetData(sepdir, seldir, "Yaxis", y, yerr, itemp, 0);

               outdata << *y << "\t" << *yerr << "\t" << *yerr << endl;
	    }
	 }
         // Getting separate tank calculation (risetime, risetimeearly, risetimelate)
         if( directiveaffil[seldir[0]] == "tankvem" )
         {
            *y = 10.e+31;
            *itemp = 0;

            if( sdtimesec.size() != 0 )
               sdtimesec.erase(sdtimesec.begin(),sdtimesec.end());
            if( sdtimensec.size() != 0 )
               sdtimensec.erase(sdtimensec.begin(),sdtimensec.end());

            for(int i = 0; i < acttanks.size(); i++)
            {
               if(acttanks[i].GetId() < 90000)
               {
                  sdtimesec.push_back(acttanks[i].GetTimeSecond());
                  sdtimensec.push_back(acttanks[i].GetTimeNSecond());

//                  printf("Second: %lf, NSecond: %lf, Distance: %lf\n", acttanks[i].GetTimeSecond(), acttanks[i].GetTimeNSecond(), acttanks[i].GetSPDistance());

                  if( acttanks[i].GetSPDistance() < *y )
                  {
                     *y = acttanks[i].GetSPDistance();
                     sdhotnsec = sdtimensec[i];
                  }
               }
            }

            printf("Hottest tank time: %lf\n", sdhotnsec);

	    if(k == 0) *x = 0;
	    else *y = 0;
	    if(k == 0) *xerr = 0;
            else *yerr = 0;

            for(int i = 0; i < acttanks.size(); i++)
            {
               if(acttanks[i].GetId() < 90000)
               {
                  cout << "Active tank: " << acttanks[i].GetId() << endl;

                  if( (distlimit != 0 ) && (acttanks[i].GetSPDistance() <= distlimit) )
		  {
                     if(k == 0)
                        GetData(sepdir, seldir, "Xaxis", x, xerr, itemp, i);
		     else
                        GetData(sepdir, seldir, "Yaxis", y, yerr, itemp, i);
	          }
                  else if( distlimit == 0 )
		  {
                     if(k == 0)
                        GetData(sepdir, seldir, "Xaxis", x, xerr, itemp, i);
		     else
                        GetData(sepdir, seldir, "Yaxis", y, yerr, itemp, i);
	          }
               }
            }

            if(acttanks.size() > 0)
            {
	       if(k == 0)
	       {
                  *x = *x/(double)*itemp;
                  *xerr = *xerr/(double)*itemp;

                  outdata << *x << "\t" << *xerr << "\t" << *xerr << "\t";
	       }
	       else
	       {
                  *y = *y/(double)*itemp;
                  *yerr = *yerr/(double)*itemp;

                  outdata << *y << "\t" << *yerr << "\t" << *yerr << endl;
	       }
            }
         }

         // Set ranges
         if(directivetype[seldir[0]] == 'S')
	    SetRanges2D(seldir, x, derr, y, yerr, k);
         else if(directivetype[seldir[1]] == 'S')
	    SetRanges2D(seldir, x, xerr, y, derr, k);
	 else
	    SetRanges2D(seldir, x, xerr, y, yerr, k);
      }// end of both directives (for k)

      // Select what to plot	// TODO: Making scatter plots still not correct
      if( (directivetype[seldir[0]] == 'E') || (directivetype[seldir[1]] == 'E') )
      {
         if( (directiveaffil[seldir[0]] == "showeyedata") || (directiveaffil[seldir[1]] == "showeyedata") ) // don't write out zero values for eye values
	 {
	    if( ((directiveaffil[seldir[0]] == "showeyedata") && ((int)*x == 0)) || ((directiveaffil[seldir[1]] == "showeyedata") && ((int)*y == 0)) || (nracteyes == 0) )
	       cout << "No values to plot." << endl;
	    else
	    {
	       if(directive[seldir[0]] == "shfoot")
	       {
                  grAsymmErr->SetPoint(pointcnt,*x,*y);
                  grAsymmErr->SetPointError(pointcnt,derr[0],derr[1],*yerr,*yerr);
                  cout << "The value to plot " << *x << " (+" << derr[1] << "-" << derr[0] << "), " << *y << " (" << *yerr << ")" << endl;

	          pointcnt++;
	       }
	       else if(directive[seldir[1]] == "shfoot")
	       {
                  grAsymmErr->SetPoint(pointcnt,*x,*y);
                  grAsymmErr->SetPointError(pointcnt,*xerr,*xerr,derr[0],derr[1]);
                  cout << "The value to plot " << *x << " (" << *xerr << "), " << *y << " (+" << derr[1] << "-" << derr[0] << ")" << endl;

	          pointcnt++;
	       }
	       else
	       {
                  grErr->SetPoint(pointcnt,*x,*y);
                  grErr->SetPointError(pointcnt,*xerr,*yerr);
                  cout << "The value to plot " << *x << " (" << *xerr << "), " << *y << " (" << *yerr << ")" << endl;

	          pointcnt++;
	       }
	    }
	 }
	 else
	 {
            grErr->SetPoint(pointcnt,*x,*y);
            grErr->SetPointError(pointcnt,*xerr,*yerr);
            cout << "The value to plot " << *x << " (" << *xerr << "), " << *y << " (" << *yerr << ")" << endl;

	    pointcnt++;
	 }
      }
      else
      {
         gr->SetPoint(infilenr,*x,*y);
         cout << "The value to plot " << *x << ", " << *y << endl;
      }// end of plot selection

      // Saving the final plot, when we finish with the last file
      if( infilenr+1 == evtcount )
      {
         if( (directivetype[seldir[0]] == 'E') || (directivetype[seldir[1]] == 'E') )
	 {
            // Removing any unused points
	    if(infilenr >= pointcnt)
	    {
	       if((directive[seldir[0]] == "shfoot") || (directive[seldir[1]] == "shfoot"))
	       {
	          *itemp = grAsymmErr->GetN();
		  cout << "All points = " << *itemp << ", Actual points = " << pointcnt << ": ";
	          for(int h = pointcnt; h < *itemp; h++)
	             grAsymmErr->RemovePoint(pointcnt);
	       }
	       else
	       {
	          *itemp = grErr->GetN();
		  cout << "All points = " << *itemp << ", Actual points = " << pointcnt << ": ";
	          for(int h = pointcnt; h < *itemp; h++)
	             grErr->RemovePoint(pointcnt);
	       }
	       cout << "Removing excess points." << endl;
	    }

	    if((directive[seldir[0]] == "shfoot") || (directive[seldir[1]] == "shfoot"))
	    {
               grAsymmErr->GetXaxis()->SetRangeUser(xrange[0]-0.1*(xrange[1]-xrange[0]),xrange[1]+0.1*(xrange[1]-xrange[0]));
               grAsymmErr->GetYaxis()->SetRangeUser(yrange[0]-0.1*(yrange[1]-yrange[0]),yrange[1]+0.1*(yrange[1]-yrange[0]));
	       if(directive[seldir[0]] == "shfoot")
                  *stemp = ";Depth of " + IntToStr((int)(shfootlimit*100)) + directivedesc[seldir[0]] + ";" + directivedesc[seldir[1]];
	       else if(directive[seldir[1]] == "shfoot")
                  *stemp = ";" + directivedesc[seldir[0]] + ";Depth of " + IntToStr((int)(shfootlimit*100)) + directivedesc[seldir[1]];
               grAsymmErr->SetTitle((*stemp).c_str());
	       grAsymmErr->SetMarkerStyle(20);
	       grAsymmErr->SetLineColor(kRed-5);
	       grAsymmErr->SetMarkerColor(kRed);
	       grAsymmErr->SetMarkerSize(0.8);
               grAsymmErr->Draw("AP");
	    }
	    else
	    {
               grErr->GetXaxis()->SetRangeUser(xrange[0]-0.1*(xrange[1]-xrange[0]),xrange[1]+0.1*(xrange[1]-xrange[0]));
               grErr->GetYaxis()->SetRangeUser(yrange[0]-0.1*(yrange[1]-yrange[0]),yrange[1]+0.1*(yrange[1]-yrange[0]));
               *stemp = ";" + directivedesc[seldir[0]] + ";" + directivedesc[seldir[1]];
               grErr->SetTitle((*stemp).c_str());
	       grErr->SetMarkerStyle(20);
	       grErr->SetLineColor(kRed-5);
	       grErr->SetMarkerColor(kRed);
	       grErr->SetMarkerSize(0.8);
               grErr->Draw("AP");
	    }
	 }
	 else
	 {
            gr->GetXaxis()->SetRangeUser(xrange[0]-0.1*(xrange[1]-xrange[0]),xrange[1]+0.1*(xrange[1]-xrange[0]));
            gr->GetYaxis()->SetRangeUser(yrange[0]-0.1*(yrange[1]-yrange[0]),yrange[1]+0.1*(yrange[1]-yrange[0]));
            *stemp = ";" + directivedesc[seldir[0]] + ";" + directivedesc[seldir[1]];
            gr->SetTitle((*stemp).c_str());
	    gr->SetMarkerStyle(20);
	    gr->SetMarkerColor(kRed);
	    gr->SetMarkerSize(0.6);
            gr->Draw("AP");
	 }
      
         *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_vs_" + directive[seldir[1]] + ".pdf";
         c1->SaveAs((*stemp).c_str());

         // Additional .C export sequence
         *stemp = string(BASEDIR) + "/results/" + directive[seldir[0]] + "_vs_" + directive[seldir[1]] + ".C";
         c1->SaveAs((*stemp).c_str());
      }// end of saving final plot
   }
   else
      return -1;
/*------------------------------------------------------------------------------------------*/

   delete x;
   delete xerr;
   delete y;
   delete yerr;
   delete nrpoints;
   delete stemp;
   delete itemp;
   delete dtemp;
   delete[] derr;
   delete cordist;

   return 0;
}

int AdstAnalyseTool::GetData(std::string *sepdir, int *seldir, std::string type, double *out, double *outerr, int *count, int intanknr)
{
   cout << "# Entering function AdstAnalyseTool::GetData()..." << endl;

   int *useddir;
   useddir = new int;

   // Save the relevant seldir to the useddir
   if(type == "Xaxis")
      *useddir = seldir[0];
   else if(type == "Yaxis")
      *useddir = seldir[1];
   else if(type == "Zaxis")
      *useddir = seldir[2];
   else
      return -1;

   double *x, *xerr, *dtemp1, *dtemp2, *dtemp3;
   int *itemp;

   double errvalue[2];

   x = new double;
   xerr = new double;
   dtemp1 = new double;
   dtemp2 = new double;
   dtemp3 = new double;
   itemp = new int;

   SdRecShower &sdrecshw = fRecEvent->GetSDEvent().GetSdRecShower();
   GenShower &genshw = fRecEvent->GetGenShower();

   *itemp = 0;

   if( (directive[*useddir] == "totalvemmax") )	// Maximum VEM signal
   {
      *x = acttanks[intanknr].GetTotalSignal();
      *xerr = acttanks[intanknr].GetTotalSignalError();

      if( *x > *out )
      {
         *out = *x;
	 *outerr = *xerr;
      }
      *count = 1;
   }
   else if( (directive[*useddir] == "xmaxqual") ) // Longest FD shower track (best Xmax estimation)
   {
      *x = acteyes[GetEyeLongestTrack()].GetXmax();
      *xerr = acteyes[GetEyeLongestTrack()].GetXmaxError();

      *out = *x;
      *outerr = *xerr;
   }
   else if( (directive[*useddir] == "totalvemaver") || (directive[*useddir] == "muoncountaver") || (directive[*useddir] == "distaver") || (directive[*useddir] == "fdenergy") || (directive[*useddir] == "fdcalenergy") || (directive[*useddir] == "xmax") )	// Average VEM signal, muon count, distance, energy, xmax
   {
      if(directive[*useddir] == "totalvemaver")
      {
         *x = acttanks[intanknr].GetTotalSignal();
         *xerr = acttanks[intanknr].GetTotalSignalError();
      }
      else if(directive[*useddir] == "muoncountaver")
      {
         *x = genshw.GetMuonNumber();
	 *xerr = 1.e-7;
      }
      else if(directive[*useddir] == "distaver")
      {
         *x = acttanks[intanknr].GetSPDistance();
         *xerr = acttanks[intanknr].GetSPDistanceError();
      }
      else if(directive[*useddir] == "fdenergy")
      {
         *x = acteyes[GetEyeLongestTrack()].GetEnergy();
         *xerr = acteyes[GetEyeLongestTrack()].GetEnergyError();
      }
      else if(directive[*useddir] == "fdcalenergy")
      {
         *x = acteyes[GetEyeLongestTrack()].GetEcal();
         *xerr = acteyes[GetEyeLongestTrack()].GetEcalError();
      }
      else if(directive[*useddir] == "xmax")
      {
         *x = acteyes[GetEyeLongestTrack()].GetXmax();
         *xerr = acteyes[GetEyeLongestTrack()].GetXmaxError();
      }

      *out += *x;
      *outerr += *xerr;

      *count += 1;
   }
   else if( (directive[*useddir] == "muoncount") || (directive[*useddir] == "totalvemsum") )		// Total muon count, VEM sum
   {
      if(directive[*useddir] == "totalvemsum")
      {
         *x = acttanks[intanknr].GetTotalSignal();
         *xerr = acttanks[intanknr].GetTotalSignalError();
      }
      else if(directive[*useddir] == "muoncount")
      {
         *x = genshw.GetMuonNumber();
	 *xerr = 1.e-7;
      }

      *out += *x;
      *outerr = *xerr;

      *count = 1;
   }
   else if( (directive[*useddir] == "actstation") || (directive[*useddir] == "sdenergy") || (directive[*useddir] == "S1000") || (directive[*useddir] == "curvature") || (directive[*useddir] == "acteyes") )
   {
      if(directive[*useddir] == "actstation")
      {
         *x = nracttanks;
         *xerr = 1.e-7;
      }
      else if(directive[*useddir] == "sdenergy")
      {
         *x = sdrecshw.GetEnergy();
	 *xerr = sdrecshw.GetEnergyError();
      }
      else if(directive[*useddir] == "S1000")
      {
         *x = sdrecshw.GetS1000();
	 *xerr = sdrecshw.GetS1000Error();
      }
      else if(directive[*useddir] == "curvature")
      {
         *x = sdrecshw.GetCurvature();
	 *xerr = sdrecshw.GetCurvatureError();
      }
      else if(directive[*useddir] == "acteyes")
      {
         cout << "Give me eye!" << endl;
	 for(int i = 0; i < acteyes.size(); i++)
	 {
	    if(eyevalid[i] == 1)
	       *itemp = *itemp + 1;
	 }
         *x = (double)*itemp;
	 *xerr = 1.e-7;
      }

      *out = *x;
      *outerr = *xerr;
      *count = 1;
   }
   else if( (directive[*useddir] == "shfoot") )	// Shower foot calculation
   {
//cout << "Entries: " << seltree->GetEntries() << endl;

      *x = 0;
      *xerr = 0;

      vector<FDEvent> fdevt = fRecEvent->GetFDEvents();
    
      vector<double> slantDepth;
      vector<double> profiledEdX;
      vector<double> profiledEdXerr;

      slantDepth = acteyes[GetEyeLongestTrack()].GetDepth();
      profiledEdX = acteyes[GetEyeLongestTrack()].GetEnergyDeposit();
      profiledEdXerr = acteyes[GetEyeLongestTrack()].GetEnergyDepositError();

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
      }

      *itemp = 0;

//cout << "Entries under limit:" << "\t";
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
//	    cout << "Izr.1a: P1 = (" << xfoot[i-1] << "," << yfoot[i-1] << "), P2 = (" << xfoot[i] << "," << yfoot[i] << ")" << endl;

	    *dtemp1 = (xfoot[i] - xfoot[i-1])/(yfoot[i] - yfoot[i-1]); // 1/k = (x2 - x1)/(y2 - y1)
	    *dtemp2 = shfootlimit*(yfoot[yfoot.size()-1]) - yfoot[i]; // y - y2
            *x = (*dtemp1)*(*dtemp2) + xfoot[i]; // x = (1/k)*(y - y2) + x2

//	    cout << "Izr.1b: 1/k = " << *dtemp1 << ", max = " << yfoot[yfoot.size()-1] << ", y - y2 = " << *dtemp2 << ", x = " << *x << endl;
//	    cout << "Izr.2a: P1err = (" << xfoot[i-1] << "," << yerrfoot[i-1] << "), P2err = (" << xfoot[i] << "," << yerrfoot[i] << ")" << endl;

            *dtemp1 = (xfoot[i] - xfoot[i-1])/(yfoot[i]+yerrfoot[i] - (yfoot[i-1]+yerrfoot[i-1])); // 1/kerr = (x2 - x1)/(y2err - y1err)
	    *dtemp2 = (yfoot[i]+yerrfoot[i]) - (1/(*dtemp1))*(xfoot[i]); // aerr = y2err - (y2err - y1err)/(x2 - x1)*x2
	    *dtemp3 = (1/(*dtemp1))*(*x) + (*dtemp2); // yerr = kerr*x + aerr
	    *xerr = (*dtemp3) - (((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(*x) + (yfoot[i] - ((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(xfoot[i]))); // Dy = yerr - y

//	    cout << "Izr.2b: 1/kerr = " << *dtemp1 << ", a = " << *dtemp2 << ", yerr = " << *dtemp3 << ", Dy = " << *xerr << endl;

	    *dtemp1 = ((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(*x) + (yfoot[i] - ((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(xfoot[i])); // y = k*x + a

            for(int j = i; ; j++)
	    {
	       if( yfoot[j] >= (*dtemp1)+(*xerr) )
	       {
        	  *dtemp2 = xfoot[j];
//                cout << "Upper limit is = " << xfoot[j] << endl;
     	     break;
               }
            }
            for(int j = i; ; j--)
            {
               if(j == 0)
               {
	          *dtemp3 = xfoot[j];
//	          cout << "Lower limit is = " << xfoot[j] << endl;
     	          break;
               }

	       if( yfoot[j] <= (*dtemp1)-(*xerr) )
	       {
	          *dtemp3 = xfoot[j];
//	          cout << "Lower limit is = " << xfoot[j] << endl;
	          break;
	       }
	    }
//	    *dtemp2 = ((xfoot[i] - xfoot[i-1])/(yfoot[i] - yfoot[i-1]))*((*dtemp1 + (*xerr)) - (yfoot[i] - ((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(xfoot[i]))); // x+Dx = ((y+Dy - a)/k
//	    *dtemp3 = ((xfoot[i] - xfoot[i-1])/(yfoot[i] - yfoot[i-1]))*((*dtemp1 - (*xerr)) - (yfoot[i] - ((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(xfoot[i])));// x-Dx = ((y-Dy - a)/k
//	    *xerr = (*dtemp2) - (*x); // Dx = x+Dx - x

	    cout << "Izr.2c: x = " << *x << ", y = " << *dtemp1 << ", x+Dx = " << *dtemp2 << ", x-Dx = " << *dtemp3 << ", Dx+ = " << (*dtemp2)-(*x) << ", Dx- = " << (*x)-(*dtemp3) << endl;

//            if( ((xmaxlimit != 0.) && (xmaxerr[*count] <= xmaxlimit)) || (xmaxlimit == 0.) )
//            {
               *out += *x;//xfoot[i-1];
               outerr[0] += (*x)-(*dtemp3);
               outerr[1] += (*dtemp2)-(*x);
/*            }
            else if( (xmaxlimit != 0.) && (xmaxerr[*count] > xmaxlimit) )
            {
               cout << "Shower with " << *x << " (+" << (*dtemp2)-(*x) << "-" << (*x)-(*dtemp3) << ") has too large Xmax error (" << xmaxerr[*count] << ")" << endl;
	       *out += 0;
	       outerr[0] += 0;
	       outerr[1] += 0;
      	       aboveerr--;
            }
//cout << i-1 << ": " << xfoot[i-1] << "\t" << yfoot[i-1] << "\t" << (yfoot[i-1]/(yfoot[yfoot.size()-1]))*100 << endl;
//cout << i-1 << ": " << *x << "(" << *xerr << ")" << "\t" << *dtemp1 << "\t" << (*dtemp1/(yfoot[yfoot.size()-1]))*100 << endl;
 */           break;
         }
      }

      *count += 1;
   }
   else if( (directive[*useddir] == "risetime") || (directive[*useddir] == "risetimeearly") || (directive[*useddir] == "risetimelate") )	// Risetime calculation
   {
      if( (directive[*useddir] == "risetimeearly") && (sdhotnsec > acttanks[intanknr].GetTimeNSecond()) )
      {
         *x = acttanks[intanknr].GetRiseTime();
         *xerr = acttanks[intanknr].GetRiseTimeRMS();

         cout << "Risetime early: " << *x << " (" << *xerr << ")" << endl;
    
         *out += *x;
         *outerr += *xerr;
    
         *count += 1;
      }
      else if( (directive[*useddir] == "risetimelate") && (sdhotnsec < acttanks[intanknr].GetTimeNSecond()) )
      {
         *x = acttanks[intanknr].GetRiseTime();
         *xerr = acttanks[intanknr].GetRiseTimeRMS();

         cout << "Risetime late: " << *x << " (" << *xerr << ")" << endl;
    
         *out += *x;
         *outerr += *xerr;
    
         *count += 1;
      }
      else if( directive[*useddir] == "risetime" )
      {
         *x = acttanks[intanknr].GetRiseTime();
         *xerr = acttanks[intanknr].GetRiseTimeRMS();

         cout << "Risetime: " << *x << " (" << *xerr << ")" << endl;
    
         *out += *x;
         *outerr += *xerr;
    
         *count += 1;
      }
   }

   delete x;
   delete xerr;
   delete dtemp1;
   delete dtemp2;
   delete dtemp3;
   delete itemp;

   return 0;
}

// Check for the eye with the longest FD track
int AdstAnalyseTool::GetEyeLongestTrack()
{
   cout << "# Entering function AdstAnalyseTool::GetEyeLongestTrack()..." << endl;
   double longtrack = 0.0;
   int itemp;

   itemp = -1;

   for(int i = 0; i < acteyes.size(); i++)
   {
      if( acteyes[i].GetTrackLength() > longtrack )
      {
         itemp = i;
	 longtrack = acteyes[i].GetTrackLength();
      }
   }
   
//   cout << "The longest track is " << longtrack << endl;

   return itemp;
}

// TODO: Check if it works!
void AdstAnalyseTool::SetRanges2D(int *seldir, double *inx, double *inxerr, double *iny, double *inyerr, int axis)
{
   if(axis == 0)
   {
      if( ((directiveaffil[seldir[axis]] == "showeyedata") || (directiveaffil[seldir[axis]] == "eyelong")) && (*inx != 0.) )
      {
         if( directivetype[seldir[axis]] == 'E' )
         {
            if( (*inx)-(*inxerr) < xrange[0] )
               xrange[0] = (*inx)-(*inxerr);
            if( (*inx)+(*inxerr) > xrange[1] )
               xrange[1] = (*inx)+(*inxerr);
         }
         else if( directivetype[seldir[axis]] == 'S' )
         {
            if( (*inx)-(inxerr[0]) < xrange[0] )
               xrange[0] = (*inx)-(inxerr[0]);
            if( (*inx)+(inxerr[1]) > xrange[1] )
               xrange[1] = (*inx)+(inxerr[1]);
         }
         else
         {
            if( *inx < xrange[0] )
               xrange[0] = *inx;
            if( *inx > xrange[1] )
               xrange[1] = *inx;
         }
      }
      else if( ((directiveaffil[seldir[axis]] == "showeyedata") || (directiveaffil[seldir[axis]] == "eyelong")) && (*inx == 0.) )
         cout << "Invalid point. Not counting towards x range." << endl;
      else
      {
         if( directivetype[seldir[axis]] == 'E' )
         {
            if( (*inx)-(*inxerr) < xrange[0] )
               xrange[0] = (*inx)-(*inxerr);
            if( (*inx)+(*inxerr) > xrange[1] )
               xrange[1] = (*inx)+(*inxerr);
         }
         else if( directivetype[seldir[axis]] == 'S' )
         {
            if( (*inx)-(inxerr[0]) < xrange[0] )
               xrange[0] = (*inx)-(inxerr[0]);
            if( (*inx)+(inxerr[1]) > xrange[1] )
               xrange[1] = (*inx)+(inxerr[1]);
         }
         else
         {
            if( *inx < xrange[0] )
               xrange[0] = *inx;
            if( *inx > xrange[1] )
               xrange[1] = *inx;
         }
      }
   }
   else
   {
      if( ((directiveaffil[seldir[axis]] == "showeyedata") || (directiveaffil[seldir[axis]] == "eyelong")) && (*iny != 0.) )
      {
         if( directivetype[seldir[axis]] == 'E' )
         {
            if( (*iny)-(*inyerr) < yrange[0] )
               yrange[0] = (*iny)-(*inyerr);
            if( (*iny)+(*inyerr) > yrange[1] )
               yrange[1] = (*iny)+(*inyerr);
         }
         else if( directivetype[seldir[axis]] == 'S' )
         {
            if( (*iny)-(inyerr[0]) < yrange[0] )
               yrange[0] = (*iny)-(inyerr[0]);
            if( (*iny)+(inyerr[1]) > yrange[1] )
               yrange[1] = (*iny)+(inyerr[1]);
         }
         else
         {
            if( *iny < yrange[0] )
               yrange[0] = *iny;
            if( *iny > yrange[1] )
               yrange[1] = *iny;
         }
      }
      else if( ((directiveaffil[seldir[axis]] == "showeyedata") || (directiveaffil[seldir[axis]] == "eyelong")) && (*iny == 0.) )
         cout << "Invalid point. Not counting towards y range." << endl;
      else
      {
         if( directivetype[seldir[axis]] == 'E' )
         {
            if( (*iny)-(*inyerr) < yrange[0] )
               yrange[0] = (*iny)-(*inyerr);
            if( (*iny)+(*inyerr) > yrange[1] )
               yrange[1] = (*iny)+(*inyerr);
         }
         else if( directivetype[seldir[axis]] == 'S' )
         {
            if( (*iny)-(inyerr[0]) < yrange[0] )
               yrange[0] = (*iny)-(inyerr[0]);
            if( (*iny)+(inyerr[1]) > yrange[1] )
               yrange[1] = (*iny)+(inyerr[1]);
         }
         else
         {
            if( *iny < yrange[0] )
               yrange[0] = *iny;
            if( *iny > yrange[1] )
               yrange[1] = *iny;
         }
      }
   }

   cout << "Ranges are: (" << xrange[0] << ", " << xrange[1] << ") and (" << yrange[0] << ", " << yrange[1] << ")" << endl;
}
