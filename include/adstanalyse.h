#ifndef _ADSTANALYSE_H_
#define _ADSTANALYSE_H_

#include <string>
#include <vector>
#include <fstream>

#include "TROOT.h"
#include "TDirectory.h"
#include "TSystem.h"

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
#include "TMath.h"
#include "TColor.h"

#include "OfflineInclude.h"

#define SUBPLOT true
#define ALLEYES 5

class AdstAnalyseTool {
public:
   // Constructor and destructor
   AdstAnalyseTool();
   virtual ~AdstAnalyseTool();

   // Determines if starting file is .root or .tar.gz, number of events/files
   int argtype;
   int evtcount;

   // Reconstructed event variables (from the ADST format)
   RecEvent *fRecEvent;
   RecEventFile *fFile;
   DetectorGeometry *fDetGeo;

   // Hello and instructions function
   void ReadOption();
   int readopt;
   int readevent;
   void Hello();

   // Output file for additional data analysis
   ofstream outdata;

   // Analysis type
   int analysistype;	// Type of analysis to perform (observables, profiles, more)
   int filesel;		// Single or multiple files for input
   std::string plotdirective;	// Directive of what to plot
   int nrdirs;		// Number of directives (to know what kind of plot we asked for)

   // Default ROOT objects
   TCanvas *c1;
   TGraph *gr, *grtemp;
   TGraphErrors *grErr;
   TGraphAsymmErrors *grAsymmErr;
   TGraphErrors *grErrtemp;
   TH1F *histf;
   TH2F *hist2f;

   // Functions that use and access private valriables
   void SetAnalysisFilename(char *file);
   std::string GetAnalysisFilename();

   // Retrieve the active SD tanks
   std::string GetActiveTanks();
   int nracttanks;
   std::vector<SdRecStation> acttanks;
   // Retrieve the active FD eyes
   void GetActiveEyes();
   std::string GetActiveEyeList();
   int nracteyes;
   std::vector<FdRecShower> acteyes;
   FdRecShower besteye;
   int eyevalid[ALLEYES];

   // Run the analysis directive
   int RunDirective(int infilenr);

   // Specific plot commands
   int PlotLDF(std::string *sepdir, int *seldir);
   int PlotVEM(std::string *sepdir, int *seldir);
   int PlotLongProf(std::string *sepdir, int *seldir);

   // Using directives for multiple input files
   int PrepareDirective(std::string *sepdir, int *seldir, int infilenr);

   // Sub-directives for multiple file plotting
//   int GetRisetime(int i, double nsec, int *seldir, int infilenr, int *itemp, double *output);
   int GetData(std::string *sepdir, int *seldir, std::string type, double *out, double *outerr, int *count, int intanknr);
   int GetEyeLongestTrack();

   // Database directive information
   const std::string directive[256] = {		// directive names
      "stationid","totalvem","starttime.nsec","distance","muoncount",
      "actstation","totalvemmax","totalvemaver","totalvemsum","risetime","risetimeearly","risetimelate","sdenergy","S1000","curvature","distaver","muoncountaver",
      "eyeid","energy","emenergy","xmax",
      "acteyes","fdenergy","fdcalenergy","xmaxqual","shfoot",
      "ldf","vemsig","longprof",
      "all"
   };

   const std::string directiveaffil[256] = {	// ROOT file tree names for directives
      "showtankdata","showtankdata","showtankdata","showtankdata","showtankdata",
      "showsrecdata","showtankdata","showtankdata","showtankdata","tankvem","tankvem","tankvem","showsrecdata","showsrecdata","showsrecdata","showtankdata","showtankdata",
      "showfrecdata","showeyedata","showeyedata","showeyedata",
      "showfrecdata","showeyedata","showeyedata","showeyedata","eyelong",
      "srecldfdata","tankvem","eyelong",
      "none"
   };

   const char directivetype[256] = {	// ROOT file datatypes for the selected directive (I = Integer, D = double, E = double with errors, S = structure for profiles)
      'I','D','D','D','I',
      'I','D','D','D','D','D','D','E','E','E','E','I',
      'I','D','D','E',
      'I','E','E','E','S',
      'D','D','D',
      'N'
   };

   const std::string directivedesc[256] = {	// decription of the directives to use for X and Y axis naming on plots
      "SD tank ID","Total VEM signal","Start time (ns)","Distance from shower core (m)","Simulated muon count",
      "Active stations","Total VEM signal from max. SD tank","Total VEM signal average","Total VEM signal sum from all SD tanks","Risetime average (ns)","Risetime early average (ns)","Risetime late average (ns)","SD reconstructed energy (eV)","S1000 (LDF value at 1000m)","Shower front curvature","Distance from shower core average (m)","Average simulated muon count",
      "FD eye ID","Total shower energy (eV)","EM shower energy (eV)","X_{max} (g/cm^{2})",
      "Active eyes","Total shower energy average (eV)","EM shower energy average (eV)","X_{max} (g/cm^{2})","\% shower foot (g/cm^{2})",
      "LDF profile","Tank VEM signal trace","Longitudinal profile",
      "all"
   };

   const int directiveavail[256] = {	// availability for the directive to use it in single file mode (0), multiple file mode (1) or both (2)
      0,0,0,0,2,
      1,1,1,1,1,1,1,1,1,1,1,1,
      0,0,0,2,
      1,1,1,1,1,
      0,0,0,
      0
   };

   const std::string directivetrees[256] = {		// relevant tree names connected to the directive
      "stationid","totalvem","starttime.nsec","distance","muoncount",
      "actstation","totalvem","totalvem","totalvem","vemsig","vemsig","vemsig","energy","S1000","curvature","distance","muoncount",
      "eyeid","energy","emenergy","xmax",
      "acteyes","energy","emenergy","xmax","longprof",
      "ldf","vemsig","longprof",
      "all"
   };

   // Saved collection of filenames (multiple input file mode)
   std::vector<std::string> tarnames;

   // Vectors for saving x and y values, Xmax values and shower foot values
   std::vector<double> xvalue;
   std::vector<double> yvalue;

   std::vector<double> xmaxerr;

   std::vector<double> xfoot;
   std::vector<double> yfoot;
   std::vector<double> yerrfoot;

   // Ranges for plots
   double xrange[2];
   double yrange[2];
   double zrange[2];
   void InitRanges();

   void SetRanges2D(int *seldir, double *inx, double *inxerr, double *iny, double *inyerr, int axis);

   // Limit for Xmax error and shower foot values
   double xmaxlimit;
   int slantset;
   double shfootlimit;

   // Point counter for plot
   int pointcnt;

   // temporary bool
   bool btemp;
private:
   // Saved analysis filename
   std::string adstanalysisfile;
};

#endif
