#ifndef _ADSTMVA_H_
#define _ADSTMVA_H_

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
#include "TKey.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "TLine.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/Config.h"

#include "OfflineInclude.h"

#define SUBPLOT true
#define ALLEYES 5

// Class holding the root file structure
class Observables {
public:
   Observables();
   virtual ~Observables();

   float xmax;		// depth of shower maximum
   float x0;		// GH shape parameter X0 (0 is at first interaction point
   float lambda;	// GH shape parameter lambda
   float shfoot;	// Shower foot
   float fdenergy;	// Total energy from FD
   float nrmu;		// number of muons at ground level
   float shwsize;	// Shower size - replacement for S1000
   float ldfbeta;	// beta variable of the LDF
   float curvature;	// curvature R
   float risetime;	// risetime
};

class AdstMva {
public:
   // Constructor and destructor
   AdstMva();
   virtual ~AdstMva();

   // Rewrite
   void RewriteObservables(int innr, Observables sig, Observables back, TTree *back_tree);
   void PrepareOtherTrees(int nr, int sig);
   void CreateMVAPlots(double cut, std::vector<std::string> obs);
//   void SetPlotRange(TH1F *hist1, TH1F *hist2, double *xval, double *yval, double xrangeboost, double yrangeboost);
   void SetPlotRange(TNtuple *sig, TNtuple *back, TH1F *base, double *val, double rangeboost, std::string obs, char xory);
   void SetupAxis(TH1F *hist, std::string obs);

   // Reconstructed event variables (from the ADST format)
   RecEvent *fRecEvent;
   RecEventFile *fFile;
   DetectorGeometry *fDetGeo;

   GenShower *genshw;
#ifdef OFFLINENEW
   UnivRecShower *unishw;
#endif
   SdRecShower *sdrecshw;

   // All input filenames that will be used for the MVA
   std::vector<std::string> inname;
   // Output file created for the MVA
   TFile *outfile;
   TTree *sig_tree;
//   TTree *back_tree;
   TTree *all_tree;
   std::string outname;
   bool graphical;

   // Active eyes
   std::vector<FdRecShower> acteyes;
   int GetEyeLongestTrack();

   int GetShowerFoot(int longestEye, std::vector<FDEvent> fdevt);
   double shfootlimit;
   double shfoot;

   // Boolean to say, if the reconstruction finished correctly
   bool goodrec;
};

#endif
